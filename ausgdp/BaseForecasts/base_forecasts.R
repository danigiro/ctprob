rm(list = ls(all = TRUE))
library(forecast)
library(progress)

load("./DATA.RData")

# arima or ets
model <- "arima"
# log or lev
trans <- "lev"

dir.create(file.path(".","BaseForecasts",model, trans), recursive = TRUE, showWarnings = FALSE)

end_traing <- NROW(DATA)-H-fixed_length+1
#Test
#end_traing = 2
#DATA = DATA[,1:3]
start_year <- min(time(DATA))
for(j in 0:(end_traing-1)){
  listFit <- list()
  listRes <- list()
  listBase <- list()
  listTest <- list()
  listScale <- list()
  
  pb <- progress_bar$new(format = paste0("Replication n. ", j+1, 
                                         " [:bar] :percent in :elapsed (ETA: :eta)"),
                         total = length(K), clear = FALSE, width= 60, show_after = 0)
  for(k in K){
    idK <- paste0("k", k)
    if(k == 1){
      train1 <- window(DATA, start=c(start_year,1), end=c(start_year,j+fixed_length))
      cond <- NROW(train1)%%m
      if(cond !=0){
        train1 <- rbind(matrix(NA, ncol=NCOL(DATA), nrow = m-NROW(train1)%%m), train1)
      }
      test1 <- window(DATA, start=c(start_year,j+fixed_length+1), end=c(start_year,j+fixed_length+H))
    }
    train <- ts(apply(train1, 2, function(x) rowSums(matrix(x,ncol = k,byrow = T))), 
                start=c(start_year, (j+1)), frequency = m/k)
    listTest[[idK]] <- ts(apply(test1, 2, function(x) rowSums(matrix(x,ncol = k,byrow = T))), 
                                     start=start_year+((j+fixed_length)/m), frequency = m/k)
    for(i in 1:NCOL(DATA)){
      idS <- colnames(DATA)[i]
      ts_train <- na.omit(train[,i])
      
      if(trans == "log"){
        ts_train[ts_train==0] <- min(ts_train[ts_train!=0])/2
      }

      if(model == "arima"){
        if(trans == "lev"){
          listFit[[idK]][[idS]] <- auto.arima(ts_train)
        }else if(trans == "log"){
          listFit[[idK]][[idS]] <- auto.arima(ts_train, lambda = 0, biasadj = TRUE)
        }else{
          stop("Transformation not implemented")
        }
      }else if(model == "ets"){
        if(trans == "lev"){
          listFit[[idK]][[idS]] <- ets(ts_train)
        }else if(trans == "log"){
          listFit[[idK]][[idS]] <- ets(ts_train, lambda = 0, biasadj = TRUE)
        }else{
          stop("Transformation not implemented")
        }
      }else{
        stop("Model not implemented")
      }
      listBase[[idK]][[idS]] <- unname(forecast(listFit[[idK]][[idS]], h = H/k)$mean)
      listScale[[idK]][[idS]] <- mean(abs(diff(ts_train, lag = m/k, differences = 1)), na.rm = TRUE)
      
      if(cond !=0){
        listRes[[idK]][[idS]] <- c(rep(NA,(NROW(train1)/k)-length(ts_train)), as.vector(ts_train - fitted(listFit[[idK]][[idS]])))
      }else{
        listRes[[idK]][[idS]] <- as.vector(ts_train - fitted(listFit[[idK]][[idS]]))
      }
    }
    listBase[[idK]] <- simplify2array(listBase[[idK]])
    listRes[[idK]] <- simplify2array(listRes[[idK]])
    listScale[[idK]] <- simplify2array(listScale[[idK]])
    pb$tick()
  }
  itername <- paste0("ite_", formatC(j+1, width = nchar(end_traing+1), 
                                     format = "d", flag = "0"), ".RData")
  save(listBase, listRes, listFit, listScale, listTest,
       file = file.path(".","BaseForecasts",model, trans, itername))
}
rm(list = ls(all = TRUE))

rm(list = ls(all = TRUE))
library(forecast)
library(progress)
library(Matrix)
library(FoReco)

load("./DATA.RData")

# arima ets
model <- "arima"
# log lev
trans <- "lev"

dir.create(file.path(".","OverlapRes",model, trans), recursive = TRUE, showWarnings = FALSE)

listFiles <- sort(list.files(file.path(".","BaseForecasts", model, trans), 
                             full.names = TRUE))
end_traing <- NROW(DATA)-H-fixed_length+1
#Test
#end_traing = 2
#DATA = DATA[,1:3]
start_year <- min(time(DATA))
for(j in 0:(end_traing-1)){
  load(listFiles[j+1])
  listResOver <- list()
  
  pb <- progress_bar$new(format = paste0("Replication n. ", j+1, 
                                         " [:bar] :percent in :elapsed (ETA: :eta)"),
                         total = length(K), clear = FALSE, width= 60, show_after = 0)
  for(k in K){
    idK <- paste0("k", k)
    if(k == 1){
      train1 <- window(DATA, start=c(start_year,1), end=c(start_year,j+fixed_length))
    }
    train <- apply(train1, 2, zoo::rollsum, k = k, align = "left", fill = c(NULL, NULL, NA))
    for(i in 1:NCOL(DATA)){
      ts_train <- train[,i]
      idS <- colnames(DATA)[i]
      if(k == 1){
        tmp_res <- ts_train-fitted(listFit[[idK]][[idS]])
        listResOver[[idK]][[idS]] <- as.vector(t(sapply(1:(H/k), function(x){
          out <- rep(NA, length(tmp_res))
          id <- x:(length(tmp_res)-((H/k)-x))
          out[1:length(id)] <- tmp_res[id]
          out
        })))
      }else{
        opt <- suppressWarnings(matrix(ts_train, nrow = k))
        opt[-c(1:length(ts_train))] <- NA
        opt_true <- opt
        
        if(trans == "log"){
          opt <- t(apply(opt, 1, function(x){
            x[x<=0 & !is.na(x)] <- min(x[x>0 & !is.na(x)], na.rm = TRUE)/2
            x
          }))
        }
        tsinfo <- tsp(listFit[[idK]][[idS]]$x)
        res_out <- as.vector(do.call(rbind, lapply(1:NROW(opt), function(idr){
          idna <- is.na(opt[idr,])
          if(model == "arima"){
            mod <- Arima(ts(opt[idr, !idna], start = tsinfo[1], frequency = tsinfo[3]), 
                         model = listFit[[idK]][[idS]])
          }else if(model == "ets"){
            mod <- ets(ts(opt[idr, !idna], start = tsinfo[1], frequency = tsinfo[3]), 
                       model = listFit[[idK]][[idS]], use.initial.values=TRUE)
          }
          
          true_val <- opt_true[idr, !idna]
          tmp_res <- rep(NA, length(idna))
          tmp_res[!idna] <- true_val - fitted(mod)
          t(sapply(1:(H/k), function(x){
            out <- rep(NA, length(tmp_res))
            #id <- x:(length(tmp_res)-((H/k)-x))
            id <- x:(length(tmp_res))
            out[1:length(id)] <- tmp_res[id]
            out
          }))
        })))
        
        listResOver[[idK]][[idS]] <- res_out[1:(length(ts_train)*(m/k))]
      }
    }
    listResOver[[idK]] <- simplify2array(listResOver[[idK]])
    pb$tick()
  }
  allres <- t(do.call(rbind, listResOver[paste0("k", sort(K, decreasing = TRUE))]))
  
  N <- NCOL(allres) / FoReco::thf_tools(m = K)$kt
  DN <- FoReco:::Dmat(h = N, m = K, n = NROW(allres))
  E <- matrix(DN %*% as.vector(t(allres)), nrow = N, byrow = TRUE)
  E <- na.omit(E)  
  
  Nna <- NROW(E)
  DN <- FoReco:::Dmat(h = Nna, m = K, n = NROW(allres))
  E <- matrix(t(DN) %*% as.vector(t(E)), nrow = NROW(allres), byrow = TRUE)
  
  listRes <- lapply(split(t(E), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                    Nna*m/sort(K, decreasing = TRUE))), 
                    function(x){
                      out <- matrix(x, ncol = sum(dim(C)))
                      colnames(out) <- rownames(allres)
                      out
                    })[paste0("k", sort(K))]
  
  itername <- paste0("ite_", formatC(j+1, width = nchar(end_traing+1), 
                                     format = "d", flag = "0"), ".RData")
  save(listRes,
       file = file.path(".","OverlapRes",model, trans, itername))
}
rm(list = ls(all = TRUE))
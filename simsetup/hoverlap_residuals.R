rm(list = ls(all = TRUE))
library(forecast)
library(progress)
library(Matrix)
library(FoReco)

load("./DATA.RData")
source("./R/hfitted.R")

dir.create(file.path(".","HOverlapRes"), recursive = TRUE, showWarnings = FALSE)

listFiles <- sort(list.files(file.path(".","BaseForecasts"), 
                             full.names = TRUE))
#end_traing <- (NROW(DATA)/m)-H/12-(fixed_length/m)+1 # 10 years

pb <- progress_bar$new(format = " [:bar] :percent in :elapsed (ETA: :eta)",
                       total = length(K)*sum(dim(C))*length(listFiles), 
                       clear = FALSE, width= 60, show_after = 0)
for(j in 1:length(listFiles)){
  load(listFiles[j])
  listResOver <- list()
  for(k in K){
    idK <- paste0("k", k)
    if(k == 1){
      train <-listTrain[[idK]]
    }else{
      train <- apply(listTrain$k1, 2, zoo::rollsum, k = k, align = "left", fill = c(NULL, NULL, NA))
    }
    
    for(i in 1:NCOL(train)){
      ts_train <- train[,i]
      idS <- colnames(train)[i]
      if(k == 1){
        listResOver[[idK]][[idS]] <- as.vector(t(sapply(1:(H/k), function(x){
          tmp_res <- ts_train-hfitted(listFit[[idK]][[idS]], h = x)
          out <- rep(NA, length(tmp_res))
          id <- x:(length(tmp_res)-((H/k)-x))
          out[1:length(id)] <- tmp_res[id]
          out
        })))
      }else{
        opt <- suppressWarnings(matrix(ts_train, nrow = k))
        opt[-c(1:length(ts_train))] <- NA
        opt_true <- opt
        
        tsinfo <- tsp(listFit[[idK]][[idS]]$x)
        res_out <- as.vector(do.call(rbind, lapply(1:NROW(opt), function(idr){
          idna <- is.na(opt[idr,])
          if(is(listFit[[idK]][[idS]], "Arima")){
            mod <- forecast::Arima(ts(opt[idr, !idna], start = tsinfo[1], frequency = tsinfo[3]), 
                         model = listFit[[idK]][[idS]])
          }else if(is(listFit[[idK]][[idS]], "ets")){
            mod <- forecast::ets(ts(opt[idr, !idna], start = tsinfo[1], frequency = tsinfo[3]), 
                       model = listFit[[idK]][[idS]], use.initial.values=TRUE)
          }
          
          true_val <- opt_true[idr, !idna]
          t(sapply(1:(H/k), function(x){
            tmp_res <- rep(NA, length(idna))
            tmp_res[!idna] <- true_val - hfitted(mod, h = x)
            out <- rep(NA, length(tmp_res))
            id <- x:(length(tmp_res)-((H/k)-x))
            out[1:length(id)] <- tmp_res[id]
            out
          }))
        })))
        
        listResOver[[idK]][[idS]] <- res_out[1:(length(ts_train)*(m/k))]
      }
      pb$tick()
    }
    listResOver[[idK]] <- simplify2array(listResOver[[idK]])
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
  
  itername <- basename(listFiles[j])
  save(listRes,
       file = file.path(".","HOverlapRes", itername))
}
rm(list = ls(all = TRUE))
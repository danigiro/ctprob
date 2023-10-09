rm(list = ls(all = TRUE))
library(FoReco)
library(forecast)
library(progress)
library(future.apply)
plan("multisession", workers = 4)

load("./VN525.RData")

args <- commandArgs(TRUE)

if(length(args)==0){
  # arima ets
  model <- "ets"
  # log lev
  trans <- "log"
  # ctjb
  boot <- "ctjb"
}else{
  model <- args[1]
  trans <- args[2]
  boot <- args[3]
}

# reconciliation
reco <- "base"

dir.create(file.path(".","ProbReco", model, trans, boot, reco), 
           recursive = TRUE, showWarnings = FALSE)

K <- c(1,2,3,4,6,12)
m <- 12
B <- 1000

listFiles <- sort(list.files(file.path(".","BaseForecasts", model, trans), 
                             full.names = TRUE))
pb <- progress_bar$new(format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = length(listFiles)*length(K), clear = FALSE, width= 60, show_after = 0)
for(j in 1:length(listFiles)){
  load(listFiles[j])
  base <- NULL
  Index_seq <- NULL
  if(boot == "ctjb"){ # Cross-temporal Joint Bootstrap
    Index <- base::sample(c(1:(nrow(na.omit(listRes[[paste0("k",m)]]))+1)), size = B , replace = TRUE)
    for(k in K){
      Mk <- m/k
      res <- listRes[[paste0("k",k)]]
      res[!is.na(res[,1]),] <- simplify2array(lapply(listFit[[paste0("k",k)]], function(x){
        x$x[x$x<=0 & !is.na(x$x)] <- min(x$x[x$x>0 & !is.na(x$x)], na.rm = TRUE)/2
        log(x$x)-log(x$fitted)
      }))
      res <- rbind(res, matrix(0, Mk, NCOL(res)))
      if(any(is.na(listRes[[paste0("k",k)]]))){
        Index_seq[[paste0("k",k)]] <- sapply(Index, function(x) seq(Mk*x+1, Mk*(x+1), by = 1))
      }else{
        Index_seq[[paste0("k",k)]] <- sapply(Index, function(x) seq(Mk*(x-1)+1, Mk*x, by = 1))
      }
      
      if(is.vector(Index_seq[[paste0("k",k)]])){
        Index_seq[[paste0("k",k)]] <- rbind(Index_seq[[paste0("k",k)]])
      }
      
      tmp <- future_lapply(1:B, 
                           function(i, res, fit){
                             id <- Index_seq[[paste0("k",k)]][,i]
                             sapply(names(fit), function(x){
                               unname(simulate(fit[[x]], innov = res[id, x], future = TRUE, 
                                               nsim = length(res[id, x])))
                             })
                           }, fit = listFit[[paste0("k",k)]], res = res,
                           future.packages = c("forecast"))
      
      base[[paste0("k",k)]] <- Reduce(rbind, tmp)
      pb$tick()
    }
  }
  
  listBase <- base
  
  itername <- basename(listFiles[j])
  save(listBase, listTest, listRes,
       Index_seq, B,
       file = file.path(".","ProbReco", model, trans, boot, reco, itername))
}
rm(list = ls(all = TRUE))

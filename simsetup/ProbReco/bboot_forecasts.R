rm(list = ls(all = TRUE))
library(FoReco)
library(forecast)
library(progress)
library(future.apply)

load("./DATA.RData")

args <- commandArgs(TRUE)

if(length(args)==0){
  # ctjb ctjbo
  boot <- "ctjb"
}else{
  boot <- args[1]
}

# reconciliation
reco <- "base"

dir.create(file.path(".","ProbReco",boot, reco), 
           recursive = TRUE, showWarnings = FALSE)

B <- 500

listFiles <- sort(list.files(file.path(".","BaseForecasts"), 
                             full.names = TRUE))
pb <- progress_bar$new(format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = length(listFiles), clear = FALSE, width= 60, show_after = 0)
for(j in 1:length(listFiles)){
  load(listFiles[j])
  base <- NULL
  Index_seq <- NULL
  if(boot == "ctjb"){ # Cross-temporal Joint Bootstrap
    Index <- base::sample(c(1:(nrow(na.omit(listRes[[paste0("k",m)]]))+1)), size = B , replace = TRUE)
    for(k in K){
      Mk <- m/k
      res <- listRes[[paste0("k",k)]]
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
                        unname(simulate(fit[[x]], innov = res[id, x], future = TRUE))
                      })
      }, fit = listFit[[paste0("k",k)]], res = res)
      
      base[[paste0("k",k)]] <- Reduce(rbind, tmp)
    }
  }else if(boot == "ctjbo"){ # Cross-temporal Joint Bootstrap
    load(file.path(".","OverlapRes", basename(listFiles[j])))
    Index <- base::sample(c(1:(nrow(na.omit(listRes[[paste0("k",m)]]))+1)), size = B , replace = TRUE)
    for(k in K){
      Mk <- m/k
      res <- listRes[[paste0("k",k)]]
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
                               unname(simulate(fit[[x]], innov = res[id, x], future = TRUE))
                             })
                           }, fit = listFit[[paste0("k",k)]], res = res)
      
      base[[paste0("k",k)]] <- Reduce(rbind, tmp)
    }
  }
  
  listBase <- base
  pb$tick()
  
  itername <- basename(listFiles[j])
  save(listBase, listTest, listRes,
       Index_seq, B,
       file = file.path(".","ProbReco", boot, reco, itername))
}
rm(list = ls(all = TRUE))

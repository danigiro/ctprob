rm(list = ls(all = TRUE))
library(FoReco)
library(forecast)
library(progress)
library(future.apply)

load("./DATA.RData")
source("./R/fupa.R")

args <- commandArgs(TRUE)

if(length(args)==0){
  # arima or ets
  model <- "arima"
  # log or lev
  trans <- "lev"
  # ctjb or csjb or indb or
  # ctjb0 or csjb0 or indb0
  boot <- "ctjb"
}else{
  # arima or ets
  model <- args[1]
  # log or lev
  trans <- args[2]
  # ctjb or csjb or tejb or indb
  boot <- args[3]
}

# reconciliation
reco <- "base"

dir.create(file.path(".","ProbReco", model, trans, boot, reco), 
           recursive = TRUE, showWarnings = FALSE)

B <- 1000

listFiles <- sort(list.files(file.path(".","BaseForecasts", model, trans), 
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
                               unname(simulate(fit[[x]], innov = res[id, x], future = TRUE, 
                                               nsim = length(res[id, x])))
                             })
                           }, fit = listFit[[paste0("k",k)]], res = res)
      
      base[[paste0("k",k)]] <- Reduce(rbind, tmp)
    }
  }else if(boot == "ctjbo"){ # Cross-temporal Joint Bootstrap
    load(file.path(".","OverlapRes", model, trans, basename(listFiles[j])))
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
                               unname(simulate(fit[[x]], innov = res[id, x], future = TRUE, 
                                               nsim = length(res[id, x])))
                             })
                           }, fit = listFit[[paste0("k",k)]], res = res)
      
      base[[paste0("k",k)]] <- Reduce(rbind, tmp)
    }
  }else if(boot == "indb"){ # Indipendent Bootstrap
    for(k in K){
      Mk <- m/k
      res <- listRes[[paste0("k",k)]]
      res <- rbind(res, matrix(0, Mk, NCOL(res)))
      Index_seq[[paste0("k",k)]] <- lapply(listFit[[paste0("k",k)]], function(x)
        matrix(sample(1:NROW(res), size = Mk*B, replace = TRUE), Mk, B))
      tmp <- future_lapply(1:B, 
                           function(i, res, fit){
                             sapply(names(fit), function(x){
                               id <- Index_seq[[paste0("k",k)]][[x]][,i]
                               unname(simulate(fit[[x]], innov = res[id, x], future = TRUE))
                             })
                           }, fit = listFit[[paste0("k",k)]], res = res)
      base[[paste0("k",k)]] <- apply(simplify2array(tmp),2,as.vector)
    }
  }else if(boot == "csjb"){ # Cross-sectional Joint Bootstrap
    for(k in K){
      Mk <- m/k
      res <- listRes[[paste0("k",k)]]
      res <- rbind(res, matrix(0, Mk, NCOL(res)))
      
      Index <- sample(1:NROW(res[-c(1:Mk), ]), size = B, replace = TRUE)
      Index_seq[[paste0("k",k)]] <- t(outer(Index, (1:Mk)-1, FUN = "+"))
      
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
      base[[paste0("k",k)]] <- apply(simplify2array(tmp),2,as.vector)
    }
  }
  
  listBase <- base
  pb$tick()
  
  itername <- basename(listFiles[j])
  save(listBase, listTest, listRes,
       Index_seq, B,
       file = file.path(".","ProbReco", model, trans, boot, reco, itername))
}
rm(list = ls(all = TRUE))

rm(list = ls(all = TRUE))
library(forecast)
library(progress)
library(Matrix)
library(FoReco)

load("./DATA.RData")
source("./R/hfitted.R")

dir.create(file.path(".","HfittedRes"), recursive = TRUE, showWarnings = FALSE)

listFiles <- sort(list.files(file.path(".","BaseForecasts"), 
                             full.names = TRUE))
pb <- progress_bar$new(format = " [:bar] :percent in :elapsed (ETA: :eta)",
                       total = length(K)*sum(dim(C))*length(listFiles), 
                       clear = FALSE, width= 60, show_after = 0)
for(j in 1:length(listFiles)){
  load(listFiles[j])
  listResH <- listRes
  for(k in K){
    idK <- paste0("k", k)
    for(i in 1:sum(dim(C))){
      idS <- c(rownames(C), colnames(C))[i]
      if(k < m){
        idna <- is.na(listResH[[idK]][, idS])
        xtrue <- listFit[[idK]][[idS]]$x

        lf <- lapply(1:(H/k), function(x){
          tmp_res <- xtrue-hfitted(listFit[[idK]][[idS]], h = x)
          tmp_res[rep(1:(m/k), length.out = length(tmp_res)) == x]
        })
        
        res_mat <- sapply(1:max(sapply(lf, length)), function(x) sapply(lf, function(l) ifelse(x<=length(l),l[x], NA)))
        listResH[[idK]][!idna, idS] <- as.vector(res_mat)[1:sum(!idna)]
      }
      pb$tick()
    }
  }
  allres <- t(do.call(rbind, listResH[paste0("k", sort(K, decreasing = TRUE))]))
  
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
       file = file.path(".","HfittedRes", itername))
}
rm(list = ls(all = TRUE))

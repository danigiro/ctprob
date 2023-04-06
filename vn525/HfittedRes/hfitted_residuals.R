rm(list = ls(all = TRUE))
library(forecast)
library(progress)
library(Matrix)
library(FoReco)

load("./VN525.RData")

# arima ets
model <- "ets"
# log lev
trans <- "log"

dir.create(file.path(".","HfittedRes",model, trans), recursive = TRUE, showWarnings = FALSE)

listFiles <- sort(list.files(file.path(".","BaseForecasts", model, trans), 
                             full.names = TRUE))
K <- c(1,2,3,4,6,12)
m <- 12
H <- 12
fixed_length <- 132 # 1998.01 - 2008.12 (first)
end_traing <- NROW(VNdata)-H-fixed_length+1

#Test
#end_traing = 2
#VNdata = VNdata[,1:3]
start_year <- min(time(VNdata))
for(j in 0:(end_traing-1)){
  load(listFiles[j+1])
  listResH <- listRes
  
  pb <- progress_bar$new(format = paste0("Replication n. ", j+1, 
                                         " [:bar] :percent in :elapsed (ETA: :eta)"),
                         total = length(K), clear = FALSE, width= 60, show_after = 0)
  for(k in K){
    idK <- paste0("k", k)
    for(i in 1:NCOL(VNdata)){
      idS <- colnames(VNdata)[i]
      if(k < m){
        idna <- is.na(listResH[[idK]][, idS])
        xtrue <- listFit[[idK]][[idS]]$x
        if(trans == "log"){
          xlog <- xtrue
          xlog[xlog<=0] <- min(xlog[xlog>0])/2
          listFit[[idK]][[idS]]$x <- xlog
        }
        lf <- lapply(1:(H/k), function(x){
          tmp_res <- xtrue-fitted(listFit[[idK]][[idS]], h = x)
          tmp_res[rep(1:(m/k), length.out = length(tmp_res)) == x]
        })
        
        res_mat <- sapply(1:max(sapply(lf, length)), function(x) sapply(lf, function(l) ifelse(x<=length(l),l[x], NA)))
        listResH[[idK]][!idna, idS] <- as.vector(res_mat)[1:sum(!idna)]
      }
    }
    pb$tick()
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
  
  itername <- paste0("ite_", formatC(j+1, width = nchar(end_traing+1), 
                                     format = "d", flag = "0"), ".RData")
  save(listRes,
       file = file.path(".","HfittedRes",model, trans, itername))
}
rm(list = ls(all = TRUE))

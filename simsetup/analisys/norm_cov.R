library(progressr)
library(tidyverse)
rho <- -0.8
Sigma <- diag(c(0.9, 1.8))%*%matrix(c(1, rho,
                                      rho, 1), 2, 2)%*%diag(c(0.9, 1.8))

PhiB <- c(1.3446416, -0.7392638)
PhiC <- c(0.9454175, -0.4218556)

hfbts_cov <- function(Sigma, PhiB, PhiC){
  s11 <- Sigma[1,1]
  s21 <- PhiB[1]*Sigma[1,1]
  s31 <- Sigma[2,1]
  s41 <- PhiC[1]*Sigma[2,1]
  s22 <- Sigma[1,1]*(1+PhiB[1]^2)
  s32 <- PhiB[1]*Sigma[2,1]
  s42 <- Sigma[2,1]*(1+PhiC[1]*PhiB[1])
  s33 <- Sigma[2,2]
  s34 <- PhiC[1]*Sigma[2,2]
  s44 <- Sigma[2,2]*(1+PhiC[1]^2)
  matrix(c(s11,s21,s31,s41,
           s21,s22,s32,s42,
           s31,s32,s33,s34,
           s41,s42,s34,s44), 4, byrow = TRUE)
}
S <- FoReco::ctf_tools(C = t(c(1,1)), m = 2)$ctf$Fmat
true_cov <- S%*%hfbts_cov(Sigma, PhiB, PhiC)%*%t(S)


### ----
extract_norm <- function(x, true_cov){
  loadRData <- load_list(x)
  if(any(names(loadRData) == "listFree")){
    norm_r <- lapply(names(loadRData$listFree), function(comb_x){
      cov <- extract_cov(loadRData$listFree[[comb_x]], K = K)
      tibble(comb = comb_x,
             ite = readr::parse_number(basename(x)),
             Frobenius = norm(cov - true_cov, "F"),
             One = norm(cov - true_cov, "O"),
             Maximum = norm(cov - true_cov, "M"))
    })
    Reduce(rbind, norm_r)
  }else{
    cov <- extract_cov(loadRData$listBase, K = K)
    tibble(comb = NA,
           ite = readr::parse_number(basename(x)),
           Frobenius = norm(cov - true_cov, "F"),
           One = norm(cov - true_cov, "O"),
           Maximum = norm(cov - true_cov, "M"))
  }
}

extract_cov <- function(listMat, K){
  mat <- t(do.call(rbind, listMat[paste0("k", sort(K, decreasing = TRUE))]))
  N <- NCOL(mat) / FoReco::thf_tools(m = K)$kt
  DN <- FoReco:::Dmat(h = N, m = K, n = NROW(mat))
  math <- matrix(DN %*% as.vector(t(mat)), nrow = N, byrow = TRUE)
  cov(math)
}

load_list <- function(x){
  load(x,  temp_env <- new.env())
  as.list(temp_env)
}
### ----
K <- c(2,1)
prob <- list.dirs("ProbReco", recursive = FALSE)
data <- NULL
for(i in 1:length(prob)){
  meth <- list.dirs(prob[i], recursive = FALSE)
  pb <- progress_bar$new(format = paste0("(", basename(prob[i]), ")",
                                         " [:bar] :percent in :elapsed (ETA: :eta)"),
                         total = length(meth), clear = FALSE, width= 60, show_after = 0)
  for(j in 1:length(meth)){
    complete <- lapply(sort(list.files(meth[j], full.names = TRUE)), 
                       extract_norm, true_cov = true_cov)
    data <- rbind(data, 
                  Reduce(rbind, complete) |>
                    add_column(prob = basename(prob[i]),
                               meth = basename(meth[j]), .before = 1))
    pb$tick()
  }
}
save(data, true_cov, rho, Sigma, PhiB, PhiC, file = "./analisys/norm_cov.RData")
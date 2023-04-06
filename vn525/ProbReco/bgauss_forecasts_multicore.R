rm(list = ls(all = TRUE))
library(future.apply)
library(progressr)
plan("multisession", workers = 5)

handlers("progress")

args <- commandArgs(TRUE)

if(length(args)==0){
  # arima ets
  model <- "ets"
  # log lev
  trans <- "log"
  # ctsam hbsam hsam bsam ctshr hbshr hshr bshr
  boot <- "ctsam"
  # in h
  restype <- "in"
}else{
  model <- args[1]
  trans <- args[2]
  gauss <- args[3]
  if(length(args) < 4){
    restype <- "in"
  }else{
    restype <- args[4]
  }
}

# reconciliation
reco <- "base"
if(restype == "o"){
  loadRes <- "OverlapRes"
  gaussn <- paste0(gauss, restype)
}else if(restype == "h"){
  loadRes <- "HfittedRes"
  gaussn <- paste0(gauss, restype)
}else if(restype == "oh"){
  loadRes <- "HOverlapRes"
  gaussn <- paste0(gauss, restype)
}else{
  loadRes <- NULL
  gaussn <- gauss
}

dir.create(file.path(".","ProbReco", model, trans, gaussn, reco), 
           recursive = TRUE, showWarnings = FALSE)

listFiles <- sort(list.files(file.path(".","BaseForecasts", model, trans), 
                             full.names = TRUE))

gauss_fun <- function(j, model, trans, gauss, gaussn, loadRes){
  require(FoReco)
  require(Matrix)
  require(forecast)
  load("./VN525.RData")
  K <- c(1,2,3,4,6,12)
  m <- 12
  B <- 1000
  # reconciliation
  reco <- "base"
  
  # fast shr
  shrink_estim_fast <- function(x, minT = T) {
    if (is.matrix(x) == TRUE && is.numeric(x) == FALSE) {
      stop("The data matrix must be numeric!", call. = FALSE)
    }
    
    x <- stats::na.omit(x)
    p1 <- ncol(x)
    n2 <- nrow(x)
    
    if (minT == T) {
      covm <- crossprod(x) / n2
    } else {
      covm <- stats::cov(x)
    }
    
    tar <- diag(diag(covm))
    corm <- correlateR::cov2cor(covm)
    xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
    xs <- xs[stats::complete.cases(xs), ]
    v <- (1 / (n2 * (n2 - 1))) * (crossprod(xs^2) - 1 / n2 * (crossprod(xs))^2)
    diag(v) <- 0
    corapn <- diag(p1)
    d <- (corm - corapn)^2
    lambda <- sum(v) / sum(d)
    lambda <- max(min(lambda, 1), 0)
    shrink.cov <- lambda * tar + (1 - lambda) * covm
    return(list(cov = shrink.cov, lambda = lambda))
  }
  
  load(listFiles[j])
  if(!is.null(loadRes)){
    load(file.path(".",loadRes, model, trans, basename(listFiles[j])))
  }
  if(gauss == "ctshr"){
    base <- t(do.call(rbind, listBase[paste0("k", sort(K, decreasing = TRUE))]))
    res <- t(do.call(rbind, listRes[paste0("k", sort(K, decreasing = TRUE))]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    cov <- shrink_estim_fast(E)
    lambda <- cov$lambda
    cov <- cov$cov
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }else if(gauss == "ctsam"){
    base <- t(do.call(rbind, listBase[paste0("k", sort(K, decreasing = TRUE))]))
    res <- t(do.call(rbind, listRes[paste0("k", sort(K, decreasing = TRUE))]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    
    cov <- crossprod(na.omit(E))/NROW(na.omit(E))
    lambda <- 0
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }else if(gauss == "hbsam"){
    Ksort <- sort(K, decreasing = TRUE)
    base <- t(do.call(rbind, listBase[paste0("k", Ksort)]))
    res <- t(do.call(rbind, listRes[paste0("k", Ksort)]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    
    kid <- rep(rep(Ksort, m/Ksort), sum(dim(C)))
    nid <- rep(1:sum(dim(C)), each = sum(Ksort))
    E <- E[, nid > sum(NROW(C)) & kid == 1]
    cov1 <- crossprod(na.omit(E))/NROW(na.omit(E))
    Sct <- FoReco::ctf_tools(C = C, m = m)$ctf$Fmat
    cov <- Sct%*%cov1%*%t(Sct)
    lambda <- 0
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }else if(gauss == "hbshr"){
    Ksort <- sort(K, decreasing = TRUE)
    base <- t(do.call(rbind, listBase[paste0("k", Ksort)]))
    res <- t(do.call(rbind, listRes[paste0("k", Ksort)]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    
    kid <- rep(rep(Ksort, m/Ksort), sum(dim(C)))
    nid <- rep(1:sum(dim(C)), each = sum(Ksort))
    E <- E[, nid > sum(NROW(C)) & kid == 1]
    cov <- shrink_estim_fast(E)
    lambda <- cov$lambda
    cov1 <- cov$cov
    
    Sct <- FoReco::ctf_tools(C = C, m = m)$ctf$Fmat
    cov <- Sct%*%cov1%*%t(Sct)
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }else if(gauss == "hsam"){
    Ksort <- sort(K, decreasing = TRUE)
    base <- t(do.call(rbind, listBase[paste0("k", Ksort)]))
    res <- t(do.call(rbind, listRes[paste0("k", Ksort)]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    
    kid <- rep(rep(Ksort, m/Ksort), sum(dim(C)))
    #nid <- rep(1:sum(dim(C)), each = sum(Ksort))
    E <- E[, kid == 1]
    cov1 <- crossprod(na.omit(E))/NROW(na.omit(E))
    kSte <- kronecker(Diagonal(sum(dim(C))), FoReco::ctf_tools(C = C, m = m)$thf$R)
    cov <- kSte%*%cov1%*%t(kSte)
    lambda <- 0
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }else if(gauss == "hshr"){
    Ksort <- sort(K, decreasing = TRUE)
    base <- t(do.call(rbind, listBase[paste0("k", Ksort)]))
    res <- t(do.call(rbind, listRes[paste0("k", Ksort)]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    
    kid <- rep(rep(Ksort, m/Ksort), sum(dim(C)))
    #nid <- rep(1:sum(dim(C)), each = sum(Ksort))
    E <- E[, kid == 1]
    cov <- shrink_estim_fast(E)
    lambda <- cov$lambda
    cov1 <- cov$cov
    kSte <- kronecker(Diagonal(sum(dim(C))), FoReco::ctf_tools(C = C, m = m)$thf$R)
    cov <- kSte%*%cov1%*%t(kSte)
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }else if(gauss == "bsam"){
    Ksort <- sort(K, decreasing = TRUE)
    base <- t(do.call(rbind, listBase[paste0("k", Ksort)]))
    res <- t(do.call(rbind, listRes[paste0("k", Ksort)]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    
    #kid <- rep(rep(Ksort, m/Ksort), sum(dim(C)))
    nid <- rep(1:sum(dim(C)), each = sum(Ksort))
    E <- E[, nid >NROW(C)]
    cov1 <- crossprod(na.omit(E))/NROW(na.omit(E))
    kScs <- kronecker(FoReco::ctf_tools(C = C, m = m)$hts$S, Diagonal(sum(m/Ksort)))
    cov <- kScs%*%cov1%*%t(kScs)
    lambda <- 0
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }else if(gauss == "bshr"){
    Ksort <- sort(K, decreasing = TRUE)
    base <- t(do.call(rbind, listBase[paste0("k", Ksort)]))
    res <- t(do.call(rbind, listRes[paste0("k", Ksort)]))
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    
    #kid <- rep(rep(Ksort, m/Ksort), sum(dim(C)))
    nid <- rep(1:sum(dim(C)), each = sum(Ksort))
    E <- E[, nid >NROW(C)]
    cov <- shrink_estim_fast(E)
    lambda <- cov$lambda
    cov1 <- cov$cov
    kScs <- kronecker(FoReco::ctf_tools(C = C, m = m)$hts$S, Diagonal(sum(m/Ksort)))
    cov <- kScs%*%cov1%*%t(kScs)
    
    tmp <- tryCatch({
      Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov)
    },
    error = function(cond){
      tryCatch({
        Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.01, NROW(cov)))
      },
      error = function(cond){
        tryCatch({
          Rfast::rmvnorm(B, mu = as.vector(t(base)), sigma = cov + diag(0.05, NROW(cov)))
        },
        error = function(cond){
          MASS::mvrnorm(B, mu = as.vector(t(base)), Sigma = cov)
        })
      })
    })
    
    DB <- FoReco:::Dmat(h = B, m = K, n = NROW(res))
    tmp <- matrix(t(DB) %*% as.vector(t(tmp)), nrow = NROW(res), byrow = TRUE)
    
    listBase <- lapply(split(t(tmp), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                         B*m/sort(K, decreasing = TRUE))), 
                       function(x){
                         out <- matrix(x, ncol = sum(dim(C)))
                         colnames(out) <- c(rownames(C), colnames(C))
                         out
                       })[paste0("k", sort(K))]
  }
  p()
  
  itername <- basename(listFiles[j])
  save(listBase, listTest, listRes, B, lambda,
       file = file.path(".","ProbReco", model, trans, gaussn, reco, itername))
}

xs <- 1:length(listFiles)
with_progress({
  p <- progressor(along = xs)
  y <- future_lapply(xs, gauss_fun, gauss = gauss, gaussn = gaussn, loadRes = loadRes,
                     model = model, trans = trans, future.seed=TRUE)
})

rm(list = ls(all = TRUE))
rm(list = ls(all = TRUE))
suppressPackageStartupMessages(library(FoReco))
suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(tidyverse))

load("./VN525.RData")
source("./R/pscore_fun.r")
args <- commandArgs(TRUE)
if(length(args)==0){
  # arima ets
  model <- "ets"
  # log lev
  trans <- "log"
  # ctjb ctsam hbsam hsam bsam ctshr hbshr hshr bshr
  boot <- "ctsam"
  # in h
  restype <- "in"
  # free sntz
  basen <- "free"
}else{
  model <- args[1]
  trans <- args[2]
  boot <- args[3]
  if(length(args) < 4){
    restype <- "in"
    basen <- "free"
  }else{
    restype <- args[4]
    if(length(args) < 5){
      basen <- "free"
    }else{
      basen <- args[5]
    }
  }
}

if(basen == "free"){
  bootB <- boot
}else{
  bootB <- boot
  boot <- paste0(boot, "0")
}
# reconciliation
reco <- "octf"
if(restype == "o"){
  loadRes <- "OverlapRes"
  reco <- paste0(reco, restype)
  listComb <- c("wlsh", "wlsv", "bdshr")
}else if(restype == "h"){
  loadRes <- "HfittedRes"
  reco <- paste0(reco, restype)
  listComb <- c("wlsh", "wlsv")
}else if(restype == "oh"){
  loadRes <- "HOverlapRes"
  reco <- paste0(reco, restype)
  listComb <- c("wlsh", "wlsv")
}else{
  loadRes <- NULL
  listComb <- c("ols", "struc", "wlsh", "wlsv", "bdshr")
}

dir.create(file.path(".","ProbReco", model, trans, boot, reco), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(".","ProbScore", model, trans), 
           recursive = TRUE, showWarnings = FALSE)

K <- c(1,2,3,4,6,12)
m <- 12

listFiles <- sort(list.files(file.path(".","ProbReco", model, 
                                       trans, bootB, "base"), 
                             full.names = TRUE))
df_crps <- df_es <- df_q <- NULL
probs_q <- sort(c(0, 0.005, 0.0125, 0.025, 0.05, 0.1, 0.125, 0.25, 1, 0.995, 
                  0.9875, 0.975, 0.95, 0.9, 0.875, 0.75, 0.5))

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

pb <- progress_bar$new(format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = length(listFiles), clear = FALSE, 
                       width= 60, show_after = 0)
for(j in 1:length(listFiles)){
  load(listFiles[j])
  if(!is.null(loadRes)){
    load(file.path(".",loadRes, model, trans, basename(listFiles[j])))
  }
  listFree <-listSntz <- NULL
  
  base <- t(do.call(rbind, listBase[paste0("k", sort(K, decreasing = TRUE))]))
  if(basen != "free"){
    base[base<0] <- 0
  }
  res <- t(do.call(rbind, listRes[paste0("k", sort(K, decreasing = TRUE))]))
  
  for(comb in listComb){
    free <- NULL
    
    time_free_start <- Sys.time()
    if(comb == "hshr"){
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
      
      free <- octrec(basef = base, C = C, m = m, comb = "omega", Omega = cov, 
                     keep = "recf", type = "M")
    }else if(comb == "bshr"){
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
      
      free <- octrec(basef = base, C = C, m = m, comb = "omega", Omega = cov, 
                     keep = "recf", type = "M")
    }else if(comb == "hbshr"){
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
      cov1 <- cov$cov
      
      Sct <- FoReco::ctf_tools(C = C, m = m)$ctf$Fmat
      cov <- Sct%*%cov1%*%t(Sct)
      free <- octrec(basef = base, C = C, m = m, comb = "omega", Omega = cov, 
                     keep = "recf", type = "M")
    }else if(comb == "shr"){
      res <- t(do.call(rbind, listRes[paste0("k", sort(K, decreasing = TRUE))]))
      N <- NCOL(res) / sum(K)
      DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
      E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
      cov <- shrink_estim_fast(E)
      cov <- cov$cov
      free <- octrec(basef = base, C = C, m = m, comb = "omega", Omega = cov, 
                     keep = "recf", type = "M")
    }else{
      free <- octrec(basef = base, C = C, m = m, comb = comb, res = res, 
                     keep = "recf", type = ifelse(comb %in% c("bdshr"), "S", "M"))
    }
    time_free_end <- Sys.time()
    time_free_end-time_free_start
    if(any(free < 0)){
      free0 <- free
      free0[free0<0] <- 0
      ks <- sum(m/sort(K, decreasing = TRUE)) - m
      sntz <- ctbu(free0[-c(1:NROW(C)), -c(1:(ks*B))], C = C, m = m)
      dimnames(sntz) <- dimnames(free)
    }else{
      sntz <- free
    }
    
    listFree[[comb]] <- lapply(split(t(free), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                                  B*m/sort(K, decreasing = TRUE))), 
                               function(x){
                                 out <- matrix(x, ncol = 525)
                                 colnames(out) <- colnames(VNdata)
                                 out
                               })[paste0("k", sort(K))]
    
    listSntz[[comb]] <- lapply(split(t(sntz), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                                  B*m/sort(K, decreasing = TRUE))), 
                               function(x){
                                 out <- matrix(x, ncol = 525)
                                 colnames(out) <- colnames(VNdata)
                                 out
                               })[paste0("k", sort(K))]
  }
  
  df_test <- list2tibble_1Lv(listTest) |> 
    rename(test = value)
  df_free <- list2tibble_2Lm(listFree, B = B) |>
    add_column(type = reco,
               nn = "free")
  df_sntz <- list2tibble_2Lm(listSntz, B = B) |>
    add_column(type = reco,
               nn = "sntz")
  df_join <- left_join(rbind(df_free, df_sntz), df_test, by = c("serie", "k"))
  
  df_crps_tmp <- df_join |>
    group_by(serie, k, comb, type, nn) |>
    summarise(crps = crps(test, value),
              .groups = "drop") |>
    unnest_longer(crps, indices_to = "h") |> 
    pivot_wider(names_from = c(type, comb, nn), values_from = crps, names_sep = "-") |>
    mutate(k = readr::parse_number(k)) |>
    add_column(prob = boot,
               ite = j, .before = 1)
  df_crps <- rbind(df_crps, df_crps_tmp)
  
  df_es_tmp <- df_join |>
    group_by(k, comb, type, nn) |>
    summarise(value = arrange_value(value, serie),
              test = arrange_value(test, serie),
              es = es(test, value),
              .groups = "drop") |>
    select(-c(value, test)) |>
    unnest_longer(es, indices_to = "h") |> 
    pivot_wider(names_from = c(type, comb, nn), values_from = es, names_sep = "-") |>
    mutate(k = readr::parse_number(k)) |>
    add_column(prob = boot,
               ite = j, 
               serie = "all", .before = 1)
  
  df_es2_tmp <- df_join |>
    mutate(nser = serie,
           serie = ifelse(serie %in% colnames(C), "bts", "uts")) |>
    group_by(k, comb, type, nn, serie) |>
    summarise(value = arrange_value(value, nser),
              test = arrange_value(test, nser),
              es = es(test, value),
              .groups = "drop") |>
    select(-c(value, test)) |>
    unnest_longer(es, indices_to = "h") |> 
    pivot_wider(names_from = c(type, comb, nn), values_from = es, names_sep = "-") |>
    mutate(k = readr::parse_number(k)) |>
    add_column(prob = boot,
               ite = j, .before = 1)
  df_es <- rbind(df_es, df_es_tmp, df_es2_tmp)
  
  if(j %in% c(1, length(listFiles))){
    df_q_test <- df_join |>
      select(-value) |>
      unnest_longer(test, indices_to = "h")
    
    df_q_tmp <- df_join |>
      select(-test) |>
      mutate(value =  lapply(value, asplit, MARGIN = 2)) |>
      unnest_longer(value, indices_to = "h") |>
      mutate(mean = sapply(value, function(x) mean(x)),
             var = sapply(value, function(x) max(var(x), 0.1)),
             value = lapply(value, function(x) quantile(x, probs = probs_q))) |>
      #value = lapply(value, function(x) quantile(x, probs = seq(0, 1, 0.005)))) |>
      unnest_longer(value, indices_to = "q") |>
      pivot_wider(names_from = q)
    
    df_q_tmp <- left_join(df_q_tmp, df_q_test, by = c("h", "serie", "k", "comb", "type", "nn")) |>
      add_column(prob = boot,
                 ite = j, .before = 1)|>
      mutate(k = readr::parse_number(k))
    
    df_q <- rbind(df_q, df_q_tmp)
  }
  
  
  if(j == 0){
    itername <- basename(listFiles[j])
    save(listSntz, listFree, #listOsqp, info_osqp, 
         file = file.path(".","ProbReco", model, trans, boot, reco, itername))
  }
  pb$tick()
}
save(df_q, df_es, df_crps, 
     file = file.path(".", "ProbScore", model, trans, 
                      paste0(paste("rboot", boot, reco, sep = "_"), 
                             ".RData")))
rm(list = ls(all = TRUE))
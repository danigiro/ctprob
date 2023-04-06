rm(list = ls(all = TRUE))
suppressPackageStartupMessages(library(FoReco))
suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(tidyverse))

load("./VN525.RData")
source("./R/pscore_fun.r")
source("./R/sntz.r")
args <- commandArgs(TRUE)

if(length(args)==0){
  # arima ets
  model <- "ets"
  # log lev
  trans <- "log"
  # ctjb ctsam hbsam hsam bsam ctshr hbshr hshr bshr
  boot <- "ctsam"
  # free sntz
  basen <- "free"
}else{
  model <- args[1]
  trans <- args[2]
  boot <- args[3]
  if(length(args) < 4){
    basen <- "free"
  }else{
    basen <- args[4]
  }
}

if(basen == "free"){
  bootB <- boot
}else{
  bootB <- boot
  boot <- paste0(boot, "0")
}
# reconciliation
reco <- "csbu"

dir.create(file.path(".","ProbReco", model, trans, boot, reco), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(".","ProbScore", model, trans), recursive = TRUE, showWarnings = FALSE)

K <- c(1,2,3,4,6,12)
m <- 12
listComb <- c("bu", "ols", "struc", "wls", "shr")

listFiles <- sort(list.files(file.path(".","ProbReco", model, trans, bootB, "base"), 
                             full.names = TRUE))
df_crps <- df_es <- df_q <- NULL
probs_q <- sort(c(0, 0.005, 0.0125, 0.025, 0.05, 0.1, 0.125, 0.25, 1, 0.995, 
                  0.9875, 0.975, 0.95, 0.9, 0.875, 0.75, 0.5))

pb <- progress_bar$new(format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = length(listFiles), clear = FALSE, width= 60, show_after = 0)
for(j in 1:length(listFiles)){
  load(listFiles[j])
  listFree <-listSntz <- NULL
  
  base <- listBase$k1
  if(basen != "free"){
    base[base<0] <- 0
  }
  
  for(comb in listComb){
    free <- NULL
    bts <- htsrec(basef = base, C = C, comb = comb, res = listRes$k1, keep = "recf")
    free <- ctbu(t(bts[,-c(1:NROW(C))]), C = C, m = m)
    
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
    save(listSntz, listFree, 
         file = file.path(".","ProbReco", model, trans, boot, reco, itername))
  }
  pb$tick()
}
save(df_q, df_es, df_crps, 
     file = file.path(".", "ProbScore", model, trans, 
                      paste0(paste("rboot", boot, reco, sep = "_"), 
                             ".RData")))
rm(list = ls(all = TRUE))
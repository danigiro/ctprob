rm(list = ls(all = TRUE))
library(FoReco)
library(forecast)
library(progress)
library(tidyverse)

load("./DATA.RData")
load("./DataRaw/groups_name.RData")
source("./R/pscore_fun.r")
args <- commandArgs(TRUE)

if(length(args)==0){
  # arima ets
  model <- "arima"
  # log lev
  trans <- "lev"
  # ctjb ctsam hsam ctshr hshr
  boot <- "ctsam"
  # res type
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
  boot <- paste(boot, basen)
}
# reconciliation
reco <- "oct"
if(restype == "o"){
  loadRes <- "OverlapRes"
  reco <- paste0(reco, restype)
  listComb <- c("wlsh", "wlsv", "bdshr", "shr", 
                "hbshr", "hshr", "bshr")
}else if(restype == "h"){
  loadRes <- "HfittedRes"
  reco <- paste0(reco, restype)
  listComb <- c("wlsh", "wlsv", "shr", 
                "hbshr", "hshr", "bshr")
}else if(restype == "oh"){
  loadRes <- "HOverlapRes"
  reco <- paste0(reco, restype)
  listComb <- c("wlsh", "wlsv", "shr", 
                "hbshr", "hshr", "bshr")
}else{
  loadRes <- NULL
  listComb <- c("ols", "struc", "wlsh", "wlsv", "bdshr", "shr", 
                "hbshr", "hshr", "bshr")
}

dir.create(file.path(".","ProbReco", model, trans, boot, reco), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(".","ProbScore", model, trans), recursive = TRUE, showWarnings = FALSE)

listFiles <- sort(list.files(file.path(".","ProbReco", model, trans, bootB, "base"), 
                             full.names = TRUE))
df_crps <- df_es <- df_q <- NULL
probs_q <- sort(c(0, 0.005, 0.0125, 0.025, 0.05, 0.1, 0.125, 0.25, 1, 0.995, 
                  0.9875, 0.975, 0.95, 0.9, 0.875, 0.75, 0.5))

pb <- progress_bar$new(format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = length(listFiles), clear = FALSE, width= 60, show_after = 0)
for(j in 1:length(listFiles)){
  load(listFiles[j])
  if(!is.null(loadRes)){
    load(file.path(".",loadRes, model, trans, basename(listFiles[j])))
  }
  listFree <- NULL
  
  base <- t(do.call(rbind, listBase[paste0("k", sort(K, decreasing = TRUE))]))
  if(basen != "free"){
    base[base<0] <- 0
  }
  res <- t(do.call(rbind, listRes[paste0("k", sort(K, decreasing = TRUE))]))
  
  for(comb in listComb){
    free <- octrec(basef = base, C = C, m = m, comb = comb, res = res, 
                   keep = "recf", type = ifelse(comb %in% c("bdshr", "hbshr", "hshr", "bshr"), 
                                                "S", "M"))
    
    listFree[[comb]] <- lapply(split(t(free), rep(paste0("k", sort(K, decreasing = TRUE)), 
                                                  B*m/sort(K, decreasing = TRUE))), 
                               function(x){
                                 out <- matrix(x, ncol = sum(dim(C)))
                                 colnames(out) <- colnames(DATA)
                                 out
                               })[paste0("k", sort(K))]
  }
  
  df_test <- list2tibble_1Lv(listTest) |> 
    rename(test = value)
  df_free <- list2tibble_2Lm(listFree, B = B) |>
    add_column(type = reco,
               nn = "free")
  df_join <- left_join(df_free, df_test, by = c("serie", "k"))
  
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
    filter(serie %in% names$Inc) |>
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
               serie = "inc",
               ite = j, .before = 1)
  
  df_es3_tmp <- df_join |> 
    filter(serie %in% names$Exp) |>
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
               serie = "exp",
               ite = j, .before = 1)
  df_es <- rbind(df_es, df_es_tmp, df_es2_tmp, df_es3_tmp)
  
  df_q_test <- df_join |>
    select(-value) |>
    unnest_longer(test, indices_to = "h")
  
  df_q_tmp <- df_join |>
    select(-test) |>
    mutate(value =  lapply(value, asplit, MARGIN = 2)) |>
    unnest_longer(value, indices_to = "h") |>
    mutate(mean = sapply(value, function(x) mean(x)),
           var = sapply(value, function(x) var(x)),
           value = lapply(value, function(x) quantile(x, probs = probs_q))) |>
    unnest_longer(value, indices_to = "q") |>
    pivot_wider(names_from = q)
  
  df_q_tmp <- left_join(df_q_tmp, df_q_test, by = c("h", "serie", "k", "comb", "type", "nn")) |>
    add_column(ite = j) |>
    mutate(k = readr::parse_number(k))
  
  df_q <- rbind(df_q, df_q_tmp)
  
  if(j == 1){
    itername <- basename(listFiles[j])
    save(listFree, 
         file = file.path(".","ProbReco", model, trans, boot, reco, itername))
  }
  pb$tick()
}
save(df_q, df_es, df_crps, 
     file = file.path(".", "ProbScore", model, trans, 
                      paste0(paste("rboot", boot, reco, sep = "_"), 
                             ".RData")))
rm(list = ls(all = TRUE))

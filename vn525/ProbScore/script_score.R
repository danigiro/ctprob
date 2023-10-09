library(tidyverse)
library(progress)

if(length(args)==0){
  # arima ets
  model <- "ets"
  # log lev
  trans <- "log"
}else{
  model <- args[1]
  trans <- args[2]
}

scores_file <- list.files(file.path(".","ProbScore", model, trans), full.names = TRUE)
df_es_mean <- NULL
df_crps_mean <- NULL
df_mis_mean <- NULL
pb <- progress_bar$new(format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = length(scores_file), clear = FALSE, width= 60, show_after = 0)

for(i in 1:length(scores_file)){
  load(scores_file[i])
  df_mis_mean <- rbind(df_mis_mean,
                       df_q |>
                         rename(lower = "5%", 
                                upper = "95%")|>
                         add_column(alpha = 0.1) |>
                         select(h, serie, k, comb, type, nn, alpha, test, ite, lower, upper) |>
                         mutate(ind_lower = test < lower,
                                ind_upper = test > upper,
                                value = (upper -lower) + ((2*(lower - test)*ind_lower)/alpha) + ((2*(test-upper)*ind_upper)/alpha)) |>
                         group_by(h, serie, k, comb, type, nn, alpha) |>
                         summarise(value = mean(value), .groups = "drop") |>
                         add_column(prob = unique(df_es$prob)),
                       df_q |>
                         rename(lower = "5%", 
                                upper = "95%")|>
                         add_column(alpha = 0.1) |>
                         select(h, serie, k, comb, type, nn, alpha, test, ite, lower, upper) |>
                         mutate(ind_lower = test < lower,
                                ind_upper = test > upper,
                                value = (upper -lower) + ((2*(lower - test)*ind_lower)/alpha) + ((2*(test-upper)*ind_upper)/alpha)) |>
                         group_by(serie, k, comb, type, nn, alpha) |>
                         summarise(value = mean(value), .groups = "drop") |>
                         add_column(prob = unique(df_es$prob), h = 0))
  rm(df_q)
  df_es_mean <- rbind(df_es_mean, 
                      df_es |>
                        pivot_longer(-c(prob, serie, ite, k, h)) |>
                        group_by(prob, serie, k, name, h) |>
                        summarise(value = mean(value), .groups = "drop") |>
                        separate(name, into = c("meth", "comb", "nn")),
                      df_es |>
                        pivot_longer(-c(prob, serie, ite, k, h)) |>
                        group_by(prob, serie, k, name) |>
                        summarise(value = mean(value), .groups = "drop") |>
                        separate(name, into = c("meth", "comb", "nn")) |>
                        mutate(h = 0))
  
  df_crps_mean <- rbind(df_crps_mean, 
                        df_crps |>
                          pivot_longer(-c(prob, serie, ite, k, h)) |>
                          group_by(prob, serie, k, name, h) |>
                          summarise(value = mean(value), .groups = "drop") |>
                          separate(name, into = c("meth", "comb", "nn")),
                        df_crps |>
                          pivot_longer(-c(prob, serie, ite, k, h)) |>
                          group_by(prob, serie, k, name) |>
                          summarise(value = mean(value), .groups = "drop") |>
                          separate(name, into = c("meth", "comb", "nn")) |>
                          mutate(h = 0))
  pb$tick()
}
save(df_crps_mean, df_es_mean, df_mis_mean, 
     file = file.path(".","ProbScore", paste(model, trans, "scores.RData", sep = "_")))

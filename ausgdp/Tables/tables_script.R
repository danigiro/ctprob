# Tables AusGDP ----
library(formattable)
library(kableExtra)
library(tidyverse)
load("./ProbScore/arima_lev_scores.RData")

## Paper ----
sel_mc <- c("arima lev" = "base",
            "csbu shr" = "ct$(shr_{cs}, bu_{te})$",
            "csbu wls" = "ct$(wls_{cs}, bu_{te})$",
            "octo wlsv" = "oct$_o(wlsv)$",
            "octo bdshr" = "oct$_o(bdshr)$",
            "octo shr" = "oct$_o(shr)$",
            "octo hshr" = "oct$_o(hshr)$",
            "octoh shr" = "oct$_{oh}(shr)$",
            "octoh hshr" = "oct$_{oh}(hshr)$")
sel_prob <- c("ctjb", "ctsamh", "hsamh", "ctsamoh", "hsamoh")
K <- setNames(c("$\\forall k \\\\in \\\\{4,2,1\\\\}$", 
                paste0("$k = ", c(1,2,4), "$")), c(0,1,2,4))

### CRPS ----
crps_k <- left_join(df_crps_mean, 
                    df_crps_mean |>
                      filter(meth == "arima", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "nn", "prob")) |>
                      rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  group_by(prob, k, meth, comb, h) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  filter(h == 0) |> select(-h) 
crps_0 <- left_join(df_crps_mean, 
                    df_crps_mean |>
                      filter(meth == "arima", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "nn", "prob")) |>
                      rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  group_by(prob, meth, comb, h) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  filter(h == 0) |> select(-h) |> add_column(k = 0, .before = 1)

rbind(crps_0, crps_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |> select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: sample covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/sam_crps.tex"), self_contained = FALSE)

### ES ----
es_k <- left_join(df_es_mean, 
                  df_es_mean |>
                    filter(meth == "arima", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "nn", "prob")) |>
                    rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  filter(h == 0, serie == "all") |>
  group_by(prob, k, meth, comb) |>
  summarise(value = exp(mean(log(value))), .groups = "drop")
es_0 <- left_join(df_es_mean, 
                  df_es_mean |>
                    filter(meth == "arima", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "nn", "prob")) |>
                    rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  filter(h == 0, serie == "all") |>
  group_by(prob, meth, comb) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  add_column(k = 0, .before = 1)

rbind(es_0, es_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex",
                           bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: sample covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/sam_es.tex"), self_contained = FALSE)


## Appendix ----
### Shrinkage matrix ----
sel_prob <- c("ctjb", "ctshrh", "hshrh", "ctshroh", "hshroh")

#### CRPS ----
rbind(crps_0, crps_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |> select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: shrinkage covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/shr_crps.tex"), self_contained = FALSE)

#### ES ----
rbind(es_0, es_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex",
                           bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: shrinkage covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/shr_es.tex"), self_contained = FALSE)

### Complete set ----
sel_mc <- c("arima lev" = "base",
            "csbu shr" = "ct$(shr_{cs}, bu_{te})$",
            "csbu wls" = "ct$(wls_{cs}, bu_{te})$",
            "oct wlsv" = "oct$(wlsv)$",
            "oct bdshr" = "oct$(bdshr)$",
            "oct shr" = "oct$(shr)$",
            "oct hshr" = "oct$(hshr)$",
            "octo wlsv" = "oct$_o(wlsv)$",
            "octo bdshr" = "oct$_o(bdshr)$",
            "octo shr" = "oct$_o(shr)$",
            "octo hshr" = "oct$_o(hshr)$",
            "octoh shr" = "oct$_{oh}(shr)$",
            "octoh hshr" = "oct$_{oh}(hshr)$")

#### Sample matrix ----
sel_prob <- c("ctjb", "ctsamh", "hsamh", "ctsamoh", "hsamoh")

##### CRPS ----
rbind(crps_0, crps_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |> select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: sample covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/sam_crps_com.tex"), self_contained = FALSE)

##### ES ----
rbind(es_0, es_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex",
                           bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: sample covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/sam_es_com.tex"), self_contained = FALSE)


#### Shrinkage matrix ----
sel_prob <- c("ctjb", "ctshrh", "hshrh", "ctshroh", "hshroh")

##### CRPS ----
rbind(crps_0, crps_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |> select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: shrinkage covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/shr_crps_com.tex"), self_contained = FALSE)

##### ES ----
rbind(es_0, es_k) |>
  filter(prob %in% sel_prob, paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(k, meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob, k) |> mutate(ming = min(value)) |> ungroup() |>
  group_by(k) |> mutate(mink = min(value)) |> ungroup() |>
  mutate(bold = value == ming,
         blue = value == mink,
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex",
                           bold = bold, 
                           color = ifelse(red, "red", ifelse(blue, "blue", "black"))))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming, -mink) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "cccccc",
      col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
                    "\\multirow{-4}{*}{Bootstrap}", "G", "H", "G", "H"),
      caption = paste0("ES - selected"), 
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[4], start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("", "", "Multi-step residuals" = 2, "Overlapping and\nmulti-step residuals" = 2), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "Gaussian frameworks: shrinkage covariance matrix" = 4), 
                   escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 5), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(2:6, width = "1.5cm") |> column_spec(1, width = "2.5cm") |>
  save_kable(paste0("./Tables/shr_es_com.tex"), self_contained = FALSE)
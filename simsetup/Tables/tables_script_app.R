# Tables AusGDP ----
library(formattable)
library(kableExtra)
library(tidyverse)
load("./ProbScore/scores.RData")
load("./analisys/norm_cov.RData")

## Paper ----
sel_mc <- c("base NA" = "base",
            "csbu bu" = "ct$(bu)$",
            "csbu shr" = "ct$(shr_{cs}, bu_{te})$",
            "tebu wlsv" = "ct$(wlsv_{te}, bu_{cs})$",
            "oct wlsv" = "oct$(wlsv)$",
            "oct bdshr" = "oct$(bdshr)$",
            "oct shr" = "oct$(shr)$",
            "oct bshr" = "oct$(bshr)$",
            "oct hshr" = "oct$(hshr)$",
            "oct hbshr" = "oct$(hbshr)$",
            "octh shr" = "oct$_h(shr)$",
            "octh bshr" = "oct$_h(bshr)$",
            "octh hshr" = "oct$_h(hshr)$",
            "octh hbshr" = "oct$_h(hbshr)$")
sel_prob <- c("ctjb", "ctsam", "bsam", "hsam", "hbsam", "ctsamh", "bsamh", "hsamh", "hbsamh")
K <- setNames(c("$\\forall k \\\\in \\\\{2,1\\\\}$", 
                paste0("$k = ", c(1,2), "$")), c(0,1,2))

### Norm cov ----
data |>
  group_by(prob, meth, comb) |>
  summarise(value = mean(Frobenius), .groups = "drop") |>
  filter(prob %in% sel_prob,
         paste(meth, comb) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb), names(sel_mc), ordered = TRUE),
         meth2 = recode(meth2, !!!sel_mc),
         prob = factor(prob, sel_prob, ordered = TRUE)) |>
  arrange(meth2, prob) |>
  select(-meth, -comb) |>
  group_by(prob) |>
  mutate(ming = min(value)) |>
  ungroup() |>
  mutate(bold = value == ming,
         blue = value == min(value),
         red = value > 1,
         value = cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, 
                           color = ifelse(blue, "blue", "black")))|>
  ungroup() |>
  select(-bold, -blue, -red, -ming) |>
  pivot_wider(names_from = c(prob), names_sep = "-") |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "c", 
      col.names = c("", "",
                    "G", "B", "H", "HB", "G", "B", "H", "HB"),
      escape = FALSE) |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "In-sample residuals" = 4, "Multi-step residuals" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "\\\\makecell[c]{Gaussian approach: sample covariance matrix}" = 8), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 9), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(1, width = "2.5cm", latex_valign = "m") |>
  save_kable(paste0("./Tables/norm_cov_app.tex"), self_contained = FALSE)

### CRPS ----
crps_k <- left_join(df_crps_mean, 
                    df_crps_mean |>
                      filter(meth == "base", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "nn", "sd", "prob")) |>
                      rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  group_by(prob, k, meth, comb, h) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  filter(h == 0) |> select(-h) 
crps_0 <- left_join(df_crps_mean, 
                    df_crps_mean |>
                      filter(meth == "base", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "nn", "sd","prob")) |>
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
      align = "c", 
      col.names = c("", "",
                    "G", "B", "H", "HB", "G", "B", "H", "HB"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "In-sample residuals" = 4, "Multi-step residuals" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "\\\\makecell[c]{Gaussian approach: sample covariance matrix}" = 8), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 9), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(1, width = "2.5cm", latex_valign = "m") |>
  save_kable(paste0("./Tables/sam_crps_app.tex"), self_contained = FALSE)

### ES ----
es_k <- left_join(df_es_mean, 
                  df_es_mean |>
                    filter(meth == "base", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "nn", "sd", "prob")) |>
                    rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  filter(h == 0, serie == "all") |>
  group_by(prob, k, meth, comb) |>
  summarise(value = exp(mean(log(value))), .groups = "drop")
es_0 <- left_join(df_es_mean, 
                  df_es_mean |>
                    filter(meth == "base", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "nn", "sd", "prob")) |>
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
  arrange(k) |> select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "c", 
      col.names = c("", "",
                    "G", "B", "H", "HB", "G", "B", "H", "HB"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "In-sample residuals" = 4, "Multi-step residuals" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "\\\\makecell[c]{Gaussian approach: sample covariance matrix}" = 8), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 9), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(1, width = "2.5cm", latex_valign = "m") |>
  save_kable(paste0("./Tables/sam_es_app.tex"), self_contained = FALSE)

## Appendix ----
### Shrinkage matrix ----
sel_prob <- c("ctshr", "bshr", "hshr", "hbshr", "ctshrh", "bshrh", "hshrh", "hbshrh", "ctjb")

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
      align = "c", 
      col.names = c("", "",
                    "G", "B", "H", "HB", "G", "B", "H", "HB"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "In-sample residuals" = 4, "Multi-step residuals" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "\\\\makecell[c]{Gaussian approach: shrinkage covariance matrix}" = 8), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 9), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(1, width = "2.5cm", latex_valign = "m") |>
  save_kable(paste0("./Tables/shr_crps_app.tex"), self_contained = FALSE)


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
  arrange(k) |> select(-k) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "c", 
      col.names = c("", "",
                    "G", "B", "H", "HB", "G", "B", "H", "HB"),
      escape = FALSE) |>
  pack_rows(K[1], start_row = 1, end_row = length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[2], start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  pack_rows(K[3], start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), 
            indent = FALSE, escape = FALSE, latex_align = "c") |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "In-sample residuals" = 4, "Multi-step residuals" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "", "\\\\makecell[c]{Gaussian approach: shrinkage covariance matrix}" = 8), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 9), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(1, width = "2.5cm", latex_valign = "m") |>
  save_kable(paste0("./Tables/shr_es_app.tex"), self_contained = FALSE)
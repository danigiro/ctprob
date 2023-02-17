# Tables VN525 ----
library(formattable)
library(kableExtra)
library(tidyverse)
load("./ProbScore/arima_lev_scores.RData")

## Paper ----
sel_mc <- c("ets log NA" = "base",
            "csbu bu sntz" = "ct$(bu)$",
            "csbu shr sntz" = "ct$(shr_{cs}, bu_{te})$",
            "tebu wlsv sntz" = "ct$(wlsv_{te}, bu_{cs})$",
            "octf ols sntz" = "oct$(ols)$",
            "octf struc sntz" = "oct$(struc)$",
            "octf wlsv sntz" = "oct$(wlsv)$",
            "octf bdshr sntz" = "oct$(bdshr)$",
            "octh hbshr sntz" = "oct$_h(hbshr)$",
            "octh bshr sntz" = "oct$_h(bshr)$",
            "octh hshr sntz" = "oct$_h(hshr)$",
            "octh shr sntz" = "oct$_h(shr)$")
sel_prob <- c("ctjb", "ctsamh0", "bsamh0", "hsamh0", "hbsamh0")
K <- setNames(c("$\\\\forall k \\\\in \\\\{12,6,4,3,2,1\\\\}$", 
                paste0("$k = ", c(1,2,3,4,6,12), "$")), c(0,1,2,3,4,6,12))

### CRPS ----
crps_k <- left_join(df_crps_mean |>
                      filter(prob %in% sel_prob,
                             paste(meth, comb, nn) %in% names(sel_mc)),
                    df_crps_mean |>
                      filter(meth == "ets", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "prob", "nn")) |>
                      rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  group_by(prob, k, meth,nn,  comb, h) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  filter(h == 0) |> select(-h) 
crps_0 <- left_join(df_crps_mean|>
                      filter(prob %in% sel_prob, 
                             paste(meth, comb, nn) %in% names(sel_mc)),
                    df_crps_mean |>
                      filter(meth == "ets", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "nn","prob")) |>
                      rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  group_by(prob, meth, comb, h, nn) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  filter(h == 0) |> select(-h) |> add_column(k = 0, .before = 1)

options(knitr.kable.NA = '') 
rbind(crps_0, crps_k) |>
  #filter(prob %in% sel_prob, paste(meth, comb, nn) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb, nn), names(sel_mc), ordered = TRUE),
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
  mutate(colgk = recode(k, !!!setNames(c(0,0,0,0,1,1,1), c(0,2,4,12,1,3,6))),
         colgk2 = recode(k, !!!setNames(c(1,2,3,4,1,2,3), c(0,2,4,12,1,3,6)))) |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  pivot_wider(names_from = c(prob, colgk), names_sep = "-") |>
  select(-colgk2, -nn) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      #align = "c|ccccc|ccccc",
      align = "c",
      col.names = c("\\multicolumn{1}{c}{}", 
                    "", "G", "B","H", "\\multicolumn{1}{c}{HB}", "", "G", "B","H", "HB"),
      #col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
      #              "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB", "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB"),
      escape = FALSE)  |>
  footnote(general = paste0("\\\\rule{0pt}{1.75em}\\\\makecell[l]{$^\\\\ast$",
                            "The Gaussian method employs a sample covariance matrix and includes ",
                            "four techniques (G, B, H, HB)\\\\\\\\ with multi-step residuals.",
                            "}"), 
           escape = FALSE, general_title = "")|>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[1], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[2], "}"), 
            start_row = 1, end_row = length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[3], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[4], "}"), 
            start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[5], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[6], "}"), 
            start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[7], "}} & \\\\multicolumn{5}{c}{"), 
            start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  #add_header_above(c("", "", "Multi-step residuals" = 4, "", "Multi-step residuals" = 4), 
  #                 escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4, 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 10), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(6, border_right = T) |> 
  column_spec(2, border_left = T) |> 
  save_kable(paste0("./Tables/sam_crps_more.tex"), self_contained = FALSE)

### ES ----
es_k <- left_join(df_es_mean |>
                    filter(prob %in% sel_prob, 
                           paste(meth, comb, nn) %in% names(sel_mc)),
                  df_es_mean |>
                    filter(meth == "ets", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "prob", "nn")) |>
                    rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  filter(h == 0, serie == "all") |>
  group_by(prob, k, nn, meth, comb) |>
  summarise(value = exp(mean(log(value))), .groups = "drop")
es_0 <- left_join(df_es_mean |>
                    filter(prob %in% sel_prob, 
                           paste(meth, comb, nn) %in% names(sel_mc)),
                  df_es_mean |>
                    filter(meth == "ets", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "prob", "nn")) |>
                    rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  filter(h == 0, serie == "all") |>
  group_by(prob, meth,nn,  comb) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  add_column(k = 0, .before = 1)

rbind(es_0, es_k) |>
  filter(prob %in% sel_prob, paste(meth, comb, nn) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb, nn), names(sel_mc), ordered = TRUE),
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
  mutate(colgk = recode(k, !!!setNames(c(0,0,0,0,1,1,1), c(0,2,4,12,1,3,6))),
         colgk2 = recode(k, !!!setNames(c(1,2,3,4,1,2,3), c(0,2,4,12,1,3,6)))) |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  pivot_wider(names_from = c(prob, colgk), names_sep = "-") |>
  select(-colgk2, -nn) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      #align = "c|ccccc|ccccc",
      align = "c",
      col.names = c("\\multicolumn{1}{c}{}", 
                    "", "G", "B","H", "\\multicolumn{1}{c}{HB}", "", "G", "B","H", "HB"),
      #col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
      #              "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB", "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB"),
      escape = FALSE)  |>
  footnote(general = paste0("\\\\rule{0pt}{1.75em}\\\\makecell[l]{$^\\\\ast$",
                            "The Gaussian method employs a sample covariance matrix and includes ",
                            "four techniques (G, B, H, HB)\\\\\\\\ with multi-step residuals.",
                            "}"), 
           escape = FALSE, general_title = "")|>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[1], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[2], "}"), 
            start_row = 1, end_row = length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[3], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[4], "}"), 
            start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[5], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[6], "}"), 
            start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[7], "}} & \\\\multicolumn{5}{c}{"), 
            start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  #add_header_above(c("", "", "Multi-step residuals" = 4, "", "Multi-step residuals" = 4), 
  #                 escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4, 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 10), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(6, border_right = T) |> 
  column_spec(2, border_left = T) |> 
  save_kable(paste0("./Tables/sam_es_more.tex"), self_contained = FALSE)

## Appendix Shrinkage matrix ----
sel_prob <- c("ctjb", "ctshrh0", "bshrh0", "hshrh0", "hbshrh0")
K <- setNames(c("$\\\\forall k \\\\in \\\\{12,6,4,3,2,1\\\\}$", 
                paste0("$k = ", c(1,2,3,4,6,12), "$")), c(0,1,2,3,4,6,12))

### CRPS ----
crps_k <- left_join(df_crps_mean |>
                      filter(prob %in% sel_prob,
                             paste(meth, comb, nn) %in% names(sel_mc)),
                    df_crps_mean |>
                      filter(meth == "ets", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "prob", "nn")) |>
                      rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  group_by(prob, k, meth,nn,  comb, h) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  filter(h == 0) |> select(-h) 
crps_0 <- left_join(df_crps_mean|>
                      filter(prob %in% sel_prob, 
                             paste(meth, comb, nn) %in% names(sel_mc)),
                    df_crps_mean |>
                      filter(meth == "ets", prob == sel_prob[1]) |>
                      select(-c("meth", "comb", "nn","prob")) |>
                      rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  group_by(prob, meth, comb, h, nn) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  filter(h == 0) |> select(-h) |> add_column(k = 0, .before = 1)

options(knitr.kable.NA = '') 
rbind(crps_0, crps_k) |>
  #filter(prob %in% sel_prob, paste(meth, comb, nn) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb, nn), names(sel_mc), ordered = TRUE),
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
  mutate(colgk = recode(k, !!!setNames(c(0,0,0,0,1,1,1), c(0,2,4,12,1,3,6))),
         colgk2 = recode(k, !!!setNames(c(1,2,3,4,1,2,3), c(0,2,4,12,1,3,6)))) |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  pivot_wider(names_from = c(prob, colgk), names_sep = "-") |>
  select(-colgk2, -nn) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      #align = "c|ccccc|ccccc",
      align = "c",
      col.names = c("\\multicolumn{1}{c}{}", 
                    "", "G", "B","H", "\\multicolumn{1}{c}{HB}", "", "G", "B","H", "HB"),
      #col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
      #              "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB", "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB"),
      escape = FALSE)  |>
  footnote(general = paste0("\\\\rule{0pt}{1.75em}\\\\makecell[l]{$^\\\\ast$",
                            "The Gaussian method employs a shrikage covariance matrix and includes ",
                            "four techniques (G, B, H, HB)\\\\\\\\ with multi-step residuals.",
                            "}"), 
           escape = FALSE, general_title = "")|>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[1], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[2], "}"), 
            start_row = 1, end_row = length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[3], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[4], "}"), 
            start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[5], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[6], "}"), 
            start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[7], "}} & \\\\multicolumn{5}{c}{"), 
            start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  #add_header_above(c("", "", "Multi-step residuals" = 4, "", "Multi-step residuals" = 4), 
  #                 escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4, 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 10), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(6, border_right = T) |> 
  column_spec(2, border_left = T) |> 
  save_kable(paste0("./Tables/shr_crps_more.tex"), self_contained = FALSE)

### ES ----
es_k <- left_join(df_es_mean |>
                    filter(prob %in% sel_prob, 
                           paste(meth, comb, nn) %in% names(sel_mc)),
                  df_es_mean |>
                    filter(meth == "ets", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "prob", "nn")) |>
                    rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  filter(h == 0, serie == "all") |>
  group_by(prob, k, nn, meth, comb) |>
  summarise(value = exp(mean(log(value))), .groups = "drop")
es_0 <- left_join(df_es_mean |>
                    filter(prob %in% sel_prob, 
                           paste(meth, comb, nn) %in% names(sel_mc)),
                  df_es_mean |>
                    filter(meth == "ets", prob == sel_prob[1]) |>
                    select(-c("meth", "comb", "prob", "nn")) |>
                    rename(base = value), by = c("serie", "k", "h")) |>
  mutate(value = value/base) |>
  filter(h == 0, serie == "all") |>
  group_by(prob, meth,nn,  comb) |>
  summarise(value = exp(mean(log(value))), .groups = "drop") |>
  add_column(k = 0, .before = 1)

rbind(es_0, es_k) |>
  filter(prob %in% sel_prob, paste(meth, comb, nn) %in% names(sel_mc)) |>
  mutate(meth2 = factor(paste(meth, comb, nn), names(sel_mc), ordered = TRUE),
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
  mutate(colgk = recode(k, !!!setNames(c(0,0,0,0,1,1,1), c(0,2,4,12,1,3,6))),
         colgk2 = recode(k, !!!setNames(c(1,2,3,4,1,2,3), c(0,2,4,12,1,3,6)))) |>
  mutate(k = factor(k, 0:12, ordered = TRUE)) |>
  arrange(k) |>
  select(-k) |>
  pivot_wider(names_from = c(prob, colgk), names_sep = "-") |>
  select(-colgk2, -nn) |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      #align = "c|ccccc|ccccc",
      align = "c",
      col.names = c("\\multicolumn{1}{c}{}", 
                    "", "G", "B","H", "\\multicolumn{1}{c}{HB}", "", "G", "B","H", "HB"),
      #col.names = c("\\multirow{-5}{*}{\\parbox{2cm}{\\centering\\textbf{Reconciliation\\\\approach}}}", 
      #              "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB", "\\multirow{-4}{*}{ctjb}", "G", "B","H", "HB"),
      escape = FALSE)  |>
  footnote(general = paste0("\\\\rule{0pt}{1.75em}\\\\makecell[l]{$^\\\\ast$",
                            "The Gaussian method employs a shrikage covariance matrix and includes ",
                            "four techniques (G, B, H, HB)\\\\\\\\ with multi-step residuals.",
                            "}"), 
           escape = FALSE, general_title = "")|>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[1], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[2], "}"), 
            start_row = 1, end_row = length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[3], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[4], "}"), 
            start_row = length(sel_mc)+1, end_row = 2*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[5], "}} & \\\\multicolumn{5}{c}{\\\\textbf{", K[6], "}"), 
            start_row = 2*length(sel_mc)+1, end_row = 3*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  pack_rows(paste0("} & \\multicolumn{5}{c}{\\\\textbf{", K[7], "}} & \\\\multicolumn{5}{c}{"), 
            start_row = 3*length(sel_mc)+1, end_row = 4*length(sel_mc), colnum = 1,
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  #add_header_above(c("", "", "Multi-step residuals" = 4, "", "Multi-step residuals" = 4), 
  #                 escape = TRUE, line_sep = 0, line = FALSE) |>
  add_header_above(c("\\\\makecell[c]{\\\\bfseries Reconciliation\\\\\\\\\\\\bfseries approach}", 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4, 
                     "ctjb", "\\\\makecell[c]{Gaussian approach\\\\textsuperscript{*}}" = 4), 
                   escape = FALSE, line_sep = 0, line = FALSE) |>
  add_header_above(c("", "Base forecasts' sample approach" = 10), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  column_spec(6, border_right = T) |> 
  column_spec(2, border_left = T) |> 
  save_kable(paste0("./Tables/shr_es_more.tex"), self_contained = FALSE)

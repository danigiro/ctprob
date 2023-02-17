# Figures AusGDP ----
library(tidyverse)

load("./ProbScore/ets_log_scores.RData")

nemenyi_fun <- function(data){
  nemenyi <- tsutils::nemenyi(data, plottype = "none")
  df_plot <- full_join(as_tibble(nemenyi$means, rownames = "name"), 
                       full_join(rename(as_tibble(nemenyi$means-nemenyi$cd/2, rownames = "name"), "l" = "value"),
                                 rename(as_tibble(nemenyi$means+nemenyi$cd/2, rownames = "name"), "u" = "value"), 
                                 by = "name"), by = "name") |>
    arrange(value) |>
    mutate(#name = gsub(" ", "", name),
      name = paste0(name, " - ", format(round(value, 2), width = 5, nsmall = 2))) |>
    add_column(fpval = nemenyi$fpval,
               fH = nemenyi$fH)
  df_plot$col <- df_plot$l <= df_plot$u[1]
  
  as_tibble(df_plot)
}

## Paper nemenyi ----
sel_mc <- c("ets log NA" = "base",
            "csbu bu sntz" = "ct(bu)",
            "csbu shr sntz" = "ct(shr[cs],bu[te])",
            "tebu wlsv sntz" = "ct(wlsv[te],bu[cs])",
            "octf ols sntz" = "oct(ols)",
            "octf struc sntz" = "oct(struc)",
            "octf wlsv sntz" = "oct(wlsv)",
            "octf bdshr sntz" = "oct(bdshr)",
            "octh hbshr sntz" = "oct[h](hbshr)",
            "octh hshr sntz" = "oct[h](hshr)",
            "octh bshr sntz" = "oct[h](bshr)",
            "octh shr sntz" = "oct[h](shr)")
sel_prob <- c("ctjb", "hbsamh0")

nem <- df_crps_mean |>
  filter(h == 0, 
         paste(meth, comb, nn) %in% names(sel_mc),
         prob %in% sel_prob) |>
  mutate(name = paste(meth, comb, nn),
         name = recode(name, !!!sel_mc)) |>
  select(prob, serie, name, value, k) |>
  group_by(prob) |>
  nest() |>
  mutate(data = map(data, pivot_wider, names_from = name),
         data = map(data, function(x) x[, -c(1:2)]),
         data = map(data, nemenyi_fun)) |>
  unnest(cols = c(data)) |>
  ungroup() |>
  arrange(value) |>
  mutate(name = factor(name, unique(name), ordered = TRUE),
         pch_name = str_detect(name, "base"))

nemk <- df_crps_mean |>
  filter(h == 0, 
         paste(meth, comb, nn) %in% names(sel_mc),
         prob %in% sel_prob) |>
  mutate(name = paste(meth, comb, nn),
         name = recode(name, !!!sel_mc)) |>
  select(prob, serie, name, value, k) |>
  group_by(prob, k) |>
  nest() |>
  mutate(data = map(data, pivot_wider, names_from = name),
         data = map(data, function(x) x[, -c(1)]),
         data = map(data, nemenyi_fun)) |>
  unnest(cols = c(data)) |>
  ungroup() |>
  arrange(value) |>
  mutate(name = factor(name, unique(name), ordered = TRUE),
         pch_name = str_detect(name, "base"))

ctjb <- rbind(nemk |>
                filter(prob == "ctjb") |>
                mutate(facet = paste0("k==",k)),
              nem |>
                filter(prob == "ctjb") |>
                mutate(k = 0,
                       facet = 'k %in% group("{",list(12,6,4,3,2,1),"}")')) |>
  arrange(value) |>
  mutate(name = factor(name, unique(name), ordered = TRUE)) |>
  arrange(k) |>
  mutate(facet = factor(facet, unique(facet), ordered = TRUE)) |>
  ggplot() + 
  geom_rect(aes(xmin=l, xmax=u, fill = col), ymin=-Inf, ymax=Inf, alpha = 0.2, 
            data = function(x) summarise(group_by(x, facet), l = min(l), col = TRUE,
                                         u = min(u), .groups = "drop"))+
  geom_segment(aes(x = l, xend = u, yend = name, y = name)) + 
  geom_point(aes(x = l, y = name), pch = "|", size = 2) + 
  geom_point(aes(x = u, y = name), pch = "|", size = 2) + 
  geom_point(aes(x = value, fill = col, y = name, pch = pch_name), size = 3) +
  geom_label(data = function(x) select(x, facet, fpval) |>
               mutate(text = paste0("Friedman test p-value ", ifelse(fpval<0.001, " < 0.001", round(fpval, 3)))),
             aes(x = Inf, y = -Inf, label = text), vjust = "inward", hjust = "inward", size = 2.5,  label.size = NA) + 
  scale_shape_manual(values=c(21, 24))+
  facet_wrap(.~facet, ncol = 2, scales = "free", 
             labeller = label_parsed)+
  labs(y = NULL, x = NULL, subtitle = "Cross-temporal Join Bootstrap approach\n") + 
  theme_minimal()+
  scale_y_discrete(labels = scales::label_parse())+
  theme(legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        legend.margin = margin())

ggsave("./Figures/ctjb_more.pdf", ctjb,
       width = 7.5,
       height = 10)

hsamh <- rbind(nemk |>
                 filter(prob == "hbsamh0") |>
                 mutate(facet = paste0("k==",k)),
               nem |>
                 filter(prob == "hbsamh0") |>
                 mutate(k = 0,
                        facet = 'k %in% group("{",list(12,6,4,3,2,1),"}")')) |>
  arrange(value) |>
  mutate(name = factor(name, unique(name), ordered = TRUE)) |>
  arrange(k) |>
  mutate(facet = factor(facet, unique(facet), ordered = TRUE)) |>
  ggplot() + 
  geom_rect(aes(xmin=l, xmax=u, fill = col), ymin=-Inf, ymax=Inf, alpha = 0.2, 
            data = function(x) summarise(group_by(x, facet), l = min(l), col = TRUE,
                                         u = min(u), .groups = "drop"))+
  geom_segment(aes(x = l, xend = u, yend = name, y = name)) + 
  geom_point(aes(x = l, y = name), pch = "|", size = 2) + 
  geom_point(aes(x = u, y = name), pch = "|", size = 2) + 
  geom_point(aes(x = value, fill = col, y = name, pch = pch_name), size = 3) +
  geom_label(data = function(x) select(x, facet, fpval) |>
               mutate(text = paste0("Friedman test p-value ", ifelse(fpval<0.001, " < 0.001", round(fpval, 3)))),
             aes(x = Inf, y = -Inf, label = text), vjust = "inward", hjust = "inward", size = 2.5,  label.size = NA) + 
  scale_shape_manual(values=c(21, 24))+
  facet_wrap(.~facet, ncol =2, scales = "free", 
             labeller = label_parsed)+
  labs(y = NULL, x = NULL, subtitle = "Gaussian approach\n(multi-step residuals, HB)") + 
  theme_minimal()+
  scale_y_discrete(labels = scales::label_parse())+
  theme(legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        legend.margin = margin())

ggsave("./Figures/hbsamh_more.pdf", hsamh,
       width = 7.5,
       height = 10)

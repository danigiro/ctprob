library(tidyverse)
# Function ----
param_global <- function(n, m){
  kt <- FoReco::thf_tools(m = m)$kt
  ((n*kt)*((n*kt)+1))/2
}

param_hf_bts <- function(nb, m){
  ((nb*m)*((nb*m)+1))/2
}

param_hf <- function(n, m){
  ((n*m)*((n*m)+1))/2
}

param_bts <- function(nb, m){
  kt <- FoReco::thf_tools(m = m)$kt
  ((nb*kt)*((nb*kt)+1))/2
}

scientific_10 <- function(x) {
  parse(text=gsub("e\\+", "%*%10^", scales::scientific_format()(x)))
}

# Dataset grid ----
data_grid <- expand_grid(m = c(2,4,12, 24),
            n = 2*round(seq(5, 300, 10))+1)

# DATA number of parameters ----
data <- data_grid |>
  mutate(nb = lapply(n, function(x) setNames(pmax(c(x-1, trunc((x-1)*3/5), (x-1)/2), 2), c(1,2,3)))) |>
  unnest_longer(nb) |>
  unique() |>
  arrange(m, n, nb) |>
  mutate(global = mapply(param_global, n = n, m = m),
         hf_bts = mapply(param_hf_bts, nb = nb, m = m),
         hf = mapply(param_hf, n = n, m = m),
         bts = mapply(param_bts, nb = nb, m = m))

# Plot lines ----
p1 <- data |> 
  pivot_longer(-c("n", "nb", "m", "nb_id"))  |>
  mutate(m = factor(paste0("m==", m), paste0("m==", c(1:24)), ordered = TRUE),
         name = factor(name, c("global", "hf_bts", "hf", "bts"), ordered = TRUE),
         name = recode(factor(name), "global" = "Global",
                       "hf_bts" = "High frequency bottom time series",
                       "hf" = "High frequency",
                       "bts" = "Bottom time series")) |>
  filter(nb_id == 2) |>
  ggplot(aes(y = value, x = n, col = name, pch = factor(nb_id))) +
  geom_line() + 
  scale_y_continuous(labels = scientific_10) +
  facet_wrap(m~., scales = "free", ncol = 2, 
             labeller = labeller(m = label_parsed)) + 
  labs(x = "n, cross-sectional dimension", 
       y = "# of parameters") +
  scale_color_hue() +
  theme_minimal() +
  theme(legend.position = "top",
        text = element_text(size = 10),
        legend.title = element_blank())

# Plot boxplot ----
p2 <- data |> 
  pivot_longer(-c("n", "nb", "m", "nb_id"))  %>%
  left_join(filter(., name == "global") |>
              select(-name) |>
              rename(global = value)) |>
  filter(nb_id == 2, name != "global") |>
  mutate(name = factor(name, c("global", "hf_bts", "hf", "bts"), ordered = TRUE),
         name = recode(factor(name), "global" = "Global",
                       "hf_bts" = "High frequency bottom time series",
                       "hf" = "High frequency",
                       "bts" = "Bottom time series")) |>
  mutate(value = 1-value/global) |>

  ggplot(aes(y = value, x = factor(m), col = name)) +
  geom_boxplot()+ 
  labs(x = "m", 
       y = "% reduction r.t. Global") +
  scale_y_continuous(labels=scales::percent) + 
  scale_color_manual(values = scales::hue_pal()(4)[-1]) +
  theme_minimal() +
  theme(legend.position = "top",
        text = element_text(size = 10),
        legend.title = element_blank())

# Arrange p1 & p2 ----
ggpubr::ggarrange(p1, p2, ncol=2, nrow=1, widths = c(0.65, 0.35),
                  common.legend = TRUE, legend="bottom")



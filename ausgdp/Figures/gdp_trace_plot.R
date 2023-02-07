library(tidyverse)
library(FoReco)
library(zoo)

load("./DATA.RData")
gdp <- DATA[,1]
time_gdp <- as.yearqtr(time(gdp))
train <- tibble(value = DATA[40:42,1],
                time = time_gdp[40:42])

fore <- NULL
extract_forecasts <- function(prob, meth, scomb, sserie, sk, site){
  filename <- file.path(".", "ProbScore", "arima", "lev", paste0("rboot_", prob, "_", meth, ".RData"))
  obj <- mget(load(filename))$df_q
  obj |>
    filter(serie == sserie,
           k == sk, 
           ite %in% site,
           comb == scomb) |>
    add_column(prob = prob)
}

wind <- 14
train_start <- 40
idx_train <- seq(40, by = 1, length.out = wind+1)
train <- tibble(value = DATA[idx_train[length(idx_train)-c(2:0)],1],
                time = time_gdp[idx_train[length(idx_train)-c(2:0)]])
fore <- rbind(extract_forecasts(prob = "ctjb", meth = "octo", scomb = "wlsv", 
                                sserie = "Gdp", sk = 1, site = wind) |>
                add_column(time = time_gdp[max(idx_train)+c(0:3)]),
              extract_forecasts(prob = "ctjb", meth = "csbu", scomb = "wls", 
                                sserie = "Gdp", sk = 1, site = wind) |>
                add_column(time = time_gdp[max(idx_train)+c(0:3)]),
              extract_forecasts(prob = "ctjb", meth = "base", scomb = "lev", 
                                sserie = "Gdp", sk = 1, site = wind) |>
                add_column(time = time_gdp[max(idx_train)+c(0:3)]),
              extract_forecasts(prob = "hsamoh", meth = "octo", scomb = "wlsv", 
                                sserie = "Gdp", sk = 1, site = wind) |>
                add_column(time = time_gdp[max(idx_train)+c(0:3)]),
              extract_forecasts(prob = "hsamoh", meth = "csbu", scomb = "wls", 
                                sserie = "Gdp", sk = 1, site = wind) |>
                add_column(time = time_gdp[max(idx_train)+c(0:3)]),
              extract_forecasts(prob = "hsamoh", meth = "base", scomb = "lev", 
                                sserie = "Gdp", sk = 1, site = wind) |>
                add_column(time = time_gdp[max(idx_train)+c(0:3)]))


sel_mc <- c("arima lev" = "base",
            "csbu shr" = "ct(shr[cs], bu[te])",
            "csbu wls" = "ct(wls[cs], bu[te])",
            "octo wlsv" = "oct[o](wlsv)",
            "octo bdshr" = "oct[o](bdshr)",
            "octo shr" = "oct[o](shr)",
            "octo hshr" = "oct[o](hshr)",
            "octoh shr" = "oct[oh](shr)",
            "octoh hshr" = "oct[oh](hshr)")
sel_prob <- c("ctjb" = "Cross-temporal Joint Bootstrap",
              "hsamoh" = "Gaussian approach: overlapping and multi-step residuals, H")
scientific_10 <- function(x,suppress_ones=TRUE) {
  s <- parse_number(strsplit(scales::scientific(x), "\\+")[[1]][2])
  bquote(10^.(s))
}

scale <- 10000
gpdtrace <- fore |>
  mutate(name = paste(type, comb),
         name = recode(name, !!!sel_mc),
         prob = recode(prob, !!!sel_prob)) |>
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin = `10%`/scale, ymax = `90%`/scale, fill = "IC 80%"), linetype = 0, alpha = 0.15) + 
  geom_line(aes(y = test/scale, col = "obs"), lty = 2) +
  geom_point(aes(y = test/scale, col = "obs", pch = "obs")) +
  geom_line(aes(y = `50%`/scale, col = "median"))+
  geom_point(aes(y = `50%`/scale, col = "median", pch = "median")) +
  geom_line(aes(y = mean/scale, col = "mean"))+
  geom_point(aes(y = mean/scale, col = "mean", pch = "mean")) +
  geom_line(aes(y = value/scale, col = "obs"), data = train, lty = 2) + 
  scale_color_manual(values = c(4, 2, 1))+ 
  scale_fill_manual(values = c(scales::hue_pal()(3)[2])) +
  theme_minimal() +
  labs(y = bquote("Australian Gross Domestic Product (GDP)"%*%.(scientific_10(scale))~
                    "Iteration F.E. #"*.(wind)), 
       x = NULL) +
  scale_shape_manual(values = c(0:2), guide = NULL) +
  facet_grid(name~prob, labeller = labeller(name = label_parsed, 
                                            prob = label_wrap_gen(width = 40,multi_line = TRUE))) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 10),
        strip.text = element_text(size = 10))

ggsave("./Figures/gdptrace.pdf", plot = gpdtrace,
       width = 7.5,
       height = 5)
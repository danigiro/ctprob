library(tidyverse)
library(Matrix)
m <- 2
C <- t(c(1,1))
n <- 3
tools <- FoReco::ctf_tools(m = m, C = C)
Sct <- tools$ctf$Fmat
Ste <- tools$thf$R
Scs <- tools$hts$S
h <- rep(c(1,1,2), 3)
k <- rep(c(2,1,1), 3)
ns <- rep(c(1:3), each = 3)

# Covariance ----
omega_base <- matrix(1, length(h), length(h))
omega_1 <- omega_base
diag(omega_1) <- 1

# Function ----
cov2tibble <- function(cov){
  as_tibble(as.matrix(cov), rownames = "row") |>
    pivot_longer(-row, names_to = "col") |>
    mutate_if(is.character, readr::parse_number)
}

global_shr0 <- function(lambda, omega){
  cov <- lambda*diag(diag(omega)) + (1-lambda)*omega
  diag(cov) <- 2
  cov2tibble(cov)
}

hf_shr0 <- function(lambda, omega){
  omega_hf <- omega[k == 1, k == 1]
  cov1 <- lambda*diag(diag(omega_hf))+(1-lambda)*omega_hf
  cov <- kronecker(diag(n), Ste) %*% cov1 %*%t(kronecker(diag(n), Ste))
  cov[cov!=0] <- 1
  diag(cov[k == 1, k == 1]) <- 2
  cov2tibble(cov)
}

hf_bts_shr0 <- function(lambda, omega){
  omega_hf_bts <- omega[k == 1 & ns != 1, k == 1 & ns != 1]
  cov1 <- lambda*diag(diag(omega_hf_bts))+(1-lambda)*omega_hf_bts
  cov <- Sct %*% cov1 %*%t(Sct)
  cov[cov!=0] <- 1
  diag(cov[k == 1 & ns > 1, k == 1 & ns > 1]) <- 2
  cov2tibble(cov)
}

bts_shr0 <- function(lambda, omega){
  omega_bts <- omega[ns != 1, ns != 1]
  cov1 <- lambda*diag(diag(omega_bts))+(1-lambda)*omega_bts
  cov <-  kronecker(Scs, diag(m+1)) %*% cov1 %*%t(kronecker(Scs, diag(m+1)) )
  
  cov[cov!=0] <- 1
  diag(cov[ns > 1, ns > 1]) <- 2
  cov2tibble(cov)
}

# Dataset lambda = 0 ----
data0 <- tibble(lambda = 0) |>
  mutate(hf = lapply(lambda, function(x){
    mask <- omega_base
    mask[k == 1, k == 1] <- -1
    diag(mask[k == 1, k == 1]) <- 2
    cov2tibble(mask)
  }), 
  global = lapply(lambda, function(x){
    mask <- omega_base
    mask[, ] <- -1
    diag(mask) <- 2
    cov2tibble(mask)
  }), 
  hf_bts = lapply(lambda, function(x){
    mask <- omega_base
    mask[k == 1 & ns > 1, k == 1 & ns > 1] <- -1
    diag(mask[k == 1 & ns > 1, k == 1 & ns > 1]) <- 2
    cov2tibble(mask)
  }), 
  bts = lapply(lambda, function(x){
    mask <- omega_base
    mask[ns > 1, ns > 1] <- -1
    diag(mask[ns > 1, ns > 1]) <- 2
    cov2tibble(mask)
  })) |>
  pivot_longer(-lambda, values_to = "cov") |>
  unnest(cov)

# Dataset lambda = 1 ----
data <- tibble(lambda = 1) |>
  mutate(
  global = lapply(lambda, global_shr0, omega = omega_1),
  hf = lapply(lambda, hf_shr0, omega = omega_1),
  hf_bts = lapply(lambda, hf_bts_shr0, omega = omega_1),
  bts = lapply(lambda, bts_shr0, omega = omega_1)) |>
  pivot_longer(-lambda, values_to = "cov") |>
  unnest(cov)

rbind(data0, data) |>
  mutate(lambda = paste0("lambda==", lambda),
         name = factor(name, c("global", "hf_bts", "hf", "bts"), ordered = TRUE),
         name = recode(factor(name), "global" = "Global",
                       "hf_bts" = "High frequency bottom time series",
                       "hf" = "High frequency",
                       "bts" = "Bottom time series")) |>
  ggplot() + 
  geom_tile(aes(x = col, y = row, fill = factor(value)), color = "black")+
  #geom_tile(aes(x = col, y = row, col = factor(mask)), size = 1, alpha = 0)+
  #scale_fill_gradient2()+
  #geom_text(aes(x = col, y = row, label = round(value, 2)), size = 2.5) +
  scale_y_reverse()+
  labs(y = NULL, x = NULL) + 
  facet_grid(lambda~name, labeller = labeller(lambda = label_parsed, .multi_line = TRUE,
                                              name = label_wrap_gen(width = 20))) +
  scale_fill_manual(values = c("#bfdbf7", "white", "#6096ba", "#022b3a")) + 
  coord_fixed(expand = 0.05)+
  theme_void() + 
  theme(#axis.text.y.left = element_text(size = 8),
    #axis.ticks.length.y.left = unit(0.5, "lines"),
    strip.text.y = element_text(angle = -90),
    legend.position = "none",
    legend.title = element_blank(),
    plot.margin = margin())


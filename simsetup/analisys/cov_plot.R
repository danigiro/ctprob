library(tidyverse)
library(Matrix)

PhiB <- c(1.3446416, -0.7392638)
PhiC <- c(0.9454175, -0.4218556)

hfbts_cov <- function(Sigma, PhiB, PhiC){
  s11 <- Sigma[1,1]
  s21 <- PhiB[1]*Sigma[1,1]
  s31 <- Sigma[2,1]
  s41 <- PhiC[1]*Sigma[2,1]
  s22 <- Sigma[1,1]*(1+PhiB[1]^2)
  s32 <- PhiB[1]*Sigma[2,1]
  s42 <- Sigma[2,1]*(1+PhiC[1]*PhiB[1])
  s33 <- Sigma[2,2]
  s34 <- PhiC[1]*Sigma[2,2]
  s44 <- Sigma[2,2]*(1+PhiC[1]^2)
  matrix(c(s11,s21,s31,s41,
           s21,s22,s32,s42,
           s31,s32,s33,s34,
           s41,s42,s34,s44), 4, byrow = TRUE)
}
S <- FoReco::ctf_tools(C = t(c(1,1)), m = 2)$ctf$Fmat

Sigma <- diag(c(0.9, 1.8))%*%matrix(c(1, -0.8,
                                      -0.8, 1), 2, 2)%*%diag(c(0.9, 1.8))
cov <- S%*%hfbts_cov(Sigma = Sigma, PhiB = PhiB, PhiC = PhiC)%*%t(S)
Sigma2 <- diag(c(0.9, 1.8))%*%matrix(c(1, 0.6,
                                       0.6, 1), 2, 2)%*%diag(c(0.9, 1.8))
cov2 <- S%*%hfbts_cov(Sigma = Sigma2, PhiB = PhiB, PhiC = PhiC)%*%t(S)

cov2tibble <- function(cov){
  as_tibble(as.matrix(cov), rownames = "row") |>
    pivot_longer(-row, names_to = "col") |>
    mutate_if(is.character, readr::parse_number)
}

p1 <- cov2tibble(cov) |>
        add_column(rho = "Covariance matrix")|>
  ggplot(aes(x = col, y = row, fill = value)) + 
  geom_tile(col = "black")+
  scale_fill_gradient2()+
  geom_text(aes(label = round(value, 2)), size = 2.5) +
  scale_y_reverse(breaks = seq(1:9), labels = function(x) paste0(rep(LETTERS[1:3], each = 3), ", k = ", rep(c(2, 1, 1), 3), ", h = ", rep(c(1, 1, 2), 3))[x])+
  labs(y = NULL, x = NULL) + 
  facet_grid(.~rho) +
  coord_fixed(expand = 0.05)+
  theme_void() + 
  theme(axis.text.y.left = element_text(size = 8),
        #axis.ticks.length.y.left = unit(0.5, "lines"),
        legend.title = element_blank(),
        legend.margin = margin(0, 10, 0, 0),
        plot.title = element_text(hjust = 0.5, face = "bold", margin=margin(0,0,10,0)))

p2 <- cov2tibble(cov2cor(as.matrix(cov))) |>
  add_column(rho = "Correlation matrix") |>
  ggplot(aes(x = col, y = row, fill = value)) + 
  geom_tile(col = "black")+
  scale_fill_gradient2()+
  geom_text(aes(label = round(value, 2)), size = 2.5) +
  scale_y_reverse(breaks = seq(1:9), labels = function(x) paste0(rep(LETTERS[1:3], each = 3), ", k = ", rep(c(2, 1, 1), 3), ", h = ", rep(c(1, 1, 2), 3))[x])+
  labs(y = NULL, x = NULL) + 
  facet_grid(.~rho) +
  coord_fixed(expand = 0.05)+
  theme_void() + 
  theme(axis.text.y.left = element_text(size = 8),
        #axis.ticks.length.y.left = unit(0.5, "lines"),
        legend.title = element_blank(),
        legend.margin = margin(0, 10, 0, 0),
        plot.title = element_text(hjust = 0.5, face = "bold", margin=margin(0,0,10,0)))

ggsave("./analisys/covcor.pdf",
       gridExtra::grid.arrange(p1, p2, ncol = 2),
       width = 10,
       height = 5)

load("./analisys/base_cov_plot.RData")

c1 <- rbind(ctjb_base, #ctsam_base, 
            ctsamh_base) |>
               filter(type == "cov") |>
               mutate(type = recode(type, "cor" = "Correlation matrix",
                                    "cov" = "Covariance matrix"),
                      prob = recode(prob, "ctjb" = "Bootstrap", 
                                    "ctsam" = "Gaussian with\nin-sample residuals", 
                                    "ctsamh" = "Gaussian with\nmulti-step residuals")) |>
               ggplot(aes(x = row, y = col, fill = value)) +
               geom_tile(col = "black") +
               scale_fill_gradient2() +
               scale_y_reverse() +
               geom_text(aes(label = round(value, 1)), size = 2.5) +
               theme_void() +
               facet_grid(prob ~ type, switch="y") +
               coord_fixed()+
               theme(legend.position = "bottom",
                     legend.title = element_blank(),
                     strip.text.y = element_text(angle = -90))
c2 <- rbind(ctjb_base, #ctsam_base, 
            ctsamh_base) |>
               filter(type == "cor") |>
               mutate(type = recode(type, "cor" = "Correlation matrix",
                                    "cov" = "Covariance matrix"),
                      prob = recode(prob, "ctjb" = "Bootstrap", 
                                    "ctsam" = "Gaussian with\nin-sample residuals", 
                                    "ctsamh" = "Gaussian with\nmulti-step residuals")) |>
               ggplot(aes(x = row, y = col, fill = value)) +
               geom_tile(col = "black") +
               scale_fill_gradient2() +
               scale_y_reverse() +
               geom_text(aes(label = round(value, 1)), size = 2.5) +
               theme_void() +
               facet_grid(prob ~ type) +
               coord_fixed()+
               theme(legend.position = "bottom",
                     legend.title = element_blank(),
                     strip.text.y = element_text(angle = -90))

ggsave("./analisys/base_cov.pdf",
       gridExtra::grid.arrange(c1, c2, ncol = 2),
       width = 7,
       height = 7)


c1_app <- rbind(ctsam_base) |>
  filter(type == "cov") |>
  mutate(type = recode(type, "cor" = "Correlation matrix",
                       "cov" = "Covariance matrix"),
         prob = recode(prob, "ctjb" = "Bootstrap", 
                       "ctsam" = "Gaussian with\nin-sample residuals", 
                       "ctsamh" = "Gaussian with\nmulti-step residuals")) |>
  ggplot(aes(x = row, y = col, fill = value)) +
  geom_tile(col = "black") +
  scale_fill_gradient2() +
  scale_y_reverse() +
  geom_text(aes(label = round(value, 1)), size = 2.5) +
  theme_void() +
  facet_grid(prob ~ type, switch="y") +
  coord_fixed()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = -90))
c2_app <- rbind(ctsam_base) |>
  filter(type == "cor") |>
  mutate(type = recode(type, "cor" = "Correlation matrix",
                       "cov" = "Covariance matrix"),
         prob = recode(prob, "ctjb" = "Bootstrap", 
                       "ctsam" = "Gaussian with\nin-sample residuals", 
                       "ctsamh" = "Gaussian with\nmulti-step residuals")) |>
  ggplot(aes(x = row, y = col, fill = value)) +
  geom_tile(col = "black") +
  scale_fill_gradient2() +
  scale_y_reverse() +
  geom_text(aes(label = round(value, 1)), size = 2.5) +
  theme_void() +
  facet_grid(prob ~ type) +
  coord_fixed()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = -90))

ggsave("./analisys/base_cov_app.pdf",
       gridExtra::grid.arrange(c1_app, c2_app, ncol = 2),
       width = 7,
       height = 4)

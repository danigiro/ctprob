library(future.apply)
library(progressr)
handlers("progress")
plan("multisession", workers = 10)

data <- tsibbledata::pelt
hare <- ts(data$Hare, start = data$Year[1])
lynx <- ts(data$Lynx, start = data$Year[1])
Bcoef <- arima(lynx, c(2, 0, 0))$coef[1:2]
Ccoef <- arima(hare, c(2, 0, 0))$coef[1:2]

init <- 200
train <- 1000
h <- 2
N <- init + train + h
rho <- -0.8
Sigma <- diag(c(0.9, 1.8))%*%matrix(c(1, rho,
                                      rho, 1), 2, 2)%*%diag(c(0.9, 1.8))

dir.create(file.path(".","BaseForecasts"), recursive = TRUE, showWarnings = FALSE)

simulation <- function(id, Sigma, N, init){
  E <- MASS::mvrnorm(n = N, rep(0, 2), Sigma = Sigma)
  
  B <- arima.sim(list(ar = Bcoef), n = N, innov = E[,1])
  C <- arima.sim(list(ar = Ccoef), n = N, innov = E[,2])
  A <- B + C
  listObs <- NULL
  listObs$k1 <- ts(cbind(A, B, C)[-c(1:init),], frequency = 2)
  
  listObs$k2 <- ts(apply(listObs$k1, 2, function(x){
    tapply(x, rep(1:(length(x)/2), each = 2), sum)
  }), frequency = 1)
  m <- 2
  sim <- mapply(function(data, k){
    require(forecast)
    Fit <- list()
    Res <- list()
    Base <- list()
    Test <- list()
    Scale <- list()
    
    data <- ts(data, frequency = m/k)
    train <- window(data, end = trunc(tsp(data)[2])-k/m)
    Test <- window(data, start = trunc(tsp(data)[2]))
    Fit <- lapply(train, auto.arima,
                  start.p = 15, start.q = 2, start.P = 0, start.Q = 0,
                  seasonal = FALSE,
                  d = 0, D = 0, max.p = 20,
                  max.q = 0, max.P = 0, max.Q = 0, max.order = 20,
                  max.d = 0, max.D = 0, allowdrift = FALSE, allowmean = FALSE)
    Res <- sapply(Fit, function(x) x$x - fitted(x))
    Base <- sapply(Fit, function(x) forecast(x, h = m/k)$mean)
    Scale <- sapply(Fit, function(x) mean(abs(diff(x$x, lag = m/k, differences = 1)), na.rm = TRUE))
    
    return(list(train = train, 
                test = Test,
                fit = Fit,
                res = Res, 
                base = Base, 
                scale = Scale))
  }, data = listObs, k = c(1, 2), SIMPLIFY = FALSE)
  
  listFit <- lapply(sim, function(x) x$fit)
  listBase <- lapply(sim, function(x) x$base)
  listRes <- lapply(sim, function(x) x$res)
  listTest <- lapply(sim, function(x) x$test)
  listTrain <- lapply(sim, function(x) x$train)
  listScale <- lapply(sim, function(x) x$scale)
  
  itername <- paste0("ite_", formatC(id, width = nchar(1000+1), 
                                     format = "d", flag = "0"), ".RData")
  save(listFit, listBase, listRes, listTest, listTrain, listScale,
       file = file.path(".","BaseForecasts", itername))
  
  p()
}

idsim <- 1:500
seed_list <- future_lapply(seq_along(idsim), FUN = function(x) .Random.seed,
                           future.chunk.size = Inf, future.seed = 93L)
with_progress({
  p <- progressor(along = idsim)
  out <- future.apply::future_lapply(idsim, simulation, N = N, init = init, Sigma = Sigma, 
                                     future.seed = seed_list)
})

C <- t(c(1,1))
rownames(C) <-  "A"
colnames(C) <- c("B", "C")
K <- c(1, 2)
m <- 2
H <- 2
save(C, K, H, m, file = "DATA.RData")








hfitted <- function(object, ...) {
  UseMethod("hfitted")
}

hfitted.Arima <- function(object, h, ...) {
  if(h == 1){
    return(object$fitted)
  }
  y <- object$fitted+residuals(object, "innovation")
  yx <- residuals(object, "regression")
  # Get fitted model
  mod <- object$model
  # Reset model to initial state
  mod <-  stats::makeARIMA(mod$phi, mod$theta, mod$Delta)
  # Calculate regression component
  xm <- y-yx
  # mod_seq <- mod
  fits <- rep_len(NA_real_, length(y))
  
  start <- length(mod$Delta) + 1
  end <- length(yx) - h
  idx <- if(start > end) integer(0L) else start:end
  for(i in idx) {
    fc_mod <- attr(stats::KalmanRun(yx[seq_len(i)], mod, update = TRUE), "mod")
    fits[i + h] <- stats::KalmanForecast(h, fc_mod)$pred[h] + xm[i+h]
  }
  tsp(fits) <- tsp(object$x)
  fits
}

hfitted.ets <- function(object, h, ...) {
  if(h == 1){
    return(object$fitted)
  }else{
    forecast:::hfitted(object, h = h)
  }
  # errortype <- object$components[1]
  # trendtype <- object$components[2]
  # seasontype <- object$components[3]
  # damped <- object$components[4]
  # fc_class <- if (errortype == "A" && trendtype %in% c("A", "N") && seasontype %in% c("N", "A")) {
  #   fable:::ets_fc_class1
  # } else if (errortype == "M" && trendtype %in% c("A", "N") && seasontype %in% c("N", "A")) {
  #   fable:::ets_fc_class2
  # } else if (errortype == "M" && trendtype != "M" && seasontype == "M") {
  #   fable:::ets_fc_class3
  # } else {
  #   abort(sprintf("Multi-step fits for %s%s%s%s ETS models is not supported."),
  #         errortype, trendtype, if(damped) "d" else "", seasontype)
  # }
  # 
  # n <- nrow(object$states)-1
  # fits <- rep_len(NA_real_, n)
  # for(i in seq_len(n-h+1)) {
  #   fits[i + h - 1] <- fc_class(
  #     h = h,
  #     last.state = as.numeric(object$states[i, NCOL(object$states)]),
  #     trendtype, seasontype, damped, object$m, object$sigma2,
  #     object$par
  #   )$mu[h]
  # }
  # tsp(fits) <- tsp(object$x)
  # fits
}

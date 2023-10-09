#!/bin/bash

# Base forecast model (ets or arima)
model="ets" 

# Base forecast transformation:
# - log: logarithm transformation
# - lev: no transformation
trans="log" 

echo "Point base forecasts"
Rscript ./BaseForecasts/base_forecasts.R $model $trans 

echo "h-step residuals"
Rscript ./HfittedRes/hfitted_residuals.R $model $trans 
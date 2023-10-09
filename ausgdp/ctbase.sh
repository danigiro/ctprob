#!/bin/bash

# Base forecast model (ets or arima)
model="arima"

# Base forecast transformation:
# - log: logarithm transformation
# - lev: no transformation
trans="lev"

echo "Point base forecasts"
Rscript ./BaseForecasts/base_forecasts.R $model $trans 

echo "h-step residuals"
Rscript ./HfittedRes/hfitted_residuals.R $model $trans 

echo "1-step overlapping residuals"
Rscript ./OverlapRes/overlap_residuals.R $model $trans 

echo "h-step overlapping residuals"
Rscript ./HOverlapRes/hoverlap_residuals.R $model $trans 
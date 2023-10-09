#!/bin/bash

# Base forecast model (ets or arima)
model="arima"

# Base forecast transformation:
# - log: logarithm transformation
# - lev: no transformation
trans="lev"

# Cross-temporal join bootstrap
meth_list="ctjb"

# Base forecasts
do_base=true        # samples
do_base_score=true  # scores

# Reconciled forecasts
do_oct=true         # ct reconciliation (1-step residuals)
do_octo=true        # ct reconciliation (1-step overlapping residuals)
do_octh=true        # ct reconciliation (h-step residuals)
do_octoh=true       # ct reconciliation (h-step overlapping residuals)
do_ctbu=true        # cs reconciliation + te bottom-up

for meth in $meth_list
do
  echo $model $trans $meth
  
  # Base forecasts samples
  if [ "do_base" = true ]
  then
    echo "base"
    Rscript ./ProbReco/bboot_forecasts.R $model $trans $meth
  fi
  
  # Base forecasts scores
  if [ "$do_base_score" = true ]
  then
    echo "base score"
    Rscript ./ProbReco/bboot_scores.R $model $trans $meth
  fi
  
  # Reconciled forecasts samples and scores
  ## ct reconciliation (1-step residuals)
  if [ "$do_oct" = true ]
  then
    echo "oct"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth
  fi
  
  ## ct reconciliation (1-step overlapping residuals)
  if [ "$do_octo" = true ]
  then
    echo "octo"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth "o"
  fi
  
  ## ct reconciliation (h-step residuals)
  if [ "$do_octh" = true ]
  then
    echo "octh"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth "h"
  fi
  
  ## ct reconciliation (h-step overlapping residuals)
  if [ "$do_octoh" = true ]
  then
    echo "octoh"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth "oh"
  fi
  
  ## cs reconciliation + te bottom-up
  if [ "$do_ctbu" = true ]
  then
    echo "csbu"
    Rscript ./ProbReco/rboot_scores_ctbu_cs.R $model $trans $meth
  fi
  
  echo $model $trans $meth
done

exec bash
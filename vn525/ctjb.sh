#!/bin/bash

# Base forecast model (ets or arima)
model="ets" 

# Base forecast transformation:
# - log: logarithm transformation
# - lev: no transformation
trans="log" 

# Cross-temporal join bootstrap
meth_list="ctjb"

# Non-negative forecast reconciliation (sntz or free)
nn="sntz"

# Base forecasts
do_base=true        # samples
do_base_score=true  # scores

# Reconciled forecasts 
do_oct=true         # ct reconciliation (h-step residuals)
do_octf=true        # ct reconciliation (1-step residuals)
do_csbu=true        # cs reconciliation + te bottom-up
do_tebu=true        # te reconciliation + cs bottom-up

for meth in $meth_list
do
  echo $meth
  
  # Base forecasts samples
  if [ "$do_base" = true ]
  then
    echo "base"
    Rscript ./ProbReco/bboot_forecasts.R $model $trans $meth
  fi
  
  # Base forecasts scores
  if [ "$do_base_score" = true ]
  then
    echo "score"
    Rscript ./ProbReco/bboot_scores.R $model $trans $meth
  fi
  
  # Reconciled forecasts samples and scores
  ## cs reconciliation + te bottom-up
  if [ "$do_csbu" = true ]
  then
    echo "csbu"
    Rscript ./ProbReco/rboot_scores_ctbu_cs.R $model $trans $meth $nn
  fi
  
  ## te reconciliation + cs bottom-up
  if [ "$do_tebu" = true ]
  then
    echo "tebu"
    Rscript ./ProbReco/rboot_scores_ctbu_te.R $model $trans $meth $nn
  fi
  
  ## ct reconciliation (h-step residuals)
  if [ "$do_oct" = true ]
  then
    echo "oct"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth "h" $nn
  fi
  
  ## ct reconciliation (1-step residuals)
  if [ "$do_octf" = true ]
  then
    echo "octf"
    Rscript ./ProbReco/rboot_scores_octf.R $model $trans $meth "in" $nn
  fi
done

exec bash
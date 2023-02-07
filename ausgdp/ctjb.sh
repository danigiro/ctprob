#!/bin/bash
model="arima"
trans="lev"
meth_list="ctjb ctjbo"
do_base_gauss=false
do_base_boot=true
do_base_score=true
do_oct=true
do_octo=true
do_octh=true
do_octoh=true
do_cs=true
do_ctbu=true
do_te=true
do_tebu=false

for meth in $meth_list
do
  echo $model $trans $meth
  
  if [ "$do_base_boot" = true ]
  then
    echo "base"
    Rscript ./ProbReco/bboot_forecasts.R $model $trans $meth
  fi
  
  if [ "$do_base_gauss" = true ]
  then
    echo "base"
    Rscript ./ProbReco/bgauss_forecasts.R $model $trans $meth
  fi
  
  if [ "$do_base_score" = true ]
  then
    echo "base score"
    Rscript ./ProbReco/bboot_scores.R $model $trans $meth
  fi
  
  if [ "$do_oct" = true ]
  then
    echo "oct"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth
  fi
  
  if [ "$do_octo" = true ]
  then
    echo "octo"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth "o"
  fi
  
  if [ "$do_octh" = true ]
  then
    echo "octh"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth "h"
  fi
  
  if [ "$do_octoh" = true ]
  then
    echo "octoh"
    Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth "oh"
  fi
  
  if [ "$do_cs" = true ]
  then
    echo "cs"
    Rscript ./ProbReco/rboot_scores_cs.R $model $trans $meth
  fi
  
  if [ "$do_ctbu" = true ]
  then
    echo "csbu"
    Rscript ./ProbReco/rboot_scores_ctbu_cs.R $model $trans $meth
  fi
  
  if [ "$do_te" = true ]
  then
    echo "te"
    Rscript ./ProbReco/rboot_scores_te.R $model $trans $meth
  fi
  
  if [ "$do_tebu" = true ]
  then
    echo "tebu"
    Rscript ./ProbReco/rboot_scores_ctbu_te.R $model $trans $meth
  fi
  
  echo $model $trans $meth
done

exec bash
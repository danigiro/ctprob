#!/bin/bash
meth_list="ctshr hbshr hshr bshr"
# Gaussian incoherent covariance matrix
## ctsam: global covariance matrix
## hbsam: high frequency bottom time series covariance matrix
## hsam: high frequency covariance matrix
## bsam: bottom time series covariance matrix
res_list="in o h oh"
# Residuals for the incoherent covariance matrix
## in: in-sample residuals
## o: overlapping residuals
## h: multi-step residuals
## oh: multi-step and overlapping residuals
do_base_gauss=true  # base forecasts: true or false
do_base_score=true  # base forecasts accuracy: true or false
do_oct=true         # optimal ct reco (in-sample): true or false
do_octo=true        # optimal ct reco (overlapping): true or false
do_octh=true        # optimal ct reco (multi-sample): true or false
do_octoh=true       # optimal ct reco (multi-sample + overlapping): true or false
do_cs=true          # optimal cs reco: true or false
do_ctbu=true        # optimal cs reco + te bu: true or false
do_te=true          # optimal te reco: true or false
do_tebu=true        # optimal te reco + cs bu: true or false

for meth in $meth_list
do
  for res in $res_list
  do
    if [ $res == "in" ]
    then
      res=""
    fi
    
    echo $meth$res
    
    if [ "$do_base_gauss" = true ]
    then
      echo "base"
      Rscript ./ProbReco/bgauss_forecasts.R $meth $res
    fi
    
    if [ "$do_base_score" = true ]
    then
      echo "base score"
      Rscript ./ProbReco/bboot_scores.R $meth$res
    fi
    
    if [ "$do_oct" = true ]
    then
      echo "oct"
      Rscript ./ProbReco/rboot_scores_oct.R $meth$res
    fi
    
    if [ "$do_octo" = true ]
    then
      echo "octo"
      Rscript ./ProbReco/rboot_scores_oct.R $meth$res "o"
    fi
    
    if [ "$do_octh" = true ]
    then
      echo "octh"
      Rscript ./ProbReco/rboot_scores_oct.R $meth$res "h"
    fi
    
    if [ "$do_octoh" = true ]
    then
      echo "octoh"
      Rscript ./ProbReco/rboot_scores_oct.R $meth$res "oh"
    fi
    
    if [ "$do_cs" = true ]
    then
      echo "cs"
      Rscript ./ProbReco/rboot_scores_cs.R $meth$res
    fi
    
    if [ "$do_ctbu" = true ]
    then
      echo "csbu"
      Rscript ./ProbReco/rboot_scores_ctbu_cs.R $meth$res
    fi
    
    if [ "$do_te" = true ]
    then
      echo "te"
      Rscript ./ProbReco/rboot_scores_te.R $meth$res
    fi
    
    if [ "$do_tebu" = true ]
    then
      echo "tebu"
      Rscript ./ProbReco/rboot_scores_ctbu_te.R $meth$res
    fi
  done
  echo $meth
done

exec bash
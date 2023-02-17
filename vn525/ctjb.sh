#!/bin/bash
model="ets"
trans="log"
meth_list="ctjb"
res_list="in" # in
nn="free" # sntz or free
do_base=true
do_base_score=true
do_oct=true
do_octf=true
do_csbu=true
do_tebu=true

for meth in $meth_list
do
  for res in $res_list
  do
    if [ $res == "in" ]
    then
      res=""
    fi
    
    echo $meth$res
    
    if [ "$do_base" = true ]
    then
      echo "base"
      Rscript ./ProbReco/bboot_forecasts.R $model $trans $meth $res
    fi
    
    if [ "$do_base_score" = true ]
    then
      echo "score"
      Rscript ./ProbReco/bboot_scores.R $model $trans $meth$res
    fi
    
    if [ "$do_csbu" = true ]
    then
      echo "csbu"
      Rscript ./ProbReco/rboot_scores_ctbu_cs.R $model $trans $meth$res $nn
    fi
    
    if [ "$do_tebu" = true ]
    then
      echo "tebu"
      Rscript ./ProbReco/rboot_scores_ctbu_te.R $model $trans $meth$res $nn
    fi
    
    if [ "$do_oct" = true ]
    then
      echo "oct"
      Rscript ./ProbReco/rboot_scores_oct.R $model $trans $meth$res "h" $nn
    fi
    
    if [ "$do_octf" = true ]
    then
      echo "octf"
      Rscript ./ProbReco/rboot_scores_octf.R $model $trans $meth$res "in" $nn
    fi
    
  done
  echo $meth
done

exec bash
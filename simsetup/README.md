# Simulation Project

#### Lunch order:

1. simulation.R: create simulations and base forecasts
2. residuals:
  - hfitted_residuals.R: produce the multi-step residuals according to the base forecasts.
  - overlap_residuals.R: produce the overlapping residuals according to the base forecasts.
  - hoverlap_residuals.R: produce the multi-step overlapping residuals according to the base forecasts.
3. Bootstrap/gaussian forecasts and accuracy:
  - ctjb.sh: bash file to produce bootstrap base and reconciled forecasts 
  - ctsam.sh: bash file to produce gaussian (sample covariance matrix) base and reconciled forecasts
  - ctshr.sh: bash file to produce bootstrap (shrinkage covariance matrix) base and reconciled forecasts 
4. ProbScore/script_score.R: script to summarise the accuracy scores
5. Tables/tables_script.R: script to produce the paper tables


#### ProbReco folder
Script to generate base and reconciled forecasts (used by bash files)

Base forecasts and accuracy:

- bboot_forecasts.R: bootstrap base forecasts
- bgauss_forecasts.R: gaussian base forecasts
- bboot_scores.R: bootstrap/gaussian base accuracy

Bootstrap/gaussian reconciled forecasts and accuracy: 

- rboot_scores_cs.R: cross-sectional reconciliation
- rboot_scores_te.R: temporal reconciliation
- rboot_scores_ctbu_cs.R: cross-sectional reconciliation + temporal bu
- rboot_scores_ctbu_te.R: temporal reconciliation + cross-sectional bu
- rboot_scores_oct.R: optimal cross-temporal reconciliation

#### analisys folder
Script to evaluate the covariance matrix and make some plot presneted in the paper
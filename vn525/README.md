# Australian Tourism Demand (VN525)

The different R scripts have been collected into sh files. In order, we have
1.  `ctbase.sh`: point base forecasts;
2.  `ctgauss.sh`: probabilistic base and reconciled forecasts with the Gaussian approach;
3.  `ctjb.sh`: probabilistic base and reconciled forecasts with the ctjb approach;
4.  `Rscript ./ProbScore/scrip_score.R model transformation`: accuracy scores (see `ctbase.sh` for the description of `model` and `transformation`).

## Directories

-   `BaseForecasts`: point base forecasts scripts and results;
-   `HfittedRes`: h-step residuals scripts and results;
-   `ProbReco`: probabilistic base and reconciled forecasts scripts and results;
-   `ProbScore`: probabilistic base and reconciled accuracy scores scripts and results;
-   `Figures`: Figures paper;
-   `Tables`: Tables paper and online appendix.

## Data

-   `VN525.RData`: Australian Tourism Demand RData (228 monthly observations of Visitor Nights (VN) from January 1998 to December 2016, 525 time series)

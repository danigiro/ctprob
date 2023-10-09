# Australian GDP (AusGDP)

The different `R` scripts have been collected into `sh` files. In order, we have
0.  `Rscript ./data_preparation.R`: data cleaning;
1.  `ctbase.sh`: point base forecasts;
2.  `ctgauss.sh`: probabilistic base and reconciled forecasts with the Gaussian approach;
3.  `ctjb.sh`: probabilistic base and reconciled forecasts with the ctjb approach;
4.  `Rscript ./ProbScore/scrip_score.R model transformation`: accuracy scores (see `ctbase.sh` for the description of `model` and `transformation`).

## Directories

-   `BaseForecasts`: point base forecasts scripts and results;
-   `HfittedRes`: h-step residuals scripts and results;
-   `HOverlapRes`: h-step overlapping residuals scripts and results;
-   `OverlapRes`: 1-step overlapping residuals scripts and results;
-   `ProbReco`: probabilistic base and reconciled forecasts scripts and results;
-   `ProbScore`: probabilistic base and reconciled accuracy scores scripts and results;
-   `Figures`: Figures paper;
-   `Tables`: Tables paper and online appendix.

## Data

-   `DataRaw`: raw data for the Australian Quarterly National Accounts (QNA) dataset (1984:Q4 – 2018:Q1)

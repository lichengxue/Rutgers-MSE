# model_runs

_This document holds all information on production level model runs_

## 06/17/2026 model runs

Ran `parallelized_historic_run.R` for the temperature optima and trend configurations.

1. Temperature trend: $0 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1}$
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.05)$
3. Temperature optimum: $0 ^{\circ}\mathrm{C}$

Code is as follows (Lines 451 onward of `code/parallelized_historic_run.R`)
```{r}
ecov_proj_error <- rnorm(n=n,mean=0, sd=0.5) # Error around the rising trend
Ecov_re[,] <- 0.04*seq(1,n) + ecov_proj_error
```

Resulting in a trend of $0 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1} + \epsilon_{y}$.
The results are stored in `models/sensitivity-analysis/2026-06-17_09-50-46`.

**PLEASE NOTE THAT THIS MODEL RUN IS UNUSABLE BECAUSE OF INCORRECT TEMPERATURE PROJECTION 
CONFIGURATION**

## 06/18/2026 model runs

As you might know, we fixed a major error in temperature projections for the model estimation portion 
of the code. The fix is now implemented in `code/parallelized_historic_run.R`.

Ran `code/2026_06_18_copy_em_sensitivity_analysis.R` (a sequential version of the parallelized code) 
for the temperature optima and trend configurations.

1. Temperature trend: $0 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1}$
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.05)$
3. Temperature optimum: $0 ^{\circ}\mathrm{C}$

This seems to work now!
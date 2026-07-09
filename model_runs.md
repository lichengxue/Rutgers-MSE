# model_runs

_This document holds all information on production level model runs_

## 06/17/2026 model runs

Ran `parallelized_historic_run.R` for the temperature optima and trend configurations.

1. Temperature trend: $0 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1}$
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.5)$
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
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.5)$
3. Temperature optimum: $0 ^{\circ}\mathrm{C}$

This seems to work now!


## 06/24/2026 model runs

### Run 1 - 

Ran `code/parallelized_historic_run.R` with the following configurations.

1. Temperature trend: $0.04 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1}$
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.05)$
3. Temperature optimum: $-1.5 ^{\circ}\mathrm{C}$

Used following configurations

| $\sigma_{\text{NAA}}$ 	| $t_{\text{g}}$ 	| $w_{\text{opt}}$ 	|
|-----------------------	|----------------	|------------------	|
|  0.2                    | 6               | 2                 |

Model objects, logs, outputs, and visualizations saved to the following location: `models/sensitivity_analysis/2026-06-24_10-10-02/`

### Run 2 - 

Ran `code/parallelized_historic_run.R` with the following configurations.

1. Temperature trend: $0.04 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1}$
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.05)$
3. Temperature optimum: $-1.5 ^{\circ}\mathrm{C}$

Used following configurations

| $\sigma_{\text{NAA}}$ 	| $t_{\text{g}}$ 	| $w_{\text{opt}}$ 	|
|-----------------------	|----------------	|------------------	|
|  0.2                    | 6               | 1                 |



## 06/25/2026 model runs

### Run 1 

Ran `code/parallelized_historic_run.R` with the following configurations.

1. Temperature trend: $0.04 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1} + \epsilon_{y}$
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.05)$
3. Temperature optimum (Standardized): $0 ^{\circ}\mathrm{C}$

Used following configurations

| $\sigma_{\text{NAA}}$ 	| $t_{\text{g}}$ 	| $w_{\text{opt}}$ 	|
|-----------------------	|----------------	|------------------	|
|  0.2                    | 6               | 1,2               |

_Ran for 100 iterations using 6 parallel cores_

Model objects, logs, outputs, and visualizations saved to the 
following location: `models/sensitivity_analysis/2026-06-25_01-53-08`

*Note*: Interrupted the run around the 50th iteration of the first configuration.

### Run 2

Ran `code/parallelized_historic_run.R` with the following configurations.

1. Temperature trend: $0.04 ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1} + \epsilon_{y}$
2. Error term: $\epsilon_{y} \sim \mathcal{N}(0, 0.05)$
3. Temperature optimum (Standardized): $0 ^{\circ}\mathrm{C}$

Used following configurations

| $\sigma_{\text{NAA}}$ 	| $t_{\text{g}}$ 	| $w_{\text{opt}}$ 	|
|-----------------------	|----------------	|------------------	|
|  0.2                    | 6               | 2                 |

_Ran for 100 iterations using 6 parallel cores_

Model objects, logs, outputs, and visualizations saved to the 
following location: `models/sensitivity_analysis/2026-06-25_09-18-13`.
All 100 iterations were run and then visualized. Looks OK. 

## 07/08/2026 model runs [Still running]

Ran `code/parallelized_historic_run.R` with the following configurations.

| $\sigma_{\text{NAA}}$ 	| $t_{\text{g}}$ 	| $t$ ($^{\circ}\mathrm{C}\,\mathrm{yr}^{-1}$) 	| $\epsilon_{y}$ 	| $t_{\text{opt}}$ ($^{\circ}\mathrm{C}$) 	| $w_{\text{opt}}$ 	|
|-----------------------	|----------------	|----------------------------------------------	|----------------	|-----------------------------------------	|------------------	|
| 0.2                   	| 6              	| 0.0, 0.04, 0.1                               	| 0.75           	| -1.5, 0, 2.5                            	| 2                	|


_Ran for 120 iterations using 8 parallel cores_

Model objects, logs, outputs, and visualizations saved to the 
following location: `models/sensitivity_analysis/2026-07-08_17-45-13`.

There seems to be some convergence issues with the 0.0 trend. Investigate as to why that might be.


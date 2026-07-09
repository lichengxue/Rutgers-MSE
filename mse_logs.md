# MSE development logs

## 08/22/2025 log

Issues related to `bsb_mse_example_01.R` have been fixed. There was an error when 
running the MSE loop (`loop_through_fn`) resulting in non-convergence. This was due to 
not using all the years from 1989 to 2023 to inform the estimation model. 

## 08/25/2025 log

To store all model runs, I'm creating a Report output directory. The file structure 
will be as follows. This will allow for better logging of model settings and 
model error and warning logs for all MSE runs.

Ran a version of `bsb_mse_example_01.R` with `MSE_years` set to _20_. This gave rise 
to a number of warnings shown below.

>Warning messages:
>1: In fit_wham(em_input, do.retro = do.retro, do.osa = do.osa,  ... : 
>** Error during model fit. **
>Check for unidentifiable parameters.
>
>system is computationally singular: reciprocal condition number = 1.9425e-20
>
>2: In project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) :
>  Difference between projection model nll and base model nll is 1.24568918942396
>3: In fit_wham(em_input, do.retro = do.retro, do.osa = do.osa,  ... : 
>** Error during model fit. **
>Check for unidentifiable parameters.
>
>system is computationally singular: reciprocal condition number = 6.1186e-20
>
>4: In project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) :
>  Difference between projection model nll and base model nll is 0.534905132018594
>5: In sqrt(diag(cov)) : NaNs produced
>6: In sqrt(diag(object$cov.fixed)) : NaNs produced
>7: In project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) :
>  Difference between projection model nll and base model nll is -1.51824730055841
>8: In sqrt(diag(cov)) : NaNs produced
>9: In sqrt(diag(object$cov.fixed)) : NaNs produced
>10: In sqrt(diag(object$cov.fixed)) : NaNs produced
>11: In sqrt(as.numeric(object$diag.cov.random)) : NaNs produced
>12: In sqrt(diag(object$cov.fixed)) : NaNs produced
>13: In sqrt(as.numeric(object$diag.cov.random)) : NaNs produced
>14: In sqrt(diag(object$cov.fixed)) : NaNs produced
>15: In sqrt(as.numeric(object$diag.cov.random)) : NaNs produced
>16: In sqrt(diag(object$cov.fixed)) : NaNs produced
>17: In sqrt(as.numeric(object$diag.cov.random)) : NaNs produced
>18: In project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) :
>  Difference between projection model nll and base model nll is -0.0448178645774533
>19: In fit_wham(em_input, do.retro = do.retro, do.osa = do.osa,  ... : 
>** Error during model fit. **
>Check for unidentifiable parameters.
>
>system is computationally singular: reciprocal condition number = 1.57802e-16

## 08/27/2025 log

Created `bsb_mse_example_02.R` - This is for trying out a different harvest control rule on the stock

New HCR is a hockey stick that fishes at a maximum of 80% SPR and goes down to 0.01%

### Meeting with Cheng

- MSE is showing warnings and errors once it's run for 20 years
- How to save all warnings and errors from a run to a textfile?
- Why isn't the MSE working beyond 3 years at the end?

Potential fixes for non-convergence issues

- Changing the years-at-age for Year 1 in the estimation model to _equilibrium_
  - ```r
      NAA_re$N1_model[] = "equilibrium"
    ```
- Increasing `sigma_vals` for the operating model to 0.75 to reflect the real stock assessment
- Increase `prior_sigma` of the estimation model's movement 
- Change whether we want to estimate movement with the MSE function (`loop_through_fn`)
- Change the operating model's population dynamics model from a state-space model with random effects to 
just a state-space model. This needs to be reflected in the estimation model as well. 
    - ```r
    sigma <- "rec"
    re_cor <- "iid"
    ini.opt <- "age-specific-fe"
    sigma_vals <- array(0.2, dim = c(n_stocks, n_regions, n_ages)) # NAA survival sigma
    sigma_vals[, , 1] <- 0.75 # Recruitment sigma
    # For the estimation model
    NAA_re$sigma = "rec"
      ```

## 08/28/2025 log

Added functionality for plotting estimated SSB vs. true SSB for MSE (Can be found in 
`bsb_mse_example_01.R`)

Still to be done

1. Make a file and analysis output pipeline that doesn't overwrite the "Report" folder
2. Run multiple estimation models
3. Run multiple HCRs
4. Figure out what other important plots are

## 09/03/2025 log

Improved the functionality for plotting estimated SSB vs. true SSB for MSE (Can be
found in `bsb_mse_example_01.R`)

Things to be done next

1. Make a file and analysis output pipeline that doesn't overwrite the "Report" folder
2. Make a version of the model where the estimation model is configured exactly the same as the OM
3. Make a version of the operating model with basically no errors so that we can check it's behavior

> Why is this important? 
> Because it will teach me what each error introduced to the system will do. This will give me 
a better understanding of the OM.

## 09/11/2025 log

- Making improvements to `bsb_mse_example_01.R` to store run time information as well as generate
separate folders for Reports.
- Made a script named `bsb_mse_example_03.R` that will be specifying the estimation model the same as 
the operating model
- Making the following changes to the estimation model and looking at different outcomes
  - `move_em$prior_sigma` set to 1 (from 0.2) - Not greater difference other than 
  
## 09/17/2025 log

- Held weekly meeting and set goals and expectations on deliverable for the MSE

## 09/19/2025 log

- Completed the notes for `bsb_mse_example_01.R` @ `mse_notebook.md`

## 10/06/2025 log

- Uninstalling `whamMSE` and `wham`

## 10/07/2025 log

- Meeting with Cheng
- Use this for plotting wham fits: `wham::plot_wham_output()`


## 10/22/2025 log

- Uninstalled and reinstalled `wham` and `whamMSE` packages
  - Used the following commands
    - `remove.packages("wham",whamMSE")`
    - `remotes::install_github("timjmiller/wham@lab")`
    - `remotes::install_github("lichengxue/whamMSE@Projection-MSE")`
  - *IMPORTANT* - I skipped installing any of the suggested package updates 
  - After this update, I can't run `bsb_mse_example_01.R` anymore
    - Getting an error about how proj.wham doesn't have enough years

> Error in project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) : 
>  
> ** Error setting up projections: **
> proj.opts$avg.yrs is not a subset of model years.


- Getting an error with bsb_mse_env_example_01.R

> Retro Peel: 7
> 
> --Ecov--------------------------------------------------------------------------------------------------------------------------------
> one or more ecov does not start by model year 1 - max(lag). Padding ecov... 
> Please check that the environmental covariates have been loaded and interpreted correctly.
> 
>      Model years: 1989 to 2014
>      Ecov years: 1988 to 2014
>
>    -------------------------------------------------------------------------------------------------------------------------------------
> 
> 
> Error in TMB::MakeADFun(temp$data, temp$par, DLL = "wham", random = temp$random,  : 
>   A map factor length must equal parameter length


## 11/04/2025 log

- Running into issues with `bsb_om_em_mse.R`
- The problem is with the final chunk of code `loop_through_fn`

> Error in project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) : 
>
> ** Error setting up projections: **
> proj.opts$avg.yrs is not a subset of model years.

- This issue now appears across all the code that uses `loop_through_fn`
- Might be related to an updated in the WHAM package (https://timjmiller.github.io/wham/reference/project_wham.html)
- Good example code for comparing model performance - https://github.com/lichengxue/whamMSE/blob/UMassD-MSE/GBK_Example/Update_MSE_code.R#L474

- *IMPORTANT: Uninstalled the whamMSE@Projection-MSE branch and reinstalled the regular WHAM*

`{r}
remotes::install_github("lichengxue/whamMSE",dependencies = FALSE)
`

## 11/05/2025 log

- Cheng reported that the bug referred to on 11/04/2025 has been fixed
- Attempting to reinstall `whamMSE@Projection-MSE` - *Successful*
- Ran `bsb_mse_example_01.R` after reinstallation - *Successful*
- Ran `bsb_om_em_mse.R` after reinstallation - *Successful*

## 11/11/2025 log

- Working through `bsb_om.R` and `bsb_om_em_mse.R` and making notes
- I think I understand how environmental covariates and linked and modeled in whamMSE

## 11/12/2025 log

- Trying to understand the four surveys in `bsb_om.R`.
- Put some thought into how to structure code for simulating the operating model sevearl hundred times

## 11/19/2025 log

- Setting up Annotate2 server (accessed on annotate2.sebs.rutgers.edu) to run simulations

## 12/10/2025 log

- Uninstalling `whamMSE` (Projection-MSE branch) and `wham`.
- Installing `wham` from Cheng's development branch for Guassian temperature linkage (`remotes::install_github("lichengxue/wham@Gussian_Rec")`)
- Used the following commands
  - `remove.packages("wham","whamMSE")`
  - `remotes::install_github("lichengxue/wham@Gussian_Rec")`
  - `remotes::install_github("lichengxue/whamMSE@Projection-MSE")`
- Renamed `Test_code.R` to `bsb_gaussian_test_code.R` - This runs the gaussian link
- Running into issues with `bsb_gaussian_test_code.R` after this update

> > temp <- wham::prepare_wham_input(asap)
>
> --Creating input---------------------------------------------------------------------------------------------------------------------
> 
> Error in set_basic_info(input, basic_info) : 
>   'list' object cannot be coerced to type 'integer'

## 12/16/2026 log

- Meeting with John and Cheng
- Reduce the spatial complexity of the model from two regions to one region
- **IMPORTANT**: Now we are going for a single region and single stock and testing for
the guassian function.

## 01/13/2026 log

- Previous issues with the 'Gaussian branch' appears to be resolved
- Uninstalled and reinstalled the `wham` and `whamMSE` branches
- Used the following commands
  - `remove.packages("wham")`
  - `remove.packages("whamMSE")`
  - `remotes::install_github("lichengxue/wham@Gussian_Rec")`
  - `remotes::install_github("lichengxue/whamMSE@Projection-MSE")`
- Ran the original code for `bsb_om_em_mse.R`
- Gave the follow error/warning: _The code still ran_

> Warning message:
> In fit_wham(em_input, do.retro = do.retro, do.osa = do.osa, do.brps = TRUE,  :
>   
> ** Error during model fit. **
> Check for unidentifiable parameters.
> 
> Lapack routine dgesv: system is exactly singular: U[184,184] = 0

- Output is saved in `Rutgers-MSE/models/mod_2026_01_13.RDS`

- **IMPORTANT** - Wrote a single region version of the base model for future development

**IMPORTANT** - <span style="color:red"> *ANY OPERATING MODEL THAT GIVES THIS ERROR SHOULD BE DISCARDED*</span>


## 01/29/2026 log

- Uninstalled and reinstalled `wham` and `whamMSE`
- Pulled the latest `Gussian_Code_Final.R` code from Cheng's commit on the master branch
- Renamed `Gussian_Code_Final.R` to 

## 02/10/2026 log

- Meeting with Cheng on 02/10/2026
- Decided to do Sensitivity analysis on the MSE (Gaussian form)
- Going to vary the following parameters
  1. Sigma for recruitment and NAA - `vals[]` - 0.2 [Realistic for NE stocks], 0.5, 1 [Very high]
  2. MSE gap - every 3 years, every 6 years
  3. Width of gaussian relationship - `input_Ecov$par$log_width_rec` - 0.1 [Very stringent], 1, 5 [Very forgiving to non optimal temperatures]
  4. Variability of temperature - For now, center around zero and increase stochasticity - 0.5, 1, 3

- Double check the black sea bass report for the magnitude of the recruitment, catch, and SSB
- This is in order to match up our estimates in scale to reality


## 02/16/2026 log

- Cleaning up the .git project folder.
- Moving all .R files into `scratch/` folder
- All temporary development will be done there from now on
- Current development is been done on the `Gussian_Code_Final_15yearMSE.R`. _Need to be renamed_
- Sensitivity analysis in progress - Current working file [`code/em_sensitivity_analysis.R`](code/em_sensitivity_analysis.R)

## 02/17/2026 log

- Done cleaning up the .git project folder
- Now there are three new folders including 
  - [`scratch/`](scratch/)(for developmental/tutorial/past .R files that are currently tracked)
  - [`code/`](code/)(for currently developing .R files)
  - [`git_untracked/`] - Only present on the local machine
- Moved `FAA.RDS` and `WAA.RDS` to the [`data/`](data/) folder

## 02/18/2026 log

- Development switched to the `dev` branch
- Changing [`code/em_sensitivity_analysis.R`](code/sensitivity_analysis.R) to be more flexible in the future
- Revert back to the following commit if something bad happens or break: [d66ff06](https://github.com/lichengxue/Rutgers-MSE/commit/d66ff0615ae8b0df63daf3e93f447c2cd3ed4067)

## 03/05/2026 log

- Development undergoing in `dev` branch
- Writing methodology 

## 03/18/2026 log

- Uninstalling `whamMSE` (Projection-MSE branch)
- Used the following commands
  - `remove.packages("whamMSE")`
  - `remotes::install_github("lichengxue/whamMSE@Projection-MSE")`
- Then incorporated new code for an estimation model that has a fixed temperature-recruitment relationship.
- The model cannot _estimate_ the parameters related to the gaussian link by itself ($T_{opt}$, $w_{\text{opt}}$)
- So as of right now, we estimate it and the model fits everything else around it. We can think of this as an *exceptionally strong prior* (_I guess...._)

- Ran `code/em_sensitivity_analysis.R` for only the first 6 configurations of the sensitivity settings to 
check if the package was broken.

## 03/23/2026 log

- Changing some parts of `code/em_sensitivity_analysis.R` for reproducibility and completeness. Listing them below.

*Major change is that we are changing the MSE years from 15 to 18*

1. Line 94 - `{r} n_feedback_years` changed from 15 to 18
2. Line 116 - `ecov$year <- c(north_bt[,"year"], 2023:(2025+12))` changed to `ecov$year <- c(north_bt[,"year"], 2023:(2025+15))`
3. Line 136 - `r MSE_years` changed from 15 to 18 (_Maybe these two variables can be the same thing. Check back later_)
4. Line 163 - `r for (i in 34:(36+12)) {` changed to `r for (i in 34:(36+15)) {`
4. Line 310 - `r index_Neff <- rbind(index_Neff, index_Neff[rep(33,15), , drop = FALSE])` changed to `r index_Neff <- rbind(index_Neff, index_Neff[rep(33,MSE_years), , drop = FALSE])`
5. Line 312 - `r for (i in 34:(36+12)) {` changed to `r for (i in 34:(36+15)) {`
5. Line 330 - `r catch_Neff <- rbind(catch_Neff, catch_Neff[rep(33,15), , drop = FALSE])` changed to `r catch_Neff <- rbind(catch_Neff, catch_Neff[rep(33,MSE_years), , drop = FALSE])`
5. Line 335 - `r for (i in 34:(36+12)) {` changed to `for (i in 34:(36+15)) {`
6. Line 466 - `r last.year <- 2024+12` changed to `r last.year <- 2024+15`

## 03/24/2026 log

*IMPORTANT*: Run time for just the first 4 configurations was 2 hours and 36 minutes and 8 seconds.
This might need a server to run properly.


## 03/26/2026 log

- New `wham` branch is up!
- Use the following updated instructions to uninstall and reinstall `wham` and `whamMSE` from now on.
  - `remove.packages("wham")`
  - `remove.packages("whamMSE")`
  - `remotes::install_github("lichengxue/wham@Gaussian_test")`
  - `remotes::install_github("lichengxue/whamMSE@Projection-MSE")`

## 03/31/2026 log

Meeting with Cheng. Some questions I have for him. Meeting was recorded on Zoom for references.
Check Documents/Zoom folder for the recording.

1. Remove lines 377-379 from the code (Setting random effects for the Ecov)
2. Need to make some modifications to Line 389, 390, because we haven't defined the future
3. Insert a new line after Line 390 -> input_Ecov$par$Ecov_rec[65:XX] <- Ecov_re

### Increasing flexibility of assessment/future projection years
Editing `code/em_sensitivity_analysis.R`
_Made a commit before the major changes - Hash is bdf64d1..02f7fee in the `dev` branch_

1. Line 49 - New line. Inserted `n_feedback_years <- 30`. Everything stems from here.
2. Line 125 - Removed `n_feedback_years <- 30`. It's now moved to Line 49. 
3. Line 149 (Previously 116) - Previously `ecov$year <- c(north_bt[,"year"], 2023:(2025+15))`
Now, I've changed it to the following code chunk
```{r}
      ecov_final_year <- max(north_bt[,"year"])
      projection_years <- seq(ecov_final_year+1, ecov_final_year+n_feedback_years)
      ecov$year <- c(north_bt[,"year"], projection_years)
```
Basically I've taken the last year from the temperature time series and made a vector for all the years 
that we are projecting to and updated `ecov$year` using that.
4. Line 171 - Changed `MSE_years  <- 18` to `MSE_years  <- n_feedback_years`
5. Line 180 - Changed `user_maturity[, i, ] <- OMa$input$data$mature[1, 33, , drop = FALSE]` to `user_maturity[, i, ] <- OMa$input$data$mature[1, hist_years, , drop = FALSE]`
6. Line 185 - Changed `user_waa$waa <- array(NA, dim = c(5, 33 + n_feedback_years, 8))` to `user_waa$waa <- array(NA, dim = c(5, hist_years + n_feedback_years, 8))`
7. Line 186 - changed `user_waa$waa[, 1:33, ] <- OMa$input$data$waa[c(1,2,5,6,9), ,]` to `user_waa$waa[, 1:hist_years, ] <- OMa$input$data$waa[c(1,2,5,6,9), ,]`
8. Lines 187, 188
```{r}
      for (i in 34:(36+15)) {
        user_waa$waa[, i, ] <- OMa$input$data$waa[c(1,2,5,6,9), 33, ]
      }
```
**Changed to **
```{r}
      for (i in om_future_start_index:om_future_end_index) {
        user_waa$waa[, i, ] <- OMa$input$data$waa[c(1,2,5,6,9), hist_years, ]
      }
```
9. Line 242: Changed `F_info$F[1:33,] <- OMa$rep$Fbar[, 1:2]` to `F_info$F[1:hist_years,] <- OMa$rep$Fbar[, 1:2]`

**Stopped at Line 332 for today**
10. Lines 332 to 375
```{r}

      # index sigma and Neff
      input_Ecov$data$agg_index_sigma[1:33,] <- OMa$input$data$agg_index_sigma[,1:2]
      input_Ecov$data$use_indices[1:33,]     <- OMa$input$data$use_indices[,1:2]
      input_Ecov$data$use_index_paa[1:33,]   <- OMa$input$data$use_index_paa[,1:2]

      for (i in 34:(36+15)) {
        input_Ecov$data$agg_index_sigma[i,] <- OMa$input$data$agg_index_sigma[33,1:2, drop = FALSE]
      }

      idx1 <- which(asap[[1]]$dat$use_index == 1)

      Neff1 <- do.call(cbind, lapply(idx1, function(i)
        asap[[1]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
      index_Neff <- Neff1
      index_Neff <- rbind(index_Neff, index_Neff[rep(33,MSE_years), , drop = FALSE])
      input_Ecov$data$index_Neff <- index_Neff

      input_Ecov <- whamMSE::update_input_index_info(
        input_Ecov,
        agg_index_sigma = input_Ecov$data$agg_index_sigma,
        index_Neff      = input_Ecov$data$index_Neff
      )

      # catch sigma & Neff
      input_Ecov$data$agg_catch_sigma[1:33,] <- OMa$input$data$agg_catch_sigma[,1:2]
      input_Ecov$data$use_agg_catch[1:33,]   <- OMa$input$data$use_agg_catch[,1:2]
      input_Ecov$data$use_catch_paa[1:33,]   <- OMa$input$data$use_catch_paa[,1:2]

      for (i in 34:(36+15)) {
        input_Ecov$data$agg_catch_sigma[i,] <- OMa$input$data$agg_catch_sigma[33,1:2]
      }

      Neff1 <- asap[[1]]$dat$catch_Neff
      catch_Neff <- cbind(Neff1)
      catch_Neff <- rbind(catch_Neff, catch_Neff[rep(33,MSE_years), , drop = FALSE])

      input_Ecov$data$catch_Neff <- catch_Neff

      input_Ecov <- update_input_catch_info(
        input_Ecov,
        agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
        catch_Neff      = input_Ecov$data$catch_Neff
      )

```
Changed to 
```{r}

      # index sigma and Neff
      input_Ecov$data$agg_index_sigma[1:hist_years,] <- OMa$input$data$agg_index_sigma[,1:2]
      input_Ecov$data$use_indices[1:hist_years,]     <- OMa$input$data$use_indices[,1:2]
      input_Ecov$data$use_index_paa[1:hist_years,]   <- OMa$input$data$use_index_paa[,1:2]

      for (i in om_future_start_index:om_future_end_index) {
        input_Ecov$data$agg_index_sigma[i,] <- OMa$input$data$agg_index_sigma[hist_years,1:2, drop = FALSE]
      }

      idx1 <- which(asap[[1]]$dat$use_index == 1)

      Neff1 <- do.call(cbind, lapply(idx1, function(i)
        asap[[1]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
      index_Neff <- Neff1
      index_Neff <- rbind(index_Neff, index_Neff[rep(hist_years,MSE_years), , drop = FALSE])
      input_Ecov$data$index_Neff <- index_Neff

      input_Ecov <- whamMSE::update_input_index_info(
        input_Ecov,
        agg_index_sigma = input_Ecov$data$agg_index_sigma,
        index_Neff      = input_Ecov$data$index_Neff
      )

      # catch sigma & Neff
      input_Ecov$data$agg_catch_sigma[1:hist_years,] <- OMa$input$data$agg_catch_sigma[,1:2]
      input_Ecov$data$use_agg_catch[1:hist_years,]   <- OMa$input$data$use_agg_catch[,1:2]
      input_Ecov$data$use_catch_paa[1:hist_years,]   <- OMa$input$data$use_catch_paa[,1:2]

      for (i in om_future_start_index:om_future_end_index) {
        input_Ecov$data$agg_catch_sigma[i,] <- OMa$input$data$agg_catch_sigma[hist_years,1:2]
      }

      Neff1 <- asap[[1]]$dat$catch_Neff
      catch_Neff <- cbind(Neff1)
      catch_Neff <- rbind(catch_Neff, catch_Neff[rep(hist_years,MSE_years), , drop = FALSE])

      input_Ecov$data$catch_Neff <- catch_Neff

      input_Ecov <- update_input_catch_info(
        input_Ecov,
        agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
        catch_Neff      = input_Ecov$data$catch_Neff
      )

```
11. Line 446 changed from `om_ecov$parList$F_pars[34:(36+15),] = om_with_data$rep$log_SPR_FXSPR_static` to 
`om_ecov$parList$F_pars[om_future_start_index:om_future_end_index,] = om_with_data$rep$log_SPR_FXSPR_static`
12. Line 502 changed from `last.year <- 2024+15` to `last.year <- terminal.year+n_feedback_years`. **Check this again. Not sure about it!!!**


**IMPORTANT: Might need to move the tryCatch from the larger loop to the seed loop**


## 04/06/2026 log

We also need to make a single model + multiple iteration seed run code available here

Change made to `code/em_sensitivity_analysis_30_years.R`. Changed mean recruitment in operating model from exp(10.5) to exp(12)

John also wanted recruitment deviates - mod1$em_full[[1]]$rep$NAA_devs
For the operating model - mod1$om$rep$NAA_devs_1

## 04/13/2026 log

Initiating a model run with the following parameters. Saved to `2026-04-13_02-13-32`

1. iterations = 10
2. years =15
3. $\sigma_{NAA}$ = 0.2
4. $t_g$ = {3 years, 6 years}
5. $w_\text{opt}$ = 0.5

*"Total execution time was 5 hours and 19 minutes and 42 seconds"*


## 04/14/2026 log

Regarding issues with convergence when setting NAA random effects process error to very low values.
What happens with the convergence issues for Model 3 even though the parameters are fixed?
The estimator is having trouble getting to such a low boundary with the other parameters.

Link to some example code for parallelization and visualizing multiple models - 
https://lichengxue.github.io/SPASAM.MSE/Performance-Analysis-Tools.html


## 04/21/2026 log

To do list for 04/21/2026

1. Set up some temperature scenarios for 

Working on the following files - `scratch/parallel_computations.R`, `scratch/environmental_projections`

## 04/28/2026 log

1. Make the plots for the true and estimated values between the models

## 05/16/2026 log

Installing `SPASAM.MSE` - devtools::install_github(
  "lichengxue/SPASAM.MSE",
  dependencies = FALSE
)
Added 'do.brps=TRUE', to code/parallelized_code_sensitivity_run.R in the `loop_through_fn` function

NOTE: Cheng said to include code to make sure that random seeds are not generated

Need to talk to John about the following
1. Duration of the length of the feedback period (longer is better: 30-50 years)
2. Historical temperature trends and future temperature scenarios

## 05/19/2026 to current log

Ran the model for an estimation period of 15 years at 20 iterations - Tooks 10 minutes 41 seconds
Ran the model for an estimation period of 50 years at 8 iterations (2026-05-21_09-20-25) - Took 39 minutes and 0 seconds

Made a new file to introduce temperature trends and set new hypothetical temperature optima

Set the following
1. `input_Ecov$par$Topt_rec` = 1.5, Width was set at 1
2. Removed the cosine varying random effects that we set in initial models
3. 


## 05/22/2026 log

- Use the following updated instructions to uninstall and reinstall `wham`, `whamMSE`, and `SPASAM.MSE`
from now on.
  - `remove.packages("wham")`
  - `remove.packages("whamMSE")`
  - `remove.packages("SPASAM.MSE")`
  - `remotes::install_github("lichengxue/wham@Gaussian_test")`
  - `remotes::install_github("lichengxue/whamMSE@Projection-MSE")`
  - `devtools::install_github("lichengxue/SPASAM.MSE", dependencies = FALSE)`
  
  
## 05/26/2026 log

*IMPORTANT* - Caught something important in the code. Estimation model doesn't set the 
PercentFXSPR to 75
Lines 516 of code/parallelized_historic_run.R and similar in code/parallelized_code_sensitivity_run.R is a new addition
`R om_ecov$input$percentFXSPR <- 75 # This is new [2026-05-26]`

## 06/11/2026 log

Working on Stephen's problem. Ran different configurations of this model varying the following parameters.  
The problem is that there isn't enough variation between recruitment (top panel in the plots).

1. $\beta_R$ - `Ecov_beta_R`
2. $\text{log}NAA_{\sigma}$ - `log_NAA_sigma`
3. $\sigma$ - `sigma_vals`

An output plot resulting from different configurations is shown below.

![alt text](https://github.com/lichengxue/Rutgers-MSE/blob/dev/plots/potts-analysis/potts_outputs.png)


## 06/18/2026 log

Found a major error in the projection code where the projected temperatures are not applied and 
future temperatures are set to zero. This has been fixed but not implemented in 
`code/parallelized_historic_run.R` yet.

To be done (Immediately):

1. Need to incorporate metadata on `iterations` to the log files [x]
2. Need to calculate convergence rates at a global scale as well as for each model type and configuration []

To be done (Down the road):

1. Ability to merge multiple model runs while accounting for duplicate seeds so that a total of at least 
100 seed outputs can be prepared for final analysis and manuscript submission.

## 06/23/2026 log

Reinstalled `SPASAM.MSE` package to check whether the `plot_mse_output()` function works now.
*_It works!_*

Moved the following files from `code/` to `scratch/`

- `2026_04_07_copy_em_sensitivity_analysis_30_years.R`
- `2026_04_21_copy_em_sensitivity_analysis_30_years.R`
- `2026_06_18_copy_em_sensitivity_analysis_30_years.R`
- `em_sensitivity_analysis_30_years.R`
- `em_sensitivity_analysis.R`


## 06/24/2026 log

Ran several models. Check `model_runs.md` for detailed notes on those.

*IMPORTANT*: Fixed something in the code at `code/parallelized_historic_runs.R`. Future temperature projections initiated 
at 0 (historical mean) instead of the last recorded temperature. This is fixed now. The image below shows the result of the fix 
in detail.

<p align="center">
  <img width="60%" src="https://github.com/lichengxue/Rutgers-MSE/blob/dev/images/ecov_re_fix_for_github_notes.png" />
  <br>
  Figure 1 - Fix implemented on how future temperature projections are incorporated into the timeseries
</p>

Line 472 in commit -  https://github.com/lichengxue/Rutgers-MSE/commit/02f8db450e54ecb33a299760573dea5da0620855

## 06/25/2026 log

Ran several models. Check `model_runs.md` for detailed notes on those.
Discussed with John about simplifying some of the temperature optima and the trends.


## 07/04/2026 log [Server log]

** WRITING THIS FROM THE SERVER SIDE **

Setting up the server to run the code. Use `renv` package to maintain packages 
and dependencies in a local/contained folder.

Installed `TMB` as a dependency.

Getting an error running `code/parallelized_historic_run.R`. Similar to the one I'm getting on the Windows Annotate2 Server.

Seems to be an issue related to `TMB` package.

```{r}
Error: package or namespace load failed for ‘TMB’:
 .onLoad failed in loadNamespace() for 'TMB', details:
  call: dyn.load(file, DLLpath = DLLpath, ...)
  error: unable to load shared object '/Users/jeewanthabandara/Library/Caches/org.R-project.R/R/renv/cache/v5/macos/R-4.6/aarch64-apple-darwin23/TMB/1.9.21/62a5714b9af765a4c286f4e447d44cea/TMB/libs/TMB.so':
  dlopen(/Users/jeewanthabandara/Library/Caches/org.R-project.R/R/renv/cache/v5/macos/R-4.6/aarch64-apple-darwin23/TMB/1.9.21/62a5714b9af765a4c286f4e447d44cea/TMB/libs/TMB.so, 0x0006): symbol not found in flat namespace '_omp_get_max_threads'
```

## 07/07/2026 log

The issue regarding running the code on a server seems to be related to how `renv` 
package masks several C++ libraries from the `TMB` package. The issue is that 
`TMB` gets installed, but it's corrupted. I initially thought that this was related 
to an issue with `gfortran` library. I'm listing all the steps/fixes I tried below.

1. Set up an `renv` environment inside a fresh project on the server.
2. Installed `TMB`, `wham`, `whamMSE`, and `SPASAM.MSE` into this environment.
3. Ran `code/parallelized_historic_run.R`. Gave the following error regarding `OpenMP` linking.
```{r}
Error: package or namespace load failed for ‘TMB’:
 .onLoad failed in loadNamespace() for 'TMB', details:
  call: dyn.load(file, DLLpath = DLLpath, ...)
  error: unable to load shared object '/Users/jeewanthabandara/Library/Caches/org.R-project.R/R/renv/cache/v5/macos/R-4.6/aarch64-apple-darwin23/TMB/1.9.21/62a5714b9af765a4c286f4e447d44cea/TMB/libs/TMB.so':
  dlopen(/Users/jeewanthabandara/Library/Caches/org.R-project.R/R/renv/cache/v5/macos/R-4.6/aarch64-apple-darwin23/TMB/1.9.21/62a5714b9af765a4c286f4e447d44cea/TMB/libs/TMB.so, 0x0006): symbol not found in flat namespace '_omp_get_max_threads'
```
4. 

**NOTE**: The quick fix is simple. Simply installed `TMB`, `wham`, `whamMSE`, and `SPASAM.MSE` 
outside of `renv` (Directly to base R). This removes all the advantages of `renv`, but it works.


## 07/08/2026 log

Future scenarios for temperature projections has been updated. The following changes are now in effect.

1. Temperature trends ($t$) are now as follows: $\[0.04,0.102,0.0\] ^{\circ}\mathrm{C}\,\mathrm{yr}^{-1} + \epsilon_{y}$
2. Stochasticity (as determined by error term $\epsilon_y$) of all temperature trends are now same as the stochasticity in the last 10 years of the historical trend
3. Number of hypothetical temperature optima has been reduced. Now as follows: $\[0,-1.5, 2.5\] ^{\circ}\mathrm{C}$

The following table shows how the scenarios were simplified.


| Before 	| After 	|
|--------	|-------	|
| <img width="350" src="https://github.com/lichengxue/Rutgers-MSE/blob/dev/images/temperature_trends_and_optima_gathered.png" />       	| <img width="350" src="https://github.com/lichengxue/Rutgers-MSE/blob/dev/images/temperature_trends_and_optima_gathered_3.png" />             	|




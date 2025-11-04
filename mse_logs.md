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
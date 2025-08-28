# MSE notebook

*For taking any and all notes regarding MSE for BSB. Posts arranged from latest to oldest*

- Notes by and (mostly for the benefit) of RMWJ Bandara

## `bsb_mse_example_01.R`

This is a brief explanation of the example MSE used for black sea bass.
The model sets up a 3 year long MSE with assessments every 3 years. This is 
preceded by 35 years of historical data (1989 to 2023) after which the MSE 
takes over. The population is structured as two meta-stocks that have movement 
between each other. 

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
````

## 08/28/2025 log

Adding functionality for plotting estimated SSB vs. true SSB for MSE (Can be found in 
`bsb_mse_example_01.R`)


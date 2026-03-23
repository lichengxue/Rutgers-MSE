# -----------------------------------------------------------------------------
#
#                 Management Strategy Evaluation (MSE)
#           Black Sea Bass (BSB) with Environmental Drivers
#                       SENSITIVITY ANALYSIS
#
#
# Author(s): RMWJ Bandara, Chengxue Li
# Date: 2026/03/02
# Runtime environment: MacOS Sequoia 15.5 on M-chip Macbook Pro, R version 4.4.1
#
#
#
# -----------------------------------------------------------------------------

library(wham)
library(whamMSE)
library(dplyr)
library(here)
library(ggplot2)
library(gridExtra)
library(beepr)
library(readr)

here::here()

# NOTE: This code is set to run on a server as well as a personal computer (MacOS for now)

#### RUN TIMES ####
# Benchmarking run time and info for save files
run_start_time <- Sys.time()

#### RUN ENVIRONMENT ####
# Options for run environment
run_env_opts <- c("local","annotate2","amarel")
run_env <- run_env_opts[1]

#### MODEL RUN SETTINGS ####
# Set iterations, a base random seed, and then generate seeds for each MSE run
# Set a model name
# NOTE: NOT USING THESE SETTINGS FOR THE SENSITIVITY ANALYSIS
iterations <- 1
base_random_seed <- 853
set.seed(base_random_seed)
mse_random_seeds <- as.integer(floor(runif(iterations, min=0, max=1000)))
model_name <- "BSB Ecov"

#### SENSITIVITY ANALYSIS SETTINGS ####
# We will be looking at model variation from the following
# Process error for recruitment and NAA - vals - 0.2, 0.5, 1
# MSE gaps - 3 years, 6 years
# Width of gaussian relationship - `input_Ecov$par$log_width_rec` - 0.1 [Very stringent], 1, 5 [Very forgiving to non optimal temperatures]
# Variability of temperature - Mean centered around zero -  increasing stochasticity - 0.5, 1, 3
# NOTE: Variability of temperature is removed for now. Integration possible in the future

proc_error_v <- c(0.2, 0.5, 1) # Process error for NAA random effects
mse_gaps_v <- c(3,6) # Time between assessments for MSE
gauss_width_v <- c(0.1, 1, 5) # Width of the gaussian relationship (Wider = Less sensitive to optimal temperature)

# Use `crossing` function from tidyr to create a dataframe of all the settings in the
# sensitivity analysis
sens_analysis_settings <- tidyr::crossing(proc_error= proc_error_v, 
                                          mse_gap = mse_gaps_v, 
                                          gauss_width=gauss_width_v)

sens_analysis_settings <- sens_analysis_settings %>% mutate(nid=row_number(), .before=1)

#### IF THIS IS A TEST RUN - YOU ONLY WANT TO RUN A COUPLE OF ROWS!!! ####
# NOTE: ONLY USE FIVE ROWS BECAUSE THIS IS A TEST RUN
sens_analysis_settings <- sens_analysis_settings %>% head(n=1)

#### EXTERNAL FUNCTIONS ####
# Source reusable functions from `functions/reusable_functions.R`
source(here("functions","reusable_functions.R"))

#### DO THE SENSITIVITY ANALYSIS OR NOT? ####
# WARNING: THIS CAN TAKE A LOT OF TIME
DO_ANALYSIS = TRUE
if(DO_ANALYSIS){
  for(i in 1:nrow(sens_analysis_settings)) {
    # Envelope in a tryCatch so that we just skip over stuff that doesn't work
    tryCatch({
      skip_to_next <- FALSE # A flag we use to deal with convergence errors
      # Get the settings for this sensitivity run
      row <- sens_analysis_settings[i,]
      proc_error <- row$proc_error
      mse_gap <- row$mse_gap
      gauss_width <- row$gauss_width
      run_id <- row$nid
      print(paste("Running sensitivity analysis Run:",run_id, sep=" "))
      
      # NOTE: OM model ends in 2021 and EM starts in 2022
      # But the bottom temperature dataset goes until 2022
      n_feedback_years <- 15
      
      OMa  <- readRDS("models/OM_base.RDS")
      asap <- read_asap3_dat("data/north.dat")
      
      #### 1. ECOV DATA ####
      north_bt <- read.csv("data/bsb_bt_temp_nmab_1959-2022.csv")
      
      ecov <- list(label = c("North_BT"))
      
      # original means
      ecov$mean <- cbind(north_bt[,"mean"])
      
      # center each column
      ecov$mean <- t(t(ecov$mean) - apply(ecov$mean, 2, mean))
      
      # extend to include 3 feedback years (2023–2025)
      ecov$mean <- rbind(ecov$mean, matrix(0, n_feedback_years, ncol(ecov$mean)))
      
      # for this test, force north BT deviations to 0 so all variation comes from Ecov_re
      # ecov$mean[,1] <- 0
      
      ecov$year <- c(north_bt[,"year"], 2023:(2025+12))
      
      plot(ecov$year, ecov$mean[,1], type = "l", main = "Ecov mean (North BT)", xlab = "Year")
      
      # obs error for Ecov
      ecov$logsigma           <- "est_1"
      ecov$use_obs            <- matrix(1, nrow(ecov$mean), ncol(ecov$mean))
      ecov$process_model      <- "ar1"
      ecov$process_mean_vals  <- apply(ecov$mean, 2, mean)
      
      # keep linear Ecov–R links defined, but we will zero out Ecov_beta_R later
      ecov$recruitment_how <- matrix("controlling-lag-0-linear")
      
      #### 2. BIOLOGICAL/MSE FRAME ####
      n_stocks  <- 1
      n_regions <- 1
      n_ages    <- 8
      
      year_start <- 1989
      year_end   <- 2021
      MSE_years  <- 15
      
      hist_years <- length(year_start:year_end)
      total_years <- length(year_start:(year_end + MSE_years))
      
      # maturity
      user_maturity <- array(NA, dim = c(n_stocks, total_years, n_ages))
      user_maturity[, 1:hist_years, ] <- OMa$input$data$mature[1,,]
      for (i in (hist_years+1):(hist_years+MSE_years)) {
        user_maturity[, i, ] <- OMa$input$data$mature[1, 33, , drop = FALSE]
      }
      
      
      user_waa <- list()
      user_waa$waa <- array(NA, dim = c(5, 33 + n_feedback_years, 8))
      user_waa$waa[, 1:33, ] <- OMa$input$data$waa[c(1,2,5,6,9), ,]
      for (i in 34:(36+12)) {
        user_waa$waa[, i, ] <- OMa$input$data$waa[c(1,2,5,6,9), 33, ]
      }
      user_waa$waa_pointer_fleets   <- 1:2
      user_waa$waa_pointer_indices  <- 3:4
      user_waa$waa_pointer_totcatch <- 5
      user_waa$waa_pointer_ssb      <- 5
      user_waa$waa_pointer_M        <- 5
      
      fracyr_spawn <- 0.5
      
      #### 3.CATCH INFO ####
      catch_info <- list(
        catch_cv       = c(asap[[1]]$dat$catch_cv[1,]),
        catch_Neff     = c(asap[[1]]$dat$catch_Neff[1,]),
        use_agg_catch  = 1,
        use_catch_paa  = 1
      )
      
      fracyr_indices <- c(0.4166667, 0.2500000, 0.4166667, 0.2500000)
      
      index_info <- list(
        index_cv        = rep(0.5, 2),
        index_Neff      = rep(25, 2),
        fracyr_indices  = fracyr_indices[1:2],
        q               = OMa$rep$q[1,1:2],
        use_indices     = rep(1, 2),
        use_index_paa   = rep(1, 2),
        units_indices   = rep(2, 2),
        units_index_paa = rep(2, 2)
      )
      
      #### 4. GENERATE BASIC INFO ####
      info <- whamMSE::generate_basic_info(
        n_stocks       = 1,
        n_regions      = 1,
        n_indices      = 2,
        n_fleets       = 2,
        n_seasons      = 1,  
        base.years     = year_start:year_end,
        n_feedback_years = MSE_years,
        n_ages         = 8,
        catch_info     = catch_info,
        index_info     = index_info,
        user_waa       = user_waa,
        user_maturity  = user_maturity,
        fracyr_spawn   = fracyr_spawn
      )
      
      basic_info     <- info$basic_info
      catch_info_use <- info$catch_info
      index_info_use <- info$index_info
      F_info         <- info$F
      F_info$F[1:33,] <- OMa$rep$Fbar[, 1:2]
      
      #### 5. SELECTIVITY/M/NAA_re ####
      sel <- list(n_selblocks = 4,
                  model = c("age-specific","age-specific", "age-specific", "age-specific"))
      sel$initial_pars <- list(
        rep(c(0.5,1),c(3,5)), #north comm
        rep(c(0.5,1),c(6,2)), #north rec
        rep(c(0.5,1,1),c(1,1,6)), #north rec cpa
        rep(c(0.5,1),c(4,4)) #north vast
      )
      sel$fix_pars <- list(
        4:8, #north comm
        7:8, #north rec
        2:8, #north rec cpa
        5:8 #north vast
      )
      
      M <- list(
        model        = "constant",
        initial_means = array(0.4, dim = c(n_stocks, n_regions, n_ages))
      )
      
      vals = exp(OMa$parList$log_NAA_sigma)[1,1,]
      
      ##### SENSITIVITY ANALYSIS POINT 1 ####
      vals[] = proc_error
      
      sigma_vals <- array(vals,
                          dim = c(n_stocks, n_regions, n_ages))
      
      NAA_re <- list(
        recruit_model = 2,
        recruit_pars  = exp(10.5), # fixed mean recruit by stock
        sigma_vals    = sigma_vals,
        sigma         = list("rec+1"),
        cor           = list("2dar1"),
        N1_model      = rep("age-specific-fe", 1)
      )
      
      #### 6. BUILD input_Ecov ####
      
      input_Ecov <- prepare_wham_input(
        basic_info  = basic_info,
        selectivity = sel,
        M           = M,
        NAA_re      = NAA_re,
        ecov        = ecov,
        catch_info  = catch_info_use,
        index_info  = index_info_use,
        F           = F_info
      )
      
      # update WAA pointers
      waa_info   <- info$par_inputs$user_waa
      input_Ecov <- update_waa(input_Ecov, waa_info = waa_info)
      
      # Use OM Estimations for Ecov process pars
      input_Ecov$par$Ecov_process_pars <- OMa$parList$Ecov_process_pars[,1, drop = FALSE] # just north
      input_Ecov$par$Ecov_beta_R       <- OMa$parList$Ecov_beta_R[1,1, ,drop = FALSE]
      
      # For this Gaussian test, zero out the linear Ecov–R effect:
      if (!is.null(input_Ecov$par$Ecov_beta_R)) {
        input_Ecov$par$Ecov_beta_R[] <- 0
      }
      
      #### 7. SET GAUSSIAN RELATIONSHIP WITH TEMPERATURE ####
      input_Ecov$data$use_gauss_T_rec <- 1L         # turn ON Gaussian link
      input_Ecov$data$Ecov_rec_T_col  <- 0L         # first Ecov column = North_BT
      
      temp_col <- input_Ecov$data$Ecov_rec_T_col + 1L
      temp_vec <- input_Ecov$data$Ecov_obs[, temp_col]
      
      # WE need to specify clearly
      input_Ecov$par$Topt_rec      <- 0.0          # peak at 0
      #### SENSITIVITY ANALYSIS POINT 2 ####
      input_Ecov$par$log_width_rec <- log(gauss_width)     # width = 1
      n_stocks <- input_Ecov$data$n_stocks
      input_Ecov$par$beta_T_rec    <- rep(1, n_stocks)
      
      # FIX map for Gaussian T-recruit params 
      for(nm in c("Topt_rec","log_width_rec","beta_T_rec")) {
        if(!is.null(input_Ecov$par[[nm]])) {
          input_Ecov$map[[nm]] <- factor(rep(NA, length(input_Ecov$par[[nm]])))
        }
      }
      
      #### 8. FIX N1, CATCH, AND INDEX[SURVEY] INFO ####
      tmp = OMa$parList$log_NAA[1, 1, 1, ]
      input_Ecov$par$log_N1 <- array(tmp, dim = c(1, 1, length(tmp)))
      
      # index sigma and Neff
      input_Ecov$data$agg_index_sigma[1:33,] <- OMa$input$data$agg_index_sigma[,1:2]
      input_Ecov$data$use_indices[1:33,]     <- OMa$input$data$use_indices[,1:2]
      input_Ecov$data$use_index_paa[1:33,]   <- OMa$input$data$use_index_paa[,1:2]
      
      for (i in 34:(36+12)) {
        input_Ecov$data$agg_index_sigma[i,] <- OMa$input$data$agg_index_sigma[33,1:2, drop = FALSE]
      }
      
      idx1 <- which(asap[[1]]$dat$use_index == 1)
      
      Neff1 <- do.call(cbind, lapply(idx1, function(i)
        asap[[1]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
      index_Neff <- Neff1
      index_Neff <- rbind(index_Neff, index_Neff[rep(33,15), , drop = FALSE])
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
      
      for (i in 34:(36+12)) {
        input_Ecov$data$agg_catch_sigma[i,] <- OMa$input$data$agg_catch_sigma[33,1:2]
      }
      
      Neff1 <- asap[[1]]$dat$catch_Neff
      catch_Neff <- cbind(Neff1)
      catch_Neff <- rbind(catch_Neff, catch_Neff[rep(33,15), , drop = FALSE])
      
      input_Ecov$data$catch_Neff <- catch_Neff
      
      input_Ecov <- update_input_catch_info(
        input_Ecov,
        agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
        catch_Neff      = input_Ecov$data$catch_Neff
      )
      
      #### 9. FORCE Ecov_re pattern: -2 to +2, 1989-2023 ####
      Ecov_re <- OMa$parList$Ecov_re[,1, drop = FALSE]
      
      # Test Only Start Here
      n <- length(Ecov_re)
      
      Ecov_re[,] <- cos(seq(pi, 3*pi, length.out = n))
      Ecov_re
      plot(Ecov_re, type = "l")
      # End here
      
      ny      <- nrow(Ecov_re)
      
      yrs <- ecov$year      # 1959–2025, length 67
      
      # Year indices for 1989–2023 (should be 35 yrs)
      idx <- which(yrs >= 1989 & yrs <= 2023)
      
      # Start with zeros, then fill desired segment
      input_Ecov$par$Ecov_re[1:64,] <- Ecov_re
      
      # Do not simulate Ecov_re – use the values above
      input_Ecov$data$do_simulate_Ecov_re <- 0
      
      # Remove Ecov_re from the list of random effects TMB will estimate
      if ("Ecov_re" %in% input_Ecov$random) {
        input_Ecov$random <- input_Ecov$random[input_Ecov$random != "Ecov_re"]
      }
      
      # par <- input_Ecov$par
      # map <- input_Ecov$map
      # 
      # cat("In par not in map:\n"); print(setdiff(names(par), names(map)))
      # cat("In map not in par:\n"); print(setdiff(names(map), names(par)))
      # 
      # shared <- intersect(names(par), names(map))
      # bad_len <- shared[sapply(shared, function(nm) length(par[[nm]]) != length(map[[nm]]))]
      # if(length(bad_len)) {
      #   print(data.frame(
      #     name = bad_len,
      #     par_len = sapply(bad_len, function(nm) length(par[[nm]])),
      #     map_len = sapply(bad_len, function(nm) length(map[[nm]]))
      #   ))
      # }
      
      #### 10. BUILD OM, REMOVE Ecov_re FROM UNFITTED OM, PLUG IN OUR NAA rho ####
      unfitted_om <- fit_wham(input_Ecov, do.fit = FALSE, do.brps = FALSE,
                              MakeADFun.silent = TRUE)
      
      if ("Ecov_re" %in% unfitted_om$input$random) {
        input_Ecov$random <- unfitted_om$input$random[
          unfitted_om$input$random != "Ecov_re"
        ]
      }
      
      input_Ecov$par$trans_NAA_rho <- OMa$parList$trans_NAA_rho[1,1,,drop = FALSE]
      
      random <- input_Ecov$random
      input_Ecov$random <- NULL
      
      om_ecov <- fit_wham(input_Ecov, do.fit = FALSE, do.brps = TRUE,
                          MakeADFun.silent = TRUE)
      
      # NOTE: OPERATING MODEL IS NOT BEING SAVED
      # saveRDS(om_ecov, file = "om_ecov.rds")
      
      om_with_data <- update_om_fn(om_ecov, seed = 123, random = random)
      
      om_with_data$input$data$agg_catch #simulated catch
      om_with_data$input$data$agg_indices #simulated index
      
      #### 11. DIAGNOSTICS: SHOW GAUSSIAN T-R RELATIONSHIP ####
      # Temperature series (North BT)
      T_series <- om_with_data$rep$Ecov_x[, 1]
      
      # Recruitment (stock 1, region 1, age 1)
      rec1 <- om_with_data$rep$NAA[1,1,,1]
      nyr  <- length(rec1)
      yrs_rec <- year_start:(year_start + nyr - 1)
      
      # par(mfrow = c(2,1), mar = c(4,4,2,1))
      
      plot(yrs, T_series, type = "l", xlab = "Year", ylab = "Ecov_x (North BT)",
           main = "Ecov_x (North BT) with -2 → +2 ramp, 1989–2023")
      abline(v = c(1989, 2023), lty = 2, col = "grey")
      
      plot(yrs_rec, rec1, type = "l", xlab = "Year", ylab = "Recruitment (stock 1)",
           main = "Recruitment vs Gaussian temperature effect")
      
      par(mfrow = c(1,1))
      
      # Optional: overlay theoretical Gaussian scalar (rescaled)
      Topt  <- om_with_data$parList$Topt_rec
      width <- exp(om_with_data$parList$log_width_rec)
      
      g_scalar <- exp(-0.5 * ((T_series - Topt) / width)^2)
      
      # Rescale for visualization
      g_scaled <- max(rec1, na.rm = TRUE) * g_scalar / max(g_scalar, na.rm = TRUE)
      
      lines(yrs, g_scaled, col = "red", lwd = 2)
      legend("topright", legend = c("Recruitment (stock 1)", "Scaled Gaussian(T)"),
             col = c("black","red"), lty = 1, bty = "n")
      
      
      #### 12. SET HARVEST CONTROL RULE (HCR) ####
      # Specify the Harvest Control Rule (HCR)
      hcr <- list()
      hcr$hcr.type <- 1 # FXSPR - Fishing pressure to keep the SPR at a certain percentage
      hcr$hcr.opts <- list(use_FXSPR = TRUE, percentFXSPR = 75) # Apply F at 75% unfished SPR
      
      #### SENSITIVITY ANALYSIS POINT 3 ####
      assess.interval <- mse_gap # Assessments occur every 3 years
      base.years <- year_start:year_end
      terminal.year <- tail(base.years, 1)
      last.year <- 2024+12
      assess.years <- seq(terminal.year, last.year - assess.interval, by = assess.interval)
      
      # Remember to use this config for EM NAA re
      NAA_re_em = NAA_re
      NAA_re_em$N1_model[] = "equilibrium"
      
      # You can use 'est_1' to let EM estimate obs error for Ecov
      ecov_em <- ecov
      ecov_em$logsigma <- 'est_1'
      
      #### 13. SETUP EMs & GATHER THEIR RESULTS ####
      # Execute the MSE loop for one realization
      ##### MODEL 1 - ASSUMES A NO LAG LINEAR RELATIONSHIP BETWEEN  #####
      mod1 <- loop_through_fn(
        om = om_with_data,
        em_info = info,
        random = random,
        sel_em = sel,
        M_em = M,
        ecov_em = ecov_em,
        NAA_re_em = NAA_re_em,
        em.opt = list(
          separate.em = FALSE,
          separate.em.type = 3,
          do.move = TRUE,
          est.move = TRUE
        ),
        update_catch_info = list(agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
                                 catch_Neff = input_Ecov$data$catch_Neff),
        update_index_info = list(agg_index_sigma = input_Ecov$data$agg_index_sigma,
                                 index_Neff = input_Ecov$data$index_Neff),
        assess_years = assess.years,
        assess_interval = assess.interval,
        base_years = base.years,
        year.use = length(base.years),
        add.years = TRUE,
        seed = 123,
        hcr = hcr,
        save.last.em = TRUE # If True, will save all EM information from every iteration, file size can be large, but you can only plot the EM output (using plot_wham_output function) when TRUE...
      )
      
      ##### MODEL 2 - REMOVES THE ENVIRONMENTAL LINKAGE ##### 
      # Execute the MSE loop for one realization
      mod2 <- loop_through_fn(
        om = om_with_data,
        em_info = info,
        random = random,
        sel_em = sel,
        M_em = M,
        # ecov_em = ecov_em,
        NAA_re_em = NAA_re_em,
        em.opt = list(
          separate.em = FALSE,
          separate.em.type = 3,
          do.move = TRUE,
          est.move = TRUE
        ),
        update_catch_info = list(agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
                                 catch_Neff = input_Ecov$data$catch_Neff),
        update_index_info = list(agg_index_sigma = input_Ecov$data$agg_index_sigma,
                                 index_Neff = input_Ecov$data$index_Neff),
        assess_years = assess.years,
        assess_interval = assess.interval,
        base_years = base.years,
        year.use = length(base.years),
        add.years = TRUE,
        seed = 123,
        hcr = hcr,
        save.last.em = TRUE # If True, will save all EM information from every iteration, file size can be large, but you can only plot the EM output (using plot_wham_output function) when TRUE...
      )
      
      ##### MODEL 3 - NO RELATIONSHIP BETWEEN TEMPERATURE AND RECRUITMENT BUT THE OPTION IS THERE #####
      # The model will try to estimate the timeseries of the temperature itself (mean, stdev, autocorrelation)
      ecov_em1 <- ecov_em
      ecov_em1$recruitment_how[] = "none"
      # Execute the MSE loop for one realization
      mod3 <- loop_through_fn(
        om = om_with_data,
        em_info = info,
        random = random,
        sel_em = sel,
        M_em = M,
        ecov_em = ecov_em1,
        NAA_re_em = NAA_re_em,
        em.opt = list(
          separate.em = FALSE,
          separate.em.type = 3,
          do.move = TRUE,
          est.move = TRUE
        ),
        update_catch_info = list(agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
                                 catch_Neff = input_Ecov$data$catch_Neff),
        update_index_info = list(agg_index_sigma = input_Ecov$data$agg_index_sigma,
                                 index_Neff = input_Ecov$data$index_Neff),
        assess_years = assess.years,
        assess_interval = assess.interval,
        base_years = base.years,
        year.use = length(base.years),
        add.years = TRUE,
        seed = 123,
        hcr = hcr,
        save.last.em = TRUE # If True, will save all EM information from every iteration, file size can be large, but you can only plot the EM output (using plot_wham_output function) when TRUE...
      )
      
      ##### MODEL 4 - Assumes Gaussian relationship between temperature and recruitment #####
      # IMPORTANT: As of now (03/18/2026) we don't estimate it directly. We provide direct values for the relationship (Topt, Width)
      mod4 <- loop_through_fn(
        om = om_with_data,
        em_info = info,
        random = random,
        sel_em = sel,
        M_em = M,
        ecov_em = ecov_em1, # remember to set no link between ecov and rec here!
        NAA_re_em = NAA_re_em,
        gauss_rec_em = list( # This will create the linkage
          use = TRUE,
          Ecov_rec_T_col = 1,   # R index, first Ecov column
          Topt_rec = 0.0,
          width_rec = 0.5,
          beta_T_rec = 1,
          estimate = TRUE      # fixed at these values
        ),
        em.opt = list(
          separate.em = FALSE,
          separate.em.type = 3,
          do.move = FALSE,
          est.move = FALSE
        ),
        update_catch_info = list(
          agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
          catch_Neff = input_Ecov$data$catch_Neff
        ),
        update_index_info = list(
          agg_index_sigma = input_Ecov$data$agg_index_sigma,
          index_Neff = input_Ecov$data$index_Neff
        ),
        assess_years = assess.years,
        assess_interval = assess.interval,
        base_years = base.years,
        year.use = length(base.years),
        add.years = TRUE,
        seed = 123,
        hcr = hcr,
        save.last.em = TRUE
      )
      
      #### SAVE MODEL ####
      # FILE/FOLDER STRUCTURE FOR MODEL SAVING: models/sensitivity_analysis
      SAVE_MODEL <- TRUE
      if(SAVE_MODEL){
        # Create a folder for saving all the data and run information
        # We will use the date+time from the start time of the code
        # Format run_start_time as a posixDate object
        folder_name <- format(run_start_time, "%Y-%m-%d_%H-%M-%S")
        folder_path <- here("models","sensitivity_analysis",folder_name,"models")
        # Create a folder. Suppress warnings and allow recursive folders to be created
        dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
        saveRDS(mod1, here(folder_path,paste("sens_run_",run_id,"_mod_1",".RDS",sep="")))
        saveRDS(mod2, here(folder_path,paste("sens_run_",run_id,"_mod_2",".RDS",sep="")))
        saveRDS(mod3, here(folder_path,paste("sens_run_",run_id,"_mod_3",".RDS",sep="")))
        saveRDS(mod4, here(folder_path,paste("sens_run_",run_id,"_mod_4",".RDS",sep="")))
        # for(iter in seq(iterations)){
        #   saveRDS(model_list[iter], here(folder_path,paste("model_run_",iter,".RDS",sep="")))
        # }
        print(paste("Models saved for Run",run_id,sep=" "))
      }
      
      #### SAVE PLOTS ####
      if(SAVE_MODEL){
        folder_path <- here("models","sensitivity_analysis",folder_name,"plots")
        dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
        ##### OPERATING MODEL - ABUNDANCE[SSB] #####
        plot(mod1$om$rep$SSB, type = "l", col = "red")
        lines(mod2$om$rep$SSB, type = "l", col = "blue")
        
        # Getting these into a nice tidyverse format
        mod1_om_ssb <- mod1$om$rep$SSB
        mod1_om_ssb <- cbind(mod1_om_ssb,model='Model 1')
        mod1_om_ssb <- cbind(ID = as.integer(1:nrow(mod1_om_ssb)), mod1_om_ssb)
        mod2_om_ssb <- mod2$om$rep$SSB
        mod2_om_ssb <- cbind(mod2_om_ssb,model='Model 2')
        mod2_om_ssb <- cbind(ID = as.integer(1:nrow(mod2_om_ssb)), mod2_om_ssb)
        mod3_om_ssb <- mod3$om$rep$SSB
        mod3_om_ssb <- cbind(mod3_om_ssb,model='Model 3')
        mod3_om_ssb <- cbind(ID = as.integer(1:nrow(mod3_om_ssb)), mod3_om_ssb)
        mod4_om_ssb <- mod4$om$rep$SSB
        mod4_om_ssb <- cbind(mod4_om_ssb,model='Model 4')
        mod4_om_ssb <- cbind(ID = as.integer(1:nrow(mod4_om_ssb)), mod4_om_ssb)
        
        
        mod_om_ssb <- rbind(mod1_om_ssb, mod2_om_ssb, mod3_om_ssb, mod4_om_ssb)
        mod_om_ssb_df <- tibble::as_tibble(mod_om_ssb)
        mod_om_ssb_df <- mod_om_ssb_df %>% rename(year=ID) %>% mutate(year=as.integer(year)+1988)
        mod_om_ssb_df <- mod_om_ssb_df %>% rename(SSB=V2) %>% mutate(SSB=as.numeric(SSB))
        
        # Plot via ggplot
        om_ssb_plot_1 <- ggplot(mod_om_ssb_df, aes(year, SSB, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + facet_wrap(~model, nrow=2) + 
          scale_x_continuous(breaks=seq(1985,2040,5)) + 
          labs(color="Model", title=paste("SSB in OM: ","sigma_naa=",proc_error, ", mse_gap=",mse_gap, ", gaussian_width=",gauss_width, sep=" ")) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
        
        # Plot via ggplot
        om_ssb_plot_2 <- ggplot(mod_om_ssb_df, aes(year, SSB, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + 
          labs(color="Model") + 
          theme_bw()
        
        #### OPERATING MODEL - PREDICTED CATCH ####
        mod1_om_pred_catch <- data.frame(PRED_CATCH=mod1$om$rep$pred_catch[,1])
        mod1_om_pred_catch <- cbind(mod1_om_pred_catch,model='Model 1')
        mod1_om_pred_catch <- cbind(ID = as.integer(1:nrow(mod1_om_pred_catch)), mod1_om_pred_catch)
        mod2_om_pred_catch <- data.frame(PRED_CATCH=mod2$om$rep$pred_catch[,1])
        mod2_om_pred_catch <- cbind(mod2_om_pred_catch,model='Model 2')
        mod2_om_pred_catch <- cbind(ID = as.integer(1:nrow(mod2_om_pred_catch)), mod2_om_pred_catch)
        mod3_om_pred_catch <- data.frame(PRED_CATCH=mod3$om$rep$pred_catch[,1])
        mod3_om_pred_catch <- cbind(mod3_om_pred_catch,model='Model 3')
        mod3_om_pred_catch <- cbind(ID = as.integer(1:nrow(mod3_om_pred_catch)), mod3_om_pred_catch)
        mod4_om_pred_catch <- data.frame(PRED_CATCH=mod4$om$rep$pred_catch[,1])
        mod4_om_pred_catch <- cbind(mod4_om_pred_catch,model='Model 4')
        mod4_om_pred_catch <- cbind(ID = as.integer(1:nrow(mod4_om_pred_catch)), mod4_om_pred_catch)
        
        mod_om_pred_catch <- rbind(mod1_om_pred_catch, mod2_om_pred_catch, 
                                   mod3_om_pred_catch, mod4_om_pred_catch)
        mod_om_pred_catch_df <- tibble::as_tibble(mod_om_pred_catch)
        mod_om_pred_catch_df <- mod_om_pred_catch_df %>% rename(year=ID) %>% mutate(year=as.integer(year)+1988)
        
        # Plot via ggplot
        om_pred_catch_plot_1 <- ggplot(mod_om_pred_catch_df, aes(year, PRED_CATCH, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + facet_wrap(~model, nrow=2) + 
          scale_x_continuous(breaks=seq(1985,2040,5)) + 
          labs(color="Model", title=paste("PRED. CATCH in OM: ","sigma_naa=",proc_error, ", mse_gap=",mse_gap, ", gaussian_width=",gauss_width, sep=" ")) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
        
        # Plot via ggplot
        om_pred_catch_plot_2 <- ggplot(mod_om_pred_catch_df, aes(year, PRED_CATCH, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + 
          labs(color="Model") + 
          theme_bw()
        
        #### ESTIMATION MODEL - ABUNDANCE [SSB] ####
        # lines(mod1$em_full[[1]]$rep$SSB, col = "red")
        mod1_em_ssb <- data.frame(SSB=mod1$em_full[[1]]$rep$SSB)
        mod1_em_ssb <- cbind(mod1_em_ssb,model='Model 1')
        mod1_em_ssb <- cbind(ID = as.integer(1:nrow(mod1_em_ssb)), mod1_em_ssb)
        mod2_em_ssb <- data.frame(SSB=mod2$em_full[[1]]$rep$SSB)
        mod2_em_ssb <- cbind(mod2_em_ssb,model='Model 2')
        mod2_em_ssb <- cbind(ID = as.integer(1:nrow(mod2_em_ssb)), mod2_em_ssb)
        mod3_em_ssb <- data.frame(SSB=mod3$em_full[[1]]$rep$SSB)
        mod3_em_ssb <- cbind(mod3_em_ssb,model='Model 3')
        mod3_em_ssb <- cbind(ID = as.integer(1:nrow(mod3_em_ssb)), mod3_em_ssb)
        mod4_em_ssb <- data.frame(SSB=mod4$em_full[[1]]$rep$SSB)
        mod4_em_ssb <- cbind(mod4_em_ssb,model='Model 4')
        mod4_em_ssb <- cbind(ID = as.integer(1:nrow(mod4_em_ssb)), mod4_em_ssb)
        
        
        
        mod_em_ssb <- rbind(mod1_em_ssb, mod2_em_ssb, mod3_em_ssb, mod4_em_ssb)
        mod_em_ssb_df <- tibble::as_tibble(mod_em_ssb)
        mod_em_ssb_df <- mod_em_ssb_df %>% rename(year=ID) %>% mutate(year=as.integer(year)+1988)
        
        # Plot via ggplot
        em_ssb_plot_1 <- ggplot(mod_em_ssb_df, aes(year, SSB, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + facet_wrap(~model, nrow=2) + 
          scale_x_continuous(breaks=seq(1985,2040,5)) + 
          labs(color="Model", title=paste("SSB in EM: ","sigma_naa=",proc_error, ", mse_gap=",mse_gap, ", gaussian_width=",gauss_width, sep=" ")) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
        
        # Plot via ggplot
        em_ssb_plot_2 <- ggplot(mod_em_ssb_df, aes(year, SSB, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + 
          labs(color="Model") + 
          theme_bw()
        
        #### NAA - Age class 1 EM vs. OM  ####
        # NAA - Age class 1
        # EM vs. OM
        mod1_om_NAA <- data.frame(NAA_1=mod1$om$rep$NAA[,,,1])
        mod1_om_NAA <- cbind(mod1_om_NAA,model='Model 1 - OM')
        mod1_om_NAA <- cbind(ID = as.integer(1:nrow(mod1_om_NAA)), mod1_om_NAA)
        mod1_em_NAA <- data.frame(NAA_1=mod1$em_full[[1]]$rep$NAA[,,,1])
        mod1_em_NAA <- cbind(mod1_em_NAA,model='Model 1 - EM')
        mod1_em_NAA <- cbind(ID = as.integer(1:nrow(mod1_em_NAA)), mod1_em_NAA)
        mod2_em_NAA <- data.frame(NAA_1=mod2$em_full[[1]]$rep$NAA[,,,1])
        mod2_em_NAA <- cbind(mod2_em_NAA,model='Model 2 - EM')
        mod2_em_NAA <- cbind(ID = as.integer(1:nrow(mod2_em_NAA)), mod2_em_NAA)
        mod3_em_NAA <- data.frame(NAA_1=mod3$em_full[[1]]$rep$NAA[,,,1])
        mod3_em_NAA <- cbind(mod3_em_NAA,model='Model 3 - EM')
        mod3_em_NAA <- cbind(ID = as.integer(1:nrow(mod3_em_NAA)), mod3_em_NAA)
        mod4_em_NAA <- data.frame(NAA_1=mod4$em_full[[1]]$rep$NAA[,,,1])
        mod4_em_NAA <- cbind(mod3_em_NAA,model='Model 4 - EM')
        mod4_em_NAA <- cbind(ID = as.integer(1:nrow(mod4_em_NAA)), mod4_em_NAA)
        
        mod_NAA <- rbind(mod1_om_NAA, mod1_em_NAA, mod2_em_NAA, mod3_em_NAA, mod4_em_NAA)
        mod_NAA_df <- tibble::as_tibble(mod_NAA)
        mod_NAA_df <- mod_NAA_df %>% rename(year=ID) %>% mutate(year=as.integer(year)+1988)
        
        # Plot via ggplot
        mod_NAA_plot_1 <- ggplot(mod_NAA_df, aes(year, NAA_1, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + facet_wrap(~model, nrow=2) + 
          scale_x_continuous(breaks=seq(1985,2040,5)) + 
          labs(color="Model", title=paste("NAA 1 in OM & EM: ","sigma_naa=",proc_error, ", mse_gap=",mse_gap, ", gaussian_width=",gauss_width, sep=" ")) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
        
        # Plot via ggplot
        mod_NAA_plot_2 <- ggplot(mod_NAA_df, aes(year, NAA_1, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + 
          scale_x_continuous(breaks=seq(1985,2040,5)) + 
          labs(color="Model") + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
        
        #### FISHING PRESSURE ####
        mod1_om_fbar <- data.frame(fbar=mod1$om$rep$Fbar[,1])
        mod1_om_fbar <- cbind(mod1_om_fbar,model='Model 1 - OM')
        mod1_om_fbar <- cbind(ID = as.integer(1:nrow(mod1_om_fbar)), mod1_om_fbar)
        mod1_em_fbar <- data.frame(fbar=mod1$em_full[[1]]$rep$Fbar[,1])
        mod1_em_fbar <- cbind(mod1_em_fbar,model='Model 1 - EM')
        mod1_em_fbar <- cbind(ID = as.integer(1:nrow(mod1_em_fbar)), mod1_em_fbar)
        mod2_em_fbar <- data.frame(fbar=mod2$em_full[[1]]$rep$Fbar[,1])
        mod2_em_fbar <- cbind(mod2_em_fbar,model='Model 2 - EM')
        mod2_em_fbar <- cbind(ID = as.integer(1:nrow(mod2_em_fbar)), mod2_em_fbar)
        mod3_em_fbar <- data.frame(fbar=mod3$em_full[[1]]$rep$Fbar[,1])
        mod3_em_fbar <- cbind(mod3_em_fbar,model='Model 3 - EM')
        mod3_em_fbar <- cbind(ID = as.integer(1:nrow(mod3_em_fbar)), mod3_em_fbar)
        mod4_em_fbar <- data.frame(fbar=mod4$em_full[[1]]$rep$Fbar[,1])
        mod4_em_fbar <- cbind(mod4_em_fbar,model='Model 4 - EM')
        mod4_em_fbar <- cbind(ID = as.integer(1:nrow(mod4_em_fbar)), mod4_em_fbar)
        
        mod_fbar <- rbind(mod1_om_fbar, mod1_em_fbar, mod2_em_fbar, mod3_em_fbar, mod4_em_fbar)
        mod_fbar_df <- tibble::as_tibble(mod_fbar)
        mod_fbar_df <- mod_fbar_df %>% rename(year=ID) %>% mutate(year=as.integer(year)+1988)
        
        # Plot via ggplot
        mod_fbar_plot_1 <- ggplot(mod_fbar_df, aes(year, fbar, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + facet_wrap(~model, nrow=2) + 
          scale_x_continuous(breaks=seq(1985,2040,5)) + 
          labs(color="Model", title=paste("Fbar in OM & EM: ","sigma_naa=",proc_error, ", mse_gap=",mse_gap, ", gaussian_width=",gauss_width, sep=" ")) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
        
        # Plot via ggplot
        mod_fbar_plot_2 <- ggplot(mod_fbar_df, aes(year, fbar, color=as.factor(model))) + 
          geom_line(alpha=0.5, linewidth=1) + 
          scale_x_continuous(breaks=seq(1985,2040,5)) + 
          labs(color="Model") + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
        
        #### PARAMETER ESTIMATES AND DIFFERENCES - STILL TO BE DONE! ####
        
        # Arrange and save the plots to a multi-page PDF (4 plots per page: 2 rows, 2 columns)
        # pdf(here(folder_path,paste("om_ssb_plot_run_",run_id,".pdf",sep="")), width = 12, height = 8)
        # pdf_path <- here(folder_path, "run_id.pdf")
        #  pdf(pdf_path, width=12, height=8)
        #  marrangeGrob(test_plots_list, nrow = 2, ncol = 1)
        #  dev.off()
        #  
        # Append the 
        diagnostics_plots <- list(om_ssb_plot_1, om_ssb_plot_2, 
                                  om_pred_catch_plot_1, om_pred_catch_plot_2,
                                  em_ssb_plot_1, em_ssb_plot_2,
                                  mod_NAA_plot_1, mod_NAA_plot_2,
                                  mod_fbar_plot_1, mod_fbar_plot_2)
        
        multi.page <- ggpubr::ggarrange(plotlist = diagnostics_plots, nrow = 2, ncol = 1) 
        ggpubr::ggexport(multi.page, filename = pdf_path <- here(folder_path, paste("sens_run_",run_id,"_plots.pdf",sep="")))
        print(paste("Done with NID:",run_id))
        beepr::beep(1)
      } # PLOTTING IF ENDS
    },
    error=function(e){
      skip_to_next <<- TRUE
    })
    if(skip_to_next){
      next
    }
    ##### SAVE MODEL RUN SETTINGS  #####
  } # FOR LOOP ENDS
} # IF ENDS

model_csv <- c("mod_1", "mod_2", "mod_3", "mod_4")
sens_analysis_to_csv <- tidyr::crossing(proc_error= proc_error_v, 
                                        mse_gap = mse_gaps_v, 
                                        gauss_width=gauss_width_v,
                                        model=model_csv)

sens_analysis_to_csv <- sens_analysis_to_csv %>% mutate(run_id=rep(1:18, each=4), .before=1)
sens_analysis_to_csv <- sens_analysis_to_csv %>% mutate(model_path=here(folder_path,paste("sens_run_",run_id,"_",model,".RDS",sep="")))
folder_path <- here("models","sensitivity_analysis",folder_name)
write_csv(sens_analysis_to_csv, here(folder_path, "all_model_settings.csv"))

print("Done with the whole sensitivity run")
run_end_time <- Sys.time()
run_total_time_hours <- floor(as.numeric(difftime(run_end_time, run_start_time, units=c("hours"))))
run_total_time_mins <- floor(as.numeric(difftime(run_end_time, run_start_time, units=c("mins"))) %% 60)
run_total_time_secs <- floor(as.numeric(difftime(run_end_time, run_start_time, units=c("secs"))) %% 60)
print(paste("Total execution time was", run_total_time_hours, "hours and", 
            run_total_time_mins, "minutes and", run_total_time_secs, 
            "seconds", sep=" "))
beepr::beep(2)


# Gather the model objects to a single list and name them
mod_list <- list(mod1, mod2, mod3, mod4)
names(mod_list) <- c("Model 1", "Model 2", "Model 3", "Model 4")
mod_list_2 <- list(mod1, mod2, mod3, mod4)
names(mod_list_2) <- c("Model 1", "Model 2", "Model 3", "Model 4")

gather_mod_list <- list(mod_list, mod_list_2)
names(gather_mod_list) <- c("List 1", "List 2")

# Access the 

#### DISCARD EVERYTHING BELOW ####

#### EXPERIMENTAL MINI CHUNK HERE ####
# Get models into a list
# Create a list of model objects
mod_list <- list(mod_names=c("Model 1", "Model 2", "Model 3"))
mod_list$models <- list(mod1, mod2, mod3)

#### OPERATING MODEL - ABUNDANCE[SSB] ####
plot(mod1$om$rep$SSB, type = "l", col = "red")
lines(mod2$om$rep$SSB, type = "l", col = "blue")

# Getting these into a nice tidyverse format
mod1_om_ssb <- mod1$om$rep$SSB
mod1_om_ssb <- cbind(mod1_om_ssb,model='Model 1')
mod1_om_ssb <- cbind(ID = as.integer(1:nrow(mod1_om_ssb)), mod1_om_ssb)
mod2_om_ssb <- mod2$om$rep$SSB
mod2_om_ssb <- cbind(mod2_om_ssb,model='Model 2')
mod2_om_ssb <- cbind(ID = as.integer(1:nrow(mod2_om_ssb)), mod2_om_ssb)
mod3_om_ssb <- mod3$om$rep$SSB
mod3_om_ssb <- cbind(mod3_om_ssb,model='Model 3')
mod3_om_ssb <- cbind(ID = as.integer(1:nrow(mod3_om_ssb)), mod3_om_ssb)


mod_om_ssb <- rbind(mod1_om_ssb, mod2_om_ssb, mod3_om_ssb)
mod_om_ssb_df <- tibble::as_tibble(mod_om_ssb)
mod_om_ssb_df <- mod_om_ssb_df %>% rename(year=ID) %>% mutate(year=as.integer(year)+1988)
mod_om_ssb_df <- mod_om_ssb_df %>% rename(SSB=V2) %>% mutate(SSB=as.numeric(SSB))

# Plot via ggplot
ggplot(mod_om_ssb_df, aes(year, SSB, color=as.factor(model))) + 
  geom_line(alpha=0.5, linewidth=1) + facet_wrap(~model, nrow=2) + 
  scale_x_continuous(breaks=seq(1985,2040,5)) + 
  labs(color="Model") + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))

# Plot via ggplot
ggplot(mod_om_ssb_df, aes(year, SSB, color=as.factor(model))) + 
  geom_line(alpha=0.5, linewidth=1) + 
  labs(color="Model") + 
  theme_bw()


#### PREDICTED CATCH -  ####
# Fleet 1
plot(mod1$om$rep$pred_catch[,1], type = "l", col = "red")
lines(mod2$om$rep$pred_catch[,1], type = "l", col = "blue")

# Fleet 2
plot(mod1$om$rep$pred_catch[,2], type = "l", col = "red")
lines(mod2$om$rep$pred_catch[,2], type = "l", col = "blue")

#### ESTIMATION MODEL - ABUNDANCE[SSB] ####
lines(mod1$em_full[[1]]$rep$SSB, col = "red")
lines(mod2$em_full[[1]]$rep$SSB, col = "blue")
lines(mod3$em_full[[1]]$rep$SSB, col = "purple")


#### DIFFERENCE BETWEEN ESTIMATION MODELS ####
mod1$em_full[[1]]$rep$SSB - mod2$em_full[[1]]$rep$SSB

mod2$em_full[[1]]$parList$mean_rec_pars - mod1$em_full[[1]]$parList$mean_rec_pars
mod2$em_full[[1]]$parList$Ecov_beta_R - mod1$em_full[[1]]$parList$Ecov_beta_R

# NAA - Age class 1
# EM vs. OM
plot(mod1$om$rep$NAA[,,,1], type = "l")
lines(mod1$em_full[[1]]$rep$NAA[,,,1], col = "red")

# Fishing pressure
plot(mod1$om$rep$Fbar[,1], type = "l")
lines(mod1$em_full[[1]]$rep$Fbar[,1], col = "red")
lines(mod2$em_full[[1]]$rep$Fbar[,1], col = "blue")
lines(mod3$em_full[[1]]$rep$Fbar[,1], col = "green")

plot(mod1$om$rep$Fbar[,2], type = "l")
lines(mod1$em_full[[1]]$rep$Fbar[,2], col = "red")
lines(mod2$em_full[[1]]$rep$Fbar[,2], col = "blue")
lines(mod3$em_full[[1]]$rep$Fbar[,2], col = "green")

plot(mod1$om$rep$Fbar[,3], type = "l")
lines(mod1$em_full[[1]]$rep$Fbar[,3], col = "red")
lines(mod2$em_full[[1]]$rep$Fbar[,3], col = "blue")
lines(mod3$em_full[[1]]$rep$Fbar[,3], col = "green")


plot(mod1$om$rep$Fbar[,4], type = "l")
lines(mod1$em_full[[1]]$rep$Fbar[,4], col = "red")
lines(mod2$em_full[[1]]$rep$Fbar[,4], col = "blue")
lines(mod3$em_full[[1]]$rep$Fbar[,4], col = "green")

# Predicted vs. estimated catch
plot(mod1$om$rep$pred_catch[,1], type = "l")
lines(mod1$em_full[[1]]$rep$pred_catch[,1], col = "red")

mod1$em_input[[1]]$par$Ecov_beta_R
mod2$em_input[[1]]$par$Ecov_beta_R

mod1$em_input[[1]]$par$Ecov_process_pars
mod2$em_input[[1]]$par$Ecov_process_pars

mod1$em_full[[1]]$parList$mean_rec_pars - mod2$em_full[[1]]$parList$mean_rec_pars
mod1$em_full[[1]]$parList$logit_q - mod2$em_full[[1]]$parList$logit_q

mod1$em_full[[1]]$sdrep
mod2$em_full[[1]]$sdrep

mod1$em_full[[1]]$parList$Ecov_beta_R
mod2$em_full[[1]]$parList$Ecov_beta_R

# Population correlation coefficient
mod1$om$parList$trans_NAA_rho[,,1]
mod1$em_full[[1]]$parList$trans_NAA_rho
mod2$em_full[[1]]$parList$trans_NAA_rho
mod3$em_full[[1]]$parList$trans_NAA_rho

# Is recruitment connected to ecovariate?
mod1$em_full[[1]]$input$data$Ecov_how_R
mod2$em_full[[1]]$input$data$Ecov_how_R
mod3$em_full[[1]]$input$data$Ecov_how_R

# EM - Negative log likelihood 
mod1$em_full[[1]]$rep$nll
mod2$em_full[[1]]$rep$nll
mod3$em_full[[1]]$rep$nll

plots_list <- lapply(1:10, function(i) {
  ggplot(mtcars, aes(wt, mpg)) +
    geom_point() +
    ggtitle(paste("Plot", i))
})

# How to add plots to an existing list
sample_plot <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point(color="red") +
  ggtitle(paste("Plot", i))
sample_plot_2 <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point(color="yellow") +
  ggtitle(paste("Plot", i))

plots_list_x <- list(sample_plot, sample_plot_2)
plots_list_x <- append(plots_list_x, list(sample_plot))

# Arrange and save the plots to a multi-page PDF (4 plots per page: 2 rows, 2 columns)
pdf(here("plots","multi_page_gridExtra2.pdf"), width = 12, height = 8)
marrangeGrob(plots_list_x, nrow = 2, ncol = 2)
dev.off()


#### DATAFRAME FOR STORING PARAMETER ESTIMATES AND DIFFERENCES ####
# Required columns
# 1. RUN_ID
# 2. proc_error
# 3. mse_gap
# 4. gauss_width
# 5. parameter
# 6. OM or EM
# 7. diff or raw_value
# 8. desc - 'Show the equation here' or the model number
# 9. value

list_1 <- list(run_id=1, proc_error=0.1, mse_gap=3, gauss_width=0.1,
               parameter="SSB",om_em="OM",
               diff_or_raw="raw",desc="Model 1", value=2.5e10)

list_2 <- list(run_id=2, proc_error=0.1, mse_gap=6, gauss_width=0.1,
               parameter="SSB",om_em="OM",
               diff_or_raw="diff",desc="Model 1 - Model 2", value=1.3e10)

list_3 <- list(run_id=3, proc_error=0.1, mse_gap=rep(5,30), gauss_width=0.1,
               parameter="SSB",om_em="OM",
               diff_or_raw="raw",desc="Model 1 - Model 2", value=1.3e10)
# NOTE: `as.data.frame(list_3)` works

all_lists <- list(list_1, list_2)

param_df <- do.call(rbind.data.frame, all_lists)

# Binding lists of unequal length

all_lists_2 <- list(list_1, list_2, list_3)

# The `rbind.data.frame` call doesn't work on unequal lists
param_df_2 <- do.call(rbind.data.frame, all_lists_2)

# But data.table::rbindlist does!
data.table::rbindlist(list(list_1, list_2, list_3), fill = TRUE)

View(data.table::rbindlist(list(list_1, list_2, list_3), fill = TRUE))

vec <- 1:18

# Repeat each element of the vector three times
result <- rep(vec, each = 3) 

# Print the result
print(result)


# Handling errors

safe_log <- function(x) {
  tryCatch(
    expr = {
      log(x)
    },
    error = function(e) {
      message(paste("An error occurred in the code:", conditionMessage(e)))
      return(NA) # Return NA in case of an error
    }
  )
}

# Example usage:
print(safe_log(10))
print(safe_log("a")) # This will trigger the error handler

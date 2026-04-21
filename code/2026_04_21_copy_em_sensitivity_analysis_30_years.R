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
# NOTE: Introducing parallelization here
#
# -----------------------------------------------------------------------------

# devtools::install_github("lichengxue/wham@Gaussian_test")
# devtools::install_github("lichengxue/whamMSE@Projection-MSE")
library(wham)
library(whamMSE)
library(dplyr)
library(here)
library(ggplot2)
library(gridExtra)
library(beepr)
library(readr)
library(foreach)
library(doParallel)

here::here()

# NOTE: This code is set to run on a server as well as a personal computer (MacOS for now)

#### RUN TIMES ####
# Benchmarking run time and info for save files
run_start_time <- Sys.time()

#### RUN ENVIRONMENT ####
# Options for run environment
run_env_opts <- c("local","annotate2","amarel")
run_env <- run_env_opts[1]

#### PARALLEL COMPUTATIONS ####
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

#### MODEL RUN SETTINGS ####
# Set iterations, a base random seed, and then generate seeds for each MSE run
# Set a model name
# NOTE: NOT USING THESE SETTINGS FOR THE SENSITIVITY ANALYSIS
iterations <- 10
base_random_seed <- 853
set.seed(base_random_seed)
mse_random_seeds <- as.integer(floor(runif(iterations, min=0, max=1000)))
model_name <- "BSB Ecov"
n_feedback_years <- 15

# Make the seeds into a dataframe and save
random_seeds_df <- data.frame(n_seed=mse_random_seeds) %>% mutate(nid=row_number(), .before=1)

#### SENSITIVITY ANALYSIS SETTINGS ####
# We will be looking at model variation from the following
# Process error for recruitment and NAA - vals - 0.2, 0.5, 1
# MSE gaps - 3 years, 6 years
# Width of gaussian relationship - `input_Ecov$par$log_width_rec` - 0.1 [Very stringent], 1, 5 [Very forgiving to non optimal temperatures]
# Variability of temperature - Mean centered around zero -  increasing stochasticity - 0.5, 1, 3
# NOTE: Variability of temperature is removed for now. Integration possible in the future

proc_error_v <- c(0.2) # Process error for NAA random effects. We keep a low value here to see that allows us to see difference in model performance
mse_gaps_v <- c(3,6) # Time between assessments for MSE
gauss_width_v <- c(0.5) # Width of the gaussian relationship (Wider = Less sensitive to optimal temperature)
total_comb_no <- length(proc_error_v)*length(mse_gaps_v)*length(gauss_width_v)

# Use `crossing` function from tidyr to create a dataframe of all the settings in the
# sensitivity analysis
sens_analysis_settings <- tidyr::crossing(proc_error= proc_error_v,
                                          mse_gap = mse_gaps_v,
                                          gauss_width=gauss_width_v)

sens_analysis_settings <- sens_analysis_settings %>% mutate(nid=row_number(), .before=1)

#### IF THIS IS A TEST RUN - YOU ONLY WANT TO RUN A COUPLE OF ROWS!!! ####
# NOTE: ONLY USE X NO. OF ROWS BECAUSE THIS IS A TEST RUN
# sens_analysis_settings <- sens_analysis_settings %>% head(n=1)

#### EXTERNAL FUNCTIONS ####
# Source reusable functions from `functions/reusable_functions.R`
source(here("functions","reusable_functions.R"))

SAVE_MODEL <- TRUE
if(SAVE_MODEL){
  # Create a folder for saving all the data and run information
  # We will use the date+time from the start time of the code
  # Format run_start_time as a posixDate object
  folder_name <- format(run_start_time, "%Y-%m-%d_%H-%M-%S")
  folder_path <- here("models","sensitivity_analysis",folder_name,"models")
  # Create a folder. Suppress warnings and allow recursive folders to be created
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  # Folder for plots
  folder_path_plots <- here("models","sensitivity_analysis",folder_name,"plots")
  dir.create(folder_path_plots, recursive = TRUE, showWarnings = FALSE)
  # Folder for diagnostics
  folder_path_diagnostics <- here("models","sensitivity_analysis",folder_name,"diagnostics")
  dir.create(folder_path_diagnostics, recursive = TRUE, showWarnings = FALSE)
}

#### DO THE SENSITIVITY ANALYSIS OR NOT? ####
# WARNING: THIS CAN TAKE A LOT OF TIME
DO_ANALYSIS = TRUE

##### IMPORTANT ####
# seed = base_random_seed
# seed = 123
##### You should make sure seed is used consistently for each EM, not a big problem here, but when you run 100 realizations...

if(DO_ANALYSIS){
  for(j in 1:nrow(sens_analysis_settings)) {
    # Envelope in a tryCatch so that we just skip over stuff that doesn't work
    foreach(k = 1:nrow(random_seeds_df)) %dopar% {
      tryCatch({
        skip_to_next <- FALSE # A flag we use to deal with convergence errors
        # Get the settings for this sensitivity run
        row <- sens_analysis_settings[j,]
        proc_error <- row$proc_error
        mse_gap <- row$mse_gap
        gauss_width <- row$gauss_width
        run_id <- row$nid
        print(paste("Running sensitivity analysis Run:",run_id, sep=" "))
        
        # Get and set the random seed
        seed_row <- random_seeds_df[k,]
        iter_id <- seed_row$nid
        seed <- seed_row$n_seed
        print(paste("Seed is now: ",seed))
        
        # NOTE: OM model ends in 2021 and EM starts in 2022
        # But the bottom temperature dataset goes until 2022
        # n_feedback_years <- 30
        
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
        # Get the final year
        ecov_final_year <- max(north_bt[,"year"])
        projection_years <- seq(ecov_final_year+1, ecov_final_year+n_feedback_years)
        ecov$year <- c(north_bt[,"year"], projection_years)
        
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
        MSE_years  <- n_feedback_years
        
        hist_years <- length(year_start:year_end)
        total_years <- length(year_start:(year_end + MSE_years))
        om_future_start_index <- hist_years+1
        om_future_end_index <- hist_years+n_feedback_years
        
        # maturity
        user_maturity <- array(NA, dim = c(n_stocks, total_years, n_ages))
        user_maturity[, 1:hist_years, ] <- OMa$input$data$mature[1,,]
        for (i in (hist_years+1):(hist_years+MSE_years)) {
          user_maturity[, i, ] <- OMa$input$data$mature[1, hist_years, , drop = FALSE]
        }
        
        
        user_waa <- list()
        user_waa$waa <- array(NA, dim = c(5, hist_years + n_feedback_years, 8))
        user_waa$waa[, 1:hist_years, ] <- OMa$input$data$waa[c(1,2,5,6,9), ,]
        for (i in om_future_start_index:om_future_end_index) {
          user_waa$waa[, i, ] <- OMa$input$data$waa[c(1,2,5,6,9), hist_years, ]
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
        F_info$F[1:hist_years,] <- OMa$rep$Fbar[, 1:2]
        
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
          recruit_pars  = exp(12), # fixed mean recruit by stock
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
        
        # IMPORTANT! #
        # We can also set 75% of F40% to be F default values in the feedback! So you have another BASELINE
        om_ecov$input$par$F_pars[om_future_start_index:om_future_end_index,] = log(exp(c(-1.167205, -1.167205))*0.75)
        om_with_data <- update_om_fn(om_ecov, seed = seed, random = random)
        # om_ecov$parList$F_pars[om_future_start_index:om_future_end_index,] = om_with_data$rep$log_SPR_FXSPR_static
        
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
        
        # Realized recruitment / age-1 abundance
        rec1_obs <- om_with_data$rep$NAA[1, 1, , 1]
        
        # Model-predicted recruitment / age-1 abundance (deterministic center)
        rec1_pred <- om_with_data$rep$pred_NAA[1, 1, , 1]
        
        nyr <- length(rec1_obs)
        yrs_rec <- year_start:(year_start + nyr - 1)
        
        plot(yrs_rec, rec1_obs, type = "l", lwd = 1,
             xlab = "Year", ylab = "Recruitment (stock 1, region 1, age 1)",
             main = "Realized vs model-predicted recruitment")
        
        lines(yrs_rec, rec1_pred, col = "red", lwd = 2)
        
        legend("topright",
               legend = c("Realized recruitment", "Predicted recruitment"),
               col = c("black", "red"), lty = 1, bty = "n")
        
        # Convert these into nicer ggplots
        # Temperature over time
        
        
        #### 12. SET HARVEST CONTROL RULE (HCR) ####
        # Specify the Harvest Control Rule (HCR)
        hcr <- list()
        hcr$hcr.type <- 1 # FXSPR - Fishing pressure to keep the SPR at a certain percentage
        hcr$hcr.opts <- list(use_FXSPR = TRUE, percentFXSPR = 75) # Apply F at 75% unfished SPR
        
        #### SENSITIVITY ANALYSIS POINT 3 ####
        assess.interval <- mse_gap # Assessments occur every 3 or 6 years. Depends on the configuration
        base.years <- year_start:year_end
        terminal.year <- tail(base.years, 1)
        last.year <- terminal.year+n_feedback_years
        assess.years <- seq(terminal.year, last.year - assess.interval, by = assess.interval)
        
        # ggplot assess.years component
        assess_years_lines <- geom_vline(xintercept=assess.years, linetype="dashed", color="grey", alpha=0.75)
        
        # Remember to use this config for EM NAA re
        NAA_re_em = NAA_re
        NAA_re_em$N1_model[] = "equilibrium"
        
        # You can use 'est_1' to let EM estimate obs error for Ecov
        ecov_em <- ecov
        ecov_em$logsigma <- 'est_1'
        
        #### 13. SETUP EMs & GATHER THEIR RESULTS ####
        # Create empty objects for models 1,2, and 3
        mod1 <- list()
        mod2 <- list()
        mod3 <- list()
        # Execute the MSE loop for one realization
        ##### MODEL 1 - ASSUMES A NO LAG LINEAR RELATIONSHIP BETWEEN TEMP AND RECRUITMENT #####
        mod1_flag <- FALSE
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
          seed = seed,
          hcr = hcr,
          save.last.em = TRUE # If True, will save all EM information from every iteration, file size can be large, but you can only plot the EM output (using plot_wham_output function) when TRUE...
        )
        mod1_flag <- TRUE
        print("Done with linear model!")
        
        
        ##### NOTE: ALL CODE REMOVED FOR HERE FOR PARALLELIZATION TEST #####
        ##### MODEL 2 - REMOVES THE ENVIRONMENTAL LINKAGE #####
        # Execute the MSE loop for one realization
        # mod2_flag <- FALSE
        # mod2 <- loop_through_fn(
        #   om = om_with_data,
        #   em_info = info,
        #   random = random,
        #   sel_em = sel,
        #   M_em = M,
        #   # ecov_em = ecov_em,
        #   NAA_re_em = NAA_re_em,
        #   em.opt = list(
        #     separate.em = FALSE,
        #     separate.em.type = 3,
        #     do.move = TRUE,
        #     est.move = TRUE
        #   ),
        #   update_catch_info = list(agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
        #                            catch_Neff = input_Ecov$data$catch_Neff),
        #   update_index_info = list(agg_index_sigma = input_Ecov$data$agg_index_sigma,
        #                            index_Neff = input_Ecov$data$index_Neff),
        #   assess_years = assess.years,
        #   assess_interval = assess.interval,
        #   base_years = base.years,
        #   year.use = length(base.years),
        #   add.years = TRUE,
        #   seed = seed,
        #   hcr = hcr,
        #   save.last.em = TRUE # If True, will save all EM information from every iteration, file size can be large, but you can only plot the EM output (using plot_wham_output function) when TRUE...
        # )
        # mod2_flag <- TRUE
        # print("Done with no-link model!")
        # 
        # ##### MODEL 3 - Assumes Gaussian relationship between temperature and recruitment #####
        # # IMPORTANT: As of now (03/18/2026) we don't estimate it directly. We fix the values for the relationship (Topt, Width)
        # mod3_flag <- FALSE
        # ecov_em1 <- ecov_em
        # ecov_em1$recruitment_how[] = "none"
        # mod3 <- loop_through_fn(
        #   om = om_with_data,
        #   em_info = info,
        #   random = random,
        #   sel_em = sel,
        #   M_em = M,
        #   ecov_em = ecov_em1, # remember to set no link between ecov and rec here!
        #   NAA_re_em = NAA_re_em,
        #   gauss_rec_em = list( # This will create the linkage
        #     use = TRUE,
        #     Ecov_rec_T_col = 1,   # R index, first Ecov column
        #     Topt_rec = 0.0,
        #     width_rec = gauss_width, # Should change to gauss_width!
        #     beta_T_rec = 1,
        #     estimate = FALSE      # fixed at these values
        #     #### IMPORTANT ! ####
        #     # We shouldn't estimate this!~
        #     ####
        #   ),
        #   em.opt = list(
        #     separate.em = FALSE,
        #     separate.em.type = 3,
        #     do.move = FALSE,
        #     est.move = FALSE
        #   ),
        #   update_catch_info = list(
        #     agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
        #     catch_Neff = input_Ecov$data$catch_Neff
        #   ),
        #   update_index_info = list(
        #     agg_index_sigma = input_Ecov$data$agg_index_sigma,
        #     index_Neff = input_Ecov$data$index_Neff
        #   ),
        #   assess_years = assess.years,
        #   assess_interval = assess.interval,
        #   base_years = base.years,
        #   year.use = length(base.years),
        #   add.years = TRUE,
        #   seed = seed,
        #   hcr = hcr,
        #   save.last.em = TRUE
        # )
        # mod3_flag <- TRUE
        # print("Done with gaussian link model!")
        # 
        
        #### SAVE MODEL RUNS ####
        # FILE/FOLDER STRUCTURE FOR MODEL SAVING: models/sensitivity_analysis
        SAVE_MODEL <- TRUE
        if(SAVE_MODEL){
          # # Create a folder for saving all the data and run information
          # # We will use the date+time from the start time of the code
          # # Format run_start_time as a posixDate object
          # folder_name <- format(run_start_time, "%Y-%m-%d_%H-%M-%S")
          # folder_path <- here("models","sensitivity_analysis",folder_name,"models")
          # # Create a folder. Suppress warnings and allow recursive folders to be created
          # dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
          saveRDS(mod1, here(folder_path,paste("sens_run_",run_id,"iter_id_",iter_id,"mod_1",".RDS",sep="")))
          saveRDS(mod2, here(folder_path,paste("sens_run_",run_id,"iter_id_",iter_id,"mod_2",".RDS",sep="")))
          saveRDS(mod3, here(folder_path,paste("sens_run_",run_id,"iter_id_",iter_id,"mod_3",".RDS",sep="")))
          # for(iter in seq(iterations)){
          #   saveRDS(model_list[iter], here(folder_path,paste("model_run_",iter,".RDS",sep="")))
          # }
          # Write in model results
          write_csv(em_ssb_dif_table, here(folder_path_diagnostics,paste("diagnostics_ssb_diff_",run_id,"iter_id",iter_id,".csv",sep="")))
          write_csv(rec_par_df, here(folder_path_diagnostics,paste("diagnostics_rec_par_diff_",run_id,"iter_id",iter_id,".csv",sep="")))
          write_csv(ecov_beta_df, here(folder_path_diagnostics,paste("ecov_beta_diff_",run_id,"iter_id",iter_id,".csv",sep="")))
          print(paste("Models saved for Run",run_id,"and seed: ",seed,sep=" "))
          beepr::beep(3)
        }
      },
      error=function(e){
        if(mod1_flag==FALSE){
          print("Model 1 failed to converge")
        } else if(mod2_flag==FALSE){
          print("Model 2 failed to converge")
        } else if(mod3_flag==FALSE){
          print("Model 3 failed to converge")
        } else{
          print("Some other error!")
        }
        # Try to save what we can of these models
        SAVE_MODEL <- TRUE
        if(SAVE_MODEL){
          # # Create a folder for saving all the data and run information
          # # We will use the date+time from the start time of the code
          # # Format run_start_time as a posixDate object
          # folder_name <- format(run_start_time, "%Y-%m-%d_%H-%M-%S")
          # folder_path <- here("models","sensitivity_analysis",folder_name,"models")
          # # Create a folder. Suppress warnings and allow recursive folders to be created
          # dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
          saveRDS(mod1, here(folder_path,paste("sens_run_",run_id,"_iter_id_",iter_id,"_mod_1",".RDS",sep="")))
          saveRDS(mod2, here(folder_path,paste("sens_run_",run_id,"_iter_id_",iter_id,"_mod_2",".RDS",sep="")))
          saveRDS(mod3, here(folder_path,paste("sens_run_",run_id,"_iter_id_",iter_id,"_mod_3",".RDS",sep="")))
          print(paste("Models saved for Run",run_id,"and seed: ",seed,sep=" "))
          print("IMPORTANT: This run ran into a critical error, so watchout for incomplete model objects!")
          beepr::beep(3)
        }
        skip_to_next <<- TRUE
      })
      if(skip_to_next){
        next
      }
    } # FOREACH LOOP - PARALLELIZATION ENDS
    print(paste("All seed runs ended for sensitivity analysis setting: ",run_id))
  } # FOR LOOP FOR SENSITIVITY SETTING ENDS
} # IF ENDS

model_csv <- c("mod_1", "mod_2", "mod_3")
sens_analysis_to_csv <- tidyr::crossing(proc_error= proc_error_v,
                                        mse_gap = mse_gaps_v,
                                        gauss_width=gauss_width_v,
                                        model=model_csv)

sens_analysis_to_csv <- sens_analysis_to_csv %>% mutate(run_id=rep(1:total_comb_no, each=3), .before=1)
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

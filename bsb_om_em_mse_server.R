#
#           Management Strategy Evaluation (MSE) Example for
#           Black Sea Bass (BSB) with Environmental Drivers
#                   BSB MSE - Environmental Drivers
#
#
# Author(s): Chengxue Li [Original], Jeewantha Bandara [This version]
# Date: 2025/10/25
# Runtime environment: MacOS Sequoia 15.5 on M-chip Macbook Pro, R version 4.4.1

library(wham)
library(whamMSE)
library(dplyr)
library(here)

# NOTE: This code is set to run on a server as well as a personal computer

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
iterations <- 2
base_random_seed <- 853
set.seed(base_random_seed)
mse_random_seeds <- as.integer(floor(runif(iterations, min=0, max=1000)))
model_name <- "BSB Ecov"

#### EXTERNAL FUNCTIONS ####
# Source reusable functions from `functions/reusable_functions.R`
source(here("functions","reusable_functions.R"))

# suppose we have feedback years = 3
n_feedback_years = 3

#### READ IN SAVED OPERATING MODEL OBJECT ####
OMa <- readRDS(here("models","OM_base.RDS"))
# Read in the asap models
asap <- wham::read_asap3_dat(here("data",c("north.dat","south.dat")))
temp <- wham::prepare_wham_input(asap)

#### SET ECOVARIATE DATA ####
# 1. Read in bottom temperatures for north and south
# 2. Mean and center the values
# 3. 
north_bt <- read.csv(here("data","bsb_bt_temp_nmab_1959-2022.csv"))
south_bt <- read.csv(here("data","bsb_bt_temp_smab_1959-2022.csv"))
ecov <- list(label = c("North_BT","South_BT"))
ecov$mean <- cbind(north_bt[,'mean'], south_bt[,'mean'])
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean,2,mean))
ecov$mean <- rbind(ecov$mean,matrix(0, n_feedback_years, ncol(ecov$mean)))
# Now we assume the future 3 years temp will be the same the temp in 2022
plot(ecov$mean[,1], type = "l") # remember this is deviations not the actual values
plot(ecov$mean[,2], type = "l") # remember this is deviations not the actual values

# You have 3 options to model obs error of Ecov time series
# 1. fix obs. error to be close to 0
ecov$logsigma <- 'est_1'  # fix obs. error to be close to 0

# 2. Use original obs error from the database
# ecov$logsigma <- log(cbind(north_bt[, "se"], south_bt[, "se"]))
# ecov$logsigma <- rbind(ecov$logsigma, ecov$logsigma[rep(33, 3), , drop = FALSE])

# 3. Use new but relatively low and time-invariant obs error
# ecov$logsigma <- matrix(log(0.2), 67, 2)

ecov$year <- c(north_bt[,'year'],2023:2025)
ecov$use_obs <- matrix(1, NROW(ecov$mean),NCOL(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- apply(ecov$mean, 2, mean) # Mean them by the column (the two subunits)
ecov$recruitment_how <- matrix(c("controlling-lag-0-linear","none","none","controlling-lag-0-linear"), 2,2)

#### BASIC INPUT DATA ####
# Now Let's work on data input
n_stocks = 2
n_regions = 2
n_ages = 8

# Define MSE timeframe
year_start <- 1989
year_end   <- 2021
MSE_years  <- 3 # This will be 2023,2024,2025
hist_years = length(year_start:year_end)
total_years = length(year_start:(year_end+3))
user_maturity <- array(NA, dim = c(n_stocks,total_years,n_ages))
user_maturity[,1:hist_years,] <- OMa$input$data$mature
for (i in (hist_years+1):(hist_years+3)) user_maturity[,i,] <- OMa$input$data$mature[, 33,, drop=FALSE] # Assume MAA for the next 3 years is the same as MAA for the terminal year (i.e, 2021)

# Check WAA Pointers 
# Trying to figure out what the `10` here represents
# The commercial and recreational fleets for north and south are 1,2,3,4
# The VAST and recreational CPA surveys for north and south are 5,6,7,8
# Spawning stock biomass WAA are 9,10
dim(OMa$input$data$waa) # 10 33  8
OMa$input$data$waa_pointer_fleets
OMa$input$data$waa_pointer_indices
OMa$input$data$waa_pointer_ssb
OMa$input$data$waa_pointer_M = OMa$input$data$waa_pointer_ssb

# Set WAA for projection years. We use `user_waa` for this purpose
user_waa <- list()
user_waa$waa <- array(NA, dim = c(10,33+n_feedback_years,8))
user_waa$waa[,1:33,] <- OMa$input$data$waa
for (i in 34:36) user_waa$waa[,i,] <- OMa$input$data$waa[, 33,, drop=FALSE] # Assume WAA for the next 3 years is the same as WAA for the terminal year (i.e. 2021)
user_waa$waa_pointer_fleets <- OMa$input$data$waa_pointer_fleets
user_waa$waa_pointer_indices <- OMa$input$data$waa_pointer_indices
user_waa$waa_pointer_totcatch <- OMa$input$data$waa_pointer_ssb
user_waa$waa_pointer_ssb <- OMa$input$data$waa_pointer_ssb
user_waa$waa_pointer_M <- OMa$input$data$waa_pointer_M

# Read in when they spawn from the `asap` model
fracyr_spawn = asap[[1]]$dat$fracyr_spawn

# Set catch information for the estimation model
catch_info <- list(
  catch_cv = c(asap[[1]]$dat$catch_cv[1,],asap[[2]]$dat$catch_cv[1,]), # just use the value from first year, we can make the CV to be yearly-varying later on 
  catch_Neff = c(asap[[1]]$dat$catch_Neff[1,],asap[[2]]$dat$catch_Neff[1,]), # just use the value from first year, we can make the CV to be yearly-varying later on 
  use_agg_catch = 1, # if 1, mean use annual total catch data across all years, will apply to all fleets 
  use_catch_paa = 1 # if 1, mean use annual proportion at age catch data across all years, will apply to all fleets 
)

# Need to back calculate fraction_index:
# This is an issue with how fractional years/indices are calculated between WHAM and ASAP
fracyr_seasons <- OMa$input$data$fracyr_seasons
fracyr_seasons
fracyr_indices <- OMa$input$data$fracyr_indices[1,]
fracyr_indices

fracyr_seasons <- OMa$input$data$fracyr_seasons
fracyr_indices <- OMa$input$data$fracyr_indices[1,]
index_seasons <- OMa$input$data$index_seasons

# Call in the function from `reusable_functions.R`
out <- reconcile_index_timing(fracyr_seasons, fracyr_indices, index_seasons,
                              boundary_rule = "next_season")

# This is a sanity check
# Sanity check: equals fracyr_indices (≈ within 1e-8)
fracyr_indices = out$fracyr_rebuilt

# Create a list of survey indices information
index_info <- list(
  index_cv = rep(0.5,4), # just make up some values, we will revisit and modify later
  index_Neff = rep(25,4), # just make up some values, we will revisit and modify later
  fracyr_indices = fracyr_indices, 
  q = OMa$rep$q[1,], # q should be estiamted from OMa; q is survey catchability
  use_indices = rep(1,4), # Assume I use all annual total index data 
  use_index_paa = rep(1,4), # Assume I use all annual proportion at age index data 
  units_indices = rep(2,4), 
  units_index_paa = rep(2,4) 
)

#### GENERATE OPERATING MODEL  ####
# Generate basic information for the WHAM OM
info <- whamMSE::generate_basic_info(n_stocks = 2, n_regions = 2, n_indices = 4, n_fleets = 4, n_seasons = 11,
                            base.years = year_start:year_end, n_feedback_years = MSE_years, n_ages = 8,
                            catch_info = catch_info, index_info = index_info, user_waa = user_waa, 
                            user_maturity = user_maturity, fracyr_spawn = fracyr_spawn, fracyr_seasons = fracyr_seasons)

# Get the basic information to generate the operating and estimation model
basic_info <- info$basic_info
catch_info_use <- info$catch_info
index_info_use <- info$index_info
# For 'F', we need to get the actual estimated values from the fitted model
F_info <- info$F 
F_info$F[1:33,] <- OMa$rep$Fbar[,1:4] # we must use newly estimated F values from OMa!

# Selectivity
# For selectivity, we are going to set some basic values and fix certain parameters so they
# won't be estimated by the model
sel <- list(n_selblocks = 8,
            model = rep(c("age-specific","logistic","age-specific"),
                        c(2,2,4)))
sel$initial_pars <- list(
  rep(c(0.5,1),c(3,5)), #north comm
  rep(c(0.5,1),c(6,2)), #north rec
  c(5,1), #south comm - South commercial and recreational are logistic selectivity
  c(5,1),	#south rec
  rep(c(0.5,1,1),c(1,1,6)), #north rec cpa 
  rep(c(0.5,1),c(4,4)), #north vast
  rep(c(0.5,1,1),c(2,4,2)), #south rec cpa
  rep(c(0.5,1),c(1,7)) #south vast
)
sel$fix_pars <- list(
  4:8, #north comm
  7:8, #north rec
  NULL, #south comm
  NULL, #south rec
  2:8, #north rec cpa
  5:8, #north vast
  3:8, #south rec cpa
  2:8 #south vast
)
# There's no time varying selectivity. Thereby no random effects
sel$re = rep("none",8) # assume no time varying selectivity

# Catchability - 
# double check q
OMa$parList$logit_q
OMa$input$data$q_lower # 0
OMa$input$data$q_upper # 1000
gen.invlogit(OMa$parList$logit_q,low = 0, upp = 1000)

# Now we are going to use the selectivity values from the fitted model and use them
# as initial values for the selectivity parameters
# selectivity used estimated Sel vals from OMa!
sel_vals = gen.invlogit(OMa$parList$logit_selpars,low = 0, upp = 1)
sel_vals
sel$initial_pars <- list(
  sel_vals[1,1:8], #north comm
  sel_vals[2,1:8], #north rec
  sel_vals[3,9:10], #south comm - Because of logistic selectivity
  sel_vals[4,9:10],	#south rec - Because of logistic selectivity
  sel_vals[5,1:8], #north rec cpa
  sel_vals[6,1:8], #north vast
  sel_vals[7,1:8], #south rec cpa
  sel_vals[8,1:8] #south vast
)

# Set mortality for the model
M <- list(model = "constant", initial_means = array(0.4, dim = c(n_stocks, n_regions, n_ages)))
# Assume constant M here (But be careful if you have M re later)

# Set process errors for Numbers at age (NAA)
# Use the process errors from the fitted model
sigma_vals <- array(NA, dim = c(n_stocks, n_regions, n_ages))
sigma_vals <- exp(OMa$parList$log_NAA_sigma)

# This is the NAA random effects
NAA_re = list(recruit_model = 2,
              recruit_pars = list(exp(9.032686), exp(9.753513)), # You need to convert to natural scale
              sigma_vals = sigma_vals,
              sigma = list("rec+1","rec+1"), 
              cor = list("2dar1","2dar1"), # I use 2dar1, it seems the EM is working well
              N1_model = rep("age-specific-fe",2))

# Below NAA re config are from previous code
# set fixed values for NAA re for stock 1 in south on Jan 1
# This is because of natal homing
NAA_re$sigma_vals <- array(1, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_vals[1,2,2:temp$data$n_ages] <- 0.05 #2+ fixed to SCAA-ish for North fish in south on Jan 1.
NAA_re$sigma_map <- array(NA, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_map[1,1:2,1] <- 1
NAA_re$sigma_map[2,2,1] <- 3
NAA_re$sigma_map[1,1,2:temp$data$n_ages] <- 2
NAA_re$sigma_map[2,2,2:temp$data$n_ages] <- 4

#turn off estimation of AR1 cor parameters for north population in the south on Jan 1
# Because there's no recruitment that happens for the north popuations when they are in the southern region
x <- array(NA, dim = dim(temp$par$trans_NAA_rho))
x[1,1,1:2] <- 1:2
x[2,,1] <- 3
x[2,,2] <- 4
NAA_re$cor_map <- x

x <- array(NA, dim = dim(temp$par$trans_NAA_rho))
x[1,1,1:3] <- 1:3
x[2,2,1:3] <- 4:6
NAA_re$cor_map <- x
# input_1$map$trans_NAA_rho <- factor(x)


# Setup the basic wham model for ecovariate modeling
# We then keep adding stuff like WAA, and initial ecovariate parameters to it
input_Ecov <- wham::prepare_wham_input(basic_info = basic_info, 
                                          selectivity = sel, 
                                          M = M, NAA_re = NAA_re, ecov = ecov, 
                                          catch_info = catch_info_use, 
                                          index_info = index_info_use, 
                                          F = F_info)

# I forgot to add update waa pointer in the previous version
# IF I don't specify the pointer, then only the first WAA will be used 
# for all fleets, surveys, and stock. 
waa_info <- info$par_inputs$user_waa
input_Ecov <- update_waa(input_Ecov, waa_info = waa_info)

# Initialize true parameter values from OM
input_Ecov$par$Ecov_process_pars <- OMa$parList$Ecov_process_pars
input_Ecov$par$Ecov_beta_R <- OMa$parList$Ecov_beta_R

# Initialize numbers at age
input_Ecov$par$log_N1 <- OMa$parList$log_NAA[,,1,]
input_Ecov$par$log_N1[2,1,] = 0

# Copy the same aggrgate index process errors from the fitted model
input_Ecov$data$agg_index_sigma[1:33,] <- OMa$input$data$agg_index_sigma
input_Ecov$data$use_indices[1:33,] <- OMa$input$data$use_indices
input_Ecov$data$use_index_paa[1:33,] <- OMa$input$data$use_index_paa

# Use the last years process errors for aggregate indices for future projections
for (i in 34:36) input_Ecov$data$agg_index_sigma[i,] = OMa$input$data$agg_index_sigma[33,]

# Use the effective sampling size from the fitted model
input_Ecov$data$index_Neff[1:33,] <- OMa$input$data$index_Neff

idx1 <- which(asap[[1]]$dat$use_index == 1)
idx2 <- which(asap[[2]]$dat$use_index == 1)

# Setting up effective sample size for the future
Neff1 <- do.call(cbind, lapply(idx1, function(i) asap[[1]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
Neff2 <- do.call(cbind, lapply(idx2, function(i) asap[[2]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))

index_Neff <- cbind(Neff1, Neff2)

# append 3 new rows, each equal to row 33. This is for the future projection years
stopifnot(nrow(index_Neff) >= 33)
index_Neff <- rbind(index_Neff, index_Neff[rep(33, 3), , drop = FALSE])

# Set the effective sampling size in `input_Ecov`
input_Ecov$data$index_Neff <- index_Neff

# Must have this!: Update the model index information with this information
input_Ecov <- whamMSE::update_input_index_info(input_Ecov,
                                      agg_index_sigma = input_Ecov$data$agg_index_sigma, 
                                      index_Neff = input_Ecov$data$index_Neff)

# Catch data. Same story as the index data
# Do you want the model to know the total catch for a year?
input_Ecov$data$agg_catch_sigma[1:33,] <- OMa$input$data$agg_catch_sigma
input_Ecov$data$use_agg_catch[1:33,] <- OMa$input$data$use_agg_catch
input_Ecov$data$use_catch_paa[1:33,] <- OMa$input$data$use_catch_paa
# paa - Proportion at age

# Use the last year to set values from 2022 to 2025
for (i in 34:36) input_Ecov$data$agg_catch_sigma[i,] = OMa$input$data$agg_catch_sigma[33,]

input_Ecov$data$catch_Neff[1:33,] <- OMa$input$data$catch_Neff

Neff1 <- asap[[1]]$dat$catch_Neff
Neff2 <- asap[[2]]$dat$catch_Neff

catch_Neff <- cbind(Neff1, Neff2)

input_Ecov$data$catch_Neff[1:33,] <- catch_Neff

# append 3 new rows, each equal to row 33
stopifnot(nrow(catch_Neff) >= 33)
catch_Neff <- rbind(catch_Neff, catch_Neff[rep(33, 3), , drop = FALSE])

# Must have this!: Update `input_Ecov` from using this catch effective size information
input_Ecov <- update_input_catch_info(input_Ecov,
                                      agg_catch_sigma = input_Ecov$data$agg_catch_sigma, 
                                      catch_Neff = input_Ecov$data$catch_Neff)

# Check for correct dimensions
length(OMa$parList$Ecov_re[,1])
length(OMa$rep$Ecov_x[,1])

# Copy the ecovariate random effects from the fitted model
Ecov_re <- OMa$parList$Ecov_re

input_Ecov$par$Ecov_re[1:nrow(Ecov_re),] <- Ecov_re
input_Ecov$data$do_simulate_Ecov_re <- 0 # Turn off the random generation of random effects

# Generate the unfitted wham operating model object: `unfitted om`
unfitted_om <- fit_wham(input_Ecov, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

# Remove Ecov_re from the list of random effects TMB will estimate
if ("Ecov_re" %in% unfitted_om$input$random) {
  input_Ecov$random <- unfitted_om$input$random[unfitted_om$input$random != "Ecov_re"]
} # we want to keep Ecov time series always the same as original, otherwise the TMB will generate a "FAKE" ecov time series by a random seed

# Important! We need to plug in AR1 coefficients for NAA here!
# Copy them using the fitted model
# We need to update 2dAR1 pars for NAA 
input_Ecov$par$trans_NAA_rho = OMa$parList$trans_NAA_rho

# Prepare OM object
random <- input_Ecov$random
input_Ecov$random <- NULL
# Fit the operating model for ecovariate MSE after removing the estimation for ecovariate random effects
# IMPORTANT STEP HERE
om_ecov <- fit_wham(input_Ecov, do.fit = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)
saveRDS(om_ecov, file = "om_ecov.rds")

# Update the fishing mortality of the operating model
# IMPORTANT: THIS IS THE FINAL OPERATING MODEL WE WILL FEED INTO THE MSE
om_with_data <- update_om_fn(om_ecov, seed = 123, random = random)

# Check SSB and Ecov in one realization
data = om_with_data
plot(data$rep$SSB[,1],type = "l")
plot(data$rep$SSB[,2],type = "l")

plot(data$rep$Ecov_x[,1],type = "l")
plot(data$rep$Ecov_x[,2],type = "l")

plot(ecov$mean[,1],type = "l")
lines(data$input$data$Ecov_obs[,1],col = "red")

plot(ecov$mean[,2],type = "l")
lines(data$input$data$Ecov_obs[,2],col = "red")

#### ESTIMATION MODEL ####

# Set the movement model for the estimation model
seasons = c(rep(1,5),2,rep(1,5))/12
move = list(stock_move = c(TRUE,FALSE), separable = TRUE) #north moves, south doesn't
move$must_move = array(0,dim = c(2,length(seasons),2))	
#if north stock in region 2 (south) must move back to region 1 (north) at the end of interval 5 right before spawning
move$must_move[1,5,2] <- 1 
move$can_move = array(0, dim = c(2,length(seasons),2,2))
move$can_move[1,c(1:4),2,1] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,c(7:11),1,2] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,5,2,] <- 1 #north stock can (and must) move in last season prior to spawning back to north 
mus <- array(0, dim = c(2,length(seasons),2,1))
mus[1,1:11,1,1] <- 0.02214863 #see here("2023.RT.Runs","transform_SS_move_rates.R") for how these numbers are derived.
mus[1,1:11,2,1] <- 0.3130358
move$mean_vals <- mus # Set prior values for movement
move$mean_model = matrix("stock_constant", 2,1)
#prior distribution on movement parameters 
move$use_prior <- array(0, dim = c(2,length(seasons),2,1))
move$use_prior[1,1,1,1] <- 1
move$use_prior[1,1,2,1] <- 1
move$prior_sigma <- array(0, dim = c(2,length(seasons),2,1))
move$prior_sigma[1,1,1,1] <- 0.2
move$prior_sigma[1,1,2,1] <- 0.2

# Configure the Estimation Model (EM) for use within the loop
# Add movement prior!
move_em <- move
n_seasons <- length(fracyr_seasons)
move_em$use_prior <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))
move_em$use_prior[1, 1, , ] <- 1
move_em$prior_sigma <- array(0.2, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))

n_seasons = 11
# Specify the Harvest Control Rule (HCR)
hcr <- list()
hcr$hcr.type <- 1 # FXSPR - Fishing pressure to keep the SPR at a certain percentage
hcr$hcr.opts <- list(use_FXSPR = TRUE, percentFXSPR = 75) # Apply F at 75% unfished SPR

assess.interval <- 3 # Assessments occur every 3 years
base.years <- year_start:year_end
terminal.year <- tail(base.years, 1)
last.year <- max(data$years) # 2024 is the last year in feedback
assess.years <- seq(terminal.year, last.year - assess.interval, by = assess.interval)

# Remember to use this config for EM NAAre
NAA_re_em = NAA_re
NAA_re_em$N1_model[] = "equilibrium"

# You can use 'est_1' to let EM estimate obs error for Ecov
# Here we are going to use the same environment link for the estimation model as the opearting model
# We will also keep the same observation error for the ecovariate
ecov_em <- ecov
ecov_em$logsigma <- 'est_1' 

#### RUN MSE ####

model_list <- list()

for(iter in seq(iterations)){
  print(iter)
  # Execute the MSE loop for one realization
  mod <- loop_through_fn(
    om = om_with_data,
    em_info = info,
    random = random,
    sel_em = sel,
    M_em = M,
    NAA_re_em = NAA_re_em,
    move_em = move_em,
    ecov_em = ecov_em,
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
  # Give the model a name
  mod$om$model_name <- paste(model_name, "Run", iter, sep=" ")
  # Add to the series of models
  model_list <- append(model_list, list(mod))
}


#### SAVE MODEL ####
SAVE_MODEL <- TRUE
if(SAVE_MODEL){
  # Create a folder for saving all the data and run information
  # We will use the date+time from the start time of the code
  # Format run_start_time as a posixDate object
  folder_name <- format(run_start_time, "%Y-%m-%d_%H-%M-%S")
  folder_path <- here("models","model_runs",folder_name)
  # Create a folder. Suppress warnings and allow recursive folders to be created
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  for(iter in seq(iterations)){
    saveRDS(model_list[iter], here(folder_path,paste("model_run_",iter,".RDS",sep="")))
  }
  print("Models saved!")
}

#### PLOT DIAGNOSTICS ####

##### Model names #####
model_names <- lapply(model_list, function(model) model$om$model_name)

##### Ecov #####
# Ecov observed in EM
ecov_em_obs <- lapply(model_list, function(model) model$em_full[[1]]$input$data$Ecov_obs)
# Ecov true in EM
ecov_em_true <- lapply(model_list, function(model) model$em_full[[1]]$rep$Ecov_x)

# Ecov observed in OM
ecov_om_obs <- lapply(model_list, function(model) model$om$input$data$Ecov_obs)
# Ecov true in OM
ecov_om_true <- lapply(model_list, function(model) model$om$rep$Ecov_x)


# Plot these nicely in ggplot
ecov_om_true_tidy <- data.table::as.data.table(ecov_om_true)

plot(mod$em_full[[1]]$input$data$Ecov_obs, type = "l", col = "orange",
     ylab = "Ecov", xlab = "Year", lwd = 2,
     main = "EM observation (input)")


# saveRDS(mod, here("models","whamMSE_models","mod_2025_11_20.RDS"))

# SSB in OM
plot(mod$om$rep$SSB[,1], type = "l")
plot(mod$om$rep$SSB[,2], type = "l")

# Catch in OM
plot(mod$om$rep$pred_catch[,1], type = "l")
plot(mod$om$rep$pred_catch[,2], type = "l")

mod$converge_list # = 2 means converged! 

plot_wham_output(mod$em_full[[1]], out.type = "html")


# Plotting diagnostics
# Copied from https://github.com/lichengxue/whamMSE/blob/UMassD-MSE/GBK_Example/Update_MSE_code.R#L474


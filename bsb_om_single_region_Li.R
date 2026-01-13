# -----------------------------------------------------------------------------
#
#                         OPERATING MODEL FOR 
#       BLACK SEA BASS WITH ENVIRONMENTAL COVARIATES FOR SINGLE REGION
#
#
# Author(s): Chengxue Li [Original], Jeewantha Bandara [This version]
# Date: 2026/01/13
# Runtime environment: MacOS Sequoia 15.5 on M-chip Macbook Pro, R version 4.4.1
#
#
#
# -----------------------------------------------------------------------------


library(wham)
library(here)
setwd("C:/Users/liche/Desktop/Rutgers-MSE")

#### READ IN ASAP MODEL ####
asap <- wham::read_asap3_dat("data/north.dat")
temp <- wham::prepare_wham_input(asap)

#### BASIC INFO ####

basic_info <- list(region_names = "North", 
                   stock_names = paste0("BSB_", "North")) 

seasons = c(rep(1,5),2,rep(1,5))/12
basic_info$fracyr_seasons <- seasons

#average recruitment over years 2000+ for SSB40 BRPs
basic_info$XSPR_R_avg_yrs <- which(temp$years>1999)
basic_info$XSPR_R_opt <- 2 #use average of recruitments (random effects), not expected/predicted given last time step

#### ENVIRONMENTAL COVARIATE ####
north_bt <- read.csv("data/bsb_bt_temp_nmab_1959-2022.csv")

ecov <- list(label = c("North_BT"))
ecov$mean <- cbind(north_bt[,'mean'])
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean,2,mean))
ecov$logsigma <- log(north_bt[,'se'])
ecov$year <- north_bt[,'year']
ecov$use_obs <- matrix(1, NROW(ecov$mean),NCOL(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- apply(ecov$mean, 2, mean)

ecov_0 <- ecov 
ecov_0$recruitment_how <- matrix("none", 1,1)
ecov_1 <- ecov_2 <- ecov_0
ecov_1$recruitment_how[1,1] <- "controlling-lag-0-linear" #north

# Check ecov_3 again
ecov_3$recruitment_how

# NAA_re - A state-space model and the random effects are 2nd degree 
# (both the year and age/parameter) first order auto-correlation (2dar1)
NAA_re = list(sigma = list("rec+1"), cor = list("2dar1"), 
              N1_model = rep("equilibrium",1))

#Wasn't decoupled in RT
#NAA_re$decouple_recruitment <- FALSE

# set fixed values for NAA re for stock 1 in south on Jan 1
NAA_re$sigma_vals <- array(1, dim = c(1,1,temp$data$n_ages))
NAA_re$sigma_vals[1,1,1] <- 0.5 # Rec sigma
NAA_re$sigma_vals[1,1,2:temp$data$n_ages] <- 0.2 # NAA sigma

#### SELECTIVITY ####

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

#### CATCH INFORMATION ####
# Effective sampling size
catch_Neff <- temp$data$catch_Neff
catch_Neff[] <- 1000
catch_info <- list(catch_Neff = catch_Neff)
x <- temp$data$selblock_pointer_fleets
x[] <- rep(1:2, each = NROW(x))
catch_info$selblock_pointer_fleets = x

#### SURVEY INDEX INFORMATION ####
index_Neff <- temp$data$index_Neff
index_Neff[] <- 1000
index_info <- list(index_Neff = index_Neff)
x <- temp$data$selblock_pointer_indices
x[] <- rep(3:4, each = NROW(x))
index_info$selblock_pointer_indices = x

#############################################
# M random effects just on age 1
# Note: I don't think this is used anywhere in the model
M_re_opt    = list(M_model = "constant", M_mean = 0.4)

#### Some candidate models ####
# NOTE: WE ARE SAVING THE ecov_3 model as the base model
##### Ecov-1 model #####
input_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, 
                              basic_info = basic_info, M = M_re_opt,
                              ecov = ecov_1, catch_info = catch_info,
                              index_info = index_info)

input_1$fleet_names <- paste0(rep(c("North_"),each = 1), temp$fleet_names)
input_1$index_names <- paste0(rep(c("North_"),each = 1), temp$index_names)

fit_1 <- fit_wham(input_1, do.fit = FALSE,
                  do.brps = FALSE, do.sdrep = TRUE, 
                  do.osa = FALSE, do.retro = FALSE)

temp <- fit_1$input
temp$par <- fit_1$parList
temp$par$trans_NAA_rho


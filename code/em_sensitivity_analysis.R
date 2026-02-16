# -----------------------------------------------------------------------------
#     BSB MSE with Environmental Drivers - Gaussian Recruitment Test
# -----------------------------------------------------------------------------
# install.packages("remotes")
# remotes::install_github("lichengxue/wham@Gussian_Rec")

# setwd("C:/Users/liche/Desktop/Rutgers-MSE")
library(wham)
library(whamMSE)
library(dplyr)
library(here)

here::here()

# ----------------------------------
# Basic setup / read in objects
# ----------------------------------
n_feedback_years <- 15

OMa  <- readRDS("models/OM_base.RDS")
asap <- read_asap3_dat("data/north.dat")

# ----------------------------------
# 1. Ecov data
# ----------------------------------
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

# ----------------------------------
# 2. Biological / MSE time frame
# ----------------------------------
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

# ----------------------------------
# 3. Catch info
# ----------------------------------
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

# ----------------------------------
# 5. Generate basic_info (whamMSE)
# ----------------------------------
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

# ----------------------------------
# 6. Selectivity / M / NAA_re
# ----------------------------------

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

# lower sigma for recruitment and NAA t0 0.005
vals[] = 0.1

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

# ----------------------------------
# 7. Build input_Ecov
# ----------------------------------
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

# ----------------------------------
# 8. Gaussian T–recruit settings
# ----------------------------------
input_Ecov$data$use_gauss_T_rec <- 1L         # turn ON Gaussian link
input_Ecov$data$Ecov_rec_T_col  <- 0L         # first Ecov column = North_BT

temp_col <- input_Ecov$data$Ecov_rec_T_col + 1L
temp_vec <- input_Ecov$data$Ecov_obs[, temp_col]

# WE need to specify clearly
input_Ecov$par$Topt_rec      <- 0.0          # peak at 0
input_Ecov$par$log_width_rec <- log(1.0)     # width = 1
input_Ecov$par$log_width_rec <- log(0.5)     # width = 1
n_stocks <- input_Ecov$data$n_stocks
input_Ecov$par$beta_T_rec    <- rep(1, n_stocks)

# ---- FIX map for Gaussian T-recruit params ----
for(nm in c("Topt_rec","log_width_rec","beta_T_rec")) {
  if(!is.null(input_Ecov$par[[nm]])) {
    input_Ecov$map[[nm]] <- factor(rep(NA, length(input_Ecov$par[[nm]])))
  }
}

# ----------------------------------
# 9. Fix N1, catch & index info
# ----------------------------------
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

# ----------------------------------
# 10. Force Ecov_re pattern: -2 to +2, 1989–2023
# ----------------------------------
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

# ----------------------------------
# 11. Build OM, remove Ecov_re from random, plug NAA rho
# ----------------------------------
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
saveRDS(om_ecov, file = "om_ecov.rds")

om_with_data <- update_om_fn(om_ecov, seed = 123, random = random)

om_with_data$input$data$agg_catch #simulated catch
om_with_data$input$data$agg_indices #simulated index

# ----------------------------------
# 12. Diagnostics: show Gaussian T–R relationship
# ----------------------------------
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



# Specify the Harvest Control Rule (HCR)
hcr <- list()
hcr$hcr.type <- 1 # FXSPR - Fishing pressure to keep the SPR at a certain percentage
hcr$hcr.opts <- list(use_FXSPR = TRUE, percentFXSPR = 75) # Apply F at 75% unfished SPR

assess.interval <- 3 # Assessments occur every 3 years
base.years <- year_start:year_end
terminal.year <- tail(base.years, 1)
last.year <- 2024+12
assess.years <- seq(terminal.year, last.year - assess.interval, by = assess.interval)

# Remember to use this config for EM NAAre
NAA_re_em = NAA_re
NAA_re_em$N1_model[] = "equilibrium"

# You can use 'est_1' to let EM estimate obs error for Ecov
ecov_em <- ecov
ecov_em$logsigma <- 'est_1'


# Execute the MSE loop for one realization
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

plot(mod1$om$rep$SSB, type = "l", col = "red")
lines(mod2$om$rep$SSB, type = "l", col = "blue")

plot(mod1$om$rep$pred_catch[,1], type = "l", col = "red")
lines(mod2$om$rep$pred_catch[,1], type = "l", col = "blue")

plot(mod1$om$rep$pred_catch[,2], type = "l", col = "red")
lines(mod2$om$rep$pred_catch[,2], type = "l", col = "blue")

lines(mod1$em_full[[1]]$rep$SSB, col = "red")
lines(mod2$em_full[[1]]$rep$SSB, col = "blue")
lines(mod3$em_full[[1]]$rep$SSB, col = "purple")

mod1$em_full[[1]]$rep$SSB - mod2$em_full[[1]]$rep$SSB

mod2$em_full[[1]]$parList$mean_rec_pars - mod1$em_full[[1]]$parList$mean_rec_pars
mod2$em_full[[1]]$parList$Ecov_beta_R - mod1$em_full[[1]]$parList$Ecov_beta_R

plot(mod$om$rep$NAA[,,,1], type = "l")
lines(mod$em_full[[1]]$rep$NAA[,,,1], col = "red")

plot(mod$om$rep$Fbar[,1], type = "l")
lines(mod$em_full[[1]]$rep$Fbar[,1], col = "red")

plot(mod$om$rep$pred_catch[,1], type = "l")
lines(mod$em_full[[1]]$rep$pred_catch[,1], col = "red")

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

mod$om$parList$trans_NAA_rho[,,1]
mod1$em_full[[1]]$parList$trans_NAA_rho
mod2$em_full[[1]]$parList$trans_NAA_rho
mod3$em_full[[1]]$parList$trans_NAA_rho


mod1$em_full[[1]]$input$data$Ecov_how_R
mod2$em_full[[1]]$input$data$Ecov_how_R
mod3$em_full[[1]]$input$data$Ecov_how_R


mod1$em_full[[1]]$rep$nll
mod2$em_full[[1]]$rep$nll
mod1$em_full[[1]]$rep$nll

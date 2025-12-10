# -----------------------------------------------------------------------------
#     BSB MSE with Environmental Drivers - Gaussian Recruitment Test
# -----------------------------------------------------------------------------

library(wham)
library(whamMSE)
library(dplyr)
library(here)

# ----------------------------------
# Basic setup / read in objects
# ----------------------------------
n_feedback_years <- 3

OMa  <- readRDS(here("models","OM_base.RDS"))
asap <- wham::read_asap3_dat(here("data",c("north.dat","south.dat")))
temp <- wham::prepare_wham_input(asap)

# ----------------------------------
# 1. Ecov data
# ----------------------------------
north_bt <- read.csv(here("data","bsb_bt_temp_nmab_1959-2022.csv"))
south_bt <- read.csv(here("data","bsb_bt_temp_smab_1959-2022.csv"))

ecov <- list(label = c("North_BT","South_BT"))

# original means
ecov$mean <- cbind(north_bt[,"mean"], south_bt[,"mean"])

# center each column
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean, 2, mean))

# extend to include 3 feedback years (2023–2025)
ecov$mean <- rbind(ecov$mean, matrix(0, n_feedback_years, ncol(ecov$mean)))

# for this test, force north BT deviations to 0 so all variation comes from Ecov_re
ecov$mean[,1] <- 0

ecov$year <- c(north_bt[,"year"], 2023:2025)

plot(ecov$year, ecov$mean[,1], type = "l", main = "Ecov mean (North BT)", xlab = "Year")

# obs error for Ecov
ecov$logsigma           <- "est_1"
ecov$use_obs            <- matrix(1, nrow(ecov$mean), ncol(ecov$mean))
ecov$process_model      <- "ar1"
ecov$process_mean_vals  <- apply(ecov$mean, 2, mean)

# keep linear Ecov–R links defined, but we will zero out Ecov_beta_R later
ecov$recruitment_how <- matrix(
  c("controlling-lag-0-linear","none",
    "none","controlling-lag-0-linear"),
  nrow = 2, byrow = TRUE
)

# ----------------------------------
# 2. Biological / MSE time frame
# ----------------------------------
n_stocks  <- 2
n_regions <- 2
n_ages    <- 8

year_start <- 1989
year_end   <- 2021
MSE_years  <- 3

hist_years <- length(year_start:year_end)
total_years <- length(year_start:(year_end + MSE_years))

# maturity
user_maturity <- array(NA, dim = c(n_stocks, total_years, n_ages))
user_maturity[, 1:hist_years, ] <- OMa$input$data$mature
for (i in (hist_years+1):(hist_years+MSE_years)) {
  user_maturity[, i, ] <- OMa$input$data$mature[, 33, , drop = FALSE]
}

# WAA and pointers
OMa$input$data$waa_pointer_M <- OMa$input$data$waa_pointer_ssb

user_waa <- list()
user_waa$waa <- array(NA, dim = c(10, 33 + n_feedback_years, 8))
user_waa$waa[, 1:33, ] <- OMa$input$data$waa
for (i in 34:36) {
  user_waa$waa[, i, ] <- OMa$input$data$waa[, 33, , drop = FALSE]
}
user_waa$waa_pointer_fleets   <- OMa$input$data$waa_pointer_fleets
user_waa$waa_pointer_indices  <- OMa$input$data$waa_pointer_indices
user_waa$waa_pointer_totcatch <- OMa$input$data$waa_pointer_ssb
user_waa$waa_pointer_ssb      <- OMa$input$data$waa_pointer_ssb
user_waa$waa_pointer_M        <- OMa$input$data$waa_pointer_M

fracyr_spawn <- asap[[1]]$dat$fracyr_spawn

# ----------------------------------
# 3. Catch info
# ----------------------------------
catch_info <- list(
  catch_cv       = c(asap[[1]]$dat$catch_cv[1,], asap[[2]]$dat$catch_cv[1,]),
  catch_Neff     = c(asap[[1]]$dat$catch_Neff[1,], asap[[2]]$dat$catch_Neff[1,]),
  use_agg_catch  = 1,
  use_catch_paa  = 1
)

# ----------------------------------
# 4. Index timing & info
# ----------------------------------
fracyr_seasons  <- OMa$input$data$fracyr_seasons
fracyr_indices  <- OMa$input$data$fracyr_indices[1,]
index_seasons   <- OMa$input$data$index_seasons

reconcile_index_timing <- function(fracyr_seasons, fracyr_indices, index_seasons,
                                   boundary_rule = c("next_season", "this_season")) {
  boundary_rule <- match.arg(boundary_rule)
  n_seasons <- length(fracyr_seasons)
  stopifnot(length(index_seasons) == length(fracyr_indices))
  if (abs(sum(fracyr_seasons) - 1) > 1e-8) stop("fracyr_seasons must sum to 1.")
  if (any(index_seasons < 1 | index_seasons > n_seasons)) stop("index_seasons out of range.")
  
  starts <- c(0, cumsum(fracyr_seasons))
  eps <- sqrt(.Machine$double.eps)
  
  x <- pmin(pmax(fracyr_indices, 0), 1 - eps)
  
  x_nudge <- if (boundary_rule == "next_season") {
    y <- x + eps
    y[x == 0] <- eps
    y
  } else {
    pmax(x - eps, 0)
  }
  implied_season <- findInterval(x_nudge, starts, rightmost.closed = FALSE)
  implied_season[implied_season < 1] <- 1
  implied_season[implied_season > n_seasons] <- n_seasons
  
  mismatch <- which(implied_season != index_seasons)
  if (length(mismatch)) {
    message("Note: ", length(mismatch),
            " index/indices have implied_season != index_seasons: [",
            paste(mismatch, collapse = ", "), "]. Using index_seasons as authoritative.")
  }
  
  season_len <- fracyr_seasons[index_seasons]
  offset <- x - starts[index_seasons]
  offset[abs(offset) < eps] <- 0
  offset <- pmin(pmax(offset, 0), season_len - eps)
  
  prop_within <- ifelse(season_len > 0, offset / season_len, 0)
  
  fracyr_rebuilt <- starts[index_seasons] + offset
  fracyr_rebuilt <- pmin(pmax(fracyr_rebuilt, 0), 1 - eps)
  
  list(
    starts              = starts,
    implied_season      = implied_season,
    used_season         = index_seasons,
    within_season_offset= offset,
    within_season_prop  = prop_within,
    fracyr_rebuilt      = fracyr_rebuilt
  )
}

out <- reconcile_index_timing(fracyr_seasons, fracyr_indices, index_seasons,
                              boundary_rule = "next_season")

fracyr_indices <- out$fracyr_rebuilt

index_info <- list(
  index_cv        = rep(0.5, 4),
  index_Neff      = rep(25, 4),
  fracyr_indices  = fracyr_indices,
  q               = OMa$rep$q[1,],
  use_indices     = rep(1, 4),
  use_index_paa   = rep(1, 4),
  units_indices   = rep(2, 4),
  units_index_paa = rep(2, 4)
)

# ----------------------------------
# 5. Generate basic_info (whamMSE)
# ----------------------------------
info <- whamMSE::generate_basic_info(
  n_stocks       = 2,
  n_regions      = 2,
  n_indices      = 4,
  n_fleets       = 4,
  n_seasons      = 11,
  base.years     = year_start:year_end,
  n_feedback_years = MSE_years,
  n_ages         = 8,
  catch_info     = catch_info,
  index_info     = index_info,
  user_waa       = user_waa,
  user_maturity  = user_maturity,
  fracyr_spawn   = fracyr_spawn,
  fracyr_seasons = fracyr_seasons
)

basic_info     <- info$basic_info
catch_info_use <- info$catch_info
index_info_use <- info$index_info
F_info         <- info$F
F_info$F[1:33,] <- OMa$rep$Fbar[, 1:4]

# ----------------------------------
# 6. Selectivity / M / NAA_re
# ----------------------------------
gen.invlogit <- function(eta, low, upp, s = 1) {
  low + (upp - low) * plogis(s * eta)
}

sel <- list(
  n_selblocks = 8,
  model       = rep(c("age-specific","logistic","age-specific"), c(2,2,4))
)

sel_vals <- gen.invlogit(OMa$parList$logit_selpars, low = 0, upp = 1)

sel$initial_pars <- list(
  sel_vals[1,1:8],  # north comm
  sel_vals[2,1:8],  # north rec
  sel_vals[3,9:10], # south comm
  sel_vals[4,9:10], # south rec
  sel_vals[5,1:8],  # north rec cpa
  sel_vals[6,1:8],  # north vast
  sel_vals[7,1:8],  # south rec cpa
  sel_vals[8,1:8]   # south vast
)

sel$fix_pars <- list(
  4:8,   # north comm
  7:8,   # north rec
  NULL,  # south comm
  NULL,  # south rec
  2:8,   # north rec cpa
  5:8,   # north vast
  3:8,   # south rec cpa
  2:8    # south vast
)
sel$re <- rep("none", 8)

M <- list(
  model        = "constant",
  initial_means = array(0.4, dim = c(n_stocks, n_regions, n_ages))
)

sigma_vals <- array(exp(OMa$parList$log_NAA_sigma),
                    dim = c(n_stocks, n_regions, n_ages))

NAA_re <- list(
  recruit_model = 2,
  recruit_pars  = list(exp(9.032686), exp(9.753513)), # fixed mean recruit by stock
  sigma_vals    = sigma_vals,
  sigma         = list("rec+1", "rec+1"),
  cor           = list("2dar1", "2dar1"),
  N1_model      = rep("age-specific-fe", 2)
)

# Use your previous sigma_map structure:
NAA_re$sigma_vals <- array(1, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_vals[1,2,2:temp$data$n_ages] <- 0.05
NAA_re$sigma_map <- array(NA, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_map[1,1:2,1] <- 1
NAA_re$sigma_map[2,2,1]   <- 3
NAA_re$sigma_map[1,1,2:temp$data$n_ages] <- 2
NAA_re$sigma_map[2,2,2:temp$data$n_ages] <- 4

# Turn OFF NAA random effects (so recruitment is deterministic w.r.t. T)
NAA_re$sigma_vals[] <- 0.0

x <- array(NA, dim = dim(temp$par$trans_NAA_rho))
x[1,1,1:3] <- 1:3
x[2,2,1:3] <- 4:6
NAA_re$cor_map <- x

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
input_Ecov$par$Ecov_process_pars <- OMa$parList$Ecov_process_pars
input_Ecov$par$Ecov_beta_R       <- OMa$parList$Ecov_beta_R

# For this Gaussian test, zero out the linear Ecov–R effect:
if (!is.null(input_Ecov$par$Ecov_beta_R))
  input_Ecov$par$Ecov_beta_R[] <- 0

# ----------------------------------
# 8. Gaussian T–recruit settings
# ----------------------------------
input_Ecov$data$use_gauss_T_rec <- 1L         # turn ON Gaussian link
input_Ecov$data$Ecov_rec_T_col  <- 0L         # first Ecov column = North_BT

temp_col <- input_Ecov$data$Ecov_rec_T_col + 1L
temp_vec <- input_Ecov$data$Ecov_obs[, temp_col]

input_Ecov$par$Topt_rec      <- 0.0          # peak at 0
input_Ecov$par$log_width_rec <- log(1.0)     # width = 1
n_stocks <- input_Ecov$data$n_stocks
input_Ecov$par$beta_T_rec    <- rep(1, n_stocks)

# ----------------------------------
# 9. Fix N1, catch & index info
# ----------------------------------
input_Ecov$par$log_N1 <- OMa$parList$log_NAA[,,1,]
input_Ecov$par$log_N1[2,1,] <- 0   # as you had

# index sigma and Neff
input_Ecov$data$agg_index_sigma[1:33,] <- OMa$input$data$agg_index_sigma
input_Ecov$data$use_indices[1:33,]     <- OMa$input$data$use_indices
input_Ecov$data$use_index_paa[1:33,]   <- OMa$input$data$use_index_paa

for (i in 34:36) {
  input_Ecov$data$agg_index_sigma[i,] <- OMa$input$data$agg_index_sigma[33,]
}

input_Ecov$data$index_Neff[1:33,] <- OMa$input$data$index_Neff

idx1 <- which(asap[[1]]$dat$use_index == 1)
idx2 <- which(asap[[2]]$dat$use_index == 1)

Neff1 <- do.call(cbind, lapply(idx1, function(i)
  asap[[1]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
Neff2 <- do.call(cbind, lapply(idx2, function(i)
  asap[[2]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
index_Neff <- cbind(Neff1, Neff2)
index_Neff <- rbind(index_Neff, index_Neff[rep(33,3), , drop = FALSE])
input_Ecov$data$index_Neff <- index_Neff

input_Ecov <- whamMSE::update_input_index_info(
  input_Ecov,
  agg_index_sigma = input_Ecov$data$agg_index_sigma,
  index_Neff      = input_Ecov$data$index_Neff
)

# catch sigma & Neff
input_Ecov$data$agg_catch_sigma[1:33,] <- OMa$input$data$agg_catch_sigma
input_Ecov$data$use_agg_catch[1:33,]   <- OMa$input$data$use_agg_catch
input_Ecov$data$use_catch_paa[1:33,]   <- OMa$input$data$use_catch_paa

for (i in 34:36) {
  input_Ecov$data$agg_catch_sigma[i,] <- OMa$input$data$agg_catch_sigma[33,]
}

Neff1 <- asap[[1]]$dat$catch_Neff
Neff2 <- asap[[2]]$dat$catch_Neff
catch_Neff <- cbind(Neff1, Neff2)
catch_Neff <- rbind(catch_Neff, catch_Neff[rep(33,3), , drop = FALSE])

input_Ecov$data$catch_Neff <- catch_Neff

input_Ecov <- update_input_catch_info(
  input_Ecov,
  agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
  catch_Neff      = input_Ecov$data$catch_Neff
)

# ----------------------------------
# 10. Force Ecov_re pattern: -2 to +2, 1989–2023
# ----------------------------------
Ecov_re <- OMa$parList$Ecov_re
ny      <- nrow(Ecov_re)

yrs <- ecov$year      # 1959–2025, length 67

# Year indices for 1989–2023 (should be 35 yrs)
idx <- which(yrs >= 1989 & yrs <= 2023)

# Gradient -2 → +2 over these years
grad <- seq(-2, 2, length.out = length(idx))

# Start with zeros, then fill desired segment
input_Ecov$par$Ecov_re[,] <- 0
input_Ecov$par$Ecov_re[idx, 1] <- grad    # North BT
input_Ecov$par$Ecov_re[, 2]    <- 0       # South BT flat

# Do not simulate Ecov_re – use the values above
input_Ecov$data$do_simulate_Ecov_re <- 0

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

input_Ecov$par$trans_NAA_rho <- OMa$parList$trans_NAA_rho

random <- input_Ecov$random
input_Ecov$random <- NULL

om_ecov <- fit_wham(input_Ecov, do.fit = FALSE, do.brps = TRUE,
                    MakeADFun.silent = TRUE)
saveRDS(om_ecov, file = "om_ecov.rds")

om_with_data <- update_om_fn(om_ecov, seed = 123, random = random)

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

# This script is for setting up temperature trends that we will port over to 
# code/parallelized_historic_run.R

library(dplyr)
library(here)
library(ggplot2)

here::here()

# Set the n_feedback_years
n_feedback_years <- 30

north_bt <- read.csv("data/bsb_bt_temp_nmab_1959-2022.csv")

ecov <- list(label = c("North_BT"))

# original means
ecov$mean <- cbind(north_bt[,"mean"])

# center each column
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean, 2, mean))

# extend to include the feedback years
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



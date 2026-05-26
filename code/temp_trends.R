# This script is for setting up temperature trends that we will port over to 
# code/parallelized_historic_run.R

library(dplyr)
library(here)
library(ggplot2)

here::here()

#### EXTERNAL FUNCTIONS ####
# Source reusable functions from `functions/reusable_functions.R`
source(here("functions","reusable_functions.R"))
source(here("functions","plot_themes.R"))

#### 1. ECOV DATA ####

n_feedback_years <- 30

north_bt <- read.csv("data/bsb_bt_temp_nmab_1959-2022.csv")

north_bt <- north_bt %>% select(year, mean) %>% arrange(year) # Only relevant columns
north_bt <- north_bt %>% rename(temp=mean) # Rename `mean` column
north_bt <- north_bt %>% mutate(temp_c=temp-mean(temp, na.rm=TRUE)) # Centered column
north_bt <- north_bt %>% mutate(year_n=row_number())
north_bt <- north_bt %>% mutate(period="Historical")
# Store the centering value
hist_temp_c <- mean(north_bt$temp, na.rm=TRUE)
# Last year of the historical timeseries
last_year <- north_bt %>% select(year) %>% tail(n=1) %>% pull()
# First year of the projection
first_proj_year <- last_year+1
# Last temperature
last_temp_c <- north_bt %>% select(temp_c) %>% tail(n=1) %>% pull()

#### FUTURE SCENARIOS ####

##### Last years temperature is retained for the future with some error #####
scen_0_trend <- 0.0 # 0 deg C/year
scen_0_sd <- 0.05
scen_0_error <- rnorm(n=n_feedback_years, mean=0, sd=scen_0_sd)
scen_0_proj <- scen_0_trend*seq(1,n_feedback_years) + last_temp_c  + scen_0_error

scen_0_df <- data.frame(year=seq(first_proj_year-1, first_proj_year+n_feedback_years-1), 
                        temp_c=c(last_temp_c, scen_0_proj),
                        period="Projected - 0")


##### REPORTED TREND #####
scen_1_trend <- 0.04 # 0.04 deg C/year
scen_1_sd <- 0.05
scen_1_error <- rnorm(n=n_feedback_years, mean=0, sd=scen_1_sd)
scen_1_proj <- scen_1_trend*seq(1,n_feedback_years) +last_temp_c  + scen_1_error

scen_1_df <- data.frame(year=seq(first_proj_year-1, first_proj_year+n_feedback_years-1), 
                        temp_c=c(last_temp_c, scen_1_proj),
                        period="Projected - 1")

##### TREND FROM LAST 5 YEARS #####
last_5_years <- north_bt %>% arrange(year) %>% tail(n=5)
scen_2_lm <- lm(temp_c ~ year, data=last_5_years)
scen_2_trend <- scen_2_lm$coefficients[["year"]]
scen_2_sd <- 0.5
scen_2_error <- rnorm(n=n_feedback_years, mean=0, sd=scen_2_sd)
scen_2_proj <- scen_2_trend*seq(1,n_feedback_years)+ last_temp_c + scen_2_error

scen_2_df <- data.frame(year=seq(first_proj_year-1, first_proj_year+n_feedback_years-1), 
                        temp_c=c(last_temp_c, scen_2_proj),
                        period="Projected - 2")


##### REPORTED TREND - BUT MORE STOCHASTIC #####
scen_3_trend <- 0.04 # 0.04 deg C/year
scen_3_sd <- 0.25
scen_3_error <- rnorm(n=n_feedback_years, mean=0, sd=scen_3_sd)
scen_3_proj <- scen_3_trend*seq(1,n_feedback_years) +last_temp_c  + scen_3_error

scen_3_df <- data.frame(year=seq(first_proj_year-1, first_proj_year+n_feedback_years-1), 
                        temp_c=c(last_temp_c, scen_3_proj),
                        period="Projected - 3")


# Bind all projections together
all_proj_df <- rbind(scen_0_df, scen_1_df, scen_2_df, scen_3_df)

ggplot(all_proj_df, aes(x=year, y=temp_c, color=as.factor(period))) + 
  geom_line() + 
  geom_line(data=north_bt, aes(x=year, y=temp_c), color="black") + 
  geom_vline(xintercept=last_year, color="grey", linewidth=0.5) + 
  labs(x="Year",y="Standardized temperature [deg C]", color="Temperature trend") + 
  bioshift_plot_theme_1 + theme(panel.grid.major=element_line(color="grey", linewidth=0.1))


#### HYPOTHETICAL OPTIMUMS ####

# Setting optimal temperatures
# Scenario 1 - Optimal temperature is the mean of the last 60 years
topt_1 <- 0
# Scenario 2 - Optimal temperature will be achieved in the future
topt_2 <- 3.5
# Scenario 3 - Optimal temperature will be achieved quickly in the future
topt_3 <- 2.5
# Scenario 4 - Optimal temperature was in the recent past
topt_4 <- 1.25
# All gaussian widths are set to 2 - Quite generous

# Same plot as above, but add in the different temperature optima
toptlines_df <- data.frame(
  value    = c(topt_1, topt_2, topt_3, topt_4),
  optima = c("Optimum 1", "Optimum 2", "Optimum 3", "Optimum 4")
)

ggplot(all_proj_df, aes(x=year, y=temp_c, color=as.factor(period))) + 
  geom_line() + 
  geom_line(data=north_bt, aes(x=year, y=temp_c), color="black") + 
  geom_vline(xintercept=last_year, color="grey", linewidth=0.5) + 
  geom_hline(data=toptlines_df, aes(yintercept=value, linetype=optima), alpha=0.5) + 
  labs(x="Year",y="Standardized temperature [deg C]", color="Temperature trend", linetype="Optima") + 
  bioshift_plot_theme_1 + theme(panel.grid.major=element_line(color="grey", linewidth=0.1))


scen_1_df %>% mutate(rec_response=gauss_rec(tx=temp_c, topt=topt_1, twidth=2))
scen_2_df %>% mutate(rec_response=gauss_rec(tx=temp_c, topt=topt_2, twidth=2))
scen_3_df %>% mutate(rec_response=gauss_rec(tx=temp_c, topt=topt_3, twidth=2))
scen_3_df %>% mutate(rec_response=gauss_rec(tx=temp_c, topt=topt_4, twidth=2))

# Visualize this response
temp_range <- range(rbind(north_bt %>% select(year, temp_c, period), all_proj_df)$temp_c)
temp_range <- c(floor(temp_range[1]), ceiling(temp_range[2]))

# Function to draw the template gaussian response
all_temps <- seq(temp_range[1], temp_range[2], by=0.1)

plot(gauss_rec(tx=all_temps, topt=topt_1, twidth=2), type="l")

# Trying to introduce symmetry
all_temps <- seq(-4.3,4.3, by=0.1)

plot(gauss_rec(tx=all_temps, topt=topt_1, twidth=2), type="l")

# Answer to find the limits
# [topt - 4×twidth, topt + 4×twidth] captures virtually all of the meaningful area of the function.


ecov <- list(label = c("North_BT"))

# original means
ecov$mean <- cbind(north_bt[,"mean"])

# center each column
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean, 2, mean))

# extend to include the feedback years
ecov$mean <- rbind(ecov$mean, matrix(0, n_feedback_years, ncol(ecov$mean)))




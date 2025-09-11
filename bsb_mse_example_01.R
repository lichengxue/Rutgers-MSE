# -----------------------------------------------------------------------------
#
#           Management Strategy Evaluation (MSE) Example for
#                       Black Sea Bass (BSB)
#                           BSB MSE - Mark I
#
#
# Author(s): Chengxue Li [Original], Jeewantha Bandara [This version]
# Last change: 2025/08/29
# Runtime environment: MacOS Sequoia 15.5 on M-chip Macbook Pro, R version 4.4.1
#
#
# Description:
# This script demonstrates the full workflow for a Management Strategy
# Evaluation (MSE) using the 'wham' and 'whamMSE' packages. It covers:
#   1. Setting up the Operating Model (OM) based on historical BSB data.
#   2. Configuring the Estimation Model (EM) used for in-loop assessments.
#   3. Executing the MSE simulation loop.
#   4. Visualizing the results.
#
# -----------------------------------------------------------------------------


# --- 1. INSTALLATION & SETUP ----

# Run these lines only if you need to install or update the packages.
# The file paths are examples and must be changed to your local directories.
# -----------------------------------------------------------------------------
# devtools::install_github("timjmiller/wham", dependencies=FALSE)
# devtools::install_github("lichengxue/whamMSE", dependencies=TRUE)

# Load required libraries
# -----------------------------------------------------------------------------
library(wham)
library(whamMSE)
library(tidyverse)
library(here)

here::here()
# Set working directory
# IMPORTANT: Set this to your project's working directory where the data
# files (e.g., "WAA.RDS") are located.
# -----------------------------------------------------------------------------
# setwd("/Users/jeewanthabandara/Graduate Studies/Rutgers-MSE-Setup")


# --- 2. DEFINE MODEL DIMENSIONS & PERIODS ----

# This section defines the core dimensions of the simulation, including time,
# populations, regions, fleets, and biological characteristics.
# -----------------------------------------------------------------------------

# Define the historical period for the model
years <- 1989:2023

# Define the MSE feedback (projection) period in years
MSE_years <- 10

# Define the number of seasons and the fraction of the year in each 'season'
# The larger the number of seasons, the greater the chance of likelihood errors
# Another potential fix is lowering the number of seasons - But keep this for last
n_seasons <- 11
fracyr_seasons <- c(rep(1, 5), 2, rep(1, 5)) / 12

# Define the number of 'stocks' (populations or subpopulations)
n_stocks <- 2

# Define the number of geographic regions
n_regions <- 2

# Define the total number of fishing fleets across all regions
# Black sea bass has 4 total fishing fleets in the stock assessment
# A recreational and commercial fleet for each subunit
n_fleets <- 4

# Define the total number of surveys across all regions
n_indices <- 2

# Define the number of age classes
n_ages <- 8

# Define the fraction of the year when the population spawns
# 0.5 means they spawn in the middle of the year
fracyr_spawn <- 0.5


# --- 3. PREPARE INPUT DATA ----

# This section loads or defines the biological and fishery data that will be
# used to parameterize the Operating Model.
# -----------------------------------------------------------------------------

# 3a. Weight-at-Age (WAA) ----
# Here, we assume time-invariant WAA, but the framework allows for
# time-varying data. The code reads a base WAA matrix and replicates it
# across all simulation years.
waa <- readRDS("WAA.RDS") # Make sure this file is in your working directory
user_waa <- list()
# Dimension: indices x years x ages
user_waa$waa <- array(NA, dim = c(8, length(years) + MSE_years, 8))
for (age in 1:8) {
  for (indx in 1:8) {
    user_waa$waa[indx, , age] <- rep(waa[indx, age], length(years) + MSE_years)
  }
}
# Pointers to assign the correct WAA matrix to fleets, indices, SSB, and M
user_waa$waa_pointer_fleets <- 1:4
user_waa$waa_pointer_indices <- 5:6
user_waa$waa_pointer_ssb <- 7:8 # The stock itself
user_waa$waa_pointer_M <- 7:8 # The stock itself


# 3b. Maturity-at-Age (MAA) ----
# We define time-invariant maturity-at-age for the two stocks.
user_maa <- rbind(
  c(0, 0.51, 0.95, 1, 1, 1, 1, 1),
  c(0, 0.45, 0.85, 0.95, 0.99, 1, 1, 1)
)


# 3c. Fleet and Survey Information ----
# Define which fleets/surveys operate in which regions and specify their
# associated observation error models.
# 
# Fleet Information
fleet_regions <- c(1, 1, 2, 2) # Commercial, Recreational, Commercial, Recreational
catch_Neff <- c(50, 50, 50, 50) # Effective sample size for age comps
catch_cv <- c(0.05, 0.15, 0.05, 0.15) # CV for total catch

# Survey Information
index_regions <- c(1, 2)
index_Neff <- c(50, 50) # Effective sample size for age comps
index_cv <- c(0.4, 0.4) # CV for survey index
# Survey catchability (q), from the BSB stock assessment
q <- c(7.370412e-05, 2.832577e-05)
# Fraction of the year when each survey is conducted
fracyr_indices <- c(0.25, 0.25)


# --- 4. CONFIGURE THE OPERATING MODEL (OM) ----

# This section builds the components of the Operating Model, which represents
# the "true" underlying population and fishery dynamics.
# -----------------------------------------------------------------------------

# 4a. Generate Basic Information List ----
# The `generate_basic_info` function consolidates the parameters into a
# list required for building the OM.
info <- generate_basic_info(
  base.years = years,
  n_feedback_years = MSE_years,
  n_stocks = n_stocks,
  n_regions = n_regions,
  n_indices = n_indices,
  n_fleets = n_fleets,
  n_ages = n_ages,
  n_seasons = n_seasons,
  catch_info = list(catch_cv = catch_cv, catch_Neff = catch_Neff),
  index_info = list(index_cv = index_cv, index_Neff = index_Neff, q = q, fracyr_indices = fracyr_indices),
  fracyr_seasons = fracyr_seasons,
  fracyr_spawn = 0.5,
  user_waa = user_waa,
  user_maturity = user_maa
)

# Extract the generated lists for easier access
basic_info <- info$basic_info
catch_info <- info$catch_info
index_info <- info$index_info
F_info <- info$F


# 4b. Update Historical Fishing Mortality (F) ----
# Override the default historical F with estimates from the BSB assessment.
F_est <- readRDS("F_hist.RDS") # Make sure this file is in your working directory
F_info$F[1:35, ] <- F_est


# 4c. Configure the Movement Model ----
# This section defines the movement dynamics between regions, which can be complex.
# It specifies where stocks can be, when they can move, and movement probabilities.
basic_info$NAA_where <- array(1, dim = c(n_stocks, n_regions, n_ages))
basic_info$NAA_where[1, 2, 1] <- 0 # Stock 1, age 1 cannot be in region 2 on Jan 1
basic_info$NAA_where[2, 1, ] <- 0  # Stock 2 does not move into region 1

move <- list(stock_move = c(TRUE, FALSE), separable = TRUE)
move$must_move <- array(0, dim = c(2, 11, 2))
move$must_move[1, 5, 2] <- 1
move$can_move <- array(0, dim = c(2, 11, 2, 2))
move$can_move[1, 1:4, 2, 1] <- 1
move$can_move[1, 7:11, 1, 2] <- 1
move$can_move[1, 5, 2, ] <- 1
mean_vals <- array(0, dim = c(2, length(fracyr_seasons), 2, 1))
mean_vals[1, 1:11, 1, 1] <- 0.02214863
mean_vals[1, 1:11, 2, 1] <- 0.3130358
move$mean_vals <- mean_vals
move$mean_model <- matrix("stock_constant", 2, 1)


# 4d. Configure Numbers-at-Age (NAA) Dynamics ----
# Define process errors for recruitment and survival. Here we assume recruitment
# is random around a mean (density-independent).
sigma <- "rec+1" # Full state-space model where all numbers at age are random effects
re_cor <- "iid" # iid - Independent and identically distributed covariate
ini.opt <- "age-specific-fe" # Need to figure out what this means....
sigma_vals <- array(0.2, dim = c(n_stocks, n_regions, n_ages)) # NAA survival sigma - Another potential fix is increasing this value
sigma_vals[, , 1] <- 0.75 # Recruitment sigma/error

NAA_re <- list(
  N1_model = rep(ini.opt, n_stocks),
  sigma = rep(sigma, n_stocks),
  cor = rep(re_cor, n_stocks),
  recruit_model = 2, # 2 = density-independent. Random about mean
  recruit_pars = list(16000, 21000), # Mean recruitment for each stock
  sigma_vals = sigma_vals
)


# 4e. Configure Natural Mortality (M) ----
# Assume a time-invariant, age-constant natural mortality (M).
M <- list(model = "constant", initial_means = array(0.4, dim = c(n_stocks, n_regions, n_ages)))


# 4f. Configure Gear Selectivity ----
# Specify selectivity patterns (e.g., age-specific, logistic) for each fleet and survey.
model <- c("age-specific", "age-specific", "logistic", "logistic", "age-specific", "age-specific")
initial_pars_f1 <- c(0.031, 0.421, 0.858, 1, 1, 1, 1, 1)
initial_pars_f2 <- c(0.036, 0.355, 0.573, 0.782, 0.868, 0.928, 1, 1)
initial_pars_f3 <- c(2.381, 0.402)
initial_pars_f4 <- c(2.788, 0.849)
initial_pars_s1 <- c(0.091, 0.463, 0.905, 0.962, 1, 1, 1, 1)
initial_pars_s2 <- c(0.362, 1, 1, 1, 1, 1, 1, 1)

sel <- list(
  model = model,
  initial_pars = c(
    list(initial_pars_f1), list(initial_pars_f2),
    list(initial_pars_f3), list(initial_pars_f4),
    list(initial_pars_s1), list(initial_pars_s2)
  ),
  fix_pars = list(4:8, 7:8, NULL, NULL, 5:8, 2:8)
)


# --- 5. BUILD THE OM OBJECT & PREPARE FOR MSE ----

# This section finalizes the OM setup by consolidating all components.
# -----------------------------------------------------------------------------

# 5a. Prepare Final `wham` Input ----
# This critical step gathers all configurations into a single input list.
input <- prepare_wham_input(
  basic_info = basic_info, 
  NAA_re = NAA_re,
  move = move,
  M = M,
  selectivity = sel,
  catch_info = catch_info,
  index_info = index_info,
  F = F_info
)

# 5b. Post-processing Updates ----
# Some inputs are easier to update after the main list is created.
waa_info <- info$par_inputs$user_waa
input <- update_waa(input, waa_info = waa_info)

# Set initial numbers-at-age for the first year of the historical period.
ini.NAA <- matrix(NA, n_ages, n_stocks)
ini.NAA[, 1] <- c(4531.15069, 2876.48828, 1484.19315, 692.67290, 314.73335, 142.50149, 64.41915, 52.58156)
ini.NAA[, 2] <- c(16023.94518, 10167.14058, 5683.71704, 2523.02265, 972.93721, 353.59712, 125.80803, 68.66246)
input$par$log_N1[] <- 0
for (i in 1:n_regions) {
  input$par$log_N1[i, i, ] <- log(ini.NAA[, i])
}

# 5c. Generate the Operating Model Object ----
# We use `fit_wham` with `do.fit = FALSE` to create the OM structure without
# fitting it to data. This defines our "true" simulated world.
random <- input$random # Identify which processes are random effects - It's Numbers at age
input$random <- NULL   # Nullify so inner optimization won't change simulated REs
om <- fit_wham(input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
# do.brps - Calculate and report biological reference points

# --- 6. MSE CONFIGURATION & EXECUTION ----

# This section sets up the specifics of the MSE feedback loop, including
# assessment frequency and the harvest control rule, and then runs it.
# -----------------------------------------------------------------------------

# 6a. Specify MSE Timeline ----
assess.interval <- 3 # Assessments occur every 3 years
base.years <- years
terminal.year <- tail(base.years, 1)
assess.years <- seq(terminal.year, tail(om$years, 1) - assess.interval, by = assess.interval)

# Print the designated assessment years to verify
print(assess.years)


# 6b. Run the MSE Simulation Loop ----
#
# WARNING: The following MSE loop can be computationally intensive and take a
# long time to run. It is wrapped in `if (FALSE)` to prevent accidental execution.
# Change to `if (TRUE)` to run the simulation.
#
if (TRUE) {
  
  # Generate a single realization of the OM (with one set of random effects)
  om_with_data <- update_om_fn(om, seed = 123, random = random)
  
  # Configure the Estimation Model (EM) for use within the loop
  move_em <- move
  move_em$use_prior <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))
  move_em$use_prior[1, 1, , ] <- 1
  # Another potential fix - Increase the value here to increase the confidence of the movement rate
  move_em$prior_sigma <- array(0.2, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))
  
  # Specify the Harvest Control Rule (HCR)
  hcr <- list()
  hcr$hcr.type <- 1 # FXSPR - Fishing pressure to keep the SPR at a certain percentage
  hcr$hcr.opts <- list(use_FXSPR = TRUE, percentFXSPR = 75) # Apply F at 75% unfished SPR
  
  # Potential fix to the non-convergence. This seemed to fix the non-convergence issues
  NAA_re$N1_model[] = "equilibrium"
  
  # Execute the MSE loop for one realization
  mod <- loop_through_fn(
    om = om_with_data,
    em_info = info,
    random = random,
    sel_em = sel,
    M_em = M,
    NAA_re_em = NAA_re,
    move_em = move_em,
    em.opt = list(
      separate.em = FALSE,
      separate.em.type = 3,
      do.move = TRUE,
      est.move = TRUE # Another potential fix to the non-convergence issues is to set this to FALSE
    ),
    assess_years = assess.years,
    assess_interval = assess.interval,
    base_years = base.years,
    year.use = 35,
    add.years = TRUE,
    seed = 123,
    hcr = hcr,
    save.last.em = TRUE
  )
}


# --- 7. VISUALIZE MSE RESULTS ----

# After the MSE is complete, this section shows how to use the built-in
# plotting functions to summarize and visualize the results.
# -----------------------------------------------------------------------------

# The code below is an example of how to call the plotting function. It will not
# run without a 'mods' object, which would be a list of completed MSE runs from
# the loop in Section 6.
# It is wrapped in `if (FALSE)` to prevent errors.
#

# plot_wham_output(mod$em_full[[1]], out.type = "html")



if (TRUE) {
  
  plot_mse_output(
    mods<-list(mod), # This would be a list of results, e.g., mods <- list(mod1, mod2)
    main_dir = getwd(),
    output_dir = "Report",
    output_format = "pdf", # or "html" or "png"
    width = 10, height = 7, dpi = 300,
    col.opt = "D",
    method = "mean",
    outlier.opt = NA,
    f.ymax = 2, # control y-axis scale
    new_model_names = c("FXSPR75"),
    base.model = 'FXSPR75',
    start.years = 31,
    use.n.years.first = 3,
    use.n.years.last = 3
  )
  
}

# Custom plots for comparing estimation model vs. operating model

# OM SSB vs. EM SSB

# No. of MSE years
assess.years.length <- length(assess.years)

# All SSB estimates from the EM from the final MSE
em_ssb <- mod$em_list[[assess.years.length]]$SSB[,]

om_ssb <- mod$om$rep$SSB

# Transform these matrices to a dataframe and plot them nicely alongside each other
em_ssb_df <- as.data.frame(em_ssb)
colnames(em_ssb_df) <- c("Stock 1","Stock 2")
em_ssb_df <- em_ssb_df %>% mutate(year = append(years,tail(years)+MSE_years))
em_ssb_df <- pivot_longer(em_ssb_df, cols=c("Stock 1","Stock 2"), names_to=c("Stock"), values_to=c("SSB") )
# em_ssb_df <- em_ssb_df %>% mutate(Stock = replace(Stock, Stock=="V1","Stock 1"))
# em_ssb_df <- em_ssb_df %>% mutate(Stock = replace(Stock, Stock=="V2","Stock 2"))

om_ssb
om_ssb_df <- as.data.frame(om_ssb)
colnames(om_ssb_df) <- c("Stock 1", "Stock 2")
om_ssb_df$year <- seq.int(nrow(om_ssb_df)) + (years[1]-1)
om_ssb_df <- pivot_longer(om_ssb_df, cols=c("Stock 1", "Stock 2"), names_to=c("Stock"), values_to=c("SSB"))
om_ssb_df <- om_ssb_df %>% relocate(year)

# Join the two dataframes together
em_ssb_df <- em_ssb_df %>% mutate(model="Estimation Model")
om_ssb_df <- om_ssb_df %>% mutate(model="Operating Model")

# Bind the two together
m_ssb_df <- rbind(em_ssb_df, om_ssb_df)

# Make a plot of em vs. om SSB
em_vs_om_plot <- ggplot(m_ssb_df, aes(year, SSB, color=model)) + geom_line(linewidth=1.5, alpha=0.75) + 
  facet_wrap(~Stock, nrow=2, scales="free") + 
  geom_vline(xintercept=c(assess.years), linetype="dotted", linewidth=0.5, alpha=0.5) + 
  labs(x="Year", y="SSB", color="Model") + 
  scale_x_continuous(breaks=seq(1990,2035,3)) + 
  scale_y_continuous(limits=c(0,35000)) + 
  scale_linetype_manual(values=c("dotted"), name="Assessment year") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.background = element_blank(),
        axis.line = element_line())

# Tried to make a legend entry for assessment years. But not currently working
em_vs_om_plot
  
ggsave(here("plots","em_vs_om_plot.png"), em_vs_om_plot, height=6, width=10, units=c("in"), dpi=300)



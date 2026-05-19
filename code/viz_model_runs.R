# -----------------------------------------------------------------------------
#
#                 Management Strategy Evaluation (MSE)
#           Black Sea Bass (BSB) with Environmental Drivers
#                       VISUALIZING MODEL RUNS
#
#
# Author(s): RMWJ Bandara, Chengxue Li
# Date: 2026/05/18
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
library(doParallel)

here::here()

model_run <- "2026-05-18_16-37-38"

all_model_files <- list.files(path=here("models","sensitivity_analysis",model_run,"models"), 
                          pattern="\\.RDS$", 
                          include.dirs = TRUE,
                          full.names = TRUE)

# Regex pattern for extracting the configuration, iteration, and model number
pattern <- "sens_run_(\\d+)_iter_id_(\\d+)_mod_(\\d+)\\.RDS"

# Read in all the files
# Function for parsing the list of models using the Regex pattern
# This returns a nice dataframe that contains details about the iterations, configs, and model number
parse_filename <- function(fname) {
  m <- regmatches(fname, regexec(pattern, fname))[[1]]
  if (length(m) == 0) return(NULL)
  data.frame(
    filename = fname,
    run      = as.integer(m[2]),
    iter_id  = as.integer(m[3]),
    mod      = as.integer(m[4])
  )
}

model_df <- do.call(rbind, lapply(all_model_files, parse_filename))
print(model_df)

# New column that contains Model #
model_df <- model_df %>% mutate(model_name=paste("Model-",mod,sep=""))

# Read these into a single list using `lapply`
all_models <- lapply(model_df$filename, readRDS)
# Name these models
names(all_models) <- model_df %>% select(model_name) %>% pull()


# Adopting this function from this vignette
# https://lichengxue.github.io/SPASAM.MSE/Performance-Analysis-Tools.html
plot_mse_output(all_models,
               main_dir = getwd(),
               output_dir = "Report",
               output_format = c("html"), # or html or png
               width = 10, height = 7, dpi = 300,
               col.opt = "D",
               # new_model_names = c("M1","M2","M3","M4","M5"),
               # base.model = "M1",
               # start.years = 31,
               # use.n.years.first = 5,
               # use.n.years.last = 5
               )

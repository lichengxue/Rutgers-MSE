# -----------------------------------------------------------------------------
#
#                 Management Strategy Evaluation (MSE)
#           Black Sea Bass (BSB) with Environmental Drivers
#                       VISUALIZING MODEL RUNS
#
#
# Author(s): RMWJ Bandara, Chengxue Li
# Date: 2026/05/18 (Last update on 2026/06/17)
# Runtime environment: MacOS Sequoia 15.5 on M-chip Macbook Pro, R version 4.4.1
#
#
#
# -----------------------------------------------------------------------------

library(wham)
# library(whamMSE)
library(SPASAM.MSE)
library(dplyr)
library(here)
library(ggplot2)
library(gridExtra)
library(beepr)
library(readr)
library(doParallel)
library(stringr)

here::here()

# model_run <- "2026-05-26_11-15-55"
model_run <- "2026-06-17_09-50-46"

all_model_files <- list.files(path=here("models","sensitivity_analysis",model_run,"models"), 
                          pattern="\\.RDS$", 
                          include.dirs = TRUE,
                          full.names = TRUE)

# Regex pattern for extracting the configuration, iteration, and model number
pattern <- "sens_run_(\\d+)_iter_id_(\\d+)_mod_(\\d+)\\.RDS"
error_pattern <- "sens_run_(\\d+)_iter_id_(\\d+)_error_log_mod_(\\d+)\\.RDS"

# Handle errors using the 'log.txt' file
log_file <- readLines(here("models","sensitivity_analysis",model_run,"logs","logs.txt"))
matches <- str_match(log_file, '\\[1\\] "Model error in run_id: (\\d+) and iter_id: (\\d+)"')

# This dataframe is important because it tells us which models possibly had convergence issues
# Sometimes, model objects can be a non-Null object even if they contain gibberish or
# resulted from a bad model run
aborted_runs <- data.frame(
  run  = as.integer(matches[!is.na(matches[, 1]), 2]),
  iter_id = as.integer(matches[!is.na(matches[, 1]), 3])
)

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

parse_error_filename <- function(fname) {
  m <- regmatches(fname, regexec(error_pattern, fname))[[1]]
  if (length(m) == 0) return(NULL)
  data.frame(
    filename = fname,
    run      = as.integer(m[2]),
    iter_id  = as.integer(m[3]),
    mod      = as.integer(m[4])
  )
}


valid_model_df <- do.call(rbind, lapply(all_model_files, parse_filename))
error_model_df <- do.call(rbind, lapply(all_model_files, parse_error_filename))
print(valid_model_df)
print(error_model_df)
if(!is.null(valid_model_df)){
  valid_model_df <- valid_model_df %>% mutate(success="Yes")
  
}
if(!is.null(error_model_df)){
  error_model_df <- error_model_df %>% mutate(success="No")
}

model_df <- rbind(valid_model_df, error_model_df)

# Introduce other parameters such as resulting file size in there
model_df$size_mb <- round(file.info(model_df$filename)$size / 1024^2, 2)

# An easy way of saying whether a model has converged or not
model_df <- model_df %>% mutate(success=ifelse(size_mb>0, "Yes","No"))

# New column that contains Model #
model_df <- model_df %>% mutate(model_name=paste("Model-",mod,sep=""))

# Write this to disk temporarily
write.csv(model_df, here("models","sensitivity_analysis",model_run,"model_meta_info.csv"))

# Select only the readable models
model_df <- model_df %>% filter(success=="Yes")

# Condense these into a list (run settings) of list (models) of lists (iterations/seeds)

model_list <- model_df %>%
  filter(success == "Yes") %>%
  arrange(run, mod, iter_id) %>%
  split(.$run) %>%                              # top level: run
  lapply(function(run_df) {
    run_df %>%
      split(.$mod) %>%                          # second level: mod
      lapply(function(mod_df) {
        mod_df %>%
          split(.$iter_id) %>%                  # third level: iter_id
          setNames(paste0("Model - ", mod_df$mod[1])) %>%  # all iters get the same model name
          lapply(function(row) {
            readRDS(row$filename)               # read the .RDS model object
          })
      })
  })


model_list_2 <- model_df %>%
  filter(success == "Yes") %>%
  arrange(run, iter_id, mod) %>%
  split(.$run) %>%                          # top level: run
  lapply(function(run_df) {
    run_df %>%
      split(.$iter_id) %>%                  # second level: iter_id
      lapply(function(iter_df) {
        iter_df %>%
          split(.$mod) %>%                  # third level: model
          setNames(iter_df$model_name) %>%  # named by model_name
          lapply(function(row) {
            readRDS(row$filename)
          })
      })
  })


# Create a separate object for Run 1
run_1_models <- model_list[["1"]]

length(model_list[["1"]][["1"]])


# Read these into a single list using `lapply`
all_models <- lapply(model_df$filename, readRDS)
# Name these models
names(all_models) <- model_df %>% select(model_name) %>% pull()



#### NEW CODE BLOCK ####

model_dir <- here("models","sensitivity_analysis","2026-05-26_10-21-15","models")
run_id <- 1
nsim <- 6
model_nums <- 1:3

mods <- lapply(1:nsim, function(r) {
  
  mod_list <- lapply(model_nums, function(m) {
    
    file_path <- file.path(
      model_dir,
      sprintf("sens_run_%d_iter_id_%d_mod_%d.RDS", run_id, r, m)
    )
    
    readRDS(file_path)
  })
  
  names(mod_list) <- paste0("Mod", model_nums)
  return(mod_list)
})

#### END NEW CODE BLOCK ####


# Adopting this function from this vignette
# https://lichengxue.github.io/SPASAM.MSE/Performance-Analysis-Tools.html
plot_mse_output(mods,
               main_dir = getwd(),
               output_dir = "Report-May-21",
               output_format = c("html"), # or html or png
               width = 10, height = 7, dpi = 300,
               col.opt = "D",
               # new_model_names = c("M1","M2","M3","M4","M5"),
               # base.model = "M1",
               # start.years = 31,
               # use.n.years.first = 5,
               # use.n.years.last = 5
               )

# Check for my method
plot_mse_output(model_list_2[["1"]],
                main_dir = getwd(),
                output_dir = "Report-May-26-2",
                output_format = c("html"), # or html or png
                width = 10, height = 7, dpi = 300,
                col.opt = "D",
                # new_model_names = c("M1","M2","M3","M4","M5"),
                # base.model = "M1",
                # start.years = 31,
                # use.n.years.first = 5,
                # use.n.years.last = 5
)


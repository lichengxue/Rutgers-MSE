# -----------------------------------------------------------------------------
#
#                 Management Strategy Evaluation (MSE)
#           Black Sea Bass (BSB) with Environmental Drivers
#                       VISUALIZING MODEL RUNS
#
#
# Author(s): RMWJ Bandara, Chengxue Li
# Date: 2026/05/18 (Last update on 2026/07/16)
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
library(purrr)

here::here()

# model_run <- "2026-05-26_11-15-55"
# model_run <- "2026-06-18_16-33-58"
# model_run <- "2026-06-19_12-52-15"
# model_run <- "2026-06-23_13-33-03"
# model_run <- "2026-06-24_10-10-02"
# model_run <- "2026-06-24_11-34-32"
# model_run <- "2026-06-25_09-18-13"
model_run <- "2026-07-07_13-48-04"


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
model_df <- NULL
print(valid_model_df)
print(error_model_df)
if(!is.null(valid_model_df)){
  valid_model_df <- valid_model_df %>% mutate(success="Yes")
  
}
if(!is.null(error_model_df)){
  error_model_df <- error_model_df %>% mutate(success="No")
}

# New data frame that will hold all the valid model runs
model_df <- valid_model_df

# Make sure that unconverging iterations are dropped
if(!is.null(error_model_df)){
  unconverged_iterations <- error_model_df %>% select(iter_id) %>% distinct() %>% pull()
  model_df <- model_df %>% filter(!iter_id %in% unconverged_iterations)
}

# model_df <- rbind(valid_model_df, error_model_df)

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

# Calculate the convergence rate
config_file <- read_csv(here("models","sensitivity_analysis",model_run,"model_configs.csv"))
config_iters <- config_file %>% filter(config=="iterations") %>% select(val) %>% pull()
converged_iters <- model_df %>% distinct(iter_id) %>% count() %>% pull()
convergence_rate <- converged_iters/config_iters*100
print(paste("Convergence rate was",round(convergence_rate),"% with",converged_iters,
            "converging from a total of",config_iters,"seeds", sep=" "))

# Condense these into a list (run settings) of list (models) of lists (iterations/seeds) of
# configurations (runs)
model_list <- model_df %>%
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


# Make a directory for model reports
model_reports_path <- here("models", "sensitivity_analysis", model_run, "model-reports")
dir.create(model_reports_path, recursive = TRUE, showWarnings = FALSE)

model_list_length <- length(model_list) # How many different configurations are there

for(i in 1:model_list_length){
  # Visualize model performance statistics
  plot_mse_output(model_list[[as.character(i)]],
                  main_dir = model_reports_path, #getwd() is default
                  output_dir = paste("Run-",i,sep=""),
                  output_format = c("html"), # or html or png
                  width = 10, height = 7, dpi = 300,
                  col.opt = "D",
                  new_model_names = c("Linear relationship","No relationship","Gaussian relationship"),
                  base.model = "No relationship",
                  start.years = 36, # This starts at 36. But check. 2023 is the start year of assessment
                  use.n.years.first = 5,
                  use.n.years.last = 5
                  # new_model_names = c("M1","M2","M3","M4","M5"),
                  # base.model = "M1",
  )
}


#### Inter-model comparison between models for different configurations ####

# Get the model settings
model_settings <- read_csv(here("models","sensitivity_analysis",model_run,"all_model_settings.csv"))
model_settings <- model_settings %>% group_by(run_id) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  select(-c(model, model_path))
model_settings <- model_settings %>% mutate(run_name=paste("Run-",run_id, sep=""), .after=1)


# Read in the `mse_performance_summary.csv` files for each of the runs and bind them together
model_result_summaries <- model_settings$run_name |>
  map(\(rn) {
    f <- here("models", "sensitivity_analysis", model_run, "model-reports",
              rn, "mse_performance_summary.csv")   # adjust name per above
    if (!file.exists(f)) return(NULL)
    read_csv(f) |> mutate(run_name = rn, .before=1)
  }) |>
  bind_rows()

# Introduce a new column (scenario) - Simply for naming/cleanliness purposes
model_result_summaries <- model_result_summaries %>% 
  mutate(scenario = str_replace(run_name, "Run", "Scenario"))

# The unique stats computed for the models during `plot_mse_output`
model_result_summaries %>% distinct(metric_detail)
print(model_result_summaries %>% distinct(metric), n=25)


##### Catch in the last 5 years #####
catch_last_results <- model_result_summaries %>% filter(metric=="Catch_last")


ggplot(catch_last_results, aes(x = scenario, fill = Model)) +
  geom_boxplot(
    aes(ymin = min, lower = q1, middle = median, upper = q3, ymax = max),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(x = "Model configuration", y = "Catch in the last 5 years") +
  theme_bw()


##### SSB #####

ssb_results <- model_result_summaries %>% filter(metric=="SSB" & level=="global" & period=="from_start_to_end")


ggplot(ssb_results, aes(x = scenario, fill = Model)) +
  geom_boxplot(
    aes(ymin = min, lower = q1, middle = median, upper = q3, ymax = max),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(x = "Model configuration", y = "SSB") +
  theme_bw()

##### SSB - AAV #####

ssb_aav_results <- model_result_summaries %>% filter(metric=="SSB" & metric_detail=="AAV" & level=="global")


ggplot(ssb_aav_results, aes(x = scenario, fill = Model)) +
  geom_boxplot(
    aes(ymin = min, lower = q1, middle = median, upper = q3, ymax = max),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(x = "Model configuration", y = "SSB [Average Annual Variation]") +
  theme_bw()


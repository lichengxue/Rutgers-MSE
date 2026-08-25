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


# Benchmarking run time and info to save files
viz_start_time <- Sys.time()

# model_run <- "2026-05-26_11-15-55"
# model_run <- "2026-06-18_16-33-58"
# model_run <- "2026-06-19_12-52-15"
# model_run <- "2026-06-23_13-33-03"
# model_run <- "2026-06-24_10-10-02"
# model_run <- "2026-06-24_11-34-32"
# model_run <- "2026-06-25_09-18-13"
# model_run <- "2026-07-07_13-48-04"


all_model_runs <- c("2026-08-19_16-21-49", "2026-08-19_15-45-40")

# Build full paths to each all_seeds.txt file
seed_files <- here("models","sensitivity_analysis", all_model_runs, "all_seeds.csv")

#### VERIFY SEEDS AND CONFIGURATIONS ARE IDENTICAL ####
check_files_identical <- function(filename, run_folders,
                                  subpath = c("models", "sensitivity_analysis")) {
  
  # Collapse subpath into one fixed prefix so it doesn't get recycled
  # element-wise against run_folders
  prefix <- paste(subpath, collapse = "/")
  
  # Build full paths to each file
  files <- here(prefix, run_folders, filename)
  
  # Confirm they all exist first
  if (!all(file.exists(files))) {
    stop("Missing ", filename, " in: ",
         paste(files[!file.exists(files)], collapse = ", "))
  }
  
  # Read contents of each file
  contents <- lapply(files, readLines)
  
  # Compare all to the first one
  all_identical <- all(vapply(contents[-1], identical, logical(1), contents[[1]]))
  
  if (all_identical) {
    message("✅ ", filename, " is identical across all folders.")
  } else {
    message("❌ ", filename, " differs between folders.")
    
    diffs <- which(!vapply(contents[-1], identical, logical(1), contents[[1]])) + 1
    message("Differing folders: ", paste(run_folders[diffs], collapse = ", "))
  }
  
  invisible(all_identical)
}

check_files_identical("all_seeds.csv", all_model_runs)
check_files_identical("model_configs.csv", all_model_runs)

#### VERIFY SETTINGS ARE DIFFERENT ####

cols_to_extract <- c("run_id","proc_error", "mse_gap", "gauss_width", 
                     "temp_trend", "temp_trend_stoch", "topt_rec")

extract_settings <- function(filename, run_folders, cols,
                             subpath = c("models", "sensitivity_analysis")) {
  
  prefix <- paste(subpath, collapse = "/")
  files <- here(prefix, run_folders, filename)
  
  if (!all(file.exists(files))) {
    stop("Missing ", filename, " in: ",
         paste(files[!file.exists(files)], collapse = ", "))
  }
  
  # Read each file, select the requested columns, dedupe, tag with run folder
  data_list <- Map(function(path, run) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    
    missing_cols <- setdiff(cols, names(df))
    if (length(missing_cols) > 0) {
      warning("In ", run, ", missing columns: ", paste(missing_cols, collapse = ", "))
    }
    
    df <- df[intersect(cols, names(df))]
    df <- distinct(df)
    df$run_folder <- run
    df
  }, files, run_folders)
  
  # Combine into one data frame
  settings_df <- bind_rows(data_list)
  
  # Create a clean sequential run_id, unique across all combined folders
  settings_df <- settings_df %>%
    mutate(run_id_new = row_number())
  
  settings_df
}

settings_df <- extract_settings("all_model_settings.csv", all_model_runs, cols_to_extract)
settings_df

#### List all the model files across multiple runs ####

list_model_files <- function(run_folders,
                             subpath = c("models", "sensitivity_analysis"),
                             models_subfolder = "models") {
  
  prefix <- paste(subpath, collapse = "/")
  
  file_list <- lapply(run_folders, function(run) {
    dir_path <- here(prefix, run, models_subfolder)
    
    if (!dir.exists(dir_path)) {
      warning("Missing models folder for run: ", run)
      return(NULL)
    }
    
    files <- list.files(dir_path, pattern = "\\.RDS$", full.names = TRUE)
    
    data.frame(
      run_folder = run,
      file_path = files,
      file_name = basename(files),
      stringsAsFactors = FALSE
    )
  })
  
  all_files <- bind_rows(file_list)
  
  # Parse out run_id (x), nid (y), model number (z), and error flag
  all_files <- all_files %>%
    mutate(
      is_error = str_detect(file_name, "_error_log_mod_"),
      run_id   = as.integer(str_match(file_name, "sens_run_(\\d+)_")[, 2]),
      nid      = as.integer(str_match(file_name, "iter_id_(\\d+)_")[, 2]),
      model    = as.integer(str_match(file_name, "mod_(\\d+)\\.RDS$")[, 2])
    )
  
  all_files
}

# all_model_runs <- c("2026-08-19_16-21-49", "2026-08-19_15-45-40")
model_files_df <- list_model_files(all_model_runs)
model_files_df

# Join `settings_df` and `model_files_df`
model_files_df <- model_files_df %>%
  left_join(settings_df, by = c("run_folder", "run_id"))

model_files_df

colnames(model_files_df)


# Identify nids that had at least one error, anywhere across all runs
error_nids <- model_files_df %>%
  filter(is_error) %>%
  distinct(nid)

# Drop all rows for those nids, across every run_id_new
model_files_df_clean <- model_files_df %>%
  anti_join(error_nids, by = "nid")

model_files_df_clean

nrow(model_files_df) - nrow(model_files_df_clean)
n_distinct(error_nids)   # number of distinct seeds/iterations dropped


model_list <- model_files_df_clean %>%
  filter(!is_error) %>%
  arrange(run_id_new, nid, model) %>%
  split(.$run_id_new) %>%                     # top level: run
  lapply(function(run_df) {
    run_df %>%
      split(.$nid) %>%                        # second level: iter_id/seed
      lapply(function(iter_df) {
        iter_df %>%
          split(.$model) %>%                  # third level: model
          setNames(paste0("model_", iter_df$model)) %>%  # named by model number
          lapply(function(row) {
            readRDS(row$file_path)
          })
      })
  })


viz_folder_name <- format(viz_start_time, "%Y-%m-%d_%H-%M-%S")
viz_folder_path <- here("models","gathered_viz",viz_folder_name)
# Add which file was used to the model_run_info data
# model_run_info <- model_run_info %>% add_row(model_spec="which_file", vals=file_used, .before = 1)
# Create a folder. Suppress warnings and allow recursive folders to be created
# dir.create(viz_folder_path, recursive = TRUE, showWarnings = FALSE)

# Make a directory for model reports
model_reports_path <- here("models","gathered_viz", viz_folder_name, "model-reports")
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


settings_df <- settings_df %>% mutate(run_name=paste("Run-",run_id_new, sep=""), .after=1)

# Read in the `mse_performance_summary.csv` files for each of the runs and bind them together
model_result_summaries <- settings_df$run_name |>
  map(\(rn) {
    f <- here("models", "gathered_viz", viz_folder_name, "model-reports",
              rn, "mse_performance_summary.csv")   # adjust name per above
    if (!file.exists(f)) return(NULL)
    read_csv(f) |> mutate(run_name = rn, .before=1)
  }) |>
  bind_rows()

# Introduce a new column (scenario) - Simply for naming/cleanliness purposes
model_result_summaries <- model_result_summaries %>% 
  mutate(scenario = str_replace(run_name, "Run", "Scenario"))


# Do the same for settings_df
settings_df %>% 
  mutate(scenario = str_replace(run_name, "Run", "Scenario"))


# The unique stats computed for the models during `plot_mse_output`
model_result_summaries %>% distinct(metric_detail)
print(model_result_summaries %>% distinct(metric), n=25)

# There are different details that's specifically computed for each metric
print(model_result_summaries %>% distinct(metric, metric_detail) %>% 
        arrange(metric, metric_detail), n=50)

# They are also computed for both local and global levels
# There are different details that's specifically computed for each metric
print(model_result_summaries %>% distinct(metric, metric_detail, level) %>% 
        arrange(metric, metric_detail), n=50)

# They are also computed for different time periods
print(model_result_summaries %>% distinct(metric, metric_detail, level, period) %>% 
        arrange(metric, metric_detail), n=75)



#### Catch ####

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

##### Catch in the first 5 years #####
catch_first_results <- model_result_summaries %>% filter(metric=="Catch_first")

ggplot(catch_first_results, aes(x = scenario, fill = Model)) +
  geom_boxplot(
    aes(ymin = min, lower = q1, middle = median, upper = q3, ymax = max),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(x = "Model configuration", y = "Catch in the last 5 years") +
  theme_bw()


#### SSB ####

##### SSB - Global - Start to finish #####

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




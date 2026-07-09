# This script is for setting up temperature trends that we will port over to 
# code/parallelized_historic_run.R

library(dplyr)
library(here)
library(ggplot2)
library(readr)
library(tidyr)
library(patchwork)
library(ggpmisc)

# Benchmarking run time and info for save files
run_start_time <- Sys.time()

# Set 'here'
here::here()

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
                        period="Projection - 1")


##### REPORTED TREND #####
scen_1_trend <- 0.04 # 0.04 deg C/year
scen_1_sd <- 0.05
scen_1_error <- rnorm(n=n_feedback_years, mean=0, sd=scen_1_sd)
scen_1_proj <- scen_1_trend*seq(1,n_feedback_years) +last_temp_c  + scen_1_error

scen_1_df <- data.frame(year=seq(first_proj_year-1, first_proj_year+n_feedback_years-1), 
                        temp_c=c(last_temp_c, scen_1_proj),
                        period="Projection - 2")

##### TREND FROM LAST 10 YEARS #####
last_5_years <- north_bt %>% arrange(year) %>% tail(n=10)
scen_2_lm <- lm(temp_c ~ year, data=last_5_years)
scen_2_trend <- scen_2_lm$coefficients[["year"]]
scen_2_sd <- 0.5
scen_2_error <- rnorm(n=n_feedback_years, mean=0, sd=scen_2_sd)
scen_2_proj <- scen_2_trend*seq(1,n_feedback_years)+ last_temp_c + scen_2_error

scen_2_df <- data.frame(year=seq(first_proj_year-1, first_proj_year+n_feedback_years-1), 
                        temp_c=c(last_temp_c, scen_2_proj),
                        period="Projection - 3")


##### REPORTED TREND - BUT MORE STOCHASTIC #####
scen_3_trend <- 0.04 # 0.04 deg C/year
scen_3_sd <- 0.25
scen_3_error <- rnorm(n=n_feedback_years, mean=0, sd=scen_3_sd)
scen_3_proj <- scen_3_trend*seq(1,n_feedback_years) +last_temp_c  + scen_3_error

scen_3_df <- data.frame(year=seq(first_proj_year-1, first_proj_year+n_feedback_years-1), 
                        temp_c=c(last_temp_c, scen_3_proj),
                        period="Projection - 4")


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
topt_4 <- -1.5
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


# Make nice plots for the theoretical optima
standard_twidth <- 2
# Trying to introduce symmetry

# Optima 1 and the trajectories
# all_temps <- seq(-4.3,4.3, by=0.1)
all_temps <- seq(topt_1 - 4*standard_twidth, topt_1 + 4*standard_twidth, by=0.05)
optima_1_response <- gauss_rec(tx=all_temps, topt=topt_1, twidth=standard_twidth)
optima_1_viz <- data.frame(temp=all_temps, response=optima_1_response)
ggplot(optima_1_viz, aes(temp, response)) + geom_line(linewidth=2) + 
  geom_rug(data=north_bt, aes(x=temp_c), color="black", sides="b", alpha=0.25, inherit.aes = FALSE) + 
  geom_rug(data=all_proj_df, aes(x=temp_c, color=period), inherit.aes=FALSE) + 
  labs(x="Standardized temperature (deg C)", y="Response", color="Temperature trend") + 
  bioshift_plot_theme_1 + theme(panel.grid.major=element_line(color="grey", linewidth=0.1))

# Optima 2 and the trajectories
all_temps <- seq(topt_2 - 4*standard_twidth, topt_2 + 4*standard_twidth, by=0.05)
optima_2_response <- gauss_rec(tx=all_temps, topt=topt_2, twidth=standard_twidth)
optima_2_viz <- data.frame(temp=all_temps, response=optima_2_response)
ggplot(optima_2_viz, aes(temp, response)) + geom_line(linewidth=2) + 
  geom_rug(data=north_bt, aes(x=temp_c), color="black", sides="b", alpha=0.25, inherit.aes = FALSE) + 
  geom_rug(data=all_proj_df, aes(x=temp_c, color=period), inherit.aes=FALSE) + 
  labs(x="Standardized temperature (deg C)", y="Response", color="Temperature trend") + 
  bioshift_plot_theme_1 + theme(panel.grid.major=element_line(color="grey", linewidth=0.1))


# Optima 3 and the trajectories
all_temps <- seq(topt_3 - 4*standard_twidth, topt_3 + 4*standard_twidth, by=0.05)
optima_3_response <- gauss_rec(tx=all_temps, topt=topt_3, twidth=standard_twidth)
optima_3_viz <- data.frame(temp=all_temps, response=optima_3_response)
ggplot(optima_3_viz, aes(temp, response)) + geom_line(linewidth=2) + 
  geom_rug(data=north_bt, aes(x=temp_c), color="black", sides="b", alpha=0.25, inherit.aes = FALSE) + 
  geom_rug(data=all_proj_df, aes(x=temp_c, color=period), inherit.aes=FALSE) + 
  labs(x="Standardized temperature (deg C)", y="Response", color="Temperature trend") + 
  bioshift_plot_theme_1 + theme(panel.grid.major=element_line(color="grey", linewidth=0.1))


# Optima 4 and the trajectories
all_temps <- seq(topt_4 - 4*standard_twidth, topt_4 + 4*standard_twidth, by=0.05)
optima_4_response <- gauss_rec(tx=all_temps, topt=topt_4, twidth=standard_twidth)
optima_4_viz <- data.frame(temp=all_temps, response=optima_4_response)
ggplot(optima_4_viz, aes(temp, response)) + geom_line(linewidth=2) + 
  geom_rug(data=north_bt, aes(x=temp_c), color="black", sides="b", alpha=0.25, inherit.aes = FALSE) + 
  geom_rug(data=all_proj_df, aes(x=temp_c, color=period), inherit.aes=FALSE) + 
  labs(x="Standardized temperature (deg C)", y="Response", color="Temperature trend") + 
  bioshift_plot_theme_1 + theme(panel.grid.major=element_line(color="grey", linewidth=0.1))

# Gather the response curves into a single plot
optima_1_viz <- optima_1_viz %>% mutate(optima="Optima - 1")
optima_2_viz <- optima_2_viz %>% mutate(optima="Optima - 2")
optima_3_viz <- optima_3_viz %>% mutate(optima="Optima - 3")
optima_4_viz <- optima_4_viz %>% mutate(optima="Optima - 4")
all_optima_viz <- rbind(optima_1_viz, optima_2_viz, optima_3_viz, optima_4_viz) # Gather all optima data
# all_optima_viz <- all_optima_viz %>% 
#   group_by(optima) %>% mutate(year=row_number(), .before=1) # Number the years
# The responses from each optima for each projection
# A dataframe that combines all the optima and the relevant widths
optima_gathered <- data.frame(optima=c("Optima - 1","Optima - 2","Optima - 3","Optima - 4"),
                              topt=c(topt_1, topt_2, topt_3, topt_4),
                              twidth=2)

trend_data <- data.frame(rate=c(scen_0_trend, scen_1_trend, scen_2_trend, scen_3_trend),
                         error=c(scen_0_sd, scen_1_sd, scen_2_sd, scen_3_sd),
                         period=c("Projection - 1", "Projection - 2", "Projection - 3", "Projection - 4"))

# Apply the function
g <- mapply(function(topt, tw) gauss_rec(all_proj_df$temp_c, topt, tw),
            optima_gathered$topt, optima_gathered$twidth)
colnames(g) <- optima_gathered$optima
all_proj_df_2 <- cbind(all_proj_df, g)
all_proj_df_2 <- all_proj_df_2 %>%
  pivot_longer(cols = starts_with("Optima"),
               names_to = "optima", values_to = "value")

# All the optima in a single plot
# We can make the trends_df nicer by introducing new column names and parsing them
# We are using the `ggpmisc` package here for rendering the optima information as a table
optima_gathered_2 <- optima_gathered
colnames(optima_gathered_2) <- c("Optima", "mu~(degree*C)", "sigma")
optima_plots <- ggplot(all_optima_viz, aes(temp, response, linetype=optima)) + geom_line(linewidth=1.5) + 
  geom_rug(data=north_bt, aes(x=temp_c), color="black", sides="b", alpha=0.25, inherit.aes = FALSE) + 
  geom_rug(data=all_proj_df, aes(x=temp_c, color=period), inherit.aes=FALSE) + 
  annotate(geom = "table", x = Inf, y = Inf,
           label = list(optima_gathered_2),
           hjust = 1, vjust = 1, parse=TRUE) + 
  labs(x="Standardized temperature (deg C)", y="Response", color="Temperature trend", linetype="Hypothetical Optima") + 
  bioshift_plot_theme_1 + 
  theme(panel.grid.major=element_line(color="grey", linewidth=0.1),
        legend.position="bottom")


# Break the optima into different panels and track the responses for the different trends
response_projections_plot <- ggplot(all_proj_df_2, aes(year, value, linetype=optima, color=period)) + 
  geom_line(linewidth=0.75) + 
  geom_label(data = trend_data,
             mapping = aes(x=Inf,
                           y=Inf,
                           label=paste("Trend ", round(rate,2), "±",error,"C/yr")),
             size=12, size.unit="pt",
             vjust=1, hjust=1, 
             fill = "white", alpha=0.5, label.size=NA,
             inherit.aes=FALSE) +
  facet_wrap(~period, ncol=1) + 
  labs(x="Year", y="Response", color="Temperature trend", linetype="Hypothetical Optima") + 
  bioshift_plot_theme_1 + 
  theme(panel.grid.major=element_line(color="grey", linewidth=0.1),
        legend.position = "none")

# Temperature trends/projections
temp_proj_plot <- ggplot(all_proj_df_2, aes(year, temp_c, color=period)) + 
  geom_line(linewidth=0.75) + 
  geom_line(data=north_bt, aes(year, temp_c), color="black", linewidth=0.50) + 
  geom_vline(xintercept=first_proj_year-1, linewidth=0.5, color="grey25", linetype=2) + 
  labs(x="Year", y="Temperature (deg C/yr)", color="Temperature trend") + 
  bioshift_plot_theme_1 + 
  theme(panel.grid.major=element_line(color="grey", linewidth=0.1),
        legend.position = "none")

# Make them into a nice plot using patchwork
# This requires some layout changes
layout <- "
AB
CC
" # This simply means that plots A,B are at the top, and the 2 'C's are the gathered legend
gathered_plot_1 <- (optima_plots | response_projections_plot) + guide_area() +
  plot_layout(design = layout, guides = "collect",
              heights = c(1, 0.15))


# Make them into a nice plot using patchwork
# This requires some layout changes
layout_2 <- "
AA
BC
DD
" 

# This simply means that plots A,B are at the top, and the 2 'C's are the gathered legend
gathered_plot_2 <- optima_plots + temp_proj_plot + response_projections_plot + guide_area() +
  plot_annotation(tag_levels = 'a', tag_suffix=')') + 
  plot_layout(design = layout_2, guides = "collect",
              heights = c(1, 1, 0.3)) & theme(legend.key.width = unit(2, "cm")) &
  guides(linetype = guide_legend(nrow = 1))

ggsave(here("images","temperature_trends_and_optima_gathered_2.png"), gathered_plot_2, height=10, width=10,
       units=c("in"), dpi=300)

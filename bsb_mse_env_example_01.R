# -----------------------------------------------------------------------------
#
#           Management Strategy Evaluation (MSE) Example for
#           Black Sea Bass (BSB) with Environmental Drivers
#                   BSB MSE - Environmental Drivers
#
#
# Author(s): Chengxue Li [Original], Jeewantha Bandara [This version]
# Date: 2025/08/06
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


library(wham)
library(here)
asap <- read_asap3_dat(here("data",c("north.dat","south.dat")))
temp <- prepare_wham_input(asap)

#basic_info
basic_info <- list(region_names = c("North", "South"), stock_names = paste0("BSB_", c("North", "South"))) 
seasons = c(rep(1,5),2,rep(1,5))/12
basic_info$fracyr_seasons <- seasons
#each age other than 1 (recruitment) for north stock can be in either region on Jan 1 
basic_info$NAA_where <- array(1, dim = c(2,2,8))
basic_info$NAA_where[1,2,1] = 0 #stock 1, age 1 can't be in region 2 
basic_info$NAA_where[2,1,] = 0 #stock 2, any age can't be in region 1 (stock 2 doesn't move) 
#average recruitment over years 2000+ for SSB40 BRPs
basic_info$XSPR_R_avg_yrs <- which(temp$years>1999)
basic_info$XSPR_R_opt <- 2 #use average of recruitments (random effects), not expected/predicted given last time step
#############################################
north_bt <- read.csv(here("data","bsb_bt_temp_nmab_1959-2022.csv"))
south_bt <- read.csv(here("data","bsb_bt_temp_smab_1959-2022.csv"))
# Cutting down to 1989 to 2014 for the temperature data
north_bt <- north_bt %>% filter(year>=1989)
south_bt <- south_bt %>% filter(year>=1989)
ecov <- list(label = c("North_BT","South_BT"))
ecov$mean <- cbind(north_bt[,'mean'], south_bt[,'mean'])
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean,2,mean))
ecov$logsigma <- log(cbind(north_bt[,'se'], south_bt[,'se']))
ecov$year <- north_bt[,'year']
ecov$use_obs <- matrix(1, NROW(ecov$mean),NCOL(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- apply(ecov$mean, 2, mean)

ecov_0 <- ecov 
ecov_0$recruitment_how <- matrix(c("none","none","none","none"), 2,2)
ecov_1 <- ecov_2 <- ecov_0
ecov_1$recruitment_how[1,1] <- "controlling-lag-0-linear" #north
ecov_2$recruitment_how[2,2] <- "controlling-lag-0-linear" #south
ecov_3 <- ecov_1
ecov_3$recruitment_how[2,2] <- "controlling-lag-0-linear" #both

ecov_4 <- ecov_5 <- ecov
ecov_4$M_how <- ecov_5$M_how <- array("none", dim = c(2,2,8,2)) #n_ecov, n_stocks, n_ages, n_regions
ecov_4$M_how[1,1,1,1] <- "lag-0-linear"

ecov_5$M_how[2,2,1,2] <- "lag-0-linear"
ecov_6 <- ecov_4
ecov_6$M_how[2,2,1,2] <- "lag-0-linear"

#############################################

#############################################
#NAA_re
NAA_re = list(sigma = list("rec+1","rec+1"), cor = list("2dar1","2dar1"), N1_model = rep("equilibrium",2))

#Wasn't decoupled in RT
#NAA_re$decouple_recruitment <- FALSE

# set fixed values for NAA re for stock 1 in south on Jan 1
NAA_re$sigma_vals <- array(1, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_vals[1,2,2:temp$data$n_ages] <- 0.05 #2+ fixed to SCAA-ish for North fish in south on Jan 1.
NAA_re$sigma_map <- array(NA, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_map[1,1:2,1] <- 1
NAA_re$sigma_map[2,2,1] <- 3
NAA_re$sigma_map[1,1,2:temp$data$n_ages] <- 2
NAA_re$sigma_map[2,2,2:temp$data$n_ages] <- 4


#turn off estimation of AR1 cor parameters for north population in the south on Jan 1
x <- array(NA, dim = dim(temp$par$trans_NAA_rho))
x[1,1,1:2] <- 1:2
x[2,,1] <- 3
x[2,,2] <- 4
NAA_re$cor_map <- x

x <- array(NA, dim = dim(temp$par$trans_NAA_rho))
x[1,1,1:3] <- 1:3
x[2,2,1:3] <- 4:6
NAA_re$cor_map <- x
# input_1$map$trans_NAA_rho <- factor(x)
#############################################

#############################################
#move
move = list(stock_move = c(TRUE,FALSE), separable = TRUE) #north moves, south doesn't
move$must_move = array(0,dim = c(2,length(seasons),2))	
#if north stock in region 2 (south) must move back to region 1 (north) at the end of interval 5 right before spawning
move$must_move[1,5,2] <- 1 
move$can_move = array(0, dim = c(2,length(seasons),2,2))
move$can_move[1,c(1:4),2,1] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,c(7:11),1,2] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,5,2,] <- 1 #north stock can (and must) move in last season prior to spawning back to north 
mus <- array(0, dim = c(2,length(seasons),2,1))
mus[1,1:11,1,1] <- 0.02214863 #see here("2023.RT.Runs","transform_SS_move_rates.R") for how these numbers are derived.
mus[1,1:11,2,1] <- 0.3130358
move$mean_vals <- mus 
move$mean_model = matrix("stock_constant", 2,1)
#prior distribution on movement parameters 
move$use_prior <- array(0, dim = c(2,length(seasons),2,1))
move$use_prior[1,1,1,1] <- 1
move$use_prior[1,1,2,1] <- 1
move$prior_sigma <- array(0, dim = c(2,length(seasons),2,1))
move$prior_sigma[1,1,1,1] <- 0.2
move$prior_sigma[1,1,2,1] <- 0.2
#############################################

#############################################
#selectivity
sel <- list(n_selblocks = 8,
            model = rep(c("age-specific","logistic","age-specific"),
                        c(2,2,4)))
sel$initial_pars <- list(
  rep(c(0.5,1),c(3,5)), #north comm
  rep(c(0.5,1),c(6,2)), #north rec
  c(5,1), #south comm
  c(5,1),	#south rec
  rep(c(0.5,1,1),c(1,1,6)), #north rec cpa
  rep(c(0.5,1),c(4,4)), #north vast
  rep(c(0.5,1,1),c(2,4,2)), #south rec cpa
  rep(c(0.5,1),c(1,7)) #south vast
)
sel$fix_pars <- list(
  4:8, #north comm
  7:8, #north rec
  NULL, #south comm
  NULL, #south rec
  2:8, #north rec cpa
  5:8, #north vast
  3:8, #south rec cpa
  2:8 #south vast
)
sel$re <- rep(c(
  "2dar1", "2dar1", "none", "ar1_y", "2dar1","none"), 
  c(1, 1, 2, 1, 1, 2))

#############################################

#############################################
#catch_info
catch_Neff <- temp$data$catch_Neff
catch_Neff[] <- 1000
catch_info <- list(catch_Neff = catch_Neff)
x <- temp$data$selblock_pointer_fleets
x[] <- rep(1:4, each = NROW(x))
catch_info$selblock_pointer_fleets = x
#############################################

#############################################
#index_info
index_Neff <- temp$data$index_Neff
index_Neff[] <- 1000
index_info <- list(index_Neff = index_Neff)
x <- temp$data$selblock_pointer_indices
x[] <- rep(5:8, each = NROW(x))
index_info$selblock_pointer_indices = x
#estimate obs sd scalar for Rec CPA in north and south
index_info$initial_index_sd_scale <- c(5,1,5,1)
index_info$map_index_sd_scale <- c(1,NA,2,NA)

#############################################

#############################################
#age_comp
age_comp = list(
  fleets = c("dir-mult","logistic-normal-miss0","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"), 
  indices = c("logistic-normal-miss0","dir-mult","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"))
#############################################

#############################################
#M random effects just on age 1
M_re_map <- array(NA, c(2,2,8))
M_re_map[1,1,1] <- 1
M_re_map[2,2,1] <- 2
M_sigma_vals <- diag(c(1,1))
M_sigma_map <- diag(1:2)
M_re_model <- matrix("none", 2,2)
M_re_model[1,1] <- "iid_y"
M_re_model[2,2] <- "iid_y"
M_list <- list(re_model = M_re_model, re_map = M_re_map, sigma_map = M_sigma_map, sigma_vals = M_sigma_vals)
#############################################

input_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp)
input_1$fleet_names <- paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names)
input_1$index_names <- paste0(rep(c("North_", "South_"),c(2,2)), temp$index_names)

fit_1 <- fit_wham(input_1, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

temp <- fit_1$input
temp$par <- fit_1$parList
temp$par$trans_NAA_rho[1,2,]
temp$par$trans_NAA_rho[2,1,] <- 0
x <- fit_wham(temp, do.brps = FALSE, do.fit = FALSE) #same
# fit_1 <- do_reference_points(fit_1)
# fit_1 <- do_sdreport(fit_1)

input_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_0, catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp)
# input_0$par <- fit_1$parList
# input_0$par$Ecov_beta_R[] <- 0
fit_0 <- fit_wham(input_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
#not quite as good as using better starting values
fit_0$fn(fits[[1]]$opt$par)
fit_0 <- do_retro_peels(fit_0, use.mle = FALSE)
# pkgbuild::compile_dll("c:/work/wham/wham",debug = FALSE)
# pkgload::load_all("c:/work/wham/wham")
xx <- retro(fit_0, n.peels = 1, use.mle = FALSE, MakeADFun.silent = TRUE, check.version = FALSE, save.input = TRUE)
x <- fit_peel(1, fit_0$input, n.newton = 0, MakeADFun.silent = TRUE, save.input = TRUE)
xx.mle <- retro(fit_0, n.peels = 1, use.mle = TRUE, MakeADFun.silent = TRUE, check.version = FALSE, save.input = TRUE)

temp <- fit_0$input
temp$par <- fit_0$parList
temp$map$logit_q <- factor(rep(NA, length(temp$par$logit_q)))
temp$map$trans_NAA_rho <- factor(rep(NA, length(temp$par$logit_q)))
x <- fit_peel(7, temp, n.newton = 0, MakeADFun.silent = FALSE, save.input = TRUE)

input_2 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_2, catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp)
input_2$par <- fit_0$parList
fit_2 <- fit_wham(input_2, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = TRUE)

input_3 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_3, catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp)
input_3$par <- fit_0$parList
fit_3 <- fit_wham(input_3, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = TRUE)

input_4 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_4, catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp)
input_4$par <- fit_0$parList
fit_4 <- fit_wham(input_4, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = TRUE)

input_5 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_5, catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp)
input_5$par <- fit_0$parList
fit_5 <- fit_wham(input_5, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = TRUE)

input_6 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_6, catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp)
input_6$par <- fit_0$parList
fit_6 <- fit_wham(input_6, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = TRUE)

# Save individually
saveRDS(fit_1, file = here("models", "fit_1.rds"))
saveRDS(fit_2, file = here("models", "fit_2.rds"))
saveRDS(fit_3, file = here("models", "fit_3.rds"))
saveRDS(fit_4, file = here("models", "fit_4.rds"))
saveRDS(fit_5, file = here("models", "fit_5.rds"))
saveRDS(fit_6, file = here("models", "fit_6.rds"))

sapply(0:6, function(x) {
  mod <- get(paste0("fit_",x))
  2*(mod$opt$obj + length(mod$opt$par))
})

input_M_re_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov, catch_info = catch_info,
                                   index_info = index_info, age_comp = age_comp, M = M_list)
nms <- names(input_M_re_0$par)
nms <- nms[which(nms != "M_repars")]
input_M_re_0$par[nms] <- fit_0$parList[nms]
nofit <- fit_wham(input_M_re_0, do.fit = FALSE, do.brps = FALSE)
x <- input_M_re_0$par$M_repars
x[] <- input_M_re_0$map$M_repars
fit_M_re_0 <- fit_wham(input_M_re_0, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

input_M_re_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info,
                                   index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_1$par[nms] <- fit_0$parList[nms]
fit_M_re_1 <- fit_wham(input_M_re_1, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

input_M_re_2 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_2, catch_info = catch_info,
                                   index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_2$par[nms] <- fit_0$parList[nms]
fit_M_re_2 <- fit_wham(input_M_re_2, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

input_M_re_3 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_3, catch_info = catch_info,
                                   index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_3$par[nms] <- fit_0$parList[nms]
fit_M_re_3 <- fit_wham(input_M_re_3, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

input_M_re_4 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_4, catch_info = catch_info,
                                   index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_4$par[nms] <- fit_0$parList[nms]
fit_M_re_4 <- fit_wham(input_M_re_4, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

input_M_re_5 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_5, catch_info = catch_info,
                                   index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_5$par[nms] <- fit_0$parList[nms]
fit_M_re_5 <- fit_wham(input_M_re_5, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

input_M_re_6 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_6, catch_info = catch_info,
                                   index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_6$par[nms] <- fit_0$parList[nms]
fit_M_re_6 <- fit_wham(input_M_re_6, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

sapply(0:6, function(x) {
  mod <- get(paste0("fit_M_re_",x))
  2*(mod$opt$obj + length(mod$opt$par))
})

for(i in 0:6) saveRDS(get(paste0("fit_", i)), here("results", paste0("fit_", i, ".RDS")))
for(i in 0:6) saveRDS(get(paste0("fit_M_re_", i)), here("results", paste0("fit_M_re_", i, ".RDS")))
sapply(Run34$input$map, function(x) sum(!is.na(x)))
sapply(nofit$input$map, function(x) sum(!is.na(x)))



fits <- lapply(0:6, function(i) readRDS(here("results", paste0("fit_",i,".RDS"))))

inpt <- fits[[1]]$input
inpt$par <- fits[[1]]$parList
fit_1_retro <- fit_wham(inpt, do.osa = FALSE, do.retro = FALSE)
fit_1_retro_peels <- do_retro_peels(fit_1_retro)

fit_1_retro_peels_alt <- do_retro_peels(fit_1_retro, use)

inpt$data$use_indices[tail(1:33,7),] <- 0
inpt$data$use_agg_catch[tail(1:33,7),] <- 0
inpt$data$use_catch_paa[tail(1:33,7),] <- 0
inpt$data$use_index_paa[tail(1:33,7),] <- 0
temp <- matrix(NA, 33, NCOL(inpt$par$F_pars))
temp[1:26,] <- 1:(26*NCOL(inpt$par$F_pars))
inpt$map$F_pars <- factor(temp)
temp <- fit_wham(inpt, do.osa=FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
inpt$par$F_pars[tail(1:33,7),] <- 0
temp <- fit_wham(inpt, do.osa=FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
inpt$map <- inpt$map[names(inpt$map) != "F_pars"]
inpt$random <- c(inpt$random, "F_pars")
temp <- fit_wham(inpt, do.osa=FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)

fits_M_re <- lapply(0:6, function(i) readRDS(here("results", paste0("fit_M_re_",i,".RDS"))))

#############################################
#M random effects just on age 1 for North stock only (South stock variance goes to 0)
M_re_map <- array(NA, c(2,2,8))
M_re_map[1,1,1] <- 1
M_re_map[2,2,1] <- 2
M_sigma_vals <- diag(c(1,1))
M_sigma_map <- diag(c(1,NA))
M_re_model <- matrix("none", 2,2)
M_re_model[1,1] <- "iid_y"
# M_re_model[2,2] <- "iid_y"
M_list <- list(re_model = M_re_model, re_map = M_re_map, sigma_map = M_sigma_map, sigma_vals = M_sigma_vals)
nms <- names(fits[[1]]$input$par)
nms <- nms[which(nms != "M_repars")]
fits_M_re <- lapply(0:6, function(i){
  input <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_",i)), catch_info = catch_info,
                              index_info = index_info, age_comp = age_comp, M = M_list)
  input$par[nms] <- fits[[1]]$parList[nms]
  return(fit_wham(input, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
})

fits_M_re <- lapply(fits_M_re, function(x) do_reference_points(x))
fits_M_re <- lapply(fits_M_re, do_sdreport)
fits_M_re <- lapply(fits_M_re, do_sdreport)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_better.RDS")))
fits_M_re <- lapply(fits_M_re, make_osa_residuals)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_better.RDS")))
fits_M_re <- lapply(fits_M_re, function(x) {
  y <- try(do_retro_peels(x))
  if(is.character(y)) return(x)
  else return(y)
})
sapply(fits_M_re, function(x) "peels" %in% names(x))

temp <- fits_M_re[[1]]$input
temp$data$q_upper[] <- 0.1
q <- fits_M_re[[1]]$rep$q[1,]
temp$par <- fits_M_re[[1]]$parList
temp$par$logit_q <- wham:::gen.logit(q, 0, 0.1)
temp <- fit_wham(temp, do.brps = TRUE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)
temp <- do_retro_peels(temp)
peels <- list()
for(i in 1:7){
  peels[[i]] <- try(fit_peel(i,input = temp$input, n.newton = 0))
}
x <- temp
x$peels <- peels
mohns_rho(x)$SSB


temp2 <- fits_M_re[[1]]$input
temp2$data$q_upper[] <- 0.1
q <- fits_M_re[[3]]$rep$q[1,]
temp2$par <- fits_M_re[[3]]$parList
temp2$par$logit_q <- wham:::gen.logit(q, 0, 0.1)
temp2 <- fit_wham(temp2, do.brps = TRUE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

gen.invlogit <- function(x,lo,hi,s){
  lo + (hi-lo)/(1+exp(-s*x))
}

gen.invlogit(-16.08 + 0.3, 0, 1000, 1)

sapply(fits_M_re, function(x) "peels" %in% names(x))


sapply(fits, function(x) {
  2*(x$opt$obj + length(x$opt$par))
})

sapply(fits_M_re, function(x) {
  2*(x$opt$obj + length(x$opt$par))
})

sapply(sort(grep("nll", names(nofit$rep), value = TRUE)), function(x) sum(nofit$rep[[x]]))
sapply(sort(grep("nll", names(Run34$rep), value = TRUE)), function(x) sum(Run34$rep[[x]]))

nms <- sort(unique(names(Run34$opt$par)))
x <- nofit$env$parList()
lapply(nms, function(i) x[[i]]- Run34$parList[[i]])


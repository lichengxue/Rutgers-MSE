# MSE trials
# Author - Stephen Potts
#####################################################################
  #####################################################################
  #####################################################################
  #####################################################################
  #####################################################################
  #####################################################################
  #####################################################################
  library(data.table)
  library(wham)
  library(whamMSE)
  library(here)
  source(here("scratch","potts-mse-functions.R"))
  #####################################################################
  #####################################################################
  
  #tjm "n_stock" not defined:
  n_stocks <- 1
  #### Beverton-Holt w pars a & b --> R=a*S*exp(beta*ecov)/(1+b*S)
  	recruit_pars <- rep(list(c(a=2, b=1*10^-4)), n_stocks)
  	l.NAA <- f.NAA(recruit_model=3, recruit_pars=recruit_pars)
  
  #### List of year variables
  	base.years=2001:2040
  	l.yr <- f.years(base.years=base.years, n_feedback_years=3, assess_interval=3)
  
  #### Basic info
  	l.info <- f.info(l.yr, n_ages=12)
  
  #### Ecov list
  	recruitment_how=matrix("controlling-lag-1-linear", nrow=1, ncol=1)
  	beta_R=0; mu_R=0; sig_R=0.5; phi_R=0.5 # beta, mean ecov, rec-ecov var, rec autocorrelation
  	l.ecov <- f.ecov(l.yr, recruitment_how=recruitment_how, beta_R=beta_R, mu_R=mu_R, sig_R=sig_R, phi_R=phi_R, ecov.form="logistic", seed.ecov=1)
  
  #### Mortality
  	M_naa <- c(0.4, rep(0.2, l.info$basic_info$n_ages-1)) #exp(-(2:l.info$basic_info$n_ages)*0.5))
  	l.M <- list(
  			model="constant",
  			initial_means=array(M_naa, dim=c(1, 1, 12)) # stocks, regions, ages
  			)
  
  #### OM input list
  	inputA <- prepare_wham_input(
  			basic_info=l.info$basic_info, 
  			M=l.M, 
  			NAA_re=l.NAA, 
  			ecov=l.ecov$l.ecov,
  			catch_info=l.info$catch_info, 
  			index_info=l.info$index_info, 
  			F=l.info$F
  			)
  			inputA$par$Ecov_beta_R[] <- beta_R # ecov-recruitment effect
  			inputA$par$Ecov_process_pars[] <- l.ecov$ecov_process_pars # R-ecov process parameters
  			inputA$par$Ecov_re[,1] <- l.ecov$ecov_re
  		# Tell WHAM NOT to simulate Ecov_re (use what we set above)
  			inputA$data$do_simulate_Ecov_re[] <- 0L
  		# Remove "ecov_re" from random effects so update_om_fn won't simulate it
  			random0 <- inputA$random
  			random0 <- random0[!random0 %in% "Ecov_re"]
  			inputA$random <- NULL
  		# Build OM
  		
  		#  Crank down variance of recruitment and older NAA and remove ecov effect to make sure things are working deterministicly
  	inputA$par$log_NAA_sigma[] <- log(0.2) # log(0.2) is the original value # log(0.00001)
  	inputA$par$Ecov_beta_R[] <- 0
  	om0_brps <- fit_wham(inputA, do.fit=FALSE, do.brps=TRUE, MakeADFun.silent=TRUE)
  
  
  # Set initial R to R at Fmsy and set F to Fmsy
  l.NAA$N1_pars <- array(c(mean(exp(om0_brps$rep$log_R_MSY)),mean(exp(om0_brps$rep$log_FMSY)),rep(NA,10)), dim = c(1,1,12))
  l.NAA$N1_model <- "equilibrium"
  l.NAA$sigma_vals <- array(0.2, dim = c(1,1,12))
  l.info$F$F[] <- mean(exp(om0_brps$rep$log_FMSY))
  inputA$par$log_NAA_sigma[] <- log(0.5)
  inputD <- inputA
  inputD$par$Ecov_beta_R[] <- -1
  inputD <- set_NAA(inputD, l.NAA)
  inputD <- set_F(inputD, l.info$F)
  # inputA$random <- inputB$random <- inputC$random <- inputD$random <- NULL
  inputA$random <- inputD$random <- NULL
  
  omD <- fit_wham(inputD, do.fit=FALSE, do.brps=FALSE, MakeADFun.silent=TRUE)
  
  ##### Run n.sim simulations and extract recruitment and SSB estimates
  	n.sim <- 10
  	m.rec <- m.ssb <- array(NA, dim=c(1, length(om0_brps$rep$log_F_tot),n.sim))
  		for(i in 1:n.sim) {
  			om_D <- update_om_fn(omD, random=random0, seed=i*sample(1:1e6, size=1))
  			m.rec[1,,i] <- om_D$rep$pred_NAA[,,,1]
  			m.ssb[1,,i] <- om_D$rep$SSB	
  		}
  
  
  	par(mfrow = c(2, 1))
  		ylim=c(0, max(m.rec))
  			plot(base.years, base.years, type="n", xlab="Year", ylab="Recruits", ylim=ylim, main=paste0("beta=-1"))
  				sapply(1:dim(m.rec)[3], function(x) lines(base.years, m.rec[1, 1:length(base.years),x]))
  		ylim=c(0, max(m.ssb))
  			plot(base.years, base.years, type="n", xlab="Year", ylab="SSB", ylim=ylim)
  				sapply(1:dim(m.ssb)[3], function(x) lines(base.years, m.ssb[1, 1:length(base.years),x]))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

#####################################################################
#####################################################################
#####################################################################
#####################################################################
###  FUNCTIONS FOR OM CREATION
# Author: Stephen Potts

### Generate various year vector list
f.years <- function(base.years, n_feedback_years, assess_interval) {
			l.yr <- list()
			l.yr$n_feedback_years <- n_feedback_years
			### Ecov
				l.yr$years_ecov_obs <- (min(base.years)-1):max(base.years) # ecov obs
				l.yr$years_ecov_proj <- max(base.years)+1:n_feedback_years # proj years
				l.yr$years_ecov <- min(l.yr$years_ecov_obs):max(l.yr$years_ecov_proj)
			### OM Parameters			
				l.yr$years_full <- min(base.years):(max(base.years)+n_feedback_years)
			### EM parameters
	   			l.yr$base.years=base.years
	   			l.yr$base_years=base.years
	   			l.yr$year.use   <- length(base.years)
	   			l.yr$terminal_year   <- max(base.years)
	   			l.yr$assess_interval <- assess_interval
				l.yr$assess_years <-seq(l.yr$terminal_year, l.yr$terminal_year  + n_feedback_years - assess_interval, by=assess_interval)
			return(l.yr)
		}
 
### Generate basic information list
f.info <- function(l.yr, n_ages) {
		info <- generate_basic_info(
				n_stocks=1,
				n_regions=1,
				n_indices=1,
				n_fleets =1,
				n_seasons=1,
				base.years=l.yr$base.years,
				n_feedback_years=l.yr$n_feedback_years,
				life_history="medium",
				n_ages=n_ages
				)
      	c.F <- rep(0.2, length(l.yr$base.years) + l.yr$n_feedback_years) # Constant F
		info$F$F[1:nrow(info$F$F),1] <- c.F # F for base years and proj
		return(info)
		}

### Configure NAA and recruitment parameters
# Recruitment --> Beverton-Holt = R(y+1) = a*S(y)/(1+b*S(y))
f.NAA <- function(recruit_model, recruit_pars) {
			### Variation NAA over stocks=1, regions=1, and ages=12
				sigma_naa_vals <- array(0.2, dim = c(1, 1, 12))
				sigma_naa_vals[1,1,1] <- 0.5
			### Stock-recuit relationship
				sigma_rec  <- "rec+1"
				re_cor_rec  <- "iid"
				ini.opt_rec <- "equilibrium"
			### Configure Numbers-at-Age (NAA) by stocks=1
				l.NAA <- list(
						N1_model=rep(ini.opt_rec, 1),
						recruit_pars=recruit_pars,
						sigma=rep(sigma_rec, 1),
						cor=rep(re_cor_rec, 1),
						recruit_model=recruit_model,
						sigma_vals=sigma_naa_vals
					)
			return(l.NAA)
		} 
	
### Build ecov parameters
# ecov.form defines shape of ecov_vals over years obs
f.ecov <- function(l.yr, recruitment_how, beta_R, mu_R, sig_R, phi_R, ecov.form, seed.ecov) {
			set.seed(seed.ecov)
			a.beta_R <- array(beta_R, dim=c(1, 1, 1), dimnames=list(stock="stock1", ecov="ecov_vals", poly="linear"))
   			ecov_process_pars <- matrix(c(mu_R, sig_R, phi_R), nrow=3, ncol=1, dimnames = list(c("mu", "log_sigma", "phi"), "ecov_vals"))
   			### Ecov values
	   			if(ecov.form=="sin") {
	   				c.t <- seq(-pi/4, 2*pi, length=length(l.yr$years_ecov_obs))
	   				c.val <- 2.5*sin(c.t) + 5
	   			}
	   			if(ecov.form=="logistic") {
	   				c.t <- seq(-5, 5, length=length(l.yr$years_ecov_obs))
	   				c.val <- exp(c.t)/(1+exp(c.t))
	   				c.val <- c.val-min(c.val)
	   				c.val <- 5*c.val/max(c.val) + 2.5
	   			}
	   			if(ecov.form=="linear") {
	   				c.val <- seq(2.5, 7.5, length=length(l.yr$years_ecov_obs))
	   			}
	   			ecov_vals_obs <- c.val + runif(c.val, -2.5, 2.5)
	   			ecov_vals_proj <- rep(mean(c.val), length(l.yr$years_ecov_proj))
	   			ecov_vals <- c(ecov_vals_obs, ecov_vals_proj)
	   			ecov_mat <- matrix(ecov_vals, ncol=1, dimnames=list(l.yr$years_ecov, "ecov_vals"))
   			### Ecov re
	    			ecov_re_obs <- scale(ecov_vals_obs-mu_R)
	    			ecov_re_proj <- rep(1.25, length(ecov_vals_proj))
	    			ecov_re <- c(ecov_re_obs, ecov_re_proj)
	  			l.ecov <- list(
	   					label="ecov_vals",
	   					year=l.yr$years_ecov,
	   					mean=ecov_mat,
	   					logsigma="est_1",
	   					use_obs=matrix(1, nrow=length(l.yr$years_ecov), ncol=1),
	   					process_model="ar1",
	   					recruitment_how=recruitment_how
	   					)
    			return(list(l.ecov=l.ecov, ecov_process_pars=ecov_process_pars, ecov_vals=ecov_vals, ecov_re=ecov_re))
		}


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
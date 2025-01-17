# ------------------------------------------------------------------------------
# Title: Functions for Weibull simulations to evaluate intercept calibration methods
# Author: Ryan Gan
# Date: 2024-06-17
# -----------------------------------------------------------------------------

# Libraries --------------------------------------------------------------------
library(survival)
library(tidyverse)

# Functions --------------------------------------------------------------------

# Function to solve for lambda given median time
# Can find lambda (scale param) for desired CDF proportion at a given time
log_lambda_from_median_time <- function(median_time, shape) {
  
  log( median_time / ( log(2)**(1/shape) ) )
  
}


# Weibull time random number generator ----
# NOTE: Can use rweibull too. 
weib_rng <- function(n = 100, shape, scale) {
  
  # Uniform prob
  u = runif(n, min = 0, max = 1)
  
  t = scale * (-log( 1 - u ))**(1 / shape) 
  
  return(t)
  
}
 
# Log-scale linear predictor given desired median time and shape ----
scale_lp <- function(median_time, shape, X = NULL, beta = NULL) {
  
  # log intercept based on desired median time
  intercept = log( median_time / ( log(2)**( 1/shape ) ) )
  
  # If no X matrix provided
  if ( is.null( X ) ) {
    
    log_scale = intercept
    
  } else {
    
    # if X matrix provided, construct linear term
    log_scale = intercept + as.matrix(X) %*% as.vector(beta) 
    
  }
  
  return( as.vector(log_scale) )
  
}


# Survival function for weibull ----
weib_surv <- function(t, shape, scale) {
  
  surv = exp( -(t / scale)**shape)
  
  return(surv)
  
}

# Weib hazard -----
weib_haz <- function(t, shape, scale) {
  
  haz = (shape / scale) * ( t / scale)**(shape - 1)
  
  return( haz )
  
}

# Cumulative hazard
weib_cumulhaz <- function(t, shape, scale) {
  
  cumulhaz = ( t / scale )**shape
  
  return(cumulhaz)
  
}


# Log hazard ----
weib_log_haz <- function(t, shape, scale) {
  
  log_haz = log(shape) - log(scale) + (shape - 1) * log(t / scale) 
  
  return( log_haz )
  
}

# Can find scale for desired CDF proportion at a given time ----
scale_from_cdf_solve <- function(t, shape, cdf_prop) {
  
  scale = t / ( ( -log(1 - cdf_prop) )**(1/shape) )
  
  return( scale )
  
}


# follow up times functions
# Can be used to simulate a more regular follow-up time that we see in trial
fup_uniform_mix <- function(n, end_of_study, enrollment_time, prob_outside_enrollment = 0.3) {
  # indicator if follow up outside enrollment period (add 1 so it's 1 or 2)
  lof_ind = rbinom(n, 1, prob = prob_outside_enrollment) + 1
  
  fup = runif(
    n, 
    min = c(end_of_study - enrollment_time, 1)[lof_ind], 
    max = c(end_of_study, end_of_study - enrollment_time - 1)[lof_ind]
    )
  
  return( fup )

}

# Bounded rexp 
fup_bound_rexp <- function(n, rate, upper_bound) {
  
  fup = pmin( rexp(n = n, rate = rate), upper_bound )
  
  return( fup )
  
}


# function for approx cox hr
beta_hr_weib <- function(hr, shape) {
  
  beta = -log( hr ^ ( 1 / shape) )
  
  return( beta )
  
}


# TODO: Fix to make scale and shape consistently named
# scale should be shape and vise versa here to be consistent with 
# common weibull parameterization, not survreg

# Weibull estimation functions
weib_median_time <- function(linpred, scale, intercept_bias = 0, shape_bias = 0) {
  
  shape = exp( sum(linpred) + intercept_bias )
  # Try calibrating the shape as well by adjusting for the bias in that parameter too
  scale_adj = scale + shape_bias 
  
  # Estimate of median time
  res = shape * log(2)**(1/scale_adj)
  
  return( res )
  
}

# Survival curve function
weib_surv_curv <- function(
  linpred, scale, intercept_bias = 0, shape_bias = 0, 
  time_vector = seq(0,80, by = 0.1)
  ) {
  
  # Calculate shape based on linear predictors provided
  shape = exp( sum(linpred) + intercept_bias ) 
  scale_adj = scale + shape_bias
  # Calculate survival probabilities
  surv = exp( -( time_vector / shape )**scale_adj )
  
  # Create a data frame with time and survival probabilities
  res = data.frame( time = time_vector, surv = surv )
  
  return( res  )
  
}

# Simulate treatment assignment 
trt_sim <- function(n, confounding = TRUE, seed = NULL) {
  
  if( !is.null(seed) ) { 
    set.seed(seed) 
  }
  
  # Simulate a high risk category
  hrisk = rbinom(n = n, size = 1, prob = 0.2)
  # Simulate age and scale it
  age = scale( rnorm(n = n, mean = 73, sd = 7)  )
    
  if( confounding ) {  
    
    # Simulate confounded treatment by age and hrisk
    p = as.numeric( 1 / (1 + exp( -( hrisk * log(0.5) + age * log(0.7) ) ) ) )
  
  } else {
    
    # If no confounding, simulate p for treatment assignment to be equal chance 
    p = 0.5
      
    }
  
  
  # Treatment assignment based on p
  trt = rbinom(n = n, size = 1, prob = p)
  
  # Covariates
  covariates = data.frame(trt, hrisk, age)
  
  return( covariates )
  
}


# Simulate survival
# TODO: Add a more flexible way to add covariates
surv_sim <- function(
  n = 365*2, mpfs = 34.2, weib_shape = 1, 
  max_time = 42, enroll_time = 12, cnsr_prob = 0.2,
  X = NULL, hr_vec, seed = NULL) {
  
  if( !is.null(seed) ) { 
    set.seed(seed) 
    }
  
  # Linear predictor of scale based on desired mpfs, X, and beta
  if( !is.null(X) ) {
    lp = scale_lp(
      median_time = mpfs, 
      shape = weib_shape,
      X = X,
      beta = beta_hr_weib(hr = hr_vec, shape = weib_shape)
    )
  } else {
    # If no covariate X matrix provided
    lp = scale_lp(
      median_time = mpfs, 
      shape = weib_shape,
      X = NULL
    )

  }

  # simulated true times
  true_times = weib_rng(n = n, shape = weib_shape, scale = exp(lp))
  
  # simulated fup time
  fup = fup_uniform_mix(
    n = n, end_of_study = max_time, enrollment_time = enroll_time,
    prob_outside_enrollment = cnsr_prob 
  )
  
  # truncate time based on fup
  times = pmin(true_times, fup)
  
  # event status
  event = as.integer(times < fup)

  # dataframe
  data = data.frame( true_time = true_times, time = times, event = event, fup = fup, X)
  # order
  data = data[order(data$trt), ] 
  # reindex
  rownames(data) = 1:nrow(data)
  
  return(data)
  
}

# TODO: Add covariate relationship to misclassification
# Add misclassification to survival time
# 1 - sp = fp rate, which subtracts time
# 1 - sn = fn rate, which adds time
# which group to apply it to indicator
add_misclassification <- function(
  data, sn = 1, sn_time_shift = c(log(3*30), 0.5), 
  sp = 1, sp_time_shift = c(log(2*30), 0.5), trt_group = 0
) {
  

  # new event indicator and time
  data$event_obs = data$event
  data$time_obs = data$time
  
  
  # treatment index
  trt_idx = which(data$trt %in% trt_group)
  
  assertthat::assert_that(length(trt_idx) == nrow(data[data$trt==0,]))
  # event index
  event_idx = which(data[trt_idx, "event"] == 1)
  # no event index
  no_event_idx = which(data[trt_idx, "event"] == 0)
  
  
  # Reclassify event as missed and add time
  # false negative index
  fn_idx = sample( event_idx, size = round(length(event_idx)*(1 - sn), 0) )
  data[fn_idx, "time_obs"] = data[fn_idx, "time"] + rlnorm(length(fn_idx), meanlog = sn_time_shift[1], sdlog = sn_time_shift[2])
  # If added time is still captured within follow-up period, include as event, else 0
  data[fn_idx, "event_obs"] = 0
  data[fn_idx, "event_obs"] = as.integer( data[fn_idx, "time_obs"] < data[fn_idx, "fup"] )
  
  # Reclassify non-event as event and subtract time
  fp_idx = sample( no_event_idx, size = round(length(no_event_idx)*(1 - sp), 0) )
  data[fp_idx, "event_obs"] = 1
  data[fp_idx, "time_obs"] = pmax(0.1, data[fp_idx, "time"] - rlnorm(length(fp_idx), meanlog = sp_time_shift[1], sdlog = sp_time_shift[2]))

  return(data)

}

# Add assessment schedule
# Helper functions for assessment schedule 

# TODO: Should add tests and asserts to make sure it's the right group based on the way I'm setting this up
# MAIA das
add_maia_das <- function(
    data, 
    trt_group = 1,
    time_var = "time_obs",
    das_time_chg = 24, 
    das = c(1,2)
    ) {
  
  # Find treatment indicator
  trt_idx = which(data$trt %in% trt_group)
  
  # Pull out time vector by treatment idx
  time_vector = data[trt_idx, time_var]
  
  assertthat::assert_that(
    !is.null(time_vector),
    msg = "A vector of time must be provided"
  )
  
  # If less than 2 years
  right_time = sapply(time_vector, function(x) {
    
    if( x <= das_time_chg ) {
      
      ceiling(x/das[1])*das[1]
    
    } else {
      
      ceiling(x/das[2])*das[2]
      
      } 
    }
  )
  
  left_time = sapply(right_time, function(x) {
    
    if( x <= das_time_chg + das[1] ) {
      
      pmax(0.1, x - das[1])
      
    } else {
      
      pmax(0.1, x - das[2])
      
      } 
    }
  )
  
  # add back in left and right time
  data[trt_idx, "left_time"] = left_time
  data[trt_idx, "right_time"] = right_time
  # replace time_obs with the observed time under given assessment schedule
  data[trt_idx, "time_obs"] = right_time
  
  return( data )
  
}


# RWD uniform mixture
add_rwd_das <- function(
    data, 
    trt_group = 0,
    time_var = "time_obs",
    on_cycle_prob = 0.3, cycle_time = 1, min_time = 0.1, max_time = 3.3
    ) {
  
  # Find treatment indicator
  trt_idx = which(data$trt %in% trt_group)
  
  # Pull out time vector by treatment idx
  time_vector = data[trt_idx, time_var]
  
  assertthat::assert_that(
    !is.null(time_vector),
    msg = "A vector of time must be provided"
    )
  
  # Get number of observations
  nobs = length(time_vector)
  
  # Generate two indicator if on cycle or off cycle
  left_ind = rbinom(n = nobs, size = 1, prob = on_cycle_prob)
  right_ind = rbinom(n = nobs, size = 1, prob = on_cycle_prob)
  
  # Find disease assessment times (DAS)
  right_das = sapply(right_ind, function(x) {
    
    if( x == 1 ) {
      
      cycle_time
    
      } else {
      
      runif(1, min = min_time, max = max_time)
      
      }
  
    }
    )
  
  left_das = sapply(left_ind, function(x) {
    
    if( x == 1 ) {
      
      cycle_time
      
    } else {
      
      runif(1, min = min_time, max = max_time)
      
    }
    
  }
  )
  
  # Times based on assessment schedules
  right_time = ceiling(time_vector/right_das)*right_das
  # Time before final assessment, bounded on the low end by 0.1
  left_time = pmax(0.1, right_time - left_das)
  
  # add back in left and right time
  data[trt_idx, "left_time"] = left_time
  data[trt_idx, "right_time"] = right_time
  # replace time_obs with the observed time under given assessment schedule
  data[trt_idx, "time_obs"] = right_time
  
  return( data )
  
}

# Function to add generic assessment schedule
add_das <- function(
    data, 
    trt_group = c(0,1),
    time_var = "time_obs",
    das = 1
) {
  
  # Find treatment indicator
  trt_idx = which(data$trt %in% trt_group)
  
  # Pull out time vector by treatment idx
  time_vector = data[trt_idx, time_var]
  
  assertthat::assert_that(
    !is.null(time_vector),
    msg = "A vector of time must be provided"
  )
  
  right_time = ceiling(time_vector/das)*das
    
  left_time = pmax(0.1, right_time - das)
      
  # add back in left and right time
  data[trt_idx, "left_time"] = left_time
  data[trt_idx, "right_time"] = right_time
  # replace time_obs with the observed time under given assessment schedule
  data[trt_idx, "time_obs"] = right_time
  
  return( data )
  
}

# Function to set up data structure for interval censoring
prep_interval_cnsr_time <- function(data) {
  
  # Check necessary columns are there
  assertthat::assert_that(
    all(c("left_time", "right_time") %in% colnames(data)), 
    msg = "'left_time' and/or 'right_time' columns not present in data. Run a das function first."
    )
  
  # Check there is no missing data
  assertthat::assert_that(
    all(!is.na(data[, c("left_time", "right_time", "event_obs")])), 
    msg = "Missing data in either event or time columns."
    )
  
  # Set up new variables for interval censoring
  # NOTE: Keeping original right time to demonstrate results when interval censoring isn't accounted for
  data$left_time_ic = pmax(0.1, data$left_time)
  data$right_time_ic = data$right_time
  
  # Right time when events == 0 should be set to Inf
  data[data$event_obs == 0, "right_time_ic"] = Inf
  
  return( data )

  }

# TODO: See if I can break this up 
# Simulation bias evaluation 
sim_bias_evaluation <- function(
  data, trt0_val_prop = 0.3, iteration = NULL, 
  interval_censoring = TRUE, att_wt_adj = TRUE) {
 
  # Check necessary columns are present
  assertthat::assert_that(
    any(c("time", "event", "time_obs", "event_obs") %in% colnames(data)), 
    msg = "Required true and observed time and event columns not in dataframe"
  )

  # Set up some formulas for interval censoring
  if( interval_censoring ) {
    
    assertthat::assert_that(
      any(c("left_time_ic", "right_time_ic") %in% colnames(data)), 
      msg = "Required left and right times for interval censored not in dataframe"
    )

    assertthat::assert_that(
      any(c(data[, "left_time_ic"] < data[,"right_time_ic"])), 
      msg = "left_time_ic > right_time_ic"
    )
    
    # Naive KM estimator under interval censoring option
    km_fml = as.formula("Surv(right_time, event_obs) ~ trt")
    
    # Weibull interval censored option
    weib_fml = as.formula("Surv(left_time_ic, right_time_ic, type = 'interval2') ~ trt")
    
    # Weibull interval censored for validation
    weib_val_fml = as.formula("Surv(left_time_ic, right_time_ic, type = 'interval2') ~ 1")
    
    # Weibull intercept calibration approach
    weib_intcal_fml = as.formula("Surv(left_time_ic, right_time_ic, type = 'interval2') ~ 1")
    
  } else {
    
    km_fml = as.formula("Surv(time_obs, event_obs) ~ trt")
    
    # Weibull right censored option
    weib_fml = as.formula("Surv(time_obs, event_obs) ~ trt")
    
    # Weibull right censored validation 
    weib_val_fml = as.formula("Surv(time_obs, event_obs) ~ 1")

    # Weibull intercept calibration approach
    weib_intcal_fml = as.formula("Surv(time_obs, event_obs) ~ 1")
    
  }
 
  # ATT Weighting for age and hrisk confounder
  trt_mod = glm(trt ~ age + hrisk , data = data, family = binomial(link = "logit"))
  
  # pscore
  pscore = as.numeric( predict(trt_mod, type = "response") )
  
  # ATT weight
  att_wt = data$trt + (1 - data$trt) * (pscore / (1 - pscore))
  
  # Add ATT weight to dataframe if true
  if( att_wt_adj ) {
   
     data$att_wt = att_wt
     
  } else {
    
    # otherwise just leave as 1 for all
    data$att_wt = 1
    
  }
  # TODO: Could add something like effective smaple size test here??
  # I feel like something is off here.
  
  # "True" time with att weight
  km_true = survfit(Surv(time, event) ~ trt, data = data, weights = att_wt)
  
  km_true_mpfs_trt1 = summary(km_true)$table[2, "median"]
  
  km_true_mpfs_trt0 = summary(km_true)$table[1, "median"]
  
  km_true_mpfs_diff = as.numeric( diff( summary(km_true)$table[, "median"] , lag = 1) )
  
  # Naive KM estimator
  km = survfit(km_fml, data = data)
  
  km_unadj_mpfs_trt1 = summary(km)$table[2, "median"]
  
  km_unadj_mpfs_trt0 = summary(km)$table[1, "median"]
  
  km_unadj_mpfs_diff = as.numeric( diff( summary(km)$table[, "median"] , lag = 1) )
  
  # ATT weighted KM estimator
  km_att = survfit(km_fml, data = data, weights = att_wt)
  
  km_att_mpfs_trt1 = summary(km_att)$table[2, "median"]
  
  km_att_mpfs_trt0 = summary(km_att)$table[1, "median"]
  
  km_att_mpfs_diff = as.numeric( diff( summary(km_att)$table[, "median"] , lag = 1) )
  
  # Validation sample in the trt=0 arm
  
  # Treatment index
  trt_idx = which(data$trt == 0)
  
  val_idx = sample(trt_idx, size = round(length(trt_idx) * trt0_val_prop, 0) )
  val_data = data[val_idx,]
  
  # Validation models
  weib_true = survreg(Surv(time, event) ~ 1, data = val_data, dist = "weibull", weights = att_wt)
  weib_obs = survreg(weib_val_fml, data = val_data, dist = "weibull", weights = att_wt)

  # intercept_bias
  intercept_bias = as.numeric( weib_true$coefficients["(Intercept)"] - weib_obs$coefficients["(Intercept)"] )

  # Full observed model, unweighted
  weib_unadj = survreg(weib_fml, data = data, dist = "weibull")
  
  # Estimate median time
  mpfs_trt0_unadj = weib_median_time(linpred = weib_unadj$coefficients[1], scale = 1/weib_unadj$scale, intercept_bias = 0)
  
  mpfs_trt1_unadj = weib_median_time(linpred = weib_unadj$coefficients[1:2], scale = 1/weib_unadj$scale, intercept_bias = 0)
  
  weib_unadj_mpfs_diff = mpfs_trt1_unadj - mpfs_trt0_unadj

  # Full observed model, ATT weighted for confounders
  # Treat 0 only
  weib_trt0_att = survreg(weib_intcal_fml, data = data[trt_idx, ], dist = "weibull", weights = att_wt)

  weib_trt1_att = survreg(weib_intcal_fml, data = data[-trt_idx, ], dist = "weibull", weights = att_wt)


  # Estimate median time
  mpfs_trt0_att = weib_median_time(linpred = weib_trt0_att$coefficients[1], scale = 1/weib_trt0_att$scale, intercept_bias = 0)
  
  mpfs_trt1_att = weib_median_time(linpred = weib_trt1_att$coefficients[1], scale = 1/weib_trt1_att$scale, intercept_bias = 0)
  
  weib_att_mpfs_diff = mpfs_trt1_att - mpfs_trt0_att

  # Estimated median time accounting for bias
  
  mpfs_trt0_att_intcal = weib_median_time(linpred = weib_trt0_att$coefficients[1], scale = 1/weib_trt0_att$scale, intercept_bias = intercept_bias )
  
  mpfs_trt1_att_intcal = weib_median_time(linpred = weib_trt1_att$coefficients[1], scale = 1/weib_trt1_att$scale, intercept_bias = 0)
  
  weib_att_intcal_mpfs_diff = mpfs_trt1_att_intcal - mpfs_trt0_att_intcal
  
  # Results
  mpfs_res = data.frame(
    method = c(
      "KM True", "KM unadjusted", "KM ATT", "Weibull unadjusted", 
      "Weibull ATT", "Weibull ATT Intercept Calibration"
    ), 
    
    mpfs_trt1 = c(km_true_mpfs_trt1, km_unadj_mpfs_trt1, km_att_mpfs_trt1, 
                  mpfs_trt1_unadj, mpfs_trt1_att, mpfs_trt1_att_intcal),
    
    mpfs_trt0 = c(km_true_mpfs_trt0, km_unadj_mpfs_trt0, km_att_mpfs_trt0, 
                  mpfs_trt0_unadj, mpfs_trt0_att, mpfs_trt0_att_intcal),
    
    mpfs_diff = c(km_true_mpfs_diff, km_unadj_mpfs_diff, km_att_mpfs_diff, 
                  weib_unadj_mpfs_diff, weib_att_mpfs_diff, weib_att_intcal_mpfs_diff), 
    
    iteration = ifelse(is.null(iteration), 1, iteration)
  )
  
  return( mpfs_res )
  
}

# ------------------------------------------------------------------------------
# Title: Weighting methods vs standaridation 
# Author: Ryan Gan
# Date: 2025-01-10
# Purpose: Trying to make sense when direct standardization should align with weighting
# for event rate and survival modeling.
# ------------------------------------------------------------------------------

# Summary:
# ATT weighting, where white patients are weighted to be similar to minority patients
# based on categorical age and sex, should produce similar results to direct
# standardization approaches for rates. This would also align with weighted
# Kaplan Meier method for describing survival, vs. the direct standardization
# approach in the direct standardization paper produces a conditional survival curve.
# Both approaches should work, but I like the weights since it's easier to report
# a marginal rate/survival curve + makes computing confidence intervals possible.

# TODO: 
# - Check another text if my epi standardization calcs are right.

# Packages ---------------------------------------------------------------------
library(survival)
library(survey)
library(tidyverse)

# Helper functions -------------------------------------------------------------
# invlogit function
inv_logit <- function(x) {1 / (1 + exp(-x)) }

# Effective sample size
ess <- function(w) { sum(w)^2/sum(w^2) }

# Simulated data ---------------------------------------------------------------
# Seed
set.seed(555)
# Simulate some data
N <- 10000
# age category >= 65
age65 <- rbinom(n = N, size = 1, prob = 0.3)
# effect of over 65 increases rate by 2
b_age65 <- log(2)
# sex male
sexm <- rbinom(n = N, size = 1, prob = 0.5)
# no sex difference
b_sexm <- log(1)

# race minority related to probability of being sampled/in sample
# older age and males are less likely to be within the sample
mu <- inv_logit( age65 * log(0.2) + sexm * log(0.4) )

# create vector based on prob of mu
minority <- rbinom(n = N, size = 1, prob = mu)

# minoriy should be no different
b_minority <- log(1.0)

# Base event rate (intercept) based on median time of 24 months
intercept <- log( -log(1 - 0.5) / 24 )

# log-lambda
log_lambda <- intercept +  b_age65 * age65 +  b_sexm * sexm + b_minority * minority
lambda <- exp(log_lambda)

# Simulate an event time based on exponential

full_times <- rexp(n = N, rate = lambda)

# Fixed follow-up of 48 months
fup <- 48

# binary event occurs during 48 month period
events <- as.integer( full_times <= fup )
# time observed
time <- pmin(full_times, fup)


# data
data <- data.frame(time, events, age65, sexm, minority) |> 
  mutate(
    age65 = factor(ifelse(age65 == 1, ">=65", "<65"), levels = c("<65", ">=65")),
    sexm = factor(ifelse(sexm == 1, "Male", "Female"), levels = c("Female", "Male")),
    minority = factor(ifelse(minority == 1, "Minority", "White"), levels = c("White", "Minority"))
  )


# check model
summary( survreg(Surv(time, events) ~ minority , data = data, dist = "weibull") )
# check fully parameterized model
summary( survreg(Surv(time, events) ~ minority + age65 + sexm, data = data, dist = "weibull") )

# Estimates of rates via standardization ---------------------------------------

# Crude Marginal rates
marginal_rates <- data |> 
  group_by(minority) |> 
  summarize(
    n = n(),
    total_time = sum(time), 
    events = sum(events)
  ) |> 
  ungroup() |> 
  mutate(
    rate = events / total_time,
    # rate per 100 person months
    rate_per100 = rate * 100
  )

# calculate rates for each group
strata_rates <- data |> 
  group_by(minority, age65, sexm) |> 
  summarize(
    n = n(),
    total_time = sum(time), 
    events = sum(events)
  ) |> 
  group_by(minority) |>
  mutate(
    rate = events / total_time, 
    # rate per 100 person months
    rate_per100 = rate * 100, 
    event_prop = events / n
  ) |> 
  ungroup()

print(strata_rates)

# standardized rates, using minority person time as the population to standardize white rate to
# working with a row indicator to identify minority or -minority (white) rows
minority_idx <- which(strata_rates$minority == "Minority")

# index to white rate and multiply by minority person time denomiator over sum of minority person time denominator
white_std_rates <- sum( strata_rates[-minority_idx, "rate"] * strata_rates[minority_idx, "total_time"] ) / sum(strata_rates[minority_idx, "total_time"])
minority_rates <- sum(strata_rates[minority_idx, "events"]) / sum(strata_rates[minority_idx, "total_time"])


# combine with marginal rates
marginal_rates$std_rates_per100 <- c(white_std_rates, minority_rates) * 100

# Print crude and standardize rates
print( marginal_rates )

# Model based approach ---------------------------------------------------------
# Attempt to model the rate using a poisson model with a log offset on time
# remember to log time to get the rate since poisson uses log link. Had some
# conversion issues before I remembered this :)

unadj_pois_mod <- glm(
  events ~ minority + offset(log(time)), 
  data = data, 
  family = poisson()
)

summary(unadj_pois_mod)

# unadjusted rate estimates
unadj_pois_est <- predict(
  unadj_pois_mod, 
  type = "response", 
  newdata = data.frame(minority = c("White", "Minority"), time = 100) # 1 month to get rate, or 100 month to get per 100
)


# propensity score model
pscore_mod <- glm(
  minority ~ age65 + sexm + age65 * sexm,
  data = data, 
  family = binomial()
)

# pscore estimate for each subject
pscore <- predict(pscore_mod, type = "response")

# att weight
att_wt <- ( ( pscore * minority ) / pscore ) + ( pscore * (1 - minority) / (1 - pscore) )

# Set index
idx <- which(data$minority == "Minority")

# effective sample size after weighting
data.frame(
  Minority = c("Minority", "White"), 
  N = c(length( att_wt[idx]), length( att_wt[-idx] ) ), 
  ESS = c( ess(att_wt[idx]), ess(att_wt[-idx]) )
)

# set up svy design matrix (important for SE)
att_design <- survey::svydesign(ids = ~1, weights = att_wt, data = data)

# att weighted rate estimates
att_pois_mod <- survey::svyglm(
  events ~ minority + offset(log(time)), 
  design = att_design,
  family = quasipoisson()
)

summary(att_pois_mod)
# att rate estimates
att_pois_est <- predict(
  att_pois_mod, 
  type = "response", 
  newdata = data.frame(minority = c("White", "Minority"), time = 1) # 1 month to get rate, or 100 month to get per 100
)

# Compare results via different methods
final_marginal_rates <- cbind(marginal_rates, unadj_pois_est, att_pois_est = att_pois_est[1:2] * 100)

print( final_marginal_rates )

# Note, weighting gets closer to the simulated "no effect between race groups"
# but both direct standardization and weighting work well. 
# Benefit of weighting approach might be that it's straightforward to estimate
# confidence intervals. Just be sure to appropriately account for weights 
# e.g. survey package in R

# Survival Evaluation of 'Direct Standardization' and Weighting ----------------

# My understanding is that the direct standardization method proposed in the 
# SAP is an adjusted cox model predicting at specific levels

km <- survfit(Surv(time, events) ~ minority, data = data)

plot( km , col = c("darkred", "darkblue"))
legend("topright", legend = c("White", "Minority"), col = c("darkred", "darkblue"), lty = 1)

max_time <- max(km$time)
  
km_summary <- summary(km, times = seq(from = 0, to = max_time, by = 1))


km_table <- data.frame(
  strata = km_summary$strata,
  time = km_summary$time, 
  surv = km_summary$surv,
  cdf = 1 - km_summary$surv,
  cumhaz = km_summary$cumhaz
)

km_table |> 
  group_by(strata) |> 
  do(tail(., 1))


# TODO: Would be helpful to chat with Laura and Brett about this
# I think cumulative hazard for the last time point should be equal to the 
# rate * last time period (assume closed cohort, e.g. no censoring except at final time)
( final_marginal_rates$unadj_pois_est / 100 ) * 48

# direct adjusted cox model
cox_mod <- coxph(Surv(time, events) ~ minority + age65 + sexm, data = data)

summary(cox_mod)

# Not sure how you get the marginal here..., this is my issue with the 
# direct adjustment method. Instead I think it's a conditional estimate....
cox_est_wht_lt65_f <- survfit(
  cox_mod, 
  newdata = data.frame(minority = "White", age65 = "<65", sexm = "Female")
)

cox_est_min_lt65_f <- survfit(
  cox_mod, 
  newdata = data.frame(minority = "White", age65 = "<65", sexm = "Female")
)

# Compare conditional median time for white and minority under 65 and female
# should be approximately the same the same. This is the direct adjustment method
all.equal( 
  summary(cox_est_min_lt65_f)$table["median"], 
  summary(cox_est_wht_lt65_f)$table["median"]
)

# Weighted KM
# using weights in km, with note that variance may not be the same as the more appropriate
# survey package, but I imagine it's close enough. Can also bootstrap

km_att <- survfit(Surv(time, events) ~ minority, data = data, weights = att_wt)

km_att

# Survey KM -----
# Weighting method, using survey package, this can take a bit to run for some reason...
# might need to use computational approaches to approximate the CI?
# This approach is too slow; commenting out
# Although if you want survey type variance (which we usually do, use this approach)
# but note per documentation, with se=TRUE, it's Aalen-Johansen estimate and not product limit

km_svy_att <- survey::svykm(
  Surv(time, events) ~ minority, 
  design = att_design, 
  #se = TRUE # commented out since it bugs on my positron ide
)

km_svy_att

# Note, same median estimates as weighted KM, but CIs need to be estimated

# overlap
plot(km_svy_att, col = c("darkred", "darkblue"))

# km_att_wt
# quantile(km_att_wt, ci=TRUE)
# plot(km_att_wt)
# print(km_att_wt)



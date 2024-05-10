# ------------------------------------------------------------------------------
# Title: KM estimator calculation
# Author: Ryan Gan
# Date: 2024-05-10
# Purpose: Coding up calculation of KM estimator for educational purposes
# ------------------------------------------------------------------------------


# survival library
library(survival)


# using lung dataset

# sorted unique times
times <- unique(sort(lung$time))
event_times <- sort(lung[which(lung$status == 2), "time"])
censor_times <- sort(lung[which(lung$status == 1), "time"])

N <- length(times)

# initial n at risk
n.risk.init <- nrow(lung)

# empty vector of survival
surv <- rep(1, N)
n.event <- rep(0, N)
# n.cesnor
n.censor <- rep(0, N)
# n at risk
n.risk <- rep(n.risk.init, N)


# for loop to calcuate metrics
for ( i in 1:N) {

    if (i == 1) {
        
        # event during interval
        n.event[i] <- sum(event_times <= times[i]) 
        # censor during interval
        n.censor[i] <- sum(censor_times <= times[i]) 

        # find at risk subset
        n.risk[i] <- n.risk.init

        # find surv prop
        surv[i] <- (1 - (n.event[i] / n.risk[i])) 

    } else {

        # event during interval
        n.event[i] <- sum(event_times > times[i-1] & event_times <= times[i]) 
        # censor during interval
        n.censor[i] <- sum(censor_times > times[i-1] & censor_times <= times[i])

        # find at risk subset
        n.risk[i] <- n.risk[i-1] - (n.event[i-1] + n.censor[i-1])

        # surv conditioned on previous surv (cumulative product)
        surv[i] <- (1 - (n.event[i] / n.risk[i])) * surv[i-1]


    }


}

# variance
#var = S(t)^2 * sum( di / (ni*(ni - di)) )

# SE = sqrt(var)
se <- sqrt( surv^2 * cumsum(n.event / (n.risk*(n.risk - n.event)))  )
# normal approx to conf
upper.95 <- surv^exp( (qnorm(0.975)*se)/surv/log(surv) )
lower.95 <- surv^exp( (qnorm(0.025)*se)/surv/log(surv) )

# surv dataframe
surv_df <- data.frame(times, n.risk, n.event, n.censor, surv, se, lower.95, upper.95)


# KM fit
km.fit <- survfit(Surv(time, status) ~ 1, data = lung)

# res <- summary(km.fit)

# Plot to check output of km.fit vs calculated. looks pretty good.
# plot km fit
plot(km.fit)
# plot calculated surv
lines(times, surv, add = TRUE, type = "s", col = "red")
lines(times, lower.95, add = TRUE, type = "s", col = "red", lty = 2)
lines(times, upper.95, add = TRUE, type = "s", col = "red", lty = 2)


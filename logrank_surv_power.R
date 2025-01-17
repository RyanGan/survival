# Simulation code for log-rank power calculations for given enrollment time, event rate, and drop out.

library(survival)
library(dplyr)
library(survival)
library(data.table)

LRPower <- function(n.arm.c = 800, accr.c = 102, fu.c = 12, ann.drop.c = 0.1, n.arm.t = 680, accr.t = 24, fu.t = 24, ann.drop.t = 0.025, HR = 0.76, type.I.error = 0.05, simulation.n=10000){
  ## need to simulation 1000 simulation and do logrank test 1000 times and calculate 1-type II error (if H1 is true but conclude H0)
  ## so power is the proportion of p>0.05
  simulation.p <- sapply(1:simulation.n, function(x){
    tryCatch({

      # compute input parameters from known values above
      haz.rate.c <- -log(0.5)/36  # monthly hazard rate, given median survival time
      drop.rate.c <- -log(1-ann.drop.c)/12            # monthly dropout rate (exponential), given annual percentage dropped
      drop.rate.t <- -log(1-ann.drop.t)/12
      total.time.c <- accr.c + fu.c
      
      # compute individual pt times (relevant times for each individual)
      enr.study.c <- runif(n.arm.c, 0, accr.c)           # enrollment time from STUDY START (assume uniform enrollment)
      
      fu.time.c <- total.time.c - enr.study.c           # individual fu time from PT ENROLLMENT TIME 
      drop.time.c <- rexp(n=n.arm.c, rate=drop.rate.c)  # drop time (from PT ENROLLMENT TIME), if observed
      event.time.c <- rexp(n=n.arm.c, rate=haz.rate.c)    # event time (from PT ENROLLMENT TIME), if observed
      
      # create dataframe for later operations
      arm.data.c <- tibble(drop.time.c=drop.time.c, 
                           event.time.c=event.time.c, 
                           fu.time.c=fu.time.c)
      
      arm.c <- arm.data.c %>%
        rowwise() %>%
        mutate(time=min(event.time.c, fu.time.c, drop.time.c)) %>%
        mutate(event=ifelse(event.time.c < fu.time.c & event.time.c < drop.time.c, 1, 0)) %>%
        mutate(group="control") %>%
        select(-c(drop.time.c, fu.time.c, event.time.c))   # use this line only after checking result
      
      
      # treatment arm
      haz.rate.t <- HR * haz.rate.c  # monthly hazard rate, given median survival time
      drop.rate.t <- -log(1-ann.drop.t)/12            # monthly dropout rate (exponential), given annual percentage dropped
      total.time.t <- accr.t + fu.t
      
      # compute individual pt times (relevant times for each individual)
      enr.study.t <- runif(n.arm.t, 0, accr.t)           # enrollment time from STUDY START (assume uniform enrollment)
      
      fu.time.t <- total.time.t - enr.study.t           # individual fu time from PT ENROLLMENT TIME 
      drop.time.t <- rexp(n=n.arm.t, rate=drop.rate.t)  # drop time (from PT ENROLLMENT TIME), if observed
      event.time.t <- rexp(n=n.arm.t, rate=haz.rate.t)    # event time (from PT ENROLLMENT TIME), if observed
      
      # create dataframe for later operations
      arm.data.t <- tibble(drop.time.t=drop.time.t, 
                           event.time.t=event.time.t, 
                           fu.time.t=fu.time.t)
      
      arm.t <- arm.data.t %>%
        rowwise() %>%
        mutate(time=min(event.time.t, fu.time.t, drop.time.t)) %>%
        mutate(event=ifelse(event.time.t < fu.time.t & event.time.t < drop.time.t, 1, 0)) %>%
        mutate(group="treatment") %>%
        select(-c(drop.time.t, fu.time.t, event.time.t))   # use this line only after checking result
      
      surv.data <- rbind(arm.c, arm.t)
      total.events <- sum(surv.data$event)
      surv.data$SurvObj <- with(surv.data, Surv(time, event))
      
      tmp.model <- survdiff(surv.data$SurvObj ~ surv.data$group) #  With rho = 0 (default) this is the log-rank test
      p.value <- 1 - pchisq(tmp.model$chisq, 1)
      lgr.result <- list("p.value" = p.value, "total.events" = total.events)
      return(lgr.result)},
      error=function(e){return(NA)})
  })
  #print(simulation.p)
  out.i <- simulation.p
  #print(length(which(na.omit(out.i[1,])<= 0.05))/simulation.n)
  #print(mean(as.numeric(out.i[2,])))

out.ii <- data.table("power"=length(which(na.omit(out.i[1,])<= 0.05))/simulation.n, "n.events"=mean(as.numeric(out.i[2,])))
  
return(out.ii) ## return power
}

LRPower(n.arm.c = 800, accr.c = 102, fu.c = 12, ann.drop.c = 0.1, n.arm.t = 680, accr.t = 24, fu.t = 24, ann.drop.t = 0.05, HR = 0.76, type.I.error = 0.05, 
        simulation.n=10000)
#power n.events
#1: 0.9252 677.7297

LRPower(n.arm.c = 800, accr.c = 102, fu.c = 12, ann.drop.c = 0.1, n.arm.t = 680, accr.t = 24, fu.t = 24, ann.drop.t = 0.05, HR = 0.78, type.I.error = 0.05, 
        simulation.n=10000)
#power n.events
#1: 0.8638 682.3198

LRPower(n.arm.c = 800, accr.c = 102, fu.c = 12, ann.drop.c = 0.1, n.arm.t = 680, accr.t = 24, fu.t = 24, ann.drop.t = 0.05, HR = 0.8, type.I.error = 0.05, 
        simulation.n=10000)
#power n.events
#1: 0.7929 687.7839


# number of treatment arm months to reach 80% power 

LRPower(n.arm.c = 800, accr.c = 102, fu.c = 3, ann.drop.c = 0.1, n.arm.t = 680, accr.t = 24, fu.t = 9, ann.drop.t = 0.05, HR = 0.76, type.I.error = 0.05, 
        simulation.n=1000)

LRPower(n.arm.c = 800, accr.c = 102, fu.c = 10, ann.drop.c = 0.1, n.arm.t = 680, accr.t = 24, fu.t = 16, ann.drop.t = 0.05, HR = 0.78, type.I.error = 0.05, 
        simulation.n=1000)

LRPower(n.arm.c = 800, accr.c = 102, fu.c = 12, ann.drop.c = 0.1, n.arm.t = 680, accr.t = 24, fu.t = 25, ann.drop.t = 0.05, HR = 0.8, type.I.error = 0.05, 
        simulation.n=10000)

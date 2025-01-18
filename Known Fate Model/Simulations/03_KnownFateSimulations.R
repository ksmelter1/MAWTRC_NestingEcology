#'---
#' title: Nest Success Modeling of Wild Turkeys in the Mid-Atlantic Region
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script uses simulated data to fit a Bayesian known fate model for female wild turkeys in our study
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modelling daily nest survival

################################
## Simulate Data 

#' Simulating hens.nests dataset
set.seed(123)

#' Number of nests
n_nests <- 44  

#' Create random dates for start and end dates of nesting
#' Sample between 1:28 days from start date to generate end date 
start_dates <- sample(seq.Date(as.Date('2024-05-01'), as.Date('2024-08-01'), by = 'day'), n_nests, replace = TRUE)
end_dates <- start_dates + sample(1:28, n_nests, replace = TRUE)

#' Simulate nest fate (1 for success, 0 for failure)
#' I want a 30% daily nest survival probability for phi
#nest_fate <- rbinom(n_nests, 1, 0.3)

#' Create the hens.nests dataframe
#hens.nests <- data.frame(
#  nestid = 1:n_nests,
#  start.date = start_dates,
#  end.date = end_dates,
#  nest_fate = nest_fate
#)

#' Create the hens.nests dataframe
hens.nests <- data.frame(
  nestid = 1:n_nests,
  start.date = start_dates
)

#' Check the structure of the dataset
glimpse(hens.nests)
summary(hens.nests)

#' Create a sequence of all unique dates between the min start and max end date
inc.dates <- sort(unique(c(hens.nests$start.date, hens.nests$end.date)))
inc.dates <- seq(inc.dates[1], inc.dates[length(inc.dates)], by = 1)

#' Calculate the number of encounter days
occs <- length(inc.dates)

#' Initialize variables
first <- last <- array(NA, dim = nrow(hens.nests))
surv.caps <- matrix(data = NA, nrow = nrow(hens.nests), ncol = occs)

#' Loop through each row of the dataset
for(i in 1:nrow(hens.nests)){ 
  #' Find the index of the start and end dates in the incubation dates sequence
  first[i] <- which(inc.dates == hens.nests$start.date[i]) 
  #last[i] <- which(inc.dates == hens.nests$end.date[i]) 

  for (j in first[i]:occs){  
  #' Fill in the survival caps (incubation days marked as 1)
   surv.caps[i, j]<-rbinom(1,1,0.97)
   if (surv.caps[i, j]==0) {
     break
   }
  } 
  
  last[i]<-ifelse(surv.caps[i,occs] %in% c(1),occs,which(surv.caps[i, ] == 0))

  #' If the nest failed, mark the last incubation day as 0 (failure)
  #if(hens.nests$nest_fate[i] == 0) {
  #  surv.caps[i, last[i]] <- 0
  #} 
}

###############################################
## Intercept Only Known Fate Model in JAGS

#' Heather Gaya Tutorial Known Fate Code
#' Switched to the cloglog link
modelstring.ns_null = "
model {
cloglog(phi) <- beta0 
for (i in 1:n.ind){
  for(t in (first[i]+1):last[i]){
    mu[i,t] <- phi*y[i,t-1]
    y[i,t] ~dbern(mu[i,t])
  }
}
beta0 ~ dunif(-6,6)    # Prior for Intercept
}
"

#' Information to provide JAGS
#' jd = lists of individual capture histories, first observed incubating, last observed incubating
#' ji = inits (0.5)
#' jp = model parameters
jd <- list(n.ind= nrow(hens.nests), y = surv.caps, 
           first = first, 
           last = last)

#' Initial values
ji <- function(){list(beta0 = .1)}

#' Initialize JAGS
jags_model <- jags.model(textConnection(modelstring.ns_null), 
                         data = jd,
                         n.chains = 1)

#' Sample from posterior distribution using parallel
jags_samples <- coda.samples(jags_model, c("beta0", "phi"),
                             n.iter = 40000, 
                             n.burnin = 10000,
                             n.thin = 1,
                             do.parallel = T)

MCMCtrace(jags_samples, pdf = F)

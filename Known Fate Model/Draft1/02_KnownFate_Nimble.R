
#'---
#' title: Nest Success Modeling of Wild Turkeys in the Mid-Atlantic Region
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script uses simulated data to fit a Bayesian known fate model for female wild turkeys in our study
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modeling daily nest survival



#####################
## Load Packages

library(nimble)
library(MCMCvis)


################################
## Nimble Model

nestsuccess <- nimbleCode({
  
  for (i in 1:n.ind) {
    for(t in first[i]:last[i]) {
      y[i, t] ~ dbern(mu[i, t])  
      cloglog(mu[i, t]) <- inprod(beta[1:J], X[i, 1:J])  
    }
  }
  
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
  phi ~ dbeta(1, 1)
  
})

#' Constants
Consts <- list(
  n.ind = nrow(nests.ready),   # Number of rows of nests
  X = X,                       # Covariate matrix (temperature)
  J = ncol(X),                 # Columns of X
  first = first,               # First incubation day for each nest
  last = last                  # Last incubation day for each nest
)

#' Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  y = surv.caps               # Survival data
)

#' Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J),    # Initial values for beta (length should match the number of covariates in X)
  phi = 0.1                   # Initial value for daily survival probability
)

start <- Sys.time()

#' Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = nestsuccess,          # Model code
  constants = Consts,          # Constants defined earlier
  inits = Inits,               # Initial values for the parameters
  data = Data,                 # Data (survival data)
  nburnin = 10000,             # Burn-in period
  niter = 40000,               # Number of iterations
  thin = 1                     # Thinning factor
)

#' View traceplots
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

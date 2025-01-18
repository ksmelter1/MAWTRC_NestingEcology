#'---
#' title: Nest Success Modeling of Wild Turkeys in the Mid-Atlantic Region
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#'  
#' **Purpose**: This script uses simulated data to fit a Bayesian known fate model for female wild turkeys in our study
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modeling daily nest survival

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

}

# Simulate a simple covariate (e.g., temperature) for each nest and day
set.seed(123)

# Assume temperature varies between 10°C and 30°C
temperature <- runif(occs, 10, 30)  # Simulated daily temperature

# Create a design matrix (e.g., intercept and temperature) for each nest
X <- matrix(temperature, nrow = n_nests, ncol = occs)


################################
## Nimble Model

nestsuccess <- nimbleCode({
  # Loop over all nests
  for (i in 1:n.ind) {
    # Loop over the incubation days for each nest (from first[i] to last[i])
    for(t in first[i]:last[i]) {

      # Assuming beta[1] is intercept, beta[2] is temperature effect
      cloglog(mu[i, t]) <- inprod(beta[1:J], X[i, 1:J])  
      
      # Survival for the current day
      y[i, t] ~ dbern(mu[i, t])  # Bernoulli outcome based on survival probability
    }
  }
  
  #' Priors for beta 
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
  # Prior for phi (daily survival probability)
  phi ~ dbeta(1, 1)
})

# Constants
Consts <- list(
  n.ind = nrow(hens.nests),    # Number of nests
  X = X,                       # Covariate matrix (temperature)
  first = first,               # First incubation day for each nest
  last = last                  # Last incubation day for each nest
)

# Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  y = surv.caps                 # Simulated survival data (daily survival/capture data)
)

# Initial values for parameters
Inits <- list(
  beta = rep(0, 2),            # Initial values for the betas (intercept and temperature effect)
  phi = 0.5                    # Initial value for daily survival probability
)

# Run MCMC sampling with Nimble
start <- Sys.time()

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

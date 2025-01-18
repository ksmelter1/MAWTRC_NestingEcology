
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
#' **Purpose**: This script fits a Bayesian known fate model for female wild turkeys in our study
#' 

#####################
## Load Packages ##
#####################

#' Vector of package names
packages <- c("R2jags",
              "MCMCvis",
              "mcmcplots",
              "tidyverse",
              "jagsUI",
              "stringr")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)

################
## Data Prep ##
################

#' Following Heather Gaya's tutorial 
#' https://github.com/heathergaya/JAGS-NIMBLE-Tutorials/blob/master/Known_Fate/Known_Fate_Models.Rmd
#' Load in subset of nesting data
#' Check structure 
hens.nests <- read.csv("Data Management/Csvs/NestAttempts_nodigits_KyleTest.csv") 
glimpse(hens.nests)

#' Rename index column to BirdID 
#' Check structure
colnames(hens.nests)[1] <- "BirdID"
glimpse(hens.nests)

#' Convert start and end dates from characters to date objects
#' Check date format
hens.nests$start.date <- as.Date(hens.nests$start.date, format = "%Y-%m-%d")
hens.nests$end.date <- as.Date(hens.nests$end.date, format = "%Y-%m-%d")
glimpse(hens.nests)

#' Simulate nest fate column in hens.nests df
#' In the process of linking the observations by nestid but not needed for this exercise
#' 1 = At least one egg hatched from a nest (Success)
#' 0 = Nest failed 
#' 44 nests
set.seed(123)
nest_fate <- rbinom(46, 1, 0.3)
hens.nests$nest_fate <- nest_fate
glimpse(hens.nests)

#' Remove rows with NA values in the start and end date columns 
#' Summary check to see if NA values remain
#' Check structure of new df
hens.nests.ready <- tidyr::drop_na(hens.nests)
summary(hens.nests.ready)
glimpse(hens.nests.ready)

#' Create incubation dates object
#' Occurrences object is just the length of the incubation dates
#' We need a list of every day between our start and end period
#' 147 encounter days is the output for occs
inc.dates <- sort(unique(c(hens.nests.ready$start.date, hens.nests.ready$end.date)))
inc.dates <- seq(inc.dates[1],inc.dates[length(inc.dates)], by = 1)
occs <- length(inc.dates)

#' Create Encounter History
#' Loops through each individual and extracts information
#' First = The first day the bird began incubating (Start.Date)
#' Last = The last day before termination of the nesting attempt 
#' Surv.Caps = The encounter history for each individual 
first <- last <- array(NA, dim = nrow(hens.nests.ready))
surv.caps <-  matrix(data = NA, nrow = nrow(hens.nests.ready), ncol = occs) 
for(i in 1:nrow(hens.nests.ready)){ 
  first[i] <- which(inc.dates == hens.nests.ready$start.date[i]) 
  last[i] <- which(inc.dates == hens.nests.ready$end.date[i]) 
  surv.caps[i,first[i]:last[i]] <- 1 
  if(hens.nests.ready$nest_fate[i] == 0)
  {surv.caps[i,last[i]] <- 0} 
}

#' Check work on first individual
hens.nests.ready[1,]

#' NA values for survival caps
#' Need to address this 
surv.caps[1,1:10]

####################################
## Bring in Habitat Predictors

#' Load in habitat covs
load("Data Management/RData/Nest-Site Selection/Covs/02_Covs.RData")

#' Read in captures csv
captures <- read_csv("Data Management/Csvs/raw data/captures.csv") 
captures

#' Filter data to include only hens 
captures <- captures %>%
  dplyr::filter(sex == "F") %>%
  dplyr::select(bandid, sex, age, weight, captyr)
captures

#' Create nest.data object for models
#' Scale continuous predictors
nest.data <- pa.nests.covs 

#' Merge data
nest.data <- merge(captures, nest.data, by = "bandid")

#' Add in columns
nest.data$basal <- pa.nests$stemcnt
nest.data$grssfrb <- pa.nests$grssfrb

#' (?<=_): Only match is there is an underscore immediately before the number we are trying to extract
#' (\\d{4}): Match exactly 4 digits
#' (?=_) : Only match if the four digits are followed 
nest.data <- nest.data %>%
  dplyr::mutate(nestyr = str_extract(nest_id, "(?<=_)(\\d{4})(?=_)"))
glimpse(nest.data)

nest.data$nestyr <- as.numeric(nest.data$nestyr)
nest.data$captyr <- as.numeric(nest.data$captyr)

#' Create a years since capture column
nest.data <- nest.data %>%
  dplyr::mutate(yrsincecap = nestyr-captyr)
glimpse(nest.data)

#' Select columns of interest
nest.data <- nest.data %>%
  dplyr::select(nest_id, bandid, avrgmxv, grssfrb, percwdy,
                case, Developed, Deciduous, Mixed, Evergreen, Agriculture, elev,
                primary, secondary, basal, weight, age, captyr, nestyr, yrsincecap)

#' Scale predictors
nest.data <- nest.data %>%
  dplyr::mutate(basal = scale(basal)) %>%
  dplyr::mutate(percwdy = scale(percwdy)) %>%
  dplyr::mutate(grssfrb = scale(grssfrb)) %>%
  dplyr::mutate(avgvisob = scale(avrgmxv)) %>%
  dplyr::mutate(elev = scale(elev)) %>%
  dplyr::mutate(primary = scale(primary)) %>%
  dplyr::mutate(secondary = scale(secondary)) %>%
  dplyr::mutate(weight = scale(weight)) 

str(nest.data)
glimpse(nest.data)

#' Change structure from a matrix to numeric (Scaled covariates)
nest.data <- nest.data %>%
  dplyr::mutate(percwdy = as.numeric(percwdy)) %>%
  dplyr::mutate(grssfrb = as.numeric(grssfrb)) %>%
  dplyr::mutate(avgvisob = as.numeric(avgvisob)) %>%
  dplyr::mutate(elev = as.numeric(elev)) %>%
  dplyr::mutate(primary = as.numeric(primary)) %>%
  dplyr::mutate(secondary = as.numeric(secondary)) %>%
  dplyr::mutate(basal = as.numeric(basal)) %>%
  dplyr::mutate(weight = as.numeric(weight)) %>%
  dplyr::mutate(Deciduous = as.numeric(Deciduous)) %>%
  dplyr::mutate(Mixed = as.numeric(Mixed)) %>%
  dplyr::mutate(Evergreen= as.numeric(Evergreen)) %>%
  dplyr::mutate(Agriculture = as.numeric(Agriculture)) %>%
  dplyr::mutate(Developed = as.numeric(Developed))

glimpse(nest.data)
str(nest.data)

#' Change case to numeric
nest.data$case <- as.numeric(nest.data$case)

#' Order df by nestid
nest.data <- nest.data[order(nest.data$nest_id),]

#' Change Nest_ID_V to numeric 
#' First must remove underlines
nest.data$nest_id <- gsub("_", "", nest.data$nest_id)
nest.data$nest_id<- as.numeric(nest.data$nest_id)

#' Create str_id column 
#' This allows the loop to iterate through the steps associated with each bird 
#' cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  group_by(nest_id) %>%
  mutate(str_ID=cur_group_id())

#' Check
str(nest.data)
glimpse(nest.data)
summary(nest.data)

#' Function to replace all NAs with 0 in a dataframe
#' Since we scaled the data, the mean is zero
fill_na_with_zero <- function(df) {
  df[is.na(df)] <- 0
  return(df)
}

#' Apply the function
nest.data.ready <- fill_na_with_zero(nest.data)
summary(nest.data.ready)

#' Assign juvenile as the reference level
nest.data.ready$age <- ifelse(nest.data.ready$age == "A", 1, 
                              ifelse(nest.data.ready$age == "J", 0, NA))


#' Dealing with scaling age ad hoc
#' If the bird is an adult and the years since capture is >1 assign it as an adult
#' If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nest.data.ready$age <- ifelse(nest.data.ready$age == "A" & nest.data.ready$yrsincecap >= 1, "A", nest.data.ready$age)
glimpse(nest.data.ready)

#' Check
glimpse(nest.data.ready)
summary(nest.data.ready)
str(nest.data.ready$age)

###############################################
## Intercept Only Known Fate Model in JAGS

#' Heather Gaya Tutorial Known Fate Code
#' Intercept only nest success model
#' Modeling Daily Nest Survival Probability
modelstring.ns_null = "
model {
logit(phi) <- beta0 
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
jd <- list(n.ind= nrow(hens.nests.ready), y = surv.caps, 
                                          first = first, 
                                          last = last)
ji <- function(){list(beta0 = .5)}

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

#############################################
## Known Fate Model with Habitat Covariates

#' Fine-scale habitat predictors
#' Coarse-scale habitat predictors

modelstring.ns_habitat = "
model {
  # Logistic regression for the logit of phi
  logit(phi) <- beta0 + inprod(X[i,], betas)  # Linear combination of covariates
  for (i in 1:n.ind) {
    for (t in (first[i] + 1):last[i]) {
      mu[i,t] <- phi * y[i,t-1]  # Logistic model for each individual
      y[i,t] ~ dbern(mu[i,t])    # Bernoulli distribution for survival
    }
  }
  
  # Priors
  beta0 ~ dunif(-6, 6)  # Prior for intercept
  for (j in 1:ncol(X)) {
    betas[j] ~ dnorm(0, 0.1)  # Priors for coefficients of covariates
  }
}
"
#' Model parameters
X <- cbind(
  rep(1, nrow(nest.data.ready)),  # Intercept (1)
  nest.data.ready$primary,         # Distance to Primary Road
  nest.data.ready$secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Deciduous,       # Deciduous Forest
  nest.data.ready$Agriculture,     # Agriculture
  nest.data.ready$elev,            # Elevation
  nest.data.ready$avgvisob,        # Visual Obstruction
  nest.data.ready$grssfrb,         # Percent Grass/Forb
  nest.data.ready$percwdy,         # Percent Woody
  nest.data.ready$basal           # Basal Area
)

#' Information to provide JAGS
jd <- list(
  n.ind = nrow(hens.nests.ready), 
  y = surv.caps, 
  first = first, 
  last = last,
  X = X  # Including the covariates matrix X
)

#' Initial values for the model
ji <- function() {
  list(
    beta0 = 0.5, 
    betas = rep(0, ncol(X))  # Initialize coefficients for covariates
  )
}

#' Parameters to monitor
jp <- c("beta0", "betas", "phi")

#' Fit the model using JAGS
ns_habitat <- jags.samples(
  model = modelstring.ns_null, 
  data = jd, 
  inits = ji, 
  n.chains = 1, 
  burnin = 1000, 
  n.iter = 2000
)

#############################################
## Full Known Fate Model without RE

#' Fine-scale habitat predictors
#' Coarse-scale habitat predictors
#' Weather predictors
#' Individual Covariates 

modelstring.ns_habitat = "
model {
  # Logistic regression for the logit of phi
  logit(phi) <- beta0 + inprod(X[i,], betas)  # Linear combination of covariates
  for (i in 1:n.ind) {
    for (t in (first[i] + 1):last[i]) {
      mu[i,t] <- phi * y[i,t-1]  # Logistic model for each individual
      y[i,t] ~ dbern(mu[i,t])    # Bernoulli distribution for survival
    }
  }
  
  # Priors
  beta0 ~ dunif(-6, 6)  # Prior for intercept
  for (j in 1:ncol(X)) {
    betas[j] ~ dnorm(0, 0.1)  # Priors for coefficients of covariates
  }
}
"
#' Model parameters (Need to add to this)
X <- cbind(
  rep(1, nrow(nest.data.ready)),  # Intercept (1)
  nest.data.ready$primary,         # Distance to Primary Road
  nest.data.ready$secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Deciduous,       # Deciduous Forest
  nest.data.ready$Agriculture,     # Agriculture
  nest.data.ready$elev,            # Elevation
  nest.data.ready$avgvisob,        # Visual Obstruction
  nest.data.ready$grssfrb,         # Percent Grass/Forb
  nest.data.ready$percwdy,         # Percent Woody
  nest.data.ready$basal           # Basal Area
)

#' Information to provide JAGS
jd <- list(
  n.ind = nrow(hens.nests.ready), 
  y = surv.caps, 
  first = first, 
  last = last,
  X = X  # Including the covariates matrix X
)

#' Initial values for the model
ji <- function() {
  list(
    beta0 = 0.5, 
    betas = rep(0, ncol(X))  # Initialize coefficients for covariates
  )
}

#' Parameters to monitor
jp <- c("beta0", "betas", "phi")

#' Fit the model using JAGS
ns_habitat <- jags.samples(
  model = modelstring.ns_null, 
  data = jd, 
  inits = ji, 
  n.chains = 1, 
  burnin = 1000, 
  n.iter = 2000
)

#############################################
## Full Known Fate Model with RE

#' Fine-scale habitat predictors
#' Coarse-scale habitat predictors
#' Weather predictors
#' Individual Covariates 
#' Random Effect Accounting for Individual Variation 

modelstring.ns_habitat = "
model {
  # Logistic regression for the logit of phi
  logit(phi) <- beta0 + inprod(X[i,], betas) + alpha[BirdID[i]]  # Linear combination of covariates
  for (i in 1:n.ind) {
    for (t in (first[i] + 1):last[i]) {
      mu[i,t] <- phi * y[i,t-1]  # Logistic model for each individual
      y[i,t] ~ dbern(mu[i,t])    # Bernoulli distribution for survival
    }
  }
  
  # Priors
  beta0 ~ dunif(-6, 6)  # Prior for intercept
  for (j in 1:ncol(X)) {
    betas[j] ~ dnorm(0, 0.1)  # Priors for coefficients of covariates
  }
  
   for(k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
}
"
#' Model parameters (Need to add to this)
X <- cbind(
  rep(1, nrow(nest.data.ready)),  # Intercept (1)
  nest.data.ready$primary,         # Distance to Primary Road
  nest.data.ready$secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Deciduous,       # Deciduous Forest
  nest.data.ready$Agriculture,     # Agriculture
  nest.data.ready$elev,            # Elevation
  nest.data.ready$avgvisob,        # Visual Obstruction
  nest.data.ready$grssfrb,         # Percent Grass/Forb
  nest.data.ready$percwdy,         # Percent Woody
  nest.data.ready$basal           # Basal Area
)

#' Information to provide JAGS
jd <- list(
  n.ind = nrow(hens.nests.ready), 
  y = surv.caps, 
  first = first, 
  last = last,
  X = X  # Including the covariates matrix X
)

#' Initial values for the model
ji <- function() {
  list(
    beta0 = 0.5, 
    betas = rep(0, ncol(X))  # Initialize coefficients for covariates
  )
}

#' Parameters to monitor
jp <- c("beta0", "betas", "phi")

#' Fit the model using JAGS
ns_habitat <- jags.samples(
  model = modelstring.ns_null, 
  data = jd, 
  inits = ji, 
  n.chains = 1, 
  burnin = 1000, 
  n.iter = 2000
)
################################################################################
################################################################################

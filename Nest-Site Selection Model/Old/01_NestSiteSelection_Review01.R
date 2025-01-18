#'---
#' title: Nest-site selection of female wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script creates a Bayesian conditional logistic regression model for nest-site selection in JAGs using the gathered covariates 

#####################
## Load Packages ##
####################

#' Vector of package names
packages <- c("R2jags",
              "MCMCvis")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


#' Read in RDS file with all data needed to fit the model
nest.data.ready <- readRDS("nest.data.ready.RDS")


############################
## JAGS Models

# We do not include developed in the model because it is correlated with agriculture


nest.selection.mod <- "  
 model{
  for (i in 1:I){
   case[i]~dpois(lambda[i])
   log(lambda[i])<-primary[i]*beta.primary +
   secondary[i]*beta.secondary +
   deciduous[i]*beta.deciduous_forest + 
   mixed[i]*beta.mixed_forest + 
   evergreen[i]*beta.evergreen_forest + 
   agriculture[i]*beta.agriculture + 
   elevation[i] *beta.elevation + 
   visob[i]* beta.visob +
   grassforb[i]*beta.grassforb +
   percwoody[i]*beta.percwoody + 
   basalarea[i]*beta.basalarea +
   age[i] * beta.age +
   alpha[str_ID[i]] 
  }
  
  #Priors
  beta.deciduous_forest~dnorm(0,0.0001)  #Deciduous Forest
  beta.mixed_forest~dnorm(0,0.0001)      #Mixed Forest
  beta.evergreen_forest~dnorm(0,0.0001)  #Mixed Forest
  beta.agriculture~dnorm(0,0.0001)       #Agriculture
  beta.primary~dnorm(0,0.0001)           #Primary
  beta.secondary~dnorm(0,0.0001)         #Secondary
  beta.elevation~dnorm(0,0.0001)         #Elevation
  beta.visob~dnorm(0,0.0001)             #Average Visual Obstruction
  beta.grassforb~dnorm(0,0.0001)         #Percent Grass Forb
  beta.percwoody~dnorm(0,0.0001)         #Percent Woody
  beta.age~dnorm(0,0.0001)             #Age Adult (Juvenile is the reference level)
  beta.basalarea~dnorm(0,0.0001)         #Basal Area
  
  for (k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
 }
"
#' Data list for JAGS
jags_data <- list(
  case = nest.data.ready$case,
  primary = nest.data.ready$primary,
  secondary = nest.data.ready$secondary,
  mixed = nest.data.ready$Mixed,
  evergreen = nest.data.ready$Evergreen,
  deciduous = nest.data.ready$Deciduous,
  agriculture = nest.data.ready$Agriculture,
  elevation = nest.data.ready$elev,
  visob = nest.data.ready$avgvisob,
  grassforb = nest.data.ready$grssfrb,
  percwoody = nest.data.ready$percwdy,
  age = nest.data.ready$age,
  basalarea = nest.data.ready$basal,
  str_ID = c(nest.data.ready$str_ID),
  I = nrow(nest.data.ready),
  K =length(unique(nest.data.ready$str_ID))
)

#' Initialize JAGS
jags_model <- jags.model(textConnection(nest.selection.subset.mod), 
                         data = jags_data,
                         n.chains = 3)

#' Burn-in and sampling from posterior distribution
jags_samples.nss.com1 <- coda.samples(jags_model, c("beta.deciduous_forest",
                                                    "beta.mixed_forest",
                                                    "beta.evergreen_forest",
                                                    "beta.agriculture",
                                                    "beta.primary",
                                                    "beta.secondary",
                                                    "beta.elevation",
                                                    "beta.visob",
                                                    "beta.grassforb",
                                                    "beta.percwoody",
                                                    "beta.age",
                                                    "beta.basalarea"), 
                                      n.iter = 200000, 
                                      n.burnin = 30000,
                                      n.thin = 5)
#' View traceplots
MCMCtrace(jags_samples.nss.com1, pdf=F)

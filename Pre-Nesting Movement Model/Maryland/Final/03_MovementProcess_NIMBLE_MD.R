#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis) in NIMBLE
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script fits SSF models in frequentist and Bayesian frameworks and extracts coefficient data
#' **Last Updated**: 1/25/25


################################################################################
## Load Packages 

#' Vector of package names
packages <- c("tidyverse",
              "MCMCvis",
              "matrixStats",
              "nimble")


#' Function to load a package or install it if not already installed
#' Found this function 
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


load("Data Management/RData/Pre-Nesting Movement Model/Maryland/Covariates/20250226_MD_Ready.RData") 
rm(x)


dat_2.ready <- final_filtered_random_steps %>%
  dplyr::mutate(secondary = scale(secondary)) %>%
  dplyr::mutate(primary = scale(primary)) 

glimpse(dat_2.ready)

dat_2.ready <- dat_2.ready %>% 
  dplyr::mutate(secondary = as.numeric(secondary)) %>%
  dplyr::mutate(primary = as.numeric(primary)) 

glimpse(dat_2.ready)


################################################################################
## Organize Data 

#' Check
str(dat_2.ready)

#' Change ID to numeric
dat_2.ready$id <- gsub("_", "", dat_2.ready$id)
dat_2.ready$id<- as.numeric(dat_2.ready$id)

#' Create strata column by merging BirdID and step_id
dat_2.ready$NA_ID <- paste(dat_2.ready$id,
                           dat_2.ready$step_id_,
                           sep = "_")

#' Change NA_ID to integer
dat_2.ready$NA_ID <- gsub("_", "", dat_2.ready$NA_ID)
dat_2.ready$NA_ID<- as.numeric(dat_2.ready$NA_ID)

class(dat_2.ready)
str(dat_2.ready)

#' Add numerical variable for animals:
dat_2.ready$ANIMAL_ID <- as.numeric(as.factor(dat_2.ready$id))

#' Stratum ID is given as "NA_ID" in the data; 
#' It is easier to have sequential enumeration, so let's generate a new stratum-ID variable str_ID:
d.map <- data.frame(NA_ID=unique(dat_2.ready$NA_ID),
                    str_ID=1:length(unique(dat_2.ready$NA_ID)))
dat_2.ready$str_ID <- d.map[match(dat_2.ready$NA_ID,d.map$NA_ID),"str_ID"]
dat_2.ready <- dat_2.ready[order(dat_2.ready$str_ID),] 
glimpse(dat_2.ready)

#' Function to fill NA values with 0 throughout dataframe
fill_NA_with_value <- function(df, value = 0) {
  #' Ensure that all columns are of numeric type to handle the replacement
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      #' Replace NA values with the specified value
      replace(col, is.na(col), value)
    } else {
      #' Return the column unchanged if it's not numeric
      col
    }
  })
  return(df)
}

#' Apply the function to fill NA values with 0.0001
dat_2.ready <- fill_NA_with_value(dat_2.ready)

#' Check
which(dat_2.ready$step_id_==3)
table(dat_2.ready$case_)
class(dat_2.ready)
summary(dat_2.ready)

#' Model parameters
X <- cbind(
  rep(1, nrow(dat_2.ready)),   # Intercept (1)
  dat_2.ready$primary,         # Distance to Primary Road
  dat_2.ready$secondary,       # Distance to Secondary Road
  dat_2.ready$Mixed,           # Mixed Forest
  dat_2.ready$Evergreen,       # Evergreen Forest
  dat_2.ready$Developed,       # Developed
  dat_2.ready$Agriculture,     # Agriculture
  dat_2.ready$Grassland,       # Grassland    
  dat_2.ready$Wetland)           # Water

rm(list = setdiff(ls(), c("dat_2.ready", "X")))


################################################################################
## Nimble Model

movemodel <- nimbleCode({
  for (i in 1:I){
    use[i] ~ dpois(lambda[i])  
    
    log(lambda[i]) <- inprod(beta[1:J], X[i, 1:J]) + alpha[str_ID[i]]  
  }
  
  #' Priors for beta 
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
  #' Priors for alpha 
  for (k in 1:K){
    alpha[k] ~ dnorm(0, sd = sqrt(1 / 0.000001))  
  }
})

#' Constants 
Consts <- list(
  str_ID = as.numeric(dat_2.ready$str_ID), 
  J = ncol(X),                             
  I = nrow(dat_2.ready),                   
  K = length(unique(dat_2.ready$str_ID))   
)

#' Data
Data <- list(
  X = X,                                   
  use = dat_2.ready$case_                  
)

#' Initial values 
Inits <- list(
  beta = rep(0, ncol(X)),                  
  alpha = rep(0, length(unique(dat_2.ready$str_ID)))  
)

start <- Sys.time()

#' Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = movemodel, 
  constants = Consts, 
  inits = Inits,
  data = Data,
  nburnin = 10000, 
  niter = 40000,
  thin = 1
)

end <- Sys.time()

saveRDS(nimbleMCMC_samples, "20250220_MovementRunMD.RDS")


#' Print the means and standard deviations of the posterior samples for the beta coefficients
colMeans(nimbleMCMC_samples[, 34077:34085])
colSds(nimbleMCMC_samples[, 34077:34085])
#' View traceplots
#MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 34077:34085]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", "Distance to Primary Road", "Distance to Secondary Road", 
               "Mixed Forest", "Evergreen Forest", "Developed", "Agriculture",
               "Grassland/Shrub", "Wetland")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)

saveRDS(samples_df, "20250226_samples_df_MD.RDS")

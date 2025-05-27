#'---
#' title: Nest-site selection of wild turkeys in Maryland (an SSF analysis)
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script fits a nest-site selection model for Maryland
#' **Last Updated**: 5/12/2025

################################################################################
## Load Packages 

packages <- c("matrixStats",
              "nimble",
              "MCMCvis",
              "tidyverse",
              "sf")

load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

lapply(packages, load_packages)


################################################################################
## Data Preparation 

#' Load in RData
load("Data Management/RData/Nest-Site Selection/Covs/Multi-State/Manuscript/02_Covs_AllStates_Ready.RData")

#' Obtain data from Maryland and drop geometry columns
#' Had some NA nest IDs that I removed
nest.data <- md.nests.covs %>%
  st_drop_geometry() %>%
  drop_na()
str(md.nests.covs)

#' Subset dataframe and rename columns 
nest.data <- nest.data %>%
  dplyr::select(NestID, BandID,
                Case, Developed, Deciduous, Mixed, Evergreen, Pasture, Crop,
                primary, secondary, Grassland, Water, Wetland) %>%
  dplyr::rename("Primary" = primary) %>%
  dplyr::rename("Secondary" = secondary) 

#' Switch coding to UTF-8
nest.data <- nest.data %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

#' Convert predictors to numeric
nest.data <- nest.data %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary)) %>%
  dplyr::mutate(Developed = as.numeric(Developed)) %>%
  dplyr::mutate(Deciduous = as.numeric(Deciduous)) %>%
  dplyr::mutate(Mixed = as.numeric(Mixed)) %>%
  dplyr::mutate(Evergreen = as.numeric(Evergreen)) %>%
  dplyr::mutate(Pasture = as.numeric(Pasture)) %>%
  dplyr::mutate(Crop = as.numeric(Crop)) %>%
  dplyr::mutate(Grassland = as.numeric(Grassland)) %>%
  dplyr::mutate(Wetland = as.numeric(Wetland)) %>%
  dplyr::mutate(Water = as.numeric(Water))
str(nest.data)
glimpse(nest.data)

#' Scale continous predictors
nest.data <- nest.data %>%
  dplyr::mutate(Primary = scale(Primary)) %>%
  dplyr::mutate(Secondary = scale(Secondary)) 
glimpse(nest.data)
str(nest.data)

#' Convert scaled predictors back to numeric
nest.data <- nest.data %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary))

#' Convert case to numeric
nest.data <- nest.data %>%
  st_drop_geometry() %>%
  dplyr::mutate(Case = as.numeric(Case))

#' Order nest.data by NestID
nest.data <- nest.data[order(nest.data$NestID),]

#' Change Nest_ID_V to numeric 
#' First must remove underlines
nest.data$NestID <- gsub("_", "", nest.data$NestID)
nest.data$NestID <- as.numeric(nest.data$NestID)

#' Create str_id column 
#' This allows the loop to iterate through the steps associated with each bird 
#' cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  group_by(NestID) %>%
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

#' Check
summary(nest.data.ready)
glimpse(nest.data.ready)
cor(nest.data.ready$Primary, nest.data.ready$Secondary)


################################################################################
## Nimble Model

#' Conditional logistic regression
#' Prior for stratum-specific intercept variance set to 10^6
#' Prior for betas set to 0,1

nestmodel<-nimbleCode({
  for (i in 1:I){
    use[i]~dpois(lambda[i])
    
    log(lambda[i])<-inprod(beta[1:J],X[i,1:J])+alpha[str_ID[i]]
  }
  
  #' Priors
  for(j in 1:J){
    beta[j]~dnorm(0,sd= 1)
  }
  
  for(k in 1:K){
    alpha[k]~dnorm(0,sd=sqrt(1/0.000001))
  }
}
)

#' Model parameters
X <- cbind(
  rep(1, nrow(nest.data.ready)),   # Intercept (1)
  nest.data.ready$Primary,         # Distance to Primary Road
  nest.data.ready$Secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Developed,       # Developed
  nest.data.ready$Pasture,         # Pasture
  nest.data.ready$Crop,            # Crop
  nest.data.ready$Grassland,       # Grassland/Shrub
  nest.data.ready$Wetland          # Wetland
)

Consts <- list(str_ID = as.numeric(nest.data.ready$str_ID),
               J=ncol(X),
               I = nrow(nest.data.ready),
               K =length(unique(nest.data.ready$str_ID)))

#' Data list for nimble
Data<-list(X=X,use = nest.data.ready$Case)

Inits<-list(beta=rep(0,ncol(X)),alpha=rep(0,length(unique(nest.data.ready$str_ID))))

nimbleMCMC_samples <- nimbleMCMC(code = nestmodel, 
                                 constants = Consts, 
                                 inits=Inits,
                                 data = Data,
                                 nburnin = 10000,
                                 niter = 20000,
                                 thin = 3,
                                 nchains = 1)
#' View traceplots
MCMCtrace(nimbleMCMC_samples [,101:110], pdf = F)
colMeans(nimbleMCMC_samples[, 101:110])
colSds(nimbleMCMC_samples[,  101:110])

#' Extract the posterior samples 
beta_samples <- nimbleMCMC_samples[, 101:110]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", "Distance to Primary Road", "Distance to Secondary Road",
               "Mixed Forest", "Evergreen Forest", "Developed", "Pasture", "Crop",
               "Grassland/Shrub", "Wetland")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names


################################################################################
## Check Proportion of Developed Samples less than 0

#' Count how many values in 'Developed' are less than 0
count_less_than_0 <- sum(samples_df$Developed < 0, na.rm = TRUE)

#' Total number of non-NA values in Developed
total_samples <- sum(!is.na(samples_df$Developed))

#' Calculate the proportion
prop_less_than_0 <- count_less_than_0 / total_samples

#' Print the proportion
prop_less_than_0


################################################################################
## Save RData file *InsertDate_NimbleResults

################################################################################
###############################################################################X


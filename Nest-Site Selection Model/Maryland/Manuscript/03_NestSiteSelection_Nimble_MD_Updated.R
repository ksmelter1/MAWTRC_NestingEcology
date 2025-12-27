#'---
#' title: Nest-site selection of female wild turkeys in Maryland 
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script fits a nest-site selection model for Maryland
#' **Last Updated**: 12/27/2025

################################################################################
## Load Packages 

packages <- c("matrixStats",
              "nimble",
              "MCMCvis",
              "tidyverse",
              "sf",
              "coda")

load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

lapply(packages, load_packages)


################################################################################
## Data Preparation 

# Load in RData
load("Data Management/RData/Nest-Site Selection/Covs/Multi-State/Manuscript/02_Covs_AllStates_Ready_Updated.RData")

# Obtain data from Maryland and drop geometry columns
# Had some NA nest IDs that I removed
# 1214_2024_2 had a lat/long error
nest.data <- md.nests.covs %>%
  st_drop_geometry() %>%
  drop_na() %>%
  dplyr::filter(NestID != "1214_2024_2")
str(md.nests.covs)

# MD Sample
md.sample <- read_csv("Samples/Maryland/NestingSample_MD.updated.csv")

# Convert to numeric to allow for data join
nest.data$BandID <- as.numeric(nest.data$BandID)

# Only keep observations of nest.data that exist in our sample
# Need to exclude the Atlantic ocean nests
nest.data <- right_join(nest.data, md.sample)

# Subset dataframe and rename columns 
# One potential nest is in open water and cannot be used for estimation
nest.data <- nest.data %>%
  dplyr::select(NestID, BandID,
                Case, Developed, Deciduous, Mixed, Evergreen, Pasture, Crop,
                primary, secondary, Grassland, Water, Wetland) %>%
  dplyr::rename("Primary" = primary) %>%
  dplyr::rename("Secondary" = secondary) 

# Switch coding to UTF-8
nest.data <- nest.data %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

# Convert predictors to numeric
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

# Scale continous predictors
nest.data <- nest.data %>%
  dplyr::mutate(Primary = scale(Primary)) %>%
  dplyr::mutate(Secondary = scale(Secondary)) 
glimpse(nest.data)
str(nest.data)

# Convert scaled predictors back to numeric
nest.data <- nest.data %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary))

# Convert case to numeric
nest.data <- nest.data %>%
  st_drop_geometry() %>%
  dplyr::mutate(Case = as.numeric(Case))

# Order nest.data by NestID
nest.data <- nest.data[order(nest.data$NestID),]

# Change Nest_ID_V to numeric 
# First must remove underlines
nest.data$NestID <- gsub("_", "", nest.data$NestID)
nest.data$NestID <- as.numeric(nest.data$NestID)

# Create str_id column 
# This allows the loop to iterate through the steps associated with each bird 
# cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  group_by(NestID) %>%
  mutate(str_ID=cur_group_id())

# Check
str(nest.data)
glimpse(nest.data)
summary(nest.data)

# Create nest.data.ready object
# Couldn't estimate a strata for one nest
nest.data.ready <- nest.data %>%
  drop_na()

# Check
summary(nest.data.ready)
glimpse(nest.data.ready)
cor(nest.data.ready$Primary, nest.data.ready$Secondary)


################################################################################
## Nimble Nest-Site Selection Model

# Conditional logistic regression
# Prior for stratum-specific intercept variance set to 10^6 following Muff et al. 2020
# Prior for betas set to 0,1
# Year random effect will be added to account for differences between years in 2023, 2024, and 2025

nestmodel<-nimbleCode({
  for (i in 1:I){
    use[i]~dpois(lambda[i])
    
    log(lambda[i])<-inprod(beta[1:J],X[i,1:J])+alpha[str_ID[i]]
  }
  
  # Priors
  for(j in 1:J){
    beta[j]~dnorm(0,sd= 1)
  }
  
  for(k in 1:K){
    alpha[k]~dnorm(0,sd=sqrt(1/0.000001))
  }
}
)

# Model parameters
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

# Data list for nimble
Data<-list(X=X,use = nest.data.ready$Case)

Inits<-list(beta=rep(0,ncol(X)),alpha=rep(0,length(unique(nest.data.ready$str_ID))))

nimbleMCMC_samples <- nimbleMCMC(code = nestmodel, 
                                 constants = Consts, 
                                 inits=Inits,
                                 data = Data,
                                 nburnin = 10000,
                                 niter = 30000,
                                 thin = 3,
                                 nchains = 1)
# View traceplots
MCMCtrace(nimbleMCMC_samples [,97:106], pdf = F)
colMeans(nimbleMCMC_samples[, 97:106])
colSds(nimbleMCMC_samples[,  97:106])

# Extract the posterior samples 
beta_samples <- nimbleMCMC_samples[, 97:106]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

# Create a vector of new names
new_names <- c("Intercept",
               "Distance to Primary Road",
               "Distance to Secondary Road",
               "Mixed Forest",
               "Evergreen Forest", 
               "Developed", 
               "Pasture",
               "Crop",
               "Grassland/Shrub",
               "Wetland")

# Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

################################################################################
## Output Trace Plots

# Force conversion to plain matrix
samples_matrix <- as.matrix(nimbleMCMC_samples[, 97:106])

# Convert to coda mcmc object
mcmc_samples <- mcmc(samples_matrix)

# Assign new names
colnames(mcmc_samples) <- new_names

# Generate trace plots
MCMCtrace(mcmc_samples,
          iter = 30000,
          pdf = T)


################################################################################
## Check Proportion of Developed Samples less than 0

# This was done to see if we could argue there was an effect for developed on selection
# Ended up just using the CrIs to maintain consistency in interpretation of effects

# Count how many values in 'Developed' are less than 0
count_less_than_0 <- sum(samples_df$Developed < 0, na.rm = TRUE)

# Total number of non-NA values in Developed
total_samples <- sum(!is.na(samples_df$Developed))

# Calculate the proportion
prop_less_than_0 <- count_less_than_0 / total_samples

# Print the proportion
prop_less_than_0


################################################################################
## Save RData file *InsertDate_NimbleResults

################################################################################
###############################################################################X


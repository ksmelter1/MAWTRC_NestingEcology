#'---
#' title: Nest-site selection of wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script fits a nest-site selection model for Pennsylvania
#' **Last Updated**: 5/12/2025

################################################################################
## Load Packages 

packages <- c("matrixStats",
              "nimble",
              "MCMCvis",
              "tidyverse",
              "sf",
              "coda",
              "GGally")

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
load("Data Management/RData/Nest-Site Selection/Covs/Multi-State/Manuscript/02_Covs_AllStates_Ready_Updated.RData")

#' This sample matches the Pre-Nesting and known fate model
Sample_PA <- read_csv("Samples/Pennsylvania/PA_Sample_Updated.csv")
Sample_PA

#' Drop geometry column and create unique identifier for NestID
pa.nests.covs <- sf::st_drop_geometry(pa.nests.covs)
pa.nests.covs$NestID1 <- paste(pa.nests.covs$NestID, pa.nests.covs$PlotType, sep = "_")

#' Read in basal area data
basal_summary <- read_csv("Data Management/Csvs/Processed/Covariates/Pennsylvania/BasalArea.csv")

#' Keep only samples that exist in Sample_PA in pa.nests.covs
pa.nests.covs.1 <- right_join(pa.nests.covs, Sample_PA)
pa.nests.covs.2 <- right_join(basal_summary, pa.nests.covs.1) 

#' Zero values represent sites with no trees for basal area
pa.nests.covs.2 <- pa.nests.covs.2 %>%
  replace_na(list(Basal = 0)) %>%
  dplyr::filter(Case != "NA")


#' Rename object
nest.data <- pa.nests.covs.2 
str(nest.data)
summary(nest.data)

#' Subset dataframe and rename columns 
nest.data <- nest.data %>%
  dplyr::select(NestID, BandID, , PercGrassForb, PercWoody,
                Case, Developed, Deciduous, Mixed, Evergreen, Pasture, Crop,
                primary, secondary, StemCount, Grassland, AvgVO, Water, Wetland, 
                PercFern, Basal) %>%
  dplyr::rename("Primary" = primary) %>%
  dplyr::rename("Secondary" = secondary)
nest.data

#' Switch coding to UTF-8
nest.data <- nest.data %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

#' Change variables to be scaled to numeric
nest.data <- nest.data %>%
  dplyr::mutate(PercWoody = as.numeric(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = as.numeric(PercGrassForb)) %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary)) %>%
  dplyr::mutate(PercFern = as.numeric(PercFern)) %>%
  dplyr::mutate(AvgVO = as.numeric(AvgVO)) %>%
  dplyr::mutate(StemCount = as.numeric(StemCount)) %>%
  dplyr::mutate(Basal = as.numeric(Basal)) 
glimpse(nest.data)
str(nest.data)

#' Scale continous predictor variables
nest.data <- nest.data %>%
  dplyr::mutate(StemCount = scale(StemCount)) %>%
  dplyr::mutate(PercWoody = scale( PercWoody)) %>%
  dplyr::mutate(PercGrassForb = scale(PercGrassForb)) %>%
  dplyr::mutate(Primary = scale(Primary)) %>%
  dplyr::mutate(Secondary = scale(Secondary)) %>%
  dplyr::mutate(PercFern = scale(PercFern)) %>%
  dplyr::mutate(AvgVO = scale(AvgVO)) %>%
  dplyr::mutate(Basal = scale(Basal)) 
str(nest.data)
glimpse(nest.data)

#' Convert all variables to numeric including all scaled predictors
nest.data <- nest.data %>%
  dplyr::mutate(PercWoody = as.numeric(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = as.numeric(PercGrassForb)) %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary)) %>%
  dplyr::mutate(PercFern = as.numeric(PercFern)) %>%
  dplyr::mutate(AvgVO = as.numeric(AvgVO)) %>%
  dplyr::mutate(Developed = as.numeric(Developed)) %>%
  dplyr::mutate(Pasture = as.numeric(Pasture)) %>%
  dplyr::mutate(Crop = as.numeric(Crop)) %>%
  dplyr::mutate(Deciduous = as.numeric(Deciduous)) %>%
  dplyr::mutate(Mixed = as.numeric(Mixed)) %>%
  dplyr::mutate(Evergreen = as.numeric(Evergreen)) %>%
  dplyr::mutate(Grassland = as.numeric(Grassland)) %>%
  dplyr::mutate(Wetland = as.numeric(Wetland)) %>%
  dplyr::mutate(Water = as.numeric(Water)) %>%
  dplyr::mutate(StemCount = as.numeric(StemCount)) %>%
  dplyr::mutate(Basal = as.numeric(Basal))
glimpse(nest.data)
str(nest.data)
summary(nest.data)

# Keep only nests that have at least one potential nest option within strata
nest.data <- nest.data %>%
  dplyr::group_by(NestID) %>%
  dplyr::filter(any(Case != 1) & any(Case == 0)) %>%
  dplyr::ungroup()

#' Drop geometry column and change case to numeric
#' 1 the nest was used, 0 it was available to the hen 
nest.data <- nest.data %>%
  st_drop_geometry() %>%
  dplyr::mutate(Case = as.numeric(Case))

#' Order dataframe by NestID
nest.data <- nest.data[order(nest.data$NestID),]

#' Change Nest_ID_V to numeric 
#' First must remove underlines
nest.data$NestID <- gsub("_", "", nest.data$NestID)
nest.data$NestID<- as.numeric(nest.data$NestID)

#' Create str_id column 
#' This allows the loop to iterate through the steps associated with each bird 
#' cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  dplyr::group_by(NestID) %>%
  dplyr::mutate(str_ID=cur_group_id())

#' Check data
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

#' Order by strata ID
nest.data.ready <- nest.data.ready %>%
  dplyr::arrange(str_ID)

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
    
    log(lambda[i])<-inprod(beta[1:J],X[i,1:J]) + alpha[str_ID[i]] 
  }
  
  ### Priors for betas ###
  for(j in 1:J){
    beta[j]~dnorm(0, sd= 1)
  }
  
  ### Priors for random effect ###
  for(k in 1:K){
    alpha[k]~dnorm(0,sd=sqrt(1/0.000001))
  }
}
)

#' Model parameters
X<- cbind(
  rep(1, nrow(nest.data.ready)),   # Intercept (1)
  nest.data.ready$Primary,         # Distance to Primary Road
  nest.data.ready$Secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Developed,       # Developed
  nest.data.ready$Pasture,         # Pasture
  nest.data.ready$Crop,            # Crop
  nest.data.ready$Wetland,         # Wetland
  nest.data.ready$PercGrassForb,   # Percent Grass/Forb
  nest.data.ready$PercWoody,       # Percent Woody
  nest.data.ready$AvgVO,           # Visual Obstruction
  nest.data.ready$PercFern,        # Percent Fern
  nest.data.ready$Grassland,       # Grassland/Shrub
  nest.data.ready$Basal,           # Basal Area
  nest.data.ready$StemCount        # Woody Stem Count
)

ggpairs(
  as.data.frame(X), 
  upper = list(continuous = wrap("cor", size = 5)),
  diag = list(continuous = "barDiag"),
  lower = list(continuous = "points")
)

Consts <- list(str_ID = as.numeric(nest.data.ready$str_ID),
               J=ncol(X),
               I = nrow(nest.data.ready),
               K =length(unique(nest.data.ready$str_ID)))

#' Data list for nimble
Data<-list(X=X,use = nest.data.ready$Case)

Inits<-list(beta=rep(0,ncol(X)),alpha=rep(0,length(unique(nest.data.ready$str_ID))))

start <- Sys.time()

nimbleMCMC_samples <- nimbleMCMC(code = nestmodel, 
                                 constants = Consts, 
                                 inits= Inits,
                                 data = Data,
                                 nburnin = 10000,
                                 niter = 30000,
                                 thin = 3,
                                 nchains = 1)

end <- Sys.time()

#' View traceplots
MCMCtrace(nimbleMCMC_samples[,266:281], iter = 20000, pdf = F)
colMeans(nimbleMCMC_samples[,266:281])
colSds(nimbleMCMC_samples[,266:281])

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 266:281]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", 
               "Distance to Primary Road",
               "Distance to Secondary Road",
               "Mixed Forest",
               "Evergreen Forest", 
               "Developed",
               "Pasture", 
               "Crop",
               "Wetland",
               "Percent Grass/Forb",
               "Percent Woody Vegetation", 
               "Visual Obstruction",
               "Percent Fern", 
               "Grassland/Shrub",
               "Basal Area",
               "Woody Stem Count")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)


################################################################################
## Output Trace Plots

# Force conversion to plain matrix
samples_matrix <- as.matrix(nimbleMCMC_samples[, 266:281])

# Convert to coda mcmc object
mcmc_samples <- mcmc(samples_matrix)

# Assign new names
colnames(mcmc_samples) <- new_names

# Generate trace plots
MCMCtrace(mcmc_samples,
          iter = 30000,
          pdf = T)


################################################################################
## Save RData file *InsertDate_NimbleResults

################################################################################
###############################################################################X


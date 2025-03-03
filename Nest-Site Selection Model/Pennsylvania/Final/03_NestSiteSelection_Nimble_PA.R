#'---
#' title: Nest-site selection of female wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: *InsertDate*_NimbleResults.RData
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script creates a Bayesian conditional logistic regression model for nest-site selection in JAGs using the gathered covariates 
#' **Last Updated**: 1/18/25

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

#' Load in RData
load("Data Management/RData/Nest-Site Selection/Covs/Multi-State/01_Covs_AllStates_Ready.RData")

Sample_PA <- read_csv("NestingSample_PA.csv")
pa.nests.covs.1 <- right_join(pa.nests.covs, Sample_PA)

nest.data <- pa.nests.covs.1 
str(pa.nests.covs)
summary(nest.data)

nest.data <- nest.data[-577, ]

nest.data <- nest.data %>%
  dplyr::select(NestID, BandID, , PercGrassForb, PercWoody, AvgMaxVO,
                Case, Developed, Deciduous, Mixed, Evergreen, Agriculture, 
                primary, secondary, StemCount, Grassland, AvgVO, Water, Wetland, 
                PercFern) %>%
  dplyr::rename("Primary" = primary) %>%
  dplyr::rename("Secondary" = secondary)

#' Switch coding to UTF-8
nest.data <- nest.data %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

nest.data <- nest.data %>%
  dplyr::mutate(PercWoody = as.numeric(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = as.numeric(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = as.numeric(AvgMaxVO)) %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary)) %>%
  dplyr::mutate(PercFern = as.numeric(PercFern)) %>%
  dplyr::mutate(AvgVO = as.numeric(AvgVO)) %>%
  dplyr::mutate(StemCount = as.numeric(StemCount)) 
glimpse(nest.data)
str(nest.data)

nest.data <- nest.data %>%
  dplyr::mutate(StemCount = scale(StemCount)) %>%
  dplyr::mutate(PercWoody = scale( PercWoody)) %>%
  dplyr::mutate(PercGrassForb = scale(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = scale(AvgMaxVO)) %>%
  dplyr::mutate(Primary = scale(Primary)) %>%
  dplyr::mutate(Secondary = scale(Secondary)) %>%
  dplyr::mutate(PercFern = scale(PercFern)) %>%
  dplyr::mutate(AvgVO = scale(AvgVO)) 
str(nest.data)
glimpse(nest.data)

nest.data <- nest.data %>%
  dplyr::mutate(PercWoody = as.numeric(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = as.numeric(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = as.numeric(AvgMaxVO)) %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary)) %>%
  dplyr::mutate(PercFern = as.numeric(PercFern)) %>%
  dplyr::mutate(AvgVO = as.numeric(AvgVO)) %>%
  dplyr::mutate(Developed = as.numeric(Developed)) %>%
  dplyr::mutate(Agriculture = as.numeric(Agriculture)) %>%
  dplyr::mutate(Deciduous = as.numeric(Deciduous)) %>%
  dplyr::mutate(Mixed = as.numeric(Mixed)) %>%
  dplyr::mutate(Evergreen = as.numeric(Evergreen)) %>%
  dplyr::mutate(Grassland = as.numeric(Grassland)) %>%
  dplyr::mutate(Wetland = as.numeric(Wetland)) %>%
  dplyr::mutate(Water = as.numeric(Water)) %>%
  dplyr::mutate(StemCount = as.numeric(StemCount))
glimpse(nest.data)
str(nest.data)
summary(nest.data)

nest.data <- nest.data %>%
  st_drop_geometry() %>%
  dplyr::mutate(Case = as.numeric(Case))

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
glimpse(nest.data.ready)

cor(nest.data$Primary, nest.data$Secondary)


###############################################
## Nimble Model

nestmodel<-nimbleCode({
  for (i in 1:I){
    use[i]~dpois(lambda[i])
    
    log(lambda[i])<-inprod(beta[1:J],X[i,1:J])+alpha[str_ID[i]]
  }
  
  #' Priors
  for(j in 1:J){
    beta[j]~dnorm(0,sd=sqrt(1/0.0001))
  }
  
  for(k in 1:K){
    alpha[k]~dnorm(0,sd=sqrt(1/0.000001))
  }
}
)

#' Model parameters
#' Wetland is not included in the model due to a small sample size
X<- cbind(
  rep(1, nrow(nest.data.ready)),   # Intercept (1)
  nest.data.ready$Primary,         # Distance to Primary Road
  nest.data.ready$Secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Developed,       # Developed
  nest.data.ready$Agriculture,     # Agriculture
  nest.data.ready$PercGrassForb,   # Percent Grass/Forb
  nest.data.ready$PercWoody,       # Percent Woody
  nest.data.ready$StemCount,       # Woody Stem Count
  nest.data.ready$AvgVO,           # Visual Obstruction
  nest.data.ready$PercFern,        # Percent Fern
  nest.data.ready$Grassland        # Grassland/Shrub
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
                                 inits= Inits,
                                 data = Data,
                                 nburnin = 10000,
                                 niter = 30000,
                                 thin = 3,
                                 nchains = 1)

#' View traceplots
MCMCtrace(nimbleMCMC_samples, pdf = F)
colMeans(nimbleMCMC_samples[,159:171])
colSds(nimbleMCMC_samples[,159:171])

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 159:171]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", "Distance to Primary Road", "Distance to Secondary Road",
               "Mixed Forest", "Evergreen Forest", "Developed", "Agriculture", 
               "Percent Grass/Forb", "Percent Woody Vegetation", "Woody Stem Count",
               "Visual Obstruction","Percent Fern", "Grassland/Shrub")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)


##########################################
## Save RData file *InsertDate_NimbleResults

################################################################################
################################################################################


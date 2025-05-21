#'---
#' title: Daily Nest Survival Modeling of Wild Turkeys in the Mid-Atlantic Region
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script uses derived incubation start and end dates to fit a Bayesian known fate model 
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modeling daily nest survival
#' **Last Updated**: 2/19/25


################################################################################
## Load Packages

library(nimble)
library(MCMCvis)
library(tidyverse)

load("Data Management/RData/Known Fate/New Jersey/01_Covs_Ready.RData")


################################################################################
## Data Prep- Encounter Histories

#' Convert start and end of incubation days to date objects
nests.scaled.ready$startI <- as.Date(nests.scaled.ready$startI)
nests.scaled.ready$endI <- as.Date(nests.scaled.ready$endI)

#' Create a sequence of all unique dates between the min start and max end date
nests.scaled.ready %>% group_by(NestYr) %>% summarise(mindate = min(startI),
                                                maxdate = max(endI)) %>%
  st_drop_geometry()

min(nests.scaled$startI)
max(nests.scaled$endI)
inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(as.Date("2024-04-08"), as.Date("2024-06-18"), by = 1)
inc.dates <- format(inc.dates, "%m-%d")
inc.dates

#' Calculate the number of encounter days
occs <- length(inc.dates)

#' Initialize variables
first <- last <- c()
surv.caps <- matrix(data = NA, nrow = nrow(nests.scaled.ready), ncol = occs)

#' Create encounter histories within surv.caps
#' If a nest fails it gets a zero at the end 
#' If it reaches hatch it ends with a 1
for(i in 1:nrow(nests.scaled.ready)){ 
  first[i] <- which(inc.dates == format(nests.scaled.ready$startI[i], "%m-%d")) 
  last[i] <- which(inc.dates == format(nests.scaled.ready$endI[i], "%m-%d")) 
  surv.caps[i,first[i]:last[i]] <- 1 
  if(nests.scaled.ready$NestFate[i] == 0)
  {surv.caps[i,last[i]] <- 0} 
}

#' Check work on first individual
surv.caps[1,]
first;
last;

#' Change surv.caps to a df
surv.caps = as.data.frame(surv.caps)
colnames(surv.caps) <- inc.dates
surv.caps$Year = nests.scaled.ready$NestYr-2023


################################################################################
## Bring in Weather Data

getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}

f <- apply(surv.caps[,1:72],1,getFirst)
k <- apply(surv.caps[,1:72],1,getLast)
f;k
f<k

#' Subset weather data to match sequence of dates for incubation
#' 4-19 to 9-16
#' This code works for a non-leap year

#' Getting 3 day rolling mean
temperatureC2024 <- weather.array.copy[1:12,1:365,1]
precip2024 <- weather.array.copy[1:12,1:365,2]

testTempC3Roll = data.frame(matrix(NA, nrow = 12, ncol = 72))   
testPrecip3Roll = data.frame(matrix(NA, nrow = 12, ncol = 72))
for(i in 1:nrow(surv.caps)){
  if(surv.caps$Year[i] == 1){
    for(j in (f[i]+99):(k[i]+99)){
      testTempC3Roll[i,j-99] <- mean(temperatureC2024[i,((j-2):j)])
      testPrecip3Roll[i,j-99] <- mean(precip2024[i,((j-2):j)])
    }
  }
  else{
    testTempC3Roll[i,j] <- "You messed up"
    testPrecip3Roll[i,] <- "You messed up"
  }
}


#' Check work
testTempC3Roll[1,1:72]
testPrecip3Roll[1,1:72]
range(testPrecip3Roll, na.rm = T) 
range(testTempC3Roll, na.rm = T)


#' Scale weather data
#' Had issues with the scale function so I did this without packages
weather.array3 <- array(data = NA, dim = c(nrow(testTempC3Roll), ncol(testTempC3Roll),2))
testTempC.mat <- as.matrix(testTempC3Roll)
testTempC.scale <- (testTempC.mat - mean(testTempC.mat,na.rm=T))/sd(testTempC.mat,na.rm=T)
testTempC.scale
testPrecip.mat <- as.matrix(testPrecip3Roll)
testPrecip.scale <- (testPrecip.mat - mean(testPrecip.mat,na.rm=T))/sd(testPrecip.mat,na.rm=T)
testPrecip.scale
weather.array3[,,1] <- testTempC.scale
weather.array3[,,2] <- testPrecip.scale


print(weather.array3[3,,1])
print(weather.array3[2,,2])

################################################################################
## Compile Parameters into matrix format

#' Function to fill NA values with 0 throughout dataframe
#' These parameters are scaled so we are replacing NA values with the mean
fill_NA_with_value <- function(df, value = 0) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      replace(col, is.na(col), value)
    } else {
      col
    }
  })
  return(df)
}
nests.ready <- fill_NA_with_value(nests.scaled.ready)
summary(nests.ready)

#' Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),               # Intercept (1)
  nests.ready$Primary,                     # Distance to Primary Road
  nests.ready$Secondary,                   # Distance to Secondary Road
  nests.ready$Mixed,                       # Mixed Forest
  nests.ready$Evergreen,                   # Evergreen Forest
  nests.ready$Developed,                   # Developed
  nests.ready$Pasture,                     # Pasture
  nests.ready$Crop,                        # Crop
  nests.ready$Grassland,                   # Grassland
  nests.ready$`Nest Incubation Date`,      # Nest Incubation Date
  nests.ready$IncubationConstancy,         # Incubation Constancy
  nests.ready$age
) 


################################################################################
## Nimble Model

knownfate.habitat <- nimbleCode({
  
  #### Loop over individuals 
  for (i in 1:n.ind) {
    
    #### Get the first and last days a nest was active
    for(t in (first[i]+1):last[i]) {
      
      #### mu and z change if an animal is alive or dead 
      mu[i,t] <- phi[i,t]*z[i,t-1]
      
      #### Estimating daily nest survival in matrix format
      cloglog(phi[i,t]) <- inprod(beta[1:J], X[i, 1:J]) +
        weather.beta[1]*weather.array[i,t,1] + weather.beta[2]*weather.array[i,t,2] + 
        weather.beta[3]*weather.array[i,t,1]*weather.array[i,t,2]
      
      #### Likelihood of Nest Survival (Bernoulli)
      z[i,t]~dbern(mu[i, t])
      
    }
  }
  
  #### Priors for beta coefficients
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = 1)  
  }
  
  #### Priors for weather coefficients
  for (k in 1:n.weather.cov) {
    weather.beta[k] ~ dnorm(0, sd = 1)
  }
  
})

################################################################################

#' Constants
Consts <- list(
  n.ind = nrow(nests.ready),   
  X = X,                       
  J = ncol(X),                 
  first = first,               
  last = last,  
  n.weather.cov = dim(weather.array3)[3] +1,
  weather.array = weather.array3
)
surv.caps.matrix <- as.matrix(surv.caps)
colnames(surv.caps.matrix) <- NULL

#' Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  z = surv.caps.matrix               
)

#' Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J),
  weather.beta = rep(0,Consts$n.weather.cov)
)

start <- Sys.time()

#' Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = knownfate.habitat,          
  constants = Consts,          
  inits = Inits,              
  data = Data,                 
  nburnin = 3000,             
  niter = 20000,
  thin = 1
)
summary(nimbleMCMC_samples)
v = colMeans(nimbleMCMC_samples)
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:15]

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
               "Grassland/Shrub",
               "Nest Incubation Date",
               "Incubation Constancy",
               "Age Class",
               "Minimum Temperature",
               "Precipitation",
               "Precipitation * Minimimum Temperature")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)


################################################################################
## Data Prep for Beta Plot

#' Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

#' View the reshaped data
head(samples_long)

#' Calculate Bayesian credible intervals
credible_intervals <- samples_long %>%
  group_by(parameter) %>%
  summarise(
    lower = quantile(estimate, 0.05),   
    upper = quantile(estimate, 0.95),   
    .groups = 'drop'
  )

#' Calculate posterior mean for each parameter
#' Calculate 95% credible intervals using quantiles
mean_estimates <- samples_long %>%
  dplyr::group_by(parameter) %>%
  summarise(
    mean_estimate = mean(estimate),  
    lower = quantile(estimate, 0.05),  
    upper = quantile(estimate, 0.95),  
    .groups = 'drop'
  ) %>%
  dplyr::filter(parameter != "Intercept")
mean_estimates

mean_estimates <- mean_estimates %>%
  mutate(Scale = case_when(
    parameter %in% c("Percent Grass/Forb", 
                     "Percent Woody Vegetation", 
                     "Aerial Visual Obstruction",
                     "Horizontal Visual Obstruction",
                     "Percent Fern") ~ "Nest",
    parameter %in% c( "Grassland/Shrub",
                      "Mixed Forest",
                      "Evergreen Forest",
                      "Developed",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Pasture",
                      "Crop"
                      ) ~ "Landscape",
    parameter %in% c("Minimum Temperature",
                    "Precipitation",
                    "Precipitation * Minimimum Temperature"
                    ) ~ "Weather",
    TRUE ~ "Individual"
  ))

mean_estimates

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c(
                                     "Minimum Temperature",
                                     "Precipitation",
                                     "Precipitation * Minimimum Temperature",
                                     "Grassland/Shrub",
                                     "Mixed Forest",
                                     "Evergreen Forest",
                                     "Developed",
                                     "Pasture",
                                     "Crop",
                                     "Distance to Primary Road",
                                     "Distance to Secondary Road",
                                     "Nest Incubation Date",
                                     "Incubation Constancy",
                                     "Age Class"
                                     )))
mean_estimates

################################################################################
## Beta Plot 

#' Beta estimates and associated 90% credible intervals 
p2.betas <- ggplot(mean_estimates, 
                   aes(x = parameter, 
                       y = mean_estimate, 
                       color = Scale, 
                       shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Beta Estimate") +
  theme_minimal() + 
  coord_flip()+
  scale_color_manual(values = c("Individual" = "#DA6509",
                                "Nest" = "#A44200",
                                "Landscape" = "#D65F5F",
                                "Weather" = "#197278")) +  
  scale_shape_manual(values = c("Individual" = 15, 
                                "Nest" = 17,
                                "Landscape" = 16,
                                "Weather" = 18)) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 12))
    
  )
p2.betas

# Save multiple objects in a single RData file
save(p2.betas, 
     samples_df,
     mean_estimates,
     file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Known Fate/kf.betas.habitat.nj.RData", overwrite = T)

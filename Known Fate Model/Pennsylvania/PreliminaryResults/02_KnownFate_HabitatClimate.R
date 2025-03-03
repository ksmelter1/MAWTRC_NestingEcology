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

load("Data Management/RData/Known Fate/PreliminaryResults/20250219_PreliminaryResults_Habitat_Climate.RData")


################################################################################
## Data Prep- Encounter Histories

#' Convert start and end of incubation days to date objects
nests.scaled$startI <- as.Date(nests.scaled$startI)
nests.scaled$endI <- as.Date(nests.scaled$endI)

#' Create a sequence of all unique dates between the min start and max end date
nests.scaled %>% group_by(NestYr) %>% summarise(mindate = min(startI),
                                                maxdate = max(endI))

min(nests.scaled$startI)
max(nests.scaled$endI)
inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(as.Date("2022-04-19"), as.Date("2022-09-16"), by = 1)
inc.dates <- format(inc.dates, "%m-%d")
inc.dates

#' Calculate the number of encounter days
occs <- length(inc.dates)

#' Initialize variables
first <- last <- c()
surv.caps <- matrix(data = NA, nrow = nrow(nests.scaled), ncol = occs)

#' Create encounter histories within surv.caps
#' If a nest fails it gets a zero at the end 
#' If it reaches hatch it ends with a 1
for(i in 1:nrow(nests.scaled)){ 
  first[i] <- which(inc.dates == format(nests.scaled$startI[i], "%m-%d")) 
  last[i] <- which(inc.dates == format(nests.scaled$endI[i], "%m-%d")) 
  surv.caps[i,first[i]:last[i]] <- 1 
  if(nests.scaled$NestFate[i] == 0)
  {surv.caps[i,last[i]] <- 0} 
}

#' Check work on first individual
surv.caps[11,]
first;

#' Change surv.caps to a df
surv.caps = as.data.frame(surv.caps)
colnames(surv.caps) <- inc.dates
surv.caps$Year = nests.scaled$NestYr-2021



################################################################################
## Bring in Weather Data

# FEB: bind them together in an array
n.weather.cov<-2
weather.array <- array(NA, dim = c(nrow(nests.scaled),ncol(tmin),n.weather.cov))
str(weather.array)
weather.array[1:169,1:730,1] <- as.matrix(tmin)
weather.array[1:169,1:730,2] <- as.matrix(precip)
weather.array[,,1]

#' Check
weather.array[1:10,109:130,2]   # 109 is the julian day for April 19th 
                                # 259 is the julian day for September 16

#FEB: here is how these dimensions work out
# weather.array[1,,] # Nest 1, all days, all cov
# weather.array[,1,] # All nests, day 1, all cov
# weather.array[,,1] # All nests, all days, cov 1
# weather.array[1,1,] # Nest 1, day 1, all cov
# weather.array[1,,1] # Nest 1, all days, cov 1
# weather.array[,1,1] # All nests, day 1, cov 1 

getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}

f <- apply(surv.caps[,1:151],1,getFirst)
k <- apply(surv.caps[,1:151],1,getLast)
f;k
f<k

#' Subset weather data to match sequence of dates for incubation
#' 4-19 to 9-16
#' This code works for a non-leap year
temperatureC2022 <- weather.array[1:169,109:259,1]
temperatureC2023 <- weather.array[1:169,(109+365):(259+365),1]
precip2022 <- weather.array[1:169,109:259,2]
precip2023 <- weather.array[1:169,(109+365):(259+365),2]

testTempC = data.frame(matrix(NA, nrow = 169, ncol = 151))   
testPrecip = data.frame(matrix(NA, nrow = 169, ncol = 151))
for(i in 1:nrow(surv.caps)){
  if(surv.caps$Year[i] == 1){
    testTempC[i,f[i]:k[i]] <- temperatureC2022[i,f[i]:k[i]]
    testPrecip[i,f[i]:k[i]] <- precip2022[i,f[i]:k[i]]
  }
  else if(surv.caps$Year[i] == 2){
    testTempC[i,f[i]:k[i]] <- temperatureC2023[i,f[i]:k[i]]
    testPrecip[i,f[i]:k[i]] <- precip2023[i,f[i]:k[i]]
  }
  else{
    testTempC[i,] <- "You messed up" 
    testPrecip[i,] <- "You messed up"
  }
}

#' Check work
testTempC[25,1:5]
testPrecip[25,1:5]
range(testPrecip, na.rm = T) # Check this out Franny

#' Scale weather data
#' Had issues with the scale function so I did this without packages
weather.array2 <- array(data = NA, dim = c(nrow(testTempC), ncol(testTempC),2))
testTempC.mat <- as.matrix(testTempC)
testTempC.scale <- (testTempC.mat - mean(testTempC.mat,na.rm=T))/sd(testTempC.mat,na.rm=T)
testTempC.scale
testPrecip.mat <- as.matrix(testPrecip)
testPrecip.mat[is.na(testPrecip.mat)] <- 0
testPrecip.scale <- (testPrecip.mat - mean(testPrecip.mat,na.rm=T))/sd(testPrecip.mat,na.rm=T)
testPrecip.scale
weather.array2[,,1] <- testTempC.scale
weather.array2[,,2] <- testPrecip.scale

#' Check
weather.array2[25,1:10,1]


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
nests.ready <- fill_NA_with_value(nests.scaled)
summary(nests.ready)

#' Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),               # Intercept (1)
  nests.ready$AvgVO,                       #  Horizontal Visual Obstruction
  nests.ready$PercWoody,                   # Percent Woody Vegetation
  nests.ready$PercGrassForb,               # Percent Grass Forb
  nests.ready$Primary,                     # Distance to Primary Road
  nests.ready$Secondary,                   # Distance to Secondary Road
  nests.ready$Mixed,                       # Mixed Forest
  nests.ready$Evergreen,                   # Evergreen Forest
  nests.ready$Developed,                   # Developed
  nests.ready$Agriculture,                 # Agriculture
  nests.ready$Grassland,                   # Grassland
  nests.ready$PercFern,                    # Percent Fern
  nests.ready$`Nest Incubation Date`       # Nest Incubation Date
) 

GGally::ggpairs(nests.ready[,c(11:12)])
#cor.test(nests.ready$PercWoody, nests.ready$PercGrassForb)

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
        weather.beta[1]*weather.array[i,t,1] + weather.beta[2]*weather.array[i,t,2] 
      
      #### Likelihood of Nest Survival (Bernoulli)
      z[i,t]~dbern(mu[i, t])
      
    }
  }
  
  #### Priors for beta coefficients
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
  #### Priors for weather coefficients
  for (k in 1:n.weather.cov) {
    weather.beta[k] ~ dnorm(0, sd = 30)
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
  n.weather.cov = dim(weather.array2)[3],
  weather.array = weather.array2
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
               "Horizontal Visual Obstruction",
               "Percent Woody Vegetation",
               "Percent Grass/Forb",
               "Distance to Primary Road",
               "Distance to Secondary Road",
               "Mixed Forest",
               "Evergreen Forest",
               "Developed",
               "Agriculture",
               "Grassland/Shrub",
               "Percent Fern",
               "Nest Incubation Date",
               "Daily Minimum Temperature",
               "Daily Precipitation",
               "Daily Precipitation * Daily Minimimum Temperature")

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
                      "Agriculture"
                      ) ~ "Landscape",
    parameter %in% c("Daily Minimum Temperature",
                    "Daily Precipitation",
                    "Daily Precipitation * Daily Minimimum Temperature"
                    ) ~ "Climate",
    TRUE ~ "Individual"
  ))

mean_estimates

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c(
                                     "Horizontal Visual Obstruction",
                                     "Percent Fern",
                                     "Visual Obstruction",
                                     "Percent Woody Vegetation", 
                                     "Percent Grass/Forb",
                                     "Daily Minimum Temperature",
                                     "Daily Precipitation",
                                     "Daily Precipitation * Daily Minimimum Temperature",
                                     "Grassland/Shrub",
                                     "Mixed Forest",
                                     "Evergreen Forest",
                                     "Developed",
                                     "Agriculture",
                                     "Distance to Primary Road",
                                     "Distance to Secondary Road",
                                     "Nest Incubation Date"
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
                                "Climate" = "#197278")) +  
  scale_shape_manual(values = c("Individual" = 15, 
                                "Nest" = 17,
                                "Landscape" = 16,
                                "Climate" = 18)) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 12))
    
  )
p2.betas

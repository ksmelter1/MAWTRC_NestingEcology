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
library(sf)
library(GGally)
library(matrixStats)
library(tidyverse)
library(tidybayes)
library(VGAM)

load("Data Management/RData/Known Fate/Pennsylvania/PreliminaryResults/20250219_Ready.RData")


################################################################################
## Data Prep- Encounter Histories

#' Convert start and end of incubation days to date objects
nests.scaled$startI <- as.Date(nests.scaled$startI)
nests.scaled$endI <- as.Date(nests.scaled$endI)

#' Create a sequence of all unique dates between the min start and max end date
nests.scaled %>% 
  sf::st_drop_geometry() %>%
  group_by(NestYr) %>% summarise(mindate = min(startI),
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

#' n.weather.cov<-2
#' weather.array <- array(NA, dim = c(nrow(nests.scaled),ncol(tmin),n.weather.cov))
#' str(weather.array)
#' weather.array[1:169,1:730,1] <- as.matrix(tmin)
#' weather.array[1:169,1:730,2] <- as.matrix(precip)
#' weather.array[,,1]
#' 
#' #' Check
#' weather.array[1:10,109:130,2]   # 109 is the julian day for April 19th 
                                # 259 is the julian day for September 16

getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}

f <- apply(surv.caps[,1:151],1,getFirst)
k <- apply(surv.caps[,1:151],1,getLast)
f;k
f<k

#' Subset weather data to match sequence of dates for incubation
#' 4-19 to 9-16
#' This code works for a non-leap year

### Filling in daily value only
#' temperatureC2022 <- weather.array[1:158,109:259,1]
#' temperatureC2023 <- weather.array[1:158,(109+365):(259+365),1]
#' precip2022 <- weather.array[1:158,109:259,2]
#' precip2023 <- weather.array[1:158,(109+365):(259+365),2]
#' 
#' testTempC = data.frame(matrix(NA, nrow = 158, ncol = 151))   
#' testPrecip = data.frame(matrix(NA, nrow = 158, ncol = 151))
#' for(i in 1:nrow(surv.caps)){
#'   if(surv.caps$Year[i] == 1){
#'     testTempC[i,f[i]:k[i]] <- temperatureC2022[i,f[i]:k[i]]
#'     testPrecip[i,f[i]:k[i]] <- precip2022[i,f[i]:k[i]]
#'   }
#'   else if(surv.caps$Year[i] == 2){
#'     testTempC[i,f[i]:k[i]] <- temperatureC2023[i,f[i]:k[i]]
#'     testPrecip[i,f[i]:k[i]] <- precip2023[i,f[i]:k[i]]
#'   }
#'   else{
#'     testTempC[i,] <- "You messed up" 
#'     testPrecip[i,] <- "You messed up"
#'   }
#' }
#' 
#' ### Getting 3 day rolling mean
#' temperatureC2022 <- weather.array[1:158,1:365,1]
#' temperatureC2023 <- weather.array[1:158,366:730,1]
#' precip2022 <- weather.array[1:158,1:365,2]
#' precip2023 <- weather.array[1:158,366:730,2]
#' 
#' testTempC3Roll = data.frame(matrix(NA, nrow = 158, ncol = 151))   
#' testPrecip3Roll = data.frame(matrix(NA, nrow = 158, ncol = 151))
#' for(i in 1:nrow(surv.caps)){
#'   if(surv.caps$Year[i] == 1){
#'     for(j in (f[i]+108):(k[i]+108)){
#'       testTempC3Roll[i,j-108] <- mean(temperatureC2022[i,((j-2):j)])
#'       testPrecip3Roll[i,j-108] <- mean(precip2022[i,((j-2):j)])
#'     }
#'   }
#'   else if(surv.caps$Year[i] == 2){
#'     for(j in (f[i]+108):(k[i]+108)){
#'     testTempC3Roll[i,j-108] <- mean(temperatureC2023[i,((j-2):j)])
#'     testPrecip3Roll[i,j-108] <- mean(precip2023[i,((j-2):j)])
#'     }
#'   }
#'   # else{
#'   #   testTempC3Roll[i,j] <- "You messed up"
#'   #   #testPrecip[i,] <- "You messed up"
#'   # }
#' }
#' 
#' #' Check work
#' range(testPrecip3Roll, na.rm = T) 
#' 
#' 
#' 
#' #' Scale weather data
#' #' Had issues with the scale function so I did this without packages
#' weather.array2 <- array(data = NA, dim = c(nrow(testTempC3Roll), ncol(testTempC3Roll),2))
#' testTempC.mat <- as.matrix(testTempC3Roll)
#' testTempC.scale <- (testTempC.mat - mean(testTempC.mat,na.rm=T))/sd(testTempC.mat,na.rm=T)
#' testTempC.scale
#' testPrecip.mat <- as.matrix(testPrecip3Roll)
#' #testPrecip.mat[is.na(testPrecip.mat)] <- 0
#' testPrecip.scale <- (testPrecip.mat - mean(testPrecip.mat,na.rm=T))/sd(testPrecip.mat,na.rm=T)
#' testPrecip.scale
#' weather.array2[,,1] <- testTempC.scale
#' weather.array2[,,2] <- testPrecip.scale
#' 
#' #' Check
#' weather.array2[25,1:10,1]

#saveRDS(weather.array2, "Known Fate Model/Pennsylvania/Weather Arrays/weather.pa_2022_2023_3daymean.RDS")
weather.array2 <-readRDS("Known Fate Model/Pennsylvania/Weather Arrays/weather.pa_2022_2023_3daymean.RDS")

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
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
  #### Priors for weather coefficients
  for (k in 1:n.weather.cov) {
    weather.beta[k] ~ dnorm(0, sd = 30)
  }
  
})

################################################################################

#' Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),               # Intercept (1)
  nests.ready$AvgVO,                       # Visual Obstruction
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
  nests.ready$StemCount,                   # Woody Stem Count
  nests.ready$`Nest Incubation Date`,      # Nest Incubation Date
  nests.ready$IncubationConstancy,         # Incubation Constancy
  nests.ready$age                          # Age Class
) 

#' Use ggpairs to visualize correlations
ggpairs(as.data.frame(X), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))

#' Constants
Consts <- list(
  n.ind = nrow(nests.ready),   
  X = X,                       
  J = ncol(X),                 
  first = first,               
  last = last,  
  n.weather.cov = dim(weather.array2)[3]+1,
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
colMeans(nimbleMCMC_samples[,1:19])
colSds(nimbleMCMC_samples[,1:19])
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:19]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix)

#' Create a vector of new names
new_names <- c("Intercept", 
               "Visual Obstruction",
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
               "Woody Stem Count",
               "Nest Incubation Date",
               "Incubation Constancy",
               "Age Class",
               "Minimum Temperature",
               "Precipitation",
               "Precipitation * Minimum Temperature")

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
                     "Visual Obstruction",
                     "Woody Stem Count",
                     "Percent Fern") ~ "Nest",
    parameter %in% c( "Grassland/Shrub",
                      "Mixed Forest",
                      "Evergreen Forest",
                      "Developed",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Agriculture"
                      ) ~ "Landscape",
    parameter %in% c("Minimum Temperature",
                    "Precipitation",
                    "Precipitation * Minimum Temperature"
                    ) ~ "Climate",
    TRUE ~ "Individual"
  ))

mean_estimates

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c( 
                                     "Precipitation * Minimum Temperature",
                                     "Minimum Temperature",
                                     "Precipitation",
                                     "Grassland/Shrub",
                                     "Mixed Forest",
                                     "Evergreen Forest",
                                     "Developed",
                                     "Agriculture",
                                     "Distance to Primary Road",
                                     "Distance to Secondary Road",
                                     "Percent Fern",
                                     "Percent Grass/Forb",
                                     "Percent Woody Vegetation",
                                     "Visual Obstruction",
                                     "Woody Stem Count",
                                     "Incubation Constancy",
                                     "Nest Incubation Date",
                                     "Age Class"
                                     )))
mean_estimates

mean_estimates <- mean_estimates %>%
  mutate(Scale = factor(Scale, levels = c("Individual", "Nest", "Landscape", "Climate")))

################################################################################
## Beta Plot 

#' Beta estimates and associated 90% credible intervals 
p1.betas <- ggplot(mean_estimates, 
                   aes(x = parameter, 
                       y = mean_estimate, 
                       color = Scale, 
                       shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Posterior Mean") +
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
    axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.45),  
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )
p1.betas

#' Beta estimates and associated 90% credible intervals 
p1.betas <- ggplot(mean_estimates, 
                   aes(x = parameter, 
                       y = mean_estimate, 
                       color = Scale, 
                       shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Posterior Mean") +
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
    axis.title.x = element_text(margin = margin(t = 10), hjust = 0.4),
    axis.title.y = element_text(margin = margin(r = 12))
    
  )
p1.betas

#' Save multiple objects in a single RData file
#save(p1.betas, mean_estimates, file = "kf.betas.habitat.PA.RData")


################################################################################
## DNS Plots

samples_df$prediction <- map2(samples_df$`Visual Obstruction`,
                           samples_df$`Percent Woody Vegetation`,
                           samples_df$`Percent Grass/Forb`,
                           samples_df$`Percent Fern`,
                           samples_df$`Distance to Primary Road`,
                           samples_df$`Distance to Secondary Road`,
                           samples_df$`Woody Stem Count`,
                           samples_df$`Nest Incubation Date`,
                           samples_df$Precipitation,
                           samples_df$`Minimum Temperature`)
                           
                           function(cov1, 
                                    cov2, 
                                    cov3, 
                                    cov4, 
                                    cov5, 
                                    cov6,
                                    cov7, 
                                    cov8, 
                                    cov9,
                                    cov10){
                             clogloglink(beta0 + 
                                           beta1 * cov1 + 
                                           beta2 * cov2 + 
                                           beta3 * cov3 +
                                           beta4 * cov4 +
                                           beta5 * cov5 +
                                           beta6 * cov6 +
                                           beta7 * cov7 +
                                           beta8 * cov8 +
                                           beta9 * cov9 +
                                           beta10* cov10, inverse = TRUE)
                           }

sim_tbl_long <- unnest(samples_df, prediction) #note how there's now one row per sample * combination of covariates

# Code to create a plot

ggplot(sim_tbl_long, aes(x = parasite_diversity, y = prediction)) +
  stat_lineribbon(.width = c(.95, .8, .5), alpha = 0.5, point_interval = mean_qi) +
  scale_fill_manual(values = rev(c("#636363", "#bdbdbd", "#f0f0f0"))) +
  facet_wrap(vars(LPDV)) +
  theme_minimal()


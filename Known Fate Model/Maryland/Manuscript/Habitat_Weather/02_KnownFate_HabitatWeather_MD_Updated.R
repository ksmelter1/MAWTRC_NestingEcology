#---
# title: Daily Nest Survival Modeling of Female Wild Turkeys in Maryland
# authors: K. Smelter
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
#' **Purpose**: This script uses derived incubation start and end dates to fit a Bayesian known fate model 
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modeling daily nest survival
#' **Last Updated**: 12/27/25


################################################################################
## Load Packages

library(nimble)
library(MCMCvis)
library(tidyverse)
library(GGally)
library(sf)
library(matrixStats)
library(tidybayes)
library(VGAM)
library(coda)

load("Data Management/RData/Known Fate/Maryland/Manuscript/Habitat_Weather/01_Covs_Ready_Updated.RData")


################################################################################
## Sample size metrics for manuscript

nests.scaled <- nests.scaled.ready %>%
  st_drop_geometry()
colnames(nests.scaled)

# Obtain sample sizes for manuscript
# Numbers don't exactly add up because a hen could nest as a juvenile one year and adult the next
# Unique Birds 
length(unique(nests.scaled$BandID))
length(unique(nests.scaled$BandID[nests.scaled$age == "0"]))
length(unique(nests.scaled$BandID[nests.scaled$age == "1"]))

# Nesting Birds
length(nests.scaled$BandID)
length(nests.scaled$BandID[nests.scaled$age == "0"])
length(nests.scaled$BandID[nests.scaled$age == "1"])

# Ensure startI is a Date type
nests.scaled$startI <- as.Date(nests.scaled$startI)

# Subset for 2023
median_2023 <- median(nests.scaled$startI[nests.scaled$NestYr == 2023], na.rm = TRUE)

# Subset for 2024
median_2024 <- median(nests.scaled$startI[nests.scaled$NestYr == 2024], na.rm = TRUE)

# Print results
median_2023
median_2024

# Calculate average number of nesting attempts per BandID
avg_attempts <- nests.scaled %>%
  group_by(BandID) %>%
  summarise(n_attempts = n_distinct(NestID1, na.rm = TRUE)) %>%
  summarise(mean_attempts = mean(n_attempts, na.rm = TRUE))

avg_attempts

################################################################################
## Data Prep- Encounter Histories

# Convert start and end of incubation days to date objects
nests.scaled$startI <- as.Date(nests.scaled$startI)
nests.scaled$endI <- as.Date(nests.scaled$endI)

# Create a sequence of all unique dates between the min start and max end date
nests.scaled %>% group_by(NestYr) %>% summarise(mindate = min(startI),
                                                maxdate = max(endI))
inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(as.Date("2024-04-15"), as.Date("2024-07-13"), by = 1)
inc.dates <- format(inc.dates, "%m-%d")
inc.dates

# Calculate the number of encounter days
occs <- length(inc.dates)

# Initialize variables
first <- last <- c()
surv.caps <- matrix(data = NA, nrow = nrow(nests.scaled), ncol = occs)

# Create encounter histories within surv.caps
# If a nest fails it gets a zero at the end 
# If it reaches hatch it ends with a 1
for(i in 1:nrow(nests.scaled)){ 
  first[i] <- which(inc.dates == format(nests.scaled$startI[i], "%m-%d")) 
  last[i] <- which(inc.dates == format(nests.scaled$endI[i], "%m-%d")) 
  surv.caps[i,first[i]:last[i]] <- 1 
  if(nests.scaled$NestFate[i] == 0)
  {surv.caps[i,last[i]] <- 0} 
}

# Check work on first individual
surv.caps[11,]
first;

# Change surv.caps to a df
surv.caps = as.data.frame(surv.caps)
colnames(surv.caps) <- inc.dates
surv.caps$Year = nests.scaled$NestYr-2022

# Functions from Kery and Schaub book to get the first and last encounter occasion
getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}

# Apply the functions
f <- apply(surv.caps[,1:90],1,getFirst)
k <- apply(surv.caps[,1:90],1,getLast)
f;k
f<k

################################################################################
## Bring in Weather Data

# Getting 3 day rolling mean
temperatureC2023 <- weather.array.copy[1:97,1:365,1]
temperatureC2024 <- weather.array.copy[1:97,366:730,1]
precip2023 <- weather.array.copy[1:97,1:365,2]
precip2024 <- weather.array.copy[1:97,366:730,2]
 
 testTempC3Roll = data.frame(matrix(NA, nrow = 97, ncol = 90))   
 testPrecip3Roll = data.frame(matrix(NA, nrow = 97, ncol = 90))
 for(i in 1:nrow(surv.caps)){
 if(surv.caps$Year[i] == 1){
  for(j in (f[i]+104):(k[i]+104)){
 testTempC3Roll[i,j-104] <- mean(temperatureC2023[i,((j-2):j)])
 testPrecip3Roll[i,j-104] <- mean(precip2023[i,((j-2):j)])
     }
   }
 else if(surv.caps$Year[i] == 2){
      for(j in (f[i]+105):(k[i]+105)){
      testTempC3Roll[i,j-105] <- mean(temperatureC2024[i,((j-2):j)])
      testPrecip3Roll[i,j-105] <- mean(precip2024[i,((j-2):j)])
      }
    }
     else{
     testTempC3Roll[i,j] <- "You messed up"
    testPrecip3Roll[i,] <- "You messed up"
     }
  }
 
testPrecip3Roll[91,first[91]:90]
testTempC3Roll[91,first[91]:90]

# Check work
testTempC3Roll[1,1:90]
testPrecip3Roll[1,1:90]
range(testPrecip3Roll, na.rm = T) 
range(testTempC3Roll, na.rm = T)
mean(unlist(testTempC3Roll), na.rm = TRUE)
mean(unlist(testPrecip3Roll), na.rm = TRUE)
quantile(unlist(testPrecip3Roll), probs = seq(0,1,0.1), na.rm = TRUE)

sd(unlist(testTempC3Roll), na.rm = TRUE)
sd(unlist(testPrecip3Roll), na.rm = TRUE)

(0-2.89)/4.09


# Scale weather data
# Had issues with the scale function so I did this without packages
weather.array3 <- array(data = NA, dim = c(nrow(testTempC3Roll), ncol(testTempC3Roll),2))
testTempC.mat <- as.matrix(testTempC3Roll)
testTempC.scale <- (testTempC.mat - mean(testTempC.mat,na.rm=T))/sd(testTempC.mat,na.rm=T)
testTempC.scale
testPrecip.mat <- as.matrix(testPrecip3Roll)
testPrecip.scale <- (testPrecip.mat - mean(testPrecip.mat,na.rm=T))/sd(testPrecip.mat,na.rm=T)
testPrecip.scale
weather.array3[,,1] <- testTempC.scale
weather.array3[,,2] <- testPrecip.scale

# Check
print(weather.array3[3,,1])
print(weather.array3[2,,2])

weather.array3[91,first[91]:90, 1]
weather.array3[91,first[91]:90, 2]

################################################################################
## Compile Parameters into matrix format

# Check for NAs
summary(nests.scaled)

# Since no NAs no need to fill non-categorical NA values with the mean
nests.ready <- nests.scaled

# Model parameters
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
  nests.ready$Wetland,                     # Wetland
  nests.ready$`Nest Incubation Date`,      # Nest Incubation Date
  nests.ready$IncubationConstancy,         # Incubation Constancy
  nests.ready$age                          # Age 
) 

# Use ggpairs to visualize correlations
ggpairs(as.data.frame(X), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))


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

# Constants
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

# Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  z = surv.caps.matrix               
)

# Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J),
weather.beta = rep(0,Consts$n.weather.cov)
)

start <- Sys.time()

# Run MCMC sampling with Nimble
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
colMeans(nimbleMCMC_samples[,1:16])
colSds(nimbleMCMC_samples[,1:16])
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

# Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:16]

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
               "Wetland",
               "Nest Incubation Date",
               "Incubation Constancy", 
               "Age Class",
               "Minimum Temperature",
               "Precipitation",
               "Precipitation * Minimum Temperature")

# Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

# View the renamed data frame
head(samples_df)

################################################################################
## Output Trace Plots

# Force conversion to plain matrix
samples_matrix <- as.matrix(nimbleMCMC_samples[, 1:16])

# Convert to coda mcmc object
mcmc_samples <- mcmc(samples_matrix)

# Assign new names
colnames(mcmc_samples) <- new_names

# Generate trace plots
MCMCtrace(mcmc_samples,
          iter = 20000,
          pdf = T)



################################################################################
## Data Prep for Beta Plot

# Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

# View the reshaped data
head(samples_long)

# Calculate posterior mean for each parameter
# Calculate 95% credible intervals using quantiles
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
                      "Wetland",
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
                    "Precipitation * Minimum Temperature"
                    ) ~ "Climate",
    TRUE ~ "Individual"
  ))

mean_estimates

# Organize variables into levels to be displayed
# Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c(
                                     "Precipitation * Minimum Temperature",
                                     "Minimum Temperature",
                                     "Precipitation",
                                     "Wetland",
                                     "Grassland/Shrub",
                                     "Mixed Forest",
                                     "Evergreen Forest",
                                     "Developed",
                                     "Pasture",
                                     "Crop",
                                     "Distance to Primary Road",
                                     "Distance to Secondary Road",
                                     "Incubation Constancy",
                                     "Nest Incubation Date",
                                     "Age Class"
                                     )))
mean_estimates
mean_estimates <- mean_estimates %>%
  mutate(Scale = factor(Scale, levels = c("Individual", 
                                          "Nest",
                                          "Landscape",
                                          "Climate")))

################################################################################
## Beta Plot 

# Beta estimates and associated 90% credible intervals 
p2.betas <- ggplot(mean_estimates, 
                   aes(x = parameter, 
                       y = mean_estimate, 
                       color = Scale, 
                       shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = " Posterior Mean") +
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

# Save multiple objects in a single RData file
save(p2.betas, samples_df, mean_estimates, file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Known Fate/kf.betas.habitat.MD.updated.RData")


################################################################################
## DNS Plots

# For each prediction plot update the covariate within colnames
# Include the intercept and multiply the beta by the predictor
# Create plots below by updating the x-axis with each predictor
sim_tbl <- expand_grid(-1:1)
colnames(sim_tbl) <- c("Distance to Secondary Road") 


sim_tbl$prediction <- pmap(list(sim_tbl$`Distance to Secondary Road`), function(cov1){
  clogloglink(samples_df$Intercept + 
                samples_df$`Distance to Secondary Road` * cov1, 
              inverse = TRUE)
})
sim_tbl_long <- unnest(sim_tbl, prediction) 

pred.p3 <- sim_tbl_long %>%
  ggplot(aes(x = `Distance to Secondary Road`, y = prediction)) +
  stat_lineribbon(.width = 0.9, alpha = 0.8, point_interval = mean_qi) +  # Set to 0.9 for 90% credible interval
  scale_fill_manual(values = rev(c("#636363", "#bdbdbd", "#f0f0f0"))) +
  ylab("Daily Nest Survival Probability") +
  theme_minimal() +
  ggtitle("") +
  theme(
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.title.x = element_text(margin = margin(t = 15))
  )
pred.p3

################################################################################
## Visualize the Interaction between Precipitation and Minimum Temperature

library(ggplot2)
library(tidyverse)
library(nimble)

# Obtain means, sds of weather variables for back transforming 
mean(unlist(testTempC3Roll), na.rm = TRUE)
mean(unlist(testPrecip3Roll), na.rm = TRUE)
sd(unlist(testTempC3Roll), na.rm = TRUE)
sd(unlist(testPrecip3Roll), na.rm = TRUE)

# Define standardized temp values
temp_vals <- seq(-2, 2, length.out = 100)

# Define low, medium, high standardized precipitation levels
# Back-transform: standardized value * sd + mean
quantile(unlist(testPrecip3Roll), probs = seq(0,1,0.1), na.rm = TRUE)
quantile(unlist(testPrecip3Roll), probs = seq(0,1,0.05), na.rm = TRUE)
precip_vals <- c(-0.7066015, 0, 0.3202135) # 0mm, 2.88mm, 4.1991667 (0-35%, mean, 75% percentile)

# Posterior samples
beta_temp <- samples_df$`Minimum Temperature`
beta_precip <- samples_df$`Precipitation`
beta_inter <- samples_df$`Precipitation * Minimum Temperature`
inter <- samples_df$Intercept

# Container for predictions
prediction_list <- list()

# Loop over each fixed precip level
for (p in precip_vals) {
  for (t in temp_vals) {
    lp <- inter + beta_temp * t + beta_precip * p + beta_inter * t * p
    phi <- icloglog(lp)
    
    prediction_list[[length(prediction_list) + 1]] <- data.frame(
      tempscaled = t,
      temp = (t*4.68)+12.83,
      precip = p,
      phi_mean = mean(phi),
      phi_lower = quantile(phi, 0.05),
      phi_upper = quantile(phi, 0.95)
    )
  }
}

# Combine predictions and label precip levels
pred_df <- do.call(rbind, prediction_list)
pred_df$precip <- factor(pred_df$precip,
                         levels = precip_vals,
                         labels = c("Low Precip", "Medium Precip", "High Precip"))

ggplot(pred_df, aes(x = temp, y = phi_mean, color = precip, fill = precip, group = precip)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = phi_lower, ymax = phi_upper), alpha = 0.25, color = NA) +
  labs(
    x = "Minimum Temperature (°C)",
    y = "Predicted daily nest survival probability",
    color = "Precipitation Level",
    fill = "Precipitation Level"
  ) +
  theme_light() +
  theme(legend.position = 'bottom')

################################################################################
###############################################################################X
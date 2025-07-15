#---
# title: Daily Nest Survival Modeling of Wild Turkeys in the Mid-Atlantic Region
# authors: "K. Smelter
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
# **Purpose**: This script uses derived incubation start and end dates to fit a Bayesian known fate model 
# **Key Changes**: This script incorporates 2024 data into the Pennsylvania known fate 
# **Last Updated**: 7/14/25

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
library(ggpubr)

# RData file from data prep script
# Contains data to fit model
load("Data Management/RData/Known Fate/Pennsylvania/Manuscript/01_Covs_Ready.updated.RData")


################################################################################
## Data Prep- Encounter Histories

# Convert start and end of incubation days to date objects
nests.scaled$startI <- as.Date(nests.scaled$startI)
nests.scaled$endI <- as.Date(nests.scaled$endI)

# Create a sequence of all unique dates between the min start and max end date
nests.scaled %>% 
  sf::st_drop_geometry() %>%
  group_by(NestYr) %>% summarise(mindate = min(startI),
                                                maxdate = max(endI))
inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(as.Date("2022-04-12"), as.Date("2022-09-16"), by = 1)
inc.dates <- format(inc.dates, "%m-%d")
inc.dates

# Calculate the number of encounter days
occs <- length(inc.dates)

# Initialize variables
# Create matrix for encounter histories
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
surv.caps[12,]
first;
last;

# Change surv.caps to a df
surv.caps = as.data.frame(surv.caps)
colnames(surv.caps) <- inc.dates
surv.caps$Year = nests.scaled$NestYr-2021


################################################################################
## Bring in Weather Data

# Functions from Kery and Schaub to get the first and last encounter occasions
# Make sure that f is always greater than k (First encounter always comes before last encounter)
getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}
f <- apply(surv.caps[,1:158],1,getFirst)
k <- apply(surv.caps[,1:158],1,getLast)
f;k
f<k

# Store weather data within arrays for 2022, 2023, and 2024
# Rows = 260 (# of nests)
# Columns = # of days within 3 years
# Dimensions = 1, or 2 (Temp = 1, Precip = 2)
temperatureC2022 <- weather.array.copy[1:260,1:365,1]
temperatureC2023 <- weather.array.copy[1:260,366:730,1]
temperatureC2024 <- weather.array.copy[1:260, 731:1095,1]
precip2022 <- weather.array.copy[1:260,1:365,2]
precip2023 <- weather.array.copy[1:260,366:730,2]
precip2024 <- weather.array.copy[1:260, 731:1095,1]

################################################################################
## Get 3 day moving average

# Convert from matrix to dataframe
# For each row in surv.caps if the year is the first year of the study, collect the rolling mean for 2022 variables
# For each row in surv.caps if the year is the second year of the study, collect the rolling mean for 2023 variables
testTempC3Roll = data.frame(matrix(NA, nrow = 260, ncol = 158))   
testPrecip3Roll = data.frame(matrix(NA, nrow = 260, ncol = 158))
for(i in 1:nrow(surv.caps)){
  if(surv.caps$Year[i] == 1){
    for(j in (f[i]+102):(k[i]+102)){
       testTempC3Roll[i,j-102] <- mean(temperatureC2022[i,((j-2):j)])
       testPrecip3Roll[i,j-102] <- mean(precip2022[i,((j-2):j)])
     }
   }
   else if(surv.caps$Year[i] == 2){
      for(j in (f[i]+102):(k[i]+102)){
     testTempC3Roll[i,j-102] <- mean(temperatureC2023[i,((j-2):j)])
    testPrecip3Roll[i,j-102] <- mean(precip2023[i,((j-2):j)])
    }
   }
  
  else if(surv.caps$Year[i] == 3){
    for(j in (f[i]+102):(k[i]+102)){
      testTempC3Roll[i,j-102] <- mean(temperatureC2023[i,((j-2):j)])
      testPrecip3Roll[i,j-102] <- mean(precip2023[i,((j-2):j)])
    }
  }
   else{
  testTempC3Roll[i,j] <- "You messed up"
  testPrecip3Roll[i,] <- "You messed up"
}
 }

# Check work
range(testPrecip3Roll, na.rm = T) 
range(testTempC3Roll, na.rm = T)

# Scale weather data
# Scale function had problems with the NAs within the matrices
weather.array2 <- array(data = NA, dim = c(nrow(testTempC3Roll), ncol(testTempC3Roll),2))
testTempC.mat <- as.matrix(testTempC3Roll)
testTempC.scale <- (testTempC.mat - mean(testTempC.mat,na.rm=T))/sd(testTempC.mat,na.rm=T)
testTempC.scale
testPrecip.mat <- as.matrix(testPrecip3Roll)
testPrecip.scale <- (testPrecip.mat - mean(testPrecip.mat,na.rm=T))/sd(testPrecip.mat,na.rm=T)
testPrecip.scale
weather.array2[,,1] <- testTempC.scale
weather.array2[,,2] <- testPrecip.scale

# Check work
weather.array2[12,,1]
weather.array2[1,1:100,2]

# Check correlations
temp_flat <- as.vector(weather.array2[,,1])
precip_flat <- as.vector(weather.array2[,,2])

# Compute correlation, excluding NA pairs
cor_result <- cor(temp_flat, precip_flat, use = "complete.obs")

# Print correlation
cor_result


################################################################################
## Compile Parameters into matrix format

# Function to fill NA values with 0 throughout dataframe
# These parameters are scaled so we are replacing NA values with the mean
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

# Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),               # Intercept (1)
  nests.ready$AvgVO,                       # Visual Obstruction
  nests.ready$StemCount,                   # Woody Stem Count
  nests.ready$Basal,                       # Basal Area
  nests.ready$PercWoody,                   # Percent Woody Vegetation
  nests.ready$PercGrassForb,               # Percent Grass Forb
  nests.ready$Primary,                     # Distance to Primary Road
  nests.ready$Secondary,                   # Distance to Secondary Road
  nests.ready$Mixed,                       # Mixed Forest
  nests.ready$Evergreen,                   # Evergreen Forest
  nests.ready$Developed,                   # Developed
  nests.ready$Pasture,                     # Pasture
  nests.ready$Crop,                        # Crop
  nests.ready$Grassland,                   # Grassland
  nests.ready$Wetland,                     # Wetland
  nests.ready$PercFern,                    # Percent Fern
  nests.ready$`Nest Incubation Date`,      # Nest Incubation Date
  nests.ready$IncubationConstancy,         # Incubation Constancy
  nests.ready$age                          # Age Class
) 

# Use ggpairs to visualize correlations
ggpairs(as.data.frame(X), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))

# Constants
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

# Check outputs and convergence
summary(nimbleMCMC_samples)
colMeans(nimbleMCMC_samples[,1:22])
colSds(nimbleMCMC_samples[,1:22])
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

# Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:22]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix)

# Create a vector of new names
new_names <- c("Intercept", 
               "Visual Obstruction",
               "Woody Stem Count",
               "Basal Area",
               "Percent Woody Vegetation",
               "Percent Grass/Forb",
               "Distance to Primary Road",
               "Distance to Secondary Road",
               "Mixed Forest",
               "Evergreen Forest",
               "Developed",
               "Pasture",
               "Crop",
               "Grassland/Shrub",
               "Wetland",
               "Percent Fern",
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
## Data Prep for Beta Plot

# Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

# View the reshaped data
head(samples_long)

# Calculate posterior mean for each parameter
# Calculate 90% credible intervals using quantiles
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

# Assign scales to variables using case_when
# Individual
# Nest
# Landscape
# Weather
mean_estimates <- mean_estimates %>%
  mutate(Scale = case_when(
    parameter %in% c("Percent Grass/Forb", 
                     "Percent Woody Vegetation", 
                     "Visual Obstruction",
                     "Woody Stem Count",
                     "Basal Area",
                     "Percent Fern") ~ "Nest",
    parameter %in% c( "Grassland/Shrub",
                      "Mixed Forest",
                      "Evergreen Forest",
                      "Developed",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Pasture",
                      "Crop",
                      "Wetland"
                      ) ~ "Landscape",
    parameter %in% c("Minimum Temperature",
                    "Precipitation",
                    "Precipitation * Minimum Temperature"
                    ) ~ "Weather",
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
                                     "Percent Fern",
                                     "Percent Grass/Forb",
                                     "Percent Woody Vegetation",
                                     "Visual Obstruction",
                                     "Basal Area",
                                     "Woody Stem Count",
                                     "Incubation Constancy",
                                     "Nest Incubation Date",
                                     "Cumulative Distance Traveled",
                                     "Age Class"
                                     )))
mean_estimates

# Assign levels for figure 
mean_estimates <- mean_estimates %>%
  mutate(Scale = factor(Scale, levels = c("Individual", "Nest", "Landscape", "Weather")))


################################################################################
## Beta Plot--- Presentations

# Beta estimates and associated 90% credible intervals 
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
                                "Weather" = "#197278")) +  
  scale_shape_manual(values = c("Individual" = 15, 
                                "Nest" = 17,
                                "Landscape" = 16,
                                "Weather" = 18)) +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.45),  
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )
p1.betas


################################################################################
## Beta Plot--- Journal 

# Beta estimates and associated 90% credible intervals 
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
                                "Weather" = "#197278")) +  
  scale_shape_manual(values = c("Individual" = 15, 
                                "Nest" = 17,
                                "Landscape" = 16,
                                "Weather" = 18)) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), hjust = 0.4),
    axis.title.y = element_text(margin = margin(r = 12))
    
  )
p1.betas


################################################################################
## Output Data

#Save multiple objects in a single RData file
#save(p1.betas, mean_estimates, samples_df, file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Known Fate/kf.betas.habitat.PA.updated.RData", overwrite = T)


################################################################################
## DNS Plots

# For each prediction plot update the covariate within colnames
# Include the intercept and multiply the beta by the predictor
# Create plots below by updating the x-axis with each predictor

sim_tbl <- expand_grid(-1:1)
colnames(sim_tbl) <- c("Evergreen Forest") 

sim_tbl$prediction <- pmap(list(sim_tbl$`Evergreen Forest`), function(cov1){
                             clogloglink(samples_df$Intercept + 
                                         samples_df$`Evergreen Forest` * cov1, 
                                         inverse = TRUE)
                           })
sim_tbl_long <- unnest(sim_tbl, prediction) 

pred.p1 <- sim_tbl_long %>%
  ggplot(aes(x = `Visual Obstruction`, y = prediction)) +
  stat_lineribbon(.width = 0.9, alpha = 0.8, point_interval = mean_qi) +
  scale_fill_manual(values = c("#bdbdbd")) +
  ylab("Daily nest survival probability") +
  xlab("Visual obstruction") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
    legend.position = 'none'
  )
pred.p1

pred.p2 <- sim_tbl_long %>%
  ggplot(aes(x = `Nest Incubation Date`, y = prediction)) +
  stat_lineribbon(.width = 0.9, alpha = 0.8, point_interval = mean_qi) +
  scale_fill_manual(values = c("#bdbdbd")) +
  ylab("Daily nest survival probability") +
  xlab("Nest incubation date") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
    legend.position = 'none'
  )
pred.p2

pred.p3 <- sim_tbl_long %>%
  ggplot(aes(x = `Incubation Constancy`, y = prediction)) +
  stat_lineribbon(.width = 0.9, alpha = 0.8, point_interval = mean_qi) +
  scale_fill_manual(values = c("#bdbdbd")) +
  ylab("Daily nest survival probability") +
  xlab("Incubation constancy") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
    legend.position = 'none'
  )
pred.p3

pred.p4 <- sim_tbl_long %>%
  ggplot(aes(x = `Evergreen Forest`, y = prediction)) +
  stat_lineribbon(.width = 0.9, alpha = 0.8, point_interval = mean_qi) +
  scale_fill_manual(values = c("#bdbdbd")) +
  ylab("Daily nest survival probability") +
  xlab("Evergreen forest") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
    legend.position = 'none'
  )
pred.p4


################################################################################
###############################################################################X

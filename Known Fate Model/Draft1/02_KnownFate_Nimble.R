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
#' **Last Updated**: 1/23/25


#####################################################
## Load Packages

library(nimble)
library(MCMCvis)
library(tidyverse)


#####################################################
## Data Prep- Encounter Histories

#' Convert start and end of incubation days to date objects
nests.scaled$startI <- as.Date(nests.scaled$startI)
nests.scaled$endI <- as.Date(nests.scaled$endI)

#' Create a sequence of all unique dates between the min start and max end date
inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(inc.dates[1], inc.dates[length(inc.dates)], by = 1)

#' Calculate the number of encounter days
occs <- length(inc.dates)

#' Initialize variables
first <- last <- array(NA, dim = nrow(nests.scaled))
surv.caps <- matrix(data = NA, nrow = nrow(nests.scaled), ncol = occs)

#' Create encounter histories within surv.caps
#' If a nest fails it gets a zero at the end 
#' If it reaches hatch it ends with a 1
for(i in 1:nrow(nests.scaled)){ #for each individual 
  first[i] <- which(inc.dates == nests.scaled$startI[i]) 
  last[i] <- which(inc.dates == nests.scaled$endI[i]) 
  surv.caps[i,first[i]:last[i]] <- 1 
  if(nests.scaled$NestFate[i] == 0)
  {surv.caps[i,last[i]] <- 0} 
}

#' Check work on first individual
surv.caps[1,]

#####################################################################
## Attempt at creating an array

#' Max number of nests a bird had
#' nn.nests <- nests.scaled %>%
#'   dplyr::group_by(BandID) %>%
#'   summarize(n = n())
#' max(nn.nests$n)
#' 
#' #' Initialize array
#' n.hens <- length(unique(nests.scaled$BandID))
#' n.nests <- max(nn.nests$n)
#' n.days <- max(occs)
#' z <- array(NA, dim = c(n.hens, n.nests, n.days))
#' 
#' #' Fill the array
#' for (h in 1:n.hens) {
#'   ## grab single hen
#'   tmp.hen <- nests.scaled[which(nests.scaled$BandID == unique(nests.scaled$BandID)[h]), ]
#'   n.nests <- unique(tmp.hen$NestID)
#'   
#'   for (q in 1:length(n.nests)) {
#'     tmp.hen.nests <- tmp.hen[which(tmp.hen$NestID == n.nests[n]), ]
#'     n.days <- tmp.hen.nests %>%
#'       dplyr::mutate(interval = length(seq.Date(startI, endI, by = "day")))
#'     
#'     for (d in 1:n.days) {
#'       
#'       # Assign filtered data to z[h, q, d]
#'       z[h, q, d] <- nests.scaled %>%
#'         dplyr::filter(BandID == unique(nests.scaled$BandID)[h]) %>%
#'         dplyr::filter(NestID == unique(NestID)[q]) %>%
#'         dplyr::mutate(occs == length(nests$endI))
#'     }
#'   }
#' }

#' Add 0s and 1s for each column 
#' rbind for each hen and nest

##################################################
## Simulate Temperature Data

#' Simulate daily temperature values for the length of occs
set.seed(123) 
mean_temp <- 60   
sd_temp <- 10   

#' Simulate temperatures for the entire sequence of days
daily_temps <- rnorm(occs, mean = mean_temp, sd = sd_temp)

#' Create a matrix to match surv.caps dimensions
temp.values <- matrix(NA, nrow = nrow(nests), ncol = occs)

#' Populate temp.values with the simulated daily temperatures
for(i in 1:nrow(nests)) {
  temp.values[i, first[i]:last[i]] <- daily_temps[first[i]:last[i]]
}

###################################################
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
  nests.ready$`Visual Obstruction`,        # Visual Obstruction
  nests.ready$`Percent Woody Vegetation`,  # Percent Woody Vegetation
  nests.ready$`Percent Grass/Forb`,        # Percent Grass Forb
  nests.ready$`Basal Area`,                # Basal Area
  nests.ready$Primary,                     # Distance to Primary Road
  nests.ready$Secondary,                   # Distance to Secondary Road
  nests.ready$Mixed,                       # Mixed Forest
  nests.ready$Evergreen,                   # Evergreen Forest
  nests.ready$Deciduous,                   # Deciduous Forest
  nests.ready$Agriculture,                 # Agriculture
  nests.ready$Grassland,                   # Grassland
  nests.ready$`Nest Incubation Date`,      # Nest Incubation Date
  nests.ready$age                          # Age Class
)   


#################################################
## Nimble Model

knownfate <- nimbleCode({
  
  #### Loop over individuals 
  for (i in 1:n.ind) {
    
    #### Get the first and last days a nest was active
    for(t in (first[i]+1):last[i]) {
      
      #### mu and z change if an animal is alive or dead 
      mu[i,t] <- phi[i,t]*z[i,t-1]
      
      #### Estimating daily nest survival in matrix format
      cloglog(phi[i,t]) <- inprod(beta[1:J], X[i, 1:J]) 
      
      #### Likelihood of Nest Survival (Bernoulli)
      z[i,t]~dbern(mu[i, t])
      
    }
  }
  
  #### Priors for beta coefficients
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
})

#####################################################

#' Constants
Consts <- list(
  n.ind = nrow(nests.ready),   
  X = X,                       
  J = ncol(X),                 
  first = first,               
  last = last                  
)

#' Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  z = surv.caps               
)

#' Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J)   
)

start <- Sys.time()

#' Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = knownfate,          
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
beta_samples <- nimbleMCMC_samples[, 1:14]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", "Distance to Primary Road", "Distance to Secondary Road", 
               "Mixed Forest", "Evergreen Forest", "Deciduous Forest", "Agriculture", 
               "Visual Obstruction", "Percent Grass/Forb", 
               "Percent Woody Vegetation", "Basal Area", "Grassland/Shrub", "Nest Incubation Date", "Age Class")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)

############################
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
    parameter %in% c("Percent Grass/Forb", "Percent Woody Vegetation", "Visual Obstruction", 
                     "Basal Area", "Native Woody Vegetation", "Invasive Woody Vegetation") ~ "Nest-Level",
    parameter %in% c( "Grassland/Shrub",
                      "Mixed Forest",
                      "Evergreen Forest",
                      "Deciduous Forest",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Agriculture") ~ "Landscape-Level",
    TRUE ~ "Individual"
  ))

mean_estimates

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c(
                                     "Grassland/Shrub",
                                     "Mixed Forest",
                                     "Evergreen Forest",
                                     "Deciduous Forest",
                                     "Distance to Primary Road",
                                     "Distance to Secondary Road",
                                     "Agriculture",
                                     "Visual Obstruction",
                                     "Percent Woody Vegetation", 
                                     "Percent Grass/Forb",
                                     "Basal Area",
                                     "Nest Incubation Date",
                                     "Age Class"))) 
mean_estimates

#########################
## Beta Plot 

#' Beta estimates and associated 95% credible intervals 
#' Macroscale and Microscale predictors
p2.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5) +  # Points for the mean estimate 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  # Error bars for credible intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at 0
  labs(x = "Parameter", y = "Beta Estimate") +
  theme_minimal() + 
  coord_flip()+
  scale_color_manual(values = c("Individual" = "#fb8b24", "Nest-Level" = "#A44200", "Landscape-Level" = "#D65F5F")) +  # Set color for Microscale, Macroscale
  scale_shape_manual(values = c("Individual" = 15, "Nest-Level" = 17, "Landscape-Level" = 16))+  # Set shapes for Microscale, Macroscale
  theme(
    axis.title.x = element_text(margin = margin(t = 10))  # Pad the x-axis label by 10 points (~0.1 inch)
  )
p2.betas

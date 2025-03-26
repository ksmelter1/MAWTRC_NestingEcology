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
#' **Last Updated**: 3/13/25


################################################################################
## Load Packages

library(nimble)
library(MCMCvis)
library(tidyverse)
library(matrixStats)

load("Data Management/RData/Known Fate/Pennsylvania/PreliminaryResults/20250219_PreliminaryResults_Individual.RData")


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

# Subset nests.scaled to remove nests before 4-22 regardless of year
nests.scaled <- nests.scaled[format(nests.scaled$startI, "%m-%d") >= "04-22",]

inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(as.Date("2022-04-22"), as.Date("2022-09-16"), by = 1)
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

surv.caps = as.data.frame(surv.caps)
colnames(surv.caps) <- inc.dates
surv.caps$Year = nests.scaled$NestYr-2021


getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}

f <- apply(surv.caps[,1:(length(inc.dates))],1,getFirst)
k <- apply(surv.caps[,1:(length(inc.dates))],1,getLast)
f;k
f<k

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
## Known Fate- Age, Incubation Date and Behavior

#' Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),                              # Intercept (1)
  nests.ready$LPDV,                                       # LPDV
  nests.ready$ParasiteDiversity,                          # Parasite Diversity
  nests.ready$LPDV * nests.ready$ParasiteDiversity        # Interaction
) 

#' Use ggpairs to visualize correlations
ggpairs(as.data.frame(X), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))

################################################################################
## Nimble Model

knownfate.ind <- nimbleCode({
  
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

################################################################################

#' Constants
Consts <- list(
  n.ind = nrow(nests.ready),   
  X = X,                       
  J = ncol(X),                 
  first = first,               
  last = last)

surv.caps.matrix <- as.matrix(surv.caps)
colnames(surv.caps.matrix) <- NULL

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
  code = knownfate.ind,          
  constants = Consts,          
  inits = Inits,              
  data = Data,                 
  nburnin = 3000,             
  niter = 20000,               
  thin = 1                     
)

summary(nimbleMCMC_samples)
colMeans(nimbleMCMC_samples[,1:4])
colSds(nimbleMCMC_samples[,1:4])
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:4]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix)

#' Create a vector of new names
new_names <- c("Intercept", 
               "LPDV",
               "Parasite Diversity",
               "LPDV * Parasite Diversity")

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
                      "Deciduous Forest",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Agriculture",
                      "Daily Minimum Temperature",
                      "Daily Precipitation") ~ "Landscape",
    TRUE ~ "Individual"
  ))
mean_estimates

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates3 <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c(
                                     "LPDV",
                                     "Parasite Diversity",
                                     "LPDV * Parasite Diversity",
                                     "Age * Parasite Diversity"
                                   )))
mean_estimates

################################################################################
## Beta Plot 

p2.betas <- ggplot(mean_estimates, 
                   aes(x = parameter, 
                       y = mean_estimate, 
                       color = Scale, 
                       shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Posterior Mean") +
  theme_minimal() + 
  scale_color_manual(values = c("Individual" = "#DA6509",
                                "Nest" = "#A44200",
                                "Landscape" = "#D65F5F")) +  
  scale_shape_manual(values = c("Individual" = 15, 
                                "Nest" = 17,
                                "Landscape" = 16)) +  
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.text.x = element_text(angle = 45, hjust = 1),   
    axis.text.y = element_text(angle = 0)   
  )

p2.betas


# Save multiple objects in a single RData file
save(p2.betas, mean_estimates3, file = "p2.betas.disease.PA.RData")


################################################################################
## Plot Interaction

################################################################################
## Plot Interaction

#' Check data
parasite_diversity_values <- unique(nests.ready$ParasiteDiversity)
beta_mean <- colMeans(samples_df)

#' Calculate lpdv0 and lpdv1 for each level of ParasiteDiversity
lpdv0 <- icloglog(beta_mean[1] + beta_mean[2]*0 + beta_mean[3]*parasite_diversity_values + beta_mean[4]*0*parasite_diversity_values)
lpdv1 <- icloglog(beta_mean[1] + beta_mean[2]*1 + beta_mean[3]*parasite_diversity_values + beta_mean[4]*1*parasite_diversity_values)

#' Dataframe with everything needed
data <- data.frame(
  ParasiteDiversity = rep(parasite_diversity_values, 2),
  lpdv = c(lpdv0, lpdv1),
  Infection = rep(c("No LPDV Infection", "LPDV Infection"), each = length(parasite_diversity_values))
)

#' Plot Interaction
ggplot(data, aes(x = ParasiteDiversity, y = lpdv, color = Infection, group = Infection)) +
  geom_line() +   
  geom_point() +   
  scale_color_manual(values = c("No LPDV Infection" = "black", "LPDV Infection" = "red")) + 
  xlab("Parasite Diversity") +
  ylab("Daily Nest Survival Probability")+
  labs(color = "Infection Status") +
  scale_x_continuous(breaks = c(1, 2)) +  
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(r = 8))) 

################################################################################

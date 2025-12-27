#---
# title: Daily Nest Survival Modeling of Wild Turkeys in Maryland
# authors: K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
#' **Purpose**: This script uses derived incubation start and end dates to fit a Bayesian known fate model 
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modeling daily nest survival
#' **Last Updated**: 12/27/25


################################################################################
## Load Packages and data

library(nimble)
library(MCMCvis)
library(tidyverse)
library(matrixStats)

load("Data Management/RData/Known Fate/Maryland/20250718_PreliminaryResults_Disease.RData")


################################################################################
## Data Prep- Encounter Histories

md.sample <- read_csv("Samples/Maryland/NestingSample_MD.updated.csv")
nests.scaled <- right_join(nests.scaled, md.sample) %>%
  drop_na()

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

getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}

f <- apply(surv.caps[,1:90],1,getFirst)
k <- apply(surv.caps[,1:90],1,getLast)
f;k
f<k


################################################################################
## Compile Parameters into matrix format

# Check to see if NAs exist in data
summary(nests.scaled)

# Create nests.ready object
nests.ready <- nests.scaled

# Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),               # Intercept (1)
nests.ready$LPDV                           # LPDV
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
      cloglog(phi[i,t]) <- inprod(beta[1:J], X[i, 1:J]) #+
        #weather.beta[3]*weather.array[i,t,1] + weather.beta[2]*weather.array[i,t,2] 
      
      #### Likelihood of Nest Survival (Bernoulli)
      z[i,t]~dbern(mu[i, t])
      
    }
  }
  
  #### Priors for beta coefficients
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
  #### Priors for weather coefficients
 # for (k in 1:n.weather.cov) {
  #  weather.beta[k] ~ dnorm(0, sd = 0.5)
 # }
  
})

################################################################################

# Constants
Consts <- list(
  n.ind = nrow(nests.ready),   
  X = X,                       
  J = ncol(X),                 
  first = first,               
  last = last  
 # n.weather.cov = dim(weather.array2)[3] +1,
  #weather.array = weather.array2
)
surv.caps.matrix <- as.matrix(surv.caps)
colnames(surv.caps.matrix) <- NULL

# Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  z = surv.caps.matrix               
)

# Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J)
#weather.beta = rep(0,Consts$n.weather.cov)
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
colMeans(nimbleMCMC_samples[,1:2])
colSds(nimbleMCMC_samples[,1:2])
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

# Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:2]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix)

# Create a vector of new names
new_names <- c("Intercept", 
               "LPDV")

# Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

# View the renamed data frame
head(samples_df)


################################################################################
## Output Trace Plots

# Force conversion to plain matrix
samples_matrix <- as.matrix(nimbleMCMC_samples[, 1:2])

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

# Calculate Bayesian credible intervals
credible_intervals <- samples_long %>%
  group_by(parameter) %>%
  summarise(
    lower = quantile(estimate, 0.05),   
    upper = quantile(estimate, 0.95),   
    .groups = 'drop'
  )

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
                      "Agriculture"
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
                                     "LPDV"
                                     )))
mean_estimates
mean_estimates1 <- mean_estimates %>%
  mutate(Scale = factor(Scale, levels = c("Individual", "Nest", "Landscape", "Climate")))


################################################################################
## Beta Plot 

# Beta estimates and associated 90% credible intervals 
p2.betas <- ggplot(mean_estimates1, 
                   aes(x = parameter, 
                       y = mean_estimate, 
                       color = Scale, 
                       shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = " Posterior Mean") +
  theme_minimal() + 
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
save(p2.betas,
     mean_estimates1, 
     samples_df,
     file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Known Fate/p2.betas.disease.MD.updated.RData",
     overwrite = T)


################################################################################
## Check Proportion of LPDV Samples less than 0

# Count how many values in 'Developed' are less than 0
count_less_than_0 <- sum(samples_df$LPDV < 0, na.rm = TRUE)

# Total number of non-NA values in LPDV
total_samples <- sum(!is.na(samples_df$LPDV))

# Calculate the proportion
prop_less_than_0 <- count_less_than_0 / total_samples

# Print the proportion
prop_less_than_0

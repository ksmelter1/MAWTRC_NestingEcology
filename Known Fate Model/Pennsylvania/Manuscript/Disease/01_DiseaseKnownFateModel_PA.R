# ---
# title: "Modelling the Effect of Disease and Parasites on Daily Nest Survival Probability in Pennsylvania"
# authors: K. Smelter
# output: html_document
# Last Updated: 12/27/25
# ---

################################################################################
## Load Packages

library(tidyverse)
library(terra)
library(mapview)
library(sf)
library(stringr)
library(nimble)
library(MCMCvis)
library(GGally)
library(matrixStats)
library(rmarkdown)

################################################################################
## Data Prep- Organize Nests

# Nest start and end date csv
nests <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250709_NestAttempts_allbirds_PA_Ready.csv")
nests

# Nest veg csv
nests.veg <- read_csv("Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250629_CleanedNestsVeg_2022_2023_2024.csv")
nests.veg

# Filter nests.veg to only include observations that have the same nestid as nests csv
nests.veg.filtered <- nests.veg %>%
  dplyr::filter(NestID%in% nests$NestID) %>%
  dplyr::filter(PlotType == "Nest")

# Merge the filtered nests.veg with nests by nestid
nests <- dplyr::left_join(nests, nests.veg.filtered, by = "NestID") %>%
  dplyr::select(-CheckDate.x) %>%
  dplyr::rename("NestFate" = NestFate.y)

# Rename and consolidate columns 
nests <- nests %>%
  dplyr::select(NestID, 
                BandID, 
                startI, 
                endI,
                NestFate,
                Lat,
                Long) 

# Switch coding to UTF-8
nests.scaled <- nests %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))


################################################################################
## Data Prep- Individual Covariates

# Read in captures csv
captures <- read_csv("Data Management/Csvs/Raw/Captures/captures_pa.csv")
captures
 
# Filter data to include only hens
# Select and rename columns
captures <- captures %>%
   dplyr::filter(sex == "F") %>%
   dplyr::select(bandid, sex, age, captyr) %>%
   dplyr::rename("BandID" = bandid)
  captures

# Merge columns 
 nests.scaled <- merge(captures, nests.scaled, by = "BandID") 

 # (?<=_): Only match is there is an underscore immediately before the number we are trying to extract
 # (\\d{4}): Match exactly 4 digits
 # (?=_) : Only match if the four digits are followed 
 # Create NestYr column 
nests.scaled <- nests.scaled %>%
  dplyr::mutate(NestYr = str_extract(NestID, "(?<=_)(\\d{4})(?=_)")) 

# Convert to numeric 
nests.scaled$NestYr <- as.numeric(nests.scaled$NestYr)
nests.scaled$captyr <- as.numeric(nests.scaled$captyr)

# Create a years since capture column
nests.scaled <- nests.scaled %>%
 dplyr::mutate(yrsincecap = NestYr-captyr)
 glimpse(nests.scaled)

# Assign Adult as the reference level
nests.scaled$age <- ifelse(nests.scaled$age == "J", 1, 
                             ifelse(nests.scaled$age == "A", 0, NA))

# Dealing with scaling age ad hoc
# If the bird is an adult and the years since capture is >1 assign it as an adult
# If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nests.scaled$age <- ifelse(nests.scaled$age == 1 & nests.scaled$yrsincecap >= 1, 0, nests.scaled$age)
table(nests.scaled$age)
glimpse(nests.scaled)


################################################################################
## Data Prep- Disease Covariates

# Read in raw virus csv
virus <- read_csv("Data Management/Csvs/Raw/Disease/LPDV_REV/Pennsylvania/Disease_Status.csv")
virus
 
# Rename columns 
virus <- virus %>%
dplyr::rename("BandID" = bandid)
virus

# # Specify the summed columns
columns_to_sum <- c("Capillaria", "Eimeria", "Ascarids")

# Create the 'ParasiteDiversity' column by summing the specified columns for each row
virus$ParasiteDiversity <- rowSums(virus[, columns_to_sum], na.rm = TRUE)

# Left join filtered viral data and nests together
# Remove all rows that contain an NA for parasite (No other rows had NAs)
nests.scaled <- left_join(nests.scaled, virus, by = "BandID") %>%
  drop_na()

# Scale and change variables to numeric
nests.scaled$ParasiteDiversity <- scale(nests.scaled$ParasiteDiversity)
nests.scaled$ParasiteDiversity <- as.numeric(nests.scaled$ParasiteDiversity)
nests.scaled$LPDV <- as.numeric(nests.scaled$LPDV)

# Assign 1 if the NestFate "Hatched", otherwise assign 0
nests.scaled$NestFate <- ifelse(nests.scaled$NestFate == "Hatched", 1, 0)


################################################################################
## Data Prep- Encounter Histories

# Sample used in PA analysis
pa.sample <- read_csv("Samples/Pennsylvania/PA_Sample_Updated.csv")

# Convert metrics to character to allow for joining the data
pa.sample$startI <- as.character(pa.sample$startI)
pa.sample$endI <- as.character(pa.sample$endI)

# Keep only nests.scaled observations that exist within our sample
nests.scaled <- right_join(nests.scaled, pa.sample) %>%
  drop_na()

# Convert start and end of incubation days to date objects
nests.scaled$startI <- as.Date(nests.scaled$startI)
nests.scaled$endI <- as.Date(nests.scaled$endI)

# Create a sequence of all unique dates between the min start and max end date
nests.scaled %>% group_by(NestYr) %>% summarise(mindate = min(startI),
                                                maxdate = max(endI))


# Create date ranges
# Julian dates will reflect these later on
inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(as.Date("2022-04-12"), as.Date("2022-09-16"), by = 1)
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

surv.caps = as.data.frame(surv.caps)
colnames(surv.caps) <- inc.dates
surv.caps$Year = nests.scaled$NestYr-2021

# getFirst and getLast functions from Kery and Schaub (2011)
# Get the first and last nesting events 
getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}

# Apply functions
# f should always be less than K
f <- apply(surv.caps[,1:(length(inc.dates))],1,getFirst)
k <- apply(surv.caps[,1:(length(inc.dates))],1,getLast)
f;k
f<k

# Calculate the number of nesting attempts for our reference level
num_matching_rows <- nests.scaled %>%
  filter(LPDV == 0, Capillaria == 0, Ascarids == 0, Eimeria == 0) %>%
  nrow()
print(num_matching_rows)


# Create nests.ready object
nests.ready <- nests.scaled
summary(nests.ready)

# Compile parameters in matrix
# Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),                        # Intercept (1)
  nests.ready$LPDV,                                 # LPDV
  nests.ready$Capillaria,                           # Capillaria spp.
  nests.ready$Eimeria,                              # Eimeria spp.
  nests.ready$Ascarids,                             # Ascarids
  nests.ready$LPDV * nests.ready$Capillaria,        # LPDV * Capillaria
  nests.ready$LPDV * nests.ready$Eimeria,           # LPDV * Eimeria
  nests.ready$LPDV * nests.ready$Ascarids           # LPDV * Ascarids
) 

# Use ggpairs to visualize correlations
ggpairs(as.data.frame(X), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))


################################################################################
## Nimble Model 

knownfate.disease <- nimbleCode({
  
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
    beta[j] ~ dnorm(0, sd = 1)  
  }
  
})

################################################################################

# Constants
Consts <- list(
  n.ind = nrow(nests.ready),   
  X = X,                       
  J = ncol(X),                 
  first = first,               
  last = last)

# Data (daily nest survival: 1 = survived, 0 = failed)
Data <- list(
  z = surv.caps            
)

# Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J)
)

start <- Sys.time()

# Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = knownfate.disease,          
  constants = Consts,          
  inits = Inits,              
  data = Data,                 
  nburnin = 3000,             
  niter = 20000,               
  thin = 1                     
)

# Check outputs and posterior chains
summary(nimbleMCMC_samples)
colMeans(nimbleMCMC_samples[,1:8])
colSds(nimbleMCMC_samples[,1:8])
MCMCtrace(nimbleMCMC_samples[,1:8],
          iter = 20000,
          pdf = FALSE)

# Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:8]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix)

# Create a vector of new names
new_names <- c("Intercept", 
               "LPDV",
               "Capillaria sp.",
               "Eimeria sp.",
               "Ascarids",
               "LPDV * Capillaria sp.",
               "LPDV * Eimeria sp.",
               "LPDV * Ascarids")

# Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

# View the renamed data frame
head(samples_df)

################################################################################
## Output Trace Plots

# Force conversion to plain matrix
samples_matrix <- as.matrix(nimbleMCMC_samples[, 1:8])

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
# Remove intercept
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

# Create scale legend for plot
# Individual Covariates represent disease outputs
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

# Organize variables into levels to be displayed
# Filter out the intercept
mean_estimates3 <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c(
                                     "LPDV",
                                     "Capillaria sp.",
                                     "Eimeria sp.",
                                     "Ascarids",
                                     "LPDV * Capillaria sp.",
                                     "LPDV * Eimeria sp.",
                                     "LPDV * Ascarids"
                                   )))
mean_estimates


################################################################################
## Beta Plots

# Create plot that visualizes the posterior mean and credible intervals
# If the estimate and credible intervals overlap zero there was no relationship detected
# If the estimate is positive and the credible intervals don't overlap zero = positive effect
# If the estimate is negative and the credible intervals don't overlap zero = negative effect
# Coordinate Flip
p2.betas <- ggplot(mean_estimates3, 
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
  coord_flip() +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.text.x = element_text(angle = 0, hjust = 0),   
    axis.text.y = element_text(angle = 0),
    legend.position = "none"
  )

p2.betas

# Standard format 
p3.betas <- ggplot(mean_estimates3, 
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
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.text.x = element_text(angle = 45, hjust = 1),   
    axis.text.y = element_text(angle = 0),  
    legend.position = "none"
  )

p3.betas

# Save multiple objects in a single RData file
save(p2.betas, mean_estimates3, samples_df,  file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Known Fate/p2.betas.disease.PA.updated.RData", overwrite = T)

################################################################################
###############################################################################X


#'---
#' title: Nest-site selection of female wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script creates a Bayesian conditional logistic regression model for nest-site selection in JAGs using the gathered covariates 

#####################
## Load Packages ##
####################

#' Vector of package names
packages <- c("R2jags",
              "MCMCvis",
              "dplyr",
              "tidyverse")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)

########################
## Data Preparation ##
########################

#' Load in RData
load("Data Management/RData/Nest-Site Selection/Covs/02_Covs.RData")

#' Read in captures csv
captures <- read_csv("Data Management/Csvs/raw data/captures.csv") 
captures

#' Filter data to include only hens 
captures <- captures %>%
  dplyr::filter(sex == "F") %>%
  dplyr::select(bandid, sex, age, weight, captyr)
captures

#' Create nest.data object for models
#' Scale continuous predictors
nest.data <- pa.nests.covs 

#' Merge data
nest.data <- merge(captures, nest.data, by = "bandid")

#' Add in columns
nest.data$basal <- pa.nests$stemcnt
nest.data$grssfrb <- pa.nests$grssfrb

#' (?<=_): Only match is there is an underscore immediately before the number we are trying to extract
#' (\\d{4}): Match exactly 4 digits
#' (?=_) : Only match if the four digits are followed 
nest.data <- nest.data %>%
  dplyr::mutate(nestyr = str_extract(nest_id, "(?<=_)(\\d{4})(?=_)"))
glimpse(nest.data)

nest.data$nestyr <- as.numeric(nest.data$nestyr)
nest.data$captyr <- as.numeric(nest.data$captyr)

#' Create a years since capture column
nest.data <- nest.data %>%
  dplyr::mutate(yrsincecap = nestyr-captyr)
glimpse(nest.data)

#' Select columns of interest
  nest.data <- nest.data %>%
    dplyr::select(nest_id, bandid, avrgmxv, grssfrb, percwdy,
                case, Developed, Deciduous, Mixed, Evergreen, Agriculture, elev,
                primary, secondary, basal, weight, age, captyr, nestyr, yrsincecap)
  
  #' Scale predictors
  nest.data <- nest.data %>%
  dplyr::mutate(basal = scale(basal)) %>%
  dplyr::mutate(percwdy = scale(percwdy)) %>%
  dplyr::mutate(grssfrb = scale(grssfrb)) %>%
  dplyr::mutate(avgvisob = scale(avrgmxv)) %>%
  dplyr::mutate(elev = scale(elev)) %>%
  dplyr::mutate(primary = scale(primary)) %>%
  dplyr::mutate(secondary = scale(secondary)) %>%
  dplyr::mutate(weight = scale(weight)) 
  
  str(nest.data)
  glimpse(nest.data)
  
  #' Change structure from a matrix to numeric (Scaled covariates)
  nest.data <- nest.data %>%
  dplyr::mutate(percwdy = as.numeric(percwdy)) %>%
  dplyr::mutate(grssfrb = as.numeric(grssfrb)) %>%
  dplyr::mutate(avgvisob = as.numeric(avgvisob)) %>%
  dplyr::mutate(elev = as.numeric(elev)) %>%
  dplyr::mutate(primary = as.numeric(primary)) %>%
  dplyr::mutate(secondary = as.numeric(secondary)) %>%
  dplyr::mutate(basal = as.numeric(basal)) %>%
  dplyr::mutate(weight = as.numeric(weight)) %>%
  dplyr::mutate(Deciduous = as.numeric(Deciduous)) %>%
  dplyr::mutate(Mixed = as.numeric(Mixed)) %>%
  dplyr::mutate(Evergreen= as.numeric(Evergreen)) %>%
  dplyr::mutate(Agriculture = as.numeric(Agriculture)) %>%
  dplyr::mutate(Developed = as.numeric(Developed))

glimpse(nest.data)
str(nest.data)

#' Change case to numeric
nest.data$case <- as.numeric(nest.data$case)

#' Order df by nestid
nest.data <- nest.data[order(nest.data$nest_id),]

#' Change Nest_ID_V to numeric 
#' First must remove underlines
nest.data$nest_id <- gsub("_", "", nest.data$nest_id)
nest.data$nest_id<- as.numeric(nest.data$nest_id)

#' Create str_id column 
#' This allows the loop to iterate through the steps associated with each bird 
#' cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  group_by(nest_id) %>%
  mutate(str_ID=cur_group_id())

#' Check
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

#' Assign juvenile as the reference level
nest.data.ready$age <- ifelse(nest.data.ready$age == "A", 1, 
                              ifelse(nest.data.ready$age == "J", 0, NA))


#' Dealing with scaling age ad hoc
#' If the bird is an adult and the years since capture is >1 assign it as an adult
#' If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nest.data.ready$age <- ifelse(nest.data.ready$age == "A" & nest.data.ready$yrsincecap >= 1, "A", nest.data.ready$age)
glimpse(nest.data.ready)

#' Check
glimpse(nest.data.ready)
summary(nest.data.ready)
str(nest.data.ready$age)

######################################################
## JAGS Models

####################
## No Developed ##
####################
#' We keep developed out of the model because it is correlated with agriculture

nest.selection.subset.mod <- "  
 model{
  for (i in 1:I){
   case[i]~dpois(lambda[i])
   log(lambda[i])<-primary[i]*beta.primary +
   secondary[i]*beta.secondary +
   deciduous[i]*beta.deciduous_forest + 
   mixed[i]*beta.mixed_forest + 
   evergreen[i]*beta.evergreen_forest + 
   agriculture[i]*beta.agriculture + 
   elevation[i] *beta.elevation + 
   visob[i]* beta.visob +
   grassforb[i]*beta.grassforb +
   percwoody[i]*beta.percwoody + 
   basalarea[i]*beta.basalarea +
   age[i] * beta.age +
   alpha[str_ID[i]] 
  }
  
  #Priors
  beta.deciduous_forest~dnorm(0,0.0001)  #Deciduous Forest
  beta.mixed_forest~dnorm(0,0.0001)      #Mixed Forest
  beta.evergreen_forest~dnorm(0,0.0001)  #Mixed Forest
  beta.agriculture~dnorm(0,0.0001)       #Agriculture
  beta.primary~dnorm(0,0.0001)           #Primary
  beta.secondary~dnorm(0,0.0001)         #Secondary
  beta.elevation~dnorm(0,0.0001)         #Elevation
  beta.visob~dnorm(0,0.0001)             #Average Visual Obstruction
  beta.grassforb~dnorm(0,0.0001)         #Percent Grass Forb
  beta.percwoody~dnorm(0,0.0001)         #Percent Woody
  beta.age~dnorm(0,0.0001)             #Age (Adult) Juvenile is the reference level
  beta.basalarea~dnorm(0,0.0001)         #Basal Area
  
  for (k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
 }
"
#' Data list for JAGS
jags_data <- list(
  case = nest.data.ready$case,
  primary = nest.data.ready$primary,
  secondary = nest.data.ready$secondary,
  mixed = nest.data.ready$Mixed,
  evergreen = nest.data.ready$Evergreen,
  deciduous = nest.data$Deciduous,
  agriculture = nest.data.ready$Agriculture,
  elevation = nest.data.ready$elev,
  visob = nest.data.ready$avgvisob,
  grassforb = nest.data.ready$grssfrb,
  percwoody = nest.data.ready$percwdy,
  age = nest.data.ready$age,
  basalarea = nest.data.ready$basal,
  str_ID = c(nest.data.ready$str_ID),
  I = nrow(nest.data.ready),
  K =length(unique(nest.data.ready$str_ID))
)

#' Initialize JAGS
jags_model <- jags.model(textConnection(nest.selection.subset.mod), 
                         data = jags_data,
                         n.chains = 3)

#' Burn-in and sampling from posterior distribution
jags_samples.nss.com1 <- coda.samples(jags_model, c("beta.deciduous_forest",
                                           "beta.mixed_forest",
                                           "beta.evergreen_forest",
                                           "beta.agriculture",
                                           "beta.primary",
                                           "beta.secondary",
                                           "beta.elevation",
                                           "beta.visob",
                                           "beta.grassforb",
                                           "beta.percwoody",
                                           "beta.age",
                                           "beta.basalarea"), 
                             n.iter = 200000, 
                             n.burnin = 30000,
                             n.thin = 10)
#' View traceplots
MCMCtrace(jags_samples.nss.com1, pdf=F)

#' Read in RDS file of results
jags_samples.nss.com1 <- readRDS("Data Management/RData/Nest-Site Selection/ModelResults/NSS_Samples.RDS")

#' Convert mcmc.list to matrix
samples_matrix <- as.matrix(jags_samples.nss.com1)

#' Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' View the first few rows of the samples data frame
head(samples_df)

#' Rename columns
samples_df.names <- samples_df %>%
  dplyr::rename("Agriculture" = beta.agriculture) %>%
  dplyr::rename("Deciduous Forest" = beta.deciduous_forest) %>%
  dplyr::rename("Mixed Forest" = beta.mixed_forest) %>%
  dplyr::rename("Evergreen Forest" = beta.evergreen_forest) %>%
  dplyr::rename("Elevation" = beta.elevation) %>%
  dplyr::rename("Distance to Primary Road" = beta.primary) %>%
  dplyr::rename("Distance to Secondary Road" = beta.secondary) %>%
  dplyr::rename("Percent Grass Forb" = beta.grassforb) %>%
  dplyr::rename("Visual Obstruction" = beta.visob) %>%
  dplyr::rename("Percent Woody" = beta.percwoody) %>%
  dplyr::rename("Basal Area" = beta.basalarea) %>%
  dplyr::rename("Percent Fern" = beta.percfern) 
  
  
#' View the first few rows of the samples data frame
head(samples_df.names)


#' Reshape the data into long format for ggplot
samples_long <- samples_df.names %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

#' View the reshaped data
head(samples_long)

#' Calculate Bayesian credible intervals
credible_intervals <- samples_long %>%
  group_by(parameter) %>%
  summarise(
    lower = quantile(estimate, 0.025),   # 2.5th percentile
    upper = quantile(estimate, 0.975),   # 97.5th percentile
    .groups = 'drop'
  )

#' Create plot
p1.density <-ggplot(samples_long, aes(x = estimate, fill = parameter)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  facet_wrap(~parameter, scales = "free", ncol = 2) +  # Create separate panels for each parameter
  labs(x = "Parameter Estimate", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none") +
  #' Add vertical lines for the 95% credible intervals
  geom_vline(data = credible_intervals, aes(xintercept = lower), linetype = "dashed", color = "red") +
  geom_vline(data = credible_intervals, aes(xintercept = upper), linetype = "dashed", color = "red")

p1.density

#' Calculate mean estimates for each parameter
mean_estimates <- samples_long %>%
  dplyr::group_by(parameter) %>%
  summarise(
    mean_estimate = mean(estimate),  # mean of the estimates
    lower = quantile(estimate, 0.025),  # 2.5th percentile (lower bound of credible interval)
    upper = quantile(estimate, 0.975),  # 97.5th percentile (upper bound of credible interval)
    .groups = 'drop'
  )
#' Take the mean value for each estimate 
mean_estimates <- mean_estimates %>%
     Scale = case_when(
      parameter %in% c("Percent Grass Forb", "Percent Woody", "Visual Obstruction", "Percent Fern", "Basal Area") ~ "Microscale",
      parameter == "Adult", "Juvenile" ~ "Individual",  
      TRUE ~ "Macroscale")
    


#' Organize variables into levels to be displayed
mean_estimates <- mean_estimates %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Visual Obstruction",
                                       "Percent Woody", 
                                       "Percent Grass Forb",
                                       "Percent Fern",
                                       "Basal Area",
                                       "Mixed Forest",
                                       "Evergreen Forest",
                                       "Deciduous Forest",
                                       "Elevation",
                                       "Distance to Primary Road",
                                       "Distance to Secondary Road",
                                       "Agriculture",
                                       "Juvenile",
                                       "Adult")))  

#' Create the plot with mean estimates, credible intervals, color by scale, shape by scale, and faceting by scale
p2.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5) +  # Points for the mean estimate 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  # Error bars for credible intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at 0
  labs(x = "Parameter", y = "Beta Estimate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("Microscale" = "#A44200", "Macroscale" = "#D65F5F", "Individual" = "#5f0f40") +  # Set color for Microscale, Macroscale, and individual predictors
  scale_shape_manual(values = c("Microscale" = 17, "Macroscale" = 16, "Individual" = 15))) 
p2.betas

p3.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5) +  # Points for the mean estimate 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  # Error bars for credible intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +  # Horizontal line at 0
  labs(x = "Parameter", y = "Beta Estimate") +
  theme_minimal() +
  scale_color_manual(values = c("Microscale" = "#A44200", "Macroscale" = "#D65F5F", "Individual" = "#5f0f40" )) +  # Set color for Microscale and Macroscale
  scale_shape_manual(values = c("Microscale" = 17, "Macroscale" = 16, "")) +
  coord_flip() +
  theme(
    axis.title.x = element_text(margin = margin(t = 10))  # Add padding to the x-axis title (top margin)
  )
p3.betas


##########################
## Prediction Plots ##
#########################

#' Mixed Forest Prediction Plot
Mixed <- ggplot(nest.data.ready, aes(x=Mixed, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Mixed Forest", y="Probability of Use") +
  expand_limits(x=2) +
  theme_minimal() +
  ylim(0,1)
Mixed

#' Deciduous Forest Prediction Plot
Deciduous <- ggplot(nest.data.ready, aes(x= Deciduous, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Deciduous Forest", y="Probability of Use")+
  expand_limits(x=2) +
  theme_minimal() + 
  ylim(0,1)
Deciduous

#' Evergreen Forest Prediction Plot
Evergreen <- ggplot(nest.data.ready, aes(x= Evergreen, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Evergreen Forest", y="Probability of Use")+
  expand_limits(x=2) +
  ylim(0,1) +
  theme_minimal()
Evergreen

#' Developed Prediction Plot
Developed <- ggplot(nest.data.ready, aes(x= Developed, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Developed", y="Probability of Use")+
  expand_limits(x=1) +
  ylim(0,1) +
  theme_minimal()
Developed

#' Agriculture Prediction Plot
Agriculture <- ggplot(nest.data.ready, aes(x= Agriculture, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Agriculture", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=2) +
  theme_minimal()
Agriculture


#' Visual Obstruction Prediction Plot
Visob<- ggplot(nest.data.ready, aes(x= avgvisob, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x=" Visual Obstruction", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=2) +
  theme_minimal()
Visob

#' Percent Woody Prediction Plot
woody<- ggplot(nest.data, aes(x= percwdy, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x=" Percent Woody", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=2) +
  theme_minimal()
woody

#' Percent Grass Forb Prediction Plot
GrassForb<- ggplot(nest.data.ready, aes(x= grssfrb, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x=" Percent Grass/Forb", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=2) +
  theme_minimal()
GrassForb

#' Distance to Primary Road Prediction Plot
Primary<- ggplot(nest.data, aes(x= primary, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Distance to Primary Road", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=1) +
  theme_minimal()
Primary

#' Distance Secondary Road
Secondary<- ggplot(nest.data.ready, aes(x= secondary, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x=" Distance to Secondary Road", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=1) +
  theme_minimal()
Secondary

#' Distance Secondary Road
basal<- ggplot(nest.data.ready, aes(x= basal, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Basal Area", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=2) +
  theme_minimal()
basal

#' Distance Secondary Road
elev<- ggplot(nest.data.ready, aes(x= elev, y= case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  labs(x="Elevation", y="Probability of Use")+
  ylim(0,1) +
  expand_limits(x=2) +
  theme_minimal()
elev


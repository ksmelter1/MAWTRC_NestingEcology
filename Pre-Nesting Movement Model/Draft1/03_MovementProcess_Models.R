
#'---
#' title: Habitat selection of female wild turkeys during incubation (an SSF analysis)
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script fits SSF models in frequentist and Bayesian frameworks and extracts coefficient data


####################
## Load Packages ##
####################

#' Vector of package names
packages <- c("R2jags",
              "MCMCvis",
              "dplyr",
              "survival",
              "lme4",
              "rstanarm",
              "tidyr",
              "jagsUI")


#' Function to load a package or install it if not already installed
#' Found this function 
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


#' Visually checked and the projections look good
load("Data Management/RData/Individual-Specific Movement Process/RData Files/Covs_Ready.RData") 

  
dat_2.ready <- random_steps%>% dplyr::mutate(elev = scale(elev)) %>%
                 dplyr::mutate(secondary = scale(secondary)) %>%
                 dplyr::mutate(primary = scale(primary)) 


glimpse(dat_2.ready)

dat_2.ready <- dat_2.ready %>% dplyr::mutate(elev = as.numeric(elev)) %>%
                 dplyr::mutate(secondary = as.numeric(secondary)) %>%
                 dplyr::mutate(primary = as.numeric(primary)) 

glimpse(dat_2.ready)

#' Plot out data based on covariates
# plot(final.data$Shrub)
# plot(final.data$Agriculture)
# plot(final.data$Forest)
# plot(final.data$Developed)

####################
## Organize Data ##
####################

#' Check
str(dat_2.ready)

#' Change ID to numeric
dat_2.ready$id <- gsub("_", "", dat_2.ready$id)
dat_2.ready$id<- as.numeric(dat_2.ready$id)

#' Create strata column by merging BirdID and step_id
dat_2.ready$NA_ID <- paste(dat_2.ready$id, dat_2.ready$step_id_, sep = "_")

#' Change NA_ID to integer
dat_2.ready$NA_ID <- gsub("_", "", dat_2.ready$NA_ID)
dat_2.ready$NA_ID<- as.numeric(dat_2.ready$NA_ID)

class(dat_2.ready)
str(dat_2.ready)

#' Add numerical variable for animals:
dat_2.ready$ANIMAL_ID <- as.numeric(as.factor(dat_2.ready$id))

#' Stratum ID is given as "NA_ID" in the data; 
#' It is easier to have sequential enumeration, so let's generate a new stratum-ID variable str_ID:
d.map <- data.frame(NA_ID=unique(dat_2.ready$NA_ID),str_ID=1:length(unique(dat_2.ready$NA_ID)))
dat_2.ready$str_ID <- d.map[match(dat_2.ready$NA_ID,d.map$NA_ID),"str_ID"]
dat_2.ready <- dat_2.ready[order(dat_2.ready$str_ID),] 
glimpse(dat_2.ready)

#' Function to fill NA values with 0 throughout dataframe
#' Not sure if this will be the approach we take but I'm going with this for now
fill_NA_with_value <- function(df, value = 0) {
  #' Ensure that all columns are of numeric type to handle the replacement
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      #' Replace NA values with the specified value
      replace(col, is.na(col), value)
    } else {
      #' Return the column unchanged if it's not numeric
      col
    }
  })
  return(df)
}

#' Apply the function to fill NA values with 0.0001
dat_2.ready <- fill_NA_with_value(dat_2.ready)

#' Check
which(dat_2.ready$step_id_==3)
table(dat_2.ready$case_)
class(dat_2.ready)
summary(dat_2.ready)


####################
## No Developed ##
####################

#' Specify JAGS model structure
#' Removed Developed term due to its correlation with Agriculture

Ind.movement.nodev.mod <- "  
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
   alpha[str_ID[i]] 
  }
  
  #Priors
  beta.deciduous_forest~dnorm(0,0.0001) #Deciduous Forest
  beta.mixed_forest~dnorm(0,0.0001) #Mixed Forest
  beta.evergreen_forest~dnorm(0,0.0001) #Mixed Forest
  beta.agriculture~dnorm(0,0.0001) #Agriculture
  beta.primary~dnorm(0,0.0001) #Primary
  beta.secondary~dnorm(0,0.0001) #Secondary
  beta.elevation~dnorm(0,0.0001) #Elevation
  for (k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
 }
"
#' Data list for JAGS
jags_data <- list(
  case = dat_2.ready$case_,
  primary = dat_2.ready$primary,
  secondary = dat_2.ready$secondary,
  mixed = dat_2.ready$Mixed,
  evergreen = dat_2.ready$Evergreen,
  deciduous = dat_2.ready$Deciduous,
  agriculture = dat_2.ready$Agriculture,
  elevation = dat_2.ready$elev,
  str_ID = c(dat_2.ready$str_ID),
  I = nrow(dat_2.ready),
  K =length(unique(dat_2.ready$str_ID))
)

#' Initialize JAGS
jags_model <- jags.model(textConnection(Ind.movement.nodev.mod), 
                         data = jags_data,
                         n.chains = 3)

#' Burn-in and sampling from posterior distribution
Ind.movement.nodev.mod <- coda.samples(jags_model, c("beta.deciduous_forest",
                                                    "beta.mixed_forest",
                                                    "beta.evergreen_forest",
                                                    "beta.agriculture",
                                                    "beta.primary",
                                                    "beta.secondary",
                                                    "beta.elevation"), 
                                      n.iter = 200000,  n.burnin = 30000)
MCMCtrace(jags_samples, pdf=F)

#' Convert mcmc.list to matrix
samples_matrix <- as.matrix(Ind.movement.nodev.mod)

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
  dplyr::rename("Distance to Secondary Road" = beta.secondary) 

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
ggplot(samples_long, aes(x = estimate, fill = parameter)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  facet_wrap(~parameter, scales = "free", ncol = 2) +  # Create separate panels for each parameter
  labs(x = "Parameter Estimate", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none") +
  #' Add vertical lines for the 95% credible intervals
  geom_vline(data = credible_intervals, aes(xintercept = lower), linetype = "dashed", color = "red") +
  geom_vline(data = credible_intervals, aes(xintercept = upper), linetype = "dashed", color = "red")

################################################################################
################################################################################

#######################################
## Conditional Logistic Regression ##
#######################################

#' amt code
#' fit_clogit is a wrapper for the survival package

#' m1 <- final.data %>% fit_clogit(case_ ~ Developed + strata(step_id_) +strata(BirdID))
#' m2 <- final.data %>% fit_clogit(case_ ~ Agriculture + strata(step_id_) + strata(BirdID))
#' m3 <- final.data %>% fit_clogit(case_ ~ Forest + strata(step_id_) + strata(BirdID))
#' m4 <- final.data %>% fit_clogit(case_ ~ Shrub + strata(step_id_) + strata(BirdID))
#' m5<- final.data %>% fit_clogit(case_ ~ distfromroad + strata(step_id_) + strata(BirdID))
#' m6 <- final.data %>% fit_clogit(case_ ~ elev + strata(step_id_) + strata(BirdID))
#' 
#' summary(m1)
#' summary(m2)
#' summary(m3)
#' summary(m4)
#' summary(m5)
#' summary(m6)
#' 
#' ###############################################
#' ## Generalized Linear Mixed Effects Models ##
#' ###############################################
#' 
#' #' Likelihood framework using glmer package
#' m7 <- lmer(case_ ~ Developed + (1|step_id_) +(Developed|BirdID),data = final.data)
#' m8 <- lmer(case_ ~ Agriculture + (1|step_id_) +(Agriculture|BirdID),data = final.data)
#' m9 <- lmer(case_ ~ Forest + (1|step_id_) +(Forest|BirdID),data = final.data)
#' m10 <- lmer(case_ ~ Shrub+ (1|step_id_) +(Shrub|BirdID),data = final.data)
#' m11 <- lmer(case_ ~ elev + (1|step_id_) +(elev|BirdID),data = final.data)
#' m12 <- lmer(case_ ~ distfromroad + (1|step_id_) +(distfromroad|BirdID),data = final.data)
#' 
#' summary(m7)
#' summary(m8)
#' summary(m9)
#' summary(m10)
#' summary(m11)
#' summary(m12)
#' 
#' 
#' ########################################################
#' ## Bayesian Generalized Linear Mixed Effects Models ##
#' ########################################################
#' 
#' #' Bayesian glmm using rstanarm
#' #' Doesn't run
#' m13 <- rstanarm::stan_lmer(case_ ~ Developed + (1|step_id_) +(Developed|BirdID),data = final.data)
#' m14 <- rstanarm::stan_lmer(case_ ~ Agriculture + (1|step_id_) +(Agriculture|BirdID),data = final.data)
#' m15 <- rstanarm::stan_lmer(case_ ~ Forest + (1|step_id_) +(Forest|BirdID),data = final.data)
#' m16 <- rstanarm::stan_lmer(case_ ~ Shrub+ (1|step_id_) +(Shrub|BirdID),data = final.data)
#' m17 <- rstanarm::stan_lmer(case_ ~ elev + (1|step_id_) +(elev|BirdID),data = final.data)
#' m18 <- rstanarm::stan_lmer(case_ ~ distfromroad + (1|step_id_) +(distfromroad|BirdID),data = final.data)
#' 
#' summary(m13)
#' summary(m14)
#' summary(m15)
#' summary(m16)
#' summary(m17)
#' summary(m18)
#' 
#' 
#' ############################
#' ## Extract Coefficients ##
#' ############################
#' 
#' Dist2Road <- broom::tidy(m5$model) %>% 
#' filter( term== "distfromroad") %>%
#' bind_cols(confint(m5$model)) %>%
#' rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`) 
#' 
#' elev <- broom::tidy(m6$model) %>% 
#' filter(term == "elev") %>%
#' bind_cols(confint(m6$model)) %>%
#'  rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)
#' 
#' Developed <- broom::tidy(m1$model) %>% 
#'   filter(term== "Developed") %>%
#'   bind_cols(confint(m1$model)) %>%
#'   rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)
#' 
#' Forest <- broom::tidy(m3$model) %>% 
#'   filter(term == "Forest") %>%
#'   bind_cols(confint(m3$model)) %>%
#'   rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)
#' 
#' Shrub <- broom::tidy(m4$model) %>% 
#'   filter(term == "Shrub") %>%
#'   bind_cols(confint(m4$model)) %>%
#'   rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)
#' 
#' Agriculture <- broom::tidy(m2$model) %>% 
#'   filter(term == "Agriculture") %>%
#'   bind_cols(confint(m2$model)) %>%
#'   rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)
#' 
#' #' Create df for visualizations
#' ggdf <- bind_rows(Developed,Forest, Shrub,
#'                   Agriculture, Dist2Road, elev) 
#' 
#' #' Changing data to factor in order to rename names from rows using levesls
#' ggdf$term<-as.factor(ggdf$term)
#' levels(ggdf$term)[levels(ggdf$term) == 'distfromroad'] <- 'Distance from Road'
#' levels(ggdf$term)[levels(ggdf$term) == 'elev'] <- 'Elevation'
#' levels(ggdf$term)[levels(ggdf$term) == 'Agriculture'] <- 'Agriculture'
#' levels(ggdf$term)[levels(ggdf$term) == 'Developed'] <- 'Developed'
#' levels(ggdf$term)[levels(ggdf$term) == 'Forest'] <- 'Forest'
#' levels(ggdf$term)[levels(ggdf$term) == 'Shrub'] <- 'Shrub/Scrub'
#' 
#' 
#' #' Save RDS file of ggdf for plots
#' #saveRDS(ggdf, "Working_SSF_Outputs.RDS")
#' 
#' 

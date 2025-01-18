
#'---
#' title: Habitat selection of wild turkeys during incubation (an SSF analysis)
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
#' **Purpose**: This script runs SSF models in frequentist and Bayesian frameworks and extracts coefficient data


####################
## Load Packages ##
####################

require(amt)
require(dplyr)
require(survival)
require(lme4)
require(rstanarm)
require(R2jags)


#' Read in RDS file with covs 
final.data <- readRDS("Data Management/RData/Individual-Specific Movement Process/Working_SSF_Ready.RDS") %>%
  dplyr::select(BirdID, case_, step_id_, tod_end_,geometry, 
                Agriculture, 
                Developed, 
                Forest,
                elev,
                Shrub,
                distfromroad,
                t1_) 

str(final.data)


#' Plot out data based on covariates
# plot(final.data$Shrub)
# plot(final.data$Agriculture)
# plot(final.data$Forest)
# plot(final.data$Developed)

####################
## Organize Data ##
####################

#' Create strata column by merging BirdID and step_id
final.data$NA_ID <- paste(final.data$BirdID, final.data$step_id_, sep = "_")

#' Change NA_ID to integer
#First get rid of underscores
# final.data$NA_ID <- gsub("_", "", final.data$NA_ID)
# final.data$NA_ID<- as.integer(final.data$NA_ID)

###############################
## Step Selection Analysis ##
###############################

#' Fit SSFs using conditional logistic regression in the survival package 
#' Used nests=1, Stratified randomly sampled available nests=0
#' Landscape covariates, and distance from road

############
## JAGS ##
############
## Franny's Bayesian Conditional Logistic Regression
#' This is a test model to try to get these bayesian clogit models running properly
#' Model initializes
test.mod <- "  
 model{
  for (i in 1:I){
   case[i]~dpois(lambda[i])
   log(lambda[i])<-distfromroad[i] +beta.distfromroad+developed[i]*beta.developed+alpha[str_ID[i]] #Dist from road is reference level
  }
  
  #Priors
  beta.developed~dnorm(0,0.0001) #Developed
  beta.distfromroad~dnorm(0,0.0001) #Distance from Road

  for (k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
 }
"
#' Data list for JAGS
jags_data <- list(
  case = final.data$case_,
  distfromroad = final.data$distfromroad,
  developed= final.data$Developed,
  str_ID = c(final.data$NA_ID),
  I=nrow(final.data),
  K=length(unique(final.data$NA_ID))
)

#' Initialize JAGS
#' Won't run
jags_model <- jags.model(textConnection(test.mod), data = jags_data, n.chains = 1)

# Burn-in and sampling
jags_samples <- coda.samples(jags_model, c("beta.developed","beta.distfromroad"), n.iter = 10000,  n.burnin = 2000)


#######################################
## Conditional Logistic Regression ##
#######################################

#' amt code
#' fit_clogit is a wrapper for the survival package

m1 <- final.data %>% fit_clogit(case_ ~ Developed + strata(step_id_) +strata(BirdID))
m2 <- final.data %>% fit_clogit(case_ ~ Agriculture + strata(step_id_) + strata(BirdID))
m3 <- final.data %>% fit_clogit(case_ ~ Forest + strata(step_id_) + strata(BirdID))
m4 <- final.data %>% fit_clogit(case_ ~ Shrub + strata(step_id_) + strata(BirdID))
m5<- final.data %>% fit_clogit(case_ ~ distfromroad + strata(step_id_) + strata(BirdID))
m6 <- final.data %>% fit_clogit(case_ ~ elev + strata(step_id_) + strata(BirdID))

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

###############################################
## Generalized Linear Mixed Effects Models ##
###############################################

#' Likelihood framework using glmer package
m7 <- lmer(case_ ~ Developed + (1|step_id_) +(Developed|BirdID),data = final.data)
m8 <- lmer(case_ ~ Agriculture + (1|step_id_) +(Agriculture|BirdID),data = final.data)
m9 <- lmer(case_ ~ Forest + (1|step_id_) +(Forest|BirdID),data = final.data)
m10 <- lmer(case_ ~ Shrub+ (1|step_id_) +(Shrub|BirdID),data = final.data)
m11 <- lmer(case_ ~ elev + (1|step_id_) +(elev|BirdID),data = final.data)
m12 <- lmer(case_ ~ distfromroad + (1|step_id_) +(distfromroad|BirdID),data = final.data)

summary(m7)
summary(m8)
summary(m9)
summary(m10)
summary(m11)
summary(m12)


########################################################
## Bayesian Generalized Linear Mixed Effects Models ##
########################################################

#' Bayesian glmm using rstanarm
#' Doesn't run
m13 <- rstanarm::stan_lmer(case_ ~ Developed + (1|step_id_) +(Developed|BirdID),data = final.data)
m14 <- rstanarm::stan_lmer(case_ ~ Agriculture + (1|step_id_) +(Agriculture|BirdID),data = final.data)
m15 <- rstanarm::stan_lmer(case_ ~ Forest + (1|step_id_) +(Forest|BirdID),data = final.data)
m16 <- rstanarm::stan_lmer(case_ ~ Shrub+ (1|step_id_) +(Shrub|BirdID),data = final.data)
m17 <- rstanarm::stan_lmer(case_ ~ elev + (1|step_id_) +(elev|BirdID),data = final.data)
m18 <- rstanarm::stan_lmer(case_ ~ distfromroad + (1|step_id_) +(distfromroad|BirdID),data = final.data)

summary(m13)
summary(m14)
summary(m15)
summary(m16)
summary(m17)
summary(m18)


############################
## Extract Coefficients ##
############################

Dist2Road <- broom::tidy(m5$model) %>% 
filter( term== "distfromroad") %>%
bind_cols(confint(m5$model)) %>%
rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`) 

elev <- broom::tidy(m6$model) %>% 
filter(term == "elev") %>%
bind_cols(confint(m6$model)) %>%
 rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)

Developed <- broom::tidy(m1$model) %>% 
  filter(term== "Developed") %>%
  bind_cols(confint(m1$model)) %>%
  rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)

Forest <- broom::tidy(m3$model) %>% 
  filter(term == "Forest") %>%
  bind_cols(confint(m3$model)) %>%
  rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)

Shrub <- broom::tidy(m4$model) %>% 
  filter(term == "Shrub") %>%
  bind_cols(confint(m4$model)) %>%
  rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)

Agriculture <- broom::tidy(m2$model) %>% 
  filter(term == "Agriculture") %>%
  bind_cols(confint(m2$model)) %>%
  rename(est=estimate, l95=`2.5 %`, u95=`97.5 %`)

#' Create df for visualizations
ggdf <- bind_rows(Developed,Forest, Shrub,
                  Agriculture, Dist2Road, elev) 

#' Changing data to factor in order to rename names from rows using levesls
ggdf$term<-as.factor(ggdf$term)
levels(ggdf$term)[levels(ggdf$term) == 'distfromroad'] <- 'Distance from Road'
levels(ggdf$term)[levels(ggdf$term) == 'elev'] <- 'Elevation'
levels(ggdf$term)[levels(ggdf$term) == 'Agriculture'] <- 'Agriculture'
levels(ggdf$term)[levels(ggdf$term) == 'Developed'] <- 'Developed'
levels(ggdf$term)[levels(ggdf$term) == 'Forest'] <- 'Forest'
levels(ggdf$term)[levels(ggdf$term) == 'Shrub'] <- 'Shrub/Scrub'


#' Save RDS file of ggdf for plots
#saveRDS(ggdf, "Working_SSF_Outputs.RDS")



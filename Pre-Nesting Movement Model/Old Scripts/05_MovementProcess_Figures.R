
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
#' **Purpose**: This script visualizes the beta outputs and creates prediction plots from the SSF models
#'

#####################
## Load Packages ##
#####################

require(tidyr)
require(ggplot2)
require(dplyr)

##########################
## Make Plot of Betas ##
##########################

#' Read in RDS file with covs 
final.data <- readRDS("Data Management/RData/Individual-Specific Movement Process/Working_SSF_Ready.RDS") 

#' Read in RDS file of coefficient data
ggdf <- readRDS("Data Management/RData/Individual-Specific Movement Process/Working_SSF_Outputs.RDS")

#' Organize Predictor variables
ggdf$term <- factor(ggdf$term, levels =c("Agriculture", "Developed", "Distance from Road",
                                         "Elevation","Forest", "Shrub/Scrub"))
#' Plot of betas
covs <-ggdf %>% 
  ggplot(mapping = aes(x = term, y = est)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = l95, ymax = u95)
                ,position=position_dodge2(0.9,reverse = TRUE))+
  theme_classic() +
  xlab("Covariate") +
  ylab("Beta Estimate") +
  coord_flip() +
  scale_x_discrete(limits=rev)
covs  

#########################
## Prediction Plots ##
#########################

#' https://stackoverflow.com/questions/47080625/plotting-predictions-from-a-logistic-regression

#' statsmooth: Aids the eye in seeing patterns 
#' method: Smoothing method function to use
#' method.args: Format of in this case the GLM
#' se: Display confidence interval around smooth
#' fullrange: If true the smoothing line gets extended beyond the range of the plot, beyond data potentially

#' Agriculture
Agriculture.Predict <- ggplot(final.data, aes(x= Agriculture, y=case_)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Agriculture", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Agriculture.Predict

#' Deciduous Forest
Forest.Predict <- ggplot(final.data, aes(x= Forest, y=case_)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Forest", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Forest.Predict 

#' Developed
Developed.Predict <- ggplot(final.data, aes(x= Developed, y=case_)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Developed", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Developed.Predict 

#' Shrub/Scrub
Shrub.Predict <- ggplot(final.data, aes(x= Shrub, y=case_)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Shrub/Scrub", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Shrub.Predict 

#' Elevation
 Elevation.Predict <- ggplot(final.data, aes(x= elev, y=case_)) +
  stat_smooth(method="glm",
           method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Elevation", y="Probability of Use")+
  ylim(0,1)+
expand_limits(x=100)
Elevation.Predict 

#' Distance from Road
Dist.Predict <- ggplot(final.data, aes(x= distfromroad, y=case_)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Distance from Road", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Dist.Predict 


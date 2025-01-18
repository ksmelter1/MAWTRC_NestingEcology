
#'---
#' title: Nest-site selection of wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script produces coefficient and prediction plots of the model outputs
#' **Last Updated**: 1/18/2025

#####################
## Load Packages ##
####################

#' Vector of package names
packages <- c("tidyverse")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


##########################
## Data Prep for Plots

#' Load in RData
load("Data Management/RData/Nest-Site Selection/ModelResults/20250118_NimbleResults.RData")

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
    lower = quantile(estimate, 0.05),   # 2.5th percentile
    upper = quantile(estimate, 0.95),   # 97.5th percentile
    .groups = 'drop'
  )

######################
## Density Plot

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

############################
## Data Prep for Beta Plot

#' Calculate mean estimates for each parameter
#' Calculate 95% credible intervals using quantiles
mean_estimates <- samples_long %>%
  dplyr::group_by(parameter) %>%
  summarise(
    mean_estimate = mean(estimate),  # mean of the estimates
    lower = quantile(estimate, 0.05),  # 2.5th percentile (lower bound of credible interval)
    upper = quantile(estimate, 0.95),  # 97.5th percentile (upper bound of credible interval)
    .groups = 'drop'
  ) %>%
  dplyr::filter(parameter != "Intercept")
mean_estimates

mean_estimates <- mean_estimates %>%
  mutate(Scale = case_when(
    parameter %in% c("Percent Grass/Forb", "Percent Woody Vegetation", "Visual Obstruction", 
                     "Basal Area", "Native Woody Vegetation", "Invasive Woody Vegetation") ~ "Nest-Level",
    TRUE ~ "Landscape-Level"
  ))

mean_estimates

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                            levels = c("Visual Obstruction",
                                       "Percent Woody Vegetation", 
                                       "Percent Grass/Forb",
                                       "Native Woody Vegetation",
                                       "Invasive Woody Vegetation",
                                       "Basal Area",
                                       "Grassland/Shrub",
                                       "Mixed Forest",
                                       "Evergreen Forest",
                                       "Deciduous Forest",
                                       "Distance to Primary Road",
                                       "Distance to Secondary Road",
                                       "Agriculture"))) 
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
  scale_color_manual(values = c("Nest-Level" = "#A44200", "Landscape-Level" = "#D65F5F")) +  # Set color for Microscale, Macroscale
  scale_shape_manual(values = c("Nest-Level" = 17, "Landscape-Level" = 16))+  # Set shapes for Microscale, Macroscale
  theme(
    axis.title.x = element_text(margin = margin(t = 10))  # Pad the x-axis label by 10 points (~0.1 inch)
  )
  p2.betas
  
  #' Save the combined microscale plot as a PNG image
  ggsave("01_BetaEstimates_90CIs.png", 
         p2.betas, 
         path = "Data Management/Figures/Nest-Site Selection Process/Beta Outputs/",
         width = 12, 
         height = 10,
         bg = "white")

#' #########################
#' ## Prediction Plots
#' 
#' #' Distance to Primary Road
#' CIs <- matrix(quantile(nimbleMCMC_samples[[1]][,2], c(0.05, 0.5, 0.95)), nrow=1)
#' pred.seq<-seq(-1,1,0.1)
#' pred.response.log<-pred.seq%*%CIs
#' pred.response<-data.frame(exp(pred.response.log))
#' colnames(pred.response)<-c("lCI","median","uCI")
#' 
#' p3.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
#'   geom_ribbon(aes(ymin=lCI, ymax=uCI), linetype=2, alpha=0.1) +
#'   theme_minimal() +
#'   xlab("Distance to Primary Road") +
#'   ylab("Relative Probability of Use") +
#'   theme(
#'     axis.title.x = element_text(margin = margin(t = 10)), # Pad the x-axis label by 10 points (~0.1 inch)
#'     axis.title.y = element_text(margin = margin(t = 10)) 
#'     )
#' p3.predict
#' 
#' #' Agriculture
#' CIs<-matrix(quantile(jags_samples[[1]][,7],c(0.05,0.5,0.95)),nrow=1)
#' pred.seq<-seq(-1,1,0.1)
#' pred.response.log<-pred.seq%*%CIs
#' pred.response<-data.frame(exp(pred.response.log))
#' colnames(pred.response)<-c("lCI","median","uCI")
#' 
#' p8.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
#'   geom_ribbon(aes(ymin=lCI, ymax=uCI), linetype=2, alpha=0.1) +
#'   theme_minimal() +
#'   xlab("Agriculture") +
#'   ylab("Relative Probability of Use") +
#'   theme(
#'     axis.title.x = element_text(margin = margin(t = 10)), # Pad the x-axis label by 10 points (~0.1 inch)
#'     axis.title.y = element_text(margin = margin(r = 10)) 
#'   )
#' p8.predict
#' 
#' #' Visual Obstruction
#' CIs<-matrix(quantile(jags_samples[[1]][,9],c(0.05,0.5,0.95)),nrow=1)
#' pred.seq<-seq(-1,1,0.1)
#' pred.response.log<-pred.seq%*%CIs
#' pred.response<-data.frame(exp(pred.response.log))
#' colnames(pred.response)<-c("lCI","median","uCI")
#' 
#' p10.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
#'   geom_ribbon(aes(ymin=lCI, ymax=uCI), linetype=2, alpha=0.1) +
#'   theme_minimal() +
#'   xlab("Visual Obstruction") +
#'   ylab("Relative Probability of Use") +
#'   theme(
#'     axis.title.x = element_text(margin = margin(t = 10)), # Pad the x-axis label by 10 points (~0.1 inch)
#'     axis.title.y = element_text(margin = margin(r = 10)) 
#'   )
#' p10.predict
#' 
#' #' Percent Woody Vegetation
#' CIs<-matrix(quantile(jags_samples[[1]][,11],c(0.05,0.5,0.95)),nrow=1)
#' pred.seq<-seq(-1,1,0.1)
#' pred.response.log<-pred.seq%*%CIs
#' pred.response<-data.frame(exp(pred.response.log))
#' colnames(pred.response)<-c("lCI","median","uCI")
#' 
#' p12.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
#'   geom_ribbon(aes(ymin=lCI, ymax=uCI), linetype=2, alpha=0.1) +
#'   theme_minimal() +
#'   xlab("Percent Woody Vegetation") +
#'   ylab("Relative Probability of Use") +
#'   theme(
#'     axis.title.x = element_text(margin = margin(t = 10)), # Pad the x-axis label by 10 points (~0.1 inch)
#'     axis.title.y = element_text(margin = margin(r = 10)) 
#'   )
#' p12.predict
#' 
#' #' Invasive Woody Vegetation
#' CIs<-matrix(quantile(jags_samples[[1]][,13],c(0.05,0.5,0.95)),nrow=1)
#' pred.seq<-seq(-1,1,0.1)
#' pred.response.log<-pred.seq%*%CIs
#' pred.response<-data.frame(exp(pred.response.log))
#' colnames(pred.response)<-c("lCI","median","uCI")
#' 
#' p14.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
#'   geom_ribbon(aes(ymin=lCI, ymax=uCI), linetype=2, alpha=0.1) +
#'   theme_minimal() +
#'   xlab("Invasive Woody Vegetation") +
#'   ylab("Relative Probability of Use") +
#'   theme(
#'     axis.title.x = element_text(margin = margin(t = 10)), # Pad the x-axis label by 10 points (~0.1 inch)
#'     axis.title.y = element_text(margin = margin(r = 10)) 
#'   )
#' p14.predict
#' 
#' #' Grassland/Shrub
#' CIs<-matrix(quantile(jags_samples[[1]][,15],c(0.05,0.5,0.95)),nrow=1)
#' pred.seq<-seq(-1,1,0.1)
#' pred.response.log<-pred.seq%*%CIs
#' pred.response<-data.frame(exp(pred.response.log))
#' colnames(pred.response)<-c("lCI","median","uCI")
#' 
#' p16.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
#'   geom_ribbon(aes(ymin=lCI, ymax=uCI), linetype=2, alpha=0.1) +
#'   theme_minimal() +
#'   xlab("Grassland/Shrub") +
#'   ylab("Relative Probability of Use") +
#'   theme(
#'     axis.title.x = element_text(margin = margin(t = 10)), # Pad the x-axis label by 10 points (~0.1 inch)
#'     axis.title.y = element_text(margin = margin(r = 10)) 
#'   )
#' p16.predict
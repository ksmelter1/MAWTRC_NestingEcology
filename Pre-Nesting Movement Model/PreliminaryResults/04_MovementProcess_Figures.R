#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: RandomSteps_Prep(R Workspace)
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script creates beta estimate and prediction plots similar to the nest-site selection model
#' **Last Updated**: 1/20/25

####################
## Load Packages 

#' Vector of package names
packages <- c("tidyverse",
              "gridExtra")


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


##########################
## Data Prep for Plots

samples_df <- readRDS("Data Management/RData/Pre-Nesting Movement Model/RData Files/Draft3/20250203_samples_df.RDS")

#' Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

#' View the reshaped data
head(samples_long)

#' Calculate 90% Bayesian credible intervals
credible_intervals <- samples_long %>%
  group_by(parameter) %>%
  summarise(
    lower = quantile(estimate, 0.05),   
    upper = quantile(estimate, 0.95), 
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
    mean_estimate = mean(estimate),  
    lower = quantile(estimate, 0.05),  
    upper = quantile(estimate, 0.95),  
    .groups = 'drop'
  ) %>%
  dplyr::filter(parameter != "Intercept") %>%
  dplyr::filter(parameter != "Elevation")
mean_estimates

#' Assign predictors to scales
mean_estimates <- mean_estimates %>%
  dplyr::mutate(Scale = "Landscape") 
  
                                          

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c("Grassland/Shrub",
                                              "Mixed Forest",
                                              "Evergreen Forest",
                                              "Deciduous Forest",
                                              "Agriculture",
                                              "Distance to Primary Road",
                                              "Distance to Secondary Road"
                                            ))) 
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
  scale_color_manual(values = c("Landscape" = "#D65F5F")) +  # Set color for Microscale, Macroscale
  scale_shape_manual(values = c("Landscape" = 16))+  # Set shapes for Microscale, Macroscale
  theme(
    axis.title.x = element_text(margin = margin(t = 10))  # Pad the x-axis label by 10 points (~0.1 inch)
  )
p2.betas

#' Save the combined microscale plot as a PNG image
# ggsave("01_BetaEstimates_90CIs.png", 
#        p2.betas, 
#        path = "Data Management/Figures/Nest-Site Selection Process/Beta Outputs/",
#        width = 12, 
#        height = 10,
#        bg = "white")


################################################################################
## Prediction Plots

#' Distance to Primary Road
CIs <- matrix(quantile(samples_df[,1], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p1.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#D65F5F", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Distance to Primary Road") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_text(margin = margin(t = 10), vjust = 3))
p1.predict

#' Mixed Forest
CIs <- matrix(quantile(samples_df[,4], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p3.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#D65F5F", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Mixed Forest") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_text(margin = margin(t = 10), vjust = 3))
p3.predict


#' Evergreen Forest
CIs <- matrix(quantile(samples_df[,5], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p4.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#D65F5F", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Evergreen Forest") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_text(margin = margin(t = 10), vjust = 3))
p4.predict


#' Deciduous Forest
CIs <- matrix(quantile(samples_df[,6], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p5.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#D65F5F", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Deciduous Forest") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_text(margin = margin(t = 10), vjust = 3))
p5.predict


#' Agriculture
CIs <- matrix(quantile(samples_df[,7], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p6.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#D65F5F", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Agriculture") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_text(margin = margin(t = 10), vjust = 3))
p6.predict


#' Grassland/Shrub
CIs <- matrix(quantile(samples_df[,8], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p7.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#D65F5F", linetype=2, alpha=0.1) +
  theme_minimal() +
  ylab("Relative Probability of Use") +
  xlab("Grassland/Shrub") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_text(margin = margin(t = 10), vjust = 3))
p7.predict


#########################
## Create panels

#' Define a theme to remove y-axis title, labels, and ticks
remove_y_axis <- theme(
  axis.title.y = element_blank())

#' Adjust the y-axis title margin for p4.predict
p4.predict <- p4.predict + theme(
  axis.title.y = element_text(margin = margin(r = 20)))

#' Adjust the y-axis title margin for p4.predict
p1.predict <- p1.predict + theme(
  axis.title.y = element_text(margin = margin(r = 20)))

#' Adjust the y-axis title margin for p4.predict
p6.predict <- p6.predict + theme(
  axis.title.y = element_text(margin = margin(r = 20)))

#' Apply the theme to the right-hand plots
p1.predict <- p1.predict + remove_y_axis
p3.predict <- p3.predict + remove_y_axis
p5.predict <- p5.predict + remove_y_axis
p6.predict <- p6.predict + remove_y_axis
p7.predict <- p7.predict + remove_y_axis

#' Arrange the plots
grid.arrange(p1.predict, p3.predict, p4.predict, p5.predict, p6.predict, p7.predict,
             ncol = 2, nrow = 3)

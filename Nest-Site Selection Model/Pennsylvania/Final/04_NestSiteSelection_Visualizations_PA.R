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
#' **Last Updated**: 2/3/2025

################################################################################
## Load Packages 

#' Vector of package names
packages <- c("tidyverse",
              "ggpubr")


#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


################################################################################
## Data Prep for Plots

#' Load in RData
load("Data Management/RData/Nest-Site Selection/ModelResults/PreliminaryResults/20250206_NimbleResults.RData")



################################################################################
## Data Prep for Beta Plot

#' Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

#' Calculate mean estimates for each parameter
#' Calculate 90% credible intervals using quantiles
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
                     "Visual Obstruction", 
                     "Woody Stem Count",
                     "Percent Fern") ~ "Nest",
    TRUE ~ "Landscape"
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
                                       "Developed",
                                       "Agriculture",
                                       "Distance to Secondary Road",
                                       "Distance to Primary Road",
                                       "Visual Obstruction",
                                       "Woody Stem Count",
                                       "Percent Woody Vegetation", 
                                       "Percent Grass/Forb",
                                       "Percent Fern"))) 
mean_estimates


################################################################################
## Betas Plot 

mean_estimates <- mean_estimates %>%
  mutate(Scale = factor(Scale, levels = c("Nest", "Landscape")))

#' Figure for display in powerpoint presentations
p2.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  ylab("Selection Relative to Deciduous Forest") +
  theme_minimal() + 
  coord_flip() +
  scale_color_manual(values = c("Nest" = "#A44200", "Landscape" = "#D65F5F")) +  
  scale_shape_manual(values = c("Nest" = 17, "Landscape" = 16)) +  
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.38), 
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
    axis.title.y = element_blank(),
        legend.position = "none")  
p2.betas

p2.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5, stroke = 1.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  ylab("Posterior Mean") +
  xlab("Parameter") +
  theme_minimal() + 
  coord_flip() +
  scale_color_manual(values = c("Nest" = "#A44200", "Landscape" = "#D65F5F")) +  
  scale_shape_manual(values = c("Nest" = 17, "Landscape" = 16)) +  
  theme(
    axis.title.x = element_text(margin = margin(t = 10), hjust = 0.55),
    axis.title.y = element_text(margin = margin(r = 10))
  )

p2.betas


################################################################################
## Prediction Plots

#' Max Height Visual Obstruction
CIs <- matrix(quantile(samples_df[,8], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p2.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#A44200", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Aerial Visual Obstruction") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_blank())
p2.predict


#' Horizontal Visual Obstruction
CIs <- matrix(quantile(samples_df[,11], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p3.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#A44200", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Horizontal Visual Obstruction") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_blank())
p3.predict


################################################################################
## Create panels

plot_combined <- ggarrange(p2.predict, p3.predict, 
                           nrow = 1, ncol = 2) %>% 
  annotate_figure(left = text_grob("Relative Probability of Use", rot = 90))
plot_combined


################################################################################
################################################################################X
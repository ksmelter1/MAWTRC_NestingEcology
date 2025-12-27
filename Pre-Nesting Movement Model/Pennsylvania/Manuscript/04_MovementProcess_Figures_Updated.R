#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script creates beta estimate and prediction plots for the pre-nesting movement model
#' **Last Updated**: 12/27/25

################################################################################
## Load Packages and Data

# Vector of package names
packages <- c("tidyverse",
              "gridExtra",
              "stringr")

# Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

# Apply the function to each package name
lapply(packages, load_packages)

# Read in named samples from Pennsylvania pre-nesting movement model
samples_df <- readRDS("Data Management/RData/Pre-Nesting Movement Model/Pennsylvania/Model Results/20250513_samples_df.RDS")


################################################################################
## Data Prep for Beta Plot

# Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

# View the reshaped data
head(samples_long)

# Calculate mean estimates for each parameter
# Calculate 90% credible intervals using quantiles
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

# Assign predictors to scales
mean_estimates <- mean_estimates %>%
  dplyr::mutate(Scale = "Landscape") 
  
                                          
# Organize variables into levels to be displayed
# Filter out the intercept
mean_estimates1 <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c("Wetland",
                                              "Grassland/Shrub",
                                              "Mixed Forest",
                                              "Evergreen Forest",
                                              "Developed",
                                              "Crop",
                                              "Pasture",
                                              "Distance to Primary Road",
                                              "Distance to Secondary Road"
                                            ))) 

################################################################################
## Beta Plot 

# Figure for Presentations
# Beta estimates and associated 90% credible intervals 
p2.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Selection Relative to Deciduous Forest") +
  theme_minimal() + 
  coord_flip()+
  scale_color_manual(values = c("Landscape" = "#D65F5F")) +  
  scale_shape_manual(values = c("Landscape" = 16))+  
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.45),  
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )
p2.betas

# Figure for Pennsylvania publication
p2.betas <- ggplot(mean_estimates1, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5, stroke = 1.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Posterior Mean") +
  theme_minimal() + 
  coord_flip()+
  scale_color_manual(values = c("Landscape" = "#D65F5F")) +  
  scale_shape_manual(values = c("Landscape" = 16))+  
  theme(
    axis.title.x = element_text(margin = margin(t = 10), hjust = 0.45))
p2.betas

# Save as outputs as RData object
save(p2.betas, 
     mean_estimates1,
     samples_df,
     file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Pre-Nesting Movement/pre.pa.betas.updated.RData", 
     overwrite = T)


################################################################################
## Prediction Plots

# Distance to Primary Road
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

################################################################################
###############################################################################X

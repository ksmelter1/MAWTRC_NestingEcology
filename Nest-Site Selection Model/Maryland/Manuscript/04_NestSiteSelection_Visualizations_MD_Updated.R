#'---
#' title: Nest-site selection of female wild turkeys in Maryland
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script produces coefficient and prediction plots of the model outputs
#' **Last Updated**: 12/27/2025

################################################################################
## Load Packages 

# Vector of package names
packages <- c("tidyverse",
              "ggpubr")


# Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

# Apply the function to each package name
lapply(packages, load_packages)


################################################################################
## Data Prep for Plots

# Load in RData
load("Data Management/RData/Nest-Site Selection/Maryland/Model Results/Manuscript/NSS_Results_MD.RData")


################################################################################
## Data Prep for Beta Plot

# Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

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
  dplyr::filter(parameter != "Intercept") 
mean_estimates

# Assign covariate classes
mean_estimates <- mean_estimates %>%
  mutate(Scale = case_when(
    parameter %in% c("Percent Grass/Forb", "Percent Woody Vegetation", "Horizontal Visual Obstruction", 
                     "Aerial Visual Obstruction", "Percent Fern") ~ "Nest",
    TRUE ~ "Landscape"
  ))
mean_estimates

# Organize variables into levels to be displayed
# Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c("Wetland",
                                              "Grassland/Shrub",
                                              "Mixed Forest",
                                              "Evergreen Forest",
                                              "Developed",
                                              "Pasture",
                                              "Crop",
                                              "Distance to Primary Road",
                                              "Distance to Secondary Road" ))) 
mean_estimates


################################################################################
## Betas Plot 

# Figure for display in powerpoint presentations
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

# Figure for Maryland publication
p2.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5, stroke = 1.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  ylab("Selection Relative to Deciduous Forest") +
  xlab("Parameter") +
  theme_minimal() + 
  coord_flip() +
  scale_color_manual(values = c("Nest" = "#A44200", "Landscape" = "#D65F5F")) +  
  scale_shape_manual(values = c("Nest" = 17, "Landscape" = 16)) +  
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust = 0.50))
p2.betas

save(p2.betas, 
     mean_estimates,
     samples_df,
     file = "MAWTRC Nesting Ecology Manuscript/Figure Code and Data/RData/Nest-Site Selection/nss.MD.betas.updated.RData")

#################################################################################
################################################################################X
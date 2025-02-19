#'---
#' title: Nest-site selection of female wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script creates a Bayesian conditional logistic regression model for nest-site selection in JAGs using the gathered covariates 
#' **Last Updated**: 2/4/25

#' The PGC is interested in a nest-site selection model with no NLCD reference level

#####################
## Load Packages 

#' Vector of package names
packages <- c("matrixStats",
              "nimble",
              "MCMCvis",
              "tidyverse",
              "sf",
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
## Data Preparation 

#' Load in RData
load("Data Management/RData/Nest-Site Selection/Covs/Draft5/20250131_Covs.RData")


#' Create nest.data object for models
nest.data <- pa.nests.covs 
str(pa.nests.covs)


#' Select columns of interest
nest.data <- nest.data %>%
  dplyr::select(NestID, 
                BandID, 
                Case,
                PercGrassForb,
                PercWoody, 
                PercLitter,
                PercBare,
                PercBoulder,
                AvgMaxVO, 
                primary,
                secondary,
                StemCount,
                WoodyType1,
                AvgVO,
                PercFern, 
                GuardHt,
                HtWoody,
                HtFern,
                HtGrassForb,
                HtLitter, 
                HtBoulder)

#' Switch coding to UTF-8 which is a widely used character encoding system
nest.data <- nest.data %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

# Convert selected columns back to numeric
nest.data$PercWoody <- as.numeric(nest.data$PercWoody)
nest.data$PercGrassForb <- as.numeric(nest.data$PercGrassForb)
nest.data$AvgMaxVO <- as.numeric(nest.data$AvgMaxVO)
nest.data$PercFern <- as.numeric(nest.data$PercFern)
nest.data$AvgVO <- as.numeric(nest.data$AvgVO)
nest.data$StemCount <- as.numeric(nest.data$StemCount)
nest.data$GuardHt <- as.numeric(nest.data$GuardHt)
nest.data$HtWoody <- as.numeric(nest.data$HtWoody)
nest.data$HtFern <- as.numeric(nest.data$HtFern)
nest.data$HtGrassForb <- as.numeric(nest.data$HtGrassForb)
nest.data$HtLitter <- as.numeric(nest.data$HtLitter)
nest.data$HtBoulder <- as.numeric(nest.data$HtBoulder)
nest.data$PercLitter <- as.numeric(nest.data$PercLitter)
nest.data$PercBoulder <- as.numeric(nest.data$PercBoulder)
nest.data$PercBare <- as.numeric(nest.data$PercBare)
str(nest.data)
glimpse(nest.data)

# Scale relevant columns
nest.data$PercWoody <- scale(nest.data$PercWoody)
nest.data$PercGrassForb <- scale(nest.data$PercGrassForb)
nest.data$AvgMaxVO <- scale(nest.data$AvgMaxVO)
nest.data$PercFern <- scale(nest.data$PercFern)
nest.data$AvgVO <- scale(nest.data$AvgVO)
nest.data$StemCount <- scale(nest.data$StemCount)
nest.data$GuardHt <- scale(nest.data$GuardHt)
nest.data$HtWoody <- scale(nest.data$HtWoody)
nest.data$HtFern <- scale(nest.data$HtFern)
nest.data$HtGrassForb <- scale(nest.data$HtGrassForb)
nest.data$HtLitter <- scale(nest.data$HtLitter)
nest.data$HtBoulder <- scale(nest.data$HtBoulder)
nest.data$PercLitter <- scale(nest.data$PercLitter)
nest.data$PercBoulder <- scale(nest.data$PercBoulder)
nest.data$PercBare <- scale(nest.data$PercBare)
glimpse(nest.data)

# Convert selected columns back to numeric
nest.data$PercWoody <- as.numeric(nest.data$PercWoody)
nest.data$PercGrassForb <- as.numeric(nest.data$PercGrassForb)
nest.data$AvgMaxVO <- as.numeric(nest.data$AvgMaxVO)
nest.data$PercFern <- as.numeric(nest.data$PercFern)
nest.data$AvgVO <- as.numeric(nest.data$AvgVO)
nest.data$StemCount <- as.numeric(nest.data$StemCount)
nest.data$GuardHt <- as.numeric(nest.data$GuardHt)
nest.data$HtWoody <- as.numeric(nest.data$HtWoody)
nest.data$HtFern <- as.numeric(nest.data$HtFern)
nest.data$HtGrassForb <- as.numeric(nest.data$HtGrassForb)
nest.data$HtLitter <- as.numeric(nest.data$HtLitter)
nest.data$HtBoulder <- as.numeric(nest.data$HtBoulder)
nest.data$PercLitter <- as.numeric(nest.data$PercLitter)
nest.data$PercBoulder <- as.numeric(nest.data$PercBoulder)
nest.data$PercBare <- as.numeric(nest.data$PercBare)
glimpse(nest.data)


################################################################################
## Consolidate Invasive vs Native Woody Vegetation

#' If else statement to define Native vs. Invasive Woody vegetation
nest.data$WoodyType1 <- ifelse(nest.data$WoodyType1 %in% c("Non-native invasive", 
                                                           "Non-native Invasive"), 
                               yes = "Invasive",
                               no = "Native")

print(nest.data$WoodyType1)
table(nest.data$WoodyType1)

#' Create dummy variables 
nest.data$Invasive <- ifelse(nest.data$WoodyType1 == "Invasive", 1, 0)
nest.data$Native <- ifelse(nest.data$WoodyType1 == "Native", 1, 0)
glimpse(nest.data)

#' Drop geometry
#' Make case numeric
#' Remove primary and secondary road columns 
nest.data <- nest.data %>%
  st_drop_geometry() %>%
  dplyr::mutate(Case = as.numeric(Case)) %>%
  dplyr::select(-primary, -secondary)
glimpse(nest.data)

################################################################################
## Data Prep

#' Order df by nestid
nest.data <- nest.data[order(nest.data$NestID),]

#' Change Nest_ID_V to numeric 
#' First must remove underlines
nest.data$NestID <- gsub("_", "", nest.data$NestID)
nest.data$NestID <- as.numeric(nest.data$NestID)

#' Create str_id column 
#' This allows the loop to iterate through the steps associated with each bird 
#' cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  group_by(NestID) %>%
  mutate(str_ID=cur_group_id())

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
glimpse(nest.data.ready)


################################################################################
## Nimble Model

nestmodel<-nimbleCode({
  for (i in 1:I){
    use[i]~dpois(lambda[i])
    
    log(lambda[i])<-inprod(beta[1:J],X[i,1:J])+alpha[str_ID[i]]
  }
  
  #' Priors
  for(j in 1:J){
    beta[j]~dnorm(0,sd=sqrt(1/0.0001))
  }
  
  for(k in 1:K){
    alpha[k]~dnorm(0,sd=sqrt(1/0.000001))
  }
}
)

#' Model parameters
X<- cbind(
  rep(1, nrow(nest.data.ready)),   # Intercept (1)
  nest.data.ready$AvgMaxVO,        # Horizontal Visual Obstruction
  nest.data.ready$PercGrassForb,   # Percent Grass/Forb
  nest.data.ready$PercWoody,       # Percent Woody
  nest.data.ready$AvgVO,           # Aerial Visual Obstruction
  nest.data.ready$PercFern,        # Percent Fern
  nest.data.ready$StemCount,       # Woody Stem Count
  nest.data.ready$Invasive,        # Invasive Woody Vegetation, Native is the reference level
  nest.data.ready$HtWoody,         # Woody Vegetation Height 
  nest.data.ready$PercLitter,      # Percent Leaf Litter
  nest.data.ready$PercBare         # Percent Bare Ground 
)

Consts <- list(str_ID = as.numeric(nest.data.ready$str_ID),
               J=ncol(X),
               I = nrow(nest.data.ready),
               K =length(unique(nest.data.ready$str_ID)))

#' Data list for nimble
Data<-list(X=X,use = nest.data.ready$Case)

Inits<-list(beta=rep(0,ncol(X)),alpha=rep(0,length(unique(nest.data.ready$str_ID))))

nimbleMCMC_samples <- nimbleMCMC(code = nestmodel, 
                                 constants = Consts, 
                                 inits=Inits,
                                 data = Data,
                                 nburnin = 10000,
                                 niter = 20000,
                                 thin=1)

colMeans(nimbleMCMC_samples[,171:181])
colSds(nimbleMCMC_samples[,171:181])

#' View traceplots
MCMCtrace(nimbleMCMC_samples, pdf = F)

#' Extract the posterior samples for the 'beta' parameters (columns 171 to 181)
beta_samples <- nimbleMCMC_samples[, 171:181]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", 
               "Aerial Visual Obstruction", 
               "Percent Grass Forb",
               "Percent Woody Vegetation",
               "Horizontal Visual Obstruction",
               "Percent Fern",
               "Woody Stem Count",
               "Invasive Woody Vegetation",
               "Woody Vegetation Height",
               "Percent Leaf Litter",
               "Percent Bare Ground")
               
#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)

#' Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 


################################################################################
## Data Prep for Beta Plot

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

#' Assign predictors to scales
mean_estimates <- mean_estimates %>%
  dplyr::mutate(Scale = "Nest")

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c("Intercept", 
                                              "Aerial Visual Obstruction", 
                                              "Horizontal Visual Obstruction",
                                              "Percent Grass Forb",
                                              "Percent Woody Vegetation",
                                              "Percent Fern",
                                              "Woody Stem Count",
                                              "Invasive Woody Vegetation",
                                              "Woody Vegetation Height",
                                              "Percent Leaf Litter",
                                              "Percent Bare Ground"))) 
mean_estimates


################################################################################
## Betas Plot 

p2.betas <- ggplot(mean_estimates, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5, fill = "white") +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Beta Estimate") +
  theme_minimal() + 
  coord_flip() +
  scale_color_manual(values = c("Nest" = "#A44200")) +  
  scale_shape_manual(values = c("Nest" = 17)) +  
  theme(
    axis.title.x = element_text(margin = margin(t = 10), hjust = 0.42)) 
p2.betas


################################################################################
## Prediction Plots

#' Woody Stem Count
CIs <- matrix(quantile(samples_df[,7], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p1.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#A44200", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Woody Stem Count") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_blank())
p1.predict

#' Percent Fern
CIs <- matrix(quantile(samples_df[,6], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p2.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#A44200", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Percent Fern") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_blank())
p2.predict


#' Horizontal Visual Obstruction
CIs <- matrix(quantile(samples_df[,5], c(0.05, 0.5, 0.95)), nrow=1)
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

#' Aerial Visual Obstruction
CIs <- matrix(quantile(samples_df[,2], c(0.05, 0.5, 0.95)), nrow=1)
pred.seq<-seq(-1,1,0.1)
pred.response.log<-pred.seq%*%CIs
pred.response<-data.frame(exp(pred.response.log))
colnames(pred.response)<-c("lCI","median","uCI")
p4.predict <- ggplot(data=pred.response, aes(x=pred.seq, y=median))+geom_line() +
  geom_ribbon(aes(ymin=lCI, ymax=uCI), fill = "#A44200", linetype=2, alpha=0.1) +
  theme_minimal() +
  xlab("Aerial Visual Obstruction") +
  ylab("Relative Probability of Use") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_blank())
p4.predict

#' Arrange the plots into a panel
#grid.arrange(p1.predict, p2.predict, p3.predict, p4.predict, ncol = 2)

plot_combined <- ggarrange(p1.predict, p2.predict, p3.predict, p4.predict,
                           nrow = 2, ncol = 2) %>% 
  annotate_figure(left = text_grob("Relative Probability of Use", rot = 90))
plot_combined


#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script fits SSF models in frequentist and Bayesian frameworks and extracts coefficient data
#' #' **Last Updated**: 12/4/24


####################
## Load Packages 

#' Vector of package names
packages <- c("R2jags",
              "MCMCvis",
              "dplyr",
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


####################
## Organize Data 

#' Check
str(dat_2.ready)

#' Change ID to numeric
dat_2.ready$id <- gsub("_", "", dat_2.ready$id)
dat_2.ready$id<- as.numeric(dat_2.ready$id)

#' Create strata column by merging BirdID and step_id
dat_2.ready$NA_ID <- paste(dat_2.ready$id,
                           dat_2.ready$step_id_,
                           sep = "_")

#' Change NA_ID to integer
dat_2.ready$NA_ID <- gsub("_", "", dat_2.ready$NA_ID)
dat_2.ready$NA_ID<- as.numeric(dat_2.ready$NA_ID)

class(dat_2.ready)
str(dat_2.ready)

#' Add numerical variable for animals:
dat_2.ready$ANIMAL_ID <- as.numeric(as.factor(dat_2.ready$id))

#' Stratum ID is given as "NA_ID" in the data; 
#' It is easier to have sequential enumeration, so let's generate a new stratum-ID variable str_ID:
d.map <- data.frame(NA_ID=unique(dat_2.ready$NA_ID),
                    str_ID=1:length(unique(dat_2.ready$NA_ID)))
dat_2.ready$str_ID <- d.map[match(dat_2.ready$NA_ID,d.map$NA_ID),"str_ID"]
dat_2.ready <- dat_2.ready[order(dat_2.ready$str_ID),] 
glimpse(dat_2.ready)

#' Function to fill NA values with 0 throughout dataframe
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


###############################################
## JAGS Individual Movement Model

#' Written in matrix notation 
#' Does not include developed as a predictor due to the fact it is correlated with agriculture

JAGS.mod <- "  
 model{
  for (i in 1:I){
   use[i]~dpois(lambda[i])

   log(lambda[i])<-inprod(beta,X[i,])+alpha[str_ID[i]]
  }
  
  #Priors
  for(j in 1:J){
   beta[j]~dnorm(0,0.0001)
  }
  
  for(k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
  
 }
"
#' Model parameters
X <- cbind(
  rep(1, nrow(dat_2.ready)),  # Intercept (1)
  dat_2.ready$primary,         # Distance to Primary Road
  dat_2.ready$secondary,       # Distance to Secondary Road
  dat_2.ready$Mixed,           # Mixed Forest
  dat_2.ready$Evergreen,       # Evergreen Forest
  dat_2.ready$Deciduous,       # Deciduous Forest
  dat_2.ready$Agriculture,     # Agriculture
  dat_2.ready$elev,            # Elevation
  dat_2.ready$Grassland)       # Grassland    
  

#' Data prep for JAGS
jags_data <- list(use = dat_2.ready$case_,
                  X=X,
                  str_ID=as.numeric(dat_2.ready$str_ID),
                  I=nrow(dat_2.ready),
                  J=ncol(X),
                  K=length(unique(dat_2.ready$str_ID)))

#' Initialize JAGS
start <- Sys.time()
jags_model <- jags.model(textConnection(JAGS.mod), 
                         data = jags_data,
                         n.chains = 1,
                         n.adapt=2000)

#' Sample from posterior distribution using parallel
jags_samples <- coda.samples(jags_model, c("beta"),
                             n.iter = 40000, 
                             n.burnin = 10000,
                             n.thin = 1,
                             do.parallel = T)
end <- Sys.time()

colMeans(jags_samples[[1]])
colSds(jags_samples[[1]])

#' View traceplots
MCMCtrace(jags_samples, pdf = F)

#' Convert mcmc.list to matrix
samples_matrix <- as.matrix(jags_samples)

#' Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' View the first few rows of the samples data frame
head(samples_df)

#' Check the actual column names to confirm their structure
colnames(samples_df)

#' Create a vector of new names
new_names <- c("Intercept", "Distance to Primary Road", "Distance to Secondary Road", 
               "Mixed Forest", "Evergreen Forest", "Deciduous Forest", "Agriculture", 
               "Elevation", "Grassland/Shrub")

#' Rename the columns based on position
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)

save(jags_samples, samples_df, 
     "Data Management/RData/Individual-Specific Movement Process/RData Files/03_MovementProcess_ModelOutputs.RData")

################################################################################
################################################################################

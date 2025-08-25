
rm(list=ls())

setwd("C:/Users/5978386/Documents/PhD/predict_detection_zones/")


# Description
# This script uses the camera trap detection zone database from the following paper
# "Inferring camera trap detection zones for rare species using species- and camera-specific traits: a meta-level analysis"
# by Wiegers et al. (2025) in Remote Sensing in Ecology and Evolution

# It predicts effective detection radius (EDD) and effective detection angle (EDA) based on several covariates:

# body_mass: in kilograms
# ct_height: in meters
# ct_brand: character, either Browning, Reconyx, Bushnell, LIttle Acorn, or Spromise
# T: snapshot interval (see manuscript) in seconds
# data.order: character, either Aves, Carnivora, Primates, Rodentia, or Ungulata. 
# max_r_site, right-truncated distance / maximum considered distance for a study site, in meters

# Input the values relevant to your species in the 'predict.csv' file, which the script then loads and predicts a detection zone for

# Packages required
library(readr)
library(nlme)

#####################

# SCENARIO 1 (AND 3)

# estimating the detection zone for a species at a study site, where detection zones for more common
# species are available, but not for the target species (1) or for any species at a novel study site 
# for which no detection zone measurements are available at all (3)

################ EFFECTIVE DETECTION RADIUS 

# Load in the model data
LOOCV_data <- read_csv("LOOCV_data 2025.csv")
LOOCV_data <- LOOCV_data[,c("T", "site", "Scientific_name", "ct_brand", 
                            "ct_height", "body_mass", 
                            "max_r_site", "data.family", "data.genus", "data.order", 
                            "data.species","full.r")]

# Prepare data for lm()
LOOCV_data$ct_height <- as.numeric(LOOCV_data$ct_height)
LOOCV_data$data.order <- as.factor(LOOCV_data$data.order)
LOOCV_data$T[LOOCV_data$T=="contact"] <- 1500 # assign 1500 (highest non-manual value) as snapshot interval when determined manually 
LOOCV_data$T <- log(as.numeric(LOOCV_data$T))
LOOCV_data$ct_brand <- factor(LOOCV_data$ct_brand)

# Train the model
mlm_r <- lme(
  full.r ~ log(body_mass)+ ct_height + ct_brand + T + data.order + max_r_site,
  random = ~1 | site,
  data = LOOCV_data
)

# Add your own data in the predict.csv file for which we want to predict the detection zone radius for.

predict <- read_csv("predict.csv")

# Prepare data for lm()
predict$ct_height <- as.numeric(predict$ct_height)
predict$data.order <- as.factor(predict$data.order)
predict$T[predict$T=="contact"] <- 1500 # assign 1500 (highest non-manual value) as snapshot interval when determined manually 
predict$T <- log(as.numeric(predict$T))
predict$ct_brand <- factor(predict$ct_brand)

# Predict using the fitted model
predictions <- predict(
  mlm_r, 
  newdata = predict,
  level   = 0,                # level = 0 → only fixed effects
  allow.new.levels = TRUE     # allows unseen levels in random effects
)

predict$model_estimate <- predictions

# Save
predict_r <- predict

################ EFFECTIVE DETECTION ANGLE #####################

LOOCV_data_theta <- read_csv("LOOCV_data_theta.csv")
LOOCV_data_theta <- LOOCV_data_theta[,c("T", "site", "Scientific_name", "ct_brand", 
                            "ct_height", "body_mass", 
                            "max_theta_site", "data.family", "data.genus", "data.order", 
                            "data.species","full.theta")]

# Prepare data for lm()
LOOCV_data_theta$ct_height <- as.numeric(LOOCV_data_theta$ct_height)
LOOCV_data_theta$data.order <- as.factor(LOOCV_data_theta$data.order)
LOOCV_data_theta$T[LOOCV_data_theta$T=="contact"] <- 1500 # assign 1500 (highest non-manual value) as snapshot interval when determined manually 
LOOCV_data_theta$T <- log(as.numeric(LOOCV_data_theta$T))
LOOCV_data_theta$ct_brand <- factor(LOOCV_data_theta$ct_brand)

# Train the model
mlm_theta <- lme(
  full.theta ~ log(body_mass)+ ct_height + ct_brand + T + data.order + max_theta_site,
  random = ~1 | site,
  data = LOOCV_data_theta
)

# Add your own data in the predict.csv file for which we want to predict the detection zone angle for.

predict <- read_csv("predict_theta.csv")

# Prepare data for lm()
predict$ct_height <- as.numeric(predict$ct_height)
predict$data.order <- as.factor(predict$data.order)
predict$T[predict$T=="contact"] <- 1500 # assign 1500 (highest non-manual value) as snapshot interval when determined manually 
predict$T <- log(as.numeric(predict$T))
predict$ct_brand <- factor(predict$ct_brand)

# Predict using the fitted model
predictions <- predict(
  mlm_theta, 
  newdata = predict,
  level   = 0,                # level = 0 → only fixed effects
  allow.new.levels = TRUE     # allows unseen levels in random effects
)

predictions

#####################

# SCENARIO 2 ADD-ON

# estimating the detection radius for a species at a study site where detection zones for more common species
# are available, but where the practitioner also has a limited amount of measurements for the target species

# Note that this is only for the EDD, not the EDA, because the procedure did not improve for the latter 

#####################

# Load in available data

radial_measurements <- read_csv("example_measurements.csv")

# Predict radius based on available data

# Loop through species-site combinations in the provided data

full.output <- list()

for(k in 1:length(unique(radial_measurements$ID))){
  
  #k=1
  print(100*k/length(unique(radial_measurements$ID)))
  
  input <- radial_measurements[radial_measurements$ID==unique(radial_measurements$ID)[k],]
  input$ndet <- nrow(input)
  
  left = as.numeric(unique(input$left_trunc))
  right = unique(input$right_trunc)
  
  distance <- input$distance
  
  if(!(right%like%"%")){
    right= as.numeric(right)
  }
  
  ## HAZARD
  
  key = "hr"
  mod <- R.utils:::withTimeout(tryCatch(ds(distance, transect="point",key=key,adjustment = "poly",max.adjustments = 2, order=NULL,
                                           truncation=list(left=left,right=right)),
                                        error = function(e){
                                          NA
                                        }),timeout=60)
  
  #print(plot(mod,main=unique(r_df$ID)[k]))
  
  if(!is.na(mod)[1]){
    gof.test <- ddf.gof(mod$ddf)
    pvalue.hr <- gof.test$dsgof$CvM$p
    max_r <- max(mod[["ddf"]][["data"]][["distance"]])
    ndet = nrow(input)
    r.hr <- sqrt(mod$ddf$fitted[1]*mod$ddf$ds$aux$width^2)
    
  }else{
    r.hr = NA
    pvalue.hr=NA
  }
  
  ## HALF-NORMAL
  key = "hn"
  mod <- R.utils:::withTimeout(tryCatch(ds(distance, transect="point",key=key,adjustment = "herm",max.adjustments = 2, order=NULL,
                                           truncation=list(left=left,right=right)),
                                        error = function(e){
                                          NA
                                        }),timeout=60)
  
  #plot(mod,main=unique(r_df$ID)[k])
  if(!is.na(mod)[1]){
    gof.test <- ddf.gof(mod$ddf)
    pvalue.hn <- gof.test$dsgof$CvM$p
    max_r <- max(mod[["ddf"]][["data"]][["distance"]])
    ndet = nrow(input)
    r.hn <- sqrt(mod$ddf$fitted[1]*mod$ddf$ds$aux$width^2)
    
  }else{
    r.hn = NA
    pvalue.hn = NA
  }
  
  # UNIFORM
  key = "unif"
  mod <- R.utils:::withTimeout(tryCatch(ds(distance, transect="point",key=key,adjustment = "cos",max.adjustments = 2, order=NULL,
                                           truncation=list(left=left,right=right)),
                                        error = function(e){
                                          NA
                                        }),timeout=60)
  
  #plot(mod,main=unique(r_df$ID)[k])
  if(!is.na(mod)[1]){
    gof.test <- ddf.gof(mod$ddf)
    pvalue.unif <- gof.test$dsgof$CvM$p
    max_r <- max(mod[["ddf"]][["data"]][["distance"]])
    ndet = nrow(input)
    r.unif <- sqrt(mod$ddf$fitted[1]*mod$ddf$ds$aux$width^2)
    
  }else{
    r.unif = NA
    pvalue.unif=NA
  }
  
  full.output[[k]] <- cbind(input[1,],r.hr,pvalue.hr,r.hn,pvalue.hn,r.unif,pvalue.unif,max_r,ndet=ndet,right=mod$ddf$meta.data$width)
  
  
}

full.output.df <- do.call("rbind.fill",full.output)

# get best performing key
full.output.df$best.key <- colnames(full.output.df[c("pvalue.hr","pvalue.hn","pvalue.unif")])[unlist(apply(full.output.df[,c("pvalue.hr","pvalue.hn","pvalue.unif")], 1, which.max))]
full.output.df$best.key <- gsub("pvalue.","",full.output.df$best.key)

# r per species is r of best key
full.output.df$sample_estimate[full.output.df$best.key=="hr"] <- full.output.df$r.hr[full.output.df$best.key=="hr"]
full.output.df$sample_estimate[full.output.df$best.key=="hn"] <- full.output.df$r.hn[full.output.df$best.key=="hn"]
full.output.df$sample_estimate[full.output.df$best.key=="unif"] <- full.output.df$r.unif[full.output.df$best.key=="unif"]

# Determine weight of measurement data based on sample size

weights <- structure(list(ndet = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 
                        100, 125, 150, 175, 200), sample_weight = c("20", "30", "40", 
                                                                    "50", "50", "60", "60", "60", "70", "70", "70", "80", "80", "80", 
                                                                    "80")), class = "data.frame", row.names = c(NA, -15L))

# Find the weights that correspond to the sample size
full.output.df$sample_weight <- as.numeric(weights$sample_weight[match(full.output.df$ndet,weights$ndet)])/100
full.output.df$model_weight <- 1-full.output.df$sample_weight

# Join model prediction data frame to sample prediction one

full.output.df <- left_join(full.output.df,predict_r[,c("ID","model_estimate")])

# Compute weighted average of model prediction vs. own data prediction
full.output.df$joint_estimate <- full.output.df$sample_estimate*full.output.df$sample_weight +
                                    full.output.df$model_estimate*full.output.df$model_weight




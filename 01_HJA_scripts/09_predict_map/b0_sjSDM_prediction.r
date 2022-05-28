### b3_sjSDM_prediction.r

## Run on ADA
# no markdown format!

## Predictions across landscape (on raster stack of predictors)
## 1. Create full model (from best tuning)
## 2. Predict across new data from raster stack (export data at regular grid of xy coords)
## convert to rasters and stack up
## 3. Do plots


#### Read data on Ada  #####

options(echo=TRUE) # if you want see commands in output file

# CHANGE HERE FOR PATH TO REQUIRED PYTHON VERSION, if necessary
# Sys.setenv(RETICULATE_PYTHON="/PATH/TO/r-sjSDM/bin/python")

# myPaths <- .libPaths()
# myPaths <- c("/gpfs/scratch/hsp20azu/newrlib",myPaths[2],myPaths[1])
# .libPaths(myPaths)

library("sjSDM")
packageVersion("sjSDM")
sjsdmV <- packageVersion('sjSDM')		# package version

library(dplyr)
library(glue)
library(here)


varsName <- 'vars11'
date.model.run <- '202204'
abund <- "pa"
k <- 5
minocc <- 6; period <- "S1"
nstep <- 1000
	
resFolder <- here('04_Output', "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
plotFolder <- here('04_Output', "prediction_map")
sppdatapath <- here('03_format_data','otu')
	
dir.exists(resFolder)
dir.exists(plotFolder)
dir.exists(sppdatapath)
dir(resFolder)

	
## load model data 
load(here(sppdatapath, glue('forbestm_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
# load(file.path(sppdatapath, paste0("modelData_",abund,".rdata")))
## only sitenames in this set are M1 S1

## Load new data for prediction
load(file.path(plotFolder, 'rdata',  "newData_scaled_clamp.rdata")) #newData_clamp_wide.sc, xy.sites.sc, newXY.sc, allVars.sc,
## too big for github - data is on ADA and locally: save in local folder
## "processed_gis_data/r_oversize/newData_scaled.rdata"

## update device - can't load new data onto GPU - too large?? to many copies??
device <- "cpu"

## update for testing
# iter= 10
# sampling = 100

# formula.env = 'envDNN'
hidden <- list(c(50L,50L,10L), c(25L,25L,10L))

## get best tune run
res <- read.csv(file.path(resFolder, 
                             paste0("manual_tuning_sjsdm_",sjsdmV,
                                    "_", varsName, "_", k, "CV_", period, "_meanEVAL_", abund, "_min",
                                    minocc,"_nSteps", nstep, ".csv")))
head(res)
res.best <- res[which.max(res$AUC.test_mean),,drop = T]
res.best
rm(res)

# Choose pa or qp reponse data and family
if(abund == "pa") {
  Y <- otu.pa.csv
  family <- stats::binomial('probit') } else {
    if(abund ==  "qp") {
      Y <- otu.qp.csv
      family <- stats::binomial('probit') # check other family?
    } else stop("check abund")
  } 


## Need to run model again on cpu -- predict doesn't work with gpu model - too large - new data

# select X data - form globally scaled data (with newData), no validation - this is final best model
## Also filter S1 as in tuning model above (in env.vars data set)
# scale.env <- allVars.sc %>%
#   filter(SiteName %in% train.Names, period == "S1") %>% ## is the validation set (75% - select.percent)
#   dplyr::select(all_of(vars))

scale.env <- allVars.sc[allVars.sc$SiteName %in% train.Names & allVars.sc$period == "S1",vars]

head(scale.env)

any(sapply(scale.env, function(x) any(is.na(x)))) # should be FALSE - no NAs

## filter new data for prediction to same vars
newData_clamp_wide.sc <- newData_clamp_wide.sc %>%
  dplyr::select(all_of(vars))


# spatial data - just sites in model data, as above
# scale.XY <- xy.sites.sc %>%
#   filter(SiteName %in% train.Names, period == "S1") %>%
#   dplyr::select(c("UTM_E", "UTM_N"))%>%
#   as.matrix()

scale.XY <- as.matrix(xy.sites.sc[xy.sites.sc$SiteName %in% train.Names & xy.sites.sc$period == "S1", c("UTM_E", "UTM_N")])

# do model - spatial
modFull_sp <- sjSDM(
  
  Y = as.matrix(Y),
  
  env = DNN(data = scale.env, formula = ~.,
            hidden=hidden[[res.best$hidden.ind]],
            lambda = res.best$lambda.env,
            alpha = res.best$alpha.env,
            activation = res.best$acti.sp,
            dropout=res.best$drop,
            bias=T),
  
  biotic = bioticStruct(lambda=res.best$lambda.bio,
                        alpha=res.best$alpha.bio, 
                        on_diag=F, inverse = FALSE),
  
  spatial = linear(data=scale.XY, ~0+UTM_E*UTM_N,
                   lambda=res.best$lambda.sp,
                   alpha=res.best$alpha.sp),
  
  learning_rate = res.best$lr,
  iter = iter, 
  family = family, 
  sampling = sampling, # 150L, 5000L
  device = device
)

## Save model
saveRDS(modFull_sp, file.path(resFolder, paste0("s-jSDM_fullMod_sp_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rds")))

gc()

### Do prediction
# pred.sp = predict(modFull_sp, newdata = newData.sc, SP = newXY.sc)

## predict 5 times and save as list
pred_all <- lapply(1:5, function(i) {
  predict(modFull_sp, newdata = newData.sc, SP = newXY.sc)}
)

# get mean prediction and sd
pred.mn = apply(abind::abind(pred_all,along = -1L), 2:3, mean)
pred.sd = apply(abind::abind(pred_all,along = -1L), 2:3, sd)

save(pred.mn, pred.sd, file = file.path(plotFolder, 'rdata', paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata")))

### DO clamp prediction ######

## predict 5 times and save as list
pred_all.cl <- lapply(1:5, function(i) {
  predict(modFull_sp, newdata = newData_clamp_wide.sc, SP = newXY.sc)}
)

# get mean prediction and sd
pred.mn.cl = apply(abind::abind(pred_all.cl, along = -1L), 2:3, mean)
pred.sd.cl = apply(abind::abind(pred_all.cl, along = -1L), 2:3, sd)

save(pred.mn.cl, pred.sd.cl, file = file.path(plotFolder, 'rdata', paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, "_clamp", ".rdata")))



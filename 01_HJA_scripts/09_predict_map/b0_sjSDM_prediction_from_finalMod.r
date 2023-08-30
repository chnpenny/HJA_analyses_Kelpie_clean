### b3_sjSDM_prediction.r

## Run on Keele Armstrong
# no markdown format!, no glue, no here.

## Predictions across landscape (on raster stack of predictors)
## 1. Create full model (from best tuning)
## 2. Predict across new data from raster stack (export data at regular grid of xy coords)
## convert to rasters and stack up
## 3. Do plots

getwd()

#### Read data on Keele cluster  #####

options(echo=TRUE) # if you want see commands in output file

# CHANGE HERE FOR PATH TO REQUIRED PYTHON VERSION, if necessary
# Sys.setenv(RETICULATE_PYTHON="/PATH/TO/r-sjSDM/bin/python")
## OR set preferred conda env for reticulate
reticulate::use_condaenv(condaenv = "/home/ggc34/.conda/envs/r-sjsdm")

library("sjSDM")
packageVersion("sjSDM")
#sjsdmV <- packageVersion('sjSDM')		# package version
sjsdmV <- "1.0.5"

library(dplyr)
library(abind)

varsName <- 'vars11'
date.model.run <- '2023'
abund <- "pa"
k <- 5
minocc <- 6; period <- "S1"
nstep <- 1000

resFolder <- file.path("./04_Output", "sjsdm_general_outputs", paste0(varsName, "_", date.model.run))
plotFolder <- file.path("./04_Output", "prediction_map")

## Make sure folders exist
if(!all(dir.exists(resFolder),
       dir.exists(plotFolder),
       dir.exists(file.path(plotFolder, "rdata")))) stop("Check folders")

dir(resFolder)



## load final model
## final model is s-jSDM_tuned.model_S1_pa_5CV_min6_vars11_AUC_2023.RDS
## in this folder: "04_Output/sjsdm_general_outputs/vars11_2023/
modFull_sp <- readRDS(file.path(resFolder, "s-jSDM_tuned.model_S1_pa_5CV_min6_vars11_AUC_2023.RDS"))

## Load new data for prediction
load(file.path(plotFolder, 'rdata',  "newData_scaled_clamp.rdata")) 
#newData_clamp_wide.sc, xy.sites.sc, newXY.sc, allVars.sc,

## prediction grid
head(newData_clamp_wide.sc)
dim(newData_clamp_wide.sc)

## Check env vars used in model are present in prediction grid and filter to these
colnames(modFull_sp$settings$env$data)
length(colnames(modFull_sp$settings$env$data))
sum(!colnames(modFull_sp$settings$env$data) %in% names(newData_clamp_wide.sc))

## filter new data for prediction to same vars
newData_clamp_wide.sc <- newData_clamp_wide.sc %>%
 dplyr::select(all_of(colnames(modFull_sp$settings$env$data)))


### Do prediction

## predict 5 times and save as list
pred_all.cl <- lapply(1:5, function(i) {
  predict(modFull_sp, newdata = newData_clamp_wide.sc, SP = newXY.sc)}
)
 
## get mean prediction and sd
pred.mn.cl = apply(abind::abind(pred_all.cl, along = -1L), 2:3, mean)
pred.sd.cl = apply(abind::abind(pred_all.cl, along = -1L), 2:3, sd)

## Save offline -- too large for github
save(pred.mn.cl, pred.sd.cl, file = file.path(plotFolder, 'rdata', paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, "_clamp", ".rdata")))

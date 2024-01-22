#### TSNE analysis

# ..... setup

rm(list=ls())

# HJA_analyses_Kelpie_clean # is the root folder and must have a .Rproj file in it for here::here() to work.
# setwd() # set here to HJA_analyses_Kelpie_clean or use here::here()

	
library(dplyr)
library(here)
library(glue)
library(Rtsne)
library(terra)


# ..... set-names
abund = "pa"
date.model.run = '2024'
varsName = 'vars11'
minocc = 6; period = "S1"

gis_in <- here::here('03_format_data','gis','raw_gis_data')
	
outputpath = here('04_Output')
	
datapath = here('01_HJA_scripts','08_sjsdm')
gis_out = here('03_format_data','gis')
	
modFolder = file.path(outputpath, "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
resFolder = file.path(outputpath, "sjsdm_prediction_outputs", glue('{varsName}_{date.model.run}'))
plotFolder = file.path(outputpath, "prediction_map")

if(!dir.exists(resFolder)) dir.create(resFolder)
dir.exists(modFolder); dir.exists(plotFolder)

# ..... load-data

## load species AUC results for filtering
load(file.path(modFolder, "spp_test_data.rdata")) #auc_by_spp, rsq_final, tune.results

## Mean AUC per species (and other eval metrics) from 5CV
str(auc_by_spp)

## Filter species by auc
auc.filt = 0.70
# threshold for presence absence data
# tr <- 0.5
	
# how many species after AUC filter?
sum(auc_by_spp$mean > auc.filt, na.rm = T)
	
# incidence 
head(auc_by_spp$incidence)
	
# load clamp predictions (working folder - too large for github)
load(file.path("./working", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", date.model.run, "_", abund, "_clamp", ".rdata")))
# pred.mn.cl, pred.sd.cl

str(pred.mn.cl, 1)
dim(pred.mn.cl)


## filter for species performance
pred.in.cl = pred.mn.cl[,auc_by_spp$mean >= auc.filt & !is.na(auc_by_spp$mean)]
dim(pred.in.cl)
	
## load raster templates
load(file.path(gis_out, "templateRaster.rdata")) ## r.msk, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting

# convert to (new) terra package format
r.msk
r.msk <- terra::rast(r.msk)
r.msk

r.aoi.pred <- terra::rast(r.aoi.pred)


## clamped grid - load each species predicted distribution
rList <- lapply(data.frame(pred.in.cl), function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- x
  tmp
  
})
# plot(tmp)
rStack.cl = terra::rast(rList)
rStack.cl

plot(rStack.cl, 1)

## Save stacked species rasters in terra format - wrapped - off github
rstack_w <- terra::wrap(rStack.cl)
r.msk_w <- terra::wrap(r.msk)
r.aoi_w <- terra::wrap(r.aoi.pred)
save(rstack_w, r.msk_w, r.aoi_w, indNA, 
     file = file.path("./working", paste0("spp_rast_cl_", varsName, "_", date.model.run, ".rdata")))


# ..... TSNE
## Full data set
Xmat <- pred.in.cl
#r <- rast(rStack.cl) # make template

# pa version
dim(Xmat)
Xmat[1:10, 1:10]
perplexity = 50			#
	
# Max
(nrow(Xmat) - 1)/3
tsne = Rtsne::Rtsne(Xmat, dims = 2, perplexity = perplexity, theta = 0.5, pca = FALSE, num_threads = 0) # can set parallel options if using openMP
	
# plot(tsne$Y, asp = 1, pch = ".")
# str(tsne, max.level =1)
# plot(tsne$Y, asp = 1)


## put site scores into raster
makeR <- function(r, siteScores, NAs) {
  
  rSites <- rast(r)
  rSites[] <- NA
  rSites[NAs] <- siteScores
  rSites
  
}
	
rSites1 <- makeR(r.msk, tsne$Y[,1], indNA)
rSites2 <- makeR(r.msk, tsne$Y[,2], indNA)

names(rSites1) <- "TSNE1"
names(rSites2) <- "TSNE2"

plot(c(rSites1, rSites2))

# wrap and save
tsne_rast <- terra::wrap(c(rSites1, rSites2))

save(tsne, r.msk, tsne_rast, indNA, file = file.path(resFolder, "ord_tsne_res_cl_p50.rdata")) # with perp = 50

terra::writeRaster(rSites1, 
                   filename = file.path(plotFolder, 'rdata', paste0("tsne1_nopca_cl_p50_", varsName, "_", date.model.run, ".tif")), 
                   datatype = "FLT4S", overwrite = T)
terra::writeRaster(rSites2, 
                   filename = file.path(plotFolder, 'rdata', paste0("tsne2_nopca_cl_p50_", varsName, "_", date.model.run, ".tif")),
                   datatype = "FLT4S", overwrite = T)
# 
	
# pdf(file.path(plotFolder, 'plot', "tsne_scatter_cl_p50.pdf"))
#plot(tsne$Y, pch = ".")
#dev.off()
	

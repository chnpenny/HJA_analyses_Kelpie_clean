
## Make raster from predictions

rm(list=ls())

library(dplyr)
library(terra)
library(sf)

utm10N <- 32610

gis <- here::here("03_format_data","gis")
gis_in <- here::here('03_format_data','gis','raw_gis_data')
gis_out <- here::here('03_format_data','gis','processed_gis_data')
dir(gis_in)

minocc = 6; period = "S1"
varsName = 'vars11'
date.model.run = '2024'
abund = "pa"

resFolder = here::here('04_Output', "sjsdm_general_outputs", glue::glue('{varsName}_{date.model.run}'))
predFolder = here::here('04_Output', "sjsdm_prediction_outputs", glue::glue('{varsName}_{date.model.run}'))
plotFolder = here::here('04_Output', "prediction_map")
dir(resFolder)


## load species results data (with names)
load(file.path(resFolder, "spp_test_data_names.rdata")) # spp, spp.in

# get results
head(spp.in)
sum(spp.in$gt07)
mean(spp.in$mean) # 0.7954
range(spp.in$mean) # 0.7003 - 0.9938

## load stacked OTU rasters
load(file.path("working", paste0("spp_rast_cl_", varsName, "_", date.model.run, ".rdata")))

## un wrap 
r.msk <- terra::rast(r.msk_w)
rStack.cl <- terra::rast(rstack_w)

rm(r.msk_w, rstack_w)

## add auc incidence names to stack
names(rStack.cl) <- paste0(spp.in$best.name, "_", "auc=", round(spp.in$mean, 2), "_",  "prev=", round(spp.in$incidence,2))

# folder for sp tifs
sppFolder <- file.path(predFolder, "spp_tifs_cl")
if(!dir.exists(sppFolder)) dir.create(sppFolder)
filenames <- file.path(sppFolder, paste0(names(rStack.cl), "_spp_cl.tif"))

terra::writeRaster(rStack.cl, filename = filenames, datatype = "FLT4S")

# threshold for binary maps for species richness
tr <- 0.5
rStack.bin.cl <- terra::classify(rStack.cl, 
                                 rcl = matrix(c(0, tr, 0, tr, 1, 1), nrow = 2, ncol = 3, byrow = TRUE))
rStack.sum.cl <- sum(rStack.cl)
names(rStack.sum.cl) <- "sp sum"
spRich.cl <- sum(rStack.bin.cl)
names(spRich.cl) <- "sp richness"

plot(spRich.cl)
plot(rStack.sum.cl)

## Do maps of richness by groups.

## load extra layers for plotting
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)

st_area(hja.utm) / 1000000 ## in km^2

## Make species richness stack
rStack.bin.cl
spRich_order.cl <- terra::tapp(rStack.bin.cl, spp.in$order, fun = sum)
names(spRich_order.cl)

terra::writeRaster(rStack.sum.cl, 
                   filename = file.path(predFolder, "spSum_cl.tif"), 
                   datatype = "FLT4S", overwrite = T)
# "./04_Output/sjsdm_prediction_outputs/vars11_2024/spSum_cl.tif"

terra::writeRaster(spRich.cl, 
                   filename = file.path(predFolder, "spRich_all_cl.tif"), 
                   datatype = "FLT4S", overwrite = T)


# convert HJA to raster and save
# hja.r <- rStack.cl
# 
# hja.utm
# plot(r.aoi.pred, colNA = "black")
# hja.r <- rasterize(hja.utm, r.aoi.pred)
# 
# plot(hja.r, colNA = "black")
# 
# # convert to 1 and 0 for inside/outside HJA
# 
# hja.r[is.na(hja.r)] <- 0
# hja.r <- mask(hja.r, r.aoi.pred)
# plot(hja.r, colNA = "black")
# 
# save(hja.r, file = file.path(gis_out, "hja_raster.rdata"))


## write single pdf with all spp
source(file.path("01_HJA_scripts/09_predict_map/source", "plotStack.r"))

# function to add to each plot
addAll <- function(){

  plot(st_geometry(hja.utm), add = T, col = NA, border = "black")
  #plot(st_geometry(aoi.pred.sf), add = T, col = NA, border = "black")
  plot(st_geometry(xy.utm), add = T, col = "grey40", pch = 3, cex = 0.2)

}


## write single pdf with all spp
pdf(file.path(predFolder, paste0("all_spp_", varsName, "_", date.model.run, ".pdf")), width = 7, height = 7)
# this is Fig. S-individual SDMs in the supplement
plotStack(rStack.cl, fun = addAll, cex.main = 0.7)
dev.off()

## file.path("04_Output/figures", "all_spp.pdf") # old path

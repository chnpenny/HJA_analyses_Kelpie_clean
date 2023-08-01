


gis_in <- here::here('03_format_data','gis','raw_gis_data')
gis_out = here('03_format_data','gis')


load(file.path(resFolder, "ord_tsne_res_cl_p50.rdata")) # tsne, r, rSites1, rSites2, NAs,

# load raster as brick
allBrck <- brick(file.path(gis_out, "processed_gis_data/r_utm", "allStack_aoi.tif"))
load(file.path(gis_out, "processed_gis_data", "brNames.rdata")) # load brick names

# get names and name groups
names(allBrck) <- brNames
allBrck
names(allBrck)

## Load old growth 
ogsi <- raster(file.path(gis_in, "ogsi_2012_smoothed.tif"))
ogsi.utm <- projectRaster(ogsi, to = rSites1)
names(ogsi.utm) <- "ogsi"
plot(ogsi.utm)
ogsi.utm <- mask(ogsi.utm, r.msk)


## make strata by elevation
elev20 <- cut(predStack$be30, breaks = 5)
names(elev20) <- "elev20"
plot(elev20)

# reduce to covariates of interest
predStack <- stack(allBrck[[c("be30", "insideHJA","cut_msk", "DistRoad","DistStream")]], ogsi.utm, rSites1, rSites2, elev20)

# predStack <- crop(predStack, beta_irr)
# predStack <- addLayer(predStack, beta_irr)
# predStack <- mask(predStack, beta_irr)


## TSNE ~ Elevation * OGSI

load(file.path(gis_out, "processed_gis_data", "hja_raster.rdata")) ## hja.r
names(hja.r) <- "HJA"

library(terra)
set.seed(99)
pts <- terra::spatSample(rast(elev20), size = 400, method = "stratified", xy = TRUE)
head(pts)

plot(elev20)
points(pts$x, pts$y, pch = 16, col = pts$layer)

pts_extr <- raster::extract(predStack, pts[,c("x", "y")], df = TRUE)
head(pts_extr)
pts_extr$elev20 <- factor(pts_extr$elev20)

library(ggplot2)

ggplot(pts_extr, aes(x = ogsi, y = TSNE1, col = elev20))+
  geom_point()+
  scale_color_viridis_d()+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))+
  xlab("Old Growth Index")
ggsave(filename = file.path("04_Output/figures", "ogsi_tsne1_x_elev20.png"), units = "mm", height = 200, width = 200)


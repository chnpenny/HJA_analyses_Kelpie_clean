#### Irreplaceability, Ordination, OSGI plots

library(sf)
library(raster)
library(terra)
library(ggplot2)
library(tmap)
library(patchwork)

# minocc = 6; period = "S1"
varsName = 'vars11'
date.model.run = '2023'
# abund = "pa"

utm10N <- 32610

gis_in <- here::here('03_format_data','gis','raw_gis_data')
gis_out = here::here('03_format_data','gis')


resFolder = here::here('04_Output', "sjsdm_prediction_outputs", glue::glue('{varsName}_{date.model.run}'))

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)


## load TSNE data
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
# mask ogsi
ogsi.utm <- mask(ogsi.utm, allBrck$insideHJA)

## make strata by elevation
elev20 <- cut(allBrck$be30, breaks = 5)
names(elev20) <- "elev20"
plot(elev20)

# load beta irreplaceability
beta_irr <- raster(file.path(gis_out, "processed_gis_data", "r_utm", "beta_r_prob_noagg.tif"))
names(beta_irr) <- "beta"

## load spp richness
rStack.sum.cl <- raster(file.path(gis_out, "processed_gis_data/r_utm", "spSum_cl.tif"))
spRich.cl <- raster(file.path(gis_out, "processed_gis_data/r_utm", "spRich_all_cl.tif"))

rStack.sum.cl
names(spRich.cl) <- "sp.richness"

# reduce to covariates of interest
predStack <- stack(allBrck[[c("be30", "insideHJA","cut_msk", "DistRoad","DistStream")]], 
                   ogsi.utm, rSites1, rSites2, elev20, beta_irr, spRich.cl, rStack.sum.cl)

predStack

plot(predStack)
plot(predStack[[c("ogsi", "TSNE1", "TSNE2")]])
summary(values(predStack$ogsi))

# convert to terra raster
predStack_t <- rast(predStack)

ogsi_33 <- terra::focal(predStack_t$ogsi, w = 33, fun = "mean", expand = TRUE, fillValue = NA, na.rm = TRUE)
names(ogsi_33) <- "ogsi_33"

ogsi_33 <- mask(ogsi_33, rast(allBrck$insideHJA))
predStack_t <- c(predStack_t, ogsi_33)
names(predStack_t)

plot(ogsi_33)

plot(predStack_t[[c("ogsi", "ogsi_33")]])

load(file.path(gis_out, "processed_gis_data", "hja_raster.rdata")) ## hja.r
names(hja.r) <- "HJA"

set.seed(99)
pts <- terra::spatSample(predStack_t$elev20, size = 2000, method = "stratified", xy = TRUE)
head(pts)

plot(elev20)
#points(pts$x, pts$y, pch = 16, col = pts$layer)

pts_extr <- terra::extract(predStack_t, pts[,c("x", "y")], raw = FALSE)
head(pts_extr)
pts_extr$elev20 <- factor(pts_extr$elev20)
pts_extr$insideHJA <- factor(pts_extr$insideHJA, levels = c(0,1), labels = c("Outside", "Inside"))

summary(pts_extr)

## get elevation of bands
e_bands <- aggregate(pts_extr$be30, list(pts_extr$elev20), range)
str(e_bands)
(e_bands$x %/% 5) * 5

### Figure 3
cols <- terrain.colors(100, rev = TRUE) # palette for ogsi map
txt1 <- 12
txt2 <- 16

## OSGI map - panel A
tm1 <- tm_shape(predStack_t$ogsi_33)+
  tm_raster(palette = cols, title = "Old-growth \nstructural index", style = "cont", breaks = seq(0,80,20))+
  tm_shape(predStack_t$ogsi_33)+
  tm_raster(palette = grey.colors(50), alpha = 0.4, legend.show = FALSE)+
  tm_shape(hja.utm)+
  tm_borders(col = "black", lwd = 2)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"), text.size = 1)+
  tm_compass(position = c("right", "top"), size = 2)+
  tm_layout(legend.title.size = 1.5, legend.text.size = 1, frame = FALSE,
            legend.position = c("LEFT", "TOP"))

tm1
p1 <- tmap::tmap_grob(tm1) # convert to grob for patchwork

## OGSI vs sp richness by HJA - panel B
p2 <- ggplot(pts_extr, aes(x = ogsi_33, y = sp.richness, col = insideHJA))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.35)+ # 
  scale_color_manual(values = c("#f1a340", "#998ec3"), name = "HJA boundary")+
  xlab("Mean old-growth structural index within 1 km")+
  ylab("Species richness")+
  xlim(c(0,75))+
  theme(axis.text = element_text(size = txt1),
        axis.title = element_text(size = txt2),
        legend.text = element_text(size = txt1),
        legend.position='bottom')
p2

## OGSI vs TSNE by elev - Panel C
p3 <- ggplot(pts_extr, aes(x = ogsi_33, y = TSNE2, col = elev20))+
  geom_point(alpha = 0.3, size = 1, show.legend = FALSE)+
  #scale_color_manual(values = pals::kovesi.diverging_linear_bjr_30_55_c53(5), name = "Elevation bands")+
  scale_color_manual(values = pals::kovesi.diverging_bky_60_10_c30(5), name = "Elevation bands")+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.4, show.legend = FALSE)+
  xlab("Mean old-growth structural index within 1 km")+
  ylab("TSNE axis 2")+
  theme(axis.text = element_text(size = txt1),
        axis.title = element_text(size = txt2),
        legend.text = element_text(size = txt1),
        legend.position='bottom')

# get smooth x and y coords for labels
df3.smooth <- 
  ggplot_build(p3)$data[[2]] %>% 
  dplyr::group_by(group) %>%
  dplyr::filter(x == min(x))


p3 <- p3 + ggrepel::geom_label_repel(data = df3.smooth, aes(x = x, y = y, col = as.factor(group), label = group),
                                     nudge_x = -1, fontface = "bold", alpha = 0.8, show.legend = FALSE)
p3


## OGSI vs IRR by elev - Panel D
p4 <- ggplot(pts_extr, aes(x = ogsi_33, y = beta*100, col = elev20))+
  geom_point(alpha = 0.3, size = 1, show.legend = FALSE)+
  scale_color_manual(values = pals::kovesi.diverging_bky_60_10_c30(5), name = "Elevation bands")+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.4, show.legend = FALSE)+ # 
  xlab("Mean old-growth structural index within 1 km")+
  ylab("Irreplaceability index")+
  theme(axis.text = element_text(size = txt1),
        axis.title = element_text(size = txt2),
        legend.text = element_text(size = txt1),
        legend.position='bottom')

# get smooth x and y coords for labels
df4.smooth <- 
  ggplot_build(p4)$data[[2]] %>% 
  dplyr::group_by(group) %>%
  dplyr::filter(x == min(x))

p4 <- p4 + ggrepel::geom_label_repel(data = df4.smooth, aes(x = x, y = y, col = as.factor(group), label = group),
            nudge_x = -1, fontface = "bold", alpha = 0.8, show.legend = FALSE)

p4

## Put plots together
pAll <- wrap_plots(p1, p2, p3, p4) + 
  plot_annotation(tag_levels = "A")

## Export and check
ggsave(plot = pAll, filename = file.path("04_Output/figures", "Figure3.png"), units = "mm", height = 400, width = 400)
shell.exec(file.path(getwd(), "04_Output/figures", "Figure3.png"))
shell.exec(file.path(getwd(), "04_Output/figures"))

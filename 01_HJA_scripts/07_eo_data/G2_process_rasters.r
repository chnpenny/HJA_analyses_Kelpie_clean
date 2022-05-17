
#### Process and standardise all predictor rasters #####

## 1. Rasters standardised to UTM 10N (wgs), at 30 m resolution.
## 2. some rasters summarised over 100, 250, 500 or 1k radius circular windows
## 3. Distance to roads, streams, insideHJA rasterized.
## All rasters exported as tif files and rdata

## RUN ON ADA (UEA HPC cluster)
# so not markdown format

options(echo=TRUE) # if you want see commands in output file
getwd() #
# HJA_analyses_Kelpie_clean # is the root (wd) folder and must have a .Rproj file in it for here() to work.

library(raster)
library(sf)
library(exactextractr) # polygon extraction and summary for rasters
library(here)	


# wgs84 UTM 10N
utm10N = 32610
prj4.utm10 = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs" #same as above for raster (but now also accepts epsg)
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 = 26910
	
#testing local --  change this to github GIS folder
gis_in = here('03_format_data','gis',"raw_gis_data") 
gis_out = here('03_format_data','gis',"processed_gis_data")  
	
## 0. Get prelim data ######
## common extent in utm10N
load(file.path(gis_out, "commonExtent.rdata")) # ext
# make template rasters at these resolutions
r1 = raster(ext, res = c(1,1), crs = prj4.utm10)
r10 = raster(ext, res = c(10,10), crs = prj4.utm10)
r30 = raster(ext, res = c(30,30), crs = prj4.utm10)
	
### 1. Canopy cover ####

# Create canopy cover (from Lidar) - subtract highest hit from bare earth.
# Create canopy gap/cover binary layer (<4 m is a canopy gap), summarise proportion of canopy cover over 250 and 500 m

# download files from online storage
download.file("https://od.lk/d/OTJfMzAwOTQzNjRf/latlong_highesthit.tif", 
              destfile = file.path(gis_in, "lidar/latlong_highesthit.tif"),
              mode = "wb") # for windows

download.file("https://od.lk/d/OTJfMzAwOTM4MjFf/latlong_bare_earth.tif", 
              destfile = file.path(gis_in, "lidar/latlong_bare_earth.tif"),
              mode = "wb") # for windows

hh = raster(file.path(gis_in, "lidar/latlong_highesthit.tif")) # from Lidar data. Provided by Oregon State University
be = raster(file.path(gis_in, "lidar/latlong_bare_earth.tif"))
# crs: +proj=lcc +lat_0=41.75 +lon_0=-120.5 +lat_1=43 +lat_2=45.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=ft +no_defs
# at 3.014 ft resolution approx 1m
	
## set up raster parallel
# # set number of cores...
n = parallel::detectCores() - 1
	
beginCluster(n = n)
	
# convert to m
f1 = function(x) x * 0.30480
	
# change to meters
hh_m = clusterR(hh, calc, args=list(fun=f1))
be_m = clusterR(be, calc, args=list(fun=f1))
	

## canopy height
s = stack(hh_m, be_m)
f2 = function(h, b) h - b
ht_m = clusterR(s, overlay, args = list(fun = f2))
	

# set values below 0 to 0
ht_m = clusterR(ht_m, reclassify, args=list(rcl = c(-Inf, 0, 0)))
# ht_m[ht_m < 0] <- 0
ht_m
# plot(ht_m)
	
## create gap layer
ht_gt4m = clusterR(ht_m, reclassify, args=list(rcl = c(0, 4, 0, 4,Inf,1)))
# ht_gt4m <- ht_m > 4
ht_gt4m
# plot(ht_gt4m)
	

# project canopy cover/gap layer
ht_gt4_utm = projectRaster(ht_gt4m, r1, method = "ngb", filename = file.path(gis_out, "r_utm/ht_gt4_utm.tif"), 
                            datatype = "INT1U", overwrite = TRUE)
# project canopy height layer
ht_utm = projectRaster(ht_m, r1, method = "bilinear", filename = file.path(gis_out, "r_utm/ht_utm.tif"), 
                        datatype = "FLT4S", overwrite = TRUE)
ht_utm
	
be_utm = projectRaster(be_m, r1, method = "bilinear", filename = file.path(gis_out, "r_utm/be_utm.tif"), 
                        datatype = "FLT4S", overwrite = TRUE)
	
endCluster()
	

ht_gt4_utm
#plot(ht_gt4_utm)
	

## get 10m and 30m resolution raster with mean of 1m raster for canopy cover and canopy height
ht_gt4_r30 = aggregate(ht_gt4_utm, fact = 30, fun = mean, expand = FALSE,
                        filename = file.path(gis_out, "r_utm/ht_gt4_r30.tif"), datatype = "FLT4S", overwrite = TRUE)
	
# plot(ht_gt4_r30)
ht_r30 <- aggregate(ht_utm, fact = 30, fun = mean, expand = FALSE,
                    filename = file.path(gis_out, "r_utm/ht_r30.tif"), datatype = "FLT4S", overwrite = TRUE)

ht_r10 <- aggregate(ht_utm, fact = 10, fun = mean, expand = FALSE,
                    filename = file.path(gis_out, "r_utm/ht_r10.tif"), datatype = "FLT4S", overwrite = TRUE)


be_r30 <- aggregate(be_utm, fact = 30, fun = mean, expand = FALSE,
                    filename = file.path(gis_out, "r_utm/be_r30.tif"), datatype = "FLT4S", overwrite = TRUE)

be_r10 <- aggregate(be_utm, fact = 10, fun = mean, expand = FALSE,
                    filename = file.path(gis_out, "r_utm/be_r10.tif"), datatype = "FLT4S", overwrite = TRUE)
	

# Bring back from saved (if in new session)
# be_r30 <- raster(file.path(gis_out, "r_utm/be_r30.tif"))
# ht_r30 <- raster(file.path(gis_out, "r_utm/ht_r30.tif"))
# ht_gt4_r30 <- raster(file.path(gis_out, "r_utm/ht_gt4_r30.tif"))

# do focal windows of 250 and 500 m radius - but based on 30m resolution template. ie for each centroid of the
## 30m raster cells, summarise (mean) the canopy cover raster over a buffer of 250 and 500 m.


# get coordinates at 30m resolution grid
xy_r30 = raster::coordinates(r30)

# extract with polygons
head(data.frame(xy_r30, id = 1:nrow(xy_r30)))
# 
# create point sf object
pts_r30 = st_as_sf(data.frame(xy_r30, id = 1:nrow(xy_r30)),coords = c("x", "y"), crs = utm10N)
# buffer points to create polygons
pts_r30_b250 = st_buffer(pts_r30, dist = 250)
pts_r30_b250
pts_r30_b500 = st_buffer(pts_r30, dist = 500)

# subset  for testing - see time
# pts_r30_b250 <- pts_r30_b250[50000:55000,]
#plot(ht_gt4_utm)
# plot(pts_r30_b250[1,], add = T)

## Use exactextractr package - much faster than raster::extract
start <- Sys.time()
vals_gt4_r250 <- exactextractr::exact_extract(ht_gt4_utm, pts_r30_b250, fun = "mean")
end <- Sys.time()
end - start  #3.43058 hours

head(vals_gt4_r250)

start <- Sys.time()
vals_gt4_r500 <- exactextractr::exact_extract(ht_gt4_utm, pts_r30_b500, fun = "mean")
end <- Sys.time()
end - start  # 13.21485 hours

# save(vals_gt4_r250, file = file.path(gis, "vals_gt4_r250.rdata"))
# save(vals_gt4_r500, file = file.path(gis, "vals_gt4_r500.rdata"))

## put back into rasters
gt4_r250 <- r30
values(gt4_r250) <- vals_gt4_r250
writeRaster(gt4_r250, filename = file.path(gis_out, "r_utm/gt4_r30_r250.tif"), datatype = "FLT4S", overwrite = TRUE)
# gt4_r250 <- raster(file.path(gis_out, "r_utm/gt4_r30_r250.tif"))
# plot(gt4_r250)

gt4_r500 <- r30
values(gt4_r500) <- vals_gt4_r500
writeRaster(gt4_r500, filename = file.path(gis_out, "r_utm/gt4_r30_r500.tif"), datatype = "FLT4S", overwrite = TRUE)
# gt4_r500 <- raster(file.path(gis_out, "r_utm/gt4_r30_r500.tif"))

## 30 m rasters already exported. 
## gather as stack

be_ht <- stack(ht_r30, ht_gt4_r30, gt4_r250, gt4_r500)
names(be_ht)
be_ht.names <- c("ht30", "gt4_r30", "gt4_250", "gt4_500")
names(be_ht) <- be_ht.names

be_ht[[1]]

save(be_ht, be_ht.names, file = file.path(gis_out, "be_ht.rdata"))

writeRaster(be_ht, bylayer = T, 
            filename = file.path(gis_out, "r_utm/be_ht.tif"), suffix = "names", overwrite = TRUE)

# plot(stack(gt4_r250, gt4_r500))

## 2 Cut (logged) Areas ######

## bring in cut areas - shapefile - from Oregon SU
cut <- st_read(file.path(gis_in, "shape/disturbance.shp"))
cut
cut[order(cut$YEAR, na.last = F),] # OJO!!! 0 in YEAR column in arcgis come in as NAs in R/.
sum(is.na(cut$YEAR))
cut.utm <- st_transform(cut, crs = utm10N)
cut.utm
# plot(cut.utm[, c("YEAR")])
# plot(xy.utm, add = T, pch = 2, col = "black", cex = 1.5)
rm(cut)

## do proportion of cut within 1 km
hist(cut.utm$YEAR)
# 0 values are NA, and outside study area
# convert to raster and do focal
extent(cut.utm)
ext
cut.r <- rasterize(cut.utm, r30, field = "YEAR", fun = "max")
cut.r
plot(cut.r)

# make a binary cut layer - areas that have or not been logged
cut.msk <- cut.r
cut.msk
cut.msk[cut.msk > 1] <- 1
plot(cut.msk)

cut.40msk <- reclassify(cut.r, rcl = c(0,1981, 0, 1981, Inf, 1))

# convert NA to 0, as this is 0 disturbance, only important where no disturbance in whole focal window
cut.msk[is.na(cut.msk)] <- 0
cut.40msk[is.na(cut.40msk)] <- 0

plot(stack(cut.msk, cut.40msk), colNA = "black")

# give focal window values of 1, and use mean instead of giving fractional values adding up to 1, and using sum.
# Use pad = TRUE, and na.rm = T, to manage edge effects and take care of any NAs in the rasters.
## This is slower than with the default sum and fractional values.

r1k <- ifelse(focalWeight(r30, d = 1000, type = "circle")>0, 1, NA)
sum(r1k, na.rm = T)
r500m <- ifelse(focalWeight(r30, d = 500, type = "circle")>0, 1, NA)
r250m <- ifelse(focalWeight(r30, d = 250, type = "circle")>0, 1, NA)
r100m <- ifelse(focalWeight(r30, d = 100, type = "circle")>0, 1, NA)

# Give mean number of cells with cut from any year within 1k, 500, 250m circles

# function gives mean values in circular window (ie doesn't sum over whole window square, but inside circle)
f3 <- function(x, ...) sum(x, ...)/sum(!is.na(x), ...) # na.rm is carried through from focal function
# fun = mean, and na.rm = T should give same results, but doesnt' seem to like NAs in the window... 

cut.r1k <- focal(cut.msk, w = r1k, fun = f3, na.rm = TRUE, pad = TRUE) 
cut.r500 <- focal(cut.msk, w = r500m, fun = f3, na.rm = TRUE, pad = TRUE) 
cut.r250 <- focal(cut.msk, w = r250m, fun = f3, na.rm = T, pad = T)

cut40.r1k <- focal(cut.40msk, w = r1k, fun = f3, na.rm = TRUE, pad = TRUE)
cut40.r500 <- focal(cut.40msk, w = r500m, fun = f3, na.rm = TRUE, pad = TRUE)
cut40.r250 <- focal(cut.40msk, w = r250m, fun = f3, na.rm = TRUE, pad = TRUE)

plot(stack(cut.r1k, cut.r500, cut.r250, cut40.r1k, cut40.r500, cut40.r250), colNA = "black")

cutStack <- stack(cut.r, cut.msk, cut.40msk,
                  cut.r1k, cut.r500, cut.r250, 
                  cut40.r1k, cut40.r500, cut40.r250)

cut.names <- c("cut_r", "cut_msk", "cut_40msk",
               "cut_r1k", "cut_r500", "cut_r250", 
               "cut40_r1k", "cut40_r500", "cut40_r250")

names(cutStack) <- cut.names
save(cutStack, cut.names, file = file.path(gis_out, "cut_stack.rdata"))

writeRaster(cutStack, bylayer = T, 
            filename = file.path(gis_out, "r_utm/disturb.tif"), suffix = "names", overwrite = TRUE)

rm(cut.r, cut.msk, cut.40msk,
   cut.r1k, cut.r500, cut.r250, 
   cut40.r1k, cut40.r500, cut40.r250)



### 3.1 Topography - Elevation, slope, aspect, NSS, ESS, TWI #####

# get elevation raster
# be_r30 <- raster(file.path(gis_out, "r_utm/be_r30.tif"))
be_r30
## Create toppographic metrics
terr <- raster::terrain(be_r30, opt = c("slope", "aspect", "TRI"), unit= "degrees")
terr

## Make eastness and northness
Nss <- cos(terr$aspect* pi / 180) # "Northness (aspect)" # convert to radians for cos/sin functions
Ess <- sin(terr$aspect* pi / 180)  # eastness

hist(terr$aspect)
hist(Ess)

# plot(terr$aspect)
# plot(Nss)

# plot(Ess)

# install.packages("dynatopmodel")
## TWI
# dynatopmodel::upslope.area
system.time(
  twi30 <- dynatopmodel::upslope.area(be_r30, atb = T)
)
twi30
# $atb is the twi
save(twi30, file = file.path(gis_out, "twi.rdata"))
# user  system elapsed  # on laptop
# 1299.66    4.90 1424.00


#### 3.2 TPI ######

# Do TPI at three different scales
# get central cellS for each focal weight window (see above - section 2 for windows)
c250 <- ceiling(prod(dim(r250m))/2)
c500 <- ceiling(prod(dim(r500m))/2)
c1k <- ceiling(prod(dim(r1k))/2) 

r250m[c250]

plot(be_r30)

TPI250 <- focal(be_r30, w=r250m, 
                fun=function(x, ...) x[c250] - sum(x[-c250], ...)/sum(!is.na(x[-c250]), ...), 
                na.rm = T, pad=TRUE)
                #,filename = file.path(gis_out, "r_utm/tpi250.tif"), datatype = "FLT4S", overwrite = T)

TPI500 <- focal(be_r30, w=r500m, 
                fun=function(x, ...) x[c500] - sum(x[-c500], ...)/sum(!is.na(x[-c500]), ...), 
                na.rm = T, pad=TRUE)
                #,filename = file.path(gis_out, "r_utm/tpi500.tif"), datatype = "FLT4S")

TPI1k <- focal(be_r30, w=r1k, 
               fun=function(x, ...) x[c1k] - sum(x[-c1k], ...)/sum(!is.na(x[-c1k]), ...), 
               na.rm = T, pad=TRUE)
               # ,filename = file.path(gis_out, "r_utm/tpi1k.tif"), datatype = "FLT4S")

plot(stack(be_r30, TPI250, TPI500, TPI1k))

## Save all topography predictors
terr30 <- stack(be_r30, terr, Nss, Ess, twi30$atb, TPI250, TPI500, TPI1k)
terr30.names <- c("be30","tri30","slope30","aspect30","Nss30","Ess30","twi30", "tpi250", "tpi500", "tpi1k")
names(terr30) <- terr30.names

writeRaster(terr30, bylayer = T, filename = file.path(gis_out, "r_utm/terr30.tif"), 
            datatype = "FLT4S", suffix = "names", overwrite = TRUE)

save(terr30, terr30.names, file = file.path(gis_out, "terr30.rdata"))


#### 4. Lidar cover ######

# List lidar cover files
list.files(file.path(gis_in, "lidar"), "lidar_metric")

lidar.fn <- list.files(file.path(gis_in, "lidar"), "lidar_metric", full.names = T)

lidarStck <- stack(lidar.fn)
lidarStck
lidarStck[[1]]

rp <- dismo::randomPoints(lidarStck, 500)
rp_extr <- raster::extract(lidarStck, rp)
pairs(rp_extr)

# use first return data
lidar.fn <- lidar.fn[!grepl("_all", lidar.fn)]
lidar.fn

# create stack with reduced lidar cover
lidarStck <- stack(lidar.fn)
lidarStck

plot(lidarStck)

# Some extreme values.. 
plot(lidarStck$lidar_metric_mosaic_p95)
hist(lidarStck, maxpixels = 750000)

sum(values(lidarStck$lidar_metric_mosaic_p95)>100, na.rm = T)
# some extreme values...  6. remove these
sum(values(lidarStck$lidar_metric_mosaic_p25)>100, na.rm = T)
sum(values(lidarStck$lidar_metric_mosaic_rumple)>100, na.rm = T)

lidarStck$lidar_metric_mosaic_p95[lidarStck$lidar_metric_mosaic_p95 > 100] <- NA
lidarStck$lidar_metric_mosaic_p25[lidarStck$lidar_metric_mosaic_p25 > 100] <- NA
lidarStck$lidar_metric_mosaic_rumple[lidarStck$lidar_metric_mosaic_rumple > 100] <- NA

hist(lidarStck, maxpixels = 750000)
par(mfrow = c(1,1))
plot(lidarStck$lidar_metric_mosaic_rumple)

names(lidarStck)
lid.names <- c("l_Cover_2m_4m","l_Cover_2m_max","l_Cover_4m_16m","l_p25","l_p95","l_rumple")
names(lidarStck) <- lid.names

## RESAMPLE TO R30 (different origin, so crop won't get exact extent)
r30
origin(lidarStck)

lidarStck <- raster::resample(lidarStck, r30, method = "bilinear")
extent(lidarStck) == extent(r30)

save(lidarStck, lid.names, file = file.path(gis_out, "lidarStack.rdata"))
writeRaster(lidarStck, bylayer = T, filename = file.path(gis_out, "r_utm/lidar.tif"), 
            datatype = "FLT4S", suffix = "names", overwrite = TRUE)

### 5. Road stream, admin data ######

# "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb"
strm <- st_read(file.path(gis_in, "shape/rivers_large_oregonexplorer.shp"))
rds <- st_read(file.path(gis_in, "shape/roads_bioarea.shp"))
## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja_bound

strm.utm <- st_transform(strm, crs = utm10N)
rds.utm <- st_transform(rds, crs = utm10N)
hja.utm <- st_transform(hja_bound, crs = utm10N)

hja.utm
strm.utm
rds.utm

plot(strm.utm[, "STREAMS_"])
plot(rds.utm, add = T, col = "red", lwd = 2)

## Create distance rasters

# get coordinates at 30m resolution grid and # create point sf object - already done above
xy_r30 <- raster::coordinates(r30)
pts_r30 <- st_as_sf(data.frame(xy_r30, id = 1:nrow(xy_r30)),coords = c("x", "y"), crs = utm10N)

# union streams - (distance to single feature)
strm.utm
strm.union <- st_union(strm.utm)

strm.dist <- st_distance(pts_r30, strm.union)
dim(strm.dist)
head(strm.dist)

distStrm <- r30
distStrm[] <- as.numeric(strm.dist[,1])

plot(distStrm)
plot(strm.utm, add = T, col = "black")

## roads
rds.union <- st_union(rds.utm)
rds.dist <- st_distance(pts_r30, rds.union)
head(rds.dist)

distRds <- r30
distRds[] <- as.numeric(rds.dist[,1])

plot(distRds)
plot(rds.utm, add = T, col = "black")

## Do INSIDE HJA as raster
HJA <- rasterize(hja.utm, r30)
plot(HJA)

# convert to 0 and 1 for outside, inside HJA
HJA[is.na(HJA)] <- 0
plot(HJA)

plot(strm.utm[, "STREAMS_"], add = T)
plot(rds.utm, add = T, col = "black", lwd = 1)

admStack <- stack(distStrm, distRds, HJA)
admStack[[1]] # in memory

adm.names <- c("DistStream", "DistRoad","insideHJA")
names(admStack) <- adm.names

save(admStack, adm.names, file = file.path(gis_out, "admStck.rdata"))
writeRaster(admStack[[1:2]], bylayer = T, filename = file.path(gis_out, "r_utm/adm.tif"), 
            datatype = "FLT4S", suffix = "names", overwrite = TRUE)

writeRaster(admStack$insideHJA, filename = file.path(gis_out, "r_utm/adm_insideHJA.tif"), 
            datatype = "INT1U", overwrite = TRUE)

#### 6. Annual landsat metrics #####

## Download stdDev.tif, quantiles.tif, leastCloud2.tif from Google Earth Engine following G1_gee_landsat_series.js
## And store in file.path(gis_in, "gee")

## Get rasters
std <- brick(file.path(gis_in, "gee/stdDev.tif"))
qnt <- brick(file.path(gis_in, "gee/quantiles.tif"))
cld <- brick(file.path(gis_in, "gee/leastCloud2.tif"))

# set NA value
NAvalue(std) <- -9999
NAvalue(qnt) <- -9999
NAvalue(cld) <- -9999

std
names(std)
names(qnt)
names(cld)
cellStats(std$ndvi_stdDev, range)
NAvalue(cld)
cellStats(cld$LC08_045029_20180726_B2, range)

# crop to raster template
std
r30
std <- raster::crop(std, r30)
qnt <- raster::crop(qnt, r30)
cld <- raster::crop(cld, r30)
extent(r30) == extent(std)

plotRGB(cld, r = 4, g = 3, b = 2, stretch = "lin", colNA = "black")
plot(cld$LC08_045029_20180726_nbr, colNA = "black")

hist(std[["savi_stdDev"]])
sum(values(std[["savi_stdDev"]])>0.5, na.rm = T)
hist(values(std[["savi_stdDev"]])[values(std[["savi_stdDev"]])<=0.5])

plot(std[[c("ndvi_stdDev", "ndmi_stdDev", "nbr_stdDev", "savi_stdDev")]])

plot(qnt[[c("ndvi_p50", "savi_p50")]])

plot(qnt[[c("ndvi_p5", "ndvi_p50", "ndvi_p95", 
            "ndmi_p5", "ndmi_p50", "ndmi_p95",
            "nbr_p5", "nbr_p50", "nbr_p95")]])

pairs(qnt[[c("ndvi_p5", "ndvi_p50", "ndvi_p95", 
             "ndmi_p5", "ndmi_p50", "ndmi_p95",
             "nbr_p5", "nbr_p50", "nbr_p95")]], maxPixels = 500)

## Make subset stack
annualStack <- raster::stack(std[[c("ndmi_stdDev", "nbr_stdDev")]],
                             qnt[[c("ndvi_p5", "ndvi_p50", "ndvi_p95", 
                                    "ndmi_p5", "ndmi_p50", "ndmi_p95")]])

annualStack
plot(annualStack)

## get window stats over 100m, and 250m

# focal doesnt' work on stacks, therefore go through stack as a list with lapply/parApply
annualL <- as.list(annualStack)

# library(parallel)
# ncore <- detectCores()-2
# cl <- makeCluster(ncore)
# 
# clusterExport(cl, c("annualL", "f3", "r100m", "r250m", "r500m"))
# 
# clusterEvalQ(cl, {
#   
#   library(raster)
#   # rasterOptions(tmpdir = "/scratch") # make temp directory, -make rasters availbale after cluster
# })
# 
# focalL <- parLapply(cl, annualL, function(x) {
#   
#   nm <- names(x)
#   tmp100 <- focal(x, w = r100m, fun = f3, na.rm = TRUE, pad = TRUE)
#   tmp250 <- focal(x, w = r250m, fun = f3, na.rm = TRUE, pad = TRUE)
#   tmp500 <- focal(x, w = r500m, fun = f3, na.rm = TRUE, pad = TRUE)
#   
#   tmp <- stack(tmp100, tmp250, tmp500)
#   names(tmp) <- paste0(nm, c("_r100", "_r250", "_r500"))
#   tmp
#   
# })

focalL <- lapply(annualL, function(x) {
  
  nm <- names(x)
  tmp100 <- focal(x, w = r100m, fun = f3, na.rm = TRUE, pad = TRUE)
  tmp250 <- focal(x, w = r250m, fun = f3, na.rm = TRUE, pad = TRUE)
  tmp500 <- focal(x, w = r500m, fun = f3, na.rm = TRUE, pad = TRUE)
  
  tmp <- stack(tmp100, tmp250, tmp500)
  names(tmp) <- paste0(nm, c("_r100", "_r250", "_r500"))
  tmp
  #as.list(tmp)
  
})

annualFocal <- stack(unlist(lapply(focalL, as.list)))
names(annualFocal)

annualStack <- addLayer(annualFocal, cld[[c("LC08_045029_20180726_B1", "LC08_045029_20180726_B3"
                                            ,"LC08_045029_20180726_B4","LC08_045029_20180726_B5",
                                            "LC08_045029_20180726_B7","LC08_045029_20180726_B10")]])

annual.names <- names(annualStack)
save(annualStack, annual.names, file = file.path(gis_out, "annualStack.rdata"))
writeRaster(annualStack, bylayer = T, filename = file.path(gis_out, "r_utm/annual.tif"), 
            datatype = "FLT4S", suffix = "names", overwrite = TRUE)

### 7 . Temperature and precipitation #####

temp.fn <- list.files(file.path(gis_in, "temp_precip"), "\\.tif$", full.names = T)
temp.fn

tp <- stack(temp.fn)
tp
plot(tp, colNA = "black") # 0 values around edge??
cellStats(tp$Average_Annual_Maximum_Temperature__1981_2010_, range)
# check
plot(tp$Average_Annual_Maximum_Temperature__1981_2010_, 
     breaks = seq(0,18,2), col = terrain.colors(9), colNA = "black")

# seems to be 1 cell to xmin, and 1 cell to y max

# crop this
extent(tp)
tp_ext <- as.vector(extent(tp))

tp_new_ext <- extent(tp_ext + c(1000,0,0,-1000))
tp <- crop(tp, tp_new_ext)
tp

plot(tp, colNA = "black") # 0 values around edge??

## project to utm
tp_utm <- projectRaster(tp, projectExtent(tp, utm10N))
tp_utm
plot(tp_utm)

## Scale down ... maybe better direct values to each cell?? 
# x max slighlty out of range...  
tp_r30 <- raster::resample(tp_utm, r30, method = "bilinear")
tp_r30

plot(tp_r30)
names(tp_r30)
tp.names <- c("maxT_annual", "meanT_annual", "minT_annual", "precipitation_mm")
names(tp_r30) <- tp.names

save(tp.names, tp_r30, file = file.path(gis_out, "tp.rdata"))
writeRaster(tp_r30, bylayer = T, filename = file.path(gis_out, "r_utm/tp.tif"), 
            datatype = "FLT4S", suffix = "names", overwrite = TRUE)

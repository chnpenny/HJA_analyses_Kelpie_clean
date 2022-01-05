## Create new data for prediction across study area ####
```{r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/10_eo_data')
	
pacman::p_load('raster','sf','here','dplyr')
	
source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")
	
```


```{r set-names}
utm10N = 32610
	
gis_in = here('..','..','format_data','gis',"raw_gis_data") 
gis_out = here('..','..','format_data','gis',"processed_gis_data") 
	
```


```{r load-data}
## create a reduced prediction area - convex hull around (all points + HJA) + buffer
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm; xy.all.utm
	
## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja = st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound = subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm = st_transform(hja_bound, crs = utm10N)
	
## convex hull 
aoi.pred.sf = st_buffer(st_union(st_convex_hull(st_union(xy.utm)), st_convex_hull(hja.utm)), dist = 500)
	
st_write(aoi.pred.sf, here(gis_out,"s_utm","aoi_pred_sf_500.shp"), delete_layer = T)
	
plot(aoi.pred.sf, add =F, col = NA)
plot(st_geometry(xy.utm),add =T)
plot(hja.utm, add = T, col = NA)
	

# load as brick
allBrck = brick(file.path(gis_out, "r_utm","allStack.tif"))
# get names and name groups
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- allNames
allBrck
names(allBrck)
	
plot(allBrck$insideHJA)
# plot(allBrck$cut_msk)
# plot(allBrck$cut_40msk)

plot(allBrck$be30)
plot(hja.utm, add = T, col = NA)
plot(aoi.pred.sf, add =T, col = NA)
plot(xy.utm, add = T, pch = 16, col = "black")

## Try difference buffer distances
# aoi.pred.lst <- lapply(c(500, 1000, 1500), function(i){
#   st_buffer(st_union(st_convex_hull(st_union(xy.utm)), st_convex_hull(hja.utm)), dist = i)
# })
# 
# ## add to map
# plot(allBrck$be30)
# plot(st_geometry(hja.utm), add = T, col = NA)
# plot(st_geometry(xy.utm), add = T, pch = 16, col = "black")
# sapply(seq_along(aoi.pred.lst), function(x) plot(aoi.pred.lst[[x]], add = T, col = NA, border = x))

### bring in manually edited prediction area outline to replace above
aoi.pred.sf_edit <- st_read(file.path(gis_out, "s_utm/aoi_pred_sf_edit.shp"))


## Add transformed raster here
allBrck$lg_DistStream <- log(allBrck$DistStream + 0.001)
allBrck$lg_DistRoad = log(allBrck$DistRoad + 0.001)
allBrck$lg_cover2m_max = log(allBrck$l_Cover_2m_max + 0.001)
allBrck$lg_cover2m_4m = log(allBrck$l_Cover_2m_4m + 0.001)
allBrck$lg_cover4m_16m = log(allBrck$l_Cover_4m_16m + 0.001)

plot(allBrck[[c("lg_DistStream", "DistStream")]])
plot(allBrck[[c("lg_DistRoad", "DistRoad")]])

plot(allBrck[[c("lg_cover2m_max", "l_Cover_2m_max", "lg_cover2m_4m", "l_Cover_2m_4m")]])

# get names for excel description
# write.csv(names(allBrck), "clipboard", row.names = F)

origin(allBrck)

##### CROP #####
## crop to reduced aoi.pred - get extent of the convex hull around sample points and HJA
aoi.pred.ext <- adjExt(st_bbox(aoi.pred.sf_edit), outF = "Extent", d = 30) ####
# make raster template from above extent
r.aoi.pred <- raster(aoi.pred.ext, crs = utm10N, res = 30)

origin(r.aoi.pred)

# mask it with aoi (points and HJA)
r.aoi.pred[] <- 1
r.aoi.pred <- raster::mask(r.aoi.pred, st_as_sf(aoi.pred.sf_edit))
plot(r.aoi.pred, colNA = "grey")

plot(hja.utm, add = T, col = NA)
plot(aoi.pred.sf, add =T, border = "blue")
plot(aoi.pred.sf_edit, add =T, col = NA)
plot(xy.utm, add = T, pch = 16, col = "black")

## crop raster brick
allBrck <- crop(allBrck, r.aoi.pred)

origin(allBrck)
allBrck
r.aoi.pred

names(allBrck)
brNames <- names(allBrck)
# names(allBrck) <- c(allNames, "lg_DistStream", "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m", "lg_cover4m_16m")
# save allBrck names
save(brNames, file = file.path(gis_out, "brNames.rdata"))

plot(allBrck$gt4_r30)

indNA <- complete.cases(values(dropLayer(allBrck, "cut_r")))
sum(!indNA)

## add NAs
r.msk <- r.aoi.pred
r.msk[!indNA] <- NA

plot(stack(r.aoi.pred, r.msk), colNA = "black")

## UPdate index of complete cases including those from new aoi edit
indNA <- !is.na(values(r.msk))

## Mask allBrick
allBrck <- mask(allBrck, r.msk)
plot(allBrck$gt4_r30)

# save raster as single tif
writeRaster(allBrck, filename = file.path(gis_out, "r_utm/allStack_aoi.tif"), overwrite = TRUE)

# plot(r, colNA = "black")
save(r.msk, indNA, r.aoi.pred, aoi.pred.sf, file = "Hmsc_CD/oregon_ada/data/gis/templateRaster.rdata")
save(r.msk, indNA, r.aoi.pred, aoi.pred.sf, file = file.path(gis_out, "templateRaster.rdata"))
# load(file.path(gis_out, "templateRaster.rdata"))

## Scale whole data set - apart from categorical predictors
allBrck.sc <- scale(dropLayer(allBrck, c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")))
# stores scale parameters in the @data slot
allBrck.sc # in memory
# str(allBrck.sc)
inMemory(allBrck.sc[[1]])

sapply(1:nlayers(allBrck.sc), function(x) inMemory(allBrck.sc[[x]]))

## add back categorical - but bring into memory first
catRasters <- readAll(allBrck[[c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")]])
catRasters[[1]]

allBrck.sc <- addLayer(allBrck.sc, catRasters)
names(allBrck.sc)

## save scaled rasters
# save(r.msk, r.aoi.pred, indNA, allBrck.sc, file = "Hmsc_CD/oregon_ada/data/gis/predRaster_sc.rdata")
# load("Hmsc_CD/oregon_ada/data/gis/predRaster_sc.rdata")
save(r.msk, r.aoi.pred, indNA, allBrck.sc, file = file.path(gis_out, "r_oversize/predRaster_sc.rdata"))

# save scaled raster
writeRaster(allBrck.sc, filename = file.path(gis_out, "r_utm/AllStack_aoi_sc.tif"), overwrite = TRUE)

#### MAKE NEW DATA #####

# Get new data prior from unscaled version
## extract site env vars from original data set
allVars <- data.frame(xy.all.utm[, c("SiteName", "trap", "period")], raster::extract(allBrck, xy.all.utm))
head(allVars)
# allVars$geometry <- NULL

## Get new data as data.frame:
allBrck
newData <- data.frame(values(dropLayer(allBrck, "cut_r")))
# remove NAs
newData <- newData[indNA, ] # only complete cases
## change categorical to predictor values to match model
newData[, "insideHJA"] <- ifelse(newData[, "insideHJA"] == 0, "no", "yes")
table(newData[,"insideHJA"], useNA = "always")

save(allVars, newData, indNA, file = file.path(gis_out, "r_oversize/newData_unscaled.rdata"))
# load(file.path(gis_out, "r_oversize/newData_unscaled.rdata")) # r.msk, r.aoi.pred, indNA, allBrck.sc

## data frame of coordinates
newXY <- coordinates(r.msk)

newXY.sc <- scale(newXY)
str(newXY.sc)

## Load sample points
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm # and xy.all.utm - all sites
xy.all.utm

xy.sites <- st_coordinates(xy.all.utm) ## OJO some of the coordinates are repeated (different site names)
head(xy.sites)

# scale sample site coords with same parameters as complete data set
attr(newXY.sc, "scaled:center")
attr(newXY.sc, "scaled:scale")

xy.sites.sc <- cbind(
  X = (xy.sites[,"X"] - attr(newXY.sc, "scaled:center")["x"])/attr(newXY.sc, "scaled:scale")["x"],
  Y = (xy.sites[,"Y"] - attr(newXY.sc, "scaled:center")["y"])/attr(newXY.sc, "scaled:scale")["y"]
)

# Join sitenames and change names
xy.sites.sc <- data.frame(xy.sites.sc, xy.all.utm[, c("SiteName", "trap", "period")])
colnames(xy.sites.sc) <- c("UTM_E", "UTM_N", "SiteName", "trap", "period")
head(xy.sites.sc)

# fix colnames and NAs
colnames(newXY.sc) <- c("UTM_E", "UTM_N")
## filter for NAs in predictors

newXY.sc <- newXY.sc[indNA,]
dim(newXY.sc)
sum(!complete.cases(newXY.sc))

## extract site env vars from scaled data set
allVars.sc <- data.frame(xy.all.utm[, c("SiteName", "trap", "period")], raster::extract(allBrck.sc, xy.all.utm))
head(allVars.sc)

## Get new data as data.frame:
allBrck.sc
newData.sc <- data.frame(values(dropLayer(allBrck.sc, "cut_r")))
head(newData.sc)
str(newData.sc)

length(indNA) == ncell(allBrck.sc)

### remove NAs, for faster prediction and less storage - add back in for raster creation
newData.sc <- newData.sc[indNA, ] # only complete cases

nrow(newData.sc) + sum(!indNA) == ncell(allBrck)

sum(!complete.cases(newData.sc))

## change categorical to predictor values to match model
newData.sc[, "insideHJA"] <- ifelse(newData.sc[, "insideHJA"] == 0, "no", "yes")
table(newData.sc[,"insideHJA"], useNA = "always")

summary(newData.sc)

#save(newData.sc, xy.sites.sc, newXY.sc, allVars.sc, file = "Hmsc_CD/oregon_ada/data/newData_scaled.rdata")
save(newData.sc, xy.sites.sc, newXY.sc, allVars.sc, 
     file = file.path(gis_out, "r_oversize/newData_scaled.rdata"))
# load(file.path(gis_out, "r_oversize/newData_scaled.rdata"))

#### Make clamped data set #####

varsOut <- c("SiteName", "trap", "period", "cut_r", "cut_msk", "cut_40msk", "insideHJA", # "geometry"
              "aspect30","maxT_annual","meanT_annual","minT_annual","precipitation_mm")

## get range of predictor values at sample sites
sample_range <- allVars %>%
  filter(period == "S1") %>%
  select(!all_of(varsOut)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value")%>%
  group_by(predictors)%>%
  summarise(min = min(value),
            max = max(value))%>%
  arrange(predictors)

sample_range

head(newData)

newData_clamp <- newData %>%
  select(!any_of(varsOut)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value") %>%
  left_join(y = sample_range)%>%
  mutate(
    value_clamp = case_when(
    value < min ~ min,
    value > max ~ max,
    TRUE ~ value),
    value_na = case_when(
      value < min ~ NA_real_,
      value > max ~ NA_real_,
      TRUE ~ value
    ))# %>% arrange(predictors)


newData_clamp

## Change to wide format
newData_clamp_wide <- do.call(data.frame, 
                        lapply(unique(newData_clamp$predictors), 
                               function(x) newData_clamp[newData_clamp$predictors == x, "value_clamp"]))
colnames(newData_clamp_wide) <- unique(newData_clamp$predictors)
head(newData_clamp_wide)

newData_clamp_wide$insideHJA <- newData$insideHJA

save(newData_clamp_wide, indNA, file = file.path(gis_out, "r_oversize/newData_clamp.rdata"))
rm(newData_clamp, newData_clamp_wide, sample_range); gc()
#load(file.path(gis_out, "r_oversize/newData_clamp.rdata"))

## Scaled version - clamp #### 
sample_range <- allVars.sc %>%
  filter(period == "S1") %>%
  select(!all_of(varsOut)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value")%>%
  group_by(predictors)%>%
  summarise(min = min(value),
            max = max(value))%>%
  arrange(predictors)

sample_range

head(newData.sc)

newData_clamp.sc <- newData.sc %>%
  select(!any_of(varsOut)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value") %>%
  left_join(y = sample_range)%>%
  mutate(
    value_clamp = case_when(
      value < min ~ min,
      value > max ~ max,
      TRUE ~ value),
    value_na = case_when(
      value < min ~ NA_real_,
      value > max ~ NA_real_,
      TRUE ~ value
    ))# %>% arrange(predictors)


newData_clamp.sc

## Change to wide format
newData_clamp_wide.sc <- do.call(data.frame, 
                              lapply(unique(newData_clamp.sc$predictors), 
                                     function(x) newData_clamp.sc[newData_clamp.sc$predictors == x, "value_clamp"]))
colnames(newData_clamp_wide.sc) <- unique(newData_clamp.sc$predictors)
head(newData_clamp_wide.sc)
newData_clamp_wide.sc$insideHJA <- newData.sc$insideHJA


save(newData_clamp_wide.sc, xy.sites.sc, newXY.sc, allVars.sc, 
     file = file.path(gis_out, "r_oversize/newData_scaled_clamp.rdata"))

rm(newData_clamp.sc, newData_clamp_wide.sc, sample_range); gc()

#### Make raster to show geogrpahical distribution of clamped predictors

# Get NA values version in wide format
newData_clamp_na_wide <- do.call(data.frame, 
                                 lapply(unique(newData_clamp$predictors), 
                                        function(x) newData_clamp[newData_clamp$predictors == x, "value_na"]))
colnames(newData_clamp_na_wide) <- unique(newData_clamp$predictors)
head(newData_clamp_na_wide)
newData_clamp_na_wide$insideHJA <- newData$insideHJA


## make rasters
rList <- lapply(newData_clamp_na_wide, function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- x
  tmp
  
})

# plot(tmp)

rStack_clamp <- stack(rList)
names(rStack_clamp) <- colnames(newData_clamp_na_wide)
rStack_clamp

plot(rStack_clamp[[1:10]], colNA = "black")
plot(rStack_clamp[[11:22]], colNA = "black")
plot(rStack_clamp[[23:35]], colNA = "black")
plot(rStack_clamp[[35:51]], colNA = "black")

## sum to get all NAs
smStack_cl <- is.na(sum(dropLayer(rStack_clamp, "insideHJA")))
plot(smStack_cl, colNA  = "black")

save(rStack_clamp, file = file.path(gis_out, "clamp_grid_na.rdata"))

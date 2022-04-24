
```{r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/10_eo_data')
	
pacman::p_load('raster','sf','here')
	
here()
# check working directory
# HJA_analyses_Kelpie_clean # is the root (wd) folder and must have a .Rproj file in it for here() to work.

```


```{r set-names}
gis_in = here('format_data','gis',"raw_gis_data") 
gis_out = here('format_data','gis',"processed_gis_data")  
	
```


```{r load-data}
### 1. Load processed data ####

## Add all rasters to a stack, export as single multilayered tif and import as rasterbrick.


## 1 Canopy
load(file.path(gis_out, "be_ht.rdata"))
## 2 Cut
load(file.path(gis_out, "cut_stack.rdata"))
## 3. Topography
load(file.path(gis_out, "terr30.rdata"))
## 4. Lidar
load(file.path(gis_out, "lidarStack.rdata"))
## 5. Streams/Roads
load(file.path(gis_out, "admStck.rdata"))
## 6. Annual indices
load(file.path(gis_out, "annualStack.rdata"))
## 7. Temperature
load(file.path(gis_out, "tp.rdata"))
	
## Stack up 
allStck <- stack(be_ht, cutStack, terr30, lidarStck, admStack, annualStack, tp_r30)
allNames <- names(allStck)
allNames
	
# get name groups
nameList <- list(be_ht = be_ht.names, cutStack = cut.names, terr30 = terr30.names, 
                  lidarStck = lid.names, admStack = adm.names, annualStack = annual.names, tp_r30 = tp.names)

writeRaster(allStck, bylayer = F, filename = file.path(gis_out, "r_utm/allStack.tif"), overwrite = TRUE)
save(allNames, nameList, file = file.path(gis_out, "allNames.rdata"))

rm(be_ht, cutStack, terr30, lidarStck, admStack, annualStack, tp_r30,
   adm.names, annual.names, be_ht.names, cut.names, lid.names, terr30.names, tp.names, all.names)

## 2. Extract point values ####

# load as brick
allBrck <- brick(file.path(gis_out, "r_utm/allStack.tif"))
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- allNames
allBrck

write.table(data.frame(band = seq_along(allNames), predictor= allNames), 
            file = file.path(gis_out, "r_utm/allStack.txt"), row.names = F)

rm(allNames, nameList)

## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

allVars <- data.frame(SiteName = xy.utm$SiteName, extract(allBrck, xy.utm))
head(allVars)

allVars$insideHJA <- ifelse(allVars$insideHJA==0, "no", "yes")
allVars$insideHJA <- factor(allVars$insideHJA, levels = c("no", "yes"))
table(allVars$insideHJA)

str(allVars)

summary(allVars)

save(allVars, file = file.path(gis_out, "envVars.rdata"))
save(allVars, file = file.path("Hmsc_CD/oregon_ada/data", "envVars.rdata"))

# load(file.path(gis_out, "envVars.rdata"))


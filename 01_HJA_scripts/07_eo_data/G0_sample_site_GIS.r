
### Create spatial object from sample points ####


```{r setup}
rm(list=ls())
q()
	
# HJA_analyses_Kelpie_clean # is the root folder and must have a .Rproj file in it for here() to work.
	
pacman::p_load('tidyverse','sf','here','conflicted','glue','pROC', 'gridExtra','ggeffects','corrplot')
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
	
here()
	

```


```{r set-names}
# ..... UTM setting .....
# wgs84 UTM 10N
utm10N = 32610
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 = 26910
# EPSG:4269 # NAD 83
# nad83 <- 4269
	
# ....... folder structure .......
gis_in = here('03_format_data','gis',"raw_gis_data") 
gis_out = here('03_format_data','gis',"processed_gis_data") 
	
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
	
datFile = glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')
# "_uncorr.csv" was used, should be same for the purpose of this script?
otupath = here('02_Kelpie_maps',outputidxstatstabulatefolder)
	

```


```{r load-data}
otuenv = read.csv(here(otupath, datFile), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
# otuenv[1:6,1:10]
coords = unique(otuenv[,c("SiteName","UTM_E", "UTM_N")])
xy.sf = st_as_sf(coords, coords = c("UTM_E", "UTM_N"), crs = nadutm10)
	
coords_allSites = unique(otuenv[,c("SiteName","trap", "period", "UTM_E", "UTM_N")])
xy.allSites.sf = st_as_sf(coords_allSites, coords = c("UTM_E", "UTM_N"), crs = nadutm10)
	
rm(otuenv, outputidxstatstabulatefolder, datFile, primer, 
   kelpierundate, minimaprundate, samtoolsfilter, samtoolsqual, coords)
	
# transform to wgs utm to match rasters
xy.utm = st_transform(xy.sf, crs = utm10N)
xy.all.utm = st_transform(xy.allSites.sf, crs = utm10N)
rm(xy.sf, xy.allSites.sf)
	
# ... save ...
st_write(xy.utm, file.path(gis_out, "s_utm/sample_sites_utm10.shp"), delete_layer = T)
st_write(xy.utm, file.path(gis_out, "s_utm/sample_sites_utm10.kml"), delete_layer = T)
	
save(xy.utm, xy.all.utm, file = file.path(gis_out, "sample_sites.rdata"))
	

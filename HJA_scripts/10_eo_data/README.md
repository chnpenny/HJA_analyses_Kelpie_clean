### Scripts to process jSDM predictors


- G0_sample_site_GIS.r  
Extract sample site coordinates, convert to `sf` and save

- G1_gee_landsat_series.js  
Download remote sensed data from Google Earth Engine

- G2_process_rasters.r  
Reproject rasters, crop to standard study area extent and process (summarise over windows, etc)

- G3_create_raster_stack.r  
Stack all rasters and save

- G4_make_newData.r  
Create new data frame with all raster values (scaled and unscaled) for predicting models to create distribution maps

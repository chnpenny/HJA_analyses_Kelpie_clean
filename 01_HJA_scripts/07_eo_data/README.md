### Scripts to process sjSDM predictors


- G0_sample_site_GIS.r  
Extract sample site coordinates, convert to `sf` and save

- G1_gee_landsat_series.js  
Download remote sensed data from Google Earth Engine (vegetation indices - annual summaries, variation, 5, median, 95 percentiles). This script is run on Google Earth Engine. A user account is required. The bounding box of the study area needs to be uploaded to the user account assets folder before running the script. The script downloads the resulting files to the user's google drive account (in an appropriate folder). From here, they must be copied to the appropriate folder (03_format_data|gis|raw_gis_data|gee).

- G2_process_rasters.r  
Reproject rasters, crop to standard study area extent and process (summarise over windows, etc)

- G3_create_raster_stack.r  
Stack all rasters and save

- G4_make_newData.r  
Create new data frame with all raster values (scaled and unscaled) for predicting models to create distribution maps

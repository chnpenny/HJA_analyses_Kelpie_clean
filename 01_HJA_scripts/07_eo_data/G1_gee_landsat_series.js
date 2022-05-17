// GEE script for annual Landsat series 

// https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR

// GET EXTERNAL DATA
// import  aoi - bbox polygon
var reg = ee.FeatureCollection('users/cdevenish/misc/oregon_bb_pts');
print(reg);

// Check
// Map.addLayer(reg, {}, 'region', true, 0.4);
// Map.centerObject(reg);

// Define the visualization parameters.
var vis = {bands: ['B4', 'B3', 'B2'], min: -1000, max:5000};
var std_vis = {bands: ['ndvi_stdDev'], min:0 , max:0.5};

// CLOUD MASK function
// https://github.com/fitoprincipe/geetools-code-editor/wiki/Cloud-Masks

var cloud_masks = require('users/fitoprincipe/geetools:cloud_masks');
print(cloud_masks.help['landsatSR']);

var landsatSR_masks = cloud_masks.landsatSR(['cloud', 'shadow', 'snow']);


// FUNCTIONS for vegetation indices
// https://www.usgs.gov/land-resources/nli/landsat/landsat-surface-reflectance-derived-spectral-indices

// Function to calculate and add an NDVI band
//In Landsat 8, NDVI = (Band 5 – Band 4) / (Band 5 + Band 4).
var addNDVI = function(image) {
  return image.addBands(image.normalizedDifference(['B5', 'B4']).rename('ndvi')); // rename band to 'ndvi'
};
  
// Function to calculate and add an NDMI band
// In Landsat 8, NDMI = (Band 5 – Band 6) / (Band 5 + Band 6).
var addNDMI = function(image) {
return image.addBands(image.normalizedDifference(['B5', 'B6']).rename('ndmi'));
};

// NBR
// In Landsat 8, NBR = (Band 5 – Band 7) / (Band 5 + Band 7).
var addNBR = function(image) {
return image.addBands(image.normalizedDifference(['B5', 'B7']).rename('nbr'));
};
// NBR2
// In Landsat 8, NBR2 = (Band 6 – Band 7) / (Band 6 + Band 7)

// SAVI
// In Landsat 8, SAVI = ((Band 5 – Band 4) / (Band 5 + Band 4 + 0.5)) * (1.5).
var addSAVI = function(image) {
return image.addBands(image.expression('((b(5) - b(4)) / (b(5) + b(4) + 0.5)) * (1.5)').rename('savi'));
};


/*
var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B5'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
});
*/


// This function adds a band representing the image timestamp.
var addTime = function(image) {
  return image.addBands(image.metadata('system:time_start')
    // Convert milliseconds from epoch to days to aid in
    // interpretation of the following trend calculation. Maybe in days??? 
    .divide(1000 * 60 * 60 * 24));
};


//load images for composite
var L8sr = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
  .filterBounds(reg) // // just images intersecting with this
  .filterMetadata('CLOUD_COVER', 'less_than', 30) 
  .filterDate('2018-01-01','2018-12-31')
  // .map(function(image) { return image.clip(reg); }) // clip here to mask extent 
  .map(landsatSR_masks)
  //.map(addTime) // add a band with time as numeric - only for regression slope analysis
  .map(addNDVI) // Add NDVI band to image collection
  .map(addNDMI) // Add NDMI band to image collection
  .map(addNBR) // Add NBR band to image collection
  .map(addSAVI) // Add NBR band to image collection
  .select(['B1', 'B2', 'B3', 'B4','B5', 'B6', 'B7', 'B10', 'B11','ndvi', 'ndmi', 'nbr', 'savi']) // remove QA60. , 'system:time_start'
  ;
  
print('no images in collection:', L8sr.size());

// Check first
var img1 = ee.Image(L8sr.first());
print(img1);
var l8sr_med = ee.Image(L8sr.median()); // create single image composite using median

Map.centerObject(reg, 9);
// Map.addLayer(img1, vis, '1st image');
// Map.addLayer(l8sr_med, vis, 'composite - all');
// Map.addLayer(reg, {}, 'region', true, 0.4);


// get annual variation
// standard deviation across all pixels in time period, for each band
var stdDev = L8sr.reduce(ee.Reducer.stdDev());

// Get median, min, max, and 5 and 95 percentiles
var q = L8sr.reduce(ee.Reducer.percentile([0,5,50,95,100]));

// Get count
var count = L8sr.reduce(ee.Reducer.count());
var noPix = count.select("B2_count"); 

print(stdDev, 'stdDev');
// get max and min
var minD = stdDev.select(['ndvi.stdDev']).reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: reg.geometry(),
  scale: 30,
  maxPixels: 1e9
});
print(minD, 'min std Dev');

var maxD = stdDev.select(['ndvi_stdDev']).reduceRegion({
  reducer: ee.Reducer.max(),
  geometry: reg.geometry(),
  scale: 30,
  maxPixels: 1e9
});

print(maxD, 'max std Dev');

// Get top three least cloudy images over summer
var sortCloud = L8sr.filterDate('2018-07-01', '2018-08-31')
                    .sort('CLOUD_COVER');

// need to mosaic first, as still can be a small fraction of image in the 'reg'

// convert to list to pick top two
var sortCloudL = sortCloud.toList(2);

// convert back to multiband image (via ImageCollection) and make all same data type
var cld = ee.ImageCollection.fromImages(sortCloudL).toBands().toFloat();
print(cld, 'cloud 1 2 3');
print(cld.bandNames(), 'least cloudy 3 - bandNames');

// vis
var vis_cld = {bands: ['LC08_045029_20180726_B4', 'LC08_045029_20180726_B3', 'LC08_045029_20180726_B2'], min: -1000, max:5000};
var styling = {color: 'blue', fillColor: '00000000'};
var ndviParams = {min: -1, max: 1, palette: ['blue', 'white', 'green']};

// Set NA to -9999 prior to export - and clip for viewing here
stdDev = stdDev.unmask(-9999).clip(reg);
q = q.unmask(-9999);
noPix = noPix.unmask(-9999);
cld = cld.unmask(-9999).clip(reg);

//Map.addLayer(cld, vis_cld, 'least cloudy image');

Map.addLayer(cld.select(['LC08_045029_20180726_savi']), ndviParams, 'savi');
Map.addLayer(cld.select(['LC08_045029_20180726_ndvi']), ndviParams, 'ndvi');

//Map.addLayer(stdDev, std_vis, 'std');
Map.addLayer(reg.style(styling), {}, 'region', true, 0.4);

print(stdDev, 'stdDev');
print(q, 'q');


// Change to smae data type for export
q = q.toFloat(); // unify to same float

print(q.bandNames(), 'q bandNames');
print(stdDev.bandNames(), 'stdDev bandNames');


// Export


Export.image.toDrive({
  image: cld,
  description: 'leastCloud2',
  region: reg,
  scale: 30,
  crs: 'EPSG:32610', // this is  wgs84 UTM 10N
  fileFormat: 'GeoTIFF',
  maxPixels: 10000000000,
  folder: 'Oregon'
});


Export.image.toDrive({
  image: stdDev,
  description: 'stdDev',
  region: reg,
  scale: 30,
  crs: 'EPSG:32610', // this is  wgs84 UTM 10N
  fileFormat: 'GeoTIFF',
  maxPixels: 10000000000,
  folder: 'Oregon'
});


Export.image.toDrive({
  image: q,
  description: 'quantiles',
  region: reg,
  scale: 30,
  crs: 'EPSG:32610', // this is  wgs84 UTM 10N
  fileFormat: 'GeoTIFF',
  maxPixels: 10000000000,
  folder: 'Oregon'
});

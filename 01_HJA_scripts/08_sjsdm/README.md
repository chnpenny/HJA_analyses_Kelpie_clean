### Scripts to run sjSDM


- 0_sjsdm_data_prep.rmd  
  (product data for cross-validation and running the best model afterwards)

- 1_sjsdm_tuning.rmd  
  (script for 5-folds cross-validation)
  (will take very long time to run on a local computer)

- 1_sjsdm_tuning-ada.r  
  (r script format used in cluster)

- 2_sjsdm_tuned_model.rmd  
  (script to run the best model selected by cross-validation)

- 3_xAI.rmd  
  (script to carry out explanatory AI analyses)

- 4_sjsdm_data_describe.rmd  
  (script for making some descriptive graphs)

08_sjsdm.Rproj is the r project file for the scripts in this folder if rstudio is used

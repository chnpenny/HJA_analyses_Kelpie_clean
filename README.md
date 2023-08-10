# HJA_analyses_Kelpie_clean
This repo contains data, scripts and results for the paper:

Yuanheng Li, Christian Devenish, Marie I. Tosa, Mingjie Luo, Paul Greenfield, David M. Bell, Damon B. Lesmeister, Maximilian Pichler, Taal Levi, Douglas W. Yu. **Combining environmental DNA and remote sensing for efficient, fine-scale mapping of arthropod biodiversity**

The folder '01_HJA_scripts' contains all the scripts for the bioinformatic pipeline (subfolders 01 to 06), the extraction of GIS and other environmental covariates (subfolders 06 and 07), joint species distribution modelling (subfolder 08), and post-model analysis and visualization of species distributions (subfolder 09). Note that some intermediate files related to environmental variables extaction are not in the repo due to their large sizes. The majority of the scripts are run in bash and R.

The folder '02_Kelpie-maps' holds the result of the bioinformatic pipeline, (the OTU table and environmental covariates combined), and the folder '03_format_data' has some of the data of environmental covariate extraction, which is needed for distribution modelling. The folder '04_Output' contains the results of model fitting and visualizing species distributions. The folder '05_supplement' contains data and figures that are listed in the supplementary information of the paper.

Dependencies of softwares and packages used are listed in the supplementary information of the paper.

The main text figures are located [here](04_Output/figures) and the supplementary figures are located [here](05_supplement/Plots)

The raw fastq files are at NCBI Short Read Archive BioProject PRJNA869351 

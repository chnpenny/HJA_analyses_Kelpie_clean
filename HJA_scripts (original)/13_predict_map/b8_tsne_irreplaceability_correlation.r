 ## TSNE Correlation #####
 
```{r setup}
here()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/13_predict_map')
	
pacman::p_load('dplyr', 'rgdal', 'raster','here','glue','corrplot','ggplot2')
	
packageVersion('sjSDM')
	
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
```

```{r set-names}
abund = "pa"
date.model.run = '20210722'
varsName = 'vars11'
	
outputpath = here('..','..','Output')
modFolder = file.path(outputpath, "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
plotFolder = file.path(outputpath, "prediction_map")
	
```


```{r load-data}
#### Load TSNE data 
load(file.path(plotFolder, 'rdata', "ord_tsne_res_cl_p50.rdata"))
# tsne, r, rSites1, rSites2, NAs, 
	
str(tsne, max.level = 1)
	
### Load irreplaceability data ##
## Load predictor data
## Load new data for prediction and new scaled data
load(file.path(plotFolder, 'rdata', "newData_unscaled.rdata")) # allVars, newData, indNA, 
	
## Final set of VIF chosen predictors
vars11 = c("gt4_500", "cut_r1k","cut_r250","cut40_r1k","cut40_r250","be30","tri30","Nss30",
            "Ess30","twi30","tpi250","tpi1k","l_rumple","nbr_stdDev_r100","ndvi_p5_r100",
            "ndvi_p5_r500","ndvi_p50_r100","ndvi_p50_r500","ndmi_p95_r100",
            "LC08_045029_20180726_B1","LC08_045029_20180726_B5","lg_DistStream",
            "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m") # insideHJA

	
```


```{r plotting}
dim(tsne$Y)
dim(newData[,vars11])
	
mod.cor = cor(tsne$Y, newData[,vars11])
mod.cor
	
getwd()
pdf(file.path(plotFolder, 'plot', "tsne_nopca_predictors_corrplot.pdf"), width = 8, height = 8)
corrplot::corrplot(t(mod.cor), 
           is.corr = F,
           method = "ellipse",
           cl.pos = "n",
           col.lim = c(-1,1),
           title = "tsne",
           oma = c(0,0,0,0),
           #oma = c(2,2,5,1),
           mar = c(0,0,1,0),
           addCoef.col = "black",
           addCoefasPercent = T,
           number.cex = 0.5)

dev.off()
	
shell.exec(file.path(getwd(), plotFolder, 'plot', "tsne_nopca_predictors_corrplot.pdf")) # check!!!
	
## scatter plots

set.seed(45)
ind = sample(1:nrow(newData), size = 5000)
	

pdf(file.path(plotFolder, 'plot', "tsne_nopca_predictors_scatterplot.pdf"), width = 2.5, height = 30)
par(mfcol=c(length(vars11), ncol(tsne$Y)), mar = c(2,2,2,2))
sapply(vars11, function(x) smoothScatter(
  newData[ind,x], tsne$Y[ind,1], 
  xlab = "", ylab = "", main = x, cex.main = 0.5, axes = T, cex.axis = 0.4, cex = 0.1, pch = 16))
sapply(vars11, function(x) smoothScatter(
  newData[ind,x], tsne$Y[ind,2], 
  xlab = "", ylab = "", main = x, cex.main = 0.5, axes = T, cex.axis = 0.4, cex = 0.1, pch = 16))
	
dev.off()
	
shell.exec(file.path(getwd(), plotFolder, 'plot', "tsne_nopca_predictors_scatterplot.pdf"))#???
	
```



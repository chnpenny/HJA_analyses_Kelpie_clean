#### ecocopula analysis

## copied initially from oregon_ada/code/P1_ecocopula_analysis_v3.r

```{r setup}
here()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/13_predict_map')
	
pacman::p_load('dplyr', 'rgdal', 'raster','here','glue','raster','Rtsne')
	
# library(mvabund)
# library(ecoCopula)
```

```{r set-names}
abund = "pa"
date.model.run = '20210722'
varsName = 'vars11'
	
outputpath = here('..','..','Output')
	
datapath = here('..','12_sjsdm')
gis_out = gis_in = here('..','..','format_data','gis')
	
modFolder = file.path(outputpath, "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
plotFolder = file.path(outputpath, "prediction_map")
	
resFolder = here('..','..','format_data','otu')
	
#baseFolder <- "code_sjSDM/r20210716b"
```

```{r load-data}
# load model data - for species classification
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
# rm(env.vars, k, noSteps, vars, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv
	
## load species AUC resutls for filtering
load(file.path(modFolder, "sp_test_results.rdata")) # # eval.results, sp.res.test, sp.res.train
	
## Mean AUC per species (and other eval metrics)
str(sp.res.test, max.level = 1)
head(sp.res.test$auc)
	
## Filter species by auc
auc.filt = 0.70
# threshold for presence absence data
# tr <- 0.5
	
# how many species after AUC filter?
sum(sp.res.test$auc > auc.filt, na.rm = T)
	
# incidence 
incidence = colSums(otu.pa.csv)/nrow(otu.pa.csv)
	
# extrapolated predictions
load(file.path(plotFolder, 'rdata', paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata")))
# pred.mn, pred.sd, 
	
# clamp predictions
load(file.path(plotFolder, 'rdata', paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, "_clamp", ".rdata")))
# pred.mn.cl, pred.sd.cl
	

## local
# load(file.path(gis_out, "r_oversize", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) 

dim(pred.mn)

## filter for species performance
pred.in = pred.mn[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
dim(pred.in)

pred.in.cl = pred.mn.cl[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
	
## load raster templates
load(file.path(gis_out, "templateRaster.rdata")) ## r, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting
	

## clamp version
rList <- lapply(data.frame(pred.in.cl), function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- x
  tmp
  
})
# plot(tmp)
rStack.cl = stack(rList)
#names(rStack.cl) <- spp.in$best.name
rStack.cl
	
## add auc incidence names to stack
#names(rStack.cl) <- paste0(spp.in$best.name, " ", "auc=", round(spp.in$auc, 2), " ",  "prev=", round(spp.in$incidence,2))
```


```{r TSNE}
## Full data set
Xmat <- pred.in.cl
r <- raster(rStack.cl)
NAs <- indNA
	
# pa version
# Xmat <- (pred.mod[indNA2, ] >= tr)*1
dim(Xmat)
Xmat[1:10, 1:10]
	
perplexity = 50			# 100 what value do you use???
	

# ## Initial PCA dimensions
# # https://towardsdatascience.com/how-to-tune-hyperparameters-of-tsne-7c0596a18868


# Max
(nrow(Xmat) - 1)/3
	
tsne = Rtsne::Rtsne(Xmat, dims = 2, perplexity = perplexity, theta = 0.5, pca = FALSE, num_threads = 0) # don't think I'm using openMP??
	
# system.time(
#   tsne <- Rtsne(Xmat, dims = 2, theta = 0.5, partial_pca = T, num_threads = 0)
# )

# plot(tsne$Y, asp = 1, pch = ".")
# str(tsne, max.level =1)
# plot(tsne$Y, asp = 1)


# library(vegan)

# system.time(
#   pca <- rda(X = Xmat, scale = T)
# )
# 
# str(pca, max.level = 1)
# biplot(pca, pch = ".")
# screeplot(pca)
# round(cumsum(100*pca$CA$eig/sum(pca$CA$eig)),2)[1:10]
# pca$CA$eig[1:15]

# plot(pca$Ybar[,1:2], asp = 1, pch = ".")

## put site scores into raster
makeR <- function(r, siteScores, NAs) {
  
  rSites <- raster(r)
  rSites[] <- NA
  rSites[NAs] <- siteScores
  rSites
  
}
	
rSites1 <- makeR(r, tsne$Y[,1], NAs)
rSites2 <- makeR(r, tsne$Y[,2], NAs)
	
# plot(stack(rSites1, rSites2))
# 
# pcaR1 <- makeR(r, scores(pca, 1, "sites"), NAs)
# pcaR2 <- makeR(r, scores(pca, 2, "sites"), NAs)
# plot(stack(pcaR1, pcaR2))

# plot(stack(rSites1, rSites2,pcaR1, pcaR2))


# ord.stck <- stack(rSites1, rSites2,pcaR1, pcaR2)
# names(ord.stck) <- c("tsne1", "tsne2", "pca1", "pca2")
# 
# 
# pdf(file.path(plotsFolder, "ord_plots.pdf"))
# plot(ord.stck)
# dev.off()


# save(tsne, r, rSites1, rSites2, NAs, file = file.path(plotFolder, 'rdata', "ord_tsne_res_cl_p100.rdata"))
	
save(tsne, r, rSites1, rSites2, NAs, file = file.path(resFolder, "ord_tsne_res_cl_p50.rdata")) # with perp = 50

# 
# save(r, f, indNA2, file = file.path(resFolder, "rast_template_data.rdata"))

load(file.path(plotFolder, 'rdata', "ord_tsne_res_cl_p50.rdata"))
	
writeRaster(rSites1, filename = file.path(plotFolder, 'rdata', "tsne1_nopca_cl_p50.tif"), 
            datatype = "FLT4S", overwrite = T)
writeRaster(rSites2, filename = file.path(plotFolder, 'rdata', "tsne2_nopca_cl_p50.tif"), 
            datatype = "FLT4S", overwrite = T)
# 
	
pdf(file.path(plotFolder, 'plot', "tsne_scatter_cl_p50.pdf"))
plot(tsne$Y, pch = ".")
dev.off()
	
```

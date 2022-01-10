## Make raster from predictions


#### Read data on Ada  #####

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()

## trial run 
```{r setup}
here()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/13_predict_map')
	
pacman::p_load('dplyr', 'rgdal', 'raster','here','glue','corrplot','sjsdm','ggplot2')
	
source(here("source","irrAB.r"))
source(here('..', '12_sjsdm','source', 'corvif-source.r'))
	
packageVersion('sjSDM')
	

```

```{r set-names}
utm10N = 32610
	
period = "S1"
date.model.run = '20210722'			# check!!!
abund = 'pa'		
varsName = 'vars11'		
minocc = 6
cv = '5CV'		
	
outputpath = here('..','..','Output')
	
sjsdmV = '0.1.8'		# check!!!
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	
sppdatapath = here('..','..','format_data','otu')
gispath = here('..','..','format_data','gis')
	
modFolder = file.path(outputpath, "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))		# if not, create!!!
irreFolder = file.path(outputpath, "prediction_map")
	
```


```{r load-data}
# load raster templates - reduced areas
load(file.path(gispath, "templateRaster.rdata")) ## r.msk, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting
# plot(r.msk)
	
# load model result - for species classification
model.train = readRDS(here(modFolder, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_AUC_{date.model.run}.RDS')))
	
# load data for the model ........
load(here(sppdatapath, glue('forbestm_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
# !!! need to change here
if (abund == 'pa')
{
	otu.train = as.data.frame((otu.train>0)*1)
	otu.test = as.data.frame((otu.test>0)*1)
	dim(otu.test)
}
	
str(otu.train) 
names(scale.env.test); formula.env
	
s.otu.train = as.matrix(otu.train)	
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
s.otu.test = as.matrix(otu.test)		
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)
	
# otu.pa.csv, otu.qp.csv

## load species AUC resutls for filtering
load(file.path(modFolder, "sp_test_results.rdata")) # eval.results, sp.res.test, sp.res.train
	
```

```{r irreplace-ana}
########### Filter species for analysis ##########
## Mean AUC per species (and other eval metrics)
str(sp.res.test)
	
table(is.na(sp.res.test$auc))
	
## Filter species by auc
auc.filt = 0.70
table(sp.res.test$auc > auc.filt)
	
# ## extract species over AUC filter
# incidence 
incidence = colSums(otu.train)/nrow(otu.train)
	
# get species names, BOLD ID etc
spp = data.frame(species = colnames(otu.train)) %>%
  tidyr::separate(col = species, into = c("OTU", "empty", "class", "order", "family",
                                          "genus", "epithet", "BOLD", "BOLDID",
                                          "size"),
                  remove = FALSE, sep = "_", convert = TRUE) %>%  ## creates real NAs with convert = T
  mutate(best.name = case_when(is.na(epithet) & is.na(genus) & is.na(family) & is.na(order) ~ class,
                               is.na(epithet) & is.na(genus) & is.na(family) ~ order,
                               is.na(epithet) & is.na(genus) ~ family,
                               is.na(epithet) ~ genus,
                               TRUE ~ paste(genus, epithet, sep = "_")
  )) %>%
  dplyr::select(-empty)%>%
  mutate(auc = sp.res.test$auc,
         incidence = incidence,
         best.name = paste(best.name, BOLDID, sep = "_"))
	
head(spp)
	
sum(is.na(spp$best.name))
sum(grepl("NA_NA", spp$best.name))
head(spp, 30)
	
sum(is.na(spp$family))
	

### Load prediction results
# clamp predictions
load(file.path(irreFolder, 'rdata', paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, "_clamp", ".rdata")))
# pred.mn.cl, pred.sd.cl
	
dim(pred.mn.cl)
	
pred.in.cl = pred.mn.cl[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
dim(pred.in.cl)
	
## get species names 
spp.in = spp[sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc), ]
head(spp.in)
	

pc = .6			# .3 , .6
beta.r.prob.res <- irrAB(x = pred.in.cl, pc = pc, type = "total")
str(beta.r.prob.res)
	
# beta.r.prob_ct <- irrAB(x = x_in_ct, pc = 0.90, type = "total", r = rBlock)
# beta.r.prob_ct.res <- irrAB(x = x_in_ct, pc = 0.90, type = "total")

save(beta.r.prob.res, file = file.path(irreFolder, 'rdata', glue("beta_res_{pc}_noagg.rdata")))
	
beta.pix <- r.msk
beta.pix[indNA] <- beta.r.prob.res$beta

writeRaster(beta.pix, file = file.path(resFolder, "beta_r_prob_noagg.tif"))

```

```{r irre-env-correlation}
# load env data across the study area
load(here(irreFolder, 'rdata', 'newData_unscaled.rdata')) # allVars, newData, indNA,
	
## Final set of VIF chosen predictors
sel.vars11 <- c("gt4_500", "cut_r1k","cut_r250","cut40_r1k","cut40_r250","be30","tri30","Nss30",
            "Ess30","twi30","tpi250","tpi1k","l_rumple","nbr_stdDev_r100","ndvi_p5_r100",
            "ndvi_p5_r500","ndvi_p50_r100","ndvi_p50_r500","ndmi_p95_r100",
            "LC08_045029_20180726_B1","LC08_045029_20180726_B5","lg_DistStream",
            "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m") # insideHJA
	
dim(newData[,sel.vars11]) 
	
if (formula.env=='vars11') {
	evnames = c("Vegetation.4m.r500", "Logged.r1k", "Logged.r250", "Logged40.r1k", "Logged40.r250", "Elevation", 'TRI', 'Northness', 'Eastness', "TWI", 'TPI.r250', 'TPI.1k', 'Canopy.p25', "Rumple", "NBR.sd.r250", "NDVI.p5.r100", "NDVI.p5.r500", "NDVI.p50.r100",  "NDVI.p50.r500", "NDVI.p95.r250", "NDMI.p95.r100", "B1", "B5", 'Stream', 'Road', 'Canopy.2m+', "Canopy.2-4m", "Canopy.4-16m", "HJA" )
}
data.frame(names(scale.env.train), evnames)
	
# load irreplaceability data
pc = .3			# .9 , .3
load(file.path(irreFolder, 'rdata', glue("beta_res_{pc}_noagg.rdata")))
	
pc = .9			# .9 , .3
load(file.path(irreFolder, 'rdata', glue("beta_res_{pc}_noagg.rdata")))
	
#ss = sample(1:length(beta.r.prob.res$beta), 30)
#table(a3[ss]==a6[ss])
	
irre.cor = cor(beta.r.prob.res$beta, newData[,sel.vars11])
str(irre.cor)
	
a = t(irre.cor) %>% data.frame() %>% rename(cor=1) %>% tibble::add_column(name = attr(irre.cor, 'dimnames')[[2]]) %>%
		tibble::add_column(abs.cor = abs(t(irre.cor))) %>%
		arrange(desc(abs.cor)) %>% slice(1:16) %>%
		dplyr::select('cor')
	
# pdf(file.path(irreFolder, 'plot', glue("irreplace_env_{pc}_corplot.pdf")), width = 4, height = 8)
	
corrplot::corrplot(as.matrix(a), 
           is.corr = F, method = "number", cl.pos = "n",
           title = glue("irreplaceability {pc}"), oma = c(0,0,0,0), mar = c(0,0,1,0), 
           addCoef.col = "black", addCoefasPercent = F, number.cex = 0.8)
	
dev.off()
	
sel = 1:6; ss = sample(1:nrow(newData), 8000)
Sparrows = data.frame(irre=beta.r.prob.res$beta[ss], newData[ss,rownames(a)[sel]])
names(Sparrows)[-1] = paste0(names(Sparrows)[-1]," (", round(a$cor[sel],2), ')')
	
# pdf(file.path(irreFolder, 'plot', glue("irreplace_env_{pc}_ggplot.pdf")), width = 8, height = 3)
Sparrows %>%
  tidyr::gather(-irre, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = irre)) +
    facet_wrap(~ var, scales = "free") +
    geom_point(shape=1) +
    stat_smooth()
	
dev.off()
	




pacman::p_load("GGally")
	
ggpairs(Sparrows, lower=list(continuous = wrap("density", alpha = 0.3)), lower='blank')
	


mtcars %>%
  tidyr::gather(-mpg, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = mpg)) +
    facet_wrap(~ var, scales = "free") +
    geom_point() +
    stat_smooth()
	
pairs(Sparrows, verInd = c(1,1,1,1,1), horInd = 1:5,
      lower.panel = panel.cor, 
#      upper.panel=panel.smooth2,diag.panel=panel.hist,
      cex.labels=1.3 )
	
 pairs(~ Fertility + Education + Catholic, data = swiss, row1attop=FALSE,
           subset = Education < 20, main = "Swiss data, Education < 20")
```




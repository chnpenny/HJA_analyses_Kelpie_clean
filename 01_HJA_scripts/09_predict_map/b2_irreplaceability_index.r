## Irreplaceability indices

rm(list=ls())

# HJA_analyses_Kelpie_clean # is the root folder and must have a .Rproj file in it for here::here() to work.
# setwd() # set here to HJA_analyses_Kelpie_clean or use here::here()

## 
# ..... setup
here()
	
pacman::p_load('dplyr', 'terra','here','glue','corrplot','sjSDM','ggplot2')
	
source(here('01_HJA_scripts', '09_predict_map', 'source', 'irrAB.r'))
source(here('01_HJA_scripts', '08_sjsdm','source', 'corvif-source.r'))
	
packageVersion('sjSDM')
	


# ..... set-names
utm10N = 32610
	
period = "S1"
date.model.run = '2024'			# check!!!
abund = 'pa'		
varsName = 'vars11'		
minocc = 6
cv = '5CV'		
	
sjsdmV = '1.0.5'		# check!!!
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	
sppdatapath = here('03_format_data','otu')
gispath = here('03_format_data','gis')
	
resFolder = here::here("04_Output", "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))	
irreFolder = here::here("04_Output", "prediction_map")
predFolder = here::here('04_Output', "sjsdm_prediction_outputs", glue::glue('{varsName}_{date.model.run}'))
	
# ..... load-data

## load species AUC results for filtering
load(file.path(resFolder, "spp_test_data.rdata")) #auc_by_spp, rsq_final, tune.results

# ..... irreplace-analysis
########### Filter species for analysis ##########

## Mean AUC per species (and other eval metrics) from 5CV
str(auc_by_spp)
head(auc_by_spp) # contains OTU name, mean (AUC) and incidence

## Filter species by auc
auc.filt = 0.70
# threshold for presence absence data

table(auc_by_spp$mean >= auc.filt)
	
# ## extract species over AUC filter

# get species names, BOLD ID etc
spp = auc_by_spp %>%
  tidyr::separate(col = OTU, into = c("OTU", "empty", "class", "order", "family",
                                          "genus", "epithet", "BOLD", "BOLDID",
                                          "size"),
                  remove = FALSE, sep = "_", convert = TRUE) %>%  ## creates real NAs with convert = T
  mutate(best.name = case_when(is.na(epithet) & is.na(genus) & is.na(family) & is.na(order) ~ class,
                               is.na(epithet) & is.na(genus) & is.na(family) ~ order,
                               is.na(epithet) & is.na(genus) ~ family,
                               is.na(epithet) ~ genus,
                               TRUE ~ paste(genus, epithet, sep = "_"))
         ,best.name = paste(best.name, BOLDID, sep = "_")
         ) %>%
  dplyr::select(-empty)
	
head(spp)
	
sum(is.na(spp$best.name))
sum(grepl("NA_NA", spp$best.name))
head(spp, 30)
head(spp$best.name)
	
sum(is.na(spp$family))
	

### Load prediction results
# clamp predictions
# load clamp predictions (working folder - too large for github)
load(file.path("./working", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", date.model.run, "_", abund, "_clamp", ".rdata")))
# pred.mn.cl, pred.sd.cl	
dim(pred.mn.cl)
	
pred.in.cl = pred.mn.cl[,auc_by_spp$mean >= auc.filt & !is.na(auc_by_spp$mean)]
dim(pred.in.cl)
	
## get species names 
spp.in <- spp[auc_by_spp$mean >= auc.filt & !is.na(auc_by_spp$mean), ]
head(spp.in)
	

pc = .5			# .3 , .6
beta.r.prob.res <- irrAB(x = pred.in.cl, pc = pc, type = "total")
str(beta.r.prob.res)
	
# beta.r.prob_ct <- irrAB(x = x_in_ct, pc = 0.90, type = "total", r = rBlock)
# beta.r.prob_ct.res <- irrAB(x = x_in_ct, pc = 0.90, type = "total")

# save copy of raster
save(beta.r.prob.res, file = file.path("working", glue("beta_res_{pc}_{varsName}_{date.model.run}_noagg.rdata")))

## load raster templates
gis_out = here('03_format_data','gis')
load(file.path(gis_out, "templateRaster.rdata")) ## r.msk, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting

# convert to (new) terra package format
r.msk
r.msk <- terra::rast(r.msk)	

beta.pix <- r.msk
beta.pix[indNA] <- beta.r.prob.res$beta

plot(beta.pix)

terra::writeRaster(beta.pix, file = file.path(predFolder, paste0("beta_r_prob_noagg_", varsName, "_", date.model.run ,".tif")), overwrite = TRUE)


save(spp, spp.in, file = file.path(resFolder, "spp_test_data_names.rdata"))


# png(filename = file.path("04_Output/figures", "irreplaceability_pc50.png"), units = "mm", height =200, width = 300, res = 200)
# plot(beta.pix)
# dev.off()







# ..... irre-env-correlation
# load env data across the study area

# load(here(irreFolder, 'rdata', 'newData_unscaled.rdata')) # allVars, newData, indNA,
# 	
# ## Final set of VIF chosen predictors
# sel.vars11 <- c("gt4_500", "cut_r1k","cut_r250","cut40_r1k","cut40_r250","be30","tri30","Nss30",
#             "Ess30","twi30","tpi250","tpi1k","l_rumple","nbr_stdDev_r100","ndvi_p5_r100",
#             "ndvi_p5_r500","ndvi_p50_r100","ndvi_p50_r500","ndmi_p95_r100",
#             "LC08_045029_20180726_B1","LC08_045029_20180726_B5","lg_DistStream",
#             "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m") # insideHJA
# 	
# dim(newData[,sel.vars11]) 
# 	
# if (formula.env=='vars11') {
# 	evnames = c("Vegetation.4m.r500", "Logged.r1k", "Logged.r250", "Logged40.r1k", "Logged40.r250", "Elevation", 'TRI', 'Northness', 'Eastness', "TWI", 'TPI.r250', 'TPI.1k', 'Canopy.p25', "Rumple", "NBR.sd.r250", "NDVI.p5.r100", "NDVI.p5.r500", "NDVI.p50.r100",  "NDVI.p50.r500", "NDVI.p95.r250", "NDMI.p95.r100", "B1", "B5", 'Stream', 'Road', 'Canopy.2m+', "Canopy.2-4m", "Canopy.4-16m", "HJA" )
# }
# data.frame(names(scale.env.train), evnames)
# 	
# # load irreplaceability data
# pc = .3			# .9 , .3
# load(file.path(irreFolder, 'rdata', glue("beta_res_{pc}_noagg.rdata")))
# 	
# pc = .9			# .9 , .3
# load(file.path(irreFolder, 'rdata', glue("beta_res_{pc}_noagg.rdata")))
# 	
# #ss = sample(1:length(beta.r.prob.res$beta), 30)
# #table(a3[ss]==a6[ss])
# 	
# irre.cor = cor(beta.r.prob.res$beta, newData[,sel.vars11])
# str(irre.cor)
# 	
# a = t(irre.cor) %>% data.frame() %>% rename(cor=1) %>% tibble::add_column(name = attr(irre.cor, 'dimnames')[[2]]) %>%
# 		tibble::add_column(abs.cor = abs(t(irre.cor))) %>%
# 		arrange(desc(abs.cor)) %>% slice(1:16) %>%
# 		dplyr::select('cor')
# 	
# # pdf(file.path(irreFolder, 'plot', glue("irreplace_env_{pc}_corplot.pdf")), width = 4, height = 8)
# 	
# corrplot::corrplot(as.matrix(a), 
#            is.corr = F, method = "number", cl.pos = "n",
#            title = glue("irreplaceability {pc}"), oma = c(0,0,0,0), mar = c(0,0,1,0), 
#            addCoef.col = "black", addCoefasPercent = F, number.cex = 0.8)
# 	
# dev.off()
# 	
# sel = 1:6; ss = sample(1:nrow(newData), 8000)
# Sparrows = data.frame(irre=beta.r.prob.res$beta[ss], newData[ss,rownames(a)[sel]])
# names(Sparrows)[-1] = paste0(names(Sparrows)[-1]," (", round(a$cor[sel],2), ')')
# 	
# # pdf(file.path(irreFolder, 'plot', glue("irreplace_env_{pc}_ggplot.pdf")), width = 8, height = 3)
# Sparrows %>%
#   tidyr::gather(-irre, key = "var", value = "value") %>% 
#   ggplot(aes(x = value, y = irre)) +
#     facet_wrap(~ var, scales = "free") +
#     geom_point(shape=1) +
#     stat_smooth()
# 	
# dev.off()
# 	
# 
# 
# pacman::p_load("GGally")
# 	
# ggpairs(Sparrows, lower=list(continuous = wrap("density", alpha = 0.3)), lower='blank')
# 	
# pairs(Sparrows, verInd = c(1,1,1,1,1), horInd = 1:5,
#       lower.panel = panel.cor, 
# #      upper.panel=panel.smooth2,diag.panel=panel.hist,
#       cex.labels=1.3 )
# 	
# 
# 
# 
# 

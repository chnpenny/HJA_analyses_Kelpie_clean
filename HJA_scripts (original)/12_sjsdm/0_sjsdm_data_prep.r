## Data set up for sjsdm models - 5 fold cross validation, S1

## 0. Contents ######
## 1. Gets data from github with criteria below
## 2. Creates k folds and sets up tuning grid, and results table for grid runs per k


```{r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm')
	
pacman::p_load('tidyverse','here','conflicted','glue','gridExtra','ggeffects','corrplot') 
# 
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
# [1] ‘0.1.8 2020.12.23
	
source(here('source', 'corvif-source.r'))
source(here('source', 'scale-train-test.r'))
	

```


```{r set-names}
# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"
	
trap <- "M1"; period = "S1"
date.model.run = '20210722'  
abund = 'pa'		
varsName = 'vars11'		
minocc = 6		
k = 5 		# 5-folds
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
otupath = here('..','..','Kelpie_maps',outputidxstatstabulatefolder)
	
eodatapath = here('..','..','format_data','gis')
sppdatapath = here('..','..','format_data','otu')
	
modpath = here('..','..', 'Output', "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
	
sjsdmV = packageVersion('sjSDM') # package version  
	
```


```{r load-data}
# ..... load OTU data ......
otuenv = read.csv(here(otupath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
load(here(eodatapath, 'processed_gis_data', "envVars.rdata")) # allVars
	

```


```{r filter-site&spp-data}
## Select S1 data 
## Split selected data set into Training/Validation and hold out testing set
## Combine M1 and M2 data and split so that sites with both M1 and M2 samples are always in same split

# set seed because randomness will be used below
set.seed(501)
	
# get list of sites with presence of M1 and M2 samples as columns
site.chk = otuenv %>%
  filter(period == "S1") %>%
  select(c(UTM_N, UTM_E, SiteName, trap)) %>%
  tidyr::pivot_wider(names_from = trap, values_from = trap, values_fn = length) %>%
  mutate(numTrap = rowSums(select(., M1, M2), na.rm = TRUE))
site.chk
	
# shuffle order of sites (seed is set globally at start - so can be changed for different runs)
site.chk = site.chk[sample(1:nrow(site.chk)),]
	
# sum numTraps until percent training is reached 
select.percent = 0.75
num.train = round(sum(site.chk$numTrap)*select.percent) # training split is total no of samples * select %
num.train
	
site.chk$cumsum = cumsum(site.chk$numTrap)
train.Names = site.chk$SiteName[which(site.chk$cumsum <= num.train)]
test.Names = site.chk$SiteName[which(site.chk$cumsum > num.train)]
	
# check
length(train.Names); length(test.Names)
sum(site.chk$numTrap[site.chk$SiteName %in% train.Names]); sum(site.chk$numTrap[site.chk$SiteName %in% test.Names])
rm(num.train)
	

########################################################
## Choose species by minocc, taxon, etc #####

## get Species columns by M1 and M2, with minocc calculated per trap 
## choose spp at lease one trap fulfills minocc/trap
spM = otuenv %>% 
  dplyr::filter(period == "S1" & SiteName %in% train.Names) %>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>% # filter by minocc
  ungroup() %>%
  tidyr::pivot_wider(names_from = trap, values_from = nSites, values_fn = function(x) sum(x)>0) %>%
#  filter(M1+M2>=minocc)
#  rowSums(spM1[,2:3])
#  spM1[21,]
  dplyr::select(OTU)
	
spM; nrow(spM)
	

```


```{r CV-folds}
### Make cross validation folds #####

## Create splits for training data into k folds for tuning
# make fold ids - same as above to split training/test
set.seed(10)
	
cv.chk = site.chk %>%
  filter(SiteName %in% train.Names)
	
# shuffle order of sites 
cv.chk = cv.chk[sample(1:nrow(cv.chk)),]
	
# sum numTraps until percent training is reached
out = sum(cv.chk$numTrap)
by = rep(out%/%k, k); by[1:(out%%k)] = sapply(1:(out%%k), function(x) by[x]=by[x]+1)
# num of data points per fold
splits = seq(1, (out+by - out%%by), by = by)
	
cv.chk$cumsum = cumsum(cv.chk$numTrap)
splits = cumsum(by); cc=0
	
cv.chk$cumsum; splits
	
while(sum(cv.chk$cumsum %in% splits)!=k) {
	cv.chk = cv.chk[sample(1:nrow(cv.chk)),]
	cv.chk$cumsum = cumsum(cv.chk$numTrap)
	cc = cc+1; print(cc)
}
	
cv.chk$fold.id = as.numeric(cut(cv.chk$cumsum, breaks = c(1,splits), include.lowest = T, labels = 1:k, right = T))
head(cv.chk)
	
otu.folds = otuenv %>% 
  dplyr::filter(period == "S1" & SiteName %in% train.Names) %>% ## filter for sites in trainig/validation data
  left_join(cv.chk)
head(otu.folds)
	
## Check
table(otu.folds$fold.id)
table(otu.folds[,c("fold.id", "SiteName")])
	
fold.id = otu.folds$fold.id
	
rm(cc, out, by, cv.chk, splits, site.chk)
	
###  filter species 
## training / validation data set
otu.qp.train = select(otu.folds, spM$OTU)## species filter 
	
# convert to presence/absence data
otu.pa.train = otu.qp.train
otu.pa.train[otu.pa.train > 0] = 1
min(colSums(otu.pa.train))  # should be minocc 
	

```


```{r test-spp-data}
### Make Test data set ######
otu.qp.test <- otuenv %>% 
  dplyr::filter(period == "S1" & SiteName %in% test.Names) %>% ## filter for sites in test data
  dplyr::select(spM$OTU) ## Only use species selected for minocc on training
	
## some species with no presences.... 
sum(apply(otu.qp.test, 2, function(x) sum(x > 0))==0)
	
# convert to presence/absence data
otu.pa.test <- otu.qp.test
otu.pa.test[otu.pa.test > 0] = 1
min(colSums(otu.pa.test)) # - some species with no presences.... 
	
# clean up
rm( kelpierundate, minimaprundate, outputidxstatstabulatefolder, primer, samtoolsfilter, samtoolsqual)
	

```


```{r filter-env-data}
### filter(VIF) predictors to data training/test sets ######
str(allVars)       # formatted env vars
	
env.vars = otuenv %>% 
  dplyr::filter(period == "S1") %>%
  dplyr::select(trap, period, UTM_E, UTM_N, SiteName) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_")) %>%
  left_join(y = allVars, by = "SiteName") %>%
  mutate(lg_DistStream = log(DistStream + 0.001),
         lg_DistRoad = log(DistRoad + 0.001),
         lg_cover2m_max = log(l_Cover_2m_max + 0.001),
         lg_cover2m_4m = log(l_Cover_2m_4m + 0.001),
         lg_cover4m_16m = log(l_Cover_4m_16m + 0.001))
	
# str(env.vars)

all.vars = c("ht30", "gt4_r30", "gt4_250", "gt4_500", "cut_r1k", "cut_r500", "cut_r250", "cut40_r1k", "cut40_r500", 
              "cut40_r250", "be30", "tri30","slope30", "Nss30", "Ess30", "twi30", "tpi250", "tpi500", "tpi1k", "l_p25",
              "l_p95", "l_rumple", 'insideHJA', "ndmi_stdDev_r100","ndmi_stdDev_r250", "ndmi_stdDev_r500", "nbr_stdDev_r100",
              "nbr_stdDev_r250", "nbr_stdDev_r500", "ndvi_p5_r100", "ndvi_p5_r250","ndvi_p5_r500", "ndvi_p50_r100",
              "ndvi_p50_r250", "ndvi_p50_r500", "ndvi_p95_r100", "ndvi_p95_r250", "ndvi_p95_r500", "ndmi_p5_r100",
              "ndmi_p5_r250", "ndmi_p5_r500", "ndmi_p50_r100", "ndmi_p50_r250", "ndmi_p50_r500", "ndmi_p95_r100",
              "ndmi_p95_r250", "ndmi_p95_r500","LC08_045029_20180726_B1", "LC08_045029_20180726_B3", 
              "LC08_045029_20180726_B4", "LC08_045029_20180726_B5", "LC08_045029_20180726_B7",
              "LC08_045029_20180726_B10", "lg_DistStream", "lg_DistRoad", "lg_cover2m_max", 
              "lg_cover2m_4m", "lg_cover4m_16m")
# add HJA !!!
	
sum(!complete.cases(env.vars[,all.vars]))
head(env.vars[,all.vars])
	

## Reduce predictos by VIF 
dd = select(env.vars, all_of(all.vars), -'insideHJA') %>% scale() %>% data.frame() %>% add_column(select(env.vars, 'insideHJA'))
# env data need to be scaled first before VIF
varrem = 'a'; maxvif = 100
# force HJA & elevation(be30) in the selection
while  (maxvif>=8 ) {
	if (varrem!='a') { dd = select(dd, -all_of(varrem)) }
	
	vif = corvif(dd)
	vif.count = vif %>% rownames_to_column(var='var') %>% arrange(GVIF)
	vif.count =  vif.count[-c(which(vif.count$var=='be30'), which(vif.count$var=='insideHJA')), ]
	maxvif = max(vif.count$GVIF)
	varrem = vif.count$var[vif.count$GVIF==maxvif]
	
	print(c(ncol(dd), varrem))
}
	
vif
vars = rownames(vif)
varsName
	
rm(dd, varrem, maxvif, vif.count, vif)
	
## separate env.vars into train and non test data sets
env.vars.test = env.vars[env.vars$SiteName %in% test.Names, ]
env.vars.train = env.vars[env.vars$SiteName %in% train.Names, ]
	
# check names 
all(vars %in% colnames(env.vars))
	

```


```{r scale-data}
### Create scaled test & training data set for the best model (after cross-validation)

# select X data
env.train = env.vars.train[,vars]
env.test = env.vars.test[,vars]
	
# select spatial data
XY.train = env.vars.train[ , c("UTM_E", "UTM_N")]
XY.test = env.vars.test[ , c("UTM_E", "UTM_N")]
	

## scale X and spatial data 
a = scale.dd(env.train, env.test)
env.train.scale = a[1]; env.test.scale = a[2]
	
# .. spatial data
a = scale.dd(XY.train, XY.test)
XY.train.scale = a[1]; XY.test.scale = a[2]
	
# .. Y data
a = scale.dd(otu.qp.train, otu.qp.test)
otu.qp.train.scale = a[1]; otu.qp.test.scale = a[2]
	

```

```{r save-data}
# save filtered data
save(spM, test.Names, train.Names, select.percent, vars, varsName, abund, k, minocc, fold.id, env.vars.train, env.vars.test, file = here(sppdatapath, glue('filtered_info_{k}CV_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
save(otu.qp.train.scale, otu.qp.test.scale, otu.pa.train, otu.pa.test, XY.train.scale, XY.test.scale, env.train.scale, env.test.scale, varsName, fold.id,
     file = here(sppdatapath, glue('forbestm_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')) )
	
# scale have to happen only on training (4-folds) when cross-valid
save(otu.pa.train, otu.qp.train, otu.pa.test, otu.qp.test, XY.train, XY.test, env.train, env.test, varsName, fold.id,
     file = here(sppdatapath, glue('fortuning_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')) )
	
	

```

# Dec 23, 2021
# xAI on tuned model 



```{r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm')
	
pacman::p_load('tidyverse','here','conflicted','sjSDM','glue','MetricsWeighted','flashlight','gridExtra','circlize', 'colorspace')
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')                                             
	
here()
packageVersion('sjSDM')
# [1] ‘0.1.8 2020.12.23
	
source(here("source", "xAI-function.r"))
source(here("source", 'sjsdm_function.r'))
	

```

```{r set-names}
trap <- "M1"; period = "S1"
date.model.run = '20210722'
abund = 'pa'
varsName = 'vars11'
minocc = 6
cv = '5CV'; nstep=1000
	
# folder structure 
predpath = here('..','..', 'Output', 'sjsdm_prediction_outputs', glue('{varsName}_{date.model.run}'))
modpath = here('..','..', 'Output', "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
xaipath = here('..','..', 'Output', "xAI", glue('{varsName}_{date.model.run}'))
	
sjsdmV = '0.1.8'		# check!!!
	

```


```{r load-data}
# load model input data 
load(here( 'source', glue('forbestm_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
# load tuned best-params based on all metrics 
tuning = read.table(here(modpath, 'tuning', glue('best_manual_tuning_sjsdm_{sjsdmV}_{varsName}_{cv}_{period}_{abund}_min{minocc}_nSteps{nstep}.csv')), header=T, sep=',') 
	
str(tuning)
	
# data transformation
if (abund == 'pa')
{
	s.otu.train = as.matrix(otu.pa.train)
	s.otu.test = as.matrix(otu.pa.test)
	attr(s.otu.train, 'dimnames') = NULL
	attr(s.otu.test, 'dimnames') = NULL
	str(otu.pa.train)
}
	
str(s.otu.train)
names(env.test.scale)
	

```

```{r load-model}
# best-model results based on which metric
metric='AUC'
select(tuning, contains(metric))
	
model.train = readRDS(here(modpath, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{metric}_{date.model.run}.RDS')) )
	
plot(model.train$history)
	

```


```{r flashlight} 
# xAI
old=Sys.time(); newt = old
# for (i in 1:63 ) {
# for (i in 64:126 ) {
# for (i in 127:ncol(otu.train) ) {
for (i in 1:ncol(otu.pa.train) ) {
	old=Sys.time()
	print(i)
	
	# customized sjsdm prediction function
	custom_predict = function (model1, new_data) {
		newdd = select(new_data, -'otu', -'UTM_E',-'UTM_N')
		spdd = select(new_data, 'UTM_E','UTM_N')
#		print(c(dim(newdd),dim(spdd)))
		apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=spdd, type='link')) , along = -1L), 2:3, mean)[,i]
	}
	
	# variable importance on its own effect
	fl = flashlight::flashlight(model = model.train, data = data.frame(env.train.scale, XY.train.scale, otu=(otu.pa.train[,i]>0)*1), 
				    y = "otu", label = paste(i), predict_function = custom_predict, metrics = list(auc = AUC))
	imp = flashlight::light_importance(fl, m_repetitions = 6, type = "permutation", v = names(env.train.scale), 
				    seed = 55, n_max = nrow(XY.train.scale), lower_is_better = F)
	
	save(imp, file = here(xaipath, 'rdata', glue('fl_min{minocc}_{abund}_{cv}_{varsName}_spp{i}.rdata') ) )
	
	# variable importance on the effect of its interaction
	vv = most_important(imp,10)
	int = flashlight::light_interaction(fl, pairwise = F, type = 'H', v = vv,
				    grid_size = nrow(otu.pa.train), n_max = nrow(otu.pa.train), seed = 42, normalize = F)
	
	save(int, file = here(xaipath,'rdata', 'interaction', 
				     glue('fl-overall-int_min{minocc}_{abund}_{cv}_{varsName}_spp{i}.rdata') ) )
	
	newt = Sys.time()-old
	print(newt)
	rm(imp, int)
}
	

```


```{r plot-individual}
## plot flashlight results of three OTUs randomly

# select individual/interaction effect to plot
mm = c('vint','vind'); 
i=2; mm = mm[i]; mm
	
# choose 3 species based on their prevalence
taxadd = data.frame(sum = colSums(otu.pa.train>0), otu = names(otu.pa.train))
taxadd = taxadd[order(taxadd$sum, decreasing=T),]
taxadd$sum.seq = 1:nrow(taxadd)
str(taxadd)
	
spp.list = taxadd[c(sample(1:10,1), sample.int((nrow(taxadd)*.45):(nrow(taxadd)*.55),1), sample.int((nrow(taxadd)*.9):nrow(taxadd),1 )), ]  
spp.list
	
# load VI data and make plots
plot.list=list()
for (j in 1:3) {
	pp=ggplot()
	ind = which(names(otu.pa.train) == spp.list$otu[j])
	
	if (mm == 'vint') {
		load(here(xaipath,'rdata', 'interaction', 
			      glue('fl-overall-int_min{minocc}_{abund}_{cv}_{varsName}_spp{ind}.rdata') ) )
		pp = plot(int, fill = "darkred") 
	}
	if (mm == 'vind') {
		load(here(xaipath,'rdata', glue('fl_min{minocc}_{abund}_{cv}_{varsName}_spp{ind}.rdata') ))
		pp = plot(imp, fill = "darkred")
	}
	
	name = toString(strsplit(strsplit(as.character(spp.list$otu[j]),'__')[[1]][2],'_')[[1]][2:6])	# taxon info of OTU
	plot.list[[j]] = pp + ggtitle(paste0(spp.list$sum[j],', ', name)) + 
					 theme(plot.title = element_text(size = 9))
}
	 

# individual/interaction effect
if (mm=='vind') {
	pdf(here(xaipath, 'plot', glue('flashlight-{metric}_{varsName}_{abund}_{cv}_{period}_min{minocc}_{date.model.run}.pdf')), width=12, height=7)
} else if (mm=='vint') {
	pdf(here(xaipath, 'plot', glue('fl-overall-int_{varsName}_{abund}_{cv}_{period}_min{minocc}_{date.model.run}.pdf')), width=12, height=7)
}
	
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]], nrow=1)
	
dev.off()
	
```


```{r table-allspp}
## make table for plotting VI of all species

# choose interaction/individual effect
mm = c('vint','vind'); name='a'; 
i = 1; mm = mm[i]; mm
	
varx = data.frame(vardd = c(names(env.train.scale),names(XY.train.scale)), 
				  ind = 1:(ncol(env.train.scale) + ncol(XY.train.scale)) )
impdd = data.frame(out = character(), var.abs = character(), varx.abs = numeric(), value.abs = numeric(), 
				   var.pos = character(), varx.pos = numeric(), value.pos = numeric() )
str(impdd)
	
# make & save the table
for (i in 1:ncol(otu.pa.train)) {
	print(i)
	xpos = 'a'; xabs = 'a'; idd = 'a'
	
	if (mm == 'vint') {
		load(here(xaipath,'rdata', 'interaction', glue( 
					'fl-overall-int_min{minocc}_{abund}_{cv}_{varsName}_spp{i}.rdata') ) )
		loadd = int; rm(int)
		
		n1=1; iabs = which(abs(loadd$data$value)==max(abs(loadd$data$value))); xabs = loadd$data$variable[iabs]
		n1=1; xpos = most_important(loadd, n1)[n1]
		
		idd = data.frame(otu = names(otu.pa.train)[i], var.abs = xabs, 
						 varx.abs = varx$ind[varx$vardd==xabs], 
						 value.abs = loadd$data$value[loadd$data$variable==xabs], 
						 var.pos = xpos, varx.pos = varx$ind[varx$vardd==xpos], 
						 value.pos = loadd$data$value[loadd$data$variable==xpos] )
		name = glue( 'overall-int-var-fl_min{minocc}_{abund}_{cv}_{varsName}.csv' )
	}
	
	if (mm == 'vind') {
		load(here(xaipath,'rdata', glue( 'fl_min{minocc}_{abund}_{cv}_{varsName}_spp{i}.rdata') ))
		loadd = imp; rm(imp)
		
		n1=1; xabs = which(abs(loadd$data$value)==max(abs(loadd$data$value))); vabs = loadd$data$variable[xabs]
		n1=1; xpos = which(loadd$data$value==max(loadd$data$value)); vpos = loadd$data$variable[xpos]
		
		idd = data.frame(otu = names(otu.pa.train)[i], var.abs = loadd$data$variable[xabs], 
						 varx.abs = varx$ind[varx$vardd==loadd$data$variable[xabs]], 
						 value.abs = loadd$data$value[xabs], var.pos = loadd$data$variable[xpos],
						 varx.pos = varx$ind[varx$vardd==loadd$data$variable[xpos]], 
						 value.pos = loadd$data$value[xpos] )
		name = glue( 'imp-var-fl_min{minocc}_{abund}_{cv}_{varsName}.csv' )
	}
	
	impdd = rbind(impdd, idd)
	rm(idd)
	
	write.table(impdd, file=here(xaipath, 'rdata', 'sum-flash', name), row.names=F,sep=',')
	
}
	
str(impdd)
	
# load the saved table in cases that the table has been made
if (mm == 'vind') { name = glue( 'imp-var-fl_min{minocc}_{abund}_{cv}_{varsName}.csv' ) }
if (mm == 'vint') { name = glue( 'overall-int-var-fl_min{minocc}_{abund}_{cv}_{varsName}.csv' ) }
name
	
impdd = read.table(here(xaipath, 'rdata', 'sum-flash', name), header=T,sep=',')
	

```


```{r plot-circular}
## plot VI for all OTUs as a circular plot

# check info
table(impdd$otu == names(otu.pa.train))
mm				# vint, vind
	
# variables index & value which has biggest coefficient
effect_comb2 = data.frame(max_effects = impdd$varx.pos, V2 = impdd$value.pos); pp='pos'
str(effect_comb2)
	
# check incidence of some variables
table(effect_comb2$max_effects==6); table(effect_comb2$max_effects==23)
	
# some text to plot
otu.text = glue("incidence ({abund})")
version.text = ''
	 
if (formula.env=='vars11') {
	evnames = c("Vegetation.4m.r500", "Logged.r1k", "Logged.r250", "Logged40.r1k", "Logged40.r250", "Elevation", 'TRI', 'Northness', 'Eastness', "TWI", 'TPI.r250', 'TPI.1k', 'Canopy.p25', "Rumple", "NBR.sd.r250", "NDVI.p5.r100", "NDVI.p5.r500", "NDVI.p50.r100",  "NDVI.p50.r500", "NDVI.p95.r250", "NDMI.p95.r100", "B1", "B5", 'Stream', 'Road', 'Canopy.2m+', "Canopy.2-4m", "Canopy.4-16m", "HJA" )
	evn.group = c('Lidar-canopy', 'Anthropogenic', 'Anthropogenic', 'Anthropogenic', 'Anthropogenic', 'Topography', 'Topography', 'Topography', 'Topography', 'Topography', 'Topography', 'Topography', 'Lidar-canopy', 'Lidar-canopy', 'Landsat-annual', 'Landsat-annual', 'Landsat-annual', 'Landsat-annual', 'Landsat-annual', 'Landsat-annual', 'Landsat-annual', 'Landsat-bands', 'Landsat-bands', 'Topography', 'Anthropogenic', 'Lidar-canopy', 'Lidar-canopy', 'Lidar-canopy', 'Anthropogenic')
#	table(evn.group)
}
	
evnames = paste0(seq(1,length(evnames)),'=',evnames)
	
# plotting
if (mm=='vint') {cc = 'overall-int'} else if (mm=='vind' & pp=='abs') {cc = 'ind'}
	
# pdf(here(xaipath, 'plot', glue('var-imp-{cc}-flashlight-train_{varsName}_{abund}_{cv}_{period}_min{minocc}_{date.model.run}.pdf')), height=10, width=12)
	
# graph setting
layout(mat = matrix(c(1,2), nrow = 1, ncol = 2), widths = c(.815,.185) )
par(cex = .78, mar = c(2,2,1,.05))
seeds = c(0,0.77)
	
# plotting the circle (function from 'source')
cov.circle.env(version.text, evnames, evnp = evn.group, otu.text, effect_comb = effect_comb2, otu.tbl = otu.pa.train, seeds) 
	
# ...... plotting the legend 
cname = c("blue4","firebrick4","gold3","green4","magenta4")
cols = vector()
for (i in 1:length(unique(evn.group))) {
	colsi = lighten(cname[i], seq(seeds[1], seeds[2],length = table(evn.group)[i]))
	cols = append(cols, colsi)
	rm(colsi)
}
rm(cname)
	
colorord = data.frame(nn = 1:length(evnames), gg = evn.group)
colorord = colorord[order(colorord$gg),]
colorord = cbind(colorord, cc = cols)
colorord = colorord[order(colorord$nn),]
cols = as.character(colorord$cc)
rm(colorord)
	
aa = data.frame(evnames = evnames, se = 1:length(evnames))
	
show = sort(unique(effect_comb2$max_effects))
legenddd = data.frame(name = character(), min = numeric(), max = numeric() )
for (i in show) {
	aa = effect_comb2$V2[effect_comb2$max_effects==i]
	a = data.frame(name = i, min = round(min(aa),2), max = round(max(aa),2))
	legenddd = rbind(legenddd, a)
	rm(aa,a)
	}
	
plot(rep(-.1,length(show)), seq(2.65,3.4*length(show),3.4), col=cols[show], pch=18, xlab='', ylab='', axes=F, xlim=c(-.1,1.1), cex=2.2) #, ylim=c(3,85.2)
text(rep(-0.07,length(show)), seq(2.8,3.4*length(show),3.4), labels=sapply(show, function(aa) paste0(str_split(evnames[aa],'=')[[1]][1],' ', str_split(evnames[aa],'=')[[1]][2])), cex=1.1, pos=4)
	
text(rep(.4,length(show)), seq(1.25,3.4*length(show),3.4), labels=sapply(1:length(show), function(aa) paste0(sprintf('%.2f',legenddd$min[aa]),', ', sprintf('%.2f',legenddd$max[aa]))), cex=.9, pos=2,offset=-.5 )
	
dev.off()
	

```







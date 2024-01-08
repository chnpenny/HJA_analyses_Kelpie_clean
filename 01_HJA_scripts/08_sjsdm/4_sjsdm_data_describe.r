# Dec 23, 2021
# some descriptive plots

```{r setup}
rm(list=ls())
pacman::p_load('tidyverse', 'here', 'conflicted', 'glue', 'gridExtra', 'metacoder', 'sf', 'sjmisc') 
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
	

```


```{r set-names}
# basic setting
period = "S1"
date.model.run = '2023'
abund = 'pa'
varsName = 'vars11'
minocc = 6
cv = '5CV'; nstep = 1000
	
# ....... folder structure .......
eodatapath = here('03_format_data','gis')
sppdatapath = here('03_format_data','otu')
	
# bioinfo structure
outputidxstatstabulatefolder = glue("outputs_minimap2_20200929_F2308_q48_kelpie20200927_BF3BR2_vsearch97")
otupath = here('02_Kelpie_maps', outputidxstatstabulatefolder)
rm(outputidxstatstabulatefolder)
	
predpath = here('04_Output', 'sjsdm_prediction_outputs', glue('{varsName}_{date.model.run}'))
modpath = here('04_Output', "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
xaipath = here( '04_Output', "xAI", glue('{varsName}_{date.model.run}'))
despath = here(modpath, '..', 'descriptive')
	
sjsdmV = packageVersion('sjSDM')
	

```


```{r load-data}
# model input data
load(here(sppdatapath, glue('fortuning_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
# model setting info
load(here(sppdatapath, glue('filtered_info_{cv}_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
# ..... load OTU data ......
alldata = read.csv(here(otupath, 'sample_by_species_table_F2308_minimap2_20200929_kelpie20200927_FSL_qp.csv'), header = T, sep = ',', stringsAsFactors = F, na.strings = 'NA')
dim(alldata)
names(alldata)[110:125]
	
alldata1 = alldata %>% filter(period == period[1]) 
dim(alldata1); unique(alldata1$period)
	
# ... load env data ................................
load(here(eodatapath, 'processed_gis_data','envVars.rdata'))	
sort(names(allVars))
# 67
	
# ..... data seperate .....
load( here(sppdatapath, glue('modelData_{abund}.rdata')) )	
	
```


```{r heat-map}
## heat map for OTUs

# .... format data for heat map ....
(names(otu.pa.train) == names(otu.pa.test)) %>% table
	
# from otu.train
sepe = 'no'		# 'tt' train/test or 'm' m1/m2
m1 = bind_rows((otu.pa.train %>% rownames_to_column(var = 'sample.id') %>% add_column(., id = paste0(sepe,1))), 
		(otu.pa.test %>% rownames_to_column(var = 'sample.id') %>% add_column(., id = paste0(sepe,2)))) %>% 
		mutate(idfull = paste0(sample.id,'_', id))
	

# ...... after creating m1
m1 = select(m1, contains('__'), idfull) %>% 
			gather(name, presence, -idfull) %>% spread(idfull, presence) %>% 
			separate(name, c('dummy','taxa'), sep = '__', remove=F) %>% 
			separate(taxa, c('taxa','dummy'), sep = '_BOLD') %>% 
			separate(taxa, c('class','order','family','genus','species'), sep = '_') %>% select(., -dummy)
	
m1[1:5,1:4]
	

# check which taxon info has NA
table(m1$order=='NA')
table(m1$family=='NA')
	

# ............. modify NA .........................
m1$species = na_if(m1$species, '')
m1$genus = na_if(m1$genus, '')
m1$family = na_if(m1$family, '')
	
# replace NA with letter+number
a = length(which(m1$species == 'NA'))
m1$species[m1$species=='NA'] = paste0('.s',seq(1, a))
	
b = length(which(m1$genus=='NA'))
m1$genus[m1$genus=='NA'] = paste0('.g', seq((1), (b)))
	
c = length(which(m1$family=='NA'))
m1$family[m1$family=='NA'] = paste0('.f', seq((1), (c)))
	
rm(a,b,c)
m1[107,]
	
# ..... steps for heat map .....
m1h = parse_tax_data(m1, class_cols = c('class','order','family','genus','species'), named_by_rank = TRUE)
names(m1h)
m1h %>% n_obs
	
m1h$data$otu_table = calc_obs_props(m1h, data = "tax_data", cols = names(select(m1, contains(paste0('_',sepe)))) )
m1h$data$tax_table = calc_taxon_abund(m1h, data = "otu_table", cols = names(select(m1, contains(paste0('_',sepe)))) )
	
m1h$data$diff_table = compare_groups(m1h, data='tax_table', cols=names(select(m1, contains(paste0('_',sepe)))), groups=as.character(sapply(names(select(m1, contains(paste0('_',sepe)))), function(a) str_split(a, '_')[[1]][2])) )
unique(m1h$data$diff_table$taxon_id)
m1h$data$tax_data
	
m1h$data$type_abund = calc_group_mean(m1h, "tax_table", cols=names(select(m1, contains(paste0('_',sepe)))), groups=as.character(sapply(names(select(m1, contains(paste0('_',sepe)))), function(a) str_split(a, '_')[[1]][2])))
	
# choose to plot till family-level or species-level
tdraw = 'allnum'
if ( tdraw == 'family') { nodraw = c("species", "genus") }
if ( tdraw == 'allnum' ) { nodraw = ''}
	
pname = glue('heat-map_{tdraw}_{varsName}_{sepe}_{abund}_min{minocc}.pdf')
# 'allnum' plot is supplement figure 'S-HEATTREE'
# 'family' plot is panel B of the figure 'STUDYAREA'	
	
	
labtext = 'count'
	
set.seed(11)		# make sure the graph has the same layout 
m1h %>% filter_taxa(! taxon_ranks %in% nodraw) %>%
heat_tree( node_label = taxon_names,
            node_size = n_obs, 
            node_color = n_obs,
            node_size_axis_label = NULL, 
            make_edge_legend = T, node_color_axis_label = labtext, 
			overlap_avoidance = 2,
            margin_size = c(0.03, 0.03, 0.03, 0.03)
            , node_label_size_range = c(0.008,0.04)
#           , output_file = here(despath, pname)
) 
	

```


```{r descriptive-graph}
## basic descriptive graphs train/test

# make table for plotting
str(env.vars.train); train.Names
exp = rbind(env.vars.train, env.vars.test) %>% select('UTM_E', 'UTM_N', 'trap', 'SiteName') %>% 
		add_column(assign = c(rep('full training', nrow(env.vars.train)), rep('test', nrow(env.vars.test))) )
	
# check if assign is correct
sort(unique(exp$SiteName[exp$assign == 'test'])) == sort(test.Names)
sort(unique(exp$SiteName[exp$assign == 'full training'])) == sort(train.Names)
	

# bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
utm10N = 32610
hja = st_read(here(eodatapath, 'raw_gis_data', 'shape', "HJA_Boundary.shp"))
	
hja_bound = subset(hja, FP_NAME == "H.J. Andrew Experimental Forest") # just get the boundary
hja.utm = st_transform(hja_bound, crs = utm10N) # transform to UTM
hjautm = as.data.frame(st_coordinates(hja.utm))
str(hjautm)
	

# plotting ......................
hjashape = geom_polygon(data = hjautm, aes(x=X, y=Y), colour = "black", fill = NA)
	
p1 = ggplot(exp, aes(UTM_E, UTM_N),shape=assign) + theme(panel.background = element_rect(fill = 'floralwhite'), legend.key = element_rect(fill = "floralwhite"), legend.text = element_text(size = 11), legend.title = element_text(size = 11)) + geom_point(aes(colour=assign,shape=assign), size = 2) + scale_shape_manual(values=c(19,8)) + scale_colour_manual(values=c('lightgreen','chocolate1'))
	
# pdf(here(despath, glue('random-test-train_s1_min{minocc}.pdf')), height=6,width=6.5)
	
p1 + hjashape
	
dev.off()
	

```



```{r simple-ana-model}
# load model input data
load(here(sppdatapath, glue('forbestm_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	

## make table for plotting
# make taxonomy table
taxadd = data.frame(sum = colSums(otu.pa.train>0), otu = names(otu.pa.train)) %>% 
		 arrange(desc(sum)) %>% mutate(sum.seq = 1:ncol(otu.pa.train))
glimpse(taxadd)
	
# add taxonomic information 
taxadd$order = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[2])
taxadd$class = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[1])
taxadd$family = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[3])
	
taxadd = taxadd %>% left_join(., (taxadd %>% count(order)), by = 'order') %>% rename(sum.order = n)
glimpse(taxadd)
	

## add prediction(auc) data
auc.all = data.frame(otu = as.character(names(otu.pa.train)), 
					 auc.test = 0.1, 
					 auc.exp = 0.1 )
str(auc.all)
	
formula1 = varsName; formula1
	
mm = 'AUC' # choose which metric to be used
for (j in c('explain', 'test')) {
	set = j
	
		# load AUC results
		roc.dd = readRDS(here(predpath, 'rdata', glue('roc_result_{set}_{period}_{abund}_{cv}_min{minocc}_{formula1}_{mm}_{date.model.run}.RDS'))) 
		 
		# ... make long table
		if (j=='explain') {
			auc.te = data.frame( auc.exp = roc.dd$roc.allS, otu = names(roc.dd$otu) )
#			str(auc.te)
		} else if (j=='test') {
			auc.te = data.frame( auc.test = roc.dd$roc.allS, otu = names(roc.dd$otu) ) }
		
		auc.all = left_join(auc.all, auc.te, by=c('otu'='otu'), suffix=c('', glue('.{mm}')), copy=T)
} 
	
# .. after all metrics are added, sort the table
auc.all = dplyr::select(auc.all, -'auc.test',-'auc.exp')
glimpse(auc.all)
	
# ... join with taxonomy table
auc.all = auc.all %>% left_join(., select(taxadd, 'otu','order','class','family','sum','sum.order'), by = 'otu')
abc = data.frame(seq.order = letters[1:length(unique(taxadd$sum.order))], 
				 order = sort(unique(taxadd$sum.order), decreasing = T))
auc.all$oOrder = sapply(1:dim(auc.all)[1], function(x) paste(abc$seq.order[abc$order==auc.all$n[x]], '.',auc.all$order[x], '.', auc.all$n[x], sep = ''))
	

## plot of incidence vs AUC 
p.test = round(summary(lm(auc.test.AUC ~ sum, data = auc.all))$coefficients[2,4], 2)
p.exp = round(summary(lm(auc.exp.AUC ~ sum, data = auc.all))$coefficients[2,4], 2)
	
plt = ggplot(auc.all, aes(x=sum)) + geom_point(aes(y=auc.test.AUC, colour='test')) + 
#	  geom_smooth(aes(y=auc.test.AUC), method = 'lm', se = T, colour='#56B4E9') + 
	  geom_point(aes(y=auc.exp.AUC, colour='training')) + 
#	  geom_smooth(aes(y=auc.exp.AUC), method = 'lm', se = T, colour='#E69F00') + 
#	  scale_color_manual(values=c('#E69F00', '#56B4E9'), breaks = c('training', 'test'), labels = c(c(glue('training, p: {p.exp}'), glue('test, p: {p.test}')) ), name = NULL) + 
	  scale_color_manual(values=c('#E69F00', '#56B4E9'), breaks = c('training', 'test'), labels = c(c(glue('training'), glue('test')) ), name = NULL) + 
	  theme(legend.text = element_text(size = 8), legend.position = c(.91, .079)) + 
	  labs(x = 'incidence', y = 'AUC')
	

## plot of family vs AUC(test)
dd = auc.all[complete.cases(auc.all$auc.test.AUC),] %>% rename("Order"=4, 'aucT'=3)
dd$family[dd$family=='NA'] = paste0('NA',seq(1, length(which(dd$family=='NA'))))
dd$family = with(dd, reorder(family, sum.order, max))
dd$Order = with(dd, reorder(Order, sum.order, max))
p.family = round(anova(lm(aucT~family, data = dd))[1,5], 2)
	
plt2 = ggplot(dd, aes(x=family, y=aucT)) + geom_point(aes(group=Order, colour=family)) +
	   facet_grid(.~Order, scales='free') +
	   ylab('AUC (test)') + xlab(glue('family (p: {p.family})')) + 
#	   scale_shape_manual(values = seq(0,n_distinct(dd$Order))) + 
	   theme(axis.text.x=element_blank(), legend.position='top') 
#	   scale_x_discrete(labels=seq(1,n_distinct(dd$family))) + 
	

# pdf(here(predpath, 'plot', glue('family-auc_{abund}_{period}_min{minocc}_{varsName}_tuned-{mm}_{date.model.run}.pdf')), width=13, height=10)
# this is the supplement plot Figure S-AucFamily
	
plt2
	
# pdf(here(predpath, 'plot', glue('lm-incidence-auc_{abund}_{period}_min{minocc}_{varsName}_tuned-{mm}_{date.model.run}.pdf')), width=7, height=5)
# this is the supplement plot Figure S-AucIncidence
	
plt
	
dev.off()
	

```

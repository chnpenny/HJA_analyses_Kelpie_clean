# analyze final (tuned) model (tuning in cluster)



```{r setup}
rm(list=ls())
pacman::p_load('tidyverse', 'here', 'conflicted', 'sjSDM', 'glue', 'pROC', 'gridExtra') 	
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
# [1] ‘1.0.5’ 2023.07.25
	

```


```{r set-names}
period = "S1"
date.model.run = '2023'
abund = 'pa'
varsName = 'vars11'
minocc = 6
cv = '5CV'; nstep=1000
	
# ....... folder structure .......
predpath = here( '04_Output', 'sjsdm_prediction_outputs', glue('{varsName}_{date.model.run}'))
modpath = here('04_Output', "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
sppdatapath = here('03_format_data','otu')
	
sjsdmV = packageVersion('sjSDM')
	

```


```{r load-data}
abund; varsName; minocc; cv
	
# data for the model 
load(here(sppdatapath, glue('forbestm_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
if (abund == 'pa')
{
	s.otu.train = as.matrix(otu.pa.train) %>% unname
	s.otu.test = as.matrix(otu.pa.test) %>% unname
#	attr(s.otu.train, 'dimnames') 
	str(s.otu.train)
}
	
names(env.test.scale)
	
# tuned results 
tuning.dd = read.table(here(modpath, 'tuning', glue('manual_tuning_sjsdm_{sjsdmV}_{varsName}_{cv}_{period}_meanEVAL_{abund}_min{minocc}_nSteps{nstep}.csv')), header = T, sep = ',')
	


```


```{r analyze-tune}
### analyze tuing result & run best model with CPU
str(tuning.dd)
	
# some metrics may select the same hyper-param combination
which.max(tuning.dd$AUC.valid_mean); which.max(tuning.dd$cor.valid_mean)
which.max(tuning.dd$nagel.valid_mean); which.max(tuning.dd$ll.valid_mean)
which.max(tuning.dd$plr.valid_mean)
	
maxdd = select(tuning.dd, ends_with('.valid_mean'), -starts_with('tjur.'))
pp = sapply(1:ncol(maxdd), function(i) which.max(maxdd[,i])) 
# index of selected rows (best params based on different metrics)
	
maxdd[pp,]
	

## plot the tuning result 
# pdf(here(modpath, 'plot', glue('plot_tuning_sjsdm_{sjsdmV}_{period}_{cv}_{abund}_min{minocc}_{varsName}_{date.model.run}.pdf')), width = 8, height = 4.5)
	
# pdf(here(modpath, 'plot', glue('TSS-plot_tuning_sjsdm_{sjsdmV}_{period}_{cv}_{abund}_min{minocc}_{varsName}_{date.model.run}.pdf')), width = 8, height = 4.5)
	
par(cex = .75, las = 1)
# plot all the tuning grid
xcol = 'll.train_mean'
	
plot(y = tuning.dd$AUC.valid_mean, x = tuning.dd[,xcol], 
     pch = 3, col = 'blue', ylim = c(0.5,1), xlab = '', ylab = '')
points(y = tuning.dd$AUC.train_mean, x = tuning.dd[,xcol], pch = 8)
	
title(ylab = 'AUC', line = 2.6, xlab = 'log-likelihood (train)')
	
# plot best combination based on AUC and other metrics
points(y = tuning.dd[pp[1], 'AUC.valid_mean'], x = tuning.dd[pp[1],xcol], pch = 20, col = 'red')
points(y = tuning.dd[pp[2:5], 'AUC.valid_mean'], x = tuning.dd[pp[2:5],xcol], pch = 20, col = 'gray')
	
points(y = tuning.dd[pp[1], 'AUC.train_mean'], x = tuning.dd[pp[1], xcol], pch = 8, col = 'red')
points(y = tuning.dd[pp[2:5], 'AUC.train_mean'], x = tuning.dd[pp[2:5], xcol], pch = 8, col = 'gray')
	
# add vertical lines for best combinations 
abline(v = tuning.dd[pp[1],xcol], lty = 2, col = 'red')
abline(v = tuning.dd[pp[2:5],xcol], lty = 2, col = 'gray')
	
# add text for the best combinations
# need to change after re-run!!!!!
text(x = tuning.dd[pp[2:3], xcol]*1.005, y = c(.97,1), c('log-likelihood', expression("Nagelkerke's R"^2)), cex = .92)
text(x = tuning.dd[pp[4], xcol]*.94, y = 1, 'positive likelihood rate', cex = .92)
text(x = tuning.dd[pp[1], xcol]*.99, y = 1, 'AUC', cex = .92)
text(x = tuning.dd[pp[5], xcol]*.99, y = .968, 'correlation', cex = .92)
# text(x = tuning.dd[pp[6], xcol]*.99, y = .968, 'TSS', cex = .92)
	
legend('topright', pch = c(8,3), col = c('black','blue'), c('auc.train','auc.valid'), bty = 'n')
	
dev.off()
	

## model with optimal parameters
# load selected hyper-params
lambda.env = tuning.dd[pp, 'lambda.env']
alpha.env = tuning.dd[pp, 'alpha.env']
lambda.sp = tuning.dd[pp, 'lambda.sp']
alpha.sp =  tuning.dd[pp, 'alpha.sp']
lambda.bio = tuning.dd[pp, 'lambda.bio']
alpha.bio =  tuning.dd[pp, 'alpha.bio']
hidden = list(c(50L, 50L, 10L), c(25L, 25L, 10L))
hidden1 = tuning.dd[pp, 'hidden.ind']
acti = as.character(tuning.dd[pp, 'acti.sp'])
drop = tuning.dd[pp, 'drop']
lrn = tuning.dd[pp,'lr']
	
if(abund == "pa") {
	famn = stats::binomial('probit')
	samn = 5000L; devn = 'cpu'
	itern = 170L } else stop("check abund")
	
# double-check the attributes are removed from model input data
str(env.test.scale); str(s.otu.train)
	
# run best models with cpu based on all metrics
for (i in unique(match(pp, unique(pp))) ) {
#	i = pp[1]
	model.train = sjSDM(Y = s.otu.train,
	  env = DNN(data = env.train.scale, formula = ~.,
	  hidden = hidden[[hidden1[i]]], lambda = lambda.env[i], alpha = alpha.env[i], activation = acti[i], dropout = drop[i], bias = T),
	  
	  biotic = bioticStruct(lambda = lambda.bio[i], alpha = alpha.bio[i], on_diag = F, inverse = F),
	  
	  spatial = linear(data = XY.train.scale, ~0 + UTM_E*UTM_N, lambda = lambda.sp[i], alpha = alpha.sp[i]),
	  
	  device = devn, learning_rate = lrn[i], 
	  step_size = NULL, iter = itern, family = famn, sampling = samn 
	) 
	mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
#	saveRDS(model.train, here(modpath, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')) )
	rm(mm, model.train)
}
	
for (i in unique(match(pp, unique(pp))) ) {
	print(i)
	mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
	model.train = readRDS(here(modpath, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')) )
	
#	pdf(here(modpath, 'plot', glue('model-history_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.pdf')), width = 5, height = 5)
	
	par(cex = .8)
	plot(model.train$history)
	mtext(paste0(names(maxdd)[i],': ', round(select(tuning.dd[pp,], names(maxdd)[i])[i,],3)), side = 3, line = -1, adj = .95)
	mtext(paste0(names(tuning.dd)[c(2:7,9:10)], ': ', select(tuning.dd[pp,], names(tuning.dd)[c(2:7,9:10)])[i,]), side = 3, line = seq(-2,-9), adj = .95)
	
	model.train$history
	
	dev.off()
	rm(mm, model.train)
}
	
# write.table(tuning.dd[pp,], here(modpath, 'tuning', glue('best_manual_tuning_sjsdm_{sjsdmV}_{varsName}_{cv}_{period}_{abund}_min{minocc}_nSteps{nstep}.csv')), row.names = F, sep = ',')
	

```


```{r sjsdm-linear}
## !!! not necessary for the following chunks !!!
# this is the linear sjsdm of the best hyperparameters combination
# for discussion part of the manuscript

i = 1 			# for best AUC
model.train = sjSDM(Y = s.otu.train,
	  env = linear(data = env.train.scale, formula = ~.,
	  lambda = lambda.env[i], alpha = alpha.env[i] ),
	  
	  biotic = bioticStruct(lambda = lambda.bio[i], alpha = alpha.bio[i], on_diag = F, inverse = F),
	  
	  spatial = linear(data = XY.train.scale, ~0 + UTM_E*UTM_N, lambda = lambda.sp[i], alpha = alpha.sp[i]),
	  
	  device = devn, learning_rate = lrn[i], 
	  step_size = NULL, iter = itern, family = famn, sampling = samn 
	) 
mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
	
# saveRDS(model.train, here(modpath, glue('s-jSDM_linear.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')) )
	

## produce & save predicted occurrence probability of each OTU based on the linear best model
for (pred in c('explain', 'test')) {
	# explanatory AUC
	if (pred == 'explain') {
		newdd = NULL
		newsp = NULL
		otudd = s.otu.train
		set = 'explain'
	}
	# predictive AUC
	if (pred == 'test') {
		newdd = env.test.scale
		newsp = XY.test.scale
		otudd = s.otu.test
		set = 'test'
	}
	
		model1 = model.train
		otudd1 = otudd
		
		pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata = newdd, SP = newsp, type = 'link')) , along = -1L), 2:3, mean) %>% unname
		
		otudd1 = data.frame(otudd1)
		names(otudd1) = names(otu.pa.train)
		
		otudd1 = rbind(otudd1, count = (base::colSums(otudd1)>0 )*1 )
		
		# some OTUs don't occur in test set
		pred.dd = pred.dd[ , which(otudd1[nrow(otudd1),] == 1)]
		
		otudd1 = otudd1[1:nrow(pred.dd), which(otudd1[nrow(otudd1), ] == 1)]
		
		otudd.pa = (otudd1>0)*1
		
		# .... calculate AUC
		roc.dd = sapply(1:ncol(otudd1), function(j) as.numeric(pROC::roc( response = otudd.pa[,j], predictor = pred.dd[,j], direction = '<', quiet = T)$auc))
		
		auc.mean = mean(roc.dd)
		
#		saveRDS(list(pred.Y = pred.dd, otu = otudd1, roc.allS = roc.dd, auc.mean = auc.mean), here(predpath, 'rdata', glue('roc_result_{set}_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}-linear_{date.model.run}.RDS')))
		
		rm(pred.dd, model1, roc.dd, auc.mean, otudd1, otudd.pa)
}
	

## plot auc.test ~ auc.train, color ~ order, circle size ~ incidence 
# linear v.s DNN

setS = c(glue(mm,'-linear'), mm)
dd.ggplot = vector(mode = "list", length = length(setS))
plot.list = list()
	
## make taxonomy table
taxadd = data.frame(sum = colSums(otu.pa.train>0), otu = names(otu.pa.train)) %>% arrange(desc(sum)) %>%
					mutate(sum.seq = 1:ncol(otu.pa.train)) 
str(taxadd)
	
# add taxonomic information 
taxadd$class = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[1])
taxadd$order = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[2])
taxadd$family = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[3])
	
taxadd = left_join(taxadd, (taxadd %>% count(order)), by = 'order') %>% rename(sum.order = n)
	
## add prediction(auc) data
auc.all = data.frame(otu = as.character(names(otu.pa.train)), auc.test = 0.1, auc.exp = 0.1 )
str(auc.all)
	
# load saved prediction data
for (j in c('explain', 'test')) {
	for (i in setS) {
	set = j
		mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
		
		roc.dd = readRDS(here(predpath,'rdata', glue('roc_result_{set}_{period}_{abund}_{cv}_min{minocc}_{varsName}_{i}_{date.model.run}.RDS'))) 
		
		# ... make long table
		if (j=='explain') {
			auc.te = data.frame( auc.exp = roc.dd$roc.allS, otu = names(roc.dd$otu) )
#			str(auc.te)
		} else if (j=='test') {
			auc.te = data.frame( auc.test = roc.dd$roc.allS, otu = names(roc.dd$otu) ) }
		
		auc.all = left_join(auc.all, auc.te, by = 'otu', suffix = c('', glue('.{i}')), copy=T)
	} 
}
	
auc.all = select(auc.all, -'auc.test',-'auc.exp')
str(auc.all)
	
# ... join with taxonomy table
abc = data.frame(seq.order = letters[1:length(unique(taxadd$sum.order))], 
		order = sort(unique(taxadd$sum.order), decreasing = T))
auc.all = left_join(auc.all, select(taxadd, 'otu','order','class','family','sum','sum.order'),
		  by = 'otu') 
	
auc.all$oOrder = sapply(1:nrow(auc.all), function(x) paste(abc$seq.order[abc$order==auc.all$
		  sum.order[x]], '.', auc.all$order[x], '.', auc.all$sum.order[x], sep = '')) 
str(auc.all); rm(abc)
	

aa = 0 
for (j in 1:2 ) { 
	aa = aa + 1
	set = setS[j]
	
	ii = aa+1; jj = ii + 2
	
	ab = strsplit(names(auc.all)[ii], '[.]')[[1]][2]
	auc.1 = select(auc.all, all_of(ii), all_of(jj), 'sum', 'order', 'class', 'family', 'oOrder')
	
	auc.1 = auc.1 %>% rename(auc.test = 2, auc.train = 1, incidence = 3) 
	dd.ggplot[[aa]] = auc.1
	
	gp = ggplot(dd.ggplot[[aa]], aes(auc.train, auc.test)) 
	gp = gp + geom_point(aes(colour = factor(oOrder), size = incidence)) + 
	     scale_size(range = c(1, 7)) + 
	     scale_colour_manual(values = colorRampPalette(c('dodgerblue3','firebrick2','yellow'))(length(unique(auc.1$oOrder))),
	         labels = as.character(sapply(unique(auc.all$oOrder), function(cc) str_split(cc, '[.]')[[1]][2]))) + 
	     geom_smooth(method = 'lm', se = F, colour='gray' ) + 
	     geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = 'gray', size = 1.5) + 
	     theme(panel.background = element_rect(fill='snow')) + 
	     ggtitle(glue('{set}')) + 
	     geom_hline(yintercept = .5, linetype = "solid", colour = 'gray', size = .35) +
	     geom_vline(xintercept = .5, linetype = "solid", colour = 'gray', size = .35) + 
	     geom_hline(yintercept = mean(auc.1$auc.test, na.rm=T), linetype = "dashed", colour = 'red' ) +
	     geom_vline(xintercept = mean(auc.1$auc.train), linetype = "dashed", colour = 'red' ) +
	     xlim(0,1) + ylim(0,1) + 
	     labs(colour = 'Order', size = 'Incidence', x = 'AUC (train)', y = 'AUC (test)') 
	
	gp = gp + annotate(geom = "text", y = 0, x = (mean(auc.1$auc.train)+.05 ), label = glue('{round(mean(auc.1$auc.train),2)}')) + 
		 annotate(geom = "text", x = .05, y = mean(auc.1$auc.test, na.rm=T)*1.05, 
				  label = glue('{round(mean(auc.1$auc.test, na.rm=T),2)}') ) 
	if (aa == 1) {
		plot.list[[aa]] = gp + theme(legend.position='right')} else {plot.list[[aa]] = gp + theme(legend.position='none') } 
}
	
# .... plot for the supplement 
# pdf(here(predpath, 'plot', glue('linearDNN-test-train-{mm}_{varsName}_{abund}_{cv}_{period}_min{minocc}_tuned_{date.model.run}.pdf')), width = 11, height = 5.5)
	
grid.arrange(plot.list[[1]], plot.list[[2]], nrow = 1, widths = c(.55, .45))
	
dev.off()
	

```



```{r prediction}
## produce & save predicted occurrence probability of each OTU based on the best models of all metrics
for (pred in c('explain', 'test')) {
	# explanatory AUC
	if (pred == 'explain') {
		newdd = NULL
		newsp = NULL
		otudd = s.otu.train
		set = 'explain'
	}
	# predictive AUC
	if (pred == 'test') {
		newdd = env.test.scale
		newsp = XY.test.scale
		otudd = s.otu.test
		set = 'test'
	}
	
	for (i in unique(match(pp, unique(pp)))) {
		otudd1 = otudd
		mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
		model1 = readRDS(here(modpath, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')) )
		
		pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata = newdd, SP = newsp, type = 'link')) , along = -1L), 2:3, mean) %>% unname
	# To be noticed, sjSDM 1.0.1 changed the model result structure.
	# If sjSDM 1.0.1 is installed, predict function cannot be run with the model result file provided here (model1) 
	#	attr(pred.dd, 'dimnames') 
		
		otudd1 = data.frame(otudd1)
		names(otudd1) = names(otu.pa.train)
		
		otudd1 = rbind(otudd1, count = (base::colSums(otudd1)>0 )*1 )
		
		# some OTUs don't occur in test set
		pred.dd = pred.dd[ , which(otudd1[nrow(otudd1),] == 1)]
		
		otudd1 = otudd1[1:nrow(pred.dd), which(otudd1[nrow(otudd1), ] == 1)]
#		dim(otudd1)
		
		otudd.pa = (otudd1>0)*1
#		table(otudd.pa==otudd1)
		j=30
		plot(pred.dd[,j], ylim = 0:1)
		points(otudd.pa[,j], col = 'red')
		
		# .... calculate AUC
		roc.dd = sapply(1:ncol(otudd1), function(j) as.numeric(pROC::roc( response = otudd.pa[,j], predictor = pred.dd[,j], direction = '<', quiet = T)$auc))
		
		auc.mean = mean(roc.dd)
		
#		saveRDS(list(pred.Y = pred.dd, otu = otudd1, roc.allS = roc.dd, auc.mean = auc.mean), here(predpath, 'rdata', glue('roc_result_{set}_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')))
		
		rm(pred.dd, model1, roc.dd, auc.mean, mm, otudd1, otudd.pa)
	}
}
	

```


```{r make-table}
### make table for plotting

## make taxonomy table
taxadd = data.frame(sum = colSums(otu.pa.train>0), 
					otu = names(otu.pa.train)) %>%
					arrange(desc(sum)) %>%
					mutate(sum.seq = 1:ncol(otu.pa.train))
str(taxadd)
	
# add taxonomic information 
taxadd$class = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[1])
taxadd$order = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[2])
taxadd$family = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[3])
	
taxadd = left_join(taxadd, (taxadd %>% count(order)), by = 'order') %>% 
		 rename(sum.order = n)
str(taxadd)
	

## add prediction(auc) data
auc.all = data.frame(otu = as.character(names(otu.pa.train)), 
					 auc.test = 0.1, 
					 auc.exp = 0.1 )
str(auc.all)
	
formula1 = varsName
	
# names(maxdd)
	
# load saved prediction data
for (j in c('explain', 'test')) {
	set = j
	for ( i in unique(match(pp, unique(pp))) ) { 
		mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
		
		roc.dd = readRDS(here(predpath,'rdata', glue('roc_result_{set}_{period}_{abund}_{cv}_min{minocc}_{formula1}_{mm}_{date.model.run}.RDS'))) 
		
		formula = paste0(abund, '.', formula1, '-', mm) 
		# ... make long table
		if (j=='explain') {
			auc.te = data.frame( auc.exp = roc.dd$roc.allS, otu = names(roc.dd$otu) )
#			str(auc.te)
		} else if (j=='test') {
			auc.te = data.frame( auc.test = roc.dd$roc.allS, otu = names(roc.dd$otu) ) }
		
		auc.all = left_join(auc.all, auc.te, by = 'otu', suffix = c('', glue('.{formula}')), copy=T)
	} 
}
	
# .. after all metrics are added, sort the table
auc.all = select(auc.all, -'auc.test',-'auc.exp')
str(auc.all)
	
# ... join with taxonomy table
abc = data.frame(seq.order = letters[1:length(unique(taxadd$sum.order))], 
		order = sort(unique(taxadd$sum.order), decreasing = T))
auc.all = left_join(auc.all, select(taxadd, 'otu','order','class','family','sum','sum.order'),
		  by = 'otu') 
	
auc.all$oOrder = sapply(1:nrow(auc.all), function(x) paste(abc$seq.order[abc$order==auc.all$
		  sum.order[x]], '.', auc.all$order[x], '.', auc.all$sum.order[x], sep = '')) 
str(auc.all); rm(abc)
	
# table(auc.all$'auc.test.pa.vars11-AUC'>0.7)
auc.all$'auc.exp.pa.vars11-AUC' %>% mean
	
```

```{r save-metrics_result}
names(auc.all)
	
all = auc.all %>% select('otu', starts_with('auc.'))
tt = sapply(names(all)[-1], function(x) str_split(x, '[[.],-]')[[1]][c(2,5)])
tt = sapply(1:ncol(tt), function(x) tolower(paste0(tt[1,x], '.', tt[2,x])))
	
names(all) = c('otu', tt)
	
sp.res.test = all %>% select('otu', starts_with('test')) %>% rename_all(. %>% str_replace_all(.,
			  'test.', '')) %>% data.frame()
str(sp.res.test)
	
sp.res.train = all %>% select('otu', starts_with('exp')) %>% rename_all(. %>% str_replace_all(.,
			   'exp.','')) %>% data.frame()
str(sp.res.train)
	
# save(sp.res.test, sp.res.train, file = here(predpath, 'rdata', 'sp_test_results.rdata') )
	
rm(all, tt, sp.res.train, sp.res.test)
	

```


```{r boxplot-order}
## plot boxplot by order 

metric = 'AUC'		# 'AUC' , 'plr'
	
names(auc.all)
	
dd = auc.all %>% select('order', 'sum.order', 'oOrder', 'otu', 'sum', 
		ends_with(glue('-{metric}'))) %>% rename(incidence = sum)
str(dd)
	
cc = taxadd %>% count(order) %>% arrange(desc(n))
	
# set  
if (set == 'test') {
	i = which(names(dd)==glue('auc.test.{abund}.{varsName}-{metric}'))
} else if (set=='explain') {
	i = which(names(dd)==glue('auc.exp.{abund}.{varsName}-{metric}')) }
	
# pdf(here(predpath, 'plot', glue('box_order-{set}_{period}_{abund}_min{minocc}_{varsName}_{metric}_{date.model.run}.pdf')), width = 12, height = 6)
	
ggplot(dd, aes(x=order, y=dd[,i])) + geom_hline(yintercept=.75, col='gray') + 
	geom_jitter(width = 0.35, aes(colour = order, size=incidence)) + geom_boxplot(width=0.1) + 
	scale_x_discrete(limit=cc$order) + annotate(geom="text", y=-.43, x=1:length(cc$n), label=as.character(cc$n), col='red') + 
	theme_minimal() + theme(legend.position = c(.95,0.3) ) + 
	guides(color = FALSE) + annotate(geom="text", x=2, y=-.3, label=glue('AUC.{set}: {round(mean(dd[,i],na.rm=T),3)}'))
	
dev.off()
	

```


```{r correlate-auc}
## plot auc.test ~ auc.train, color ~ order, circle size ~ incidence 
# loop over all metrics
setS = sapply(1:ncol(maxdd), function(i) strsplit(names(maxdd)[i],'.valid')[[1]][1])
dd.ggplot = vector(mode = "list", length = ncol(maxdd))
plot.list = list()
	
aa = 0 
for (j in unique(match(pp, unique(pp))) ) { 
	aa = aa + 1
	set = setS[j]
	
	if (set == 'll') { set = 'log-likelihood' 
	} else if (set == 'nagel') { set = paste0("Nagelkerke","'","s",' R2')
	} else if (set == 'plr')   { set = 'positive likelihood rate' 
	} else if (set == 'cor')   { set = 'correlation'}
	
	ii = aa+1; jj = ii + n_distinct(pp)
#	print(c(names(auc.all)[ii], names(auc.all)[jj]))
	
	ab = strsplit(names(auc.all)[ii], '[.]')[[1]][2]
	auc.1 = select(auc.all, all_of(ii), all_of(jj), 'sum', 'order', 'class', 'family', 'oOrder')
	
	auc.1 = auc.1 %>% rename(auc.test = 2, auc.train = 1, incidence = 3) 
	dd.ggplot[[aa]] = auc.1
#	summary(lm(auc.test~incidence, auc.1))
	
	gp = ggplot(dd.ggplot[[aa]], aes(auc.train, auc.test)) 
	gp = gp + geom_point(aes(colour=factor(oOrder), size=incidence)) + 
	     scale_size(range = c(1, 7)) + 
	     scale_colour_manual(values = colorRampPalette(c('dodgerblue3','firebrick2','yellow'))(length(unique(auc.1$oOrder))),
	         labels = as.character(sapply(unique(auc.all$oOrder), function(cc) str_split(cc, '[.]')[[1]][2]))) + 
	     geom_smooth(method = 'lm', se = F, colour='gray' ) + 
	     geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = 'gray', size = 1.5) + 
	     theme(panel.background = element_rect(fill='snow')) + 
	     ggtitle(glue('{set}')) + 
	     geom_hline(yintercept = .5, linetype = "solid", colour = 'gray', size = .35) +
	     geom_vline(xintercept = .5, linetype = "solid", colour = 'gray', size = .35) + 
	     geom_hline(yintercept = mean(auc.1$auc.test, na.rm=T), linetype = "dashed", colour = 'red' ) +
	     geom_vline(xintercept = mean(auc.1$auc.train), linetype = "dashed", colour = 'red' ) +
	     xlim(0,1) + ylim(0,1) + 
	     labs(colour = 'Order', size = 'Incidence', x = 'AUC (train)', y = 'AUC (test)') 
	
	gp = gp + annotate(geom = "text", y = 0, x = (mean(auc.1$auc.train)+.05 ), label = glue('{round(mean(auc.1$auc.train),2)}')) + 
		 annotate(geom = "text", x = .05, y = mean(auc.1$auc.test, na.rm=T)*1.05, 
				  label = glue('{round(mean(auc.1$auc.test, na.rm=T),2)}') ) 
	if (aa==1 |aa==3 | aa==5 | aa==7) {
		plot.list[[aa]] = gp + theme(legend.position='right')} else {plot.list[[aa]] = gp + theme(legend.position='none') } 
}
	
if (n_distinct(pp) == 5) {hh = 15}; if (n_distinct(pp) == 4) {hh = 13}
	
# .... plot for only metric AUC
# nn = 'auc'; pdf(here(predpath, 'plot', glue('{nn}-test-train-AUC_{varsName}_{abund}_{cv}_{period}_min{minocc}_tuned_{date.model.run}.pdf')), width = 7.2, height = hh/2)
	
plot.list[[1]]
	
# pdf(here(predpath, 'plot', glue('test-train_{varsName}_{abund}_{cv}_{period}_min{minocc}_tuned_{date.model.run}.pdf')), width = 14, height = hh)
# this is supplement figure 'S-PERFORMANCE'
	
grid.arrange( plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]], nrow = 3, widths = c(.545,.455))
	
dev.off()
	

```


```{r AUC-prevalence}
set = 'AUC'
	
ii = 2; jj = ii + n_distinct(pp)
print(c(names(auc.all)[ii], names(auc.all)[jj]))
	
ab = strsplit(names(auc.all)[ii], '[.]')[[1]][2]
auc.1 = select(auc.all, all_of(ii), all_of(jj), 'sum', 'order', 'class', 'family', 'oOrder')
	
auc.1 = auc.1 %>% rename(auc.test = 2, auc.train = 1, incidence = 3) 
	
# linear model
m1 = lm(log(auc.1$incidence) ~ auc.1$auc.train)
m2 = lm(log(auc.1$incidence) ~ auc.1$auc.test)
summary(m2)$coefficients[2,4]
summary(m2)$r.squared
	
g1 = ggplot(auc.1, aes(auc.train, log(incidence))) + geom_point(aes(colour=factor(order))) +
		 labs(colour = 'Order', x ='AUC (train)', y = 'ln(incidence)') + geom_smooth(method='lm', se = F, color = 'gray') +
		 annotate(geom = "text", y = Inf, x = -Inf, label = bquote(paste("p-value: ",.(round(summary(m1)$coefficients[2,4],2)),', ', R^2, ': ', .(round(summary(m1)$r.squared,2)))), size = 3, vjust=1.5, hjust=-.1)
			
g2 = ggplot(auc.1, aes(auc.test, log(incidence))) + geom_point(aes(colour=factor(order))) + 
	     labs(x = 'AUC (test)', y = 'ln(incidence)')+ theme(legend.position = "none") +
	     geom_smooth(method='lm', se = F, color = 'gray') +
		 annotate(geom = "text", y = Inf, x = -Inf, label = bquote(paste("p-value: ",.(round(summary(m2)$coefficients[2,4],2)),', ', R^2, ': ', .(format(summary(m2)$r.squared, scientific = T, digits = 1)))), size = 3, vjust=1.5, hjust=-.1)
	
# pdf(here(predpath, 'plot', glue('auc-incidence_{abund}_{period}_min{minocc}_{varsName}_tuned-{metric}_{date.model.run}.pdf')), width=12, height=5.5)
	
grid.arrange(g1, g2, nrow = 1, widths = c(0.56, 0.44))
	 
dev.off()
	 

```







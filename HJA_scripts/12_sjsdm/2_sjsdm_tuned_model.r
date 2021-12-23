# last modified: Dec 23, 2021
# analyze final (tuned) model (tuning in ADA cluster)



```{r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm')
	
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','pROC', 'gridExtra','ggeffects','cowplot','plotmath','graphics', 'MuMIn') 
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
# [1] ‘0.1.8 2020.12.23
	

```


```{r set-names}
trap <- "M1"; period = "S1"
date.model.run = '20210722'
abund = 'pa'
varsName = 'vars11'
minocc = 6
cv = '5CV'; nstep=1000
	
# ....... folder structure .......
predpath = here('..','..', 'Output', 'sjsdm_prediction_outputs', glue('{varsName}_{date.model.run}'))
modpath = here('..','..', 'Output', "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
	
sjsdmV = '0.1.8'		# check!!!
	

```


```{r load-data}
abund; varsName; minocc; cv
	
# data for the model 
load(here('source', glue('forbestm_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
if (abund == 'pa')
{
	s.otu.train = as.matrix(otu.pa.train)
	s.otu.test = as.matrix(otu.pa.test)
	attr(s.otu.train, 'dimnames') = NULL
	attr(s.otu.test, 'dimnames') = NULL
}
	
str(s.otu.train)
	
names(env.test.scale)
	
# tuned results 
tuning.dd = read.table(here(modpath, 'tuning', glue('manual_tuning_sjsdm_{sjsdmV}_{varsName}_{cv}_{period}_meanEVAL_{abund}_min{minocc}_nSteps{nstep}.csv')), header=T, sep=',')
	

```


```{r analyze-tune}
### analyze tuing result & run best model with CPU

str(tuning.dd)
	
# some metrics may select the same hyper-param combination
which.max(tuning.dd$AUC.valid_mean); which.max(tuning.dd$cor.valid_mean)
which.max(tuning.dd$nagel.valid_mean); which.max(tuning.dd$ll.valid_mean)
# AUC & cor result in the same combination !
	
maxdd = select(tuning.dd, ends_with('.valid_mean'))
pp = sapply(1:ncol(maxdd), function(i) which.max(maxdd[,i])) 
# index of selected rows (best params based on different metrics)
	
maxdd[pp,]
	

## plot the tuning result 
# pdf(here(modpath, 'plot', glue('plot_tuning_sjsdm_{sjsdmV}_{period}_{cv}_{abund}_min{minocc}_{varsName}_{date.model.run}.pdf')), width=6, height=4.5)
	
par(cex=.7, las=1)
# plot all the tuning grid
xcol = 'll.train_mean'
	
plot(y = tuning.dd$AUC.valid_mean, x = tuning.dd[,xcol], 
     pch=3, col='blue', ylim = c(0.5,1), xlab='', ylab='')
points(y = tuning.dd$AUC.train_mean, x = tuning.dd[,xcol], pch=8)
	
title(ylab = 'AUC', line = 2.6, xlab = 'log-likelihood (train)')
	
# plot best combination based on AUC and other metrics
points(y = tuning.dd[pp[1], 'AUC.valid_mean'], x = tuning.dd[pp[1],xcol], pch=20, col='red')
points(y = tuning.dd[pp[2:4], 'AUC.valid_mean'], x = tuning.dd[pp[2:4],xcol], pch=20, col='gray')
	
points(y = tuning.dd[pp[1], 'AUC.train_mean'], x = tuning.dd[pp[1], xcol], pch=8, col='red')
points(y = tuning.dd[pp[2:4], 'AUC.train_mean'], x = tuning.dd[pp[2:4], xcol], pch=8, col='gray')
	
# add vertical lines for best combinations 
abline(v = tuning.dd[pp[1],xcol], lty=2, col='red')
abline(v = tuning.dd[pp[2:4],xcol], lty=2, col='gray')
	
# add text for the best combinations
# need to change after re-run!!!!!
text(x = tuning.dd[pp[2:3],xcol]*.975, y=rep(.5, 2), c('log-likelihood', expression("Nagelkerke's R"^2)), cex=.75)
text(x = tuning.dd[pp[4],xcol]*.953, y=.497, 'positive likelihood ratio', cex=.75)
text(x = tuning.dd[pp[1],xcol]*.988, y=.497, 'AUC', cex=.75)
text(x = tuning.dd[pp[5],xcol]*.96, y=.4968, ', correlation', cex=.75)
	
legend('topright',pch=c(8,3),col=c('black','blue'),c('auc.train','auc.valid'),bty='n')
	
dev.off()
	

## model with optimal parameters
# load selected hyper-params
lambda.env = tuning.dd[pp, 'lambda.env']
alpha.env = tuning.dd[pp, 'alpha.env']
lambda.sp = tuning.dd[pp, 'lambda.sp']
alpha.sp =  tuning.dd[pp, 'alpha.sp']
lambda.bio = tuning.dd[pp, 'lambda.bio']
alpha.bio =  tuning.dd[pp, 'alpha.bio']
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
hidden1 = tuning.dd[pp, 'hidden.ind']
acti = as.character(tuning.dd[pp, 'acti.sp'])
drop = tuning.dd[pp, 'drop']
lrn = tuning.dd[pp,'lr']
	
if(abund == "pa") {
	famn = stats::binomial('probit')
	samn = 5000L; devn = 'cpu'
	itern = 170L } else stop("check abund")
	
# double-check the attributes are removed from model input data
str(env.test.scale)
str(s.otu.train)
	
# run best models with cpu based on all metrics
for (i in unique(match(pp,unique(pp))) ) {
	model.train = sjSDM(Y = s.otu.train,
	  env = DNN(data = env.train.scale, formula = ~.,
	  hidden = hidden[[hidden1[i]]], lambda = lambda.env[i], alpha = alpha.env[i], activation=acti[i], dropout=drop[i], bias=T),
	  
	  biotic = bioticStruct(lambda = lambda.bio[i], alpha = alpha.bio[i], on_diag=F, inverse = FALSE),
	  
	  spatial = linear(data = XY.train.scale, ~0+UTM_E*UTM_N, lambda = lambda.sp[i], alpha = alpha.sp[i]),
	  
	  device = devn, learning_rate = lrn[i], 
	  step_size = NULL, iter = itern, family = famn, sampling = samn 
) 
	mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
	saveRDS(model.train, here(modpath, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')) )
	rm(mm)
}
	
for (i in unique(match(pp,unique(pp))) ) {
	print(i)
	mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
	model.train = readRDS(here(modpath, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')) )
	
	pdf(here(modpath, 'plot', glue('model-history_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.pdf')), width=5, height=5)
	
	par(cex=.8)
	plot(model.train$history)
	mtext(paste0(names(maxdd)[i],': ', round(select(tuning.dd[pp,],names(maxdd)[i])[i,],3)), side=3,line=-1, adj=.95)
	mtext(paste0(names(tuning.dd)[c(2:7,9:10)],': ',select(tuning.dd[pp,],names(tuning.dd)[c(2:7,9:10)])[i,]), side=3,line=seq(-2,-9), adj=.95)
	
	model.train$history
	
	dev.off()
	rm(mm, model.train)
}
	
# write.table(tuning.dd[pp,], here(modpath, 'tuning', glue('best_manual_tuning_sjsdm_{sjsdmV}_{varsName}_{cv}_{period}_{abund}_min{minocc}_nSteps{nstep}.csv')), row.names=F, sep=',')
	

```


```{r prediction}
## produce & save predicted occurrence probability of each OTU based on the best models of all metrics
for (pred in c('explain', 'test')) {
	# explanatory AUC
	if (pred=='explain') {
		newdd = NULL
		newsp = NULL
		otudd = s.otu.train
		set = 'explain'
	}
	# predictive AUC
	if (pred=='test') {
		newdd = env.test.scale
		newsp = XY.test.scale
		otudd = s.otu.test
		set = 'test'
	}
	
	for (i in unique(match(pp,unique(pp)))) {
		otudd1 = otudd
		mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
		model1 = readRDS(here(modpath, glue('s-jSDM_tuned.model_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')) )
		
		pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=newsp, type='link')) , along = -1L), 2:3, mean)
		attr(pred.dd, 'dimnames') = NULL
		
		otudd1 = data.frame(otudd1)
		names(otudd1) = names(otu.pa.train)
		
		otudd1 = rbind(otudd1, count = (base::colSums(otudd1)>0 )*1 )
		# some OTUs don't occur in test set
		
		pred.dd = pred.dd[ , which(otudd1[nrow(otudd1),] == 1)]
		
		otudd1 = otudd1[1:nrow(pred.dd), which(otudd1[nrow(otudd1), ] == 1)]
#		dim(otudd1)
		
		otudd.pa = (otudd1>0)*1
#		table(otudd.pa==otudd1)
		
		# .... calculate AUC
		roc.dd = sapply(1:ncol(otudd1), function(j) as.numeric(pROC::roc( response = otudd.pa[,j], predictor = pred.dd[,j],direction='<',quiet=T)$auc))
		
		auc.mean = mean(roc.dd)
		
		saveRDS(list(pred.Y=pred.dd, otu=otudd1, roc.allS=roc.dd, auc.mean=auc.mean), here(predpath, 'rdata', glue('roc_result_{set}_{period}_{abund}_{cv}_min{minocc}_{varsName}_{mm}_{date.model.run}.RDS')))
		
		rm(pred.dd, model1, roc.dd, auc.mean, mm, otudd1, otudd.pa)
	}
}
	

```


```{r make-table}
### make table for plotting

## make taxonomy table
taxadd = data.frame(sum = colSums(otu.pa.train>0), otu = names(otu.pa.train))
taxadd = taxadd[order(taxadd$sum, decreasing=T),]
taxadd$sum.seq = 1:nrow(taxadd)
str(taxadd)
	
# add taxonomic information 
taxadd$order = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[2])
taxadd$class = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[1])
taxadd$family = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[3])
	
taxadd = left_join(taxadd, (taxadd %>% count(order)), by=c('order'='order')) %>% rename(sum.order=n)
str(taxadd)
	

## add prediction(auc) data
auc.all = data.frame(otu = as.character(names(otu.pa.train)), 
					 auc.test = 0.1, 
					 auc.exp = 0.1 )
str(auc.all)
	
formula1 = varsName; formula1
	
#names(maxdd)
	
# load saved prediction data
for (j in c('explain', 'test')) {
	set = j
	for ( i in unique(match(pp,unique(pp))) ) { 
		mm = substring(names(maxdd)[i], 1,regexpr('.v',names(maxdd)[i])-1)
		
		roc.dd = readRDS(here(predpath,'rdata', glue('roc_result_{set}_{period}_{abund}_{cv}_min{minocc}_{formula1}_{mm}_{date.model.run}.RDS'))) 
		
		formula = paste0(abund, '.', formula1, '-', mm) 
		# ... make long table
		if (j=='explain') {
			auc.te = data.frame( auc.exp = roc.dd$roc.allS, otu = names(roc.dd$otu) )
#			str(auc.te)
		} else if (j=='test') {
			auc.te = data.frame( auc.test = roc.dd$roc.allS, otu = names(roc.dd$otu) ) }
		
		auc.all = left_join(auc.all, auc.te, by=c('otu'='otu'), suffix=c('', glue('.{formula}')), copy=T)
	} 
}
	
# .. after all metrics are added, sort the table
auc.all = dplyr::select(auc.all, -'auc.test',-'auc.exp')
str(auc.all)
	
# ... join with taxonomy table
auc.all = left_join(auc.all, select(taxadd, 'otu','order','class','family','sum','sum.order'), by=c('otu'='otu'))
abc = data.frame(seq.order=letters[1:length(unique(taxadd$sum.order))], order=sort(unique(taxadd$sum.order),decreasing=T))
auc.all$oOrder = sapply(1:dim(auc.all)[1], function(x) paste(abc$seq.order[abc$order==auc.all$n[x]],'.',auc.all$order[x],'.',auc.all$n[x], sep=''))
names(auc.all)
	
#table(auc.all$'auc.exp.pa.vars11-AUC'>0.7)
	
```


```{r boxplot-order}
## plot boxplot by order 

metric = 'AUC'		# 'AUC' , 'plr'
dd = select(auc.all, 'order', 'sum.order', 'oOrder', 'otu', 'sum', ends_with(glue('-{metric}'))) %>% rename(incidence=sum)
names(auc.all)
str(dd)
	
cc = taxadd %>% count(order)
cc = cc[order(cc$n, decreasing=T),]
	
# set = 
if (set=='test') {i = which(names(dd)==glue('auc.test.{abund}.{varsName}-{metric}'))} else if 
(set=='explain') {i = which(names(dd)==glue('auc.exp.{abund}.{varsName}-{metric}'))}
	
# pdf(here(predpath, 'plot', glue('box_order-{set}_{period}_{abund}_min{minocc}_{varsName}_{metric}_{date.model.run}.pdf')), width=12, height=6)
	
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
for (j in unique(match(pp,unique(pp))) ) { 
	aa = aa + 1
	set = setS[j]
	
	if (set == 'll') {set = 'log-likelihood'} else if 
	(set == 'nagel') {set = paste0("Nagelkerke","'","s",' R2') } else if 
	(set == 'plr') {set = 'positive likelihood rate'} else if 
	(set == 'cor') 	{set = 'correlation'}
	
	ii = aa+1; jj = ii + n_distinct(pp)
#	print(c(names(auc.all)[ii], names(auc.all)[jj]))
	
	ab = strsplit(names(auc.all)[ii],'[.]')[[1]][2]
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
	
# plot for only metric AUC
# nn='auc'; pdf(here(predpath, 'plot', glue('{nn}-test-train-AUC_{varsName}_{abund}_{cv}_{period}_min{minocc}_tuned_{date.model.run}.pdf')), width=7.2, height=hh/2)
	
plot.list[[1]]
	
# pdf(here(predpath, 'plot', glue('test-train_{varsName}_{abund}_{cv}_{period}_min{minocc}_tuned_{date.model.run}.pdf')), width=14, height=hh)
	
if (n_distinct(pp) == 4) {
	grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]], nrow=2, widths=c(.545,.455))
}
	
dev.off()
	

```


```{r posthoc-linear}
## linear regression
# ............. lme model ................
auc.all$order =  as.factor(auc.all$order)
names(auc.all)
i = 2; metric='AUC'
abund
	
dd = auc.all[which(auc.all[,i]!='NA'),]
mod1 <- lme4::lmer(dd[,i] ~ sum + (1 | order), data = dd, REML = FALSE)
mod2 <- lme4::lmer(dd[,i] ~ 1 + (1 | order), data = dd, REML = FALSE)
anova(mod1, mod2)
	
mod1 <- lme4::lmer(dd[,i] ~ sum + (1 | order), data = dd, REML = T)
mod1 <- lme4::lmer(dd[,i] ~ 1 + (1 | order), data = dd, REML = T)
r2 = MuMIn::r.squaredGLMM(mod1)
	
# pdf(here(predpath, 'plot', glue('lme_{abund}_{period}_min{minocc}_{varsName}_tuned-{metric}_{date.model.run}.pdf')), width=7, height=7)
	
p1 = ggplot(dd, aes(x = sum, y = dd[,i], group = order, colour = order)) + geom_point() + geom_line(data = cbind(dd, pred = predict(mod1)), aes(y = pred), size = 1) + labs(x = "incidence", y = strsplit(names(dd)[i],'-')[[1]][1]) + annotate(geom='text',x=(max(dd$sum)*.9), y=.1, label=paste0('R^2: ',round(r2[[1]],3))) + scale_colour_viridis_d(option = "cividis") + theme(legend.position='bottom')
p2 = ggplot(dd, aes(x = sum, y = dd[,i], group = order, colour = order)) + geom_point() + geom_line(data = cbind(dd, pred = predict(mod1)), aes(y = pred), size = 1) + labs(x = "incidence", y = paste0('AUC.',abund,'.test.',varsName)) + theme_cowplot() + facet_wrap(vars(order)) + scale_colour_viridis_d(option = "cividis") + theme(legend.position='none')
	
grid.arrange(p1,p2, nrow=2, heights=c(.65,.35))
	
p1
	
dev.off()
	

```








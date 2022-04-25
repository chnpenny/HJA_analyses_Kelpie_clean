# tuning hyper-params in cluster with GPU


# {r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm')
	
pacman::p_load('tidyverse','sjSDM','here','conflicted','glue','pROC', 'gridExtra','ggeffects','corrplot') 
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
packageVersion('sjSDM')
# [1] ‘1.0.1 2022.04.14
	
source(here("source",'scale-train-test.r'))
	



# {r set-names}
# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"
	
period = "S1"	# ???
date.model.run = '202204'	# Apr 2022  
abund = 'pa'
varsName = 'vars11'
minocc = 6
k = 5 		# 5-folds
nstep = 1000		# sjSDM iteration
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
otupath = here('..','..','Kelpie_maps',outputidxstatstabulatefolder)
	
eodatapath = here('..','..','format_data','gis')
sppdatapath = here('..','..','format_data','otu')
	
modpath = here('..','..', 'Output', "sjsdm_general_outputs", glue('{varsName}_{date.model.run}'))
	
sjsdmV = packageVersion('sjSDM')		# package version  
	



# {r load-data}
# data for tuning from rdata 
load(here(sppdatapath, glue('fortuning_data_{period}_random_min{minocc}_{date.model.run}_{varsName}.rdata')))
	
# ... OTU data need to be matrix form for sjSDM
if (abund=='pa') {
	m.otu.train = as.matrix(otu.pa.train) %>% unname
	m.otu.test = as.matrix(otu.pa.test) %>% unname
}
	
str(m.otu.train)
	



# {r tuning}
####  Set up tuning grid 
### Create scaled data sets for each fold and model running

## ... Create tuning grid (define parameters)
lambda.env = seq(0, 1, length.out = 6)
alpha.env = seq(0, 1, length.out = 6)
lambda.sp = seq(0, 1, length.out = 5)
alpha.sp =  seq(0, 1, length.out = 5)
hidden = list(c(50L, 50L, 10L), c(25L, 25L, 10L))
hidden.ind = seq_along(hidden)
acti.sp = 'relu'
drop = c(0.1,0.2) 
sample.bio = seq(0, 1, length.out = 11)
lr = c(0.001, 0.002)
	
## Make grid of priority tuning parameters, from these in sampling; then add lower priority parameters
tune.grid = expand.grid(lambda.env = lambda.env, alpha.env = alpha.env, 
                         lambda.sp = lambda.sp, alpha.sp = alpha.sp, 
                         lr = lr, hidden.ind = hidden.ind, drop = drop,
                         acti.sp = acti.sp, stringsAsFactors = F)
head(tune.grid)
	
# choose 'nstep' combinations of the priority parameters from all the possible combinations ('tune.grid')
tune.rows = sample(1:nrow(tune.grid), size = nstep, replace = F)
	
# add in lambda.bio and alpha.bio from random samples of 'sample.bio', and add results columns
tune.results = data.frame(tr = 1:nstep,
                           tune.grid[tune.rows,], 
                           lambda.bio = sample(sample.bio, size = nstep, replace = T),
                           alpha.bio = sample(sample.bio, size = nstep, replace = T),
                           loglike_sjSDM = NA, 
                           loss= NA, 
                           AUC.train = NA, AUC.valid = NA,
                           ll.train = NA, ll.valid = NA,
                           nagel.train = NA, nagel.valid = NA,
                           plr.train = NA, plr.valid = NA,
#                           tjur.train = NA, tjur.test = NA,
                           cor.train = NA, cor.valid = NA
#                           ,auc.lt5.train = NA, auc.lt5.test= NA
                           )
	
# Add in k, to keep all cross-validation runs. Average later.
tune.results = tune.results[rep(seq(nstep), k), ] %>% add_column(., k = rep(1:k, each = nstep)) %>%
				remove_rownames
head(tune.results)
	
rm(lambda.env, alpha.env, lambda.sp, alpha.sp, hidden.ind, drop, sample.bio, acti.sp)
	
if(abund == "pa") {
	family = stats::binomial('probit')
	sampling = 5000L; device = 'gpu'
#	sampling = 50L; device = 'cpu'
#	iter = 17L
	iter = 170L; otu.train = m.otu.train} else stop("check abund")
	
# ..... tuning .....
for (i in seq_len(k)) {
	## select training (4 folds) and validation (1 fold) & scale 
	t.env.train = env.train[fold.id != i, ]
	t.env.valid = env.train[fold.id == i, ]
	
	t.otu.train = otu.train[fold.id != i,]
	t.otu.valid = otu.train[fold.id == i,]
	
	t.XY.train = XY.train[fold.id != i, ]
	t.XY.valid = XY.train[fold.id == i, ]
	
	## ... scale with source code
	a = scale.dd(t.env.train, t.env.valid)
	t.env.train = a[[1]]; t.env.valid = a[[2]]
	
	a = scale.dd(t.XY.train, t.XY.valid)
	t.XY.train = a[[1]]; t.XY.valid = a[[2]]
	rm(a)
	
	for(j in seq_len(nstep)) {
#		cat("\nk", i, "tune run", j)
		print(c(i, j))
		# subset this round of tuning parameters to make easier to insert in model specs
		tr = subset(tune.results, k == i)[j, ,drop = T]
		# train the model
		model.train = sjSDM( Y = t.otu.train,
					  env = DNN(data = t.env.train, formula = ~.,
							  hidden = hidden[[tr$hidden.ind]],
							  lambda = tr$lambda.env,
							  alpha = tr$alpha.env,
							  activation = tr$acti.sp,
							  dropout = tr$drop,
							  bias = T),
		  
					  biotic = bioticStruct(lambda = tr$lambda.bio,
										  alpha = tr$alpha.bio, 
										  on_diag = F, inverse = FALSE),
		  
					  spatial = linear(data = t.XY.train, ~0 + UTM_E*UTM_N, 
									lambda = tr$lambda.sp, 
									alpha = tr$alpha.sp),
		  
					  learning_rate = tr$lr, # part of tuning grid
					  iter = iter, family = family, #
					  sampling = sampling, device = device
		)
		## SAve each model
#		saveRDS(model.train, file.path(modpath,'tuning', paste0("s-jSDM_tuning_model_", varsName, "_", k, "CV_", i, "_tr_", j, "_", abund, ".rds")))
		
		## ... Do testing and save results in data frame
		
		tune.results$loglike_sjSDM[tune.results$k == i][j] = logLik(model.train)
		tune.results$loss[tune.results$k == i][j] = model.train$history[length(model.train$history)]
		
		for (pred in c('train','valid')) {
			if (pred == 'valid') {
				newdd = t.env.valid; newsp = t.XY.valid; otudd = t.otu.valid }
			if (pred == 'train') {
				newdd = NULL; newsp = NULL; otudd = t.otu.train }
			
			# predict for all species = sites X columns
			pred.dd = apply(abind::abind(lapply(1:3, function(i) {
							predict(model.train, newdata = newdd, SP = newsp) }
					  ),along = -1L), 2:3, mean) %>% unname
			
			# convert observed to pa (if qp)
			otudd.pa = (otudd>0)*1
			# sum(colSums(otudd.pa)==0)
			
			# AUC for all species - if not present, then NA
			auc = sapply(1:dim(otudd)[2], function(i) {
				tryCatch({
				  as.numeric(pROC::roc(otudd.pa[,i], pred.dd[,i], direction = "<", quiet=T)$auc)},
				  error = function(err){ return(NA)}
				)
			})
			
			# Extra evaluation metrics
			# log-likelihood (ll), nagelkerkes' R (nagel) & positive likelihood rate (plr) for spp 
			rsq = data.frame(ll = rep(.1, length.out = ncol(pred.dd)), 
							 nagel = rep(.1, length.out = ncol(pred.dd)), 
							 plr = rep(.1, length.out = ncol(pred.dd)),
#							   tjur = rep(NA, length.out=ncol(pred.dd)),
							 cor = rep(NA, length.out = ncol(pred.dd)) )
			
			for (m in 1:ncol(pred.dd)) { 
				p = pred.dd[ ,m]; y = otudd.pa[ ,m]
				
				loglikP = sum( log( p*y + (1-p)*(1-y) ) )
				loglikN = sum( log( mean(p)*y + (1-mean(p))*(1-y) ) )
				
				rsq$nagel[m] = (1-exp(2/length(p)*(loglikN-loglikP))) / (1-exp(2/length(p)*loglikN))
				rsq$ll[m] = loglikP
				
				tppp = sum(p*y)				# true presence
				fppp = sum(p*(1-y))			# false presence
				fapp = sum((1-p)*y)			# false absence 
				tapp = sum((1-p)*(1-y))		# true absence
				
				rsq$plr[m] = tppp/(tppp+fapp)/fppp*(fppp+tapp) # get NaN if species missing at all sites. 
				
				rsq$cor[m] = suppressWarnings(cor(p, y)) # NA when no presences in species in test set
				
			}
			
			if (pred == 'train') {
				tune.results$AUC.train[tune.results$k == i][j] = mean(auc, na.rm = T)
				tune.results$ll.train[tune.results$k == i][j] = mean(rsq$ll, na.rm = T)
				tune.results$nagel.train[tune.results$k == i][j] = mean(rsq$nagel, na.rm = T)
				tune.results$plr.train[tune.results$k == i][j]  = mean(rsq$plr, na.rm = T)
#				tune.results$tjur.train[tune.results$k == i][j]  = mean(rsq$tjur, na.rm = T)
				tune.results$cor.train[tune.results$k == i][j]  = mean(rsq$cor, na.rm = T)
#				tune.results$auc.lt5.train[tune.results$k == i][j]  = sum(auc < 0.5, na.rm = T)
			}
			  
			if (pred == 'valid') {
				tune.results$AUC.valid[tune.results$k == i][j] = mean(auc, na.rm = T)
				tune.results$ll.valid[tune.results$k == i][j] = mean(rsq$ll, na.rm = T)
				tune.results$nagel.valid[tune.results$k == i][j] = mean(rsq$nagel, na.rm = T)
				tune.results$plr.valid[tune.results$k == i][j] = mean(rsq$plr, na.rm = T)
#				tune.results$tjur.test[tune.results$k == i][j]  = mean(rsq$tjur, na.rm = T)
				tune.results$cor.valid[tune.results$k == i][j]  = mean(rsq$cor, na.rm = T)
#				tune.results$auc.lt5.test[tune.results$k == i][j]  = sum(auc < 0.5, na.rm = T)
			}
			
		} # end of evaluation loop
		rm(model.train)
	}
} # end of model loop
	
head(tune.results)
	



# {r save-tuning}
## Write results to csv 
write.table(tune.results,
			file = file.path(modpath,'tuning', paste0("manual_tuning_sjsdm_",sjsdmV, "_", varsName, "_", k, 'CV_', period, "_", abund, "_min", minocc, "_nSteps", nstep, ".csv")), 
			row.names = F, sep = ',')
	
## Average AUC by tune runs 
names(tune.results)
	
tune.mean = tune.results %>%
		group_by(across(c(-loglike_sjSDM, -loss, -k, -contains(c("train", "valid"))))) %>%
		summarise(across(contains(c("train", "valid")), list(mean = mean))) %>%
		arrange(desc(AUC.valid_mean))
head(data.frame(tune.mean))
	

write.table(tune.mean,
			file = file.path(modpath, 'tuning', paste0("manual_tuning_sjsdm_",sjsdmV, 
			"_", varsName, "_", k, "CV_", period, "_meanEVAL_", abund, "_min", minocc,
			"_nSteps", nstep, ".csv")), row.names = F, sep = ',')
	






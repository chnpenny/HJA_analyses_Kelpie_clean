# from Line 172 of '1_sjsdm_tuning.rmd'


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
		print(c(i, j))
		## read-in each model
		readRDS(file.path(modpath,'tuning', paste0("s-jSDM_tuning_model_", varsName, "_", k, "CV_", i, "_tr_", j, "_", abund, ".rds")))		# model.train
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
			# Extra evaluation metrics
			# log-likelihood (ll), nagelkerkes' R (nagel) & positive likelihood rate (plr) & TSS (true skill statistic) for spp 
			rsq = data.frame(ll = rep(.1, length.out = ncol(pred.dd)), 
							 nagel = rep(.1, length.out = ncol(pred.dd)), 
							 plr = rep(.1, length.out = ncol(pred.dd)),
							 cor = rep(NA, length.out = ncol(pred.dd)),
								   tss = rep(NA, length.out = ncol(pred.dd)) )
			for (m in 1:ncol(pred.dd)) { 
				p = pred.dd[ ,m]; y = otudd.pa[ ,m]
				tppp = sum(p*y)				# true presence
				fppp = sum(p*(1-y))			# false presence
				fapp = sum((1-p)*y)			# false absence 
				tapp = sum((1-p)*(1-y))		# true absence
				rsq$tss[m] = (tppp+tapp)/(tppp+fppp+tapp+fapp) 
				
			
			
		
		
		
	

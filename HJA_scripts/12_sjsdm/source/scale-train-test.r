# the code for scaling according to training data and apply to test

scale.dd = function(dd1, dd2) {
# dd1 -> training, dd2 -> test
# .. env data - without factors
factV = colnames(dd1)[sapply(dd1, is.factor)]
	
scaling.dd1 = dplyr::select(dd1, -any_of(factV)) %>% scale()
	
# put factors back in 
dd1.scale = data.frame(scaling.dd1, dd1[,factV, drop = F])
	
# data frame of means and sd for scaling
scaler = data.frame(t(data.frame(d.mean = attr(scaling.dd1, "scaled:center"), d.sd = attr(scaling.dd1, "scaled:scale"))))
	
dd2.scale = as.data.frame(do.call(rbind, apply(dplyr::select(dd2, -any_of(factV)), 1, function(x){(x - scaler['d.mean',])/scaler['d.sd',]} ) )) %>%
tibble::add_column(dd2[,factV, drop = F])
	
dd.scale = list(dd1.scale, dd2.scale)
	
return(dd.scale)
}
	



### 5CV on combined final data set.

rm(list=ls())


# This script runs a 5-fold cross validation using the best hyperparameters 
#  from tuning stage, on full data set 
# It saves species results. Following reviewer comments to show variability on species AUC results.

## set preferred conda env for reticulate
reticulate::use_condaenv(condaenv = "C:/Users/ggc34/.conda/envs/r-sjsdm") # laptop

getwd()

library(sjSDM)
library(pROC)
library(dplyr)

packageVersion('sjSDM')
# [1] ‘1.0.5’ 2024.01.12

period = "S1"
date.model.run = '2024'
abund = 'pa'
varsName = 'vars11'
minocc = 6
cv = '5CV'
nstep=1000

## Function for scaling data
source("./01_HJA_scripts/08_sjsdm/source/scale-train-test.r")

# ....... folder structure .......
#predpath = file.path("./04_Output/sjsdm_prediction_outputs", paste0(varsName, "_", date.model.run))
modpath = file.path("04_Output", "sjsdm_general_outputs", paste0(varsName, "_", date.model.run))
sppdatapath = file.path("./03_format_data","otu")

sjsdmV = packageVersion('sjSDM')


# Get tuned results 
# results file path
res_fp <- file.path(modpath, paste0("manual_tuning_sjsdm_",
                           sjsdmV, "_", varsName, "_", cv, "_", period,
                           "_meanEVAL_", abund, "_min", minocc, 
                           "_nSteps", nstep, ".csv"))


tuning.dd = read.table(res_fp, header = T, sep = ',')
head(tuning.dd,1)

# best hyperparameters - for top AUC
# tr lambda.env alpha.env lambda.sp alpha.sp    lr hidden.ind drop acti.sp lambda.bio
# 374          0         0       0.5     0.75 0.002          1  0.2    relu        0.2
# alpha.bio AUC.train_mean ll.train_mean nagel.train_mean plr.train_mean cor.train_mean
#      0.1       0.860051     -27.10582        0.3745389       3.851173      0.4932036
# AUC.valid_mean ll.valid_mean nagel.valid_mean plr.valid_mean cor.valid_mean
#   0.6674846     -10.29037       -0.2104988       2.359365      0.2091083

## Get best results for AUC
str(tuning.dd)0
which.max(tuning.dd$AUC.valid_mean)

tr <- tuning.dd[which.max(tuning.dd$AUC.valid_mean), , drop = TRUE]
tr

## Get model data
data_path <- file.path(sppdatapath, paste0("fortuning_data_", period, 
                                           "_random_min", minocc, "_", 
                                           date.model.run, "_", 
                                           varsName, ".rdata"))
data_path
load(data_path) # env.train, otu.pa.train, XY.train 
# these have all data (as test split set to 0%)
rm(env.test, otu.pa.test, otu.qp.test, XY.test)

# combine data sets
head(env.train) # 121 samples
str(env.train)

otu.pa.train <- as.matrix(otu.pa.train)
str(otu.pa.train) # 190 OTUs

str(XY.train) # 121 samples

## set model parameters

hidden <- list(c(50L, 50L, 10L), c(25L, 25L, 10L))
hidden.ind <- seq_along(hidden)

family = stats::binomial('probit')
sampling = 5000L
device = 'gpu'
iter = 170L

## testing
# sampling <- 100L
# iter <- 50

## simple CV
k <- 5

set.seed(99)
fold.id <- sample(rep(1:k, length.out = nrow(env.train)))
table(fold.id)


## prepare mean CV results table
tune.results <- data.frame(k = 1:k,
                           loglike_sjSDM = NA, 
                           loss= NA, 
                           AUC.train = NA, AUC.valid = NA,
                           ll.train = NA, ll.valid = NA,
                           nagel.train = NA, nagel.valid = NA,
                           plr.train = NA, plr.valid = NA,
                           tjur.train = NA, tjur.test = NA,
                           cor.train = NA, cor.valid = NA,
                           tss.train =NA, tss.valid = NA)


## species results list
rsq_all <- vector(mode = "list", length = k)

# i = 1
for (i in seq_len(k)) {

  ## select training (4 folds) and validation (1 fold) & scale 
  t.env.train = env.train[fold.id != i, ]
  t.env.valid = env.train[fold.id == i, ]
  
  t.otu.train = otu.pa.train[fold.id != i,]
  t.otu.valid = otu.pa.train[fold.id == i,]
  
  t.XY.train = XY.train[fold.id != i, ]
  t.XY.valid = XY.train[fold.id == i, ]
  
  ## ... scale with source code
  a = scale.dd(t.env.train, t.env.valid)
  t.env.train = a[[1]]; t.env.valid = a[[2]]
  
  a = scale.dd(t.XY.train, t.XY.valid)
  t.XY.train = a[[1]]; t.XY.valid = a[[2]]
  rm(a)
  
  # run model
  model.train = sjSDM(Y = t.otu.train,
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
    
  ## ... Do testing and save results in data frame
  tune.results$loglike_sjSDM[tune.results$k == i] = logLik(model.train)
  tune.results$loss[tune.results$k == i] = model.train$history[length(model.train$history)]
    
  rsq_list <- list(train = NA, valid = NA)
  
    for (pred in c('train','valid')) {
      if (pred == 'valid') {
        newdd = t.env.valid; newsp = t.XY.valid; otudd = t.otu.valid }
      if (pred == 'train') {
        newdd = NULL; newsp = NULL; otudd = t.otu.train }
      
      # predict for all species = sites X columns
      pred.dd = unname(apply(abind::abind(lapply(1:3, function(i) {
        predict(model.train, newdata = newdd, SP = newsp) }
      ),along = -1L), 2:3, mean))
      
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
      rsq = data.frame(k = i,
                       ll = rep(.1, length.out = ncol(pred.dd)), 
                       nagel = rep(.1, length.out = ncol(pred.dd)), 
                       plr = rep(.1, length.out = ncol(pred.dd)),
                       tjur = rep(NA, length.out=ncol(pred.dd)),
                       cor = rep(NA, length.out = ncol(pred.dd)),
                       tss = rep(NA, length.out = ncol(pred.dd)))
      
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
        
        tjur <- base::diff(tapply(p, y, mean, na.rm = T))
        rsq$tjur[m] <- ifelse(length(tjur) > 0, tjur, NA)
        
        rsq$cor[m] = suppressWarnings(cor(p, y)) # NA when no presences in species in test set
        rsq$tss[m] = (tppp+tapp)/(tppp+fppp+tapp+fapp)
      }
      
      if (pred == 'train') {
        # mean of all OTUs
        tune.results$AUC.train[tune.results$k == i] <- mean(auc, na.rm = T)
        tune.results$ll.train[tune.results$k == i] <- mean(rsq$ll, na.rm = T)
        tune.results$nagel.train[tune.results$k == i] <- mean(rsq$nagel, na.rm = T)
        tune.results$plr.train[tune.results$k == i] <- mean(rsq$plr, na.rm = T)
        tune.results$tjur.train[tune.results$k == i] <- mean(rsq$tjur, na.rm = T)
        tune.results$cor.train[tune.results$k == i] <- mean(rsq$cor, na.rm = T)
        tune.results$tss.train[tune.results$k == i] <- mean(rsq$tss, na.rm = T)
        
        # to save species data
        colnames(rsq) <- paste0(colnames(rsq), "_", pred)
        
        ## add AUC to rsq
        rsq$AUC_train = auc
        
        ## add OTU names
        rsq$OTU <- colnames(otu.pa.train)
        
      }
      
      if (pred == 'valid') {
        
        # do mean
        tune.results$AUC.valid[tune.results$k == i] <- mean(auc, na.rm = T)
        tune.results$ll.valid[tune.results$k == i] <- mean(rsq$ll, na.rm = T)
        tune.results$nagel.valid[tune.results$k == i] <- mean(rsq$nagel, na.rm = T)
        tune.results$plr.valid[tune.results$k == i] <- mean(rsq$plr, na.rm = T)
        tune.results$tjur.test[tune.results$k == i] <- mean(rsq$tjur, na.rm = T)
        tune.results$cor.valid[tune.results$k == i] <- mean(rsq$cor, na.rm = T)
        tune.results$tss.valid[tune.results$k == i] <- mean(rsq$tss, na.rmt =T)
        
        # to save species data
        colnames(rsq) <- paste0(colnames(rsq), "_", pred)
        
        ## add AUC to rsq
        rsq$AUC.valid = auc
        
        rsq$OTU <- colnames(otu.pa.train)
        #head(rsq)
      }
      
      rsq_list[[pred]] <- rsq
      
    } # end of evaluation loop
  
  # str(rsq_list, 2)
  rsq.df <- merge(rsq_list$train, rsq_list$valid, by = "OTU", all = TRUE)
  head(rsq.df)
  
  rsq_all[[i]] <- rsq.df
  
  rm(model.train)
  
  
  } # end of CV loop

head(tune.results)

rsq_final <- do.call(rbind, rsq_all)
head(rsq_final)

colMeans(tune.results)


## add species results from test data AUC for comparison
resFolder = file.path("./04_Output/sjsdm_prediction_outputs", paste0(varsName, "_", date.model.run))

## add incidence
incidence <- data.frame(incidence = colSums(otu.pa.train), OTU = dimnames(otu.pa.train)[[2]])
rownames(incidence) <- NULL
head(incidence)

## Filter species by auc
auc.filt = 0.70

auc_by_spp <- rsq_final %>%
  group_by (OTU) %>%
  summarise(min = min(AUC.valid, na.rm =T),
            max = max(AUC.valid, na.rm = T),
            n = sum(!is.na(AUC.valid)),
            mean = mean(AUC.valid, na.rm = TRUE),
            median = median(AUC.valid, na.rm = TRUE),
            gt07 = mean >= auc.filt, 
            gt07md = median >= auc.filt,
            .groups = "drop") %>%
  dplyr::left_join(y = incidence, by = "OTU") %>%
  arrange(desc(mean))

sum(auc_by_spp$gt07) # 112

head(auc_by_spp)
table(auc_by_spp$n)
table(auc_by_spp$gt07)
table(auc_by_spp[,c("n", "gt07")], useNA = "always")
table(auc_by_spp[,c("n", "gt07md")], useNA = "always")
summary(auc_by_spp)

auc_by_spp %>%
  filter(gt07) %>%
  summarise(med_min = median(min, na.rm =TRUE),
            med_max = median(max, na.rm =TRUE),
            min = min(min, na.rm = TRUE),
            max = max(max, na.rm = TRUE))

# A tibble: 1 × 4
# med_min med_max   min   max
# 0.636   0.935    0.348   1

auc_by_spp$OTU <- factor(auc_by_spp$OTU, levels = auc_by_spp$OTU)
rsq_final$OTU <- factor(rsq_final$OTU, levels = auc_by_spp$OTU)

# add gt07 for all species to complete table
rsq_final <- rsq_final %>%
  dplyr::left_join(y = select(auc_by_spp, OTU, gt07), by = "OTU")
head(rsq_final)

library(ggplot2)

ggplot(rsq_final, aes(y = AUC.valid, x = OTU, fill = gt07))+
  scale_fill_brewer(type = "seq", palette = "Dark2")+
  geom_boxplot(linewidth = 0.2, outlier.size = 0.6, colour = "grey30", alpha = 0.2, show.legend = FALSE)+
  geom_hline(yintercept = 0.7, col = "darkred", linetype = 2, linewidth = 0.75)+
  theme(axis.text.x=element_blank())+
  ylab("AUC (test)")+
  xlab("OTUs")
ggsave(file.path(modpath, "plot", "Fig12S_spp_auc_boxplot.png"))
# "./04_Output/sjsdm_general_outputs/vars11_2024/2024/plot/spp_auc_x_incidence_boxplot.png"

ggplot(auc_by_spp, aes(y = mean, x = OTU))+
  scale_color_brewer(type = "seq", palette = "Dark2")+
  geom_errorbar(aes(ymin = min, ymax = max), col = "grey60")+
  geom_point(show.legend = FALSE, size = 0.7, col = "grey30")+
  geom_hline(yintercept = 0.7, col = "darkred", linetype = 2, linewidth = 0.75)+
  theme(axis.text.x=element_blank())+
  ylim(c(0,1))+
  ylab("AUC (test)")+
  xlab("OTU (ordered by decreasing mean AUC)")


ggplot(auc_by_spp, aes(y = mean, x = incidence, col = gt07))+
  scale_color_brewer(type = "seq", palette = "Dark2")+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin = min, ymax = max), col = "grey60", width = 0.75, size = 0.5)+
  geom_hline(yintercept = 0.7, col = "black", linetype = 2, linewidth = 0.8)+
  #theme(axis.text.x=element_blank())+
  ylim(c(0,1))+
  ylab("AUC (test)")+
  xlab("OTU incidence (across all sample sites)")

head(auc_by_spp)

save(rsq_final, auc_by_spp, tune.results, 
     file = file.path(modpath, "spp_test_data.rdata"))

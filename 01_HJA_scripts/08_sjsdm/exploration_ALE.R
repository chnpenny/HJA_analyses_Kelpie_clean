#### 
library(iml)
library(ggplot2)
library(sjSDM)

varsName <- 'vars11'
date.model.run <- '2023'
abund <- "pa"
k <- 5
minocc <- 6; period <- "S1"
nstep <- 1000
resFolder <- file.path("./04_Output", "sjsdm_general_outputs", paste0(varsName, "_", date.model.run))
modFull_sp <- readRDS(file.path(resFolder, "s-jSDM_tuned.model_S1_pa_5CV_min6_vars11_AUC_2023.RDS"))

Env = modFull_sp$settings$env$data
Sp = modFull_sp$settings$spatial$data
Y = modFull_sp$data$Y


ALE_plots = 
  sapply(1:5, function(i) {
    
    pf = function(model, newdata, K = i) {
      env = newdata[,1:ncol(Env)]
      sp = newdata[,-(1:ncol(Env))]
      pred = predict(model, newdata= env, SP = sp, type = "raw")
      return(pred[,K])
    }
    explainer = iml::Predictor$new(modFull_sp, data = cbind(Env, Sp), y = Y[,i],predict.function = pf)
    return(FeatureEffect$new(explainer, "be30"))
    
  })
res = 
  lapply(1:5, function(j) {
    df = ALE_plots[[j]]$results
    df$species = paste0("Species_", j)
    return(df)
  })
res = do.call(rbind, res)
ggplot(res, aes(x=be30, y = .value, color = species)) + geom_line() + theme_bw() + labs(y = "ALE") 



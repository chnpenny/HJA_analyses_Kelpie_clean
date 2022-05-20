#### Beta irreplaceability vs covariates


## Local
wd <- here::here()
wd 
setwd(wd)
getwd()

library(raster)

utm10N <- 32610

here()


gis_in <- here::here('03_format_data','gis','raw_gis_data')
gis_out <- here::here('03_format_data','gis','processed_gis_data')
dir(gis_in)

# load beta irreplaceability
beta_irr <- raster(file.path(gis_out, "r_utm", "beta_r_prob_noagg.tif"))
names(beta_irr) <- "beta"

# load raster as brick
allBrck <- brick(file.path(gis_out, "r_utm", "allStack_aoi.tif"))
#allBrck <- brick(file.path("J:/UEA/Oregon/gis/processed_gis_data", "r_utm", "allStack_aoi.tif"))
load(file.path(gis_out, "brNames.rdata")) # load brick names

# get names and name groups
names(allBrck) <- brNames
allBrck
names(allBrck)

## Load old growth 
ogsi <- raster(file.path(gis_in, "ogsi_2012_smoothed.tif"))
ogsi.utm <- projectRaster(ogsi, to = allBrck)
names(ogsi.utm) <- "ogsi"
plot(ogsi.utm)

# reduce to covariates of interest
predStack <- stack(allBrck[[c("be30", "insideHJA","cut_msk", "DistRoad","DistStream")]], ogsi.utm)
predStack <- crop(predStack, beta_irr)
predStack <- addLayer(predStack, beta_irr)
predStack <- mask(predStack, beta_irr)

rm(ogsi, ogsi.utm, brNames, beta_irr, allBrck)
plot(predStack)

# Sample for plots
set.seed(99)
beta.data <- data.frame(sampleRandom(predStack, size = 1000, na.rm = TRUE))
beta.data$insideHJA <- factor(beta.data$insideHJA, labels = c("Outside", "Inside"))
beta.data$cut_msk <- factor(beta.data$cut_msk, labels = c("Outside", "Inside"))

str(beta.data)
summary(beta.data)
table(beta.data$insideHJA)
table(beta.data$cut_msk)

### Do plots 

library(ggplot2)

alpha = 0.3

p1 <- ggplot(beta.data, aes(cut_msk,beta))+
  geom_boxplot(outlier.shape = NA, na.rm = TRUE)+
  geom_jitter(width = 0.1, col = "grey10", alpha = 0.1, na.rm = TRUE)+
  ylab("Irreplaceability")+
  xlab("Plantation areas")

p1

p2 <- ggplot(beta.data, aes(ogsi, beta))+
  geom_point(alpha = alpha)+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr"), se = FALSE)+
  ylab("Irreplaceability")+
  xlab("Old growth index")
p2

p3 <- ggplot(beta.data, aes(be30, beta))+
  geom_point(alpha = alpha)+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr"), se = FALSE)+
  ylab("Irreplaceability")+
  xlab("Elevation (m)")
p3

p4 <- ggplot(beta.data, aes(DistStream, beta))+
  geom_point(alpha = alpha)+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr"), se = FALSE)+
  ylab("Irreplaceability")+
  xlab("Distance to streams (m)")

p4

library(cowplot)

plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("a)", "b)", "c)", "d)"),
          label_x = 0, label_y= 0, 
          hjust = -0.5, vjust = -0.5)
ggsave(ggsave("04_Output/fig_irreplaceability_plots.png"))





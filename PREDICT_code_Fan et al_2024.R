# * BEGIN----------------------------------------

#Author:Shunxiang Fan
#Data: 2024.04.23

# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------

# * Preparation-----------------------------------------------

# clean the environment #
rm(list=ls())

# load required libraries
library(predictsFunctions)
library(StatisticalModels)
library(sf)
library(ggplot2)
library(raster)
library(BioStatR)
library(ggthemes)
library(dplyr)
library(ggeffects)
library(devtools)
library(sjPlot)
library(cowplot)
library(sjmisc)

# extra packages called in place,function from (Outhwaite et al 2022) 
source("D:/R/R_Study/Agriculture and climate change are reshaping insect biodiversity worldwide_R/Functions.R")

#Set the data folder
dataDir <- "D:/R/1_PREDICT20221104/"
outDir <- "D:/R/1_PREDICT20221104/"

# * Data selection from PREDICT-----------------------------------------------
# Read in the PREDICTS data
pred.data <- readRDS(paste0("D:/R/predicts/database.rds")) 

### organise using functions from predictsFunctions package ###
# correct sampling effort 
predicts <- CorrectSamplingEffort(pred.data)

# merge sites: this combines potential subsamples within one site
predicts <- MergeSites(predicts) 

# select animal data with Class
predicts1=droplevels(filter(predicts,Class=="Aves"| Class=="Insecta"| Class=="Arachnida"| Class=="Mammalia"| Class=="Amphibia"| Class=="Gastropoda"| Class=="Reptilia"| Class=="Diplopoda"| Class=="Chilopoda"| Class=="Entognatha"| Class=="Malacostraca") )

# remove those sites that are the problem
# animal should not have diversity metric "percent cover"
predicts1 <- predicts1[!predicts1$Diversity_metric == "percent cover", ]

# remove rows where land use or use intensity info is missing
predicts.complete <- droplevels(predicts1[(
  predicts3$Predominant_land_use!="Cannot decide"),])
predicts.complete <- droplevels(predicts.complete[(
  predicts.complete$Use_intensity!="Cannot decide"),])

# Calculate site level metrics
pred.sites.metrics <- SiteMetrics(predicts.complete, extra.cols = c("Predominant_land_use", "SSB", "SSBS", "Class","group","UN_region","Country","order","Family")) 

# remove sites with land use "Cannot decide" and "Urban" 
sites.sub<- pred.sites.metrics[!pred.sites.metrics$Predominant_land_use %in% c( "Cannot decide","Urban"), ]

# remove sites with NA in lat/long columns 
sites.sub <- sites.sub[!is.na(sites.sub$Longitude),  ]

# remove sites with NA in Total_abundance columns 
sites.sub <- sites.sub[!is.na(sites.sub$Total_abundance),  ]

# convert the PREDICTS lat/longs into spatial points
sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]
sites.sub_xy <- SpatialPoints(sites.sub_xy, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84"))

# * Extract Variables -----------------------------------------------

# ** Extract N deposition data-----------------------------------------------

#read variables raster

Nitrogen_Deposition_noy <- raster(paste0("D:/R/variables/drywet_noy_layer.tif"))#oxidized nitrogen deposition
Nitrogen_Deposition_nhx <- raster(paste0("D:/R/variables/drywet_nhx_layer.tif"))#reduce nitrogen deposition
Nitrogen_Deposition_t <- raster(paste0("D:/R/variables/drywet_total_layer.tif"))#total nitrogen deposition

# get the values for each site
sites.sub$NDnoy<- extract(Nitrogen_Deposition_noy , sites.sub_xy, na.rm = FALSE)
sites.sub$NDnhx<- extract(Nitrogen_Deposition_nhx , sites.sub_xy, na.rm = FALSE)
sites.sub$NDt<- extract(Nitrogen_Deposition_t , sites.sub_xy, na.rm = FALSE)


# ** Extract percentage of natural/semi-natural habitat (PNSH) and cropland -----------------------------------------------

#read variables raster
Percentage_nh <- raster(paste0("D:/R/variables/PNH_1km.tif"))

crop_p<- raster(paste0("D:/R/variables/CRP_1km_2005prj.tif"))

#set buffer
sites.buffer.55k <- buffer(sites.sub_xy, width = 5000, dissolve = FALSE)

#parallel processing
install.packages("parallel")
library(raster)
library(foreach)
library(doParallel)
library(Rmpi)
library(snow)
library(parallel) 

## SET PARAMETERS FOR PARALLEL
nCores <- detectCores()  # manually for non-cluster machines
cl <- parallel::makeCluster(nCores-4) # by default this uses the PSOCK mechanism as in the SNOW package
#cl <- snow::makeCluster(nCores-4)
registerDoParallel(cl)

# export data to cores
start.time <- Sys.time()

#extract cropland
mean_crop <- foreach(i = 1:length(sites.buffer.55k),.combine=rbind,.packages = c("raster","sf")) %dopar%{
  tmp=extract(crop_p,sites.buffer.55k[i,],mean,na.rm=FALSE,ID=TRUE)
  tmp <- t(tmp)
}

#extract percentage of natural/semi-natural habitat
mean_pnh <- foreach(i = 1:length(sites.buffer.55k),.combine=rbind,.packages = c("raster","sf")) %dopar%{
  tmp=extract(crop_p,sites.buffer.55k[i,],mean,na.rm=FALSE,ID=TRUE)
  tmp <- t(tmp)
}

Sys.time()-start.time
# 10.11128mins

#close parallel processing
stopCluster(cl)

#merge into sites.sub 
sites.sub$crp <- mean_crop
sites.sub$pnh <- mean_pnh

# ** Extract climate zone -----------------------------------------------
#read climate zone  raster
Climatez <- raster(paste0("D:/R/variables/Beck_KG_V1_present_0p083.tif"))

#extract climate zone for each sites
sites.sub$Climate <- extract(Climatez , sites.sub_xy, na.rm = FALSE)

### recode climatezone data
sites.sub$climatezone <- dplyr::recode(sites.sub$Climate,
                                              '1' = 'zoneA','2' = 'zoneA','3' = 'zoneA',
                                              '4' = 'zoneB','5' = 'zoneB','6' = 'zoneB','7' = 'zoneB',
                                              '8' = 'zoneC','9' = 'zoneC','11' = 'zoneC','12' = 'zoneC',
                                              '14' = 'zoneC','15' = 'zoneC','16' = 'zoneC',
                                              '18' = 'zoneD','19' = 'zoneD','22' = 'zoneD','25' = 'zoneD',
                                              '26' = 'zoneD','27' = 'zoneD',
                                              '29' = 'zoneE', '0' = '0')

# * Data save-----------------------------------------------
# where there are NAs, change to 0
sites.sub[is.na(sites.sub$NDnoy), "NDnoy"] <- 0
sites.sub[is.na(sites.sub$NDnhx), "NDnhx"] <- 0
sites.sub[is.na(sites.sub$NDt), "NDt"] <- 0
sites.sub[is.na(sites.sub$pnh), "pnh"] <- 0
sites.sub[is.na(sites.sub$crp), "crp"] <- 0

# remove any rows that have NA in the variable columns
# remove rows with NAs for any variable of interest
sites.sub <- sites.sub[!is.na(sites.sub$NDnoy), ] # 
sites.sub <- sites.sub[!is.na(sites.sub$NDnhx), ] # 
sites.sub <- sites.sub[!is.na(sites.sub$NDt), ] # 
sites.sub <- sites.sub[!is.na(sites.sub$pnh), ] # 
sites.sub <- sites.sub[!is.na(sites.sub$crp), ] # 
sites.sub <- sites.sub[!is.na(sites.sub$Climate), ] # 

#save dataset
save(sites.sub, file = paste0(outDir, "/PREDICTS_dataset_variables.rdata"))
load(file = paste0(outDir, "/PREDICTS_dataset_variables.rdata"))

# * Data transformations-----------------------------------------------

# subset columns
final.data.trans <- sites.sub[, c(1,7:9, 13:34)]

# standardize all continuous variables

final.data.trans$NDRS <- scale(final.data.trans$NDt)
final.data.trans$NDRSnoy <- scale(final.data.trans$NDnoy)
final.data.trans$NDRSnhx <- scale(final.data.trans$NDnhx)
final.data.trans$pnhRS <- scale(final.data.trans$pnh)
final.data.trans$crpRS <- scale(final.data.trans$crp)

# drop unused levels of factors
final.data.trans <- droplevels(final.data.trans)

# set land use as character variable
final.data.trans$Predominant_land_use <- as.character(final.data.trans$Predominant_land_use)

# combine secondary land uses
final.data.trans$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans[final.data.trans$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"

# set factor levels of predominant land use
final.data.trans$Predominant_land_use <- factor(final.data.trans$Predominant_land_use,
                                                levels=c("Primary vegetation","Secondary vegetation", "Cropland","Pasture","Plantation forest"))

# separate out the data where abundance column is not NA
final.data.trans <- final.data.trans[!is.na(final.data.trans$Total_abundance), ] # 

# log the abundance values
final.data.trans$logAbun <- log(final.data.trans$Total_abundance+1)

#save final dataset
save(final.data.trans, file = paste0(outDir, "/PREDICTS_dataset_inc_variables_TRANS_log.rdata"))


# * Data summary-----------------------------------------------
table(final.data.trans$Class, final.data.trans$Predominant_land_use)
table(final.data.trans$Class, final.data.trans$UN_region)
table(final.data.trans$Class, final.data.trans$climatezone)
#
length(unique(final.data.trans[!is.na(final.data.trans$logAbun) , 'SS'])) # 443
length(unique(final.data.trans[!is.na(final.data.trans$logAbun) , 'SSB'])) # 1472
length(unique(final.data.trans[!is.na(final.data.trans$logAbun) , 'SSBS'])) # 11974

# ** Plot N deposition and sites-----------------------------------------------
library(maps)
# extract the PREDICTS points
plot_data <- final.data.trans[, c("SS", "SSBS", "Longitude", "Latitude")]

# look at nsites per study
nsites <- as.matrix(table(plot_data$SS))
nsites <- cbind(rownames(nsites), nsites[,1])

# Get the lat/long per study
lon <- unlist(apply(nsites, MARGIN = 1, FUN = function(x){plot_data[plot_data$SS == x[1], 3][1]}))
lat <- unlist(apply(nsites, MARGIN = 1, FUN = function(x){plot_data[plot_data$SS == x[1], 4][1]}))

nsites <- cbind(nsites, lon, lat)

nsites <- as.data.frame(nsites)

colnames(nsites) <- c("SS", "nsites", "lon", "lat")

nsites$nsites <- as.numeric(as.character(nsites$nsites))
nsites$lon <- as.numeric(as.character(nsites$lon))
nsites$lat <- as.numeric(as.character(nsites$lat))

map.world <- map_data('world')

# transfrom n deposition raster into dataframe
world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))

test_spdf <- as(Nitrogen_Deposition_t, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
test_df$logvaule <-as.numeric(log(test_df$drywet_total_layer+1))
colnames(test_df) <- c("value", "x", "y")

#read world shp
sd_sf <- st_read("E:/continent.shp") 

library(colorBlindness)
min_value <- round(min(test_df$logvaule),1)
max_value <- max(test_df$logvaule)
max_value <- floor(max_value * 10) / 10
#plot
ggplot()+
  geom_raster(data=test_df,aes(x = x, y = y,fill=logvaule))+ 
  scale_fill_gradientn(breaks = c(min_value,4.5 ,max_value),guide = "colourbar",
                       colours =Blue2Orange12Steps)+  
  geom_sf(data=sd_sf,fill=NA,col='grey', size = 0.1)+
  theme_bw()+
  geom_point(data = nsites, aes(x = lon, y = lat, size = nsites), col = c("#0F0F0F"), alpha = 0.4) +
  theme(panel.grid.major = element_line(colour = "transparent"), 
        panel.background = element_blank(),
        legend.position = "right",legend.title = element_text(size = 12),
        text = element_text(size = 13)) +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  scale_size_continuous(range = c(0.2, 5), breaks = c(1, 10, 50, 100, 150)) +
  guides(fill=guide_colorbar(nrow=6,byrow=T,title = "log(value)"),size = guide_legend(nrow = 6, byrow = T,title = expression(paste(N[sites]))))


#save
ggsave(filename = paste0("D:/MAP_Predicts_points.pdf"),
       plot = last_plot(),
       width = 10,
       height = 6)

# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * (1)	Does nitrogen deposition have a significant negatively effect on animal species biodiversity?-----------------------------------------------

model_data <- final.data.trans[!is.na(final.data.trans$logAbun),] 

model_data <- model_data[!is.na(model_data$NDRS), ] 

# * Model N deposition-----------------------------------------------
#check random effect
SR_check1<- glmer(Species_richness ~ NDRS + (1|SS) ,
               data = model_data, family = "poisson")
SR_check2<- glmer(Species_richness ~ NDRS + (1|SS) +(1|SSB),
                 data = model_data, family = "poisson")
SR_check3<- glmer(Species_richness ~ NDRS + (1|SS) +(1|SSB)+(1|SSBS),
                 data = model_data, family = "poisson")

# check the AIC values - random structure from SR_check3 is lowest
AIC(SR_check1, SR_check2,SR_check3)

#check random effect
AB_check1<- glmer(logAbun ~ NDRS+ (1|SS) ,
               data = model_data, family = "gaussian")
AB_check2<- glmer(logAbun ~ NDRS+ (1|SS) +(1|SSB),
               data = model_data, family = "gaussian")
# check the AIC values - random structure from AB_check2 is lowest
AIC(AB_check1, AB_check2)


# ** richness model -----------------------------------------------
SR_ND1<- glmer(Species_richness ~ NDRS + (1|SS) +(1|SSB)+(1|SSBS),
               data = model_data, family = "poisson")
car::Anova(SR_ND1)
summary(SR_ND1)

# likelihood ratio test
SR_ND2<- glmer(Species_richness ~ 1 + (1|SS) +(1|SSB)+(1|SSBS),
               data = model_data, family = "poisson")


AIC(SR_ND1,SR_ND2)


# ** plot richness model----------------------------------------------

SR_ND <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "NDRS",
               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",optimizer="Nelder_Mead")

#set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

nd_sr <- expand.grid(
  NDRS=seq(from = min(SR_ND$data$NDRS),
           to = max(SR_ND$data$NDRS),length.out=100))

# back transform the predictors

nd_sr$ND <- BackTransformCentreredPredictor(
  transformedX = nd_sr$NDRS,originalX = final.data.trans$NDt)
# set richness and abundance to 0 - to be predicted
nd_sr$logAbun <- 0
nd_sr$Species_richness <- 0

# reference for % difference = primary vegetation and positive ND closest to 0
refRow <- which((nd_sr$ND==min(abs(nd_sr$ND))))

# # quantiles for presenting results
QSR <- quantile(x = SR_ND$data$NDRS,
                probs = exclQuantiles)

# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = SR_ND$model,data = nd_sr, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')


# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd_sr$NDRS  < QSR[1])),] <- NA
a.preds.tmean[which((nd_sr$NDRS  > QSR[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd_sr$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                            FUN = median,na.rm=TRUE))*100)-100
nd_sr$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_sr$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# plot
pnd_sr <- ggplot(data = nd_sr, aes(x = ND, y = PredMedian)) + 
  geom_line(size = 0.75,color ="#0072B2") +
  geom_ribbon(aes(ymin = nd_sr$PredLower, ymax = nd_sr$PredUpper), alpha = 0.2,fill ="#0072B2") +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#0072B2")) +
  #scale_colour_manual(values = c("#7fb80e","#009E73", "#0072B2","#E69F00","#D55E00")) +
  #facet_wrap(~Predominant_land_use, ncol = 2, labeller = as_labeller(c('Secondary vegetation' = "a Secondary vegetation", 'Cropland' = "b Cropland",'Pasture' = "c Pasture",'Plantation forest' = "d Plantation forest"))) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0, 2000)) +
  #scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
  scale_y_continuous( limits = c(-50, 25)) +
  ylab("Difference of species richness (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        legend.position = c(0.8, 0.9),
        legend.background = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Richness")

pnd_sr
ggsave(pnd_sr,filename = "E:/sr_nd_new.tiff",width = 5,height = 5, dpi = 600)


# ** abundance model----------------------------------------------

AB_ND1<- glmer(logAbun ~ NDRS+ (1|SS) +(1|SSB),
               data = model_data, family = "gaussian")
car::Anova(AB_ND1)
summary(AB_ND1)

# likelihood ratio test
AB_ND2<- glmer(logAbun ~ 1+ (1|SS) +(1|SSB),
               data = model_data, family = "gaussian")

AIC(AB_ND1,AB_ND2)


# ** plot abundance model----------------------------------------------
AB_ND <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
               fixedStruct = "NDRS",
               randomStruct = "(1|SS)+(1|SSB)",saveVars = c("SSBS"),optimizer="Nelder_Mead",maxIters=2e5)

#set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)
nd_ab <- expand.grid(
  NDRS=seq(from = min(AB_ND$data$NDRS),
           to = max(AB_ND$data$NDRS),length.out=100))

# back transform the predictors

nd_ab$ND <- BackTransformCentreredPredictor(
  transformedX = nd_ab$NDRS,originalX = final.data.trans$NDt)
# set richness and abundance to 0 - to be predicted
nd_ab$logAbun <- 0
nd_ab$Species_richness <- 0

# reference for % difference = primary vegetation and positive ND closest to 0
refRow <- which((nd_ab$ND==min(abs(nd_ab$ND))))

# # quantiles for presenting results
QAB <- quantile(x = AB_ND$data$NDRS,
                probs = exclQuantiles)

# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AB_ND$model,data = nd_ab, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-1

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd_ab$NDRS  < QAB[1])),] <- NA
a.preds.tmean[which((nd_ab$NDRS  > QAB[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd_ab$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                            FUN = median,na.rm=TRUE))*100)-100
nd_ab$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_ab$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# plot
pnd_ab <- ggplot(data = nd_ab, aes(x = ND, y = PredMedian)) + 
  geom_line(size = 0.75,color ="#009E73") +
  geom_ribbon(aes(ymin = nd_ab$PredLower, ymax = nd_ab$PredUpper), alpha = 0.2,fill ="#009E73") +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73")) +
  #scale_colour_manual(values = c("#7fb80e","#009E73", "#0072B2","#E69F00","#D55E00")) +
  #facet_wrap(~Predominant_land_use, ncol = 2, labeller = as_labeller(c('Secondary vegetation' = "a Secondary vegetation", 'Cropland' = "b Cropland",'Pasture' = "c Pasture",'Plantation forest' = "d Plantation forest"))) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0, 2000)) +
  #scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
  scale_y_continuous( limits = c(-75, 25)) +
  ylab("Difference of total abundance (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        legend.position = c(0.8, 0.9),
        legend.background = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Aundance")

pnd_ab

ggsave(pnd_ab,filename = "E:/ab_nd_new.tiff",width = 5,height = 5, dpi = 600)

# ** Test Chao estimator----------------------------------------------

# remove sites that don't have a ChaoR estimate
model_data2 <- model_data[!is.na(model_data$ChaoR), ] 
model_data2<-droplevels(model_data2)

# round the estimated species richness values to integers.
model_data2$ChaoR <- round(model_data2$ChaoR,0)

CH_ND <- GLMER(modelData = model_data2,responseVar = "ChaoR",fitFamily = "poisson",
               fixedStruct = "NDRS",
               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
car::Anova(CH_ND$model)
summary(CH_ND$model)

#### create tables of coefficients for SR and ChaoR models ###
tab_model(SR_ND$model, CH_ND$model, transform = NULL, file = paste0("E:/Table_SR_ChaoR_ND.html"), 
          show.icc = F, show.obs = F, show.ngroups = F)
# ** plot Chao estimator model----------------------------------------------
exclQuantiles <- c(0.025,0.975)

nd_ch <- expand.grid(
  NDRS=seq(from = min(CH_ND$data$NDRS),
           to = max(CH_ND$data$NDRS),length.out=100))

# back transform the predictors

nd_ch$ND <- BackTransformCentreredPredictor(
  transformedX = nd_ch$NDRS,originalX = final.data.trans$NDt)
# set richness and abundance to 0 - to be predicted
nd_ch$logAbun <- 0
nd_ch$Species_richness <- 0
nd_ch$ChaoR <- 0

# reference for % difference = primary vegetation and positive ND closest to 0
refRow <- which((nd_ch$ND==min(abs(nd_ch$ND))))

# # quantiles for presenting results
QSR <- quantile(x = CH_ND$data$NDRS,
                probs = exclQuantiles)

# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = CH_ND$model,data = nd_ch, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')


# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd_ch$NDRS  < QSR[1])),] <- NA
a.preds.tmean[which((nd_ch$NDRS  > QSR[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd_ch$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                            FUN = median,na.rm=TRUE))*100)-100
nd_ch$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_ch$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# plot
pnd_ch <- ggplot(data = nd_ch, aes(x = ND, y = PredMedian)) + 
  geom_line(size = 0.75,color ="#E69F00") +
  geom_ribbon(aes(ymin = nd_ch$PredLower, ymax = nd_ch$PredUpper), alpha = 0.2,fill ="#E69F00") +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#E69F00")) +
  #scale_colour_manual(values = c("#7fb80e","#009E73", "#0072B2","#E69F00","#D55E00")) +
  #facet_wrap(~Predominant_land_use, ncol = 2, labeller = as_labeller(c('Secondary vegetation' = "a Secondary vegetation", 'Cropland' = "b Cropland",'Pasture' = "c Pasture",'Plantation forest' = "d Plantation forest"))) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0, 2000)) +
  #scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
  scale_y_continuous( limits = c(-75, 25)) +
  ylab("Change in Chao estimated richness (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        legend.position = c(0.8, 0.9),
        legend.background = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Chao estimated")

ggsave(pnd_ch,filename = "E:/ch_nd.tiff",width = 5,height = 5, dpi = 600)

# ** Test NOy and NHx data-----------------------------------------------

# *** richness model---------------------------------------------  
# oxidized nitrogen deposition
SR_noy <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                fixedStruct = "NDRSnoy",
                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")

car::Anova(SR_noy$model)
summary(SR_noy$model)

# reduced nitrogen deposition
SR_nhx <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                fixedStruct = "NDRSnhx",
                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")

car::Anova(SR_nhx$model)
summary(SR_nhx$model)

#output model summary
tab_model(SR_noy$model, SR_nhx$model, transform = NULL, file = paste0("E:/Table_initial/Table_SR_noynhx.html"), 
          show.icc = F, show.obs = F, show.ngroups = F)

# *** abundance model---------------------------------------------  
# oxidized nitrogen deposition
AB_noy <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
                fixedStruct = "NDRSnoy",
                randomStruct = "(1|SS)+(1|SSB)",
                saveVars = c("SSBS"))
car::Anova(AB_noy$model)
summary(AB_noy$model)

# reduced nitrogen deposition
AB_nhx <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
                fixedStruct = "NDRSnhx",
                randomStruct = "(1|SS)+(1|SSB)",
                saveVars = c("SSBS"))
car::Anova(AB_nhx$model)
summary(AB_nhx$model)

# output model summary
tab_model(AB_noy$model, AB_nhx$model, transform = NULL, file = paste0("E:/Table_initial/Table_Ab_noynhx.html"), 
          show.icc = F, show.obs = F, show.ngroups = F)


# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * (2)	Whether nitrogen deposition effect on animal species biodiversity varied with  taxonomic groups?-----------------------------------------------

# * Model N deposition/Class-----------------------------------------------
table(model$Class)
#Amphibia    Arachnida  Birds    Chilopoda   Entognatha   Gastropoda      Insecta Malacostraca 
#343          840         2890           42          121          216         5392           21 
#Mammalia     Reptilia 
#1537          572 

#remove data with less records
model_data3<- model_data[!model_data$group %in% c("Chilopoda","Malacostraca"), ]
model_data3<- droplevels(model_data3)

#set levels
model_data3$Class <- factor(model_data3$Class, 
                            levels = c("Amphibia","Arachnida", "Birds","Entognatha","Gastropoda","Insecta","Mammalia","Reptilia"))

# ** richness model -----------------------------------------------
SR_ND_class1<- glmer(Species_richness ~ NDRS*Class + (1|SS) +(1|SSB)+(1|SSBS),
                     data = model_data3, family = "poisson",glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
car::Anova(SR_ND_class1)
summary(SR_ND_class1)


# likelihood ratio test
SR_ND_class2<- glmer(Species_richness ~ NDRS+group + (1|SS) +(1|SSB)+(1|SSBS),
                     data = model_data3, family = "poisson",glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
car::Anova(SR_ND_class2)
summary(SR_ND_class2)

anova(SR_ND_class1,SR_ND_class2)

# ** abundance model -----------------------------------------------

AB_ND_class1<- glmer(logAbun ~ NDRS*Class+ (1|SS) +(1|SSB),
                     data = model_data3, family = "gaussian")

car::Anova(AB_ND_class1)
summary(AB_ND_class1)

# likelihood ratio test
AB_ND_class2<- glmer(logAbun ~ NDRS+Class+ (1|SS) +(1|SSB),
                     data = model_data3, family = "gaussian")

anova(AB_ND_class1,AB_ND_class2)

AIC(AB_ND_class1,AB_ND_class2)

# ** plot abundance model----------------------------------------------
ggAB_ND_class<-ggeffect(AB_ND_class1,terms = c("NDRS[all]","group"),ci.lvl = 0.95)


plot(ggAB_ND_class)

summary(ggAB_ND_class$group)

ggAB_ND_class$ND <- BackTransformCentreredPredictor(
  transformedX = ggAB_ND_class$x,
  originalX = final.data.trans$NDt)

ggAB_ND_class$group <- factor(ggAB_ND_class$group, 
                              levels = c("Amphibia","Arachnida","Birds","Entognatha","Gastropoda","Insecta","Mammalia","Reptilia"))


#d16ba5, #c777b9, #ba83ca, #aa8fd8, #9a9ae1, #8aa7ec, #79b3f4, #69bff8, #52cffe, #41dfff, #46eefa, #5ffbf1
#("#a50026","#f46d43","#74add1","#313695","#4e72b8")
pggAB_ND_class <- ggplot(da=ggAB_ND_class, aes(x = ND, y = predicted)) + 
  geom_line(aes(col =group), size = 0.75) +
  geom_ribbon(aes(ymin = ggAB_ND_class$conf.low, ymax = ggAB_ND_class$conf.high, fill = group), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#C07A91","#3D9140","#80AFBF","#FFD104","#4e72b8","#f46d43","#313695","#DFC285","#41dfff"))) +
  scale_colour_manual(values = rev(c("#C07A91","#3D9140","#80AFBF","#FFD104","#4e72b8","#f46d43","#313695","#DFC285","#41dfff"))) +
  theme_bw() + 
  labs(fill = "Class", col = "Class") + 
  ylab("Total Abandacne") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  xlim(c(0, 4000)) +
  ylim(c(-10, 20)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13.5),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2),
        panel.spacing = unit(1,"cm"))

pggAB_ND_class

ggsave(pggAB_ND_class,filename = "E:/Fig_initial/pggAB_ND_class.tiff",width = 7,height = 5, dpi = 600)

# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * (2)	Whether nitrogen deposition effect on animal species biodiversity varied with land use types?-----------------------------------------------

# * Model N deposition/Land use-----------------------------------------------

# ** richness model--------------------------------------------- 
SR_ND_land1<- glmer(Species_richness ~ NDRS*Predominant_land_use + (1|SS) +(1|SSB)+(1|SSBS),
                    data = model_data, family = "poisson",glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
car::Anova(SR_ND_land1)
summary(SR_ND_land1)

# likelihood ratio test
SR_ND_land2<- glmer(Species_richness ~ NDRS+Predominant_land_use + (1|SS) +(1|SSB)+(1|SSBS),
                    data = model_data, family = "poisson",glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 100000)))
car::Anova(SR_ND_land2)
summary(SR_ND_land2)

anova(SR_ND_land1,SR_ND_land2)

# ** plot richness model--------------------------------------------- 
SR_ND_land <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                    fixedStruct = "NDRS *Predominant_land_use",
                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
car::Anova(SR_ND_land$model)
summary(SR_ND_land$model)

#set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

ndsr_land <- expand.grid(
  NDRS=seq(from = min(SR_ND_land$data$NDRS),
           to = max(SR_ND_land$data$NDRS),length.out=100),
  Predominant_land_use=factor(c("Primary vegetation", "Secondary vegetation", "Cropland", "Pasture","Plantation forest"),
                              levels = levels(SR_ND_land$data$Predominant_land_use)))

# back transform the predictors

ndsr_land$ND <- BackTransformCentreredPredictor(
  transformedX = ndsr_land$NDRS,originalX = final.data.trans$NDt)

# set richness and abundance to 0 - to be predicted
ndsr_land$logAbun <- 0
ndsr_land$Species_richness <- 0

# reference for % difference = primary vegetation and positive ND closest to 0
refRow <- which((ndsr_land$Predominant_land_use=="Primary vegetation") & (ndsr_land$ND==min(abs(ndsr_land$ND))))

# # quantiles for presenting results
QPV <- quantile(x = SR_ND_land$data$NDRS[
  SR_ND_land$data$Predominant_land_use=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = SR_ND_land$data$NDRS[
  SR_ND_land$data$Predominant_land_use=="Secondary vegetation"],
  probs = exclQuantiles)
QCR <- quantile(x = SR_ND_land$data$NDRS[
  SR_ND_land$data$Predominant_land_use=="Cropland"],
  probs = exclQuantiles)
QPA <- quantile(x = SR_ND_land$data$NDRS[
  SR_ND_land$data$Predominant_land_use=="Pasture"],
  probs = exclQuantiles)
QPF <- quantile(x = SR_ND_land$data$NDRS[
  SR_ND_land$data$Predominant_land_use=="Plantation forest"],
  probs = exclQuantiles)

# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = SR_ND_land$model,data = ndsr_land, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Primary vegetation") & (ndsr_land$NDRS < QPV[1])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Primary vegetation") & (ndsr_land$NDRS  > QPV[2])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Secondary vegetation") & (ndsr_land$NDRS  < QSV[1])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Secondary vegetation") & (ndsr_land$NDRS  > QSV[2])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Cropland") & (ndsr_land$NDRS  < QCR[1])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Cropland") & (ndsr_land$NDRS  > QCR[2])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Pasture") & (ndsr_land$NDRS  < QPA[1])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Pasture") & (ndsr_land$NDRS  > QPA[2])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Plantation forest") & (ndsr_land$NDRS  < QPF[1])),] <- NA
a.preds.tmean[which((ndsr_land$Predominant_land_use=="Plantation forest") & (ndsr_land$NDRS  > QPF[2])),] <- NA
# Get the median, upper and lower quants for the plot
ndsr_land$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                                FUN = median,na.rm=TRUE))*100)-100
ndsr_land$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                               FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
ndsr_land$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                               FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
ndsr_land$Predominant_land_use <- factor(ndsr_land$Predominant_land_use, levels = c("Primary vegetation", "Secondary vegetation", "Cropland", "Pasture","Plantation forest"))

# plot
pndsr_land <- ggplot(data = ndsr_land, aes(x = ND, y = PredMedian)) + 
  geom_line(aes(col = Predominant_land_use), size = 0.75) +
  geom_ribbon(aes(ymin = ndsr_land$PredLower, ymax = ndsr_land$PredUpper, fill = Predominant_land_use), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#01BB16","#305D50", "#0072B2","#E69F00","#B62F4B")) +
  scale_colour_manual(values = c("#01BB16","#305D50", "#0072B2","#E69F00","#B62F4B")) +
  theme_bw() + 
  scale_x_continuous(limits = c(100, 2500)) +
  scale_y_continuous( limits = c(-80, 20)) +
  ylab("Difference of species richness (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        legend.position = "right",
        legend.background = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Richness")

pndsr_land

ggsave(pndsr_land,filename = "E:/Fig_initial/sr_nd_land.tiff",width = 7,height = 5, dpi = 600)


# ** abundance model--------------------------------------------- 
AB_ND_land1<- glmer(logAbun ~ NDRS *Predominant_land_use+ (1|SS) +(1|SSB),
                    data = model_data, family = "gaussian")
car::Anova(AB_ND_land1)
summary(AB_ND_land1)

# likelihood ratio test
AB_ND_land2<- glmer(logAbun ~ NDRS +Predominant_land_use+ (1|SS) +(1|SSB),
                    data = model_data, family = "gaussian")

anova(AB_ND_land1,AB_ND_land2)

# ** plot abundance model--------------------------------------------- 

AB_ND_land <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
                    fixedStruct = "NDRS*Predominant_land_use",
                    randomStruct = "(1|SS)+(1|SSB)",
                    saveVars = c("SSBS"))
car::Anova(AB_ND_land$model)
summary(AB_ND_land$model)

#set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

nd_land <- expand.grid(
  NDRS=seq(from = min(AB_ND_land$data$NDRS),
           to = max(AB_ND_land$data$NDRS),length.out=100),
  Predominant_land_use=factor(c("Primary vegetation", "Secondary vegetation", "Cropland", "Pasture","Plantation forest"),
                              levels = levels(AB_ND_land$data$Predominant_land_use)))

# back transform the predictors

nd_land$ND <- BackTransformCentreredPredictor(
  transformedX = nd_land$NDRS,originalX = final.data.trans$NDt)
# set richness and abundance to 0 - to be predicted
nd_land$logAbun <- 0
nd_land$Species_richness <- 0

# reference for % difference = primary vegetation and positive ND closest to 0
refRow <- which((nd_land$Predominant_land_use=="Primary vegetation") & (nd_land$ND==min(abs(nd_land$ND))))


# # quantiles for presenting results
QPV <- quantile(x = AB_ND_land$data$NDRS[
  AB_ND_land$data$Predominant_land_use=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = AB_ND_land$data$NDRS[
  AB_ND_land$data$Predominant_land_use=="Secondary vegetation"],
  probs = exclQuantiles)
QCR <- quantile(x = AB_ND_land$data$NDRS[
  AB_ND_land$data$Predominant_land_use=="Cropland"],
  probs = exclQuantiles)
QPA <- quantile(x = AB_ND_land$data$NDRS[
  AB_ND_land$data$Predominant_land_use=="Pasture"],
  probs = exclQuantiles)
QPF <- quantile(x = AB_ND_land$data$NDRS[
  AB_ND_land$data$Predominant_land_use=="Plantation forest"],
  probs = exclQuantiles)

# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AB_ND_land$model,data = nd_land, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-1

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')


# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Primary vegetation") & (nd_land$NDRS < QPV[1])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Primary vegetation") & (nd_land$NDRS  > QPV[2])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Secondary vegetation") & (nd_land$NDRS  < QSV[1])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Secondary vegetation") & (nd_land$NDRS  > QSV[2])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Cropland") & (nd_land$NDRS  < QCR[1])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Cropland") & (nd_land$NDRS  > QCR[2])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Pasture") & (nd_land$NDRS  < QPA[1])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Pasture") & (nd_land$NDRS  > QPA[2])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Plantation forest") & (nd_land$NDRS  < QPF[1])),] <- NA
a.preds.tmean[which((nd_land$Predominant_land_use=="Plantation forest") & (nd_land$NDRS  > QPF[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd_land$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                              FUN = median,na.rm=TRUE))*100)-100
nd_land$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_land$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
nd_land$Predominant_land_use <- factor(nd_land$Predominant_land_use, levels = c("Primary vegetation", "Secondary vegetation", "Cropland", "Pasture","Plantation forest"))

# plot
pnd_land <- ggplot(data = nd_land, aes(x = ND, y = PredMedian)) + 
  geom_line(aes(col = Predominant_land_use), size = 0.75) +
  geom_ribbon(aes(ymin = nd_land$PredLower, ymax = nd_land$PredUpper, fill = Predominant_land_use), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#01BB16","#305D50", "#0072B2","#E69F00","#B62F4B")) +
  scale_colour_manual(values = c("#01BB16","#305D50", "#0072B2","#E69F00","#B62F4B")) +
  theme_bw() + 
  scale_x_continuous(limits = c(0, 2500)) +
  scale_y_continuous( limits = c(-100, 50)) +
  ylab("Difference of total abundance (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        #legend.position = c(0.8, 0.85),
        legend.position = "right",
        legend.background = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Abundance")
pnd_land

ggsave(pnd_land,filename = "E:/Fig_initial/ab_nd_land1.tiff",width = 8,height = 5, dpi = 600)


# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * (2)	Whether nitrogen deposition effect on animal species biodiversity varied with continental?-----------------------------------------------

# * Model N deposition/Continent-----------------------------------------------

# ** richness model--------------------------------------------- 
SR_ND_region1<- glmer(Species_richness ~ NDRS*UN_region + (1|SS) +(1|SSB)+(1|SSBS),
                      data = model_data, family = "poisson",glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
car::Anova(SR_ND_region1)
summary(SR_ND_region1)

# likelihood ratio test
SR_ND_region2<- glmer(Species_richness ~ NDRS+UN_region + (1|SS) +(1|SSB)+(1|SSBS),
                      data = model_data, family = "poisson")

anova(SR_ND_region1,SR_ND_region2)

# ** abundance model--------------------------------------------- 
AB_ND_region1<- glmer(logAbun ~ NDRS*UN_region + (1|SS) +(1|SSB),
                      data = model_data, family = "gaussian")

summary(AB_ND_region1)
car::Anova(AB_ND_region1)

# likelihood ratio test
AB_ND_region2<- glmer(logAbun ~ NDRS+UN_region + (1|SS) +(1|SSB),
                      data = model_data, family = "gaussian")

summary(AB_ND_region2)
car::Anova(AB_ND_region2)
anova(AB_ND_region1,AB_ND_region2)

# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * (2)	Whether nitrogen deposition effect on animal species biodiversity varied with climate zones?-----------------------------------------------

# * Model N deposition/Climate zone-----------------------------------------------
#remove data with few records
model_data4<- model_data[!model_data$climatezone1 %in% c("0","zoneE"), ]
model_data4<- droplevels(model_data4)


# ** richness model--------------------------------------------- 
SR_ND_climate1<- glmer(Species_richness ~ NDRS*climatezone1  + (1|SS) +(1|SSB)+(1|SSBS),
                       data = model_data4, family = "poisson",glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
car::Anova(SR_ND_climate1)
summary(SR_ND_climate1)

#likelihood ratio test
SR_ND_climate2<- glmer(Species_richness ~ NDRS+climatezone1  + (1|SS) +(1|SSB)+(1|SSBS),
                       data = model_data4, family = "poisson")

anova(SR_ND_climate2,SR_ND_climate1)

# ** abundance model--------------------------------------------- 
AB_ND_climate1<- glmer(logAbun ~ NDRS*climatezone1 + (1|SS) +(1|SSB),
                       data = model_data4, family = "gaussian")
car::Anova(AB_ND_climate1)
summary(AB_ND_climate1)

#likelihood ratio test
AB_ND_climate2<- glmer(logAbun ~ NDRS+climatezone2 + (1|SS) +(1|SSB),
                       data = model_data4, family = "gaussian")

anova(AB_ND_climate1,AB_ND_climate2)

# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * (3)	The effect of SNH in buffering the impact of N deposition on animal diversity recorded at different land-use types-----------------------------------------------

# * Model N deposition/percentage of natural/semi-natural habitat-----------------------------------------------

# ** richness model--------------------------------------------- 
SR_ND_PNH1<- glmer(Species_richness ~ NDRS*pnhRS + (1|SS) +(1|SSB)+(1|SSBS),
                   data = model_data, family = "poisson",glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 100000)))

car::Anova(SR_ND_PNH1)
summary(SR_ND_PNH1)

#likelihood ratio test
SR_ND_PNH2<- glmer(Species_richness ~ NDRS+pnhRS + (1|SS) +(1|SSB)+(1|SSBS),
                   data = model_data, family = "poisson",glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 100000)))

anova(SR_ND_PNH1,SR_ND_PNH2)
# ** plot richness model--------------------------------------------- 
SR_ND_PNH <- GLMER(modelData = model_data,responseVar = "Species_richness",
                   fitFamily = "poisson",
                   fixedStruct = "NDRS * pnhRS",
                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")


#set as different classes of PNH
#25% -1.137128, 50% -0.1354202, 75% 0.8654312, 100% 1.863949

ndsr_pnh<- expand.grid(
  NDRS=seq(from = min(SR_ND_PNH$data$NDRS),
           to = max(SR_ND_PNH$data$NDRS),
           length.out = 100),
  pnhRS=c(-1.137128,-0.1354202,0.8654312,1.863949))

# back transform the climate data range
ndsr_pnh$ND <- BackTransformCentreredPredictor(
  transformedX = ndsr_pnh$NDRS,
  originalX = final.data.trans$NDt)

# back transform PNH data range
ndsr_pnh$pnh <- round(BackTransformCentreredPredictor(
  transformedX = ndsr_pnh$pnhRS,originalX = final.data.trans$pnh)*100,0)

# set values for richness and Abundance
ndsr_pnh$logAbun <- 0
ndsr_pnh$Species_richness <- 0

# set the reference row
refRow <- which((ndsr_pnh$ND==min(abs(ndsr_pnh$ND))) & 
                  (ndsr_pnh$pnh==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)

QSR <- quantile(x = SR_ND_PNH$data$NDRS,
                probs = exclQuantiles)
# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = SR_ND_PNH$model,data = ndsr_pnh, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')


# set anything outside the desired quantiles to NA
a.preds.tmean[which((ndsr_pnh$NDRS  < QSR[1])),] <- NA
a.preds.tmean[which((ndsr_pnh$NDRS  > QSR[2])),] <- NA
# get the median and upper/lower intervals for plots
ndsr_pnh$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                               FUN = median,na.rm=TRUE))*100)-100
ndsr_pnh$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
ndsr_pnh$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
#ndsr_pnh$Predominant_land_use <- factor(ndsr_pnh$Predominant_land_use, levels = c("Primary vegetation", "Secondary vegetation", "Cropland", "Pasture","Plantation forest"))

ndsr_pnh$pnh <- factor(ndsr_pnh$pnh, levels = c("100", "75", "50", "25"))

# just take agriculture values
#ndsr_pnh <- ndsr_pnh[ndsr_pnh$Predominant_land_use %in% c("Cropland"), ]

# plot
pndsr_pnh <- ggplot(data = ndsr_pnh, aes(x = ND, y = PredMedian)) + 
  geom_line(aes(col =pnh), size = 0.75) +
  geom_ribbon(aes(ymin = ndsr_pnh$PredLower, ymax = ndsr_pnh$PredUpper, fill = pnh), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695","#4e72b8"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695","#4e72b8"))) +
  #facet_wrap(~Predominant_land_use, ncol = 2, labeller = as_labeller(c('Secondary vegetation' = "a Secondary vegetation", 'Cropland' = "b Cropland",'Pasture' = "c Pasture",'Plantation forest' = "d Plantation forest"))) + 
  theme_bw() + 
  labs(fill = "PNSH (%)", col = "PNSH (%)") + 
  ylab("Change in species richness (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  xlim(c(0, 2000)) +
  ylim(c(-60, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2)) 

pndsr_pnh

ggsave(pndsr_pnh,filename = "E:/Fig_initial/pndsr_pnsh.tiff",width = 7,height = 5, dpi = 600)

# ** abundance model--------------------------------------------- 
AB_ND_PNH1<- glmer(logAbun ~ NDRS*pnhRS + (1|SS) +(1|SSB),
                   data = model_data, family = "gaussian")

car::Anova(AB_ND_PNH1)

#likelihood ratio test
AB_ND_PNH2<- glmer(logAbun ~ NDRS+pnhRS + (1|SS) +(1|SSB),
                   data = model_data, family = "gaussian")
anova(AB_ND_PNH1,AB_ND_PNH2)


# ** plot abundance model--------------------------------------------- 
AB_ND_PNH <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
                   fixedStruct = "NDRS * pnhRS",
                   randomStruct = "(1|SS)+(1|SSB)",
                   saveVars = c("SSBS"))

### extract pnhRS values with pnh = 25%
data1 = final.data.trans[final.data.trans$pnh == max(final.data.trans$pnh[which(final.data.trans$pnh <= 0.25)]),]

data1$pnhRS

#set as different classes of PNH
#25% -1.137128, 50% -0.1354202, 75% 0.8654312, 100% 1.863949

#plot
ndab_pnh<- expand.grid(
  NDRS=seq(from = min(AB_ND_PNH$data$NDRS),
           to = max(AB_ND_PNH$data$NDRS),
           length.out = 100),
  pnhRS=c(-1.137128,-0.1354202,0.8654312,1.863949))


# back transform the climate data range
ndab_pnh$ND <- BackTransformCentreredPredictor(
  transformedX = ndab_pnh$NDRS,
  originalX = final.data.trans$NDt)

# back transform PNH data range
ndab_pnh$pnh <- round(BackTransformCentreredPredictor(
  transformedX = ndab_pnh$pnhRS,originalX = final.data.trans$pnh)*100,0)

# set values for richness and Abundance
ndab_pnh$logAbun <- 0
ndab_pnh$Species_richness <- 0

# set the reference row
refRow <- which((ndab_pnh$ND==min(abs(ndab_pnh$ND))) & 
                  (ndab_pnh$pnh==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)

QAB <- quantile(x = AB_ND_PNH$data$NDRS,
                probs = exclQuantiles)
# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AB_ND_PNH$model,data = ndab_pnh, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')


# set anything outside the desired quantiles to NA
a.preds.tmean[which((ndab_pnh$NDRS  < QAB[1])),] <- NA
a.preds.tmean[which((ndab_pnh$NDRS  > QAB[2])),] <- NA
# get the median and upper/lower intervals for plots
ndab_pnh$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                               FUN = median,na.rm=TRUE))*100)-100
ndab_pnh$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
ndab_pnh$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels of PNH
ndab_pnh$pnh <- factor(ndab_pnh$pnh, levels = c("100", "75", "50", "25"))


# plot
pndab_pnh <- ggplot(data = ndab_pnh, aes(x = ND, y = PredMedian)) + 
  geom_line(aes(col =pnh), size = 0.75) +
  geom_ribbon(aes(ymin = ndab_pnh$PredLower, ymax = ndab_pnh$PredUpper, fill = pnh), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695","#4e72b8"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695","#4e72b8")))+
  theme_bw() + 
  labs(fill = "PNSH (%)", col = "PNSH (%)") + 
  ylab("Change in total abundance (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')'))))+
  xlim(c(0, 2000)) +
  ylim(c(-80, 80)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2)) 

pndab_pnh

ggsave(pndab_pnh,filename = "E:/Fig_initial/pndab_pnsh.tiff",width = 7,height = 5, dpi = 600)


# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------

# * (4)	The effect of cropland coverage in mediating the impact of N deposition-----------------------------------------------

# * Model N deposition/percentage of cropland-----------------------------------------------

# ** richness model--------------------------------------------- 
SR_ND_crp1<- glmer(Species_richness ~ NDRS*crpRS + (1|SS) +(1|SSB)+(1|SSBS),
                   data = model_data, family = "poisson",glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 100000)))

car::Anova(SR_ND_crp1)
summary(SR_ND_crp1)

# likelihood ratio test
SR_ND_crp2<- glmer(Species_richness ~ NDRS +crpRS+ (1|SS) +(1|SSB)+(1|SSBS),
                   data = model_data, family = "poisson",glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 100000)))

anova(SR_ND_crp1,SR_ND_crp2)

# ** plot richness model--------------------------------------------- 
SR_ND_crp <- GLMER(modelData = model_data,responseVar = "Species_richness",
                   fitFamily = "poisson",
                   fixedStruct = "NDRS* crpRS",
                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")

#set different percentage of cropland
# 10% -0.3799434, 30% 0.6872444, 50% 1.753084,79% 3.279165

ncrp<- expand.grid(
  NDRS=seq(from = min(SR_ND_crp$data$NDRS),
           to = max(SR_ND_crp$data$NDRS),
           length.out = 100),
  crpRS=c(-0.3799434,0.6872444,1.753084,3.279165))

# back transform the climate data range
ncrp$ND <- BackTransformCentreredPredictor(
  transformedX = ncrp$NDRS,
  originalX = final.data.trans$NDt)

# back transform NH data range
ncrp$crp <- round(BackTransformCentreredPredictor(
  transformedX = ncrp$crpRS,originalX = final.data.trans$crp)*100,0)

# set values for richness and Abundance
ncrp$logAbun <- 0
ncrp$Species_richness <- 0

# set the reference row
refRow <- which((ncrp$ND==min(abs(ncrp$ND))) & 
                  (ncrp$crp==10))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)

QSR <- quantile(x = SR_ND_crp$data$NDRS,
                probs = exclQuantiles)
# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = SR_ND_crp$model,data = ncrp, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')


# set anything outside the desired quantiles to NA
a.preds.tmean[which((ncrp$NDRS  < QSR[1])),] <- NA
a.preds.tmean[which((ncrp$NDRS  > QSR[2])),] <- NA
# get the median and upper/lower intervals for plots
ncrp$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                           FUN = median,na.rm=TRUE))*100)-100
ncrp$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
ncrp$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
#ncrp$Predominant_land_use <- factor(ncrp$Predominant_land_use, levels = c("Primary vegetation", "Secondary vegetation", "Cropland", "Pasture","Plantation forest"))

ncrp$crp <- factor(ncrp$crp, levels = c("10", "30", "50","79"))

# just take agriculture values
#ncrp <- ncrp[ncrp$Predominant_land_use %in% c("Cropland"), ]

# plot
pncrp <- ggplot(data = ncrp, aes(x = ND, y = PredMedian)) + 
  geom_line(aes(col =crp), size = 0.75) +
  geom_ribbon(aes(ymin = ncrp$PredLower, ymax = ncrp$PredUpper, fill = crp), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#a50026","#397E90","#6CA684","#CBBA7B"))) +
  scale_colour_manual(values = rev(c("#a50026","#397E90","#6CA684","#CBBA7B"))) +
  #facet_wrap(~Predominant_land_use, ncol = 2, labeller = as_labeller(c('Secondary vegetation' = "a Secondary vegetation", 'Cropland' = "b Cropland",'Pasture' = "c Pasture",'Plantation forest' = "d Plantation forest"))) + 
  theme_bw() + 
  labs(fill = "PCR (%)", col = "PCR (%)") + 
  ylab("Change in species richness (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  xlim(c(0, 2000)) +
  ylim(c(-60, 30)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2)) 

pncrp
ggsave(pncrp,filename = "E:/Fig_initial/SR_crp_nd1.tiff",width = 7,height = 5, dpi = 600)


# ** abundance model--------------------------------------------- 
AB_ND_crp1<- glmer(logAbun ~NDRS * crpRS + (1|SS) +(1|SSB),
                   data = model_data, family = "gaussian")

car::Anova(AB_ND_crp$model)
summary(AB_ND_crp$model)

#likelihood ratio test
AB_ND_crp2<- glmer(logAbun ~NDRS + crpRS + (1|SS) +(1|SSB),
                   data = model_data, family = "gaussian")

anova(AB_ND_crp1,AB_ND_crp2)

# ** plot abundance model--------------------------------------------- 

AB_ND_crp <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
                   fixedStruct = "NDRS * crpRS",
                   randomStruct = "(1|SS)+(1|SSB)",
                   saveVars = c("SSBS"))

### extract crpRS values crp = 10%
data1 = final.data.trans[final.data.trans$crp == max(final.data.trans$crp[which(final.data.trans$crp <= 0.1)]),]

data1$crpRS

#set different percentage of cropland
# 10% -0.3799434, 30% 0.6872444, 50% 1.753084,79% 3.279165
##
ndab_crp<- expand.grid(
  NDRS=seq(from = min(AB_ND_crp$data$NDRS),
           to = max(AB_ND_crp$data$NDRS),
           length.out = 100),
  crpRS=c(-0.3799434,0.6872444,1.753084,3.279165))


# back transform the climate data range
ndab_crp$ND <- BackTransformCentreredPredictor(
  transformedX = ndab_crp$NDRS,
  originalX = final.data.trans$NDt)

# back transform NH data range
ndab_crp$crp <- round(BackTransformCentreredPredictor(
  transformedX = ncrp$crpRS,originalX = final.data.trans$crp)*100,0)

# set values for richness and Abundance
ndab_crp$logAbun <- 0
ndab_crp$Species_richness <- 0

# set the reference row
refRow <- which((ndab_crp$ND==min(abs(ndab_crp$ND))) & 
                  (ndab_crp$crp==10))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)

QAB <- quantile(x = AB_ND_crp$data$NDRS,
                probs = exclQuantiles)
# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AB_ND_crp$model,data = ndab_crp, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-1

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')


# set anything outside the desired quantiles to NA
a.preds.tmean[which((ndab_crp$NDRS  < QAB[1])),] <- NA
a.preds.tmean[which((ndab_crp$NDRS  > QAB[2])),] <- NA
# get the median and upper/lower intervals for plots
ndab_crp$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                               FUN = median,na.rm=TRUE))*100)-100
ndab_crp$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
ndab_crp$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
ndab_crp$crp <- factor(ndab_crp$crp, levels = c("10", "30", "50","79"))

# just take agriculture values
#ndab_crp <- ndab_crp[ndab_crp$Predominant_land_use %in% c("Cropland"), ]

# plotshijie shi
pndab_crp <- ggplot(data = ndab_crp, aes(x = ND, y = PredMedian)) + 
  geom_line(aes(col =crp), size = 0.75) +
  geom_ribbon(aes(ymin = ndab_crp$PredLower, ymax = ndab_crp$PredUpper, fill = crp), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#313695","#4e72b8"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#313695","#4e72b8")))+
  theme_bw() + 
  labs(fill = "PCR (%)", col = "PCR (%)") + 
  ylab("Change in total abundance (%)") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')'))))+
  xlim(c(0, 2000)) +
  ylim(c(-80, 30)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2)) 

pndab_crp

ggsave(pndab_crp,filename = "E:/Fig_initial/AB_crp_nd.tiff",width = 7,height = 5, dpi = 600)

# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * (5)	Effect of the interaction between SNH and cropland cover on N deposition impacts -----------------------------------------------
# * Model N deposition/ natural habitat/cropland-----------------------------------------------

# ** richness model--------------------------------------------- 
SR_pnh_crp1<- glmer(Species_richness ~ NDRS*pnhRS*crpRS + (1|SS) +(1|SSB)+(1|SSBS),
                    data = model_data, family = "poisson",glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
car::Anova(SR_pnh_crp1)
summary(SR_pnh_crp1)


#likelihood ratio test
SR_pnh_crp2<- glmer(Species_richness ~ NDRS+pnhRS+crpRS+NDRS:pnhRS+crpRS:pnhRS+NDRS:crpRS+ (1|SS) +(1|SSB)+(1|SSBS),
                    data = model_data, family = "poisson",glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

anova(SR_pnh_crp1,SR_pnh_crp2)

# ** plot richness model--------------------------------------------- 

#set cropland and pnh as class variables
crp_SR_ND <-ggeffect(SR_pnh_crp1,terms = c("NDRS[all]","pnhRS[-1.137128,-0.1354202,0.8654312,1.863949]","crpRS[-0.3799434,0.1520204,0.6872444,1.220082,1.753084,3.279165]"),ci.lvl = 0.95)
summary(crp_SR_ND)

##
crp_SR_ND$ND <- BackTransformCentreredPredictor(
  transformedX = crp_SR_ND$x,
  originalX = final.data.trans$NDt)

##
pcrp_SR_ND <- ggplot(da=crp_SR_ND, aes(x = ND, y = predicted)) + 
  geom_line(aes(col =group), size = 0.75) +
  geom_ribbon(aes(ymin = crp_SR_ND$conf.low, ymax = crp_SR_ND$conf.high, fill = group), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#4e72b8","#74add1","#FFCB48","#A71111")),labels = c("25","50","75","100")) +
  scale_colour_manual(values = rev(c("#4e72b8","#74add1","#FFCB48","#A71111")),labels = c("25","50","75","100")) +
  facet_wrap(~facet, ncol = 3, labeller = as_labeller(c('-0.3799434' = "10%cropland", '0.1520204' = "20%cropland",'0.6872444' = "30%cropland",'1.220082' = "40%cropland",'1.753084' = "50%cropland",'3.279165' = "79%cropland"))) + 
  theme_bw() + 
  labs(fill = "PNSH (%)", col = "PNSH (%)") + 
  ylab("Number of Species richness") +
  xlab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')')))) +
  xlim(c(0, 4000)) +
  ylim(c(0, 25)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13.5),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2),
        panel.spacing = unit(1,"cm"))

pcrp_SR_ND

ggsave(pcrp_SR_ND,filename = "E:/Fig_initial/pcrp_SR_ND_phn.tiff",width = 5,height = 5, dpi = 600)


# ** abundance model--------------------------------------------- 
AB_pnh_crp1<- glmer(logAbun ~ NDRS*pnhRS*crpRS+ (1|SS) +(1|SSB),
                    data = model_data, family = "gaussian")

car::Anova(AB_pnh_crp1)
summary(AB_pnh_crp1)

#likelihood ratio test
AB_pnh_crp2<- glmer(logAbun ~ NDRS+pnhRS+crpRS+NDRS:pnhRS+crpRS:pnhRS+NDRS:crpRS+(1|SS) +(1|SSB),
                    data = model_data, family = "gaussian")

anova(AB_pnh_crp1,AB_pnh_crp2)


# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------
# * test 0.1 ° N deposition dataset-----------------------------------------

#read 0.1° N deposition raster
oxnrdn <- raster(paste0("D:/R/variables/oxnrdn_01grid.tif"))

#extract values with sites
final.data.trans$oxnrdn<- extract(oxnrdn , final.data.trans_xy, na.rm = FALSE)
# where there are NAs, change to 0
final.data.trans[is.na(final.data.trans$oxnrdnRS), "oxnrdnRS"] <- 0

# remove any rows that have NA in the variable columns
# remove rows with NAs for any variable of interest
sites.sub <- sites.sub[!is.na(sites.sub$oxnrdnRS), ] # 

# standardize variables
final.data.trans1$oxnrdnRS <- scale(final.data.trans1$oxnrdn)

#Pearson correlation
cor.test(final.data.trans1$oxnrdn,final.data.trans1$NDt) #0.878

# ** richness model--------------------------------------------- 


test_SRmodel <- GLMER(modelData = final.data.trans1,responseVar = "Species_richness",fitFamily = "poisson",
                      fixedStruct = "oxnrdnRS",
                      randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")

car::Anova(test_SRmodel$model)
summary(test_SRmodel$model)

#output results
tab_model(test_SRmodel$model, transform = NULL, file = paste0("E:/OneDrive/!文件/PREDICT/Table_initial/Table_01SR_ND.html"),
          show.icc = F, show.obs = T, show.ngroups = F)
# ** abundance model----------------------------------------- 
test_ABmodel <- GLMER(modelData = final.data.trans1,responseVar = "logAbun",fitFamily = "gaussian",
                      fixedStruct = "oxnrdnRS",
                      randomStruct = "(1|SS)+(1|SSB)",
                      saveVars = c("SSBS"))
car::Anova(test_ABmodel$model)
summary(test_ABmodel$model)

#output results
tab_model(test_ABmodel$model, transform = NULL, file = paste0("E:/Table_initial/Table_01AB_ND.html"),
          show.icc = F, show.obs = T, show.ngroups = F)


# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------

# * N deposition values across land use et al.-----------------------------------------
library("agricolae")
library(ggsignif)

# ** anova,landuse,nd----------------------------------------------
oneway<-aov(model_data$NDt~model_data$Predominant_land_use,data = model_data)
anova(oneway)

out <- LSD.test(oneway,"model_data$Predominant_land_use",p.adj="none")
summary(out)

pndlanduse<-ggplot(model_data,aes(x=Predominant_land_use,y=NDt))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=Predominant_land_use), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  #coord_flip()+
  scale_y_continuous( limits = c(0, 3000)) +
  annotate("text",x=1,y=1500,label="e",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=2,y=1900,label="c",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=3,y=2600,label="a",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=4,y=2100,label="b",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=5,y=1900,label="d",family = "serif", fontface = "italic",size=5)+
  scale_fill_manual(values = c("#11527F","#892F1F", "#A76622","#14623B","#575982")) +
  scale_colour_manual(values = c("#11527F","#892F1F", "#A76622","#14623B","#575982")) +
  xlab("Land use types") +
  ylab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')'))))+
  scale_x_discrete(breaks = c('Primary vegetation','Secondary vegetation','Cropland','Plantation forest','Pasture'),
                   label = c('Primary\nvegetation','Secondary\nvegetation','Cropland','Plantation\nforest','Pasture'))+
  stat_summary(fun="mean", geom="point", shape=20, size=2.5, color="red", fill="red",alpha=0.7) +
  theme(aspect.ratio = 1, 
        title = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position="none",
        panel.background =element_blank(),
        panel.border = element_rect(size = 0.2, fill=NA), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))

pndlanduse 

ggsave(pndlanduse,filename = "E:/nd_landuse.tiff",width = 5,height = 5, dpi = 600)
# ** anova,Continental,nd----------------------------------------------
oneway<-aov(model_data$NDt~model_data$UN_region,data = model_data)
anova(oneway)

out <- LSD.test(oneway,"model_data$UN_region",p.adj="none")
out
summary(out)

pndregion<-ggplot(model_data,aes(x=UN_region,y=NDt))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=UN_region), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  scale_y_continuous( limits = c(0, 3000)) +
  annotate("text",x=1,y=1500,label="d",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=2,y=1900,label="b",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=3,y=1600,label="c",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=4,y=2800,label="a",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=5,y=600,label="e",family = "serif", fontface = "italic",size=5)+
  scale_fill_manual(values = c("#11527F","#892F1F", "#A76622","#14623B","#575982")) +
  scale_colour_manual(values = c("#11527F","#892F1F", "#A76622","#14623B","#575982")) +
  xlab("Continents") +
  ylab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')'))))+
  theme(aspect.ratio = 1, 
        title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position="none",
        panel.background =element_blank(),
        panel.border = element_rect(size = 0.2, fill=NA), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))

pndregion 

ggsave(pndregion,filename = "E:/Fig_initial/nd_Continent.tiff",width = 5,height = 5, dpi = 600)

# ** anova,climate zone,nd---------------------------------------------
oneway<-aov(model_data4$NDt~model_data4$climatezone1,data = model_data4)
anova(oneway)

out <- LSD.test(oneway,"model_data4$climatezone1",p.adj="none")
summary(out)

#
pndclimate<-ggplot(model_data4,aes(x=climatezone4,y=NDt))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=climatezone1), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  scale_y_continuous( limits = c(0, 3000)) +
  annotate("text",x=1,y=1500,label="e",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=2,y=1900,label="d",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=3,y=1600,label="f",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=4,y=1000,label="c",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=5,y=2500,label="b",family = "serif", fontface = "italic",size=5)+
  annotate("text",x=6,y=2200,label="a",family = "serif", fontface = "italic",size=5)+
  scale_fill_manual(values = c("#11527F","#892F1F", "#A76622","#14623B","#575982","#547282")) +
  scale_colour_manual(values = c("#11527F","#892F1F", "#A76622","#14623B","#575982","#547282")) +
  xlab("Climate Zones") +
  ylab(expression(bold(paste('Nitrogen deposition (mg N m'^-2,' yr'^-1,')'))))+
  theme(aspect.ratio = 1, 
        title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position="none",
        panel.background =element_blank(),
        panel.border = element_rect(size = 0.2, fill=NA), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))

pndclimate 

ggsave(pndclimate,filename = "E:/Fig_initial/nd_climate.tiff",width = 5,height = 5, dpi = 600)

# ** plot all----------------------------------------------

panova<-cowplot::plot_grid(pndregion,pndclimate,pndlanduse,
                           labels = c("A.", "B.", "C."), nrow = 2)

ggsave(panova,filename = "E:/Fig_initial/panova.tiff",width = 10,height = 10, dpi = 600)

# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------

# * Model result check-----------------------------------------

#Tim's roquefort function for spatial autocorrelation
SpatialAutocorrelationTest<-function(model,all.data,siteModel=TRUE){
  
  if ("data" %in% names(model)){
    model.data<-model$data
    model.data$res<-residuals(model$model)
  } else {
    model.data<-model@frame
    model.data$res<-residuals(model)
  }
  
  model.data$Longitude<-all.data$Longitude[match(model.data$SSBS,all.data$SSBS)]
  model.data$Latitude<-all.data$Latitude[match(model.data$SSBS,all.data$SSBS)]
  
  studies<-character()
  failed<-character()
  moran.i<-numeric()
  moran.p<-numeric()
  
  i=1
  for (ss in unique(model.data$SS)){
    cat(paste("\rProcessing study ",i," of ",length(unique(model.data$SS)),sep=""))
    data.sub<-droplevels(model.data[model.data$SS==ss,])
    
    if(!siteModel){
      resids<-tapply(data.sub$res,data.sub$SSBS,mean)
      long<-tapply(data.sub$Longitude,data.sub$SSBS,mean)
      lat<-tapply(data.sub$Latitude,data.sub$SSBS,mean)
      data.sub<-data.frame(res=resids,Longitude=long,Latitude=lat)
    }
    
    ds.nb<-try(dnearneigh(cbind(data.sub$Longitude,data.sub$Latitude),
                          d1=0.00000001,d2=10),silent=TRUE)
    ds.listw<-try(nb2listw(ds.nb),silent=TRUE)
    mt<-tryCatch(moran.test(data.sub$res,ds.listw),silent=TRUE,error=function(e) e, 
                 warning=function(w) w)
    
    if(class(mt)[1]=="htest"){
      if ((!is.na(mt$statistic))){
        studies<-c(studies,ss)
        moran.i<-c(moran.i,mt$statistic)
        moran.p<-c(moran.p,mt$p.value)
      } else {
        failed<-c(failed,ss)
      }
      
    } else {
      failed<-c(failed,ss)
    }
    
    i<-i+1
  }
  
  return(list(studies=studies,I=moran.i,P=moran.p,failed=failed))
  
}

# ** model list----------------------------------------------

# AB_ND , SR_ND , AB_ND_land , SR_ND_land , AB_ND_region , SR_ND_region , AB_ND_class, SR_ND_class , AB_ND_climate , SR_ND_climate
#AB_ND_PNH , SR_ND_PNH , AB_ND_crp, SR_ND_crp ,  AB_pnh_crp , SR_pnh_crp

# ** Model check SR----------------------------------------------

p1 <- plot(SR_ND$model)

## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(SR_ND$model), main = "")
qqline(resid(SR_ND$model))
p2 <- recordPlot()
invisible(dev.off())
## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(SR_ND, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$Species_richness ,fitted(SR_ND$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"), align = "v", axis = "bt",rel_heights= c(0.5, 1))
psr

psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"))

psr
#
rm(p1, p2, p3, p4)
rm(perc_auto)

# ** Model check AB----------------------------------------------

p1 <- plot(AB_ND$model)
## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(AB_ND$model), main = "")
qqline(resid(AB_ND$model))
p2 <- recordPlot()
invisible(dev.off())

## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(AB_ND, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$logAbun ,fitted(AB_ND$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

pab <-cowplot::plot_grid(p2,p3,
                         labels = c("C", "D"))
pab

rm(p1, p2, p3, p4)
rm(perc_auto)

#
cowplot::plot_grid(psr,pab,nrow = 2,rel_heights = c(1.5,1.5,1,1))

cowplot::plot_grid(p1, p2 +  theme(legend.position = "none"), legend,  nrow = 3, rel_heights = c(1,1,0.2))


ggsave(filename = "E:/Fig_initial/model_checks.pdf", height = 10, width = 10)

ggsave(filename = "E:/Fig_initial/model_checks.tiff", height = 10, width = 10)


# ** Model check SR land----------------------------------------------

p1 <- plot(SR_ND_land$model)

## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(SR_ND_land$model), main = "")
qqline(resid(SR_ND_land$model))
p2 <- recordPlot()
invisible(dev.off())
## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(SR_ND_land, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$Species_richness ,fitted(SR_ND_land$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

#psr <-cowplot::plot_grid(p2,p3,
#                         labels = c("A", "B"), align = "v", axis = "bt",rel_heights= c(0.5, 1))


psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"))

psr
#
rm(p1, p2, p3, p4)
rm(perc_auto)

# ** Model check AB land----------------------------------------------

p1 <- plot(AB_ND_land$model)
## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(AB_ND_land$model), main = "")
qqline(resid(AB_ND_land$model))
p2 <- recordPlot()
invisible(dev.off())

## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(AB_ND_land, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$logAbun ,fitted(AB_ND_land$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

pab <-cowplot::plot_grid(p2,p3,
                         labels = c("C", "D"))
pab

rm(p1, p2, p3, p4)
rm(perc_auto)

#
cowplot::plot_grid(psr,pab,nrow = 2,rel_heights = c(1.5,1.5,1,1))


ggsave(filename = "E:/Fig_initial/model_checks_land.pdf", height = 10, width = 10)

ggsave(filename = "E:/Fig_initial/model_checks_land.tiff", height = 10, width = 10)


# ** Model check SR pnh----------------------------------------------

p1 <- plot(SR_ND_PNH$model)

## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(SR_ND_PNH$model), main = "")
qqline(resid(SR_ND_PNH$model))
p2 <- recordPlot()
invisible(dev.off())
## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(SR_ND_PNH, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$Species_richness ,fitted(SR_ND_PNH$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

#psr <-cowplot::plot_grid(p2,p3,
#                         labels = c("A", "B"), align = "v", axis = "bt",rel_heights= c(0.5, 1))


psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"))

psr
#
rm(p1, p2, p3, p4)
rm(perc_auto)

# ** Model check AB pnh----------------------------------------------

p1 <- plot(AB_ND_PNH$model)
## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(AB_ND_PNH$model), main = "")
qqline(resid(AB_ND_PNH$model))
p2 <- recordPlot()
invisible(dev.off())

## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(AB_ND_PNH, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$logAbun ,fitted(AB_ND_PNH$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

pab <-cowplot::plot_grid(p2,p3,
                         labels = c("C", "D"))
pab

rm(p1, p2, p3, p4)
rm(perc_auto)

#
cowplot::plot_grid(psr,pab,nrow = 2,rel_heights = c(1.5,1.5,1,1))


ggsave(filename = "E:/Fig_initial/model_checks_pnh.pdf", height = 10, width = 10)

ggsave(filename = "E:/Fig_initial/model_checks_pnh.tiff", height = 10, width = 10)

# ** Model check SR crp----------------------------------------------

p1 <- plot(SR_ND_crp$model)

## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(SR_ND_crp$model), main = "")
qqline(resid(SR_ND_crp$model))
p2 <- recordPlot()
invisible(dev.off())
## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(SR_ND_crp, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$Species_richness ,fitted(SR_ND_crp$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

#psr <-cowplot::plot_grid(p2,p3,
#                         labels = c("A", "B"), align = "v", axis = "bt",rel_heights= c(0.5, 1))


psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"))

psr
#
rm(p1, p2, p3, p4)
rm(perc_auto)

# ** Model check AB crp----------------------------------------------

p1 <- plot(AB_ND_crp$model)
## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(AB_ND_crp$model), main = "")
qqline(resid(AB_ND_crp$model))
p2 <- recordPlot()
invisible(dev.off())

## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(AB_ND_crp, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$logAbun ,fitted(AB_ND_crp$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

pab <-cowplot::plot_grid(p2,p3,
                         labels = c("C", "D"))
pab

rm(p1, p2, p3, p4)
rm(perc_auto)

#
cowplot::plot_grid(psr,pab,nrow = 2,rel_heights = c(1.5,1.5,1,1))


ggsave(filename = "E:/Fig_initial/model_checks_crp.pdf", height = 10, width = 10)

ggsave(filename = "E:/Fig_initial/model_checks_crp.tiff", height = 10, width = 10)

# ** Model check SR crp pnh----------------------------------------------

p1 <- plot(SR_pnh_crp$model)

## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(SR_pnh_crp$model), main = "")
qqline(resid(SR_pnh_crp$model))
p2 <- recordPlot()
invisible(dev.off())
## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(SR_pnh_crp, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$Species_richness ,fitted(SR_pnh_crp$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

#psr <-cowplot::plot_grid(p2,p3,
#                         labels = c("A", "B"), align = "v", axis = "bt",rel_heights= c(0.5, 1))


psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"))

psr
#
rm(p1, p2, p3, p4)
rm(perc_auto)

# ** Model check AB crp pnh----------------------------------------------

p1 <- plot(AB_pnh_crp$model)
## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(AB_pnh_crp$model), main = "")
qqline(resid(AB_pnh_crp$model))
p2 <- recordPlot()
invisible(dev.off())

## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(AB_pnh_crp, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$logAbun ,fitted(AB_pnh_crp$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

pab <-cowplot::plot_grid(p2,p3,
                         labels = c("C", "D"))
pab

rm(p1, p2, p3, p4)
rm(perc_auto)

#
cowplot::plot_grid(psr,pab,nrow = 2,rel_heights = c(1.5,1.5,1,1))


ggsave(filename = "E:/Fig_initial/model_checks_pnh_crp.pdf", height = 10, width = 10)

ggsave(filename = "E:/Fig_initial/model_checks_pnh_crp.tiff", height = 10, width = 10)

# ** Model check SR 0.1 ----------------------------------------------

p1 <- plot(SR_ND$model)

## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(SR_ND$model), main = "")
qqline(resid(SR_ND$model))
p2 <- recordPlot()
invisible(dev.off())
## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(SR_ND, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$Species_richness ,fitted(SR_ND$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"), align = "v", axis = "bt",rel_heights= c(0.5, 1))
psr

psr <-cowplot::plot_grid(p2,p3,
                         labels = c("A", "B"))

psr
#
rm(p1, p2, p3, p4)
rm(perc_auto)

# ** Model check AB 0.1----------------------------------------------

p1 <- plot(AB_ND$model)
## 2. Normality of Residuals
pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(AB_ND$model), main = "")
qqline(resid(AB_ND$model))
p2 <- recordPlot()
invisible(dev.off())

## 3. Check for spatial autocorrelation
sa_test<-SpatialAutocorrelationTest(AB_ND, model_data)

perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100


sa_test_vals <- as.data.frame(sa_test$P)
sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sa_test_vals ) +
  geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
# 4. plot of observed vs fitted values
pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$logAbun ,fitted(AB_ND$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p4 <- recordPlot()
invisible(dev.off())
#logAbun #Species_richness

pab <-cowplot::plot_grid(p2,p3,
                         labels = c("C", "D"))
pab

rm(p1, p2, p3, p4)
rm(perc_auto)

#
cowplot::plot_grid(psr,pab,nrow = 2,rel_heights = c(1.5,1.5,1,1))

cowplot::plot_grid(p1, p2 +  theme(legend.position = "none"), legend,  nrow = 3, rel_heights = c(1,1,0.2))


ggsave(filename = "E:/Fig_initial/model_checks.pdf", height = 10, width = 10)

ggsave(filename = "E:/Fig_initial/model_checks.tiff", height = 10, width = 10)


# * ★★★★★★★★★★★★★★★★★★---------------------------------------------------------

# * END ----------------------------------------


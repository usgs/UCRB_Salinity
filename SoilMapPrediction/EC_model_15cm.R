######################
## Random Forest script that includes:
## Extraction of covariates to points
## Confustion matrix creation
## Kappa calculation
## Prediction interval creation
## Cross Validation
## Most steps parallelized
######################



# Workspace setup
# Install packages if not already installed

required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall", "quantregForest","hydroGOF")# might need snowfall
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
memory.limit(500000)
rasterOptions(maxmemory = 1e+9)
options(scipen = 999)
par(mar=c(0.3,0.3,0.3,0.3))

######## Load shapefile ##############
setwd("V:/PROJECTS/TRAVIS_NAUMAN/GIS_Archive/NCSS_Soil_Characterization_Database_2017/EC_wLIMS")## FOlder with points
shp.pts <-readOGR(".", "NCSS17_ec_LIMS_5m")
point.proj <- projection(shp.pts)
## If not prj file and you know proj, can specify by name
shp.proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##Remove unneccesary attribution
shp.pts <-(shp.pts[c(3,6,29,30,52,59,63,65,69,73,74)])
pts<-shp.pts

##Create Fields to hold merged columns
shp.pts$Sample_ID<-NA
shp.pts$LONG_dec<-NA
shp.pts$LAT_dec<-NA
shp.pts$Hzn_depth<-NA
shp.pts$Pedon_key<-NA

##ID fields to Character
shp.pts$smp_id<-as.character(shp.pts$smp_id)
shp.pts$labsampnum<-as.character(shp.pts$labsampnum)
shp.pts$user_pedon<-as.character(shp.pts$user_pedon)
shp.pts$pedon_key<-as.character(shp.pts$pedon_key)

##Merge Fields
shp.pts$LAT_dec<-coalesce(shp.pts$latitude_d,shp.pts$lat_dec)
shp.pts$LONG_dec<-coalesce(shp.pts$longitude_,shp.pts$long_dec)
shp.pts$Sample_ID<-coalesce(shp.pts$labsampnum,shp.pts$smp_id)
shp.pts$Pedon_key<-coalesce(shp.pts$pedon_key,shp.pts$user_pedon)

##As numeric
shp.pts$ec_12pre<-as.numeric(as.character(shp.pts$ec_12pre))
shp.pts$depth_t<-as.numeric(as.character(shp.pts$depth_t))
shp.pts$depth_b<-as.numeric(as.character(shp.pts$depth_b))

#calculate depth
shp.pts$Hzn_depth <- (shp.pts$depth_b-shp.pts$depth_t)/2

##Clean Attribution Again
shp.pts<-shp.pts[,c(5,10,11,12,13,14,15,16)]
pts1<-shp.pts

### Weed out points with imprecise coordinates ###
shp.pts$latnchar = nchar(abs(shp.pts$LAT_dec))
shp.pts$longnchar = nchar(abs(shp.pts$LONG_dec))
ptsc = subset(shp.pts, shp.pts$latnchar > 5 & shp.pts$longnchar > 6)
shp.pts<-ptsc

######## Load map clip boundary (if needed) ###########
setwd("V:/PROJECTS/TRAVIS_NAUMAN/GIS_Archive/CO_River_watershed_Meade")
polybound <- readOGR(".", "CO_River_watershed_Meade_ACEA")
polybound <- spTransform(polybound, shp.proj)
## Now clip points and check with visualization
shp.pts = shp.pts[polybound,]#clip by outer extent of all polybound features
plot(polybound)
plot(shp.pts, add=TRUE)

######### Grid Prep #################
## Make list of grids
setwd("O:/Models_active_work/COP_Covariates")
cov.grids <- list.files(pattern=".tif$")
## If points need to be matched up to grids ###
projgrid = raster("LFelems.tif")
## Or Make a stack of grids to extract all at once (for smaller datasets)
#cov.stack <- stack()
cov.proj <- projection(projgrid)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match rasters

## Plot to ensure alignment bw points and rasters
plot(projgrid)
plot(shp.pts, add=TRUE)

## Parallelized extract: (larger datasets)
cpus = 30
sfInit(parallel=TRUE, cpus=cpus)
sfExport("shp.pts", "cov.grids")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(cov.grids, function(i){try( raster::extract(raster(i), shp.pts) )}) 
snowfall::sfStop()
detach(package:snowfall, unload=TRUE)
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) = tools::file_path_sans_ext(basename(cov.grids))
ov.lst$DID <- seq.int(nrow(ov.lst))
shp.pts$DID <- seq.int(nrow(shp.pts))
pts.ext <- merge(as.data.frame(shp.pts),ov.lst, by="DID")

## Save points
setwd("O:/Models_active_work/UpCo/ECmodel_wLIMS_2D")
write.table(pts.ext, "cop_ncss17_ec12_2D_covarsc.txt", sep = "\t", row.names = FALSE)

## Prep for Random Forest
pts.extc <- subset(pts.ext, as.numeric(pts.ext$depth_t) <= 15 & as.numeric(pts.ext$depth_b) > 15) # subset to chosen depth
ptspred.list <- gsub(".tif","", cov.grids)# Take .tif off of the grid list to just get the variable names
pred <- "ec_12pre" ## Dependent variable
ptspred.list <- c(ptspred.list,pred) #Add dependent variable
pts.extc <- pts.extc[c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use 
pts.extc <- na.omit(pts.extc)# Remove any record with NA's (in any column - be careful)
pts.extc  <- subset(pts.extc, pts.extc$ec_12pre != "NA")
xtrain <- as.matrix(pts.extc[c(gsub(".tif","", cov.grids))])
ytrain <- c(as.matrix(pts.extc[c(pred)]))
logytrain <- log(ytrain)

############### Build quantile Random Forest
Qsoiclass <- quantregForest(x=xtrain, y=logytrain, importance=TRUE, ntree=100, keep.forest=TRUE)
#soiclass = randomForest(ec_12pre ~ ., data = ptsc, importance=TRUE, proximity=FALSE, ntree=100, keep.forest=TRUE)
soiclass = Qsoiclass
class(soiclass) = "randomForest"
soiclass## Get oob error
varImpPlot(soiclass)
############### Create Confusion matrix (for categorical models)
## Need to strip last column e.g. confusion[1:9,1:10] in rf object would be confusion[1:9,1:9]
#confmatrix = soiclass$confusion[1:9,1:9] 
#OOBkappa = Kappa.test(confmatrix, conf.level = 0.95)
#write.table(confmatrix, file = "MLRA35_ESG_oob_RFSoilGrids100m_conf_matr.txt", sep = "/t")


## Reference covar rasters to use in prediction
setwd("O:/Models_active_work/UCRB_Covariates")
rasters=stack(cov.grids)
#rasters = setMinMax(brick(rasters))
names(rasters)

## Predict onto covariate grid
setwd("O:/Models_active_work/UpCo/ECmodel_wLIMS_2D/outputs")
## Parallelized predict
beginCluster(31,type='SOCK')
predl = clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.025)),progress="text")
predh = clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.975)),progress="text")
pred = clusterR(rasters, predict, args=list(model=soiclass),progress="text")
endCluster()
#s = stack(predh,predl)
#PIwidth = overlay(s, fun=function(a,b) (a-b),progress = "text")
#varrange = as.numeric(quantile(log(ptsc$ec_12pre), probs=c(0.975))-quantile(log(ptsc$ec_12pre), probs=c(0.025)))
#PIrelwidth = overlay(s, fun=function(a,b) ((a-b)/varrange), progress = "text")
predh_bt = calc(predh, fun=function(x) (exp(x)), progress="text")#If a backtransform is needed 10^(x) or exp(x)
predl_bt = calc(predl, fun=function(x) (exp(x)), progress="text")
pred_bt = calc(pred, fun=function(x) (exp(x)), progress="text")
s_bt = stack(predh_bt,predl_bt)
PIwidth_bt = overlay(s_bt, fun=function(a,b) (a-b),progress = "text")
varrange_bt = as.numeric(quantile(pts.extc$ec_12pre, probs=c(0.975))-quantile(pts.extc$ec_12pre, probs=c(0.025)))
PIrelwidth_bt = overlay(s_bt, fun=function(a,b) ((a-b)/varrange_bt), progress = "text")
writeRaster(pred_bt, overwrite=F,filename="EC_2D_15cm_QRF_bt.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(pred, overwrite=F,filename="EC_2D_15cm_QRF.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
writeRaster(predl_bt, overwrite=F,filename="EC_2D_15cm_QRF_95PI_l_bt.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
writeRaster(predh_bt, overwrite=F,filename="EC_2D_15cm_QRF_95PI_h_bt.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(predl, overwrite=F,filename="EC_2D_15cm_QRF_95PI_l.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(predh, overwrite=F,filename="EC_2D_15cm_QRF_95PI_h.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(PIrelwidth, overwrite=F,filename="EC_2D_15cm_QRF_95PI_relwidth.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(PIwidth, overwrite=F,filename="EC_2D_15cm_QRF_95PI_width.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
writeRaster(PIwidth_bt, overwrite=F,filename="EC_2D_15cm_QRF_95PI_width_bt.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
writeRaster(PIrelwidth_bt, overwrite=F,filename="EC_2D_15cm_QRF_95PI_relwidth_bt.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
## Create lookup table (for categorical predictions)
#lookup_tab = as.data.frame(soiclass$classes)
#write.table(lookup_tab, file = "ESG_MLRA35_lookup_tab.txt", sep = "/t")


##################### Run Cross Validation ######################################
trainx = gsub(".tif","", cov.grids)
trainx = c(trainx)
trainx = pts.extc[c(trainx)]
trainy = pts.extc[c("ec_12pre")]
set.seed(41)
## rfcv runs a cross val and tests variable importance
nfolds <- 10
cv.rf = rfcv(trainx, log(trainy$ec_12pre), cv.fold = nfolds, proximity=FALSE, ntree=100, keep.forest=TRUE)
par(mar=c(5, 4, 4, 2) + 0.1)
with(cv.rf, plot(n.var, error.cv, log="x", type="o", lwd=2))
# Now put the full cross val out of the rfcv
pts.extc$cvpred = cv.rf$predicted$`40` ## the '#' at end corresponds to the number of variables included should = number of vars used
cv.RMSE = sqrt(mean((log(pts.extc$ec_12pre) - pts.extc$cvpred)^2, na.rm=TRUE))
cv.Rsquared = 1-var(log(pts.extc$ec_12pre) - pts.extc$cvpred, na.rm=TRUE)/var(log(pts.extc$ec_12pre), na.rm=TRUE)
NSE(pts.extc$cvpred,log(pts.extc$ec_12pre))


################### Manual Cross validation ################################

ptspred.listcvm <- c(ptspred.list)
pts.extcvm <- pts.ext[c(ptspred.listcvm)]
pts.extcvm <- na.omit(pts.extcvm)# Remove any record with NA's (in any column - be careful)
pts.extcvm  <- subset(pts.extcvm, pts.extcvm$ec_12pre != "NA")
nfolds <- 10
pts.extcvm$folds <- sample.int(nfolds,size =length(pts.extcvm$ec_12pre),replace=T)
formulaStringCVm <- as.formula(paste('log(ec_12pre) ~','depth','+', paste(gsub(".tif","", cov.grids), collapse="+")))
pts.extcvm$mcvpred <- "NA"
for (g in seq(nfolds)){
  traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
  rf.pcvm <- randomForest(formulaStringCVm, data=traindf, importance=FALSE, proximity=FALSE, ntree=100, keep.forest=TRUE)
  pts.extcvm$mcvpred <- ifelse(pts.extcvm$folds == g, predict(rf.pcvm, newdata=pts.extcvm),pts.extcvm$mcvpred)
  print(g)
}
pts.extcvm$mcvpred = as.numeric(pts.extcvm$mcvpred)
cvm.RMSE = sqrt(mean((pts.extcvm$ec_12pre - pts.extcvm$)^2, na.rm=TRUE))
cvm.Rsquared = 1-var(log(pts.extcvm$ec_12pre) - pts.extcvm$mcvpred, na.rm=TRUE)/var(log(pts.extcvm$ec_12pre), na.rm=TRUE)




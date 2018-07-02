######################
## Random Forest script that includes:
## Extraction of covariates to points
## Confustion matrix creation
## Kappa calculation
## Prediction interval creation
## Cross Validation
## Most steps parallelized
######################

rm(list=ls())

# Workspace setup
# Install packages if not already installed

required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall", "quantregForest")# might need snowfall
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
rasterOptions(maxmemory = 1e+09)
memory.limit(500000)
par(mar=c(0.3,0.3,0.3,0.3))


######## Load shapefile ##############
setwd("O:/Models_active_work/UpCo/Kw")## FOlder with points
   ##### Load Points ######
shp.pts <-read.table("NCSSlabpeds_161220_UpCo_kfact.txt",header=TRUE)
   #####convert to spatial dataframe######
shp.pts<-SpatialPointsDataFrame(shp.pts[,1:2],shp.pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
#crs(shp.pts)<-"+proj=longlat +datum=WGS84"
point.proj <- projection(shp.pts)
## If not prj file and you know proj, can specify by name
#shp.proj <- CRS("+proj=longlat +datum=WGS84")

######### Grid Prep #################
## Make list of grids
setwd("O:/Models_active_work/UpCo/final_covars")
cov.grids <- list.files(pattern=".tif$")
## If points need to be matched up to grids ###
projgrid = raster("redblue.tif")
## Or Make a stack of grids to extract all at once (for smaller datasets)
#cov.stack <- stack()
cov.proj <- projection(projgrid)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match rasters
#writeOGR(shp.pts, ".", "k_sw_pts_prob", driver="ESRI Shapefile")

## Plot to ensure alignment bw points and rasters
plot(projgrid)
plot(shp.pts, add=TRUE)

## Parallelized extract: (larger datasets)
cpus = 31
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
shp.pts$DID = seq.int(nrow(shp.pts))
pts = merge(as.data.frame(shp.pts),ov.lst, by="DID")

## Save points
setwd("O:/Models_active_work/UpCo/Kw")
write.table(pts, "UPCO_NCSSlabpeds_Kw_covars.txt", sep = "\t", row.names = FALSE)
#pts$swregap = as.factor(pts$swregap)
#pts$LFelems  = as.factor(as.character(pts$LFelems))
#pts$pscsmodalb = as.factor(pts$pscsmodalb)

## Prep for Random Forest
ptspred.list = gsub(".tif","", cov.grids)# Take .tif off of the grid list to just get the variable names
pred = "K_sw" ## Dependent variable
DEPTH = "DEPTH" # to use in 3D soil property mapping
ptspred.list = c(ptspred.list,pred,DEPTH) #Add dependent variable
ptsc = pts[c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use 
ptsc$DEPTH = as.numeric(as.character(ptsc$DEPTH))
ptsc = na.omit(ptsc)# Remove any record with NA's (in any column - be careful)
ptsc  = subset(ptsc, ptsc$K_sw != "NA")
xtrain = as.matrix(ptsc[c("DEPTH",gsub(".tif","", cov.grids))])
ytrain = as.matrix(ptsc[c(pred)])
#logytrain = log(ytrain)

############### Build quantile Random Forest
Qsoiclass = quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=100, keep.forest=TRUE)
#soiclass = randomForest(ec_12pre ~ ., data = ptsc, importance=TRUE, proximity=FALSE, ntree=100, keep.forest=TRUE)
soiclass = Qsoiclass
class(soiclass) = "randomForest"
soiclass## Get oob error
varImpPlot(soiclass)
############### Create Confusion matrix (for categorical models)
## Need to strip last column e.g. confusion[1:9,1:10] in rf object would be confusion[1:9,1:9]
#confmatrix = soiclass$confusion[1:9,1:9] 
#OOBkappa = Kappa.test(confmatrix, conf.level = 0.95)
#write.table(confmatrix, file = "MLRA35_ESG_oob_RFSoilGrids100m_conf_matr.txt", sep = "\t")


## Reference covar rasters to use in prediction
setwd("O:/Models_active_work/UpCo/final_covars")
#rasters=stack(list.files(getwd(),pattern=".tif$",full.names=FALSE))
DEPTH =  calc(projgrid, fun=function(x)(ifelse(x>-999,0)), progress="text") #raster to set prediction depth
names(DEPTH)<-"DEPTH"
rasters=stack(cov.grids, DEPTH)
#rasters = setMinMax(brick(rasters))
names(rasters)

## Predict onto covariate grid
setwd("O:/Models_active_work/UpCo/Kw")
## Parallelized predict
beginCluster(31,type='SOCK')
pred <- clusterR(rasters, predict, args=list(model=soiclass),progress="text")
endCluster()

###Block for outputs
setwd("O:/Models_active_work/UpCo/Kw/outputs")
writeRaster(pred, overwrite=TRUE,filename="kw_0cm_QRF.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")



################### Pedon-based Cross validation ################################
pts$pedon_key <- paste(pts$long.1,pts$lat.1,sep="")
pedons <- unique(pts$pedon_key)
pedondf <- as.data.frame(pedons)
nfolds <- 10
pedondf$folds <- sample.int(nfolds,size =length(pedons),replace=T)
ptspred.listcv <- c(ptspred.list,"pedon_key")
pts.extcv <- pts[c(ptspred.listcv)]
pts.extcv$DEPTH <- as.numeric(as.character(pts.extcv$DEPTH))
pts.extcv <- na.omit(pts.extcv)# Remove any record with NA's (in any column - be careful)
pts.extcv  <- subset(pts.extcv, pts.extcv$K_sw != "NA")
pts.extcv <- cbind(pedondf[match(pts.extcv$pedon_key, pedondf$pedons),], pts.extcv) #attaches to nasispedons
formulaStringCVp <- as.formula(paste('K_sw ~','DEPTH','+', paste(gsub(".tif","", cov.grids), collapse="+")))
pts.extcv$pcvpred <- "NA"
for (g in seq(nfolds)){
  traindf <- subset(pts.extcv, pts.extcv$folds != g)
  rf.pcv <- randomForest(formulaStringCVp, data=traindf, importance=FALSE, proximity=FALSE, ntree=100, keep.forest=TRUE)
  pts.extcv$pcvpred <- ifelse(pts.extcv$folds == g, predict(rf.pcv, newdata=pts.extcv),pts.extcv$pcvpred)
  print(g)
}
pts.extcv$pcvpred = as.numeric(pts.extcv$pcvpred)
cvp.RMSE = sqrt(mean((pts.extcv$K_sw - pts.extcv$pcvpred)^2, na.rm=TRUE))
cvp.Rsquared = 1-var(pts.extcv$K_sw - pts.extcv$pcvpred, na.rm=TRUE)/var(pts.extcv$K_sw, na.rm=TRUE)
pts.extcv5 <- subset(pts.extcv, DEPTH<5)
cvp5.RMSE = sqrt(mean((pts.extcv5$K_sw - pts.extcv5$pcvpred)^2, na.rm=TRUE))
cvp5.Rsquared = 1-var(pts.extcv5$K_sw - pts.extcv5$pcvpred, na.rm=TRUE)/var(pts.extcv5$K_sw, na.rm=TRUE)



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

required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall", "quantregForest")# might need snowfall
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
memory.limit(500000)
rasterOptions(maxmemory = 1e+9)
options(scipen = 999)
par(mar=c(0.3,0.3,0.3,0.3))


######## Get points for extraction if in table form ###########
setwd("O:/Models_active_work/UpCo/Avail_H2O_2D_CV")
pts <- read.delim("NCSS17_BD_Moisture_ttab.txt") # If in delimited file other than csv
pts<-pts[,c(2,3,6,12,13,29,30,58)]

### Weed out points with imprecise coordinates ###
pts$latnchar = nchar(abs(pts$latitude_decimal_degrees))
pts$longnchar = nchar(abs(pts$longitude_decimal_degrees))
ptsc = subset(pts, pts$latnchar > 5 & pts$longnchar > 6)

### Turn into spatial file
shp.pts <- ptsc
coordinates(shp.pts) <- ~ longitude_decimal_degrees + latitude_decimal_degrees
temp.proj <- CRS("+proj=longlat +datum=WGS84") ## specify projection
projection(shp.pts) <- temp.proj


######## Load map clip boundary (if needed) ###########
setwd("V:/PROJECTS/TRAVIS_NAUMAN/GIS_Archive/CO_River_watershed_Meade")
polybound <- readOGR(".", "CO_River_watershed_Meade_ACEA")
polybound <- spTransform(polybound, temp.proj)
## Now clip points and check with visualization
shp.pts = shp.pts[polybound,]#clip by outer extent of all polybound features
plot(polybound)
plot(shp.pts, add=TRUE)
shp.pts$depth = (shp.pts$hzn_bot-shp.pts$hzn_top)/2

######### Grid Prep #################
## Make list of grids
setwd("O:/Models_active_work/UCRB_Covariates")
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
shp.pts$DID <- seq.int(nrow(shp.pts))
pts.ext <- merge(as.data.frame(shp.pts),ov.lst, by="DID")

## Save points
setwd("O:/Models_active_work/UpCo/Avail_H2O_2D_CV")
write.table(pts.ext, "cop_ncss17_w15l2_covarsc.txt", sep = "\t", row.names = FALSE)
#pts$swregap = as.factor(pts$swregap)
#pts$LFelems  = as.factor(as.character(pts$LFelems))
#pts$pscsmodalb = as.factor(pts$pscsmodalb)

## Prep for Random Forest
pts.extc <- subset(pts.ext, as.numeric(pts.ext$hzn_top) <= 00 & as.numeric(pts.ext$hzn_bot) > 00) # subset to chosen depth
ptspred.list <- gsub(".tif","", cov.grids)# Take .tif off of the grid list to just get the variable names
pred <- "w15l2" ## Dependent variable
ptspred.list <- c(ptspred.list,pred) #Add dependent variable
pts.extc <- pts.extc[c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use 
pts.extc <- na.omit(pts.extc)# Remove any record with NA's (in any column - be careful)
pts.extc  <- subset(pts.extc, pts.extc$w15l2 != "NA")
xtrain <- as.matrix(pts.extc[c(gsub(".tif","", cov.grids))])
ytrain <- c(as.matrix(pts.extc[c(pred)]))
logytrain <- log(ytrain)
sqrtytrain<-sqrt(ytrain)

summary(pts.extc$w15l2)
############### Build quantile Random Forest
Qsoiclass <- quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=100, keep.forest=TRUE)
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
setwd("O:/Models_active_work/UpCo/BD_2D_CV/outputs")
## Parallelized predict
beginCluster(31,type='SOCK')
pred = clusterR(rasters, predict, args=list(model=soiclass),progress="text")
endCluster()
writeRaster(pred, overwrite=F,filename="db_od_00cm_QRF.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")


##################### Run Cross Validation ######################################
trainx = gsub(".tif","", cov.grids)
trainx = c(trainx)
trainx = pts.extc[c(trainx)]
trainy = pts.extc[c("w15l2")]
set.seed(41)
## rfcv runs a cross val and tests variable importance
nfolds <- 10
cv.rf = rfcv(trainx, trainy$w15l2, cv.fold = nfolds, proximity=FALSE, ntree=100, keep.forest=TRUE)
par(mar=c(5, 4, 4, 2) + 0.1)
with(cv.rf, plot(n.var, error.cv, log="x", type="o", lwd=2))
# Now put the full cross val out of the rfcv
pts.extc$cvpred = cv.rf$predicted$`40` ## the '#' at end corresponds to the number of variables included should = number of vars used
cv.RMSE = sqrt(mean(((pts.extc$w15l2) - pts.extc$cvpred)^2, na.rm=TRUE))
cv.Rsquared = 1-var((pts.extc$w15l2) - pts.extc$cvpred, na.rm=TRUE)/var((pts.extc$w15l2), na.rm=TRUE)

cv.RMSE
cv.Rsquared



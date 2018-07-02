######################
## Random Forest script that includes:
## Extraction of covariates to points
## Confustion matrix creation
## Kappa calculation
######################



# Workspace setup
# Install packages if not already installed
required.packages <- c("ggplot2", "raster", "sp", "rgdal", "plyr", "ncdf4","maptools", "randomForest", "fmsb","snow", "snowfall","parallel")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)

## Load shapefile ##
setwd("/media/tnaum/D/GIS_Archive/AIM2015/2015_Data_TOSHARE")## FOlder with points
shp.pts <-readOGR(".", "AIM15_merge_cop14")
point.proj <- projection(shp.pts)
## If not prj file and you know proj, can specify by name
#shp.proj <- CRS("+proj=longlat +datum=WGS84")

######### Grid Prep #################
## Make list of grids
setwd("G:/GIS_Archive/COP/Landsat/2014bands")
cov.grids <- list.files(pattern=".tif$")
## If points need to be matched up to grids ###
projgrid = raster("b2.tif")
## Or Make a stack of grids to extract all at once (for smaller datasets)
#cov.stack <- stack()
cov.proj <- projection(projgrid)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match aet rasters

## Plot to ensure alignment bw points and rasters
plot(projgrid)
plot(shp.pts, add=TRUE)

#Clean up for weird error in extract
shp.pts@coords <- shp.pts@coords[, 1:2]

## Loop through each raster to extract points (larger datasets)
for (g in cov.grids){
  ext.grid <- raster(g)
  shp.pts <- extract(ext.grid, shp.pts, df=TRUE, sp=TRUE)
  print(g)
}
rm(g,ext.grid)


## Save points
setwd("G:/BLM_Salinity/C_fact_RF")
write.table(as.data.frame(shp.pts), "AIM14_cover_landsat_covars.txt", sep = "\t")
pts = as.data.frame(shp.pts)
pts$satvit = ((pts$b6-ptsc$b2)/(pts$b6-pts$b2+0.9))*(1.9)-(pts$b7/2)

## Prep for Random Forest
ptspred.list = gsub(".tif","", cov.grids)# Take .tif off of the grid list to just get the variable names
pred = "BSOILCOVER" ## Dependent variable
satvi = "satvit"
ptspred.list = c(ptspred.list,pred,satvi) #Add dependent variable
ptsc = pts[c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use 
ptsc = na.omit(ptsc)# Remove any record with NA's (in any column - be careful)

## Build Random Forest
soiclass = randomForest(as.numeric(BSOILCOVER) ~ ., data = ptsc, importance=TRUE, proximity=FALSE, ntree=1000)
soiclass ## Get oob error
## Creat Confusion matrix
## Need to strip last column e.g. confusion[1:9,1:10] in rf object would be confusion[1:9,1:9]

## Reference covar rasters to use in prediction
#setwd("G:/BLM_Salinity/LS8_2016/bands")
#cov.grids <- list.files(pattern=".tif$")
#rasters=stack(list.files(getwd(),pattern=".tif$",full.names=FALSE))
#rasters=stack(cov.grids)

## Predict onto covariate grid
#setwd("G:/BLM_Salinity/C_fact_RF")
#predict(rasters, soiclass, type="response",index=2,na.rm=TRUE,progress="window",overwrite=TRUE,filename="C_fact_COriv16_frCOP14_30m_RF.tif")
## Create lookup table

##################################################################
### Stepwise regression
library(DAAG)
library(MASS)
lfit = lm(BSOILCOVER ~ ., data = ptsc)
summary(lfit)
stepbo = stepAIC(lfit, direction = "both")
summary(stepbo)

## Cross Validation
par(mar=c(5, 4, 4, 2) + 0.1)
stepbocv = CVlm(data=ptsc, form.lm=stepbo, plotit="Observed", m=10)
plot(stepbocv$BSOILCOVER~stepbocv$cvpred)
stepbocv.RMSE = sqrt(mean((stepbocv$`BSOILCOVER` - stepbocv$cvpred)^2, na.rm=TRUE))
stepbocv.Rsq = 1-var(stepbocv$`BSOILCOVER` - stepbocv$cvpred, na.rm=TRUE)/var(stepbocv$`BSOILCOVER`, na.rm=TRUE)

## Reference covar rasters to use in prediction: 2013 Landsat
setwd("/home/tnaum/data/BLM_salinity/Landsat8_2013_07_01_10_30/bands")
cov.grids <- list.files(pattern=".tif$")
#rasters=stack(list.files(getwd(),pattern=".tif$",full.names=FALSE))
rasters=stack(cov.grids)

## Predict Regression
setwd("/home/tnaum/data/BLM_salinity/Landsat8_2013_07_01_10_30/bareground_2013")
beginCluster(31,type='SOCK')
pred = clusterR(rasters, predict, args=list(model=stepbo),progress="text",filename="bg_ls8_2013_07_01_10_30.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT2S')
endCluster()


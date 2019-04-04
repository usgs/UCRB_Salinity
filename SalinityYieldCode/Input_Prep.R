#### Layer Creation for DSM SPARROW


## Load packages
required.packages <- c("raster", "sp", "rgdal","snow", "snowfall","parallel", "itertools","doParallel", "plyr", "ncdf4","maptools", "rgeos")# maybe need dplyr??
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)

##### Base rasters
ec0 <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/ec_r_0_cm_2D_QRF_bt.tif")
ec30 <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/ec_r_30_cm_2D_QRF_bt.tif")
ec60 <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/ec_r_60_cm_2D_QRF_bt.tif")
ec100 <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/ec_r_100_cm_2D_QRF_bt.tif")
## Initializew cluster for computing
setwd("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs")
beginCluster(30,type='SOCK')
# compute ave ec
ecm_fn <- function(ec0, ec30,ec60,ec100) { ind <- ec0 + ec30 + ec60 + ec100
return(ind)
}
ecm_stk <-stack(ec0, ec30, ec60, ec100)
ecm <- clusterR(ecm_stk, overlay, args=list(fun=ecm_fn), progress='text', filename="ecave_mask.tif")
#ecm <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ecave_mask.tif") # to load if restarting
mask <- raster("/home/tnaum/data/BLM_Salinity/water_mask/NLCD_water_mask.tif")

##### UCRB extent
ucrb_bnd <- readOGR("/home/tnaum/Dropbox/USGS/BLM_projects/Utah_BLM_Salinity/Huc6_boundary", "CO_River_watershed_Meade", stringsAsFactors = F)

##### Agricultural lands Sources
agland <- readOGR("/home/tnaum/data/BLM_Salinity/SPARROW/data/UCRB_Ag_Buto/sir2014_5039_UCRBAgriculture_v2.gdb", stringsAsFactors = F)
agland$IrrigClass <- paste(agland$UCRB_irrigation_status,agland$UCRB_irrigation_method,sep="-")
agIrrLand <- agland[agland$UCRB_irrigation_status=="Irrigated",]
agIrrLand$IrrigClassID <- ifelse(agIrrLand$IrrigClass == "Irrigated-Flood", 1, 2)## 1s are flooded, 2s are other types of irrigation
agIrrLandrast <- rasterize(agIrrLand, ecm,  field=agIrrLand$IrrigClassID, progress="text", datatype='INT1U', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/intermed_layers/IrrigAgrast.tif") # slow but steady
agIrrLandrast <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/intermed_layers/IrrigAgrast.tif") # To open in future runs

###### Land Source files: EC and agriculture combinations
# Highest (>90th quantile) EC non-irrigated
ec90_100s_fn <- function(ecm, agIrrLandrast) { ind <- ifelse(ecm>unname(quantile(ecm, probs=0.90,na.rm=T))&agIrrLandrast==0,1,0) #EC50q ~ 0.35 dS/m
return(ind)
}
ec90_100s_stk <-stack(ecm, agIrrLandrast)
ec90_100s <- clusterR(ec90_100s_stk, overlay, args=list(fun=ec90_100s_fn), progress='text', datatype='INT1U', filename="ec90_100s.tif")
# 75-90th quantile EC non-irrigated
ec75_90s_fn <- function(ecm, agIrrLandrast,ec90_100s) { ind <- ifelse(ecm>unname(quantile(ecm, probs=0.75,na.rm=T))&agIrrLandrast==0&ec90_100s==0,1,0) #EC50q ~ 0.35 dS/m
return(ind)
}
ec75_90s_stk <-stack(ecm, agIrrLandrast,ec90_100s)
ec75_90s <- clusterR(ec75_90s_stk, overlay, args=list(fun=ec75_90s_fn), progress='text', datatype='INT1U', filename="ec75_90s.tif")
# 50-75th quantile EC non-irrigated
ec50_75s_fn <- function(ecm, agIrrLandrast,ec90_100s,ec75_90s) { ind <- ifelse(ecm>unname(quantile(ecm, probs=0.50,na.rm=T))&agIrrLandrast==0&ec90_100s==0&ec75_90s==0,1,0) #EC50q ~ 0.35 dS/m
return(ind)
}
ec50_75s_stk <-stack(ecm, agIrrLandrast,ec90_100s,ec75_90s)
ec50_75s <- clusterR(ec50_75s_stk, overlay, args=list(fun=ec50_75s_fn), progress='text', datatype='INT1U', filename="ec50_75s.tif")
# <10th quantile EC non-irrigated
ec0_10s_fn <- function(ecm, agIrrLandrast) { ind <- ifelse(ecm<unname(quantile(ecm, probs=0.10,na.rm=T))&agIrrLandrast==0,1,0) 
return(ind)
}
ec0_10s_stk <-stack(ecm, agIrrLandrast)
ec0_10s <- clusterR(ec0_10s_stk, overlay, args=list(fun=ec0_10s_fn), progress='text', datatype='INT1U', filename="ec0_10s.tif")
# 10-25th quantile EC non-irrigated
ec10_25s_fn <- function(ecm, agIrrLandrast, ec0_10s) { ind <- ifelse(ecm<unname(quantile(ecm, probs=0.25,na.rm=T))&agIrrLandrast==0&ec0_10s==0,1,0) 
return(ind)
}
ec10_25s_stk <-stack(ecm, agIrrLandrast,ec0_10s)
ec10_25s <- clusterR(ec10_25s_stk, overlay, args=list(fun=ec10_25s_fn), progress='text', datatype='INT1U', filename="ec10_25s.tif")
# 25-50th quantile EC non-irrigated
ec25_50s_fn <- function(ecm, agIrrLandrast, ec0_10s, ec10_25s) { ind <- ifelse(ecm<unname(quantile(ecm, probs=0.50,na.rm=T))&agIrrLandrast==0&ec0_10s==0&ec10_25s==0,1,0) 
return(ind)
}
ec25_50s_stk <-stack(ecm, agIrrLandrast,ec0_10s, ec10_25s)
ec25_50s <- clusterR(ec25_50s_stk, overlay, args=list(fun=ec25_50s_fn), progress='text', datatype='INT1U', filename="ec25_50s.tif")

# Lower EC flood irrigated ag
ec0_75sF_fn <- function(ecm, agIrrLandrast) { ind <- ifelse(ecm<unname(quantile(ecm, probs=0.75,na.rm=T))&agIrrLandrast==1,1,0) #EC50q ~ 0.35 dS/m
return(ind)
}
ec0_75sF_stk <-stack(ecm, agIrrLandrast)
ec0_75sF <- clusterR(ec0_75sF_stk, overlay, args=list(fun=ec0_75sF_fn),progress='text', datatype='INT1U', filename="ec0_75sF.tif")
# High EC flood irrigated ag
ec75_100sF_fn <- function(ecm, agIrrLandrast) { ind <- ifelse(ecm>unname(quantile(ecm, probs=0.75,na.rm=T))&agIrrLandrast==1,1,0) #EC50q ~ 0.35 dS/m
return(ind)
}
ec75_100sF_stk <-stack(ecm, agIrrLandrast)
ec75_100sF <- clusterR(ec75_100sF_stk, overlay, args=list(fun=ec75_100sF_fn),progress='text', datatype='INT1U', filename="ec75_100sF.tif")
# Low EC Non-flood irrigated ag
ec0_75sN_fn <- function(ecm, agIrrLandrast) { ind <- ifelse(ecm<unname(quantile(ecm, probs=0.75,na.rm=T))&agIrrLandrast==2,1,0) #EC50q ~ 0.35 dS/m
return(ind)
}
ec0_75sN_stk <-stack(ecm, agIrrLandrast)
ec0_75sN <- clusterR(ec0_75sN_stk, overlay, args=list(fun=ec0_75sN_fn),progress='text', datatype='INT1U', filename="ec0_75sN.tif")
# High EC Non-flood irrigated ag
ec75_100sN_fn <- function(ecm, agIrrLandrast) { ind <- ifelse(ecm>unname(quantile(ecm, probs=0.75,na.rm=T))&agIrrLandrast==2,1,0) #EC50q ~ 0.35 dS/m
return(ind)
}
ec75_100sN_stk <-stack(ecm, agIrrLandrast)
ec75_100sN <- clusterR(ec75_100sN_stk, overlay, args=list(fun=ec75_100sN_fn),progress='text', datatype='INT1U', filename="ec75_100sN.tif")

endCluster()
gc()

##### Point Sources: turn to raster for summation in guage summary function
pts <- read.delim("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/UCRB_Spring_Locations_edited.txt") # If in delimited file other than csv
# Turn into spatial file
shp.pts <- pts
coordinates(shp.pts) <- ~ Primary_lon_dec + Primary_lat_dec
temp.proj <- CRS("+proj=longlat +datum=WGS84") ## specify projection
projection(shp.pts) <- temp.proj
writeOGR(shp.pts, ".", "/home/tnaum/data/BLM_Salinity/DSM_SPARROW/UCRB_springs_wgs", driver="ESRI Shapefile")
cov.proj <- projection(ecm)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match rasters
# Now rasterize
pt.sourcerast <- rasterize(shp.pts, ecm,  field=shp.pts$Tons_2010, progress="text", datatype='INT4U', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/pt_sourcerast.tif") # slow but steady


##### New BCM layer preparation
## Soil Grids bedrock within 2m probability
Rprob <- raster("/media/tnaum/Seagate Backup Plus Drive/GIS_Archive/SoilGrids250m/props_extract/BDRLOG_M_250m_ll.tif")
e <- extent(ucrb_bnd)
Rprobc <- crop(Rprob,e, progress="text")
rasterOptions(maxmemory = 1e+07, chunksize = 1e+05)
beginCluster(30,type='SOCK')
Rprobcc <- projectRaster(Rprobc,mask, progress='text') ## Still experimenting with this in parallel, it eats memory like a hog
rm(Rprobc)
endCluster()
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)
beginCluster(30,type='SOCK')
f_mask <- function(a,b) a*b
Rprobm_stk <- stack(Rprobcc,mask)
clusterR(Rprobm_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/SG250_Rprob.tif")
rm(Rprobcc)
gc()
## 30m layers: NCSS soil preds, DART
# Soil Preds: 0 cm depth predictions
sar <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/sar_r_0_cm_2D_QRF_bt.tif")
sar_stk <- stack(sar,mask)
clusterR(sar_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/sar_m.tif")
rock <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/fragvol_r_0_cm_2D_QRF_bt.tif")
rock_stk <- stack(rock,mask)
clusterR(rock_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/rock_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
fs <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/sandfine_r_0_cm_2D_QRF.tif")
fs_stk <- stack(fs,mask)
clusterR(fs_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/fs_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
ph <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/ph1to1h2o_r_0_cm_2D_QRF.tif")
ph_stk <- stack(ph,mask)
clusterR(ph_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ph_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
awc <- raster("/home/tnaum/data/BLMsoils/All_Prop_Maps/awc_r_0_cm_2D_QRF.tif")
awc_stk <- stack(awc,mask)
clusterR(awc_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/awc_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
# DART layers
elevm <- raster("/home/tnaum/Dropbox/CP DSM/UCRB_Covariates/ELEVm.tif")
elevm_stk <- stack(elevm,mask)
clusterR(elevm_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/elevm_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
ppt_ratio <- raster("/home/tnaum/Dropbox/CP DSM/UCRB_Covariates/ppt_ratio.tif")
ppt_ratio_stk <- stack(ppt_ratio,mask)
clusterR(ppt_ratio_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ppt_ratio_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
prot_index <- raster("/home/tnaum/Dropbox/CP DSM/UCRB_Covariates/PROTINDEX.tif")
prot_index_stk <- stack(prot_index,mask)
clusterR(prot_index_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/prot_index_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
slope <- raster("/home/tnaum/Dropbox/CP DSM/UCRB_Covariates/SLOPE.tif")
slope_stk <- stack(slope,mask)
clusterR(slope_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/slope_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
sness <- raster("/home/tnaum/Dropbox/CP DSM/UCRB_Covariates/SOUTHNESS.tif")
sness_stk <- stack(sness,mask)
clusterR(sness_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/sness_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
ppt <- raster("/home/tnaum/Dropbox/CP DSM/UCRB_Covariates/ppt_ann.tif")
ppt_stk <- stack(ppt,mask)
clusterR(ppt_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ppt_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))


endCluster()

## UCRB BCM variables from Miller et al., 2017
rasterOptions(maxmemory = 1e+07, chunksize = 1e+06)
beginCluster(30,type='SOCK')
excess <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/UCRB_BCM_SPARROW17/UCRB_BCM/exc1985to2012wy_meanTotalAnnual_x100i.tif")
excessc <- projectRaster(excess,mask, progress='text') ## Still experimenting with this in parallel, it eats memory like a hog
excessc_stk <- stack(excessc, mask)
clusterR(excessc_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/excess_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
deficit <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/UCRB_BCM_SPARROW17/UCRB_BCM/cwd1985to2012wy_meanTotalAnnual_x100i.tif")
deficitc <- projectRaster(deficit,mask, progress='text') ## Still experimenting with this in parallel, it eats memory like a hog
deficitc_stk <- stack(deficitc, mask)
clusterR(deficitc_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/deficit_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
melt <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/UCRB_BCM_SPARROW17/UCRB_BCM/mlt1985to2012wy_meanTotalAnnual_x100i.tif")
meltc <- projectRaster(melt,mask, progress='text') ## Still experimenting with this in parallel, it eats memory like a hog
meltc_stk <- stack(meltc, mask)
clusterR(meltc_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/melt_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))
rech <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/UCRB_BCM_SPARROW17/UCRB_BCM/rch1985to2012wy_meanTotalAnnual_x100i.tif")
rechc <- projectRaster(rech,mask, progress='text') ## Still experimenting with this in parallel, it eats memory like a hog
rechc_stk <- stack(rechc, mask)
clusterR(rechc_stk, overlay, args=list(fun=f_mask),progress='text', datatype='FLT4S', filename="/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/rech_m.tif",options=c("COMPRESS=DEFLATE", "TFW=YES"))

endCluster()

###############################################
#### Erosion risk index creation using various rasters and using quantile breaks to 
#### identify areas at high risk of contributing salts to 
#### washes in the UCRB

## Load packages
required.packages <- c("raster", "sp", "rgdal","snow", "snowfall","parallel", "itertools","doParallel")# might need snowfall
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)

## Rasters to use
bg <- raster("/home/tnaum/data/BLM_Salinity/Landsat8_2013_07_01_10_30/bareground_2013/bg_ls8_2013_07_01_10_30alb.tif")
mask <- raster("/home/tnaum/data/BLM_Salinity/water_mask/NLCD_water_mask.tif")
gap <- raster("/home/tnaum/data/UCRB_Covariates/GAP.tif")
ec <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ecave_mask.tif")
kwc <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/kw_m.tif")
flen <- raster("/home/tnaum/data/BLM_Salinity/risk_index/flength_mask.tif")
facc<- raster("/home/tnaum/data/BLM_Salinity/risk_index/CAlog10_mask.tif")
## Gap lookup table
gap_lup <- read.csv("/home/tnaum/Dropbox/USGS/BLM_projects/Utah_BLM_Salinity/GAP_analysis/Class_lookup_UCRB.csv")

## Set working Drive
setwd("/home/tnaum/data/BLM_Salinity/risk_index")

## Load up cluster for processing
cpus <- detectCores(logical=TRUE)-2
beginCluster(cpus,type='SOCK')

## Make ec and kw quantile base index dummy rasters
# Areas with EC >75 %tile
ec75q_fn <- function(ec) { ind <- ifelse(ec>unname(quantile(ec, probs=0.75,na.rm=T)),1,0)
return(ind)
}
ec75q <- clusterR(ec, calc, args=list(fun=ec75q_fn),progress='text', datatype='INT1U', filename="srisk_ec75q.tif")
# Areas with EC >50 %tile
ec50q_fn <- function(ec) { ind <- ifelse(ec>unname(quantile(ec, probs=0.50,na.rm=T)),1,0)
return(ind)
}
ec50q <- clusterR(ec, calc, args=list(fun=ec50q_fn),progress='text', datatype='INT1U', filename="srisk_ec50q.tif")
# Areas with EC >90 %tile
ec90q_fn <- function(ec) { ind <- ifelse(ec>unname(quantile(ec, probs=0.90,na.rm=T)),1,0)
return(ind)
}
ec90q <- clusterR(ec, calc, args=list(fun=ec90q_fn),progress='text', datatype='INT1U', filename="srisk_ec90q.tif")
# Areas with kw >75 %tile
kw75q_fn <- function(kwc) { ind <- ifelse(kwc>unname(quantile(kwc, probs=0.75,na.rm=T)),1,0)
return(ind)
}
kw75q <- clusterR(kwc, calc, args=list(fun=kw75q_fn),progress='text', datatype='INT1U', filename="srisk_kw75q.tif")

## GAP reclass to Macrogroup
gap_macro_mat <- as.matrix(gap_lup[c("Value","macroval")], dimnames = NULL)
colnames(gap_macro_mat) <- NULL
#gap_macro <- clusterR(gap, reclassify, args=list(rcl=gap_macro_mat, right=FALSE), progress='text')
#writeRaster(gap_macro, overwrite=F,filename="gap_macrogroup.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
gap_macro <- raster("gap_macrogroup.tif")

## Set up mask function to get rid of boundary and water pixels
f_mask <- function(a,b) a*b # a is the mask raster (1 for data pixels, NAs for others), and b is the raster to be masked

## Bare ground quantiles
## first reproject
gc()
rasterOptions(maxmemory = 1e+08,chunksize = 1e+06)
bgc <- projectRaster(bg,mask, progress='text') 
rm(bg)
gc()

rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)
bg_stack <- stack(mask,bgc)
bgm <- clusterR(bg_stack, overlay, args=list(fun=f_mask), progress='text')
# fix abnormally high and low values
f_lowcor <- function(a) {a[a<0]<-0
return(a)
}
f_highcor <- function(a) {a[a>100]<-100
return(a)
}
bg_stack <- stack(bgm)
bgm <- clusterR(bg_stack, calc, args=list(fun=f_lowcor), progress='text')
bg_stack <- stack(bgm)
bgm <- clusterR(bg_stack, calc, args=list(fun=f_highcor), progress='text')
writeRaster(bgm, overwrite=F,filename="bareground_mask.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
rm(bgc,bg_stack)
gc()
bgm <- raster("/home/tnaum/data/BLM_Salinity/risk_index/bareground_mask.tif")
## Create empty grid for use in bare ground GAP macrogroup loop
f_empty <- function(a) {a[a==1]<-0
return(a)
}
mask_stk <- stack(mask)
empty <- clusterR(mask_stk, calc, args=list(fun=f_empty), progress='text',filename="empty.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc()
empty <- raster("/home/tnaum/data/BLM_Salinity/risk_index/empty.tif")
##### Reclassify cover within GAP macrogroups by 90th percentiles
gaplist <- unique(gap_lup$macroval)
#fn_class90quan <- function(gapmacroval) {bgmask[bgm>quantile(bgm[gap_macro==gapmacroval], probs=0.9,na.rm=T)]<-1}
bcls_fn <- function(b,c) {bcls <- ifelse(c==i,b,NA)
return(bcls)
}
fn_class90quan <- function(a,bcls) {
  t <- unname(quantile(bcls, probs=0.9,na.rm=T))
  a[bcls>t]<-1
  return(a)
}
a <- raster("empty.tif")
b <- raster("bareground_mask.tif")
c <- raster("gap_macrogroup.tif")
for(i in gaplist){
  b_stk <- stack(b,c)
  bcls <- clusterR(b_stk, overlay, args=list(fun=bcls_fn), export='i', progress='text')
  a_stk <- stack(a,bcls)
  a <- clusterR(a_stk, overlay, args=list(fun=fn_class90quan), progress='text')
  rm(bcls)
  gc()
  print(i)
}
writeRaster(a, overwrite=F,filename="bgmacro90_done.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")

##### Reclassify cover within GAP macrogroups by 75th percentiles
gaplist <- unique(gap_lup$macroval)
#fn_class90quan <- function(gapmacroval) {bgmask[bgm>quantile(bgm[gap_macro==gapmacroval], probs=0.9,na.rm=T)]<-1}
bcls_fn <- function(b,c) {bcls <- ifelse(c==i,b,NA)
return(bcls)
}
fn_class75quan <- function(a,bcls) {
  t <- unname(quantile(bcls, probs=0.75,na.rm=T))
  a[bcls>t]<-1
  return(a)
}
a <- raster("empty.tif")
b <- raster("bareground_mask.tif")
c <- raster("gap_macrogroup.tif")
for(i in gaplist){
  b_stk <- stack(b,c)
  bcls <- clusterR(b_stk, overlay, args=list(fun=bcls_fn), export='i', progress='text')
  a_stk <- stack(a,bcls)
  a <- clusterR(a_stk, overlay, args=list(fun=fn_class75quan), progress='text')
  rm(bcls)
  gc()
  print(i)
}
writeRaster(a, overwrite=F,filename="bgmacro75.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
rm(a,b,c)
gc()

##### Now complex indices
## Some new rasters to include
bgmacro90d <- raster("/home/tnaum/data/BLM_Salinity/risk_index/bgmacro90_done.tif")
bgm <- raster("/home/tnaum/data/BLM_Salinity/risk_index/bareground_mask.tif")
bgmacro75d <- raster("/home/tnaum/data/BLM_Salinity/risk_index/bgmacro75.tif")
## Most conservative index
srisk_b90m40p_facc90q_ec90q_kw90q_f500m_fn <- function(bgm,bgmacro90d,ec90q,flen,facc,kwc) {
  ind <- ifelse((bgm>40|bgmacro90d==1)&ec90q==1&flen<500&facc>unname(quantile(facc, probs=0.9,na.rm=T))&kwc>unname(quantile(kwc, probs=0.9,na.rm=T)),1,NA)
  return(ind) # q90 Ksw is 0.0355, q90 of facc is 4804
}
b90m40p_facc90q_ec90q_kw90q_f500m_stk <-stack(bgm,bgmacro90d,ec90q,flen,facc,kwc)
srisk_b90m40p_facc90q_ec90q_kw90q_f500m <- clusterR(b90m40p_facc90q_ec90q_kw90q_f500m_stk, overlay, args=list(fun=srisk_b90m40p_facc90q_ec90q_kw90q_f500m_fn),progress='text',filename="srisk_ec90q_b90m40p_kw90q_facc90q_f500m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## Bareground 75q or 30%
srisk_b75m30p_fn <- function(bgm,bgmacro75d) {
  ind <- ifelse((bgm>30|bgmacro75d==1),1,NA)
  return(ind) # q90 Ksw is 0.0355, q90 of facc is 4804
}
b75m30p_stk <-stack(bgm,bgmacro75d)
srisk_b75m30p <- clusterR(b75m30p_stk, overlay, args=list(fun=srisk_b75m30p_fn),progress='text',filename="srisk_bgm75q30p.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## EC75q Bareground 75q or 30%
srisk_ec75q_b75m30p_fn <- function(bgm,bgmacro75d,ec75q) {
  ind <- ifelse((bgm>30|bgmacro75d==1)&ec75q==1,1,NA)
  return(ind) # q90 Ksw is 0.0355, q90 of facc is 4804
}
ec75q_b75m30p_stk <-stack(bgm,bgmacro75d,ec75q)
srisk_ec75q_b75m30p <- clusterR(ec75q_b75m30p_stk, overlay, args=list(fun=srisk_ec75q_b75m30p_fn),progress='text',filename="srisk_ec75q_b75m30p.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## EC75q Bareground 75q
srisk_ec75q_bgm75q_fn <- function(bgmacro75d,ec75q) {
  ind <- ifelse(bgmacro75d==1&ec75q==1,1,NA)
  return(ind) # q90 Ksw is 0.0355, q90 of facc is 4804
}
ec75q_bgm75q_stk <-stack(bgmacro75d,ec75q)
srisk_ec75q_bgm75q <- clusterR(ec75q_bgm75q_stk, overlay, args=list(fun=srisk_ec75q_bgm75q_fn),progress='text',filename="srisk_ec75q_bgm75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## EC75q Bareground 75q, kw > 75th quan
srisk_ec75q_bgm75q_kw75q_fn <- function(bgmacro75d,ec75q,kw75q) {
  ind <- ifelse(bgmacro75d==1&ec75q==1&kw75q==1,1,NA)
  return(ind) 
}
ec75q_bgm75q_kw75q_stk <-stack(bgmacro75d,ec75q,kw75q)
srisk_ec75q_bgm75q_kw75q <- clusterR(ec75q_bgm75q_kw75q_stk, overlay, args=list(fun=srisk_ec75q_bgm75q_kw75q_fn),progress='text',filename="srisk_ec75q_bgm75q_kw75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## EC75q Bareground 75q, facc > 75th quan
srisk_ec75q_bgm75q_facc75q_fn <- function(bgmacro75d,ec75q,facc) {
  ind <- ifelse(bgmacro75d==1&ec75q==1&facc>unname(quantile(facc, probs=0.75,na.rm=T)),1,NA)
  return(ind) 
}
ec75q_bgm75q_facc75q_stk <-stack(bgmacro75d,ec75q,facc)
srisk_ec75q_bgm75q_facc75q <- clusterR(ec75q_bgm75q_facc75q_stk, overlay, args=list(fun=srisk_ec75q_bgm75q_facc75q_fn),progress='text',filename="srisk_ec75q_bgm75q_facc75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## EC75q Bareground 75q, flen < 500m
srisk_ec75q_bgm75q_f500m_fn <- function(bgmacro75d,ec75q,flen) {
  ind <- ifelse(bgmacro75d==1&ec75q==1&flen<500,1,NA)
  return(ind) 
}
ec75q_bgm75q_f500m_stk <-stack(bgmacro75d,ec75q,flen)
srisk_ec75q_bgm75q_f500m <- clusterR(ec75q_bgm75q_f500m_stk, overlay, args=list(fun=srisk_ec75q_bgm75q_f500m_fn),progress='text',filename="srisk_ec75q_bgm75q_f500m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## Moderate complete index
srisk_ec75q_bgm75q_kw75q_facc75q_f500m_fn <- function(bgm,bgmacro75d,ec75q,kw75q,facc,flen) {
  ind <- ifelse(bgmacro75d==1&ec75q==1&kw75q==1&flen<500&facc>unname(quantile(facc, probs=0.75,na.rm=T)),1,NA)
  return(ind) 
}
ec75q_bgm75q_kw75q_facc75q_f500m_stk <-stack(bgm,bgmacro75d,ec75q,kw75q,facc,flen)
srisk_ec75q_bgm75q_kw75q_facc75q_f500m<- clusterR(ec75q_bgm75q_kw75q_facc75q_f500m_stk, overlay, args=list(fun=srisk_ec75q_bgm75q_kw75q_facc75q_f500m_fn),progress='text',filename="srisk_ec75q_bgm75q_kw75q_facc75q_f500m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## Moderate complete index with more liberal bare ground inclusion (30%)
srisk_ec75q_b75m30p_kw75q_facc75q_f500m_fn <- function(bgm,bgmacro75d,ec75q,kw75q,facc,flen) {
  ind <- ifelse((bgm>30|bgmacro75d==1)&ec75q==1&kw75q==1&flen<500&facc>unname(quantile(facc, probs=0.75,na.rm=T)),1,NA)
  return(ind) 
}
ec75q_b75m30p_kw75q_facc75q_f500m_stk <-stack(bgm,bgmacro75d,ec75q,kw75q,facc,flen)
srisk_ec75q_b75m30p_kw75q_facc75q_f500m <- clusterR(ec75q_b75m30p_kw75q_facc75q_f500m_stk, overlay, args=list(fun=srisk_ec75q_b75m30p_kw75q_facc75q_f500m_fn),progress='text',filename="srisk_ec75q_b75m30p_kw75q_facc75q_f500m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## Liberal complete index
srisk_ec50q_bgm75q30p_kw50q_facc50q_f1000m_fn <- function(bgm,bgmacro75d,ec50q,kwc,facc,flen) {
  ind <- ifelse((bgm>30|bgmacro75d==1)&ec50q==1&kwc>unname(quantile(kwc, probs=0.50,na.rm=T))&flen<1000&facc>unname(quantile(facc, probs=0.50,na.rm=T)),1,NA)
  return(ind) 
}
ec50q_bgm75q30p_kw50q_facc50q_f1000m_stk <-stack(bgm,bgmacro75d,ec50q,kwc,facc,flen)
srisk_ec50q_bgm75q30p_kw50q_facc50q_f1000m <- clusterR(ec50q_bgm75q30p_kw50q_facc50q_f1000m_stk, overlay, args=list(fun=srisk_ec50q_bgm75q30p_kw50q_facc50q_f1000m_fn),progress='text',filename="srisk_ec50q_bgm75q30p_kw50q_facc50q_f1000m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## Moderate index without hydro terms
srisk_ec75q_b75m30p_kw75q_fn <- function(bgm,bgmacro75d,ec75q,kw75q) {
  ind <- ifelse((bgm>30|bgmacro75d==1)&ec75q==1&kw75q==1,1,NA)
  return(ind) 
}
ec75q_b75m30p_kw75q_stk <-stack(bgm,bgmacro75d,ec75q,kw75q)
srisk_ec75q_b75m30p_kw75q <- clusterR(ec75q_b75m30p_kw75q_stk, overlay, args=list(fun=srisk_ec75q_b75m30p_kw75q_fn),progress='text',filename="srisk_ec75q_b75m30p_kw75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## Moderate index with flow accum >75th perc and no erodibility
srisk_ec75q_bgm75q30p_facc75q_fn <- function(bgm,bgmacro75d,ec75q,facc) {
  ind <- ifelse((bgm>30|bgmacro75d==1)&ec75q==1&facc>unname(quantile(facc, probs=0.75,na.rm=T)),1,NA)
  return(ind) 
}
ec75q_bgm75q30p_facc75q_stk <-stack(bgm,bgmacro75d,ec75q,facc)
srisk_ec75q_bgm75q30p_facc75q <- clusterR(ec75q_bgm75q30p_facc75q_stk, overlay, args=list(fun=srisk_ec75q_bgm75q30p_facc75q_fn),progress='text',filename="srisk_ec75q_bgm75q30p_facc75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## Moderate index with flow length < 500m and no erodibility
srisk_ec75q_bgm75q30p_f500m_fn <- function(bgm,bgmacro75d,ec75q,flen) {
  ind <- ifelse((bgm>30|bgmacro75d==1)&ec75q==1&flen<500,1,NA)
  return(ind) 
}
ec75q_bgm75q30p_f500m_stk <-stack(bgm,bgmacro75d,ec75q,flen)
srisk_ec75q_bgm75q30p_f500m <- clusterR(ec75q_bgm75q30p_f500m_stk, overlay, args=list(fun=srisk_ec75q_bgm75q30p_f500m_fn),progress='text',filename="srisk_ec75q_bgm75q30p_f500m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

## EC75q  kw > 75th quan
srisk_ec75q_kw75q_fn <- function(ec75q,kw75q) {
  ind <- ifelse(ec75q==1&kw75q==1,1,NA)
  return(ind) 
}
srisk_ec75q_kw75q_stk <-stack(ec75q,kw75q)
srisk_ec75q_kw75q <- clusterR(srisk_ec75q_kw75q_stk, overlay, args=list(fun=srisk_ec75q_kw75q_fn),progress='text',filename="srisk_ec75q_kw75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
## EC75q  facc > 75th quan
srisk_ec75q_facc75q_fn <- function(ec75q,facc) {
  ind <- ifelse(ec75q==1&facc>unname(quantile(facc, probs=0.75,na.rm=T)),1,NA)
  return(ind) 
}
srisk_ec75q_facc75q_stk <-stack(ec75q,facc)
srisk_ec75q_facc75q <- clusterR(srisk_ec75q_facc75q_stk, overlay, args=list(fun=srisk_ec75q_facc75q_fn),progress='text',filename="srisk_ec75q_facc75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
## EC75q,  flen > 75th quan
srisk_ec75q_f500m_fn <- function(ec75q,flen) {
  ind <- ifelse(ec75q==1&flen>unname(quantile(flen, probs=0.75,na.rm=T)),1,NA)
  return(ind) 
}
srisk_ec75q_f500m_stk <-stack(ec75q,flen)
srisk_ec75q_f500m <- clusterR(srisk_ec75q_facc75q_stk, overlay, args=list(fun=srisk_ec75q_facc75q_fn),progress='text',filename="srisk_ec75q_f500m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
## EC75q,  kw > 75th quan
srisk_ec75q_kw75q_fn <- function(ec,kwc) {
  ind <- ifelse(ec>unname(quantile(ec, probs=0.75,na.rm=T))&kwc>unname(quantile(kwc, probs=0.75,na.rm=T)),1,NA)
  return(ind) 
}
srisk_ec75q_kw75q_stk <-stack(ec,kwc)
srisk_ec75q_kw75q <- clusterR(srisk_ec75q_kw75q_stk, overlay, args=list(fun=srisk_ec75q_kw75q_fn),progress='text',filename="srisk_ec75q_kw75q.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc() # Flush out RAM

endCluster()





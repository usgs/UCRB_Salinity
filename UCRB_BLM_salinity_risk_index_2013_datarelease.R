###############################################
#### Script taking various rasters and using quantile breaks to 
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
bg <- raster("/home/tnaum/data/BLM_salinity/Landsat8_2013_07_01_10_30/bareground_2013/bg_ls8_2013_07_01_10_30alb.tif")
mask <- raster("/media/tnaum/D/BLM_Salinity/water_mask/NLCD_water_mask.tif")
flen <- raster("/media/tnaum/D/BLM_Salinity/FlowLength/FPL_MFD_seedsplus_1km2.tif")
gap <- raster("/media/tnaum/D/GIS_Archive/UCRB_Covariates/GAP.tif") # USGS GAP dataset
ec50q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec50q.tif")
ec75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q.tif")
ec <- raster("/home/tnaum/data/BLM_salinity/risk_index/ecave_mask.tif")
kwc <- raster("/home/tnaum/data/BLM_salinity/risk_index/kw_mask.tif")
ec90q <- raster("/home/tnaum/data/BLM_salinity/risk_index/ec90high.tif")
kw75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_kw75q.tif")
flen <- raster("/home/tnaum/data/BLM_salinity/risk_index/flength_mask.tif")
facc<- raster("/home/tnaum/data/BLM_salinity/risk_index/CAlog10_mask.tif")
## Gap lookup table
gap_lup <- read.csv("/home/tnaum/Dropbox/USGS/Utah_BLM_Salinity/GAP_analysis/Class_lookup_UCRB.csv")

## Set working Drive
setwd("/home/tnaum/data/BLM_salinity/risk_index")

## Load up cluster for processing
#cpus = detectCores(logical=TRUE)-1
beginCluster(30,type='SOCK')

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

## Create empty grid for use in bare ground GAP macrogroup loop
f_empty <- function(a) {a[a==1]<-0
return(a)
}
mask_stk <- stack(mask)
empty <- clusterR(mask_stk, calc, args=list(fun=f_empty), progress='text',filename="empty.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
gc()

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
bgmacro90d <- raster("/home/tnaum/data/BLM_salinity/risk_index/bgmacro90_done.tif")
bgm <- raster("/home/tnaum/data/BLM_salinity/risk_index/bareground_mask.tif")
bgmacro75d <- raster("/home/tnaum/data/BLM_salinity/risk_index/bgmacro75.tif")
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

endCluster()




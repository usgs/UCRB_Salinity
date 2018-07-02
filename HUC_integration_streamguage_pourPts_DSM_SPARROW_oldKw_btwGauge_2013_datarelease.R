###################################################################
### Script to
### 1) integrate upstreme HUCs to stream guage pour points
### 2) Summarize raster layers by those upstream HUCs
### 3) Correlate the raster summaries to stream guage metrics
###################################################################

## Load packages
required.packages <- c("raster", "sp", "rgdal","snow", "snowfall","parallel", "itertools","doParallel", "plyr", "ncdf4","maptools", "rgeos")# maybe need dplyr??
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)

#### Read in data sources
guages_df <- read.delim("/home/tnaum/data/BLM_salinity/UCRB_streamguage_ec_ttab.txt", stringsAsFactors=F)
hucs <- readOGR("/home/tnaum/data/BLM_salinity/HydroGDB/SIR20175009_UCRB_HydroNetwork.gdb/SIR20175009_UCRB_HydroNetwork.gdb", "sir20175009_UCRB_SPARROW_catchment_p")
strm.reaches <- readOGR("/home/tnaum/data/BLM_salinity/HydroGDB/SIR20175009_UCRB_HydroNetwork.gdb/SIR20175009_UCRB_HydroNetwork.gdb", "sir20175009_UCRB_SPARROW_network_l")
hucs$STAID <- as.numeric(as.character(hucs$STAID))
guages <- as.numeric(guages_df$guageid)
calreach.path <- "/home/tnaum/data/BLM_salinity/DSM_SPARROW/calibration_reaches/"

#### Function/Loop to create shps of incremental (between gauge) calibration reaches of each guage
huc_fn <- function(g,hucs,calreach.path){
  #for(g in guages){
  huclist <- c()
  guagelist <- c()
  starthuc <- hucs[hucs$STAID %in% g,]
  huclist <- c(huclist,starthuc$WATERID)
  tnode <- starthuc$FNODE
  while(length(tnode)>0){ # Loop to work upstream until an endpoint
    newhucs <- hucs[hucs$TNODE %in% tnode,]
    newguages <- ifelse(newhucs$STAID > 0, newhucs$STAID)
    newguages <- newguages[!is.na(newguages)] # check to see if a gauge has been encountered
    if (length(newguages)>0){
      newhucs <- subset(newhucs, !(STAID %in% newguages))
      guagelist <- c(guagelist, newguages)}
    huclist <- c(huclist, newhucs$WATERID)
    tnode <- newhucs$FNODE
  }
  guagehucs <- hucs[hucs$WATERID %in% huclist,]
  guagehucs$dissid <- 1
  guagehuc <- gUnaryUnion(guagehucs, id = guagehucs$dissid)
  filenmpoly <- paste(calreach.path,g,"hucpoly.rds",sep="")
  filenmhuclist <- paste(calreach.path,g,"huclist.rds",sep="")
  filenmguagelist <- paste(calreach.path,g,"guagelist.rds",sep="")
  saveRDS(huclist, filenmhuclist)
  saveRDS(guagehuc, filenmpoly)
  saveRDS(guagelist, filenmguagelist)
  #print(paste("Done with ", g, sep=""))
  gc()
}
## Setup up parallel list apply
snowfall::sfInit(parallel=TRUE, cpus=30) 
snowfall::sfExport("guages","hucs", "calreach.path", "huc_fn")
snowfall::sfLibrary(plyr)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(maptools)
snowfall::sfLibrary(raster)
snowfall::sfLibrary(ncdf4)
snowfall::sfLibrary(rgeos)
Sys.time()
snowfall::sfLapply(guages, function(g){huc_fn(g, hucs=hucs, calreach.path=calreach.path)})
Sys.time()
snowfall::sfStop()


### Now compute downstream load contributions for all guages after accounting for diversions
hucs$FRAC <- NULL ## Need updated version
huc_frac_tab <- read.delim("/home/tnaum/data/BLM_salinity/huc_frac_diversion.txt", stringsAsFactors = F)
hucs <- merge(hucs,huc_frac_tab, by="WATERID") # New diversion fractions provided by Matt Miller 3/6/2018
hucs_div <- subset(hucs, FRAC<1)
#hucs_div$FRAC <- hucs_div$FRAC + 0.01 # Eliminate errors from zeros.
hucdivlist <- hucs_div$WATERID
guages_df$upstrmFrac <- 1
## Loop to attribute diversion losses to upstream loads for calibration reach calcs
for(h in hucdivlist){
  starthuc <- hucs[hucs$WATERID %in% h,]
  startreach <- strm.reaches[strm.reaches$WATERID %in% h,]
  frac <- starthuc$FRAC
  startflow <- startreach$Corrected_Q_cms # F1
  flow <- startreach$Corrected_Q_cms*frac #F1,
  ## Key assumption that modeled Q doesn't account for Frac and needs to be adjusted down did have flow being *(1/frac) & the correction being flow-startflow
  startflow_correction = startflow-flow 
  newguage <- ifelse(starthuc$STAID > 0, starthuc$STAID)
  newguage <- newguage[!is.na(newguage)]
  guages_df$upstrmFrac <- ifelse(guages_df$guageid %in% newguage, frac, guages_df$upstrmFrac)
  if (length(newguage)>0){
    starthuc <- starthuc[FALSE,]
  }
  fnode <- starthuc$TNODE
  while(length(fnode)>0){
    newhuc <- hucs[hucs$FNODE %in% fnode,]
    newguage <- ifelse(newhuc$STAID > 0, newhuc$STAID)
    newguage <- newguage[!is.na(newguage)]
    if (length(newguage)>0){ ## Finding proportion of flow from headwater reaches to correct,
      wid <- newhuc$WATERID
      newreach <- strm.reaches[strm.reaches$WATERID %in% wid,]
      newflow <- newreach$Corrected_Q_cms + startflow_correction
      #flowratio <- flow/newflow ## old approach
      #newfrac <- 1-(flowratio*(1-frac)) ## old approach
      newfrac <- newflow/(newflow+startflow_correction)
      guages_df$upstrmFrac <- ifelse(guages_df$guageid == newguage, newfrac*guages_df$upstrmFrac, guages_df$upstrmFrac)
      newhuc <- newhuc[FALSE,]
    }
    fnode <- newhuc$TNODE
  }
  print(paste("Finished with ", h, sep=""))
  gc()
}
setwd('/home/tnaum/data/BLM_salinity/DSM_SPARROW')
write.table(guages_df, "UCRB_guages_with_DiversionFractions.txt", sep = "\t", row.names = FALSE)




#### Update data for calibration reach covariate and load summarization
guages_df <- subset(guages_df, !(is.na(guages_df$upstrmFrac)))
guages_df$div_adj_load_tonsyr <- guages_df$adj_mean_ds_tonyr*(1/guages_df$upstrmFrac) # Diversion correction
guages <- as.numeric(guages_df$guageid)

#### Summarize rasters and SPARROW HUC Yields by Guage HUCs and append to guages_df
huc_summary_fn <- function(g, hucs, guages_df){ # guages will be list
  hucpolypath <- paste("/home/tnaum/data/BLM_salinity/DSM_SPARROW/calibration_reaches/", g, "hucpoly.rds", sep="")
  huclistpath <- paste("/home/tnaum/data/BLM_salinity/DSM_SPARROW/calibration_reaches/", g, "huclist.rds", sep="")
  upstrguagespath <- paste("/home/tnaum/data/BLM_salinity/DSM_SPARROW/calibration_reaches/", g, "guagelist.rds", sep="")
  hucpoly <- readRDS(hucpolypath)
  huclist <- readRDS(huclistpath)
  upstrguagelist <- readRDS(upstrguagespath)
  upstrguages <- guages_df[guages_df$guageid %in% upstrguagelist,]
  calguage <- guages_df[guages_df$guageid %in% g,]
  upstr.load <- sum(upstrguages$adj_mean_ds_tonyr)
  cal.load <- calguage$div_adj_load_tonsyr - upstr.load
  hucsubset <- hucs[hucs$WATERID %in% huclist,]
  hucarea_sqm <- gArea(hucpoly)
  sqkm <- hucarea_sqm/1000000
  hucarea_sqmi <- sqkm*0.386102
  ag_load <- sum(hucsubset$Irrigated_Yield_tons_mi2/((hucsubset$SHAPE_Area/1000000)*0.386102))
  geo_load <- sum(hucsubset$Geologic_Yield_tons_mi2/((hucsubset$SHAPE_Area/1000000)*0.386102))
  ####  Rasters to summarize
  ### Source rasters
  ## Non irrigated ag
  ec0_10 <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec0_10s.tif")
  ec10_25 <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec10_25s.tif")
  ec25_50 <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec25_50s.tif")
  ec50_75 <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec50_75s.tif")
  ec75_90 <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec75_90s.tif")
  ec90_100 <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec90_100s.tif")
  ## Irrigated Ag
  # Flooded
  ec0_75F <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec0_75sF.tif")
  ec75_100F <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec75_100sF.tif")
  # Non-flooded irrigation type
  ec75_100N <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec75_100sN.tif")
  ec0_75N <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ec0_75sN.tif")
  ## Point Source springs
  sprg <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/pt_sourcerast.tif")
  ### Basin Characterization Rasters (BCMs)
  ## Main BCM variables
  ec50q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec50q.tif")
  ec75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q.tif")
  ec <- raster("/home/tnaum/data/BLM_salinity/risk_index/ecave_mask.tif")
  kw <- raster("/home/tnaum/data/BLM_salinity/risk_index/kw_mask.tif")
  ec90q <- raster("/home/tnaum/data/BLM_salinity/risk_index/ec90high.tif")
  kw75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_kw75q.tif")
  bgm <- raster("/home/tnaum/data/BLM_salinity/risk_index/bareground_mask.tif")
  flen <- raster("/home/tnaum/data/BLM_salinity/risk_index/flength_mask.tif")
  facc<- raster("/home/tnaum/data/BLM_salinity/risk_index/CAlog10_mask.tif")
  bgm75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/bgmacro75.tif")
  bgm90q <- raster("/home/tnaum/data/BLM_salinity/risk_index/bgmacro90_done.tif")
  bgm75q30p <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_bgm75q30p.tif")
  ## Risk indexes
  ec75q_kw75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_kw75q.tif")
  ec75q_facc75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_facc75q.tif")
  ec75q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_f500m.tif")
  ec75q_bgm75q30p <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_bgm75q30p.tif")
  ec90q_bgm90q40p_kw90q_facc90q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_b90m40p_facc90q_ec90q_kw90q_f500m.tif")
  ec75q_bgm75q30p_kw75q_facc75q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_bgm75q30p_kw75q_facc75q_f500m.tif")
  ec50q_bgm75q_kw75q_facc75q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec50q_bgm75_kw75q_facc75q_f500m.tif")
  ## USPED Salt Predictions
  edskg648 <- raster("/home/tnaum/data/BLM_salinity/risk_index/uspedsaltkg_m_maxmin648.tif")
  edskg648abs <- raster("/home/tnaum/data/BLM_salinity/risk_index/uspedsaltkg_m_maxmin648abs.tif")
  ## Other Soils variables
  Brock <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/SG250_Rprob.tif") # Prob of bedrock <2m
  sar <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/sar_m.tif") # Na Absorp ratio
  rock <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/rock_m.tif")# surface rock content
  fs <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/fs_vfs_m.tif") # % fine sand + vfs
  awc <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/awc_m.tif")
  ## DART variables
  elev <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/elevm_m.tif")
  ppt <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ppt_m.tif")
  pptratio <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/ppt_ratio_m.tif")
  protind <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/prot_index_m.tif")
  slp <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/slope_m.tif")
  sness <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/sness_m.tif") # raster of south vs north aspect
  ## Miller et al., (2017) BCMs as pulled from Flint and Flint (2007)
  exc <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/excess_m.tif")# excess water
  cwd <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/deficit_m.tif")# Climate water deficit
  mlt <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/melt_m.tif")# snowmelt
  rch <- raster("/home/tnaum/data/BLM_salinity/DSM_SPARROW/inputs/rech_m.tif")
  #### Guage raster prep
  hucpoly$rastfield <- 1
  e <- extent(hucpoly)
  ec0_10_hucrast <- crop(ec0_10,e)
  guagerast.path <- (paste("/home/tnaum/data/BLM_salinity/DSM_SPARROW/calibration_reach_rasters/", g, "rast.tif", sep=""))
  guagerast <- rasterize(hucpoly, ec0_10_hucrast,  field=hucpoly$rastfield, progress="text", datatype='INT1U', filename=guagerast.path) # If guage rasts don't exist
  #guagerast <- raster(guagerast.path) # for runs after the guage rasts have been made
  #### Process raster statistics for HUC
  f_mask <- function(a,b) a*b # masking function, one rast to mask and one to be masked, order irrelevant
  ### Sources
  ## EC source areas 0-10% quantiles, non ag
  ec0_10_stk <- stack(ec0_10_hucrast,guagerast)
  ec0_10_hucrast <- overlay(ec0_10_stk,fun=f_mask) #fast
  ec0_10.sqkm <- (cellStats(ec0_10_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec0_10_hucrast,ec0_10_stk)
  ## EC source areas 10-25% quantiles, non ag
  ec10_25_hucrast <- crop(ec10_25,e)
  ec10_25_stk <- stack(ec10_25_hucrast,guagerast)
  ec10_25_hucrast <- overlay(ec10_25_stk,fun=f_mask) #fast
  ec10_25.sqkm <- (cellStats(ec10_25_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec10_25_hucrast,ec10_25_stk)
  gc()
  ## EC source areas 25-50% quantiles, non ag
  ec25_50_hucrast <- crop(ec25_50,e)
  ec25_50_stk <- stack(ec25_50_hucrast,guagerast)
  ec25_50_hucrast <- overlay(ec25_50_stk,fun=f_mask) #fast
  ec25_50.sqkm <- (cellStats(ec25_50_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec25_50_hucrast,ec25_50_stk)
  ## EC source areas 50-75% quantiles, non ag
  ec50_75_hucrast <- crop(ec50_75,e)
  ec50_75_stk <- stack(ec50_75_hucrast,guagerast)
  ec50_75_hucrast <- overlay(ec50_75_stk,fun=f_mask) #fast
  ec50_75.sqkm <- (cellStats(ec50_75_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec50_75_hucrast,ec50_75_stk)
  ## EC source areas 75-90% quantiles, non ag
  ec75_90_hucrast <- crop(ec75_90,e)
  ec75_90_stk <- stack(ec75_90_hucrast,guagerast)
  ec75_90_hucrast <- overlay(ec75_90_stk,fun=f_mask) #fast
  ec75_90.sqkm <- (cellStats(ec75_90_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec75_90_hucrast,ec75_90_stk)
  ## EC source areas 90-100% quantiles, non ag
  ec90_100_hucrast <- crop(ec90_100,e)
  ec90_100_stk <- stack(ec90_100_hucrast,guagerast)
  ec90_100_hucrast <- overlay(ec90_100_stk,fun=f_mask) #fast
  ec90_100.sqkm <- (cellStats(ec90_100_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec90_100_hucrast,ec90_100_stk)
  gc()
  ## EC source areas 75-100% quantiles, flood irrigation ag
  ec75_100F_hucrast <- crop(ec75_100F,e)
  ec75_100F_stk <- stack(ec75_100F_hucrast,guagerast)
  ec75_100F_hucrast <- overlay(ec75_100F_stk,fun=f_mask) #fast
  ec75_100F.sqkm <- (cellStats(ec75_100F_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec75_100F_hucrast,ec75_100F_stk)
  ## EC source areas 0-75% quantiles, flood irrigation ag
  ec0_75F_hucrast <- crop(ec0_75F,e)
  ec0_75F_stk <- stack(ec0_75F_hucrast,guagerast)
  ec0_75F_hucrast <- overlay(ec0_75F_stk,fun=f_mask) #fast
  ec0_75F.sqkm <- (cellStats(ec0_75F_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec0_75F_hucrast,ec0_75F_stk)
  ## EC source areas 0-75% quantiles, non-flood irrigation ag
  ec0_75N_hucrast <- crop(ec0_75N,e)
  ec0_75N_stk <- stack(ec0_75N_hucrast,guagerast)
  ec0_75N_hucrast <- overlay(ec0_75N_stk,fun=f_mask) #fast
  ec0_75N.sqkm <- (cellStats(ec0_75N_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec0_75N_hucrast,ec0_75N_stk)
  ## EC source areas 75-100% quantiles, non-flood irrigation ag
  ec75_100N_hucrast <- crop(ec75_100N,e)
  ec75_100N_stk <- stack(ec75_100N_hucrast,guagerast)
  ec75_100N_hucrast <- overlay(ec75_100N_stk,fun=f_mask) #fast
  ec75_100N.sqkm <- (cellStats(ec75_100N_hucrast, stat='sum')*900)/1000000 #fast
  rm(ec75_100N_hucrast,ec75_100N_stk)
  ## Saline Spring point sources
  sprg_hucrast <- crop(sprg,e)
  sprg_stk <- stack(sprg_hucrast,guagerast)
  sprg_hucrast <- overlay(sprg_stk,fun=f_mask) #fast
  sprg.load <- cellStats(sprg_hucrast, stat='sum') #fast
  rm(sprg_hucrast,sprg_stk)
  gc()
  ### BCM variables
  ## EC >75th quantile % area
  ecq75q_hucrast <- crop(ec75q,e, progress="text") # pretty fast
  ec75q_stk <- stack(ecq75q_hucrast,guagerast)
  ec75q_hucrast <- overlay(ec75q_stk,fun=f_mask) #fast
  ec75q.pct <- (cellStats(ec75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_hucrast,ec75q_stk)
  ## EC >50th quantile % area
  ec50q_hucrast <- crop(ec50q,e)
  ec50q_stk <- stack(ec50q_hucrast,guagerast)
  ec50q_hucrast <- overlay(ec50q_stk,fun=f_mask) #fast
  ec50q.pct <- (cellStats(ec50q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec50q_hucrast,ec50q_stk)
  ## EC ave (0-30cm)
  ec_hucrast <- crop(ec,e, progress="text")
  ec_stk <- stack(ec_hucrast,guagerast)
  ec_hucrast <- overlay(ec_stk,fun=f_mask) #fast
  ec.ave <- cellStats(ec_hucrast, stat='mean') #fast
  ec.max <- cellStats(ec_hucrast, stat='max') #fast
  ec.75q <- unname(quantile(ec_hucrast, probs=0.75,na.rm=T))
  rm(ec_hucrast,ec_stk)
  ## Soil Erodibility (Kw) ave (0-30cm)
  kw_hucrast <- crop(kw,e, progress="text")
  kw_stk <- stack(kw_hucrast,guagerast)
  kw_hucrast <- overlay(kw_stk,fun=f_mask) #fast
  kw.ave <- cellStats(kw_hucrast, stat='mean') #fast
  kw.max <- cellStats(kw_hucrast, stat='max') #fast
  kw.75q <- unname(quantile(kw_hucrast, probs=0.75,na.rm=T))
  rm(kw_hucrast,kw_stk)
  gc()
  ## EC >90th quantile % area
  ec90q_hucrast <- crop(ec90q,e, progress="text")
  ec90q_stk <- stack(ec90q_hucrast,guagerast)
  ec90q_hucrast <- overlay(ec90q_stk,fun=f_mask) #fast
  ec90q.pct <- (cellStats(ec90q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec90q_hucrast,ec90q_stk)
  ## Soil Erodibility >75th quantile % area
  kw75q_hucrast <- crop(kw75q,e, progress="text")
  kw75q_stk <- stack(kw75q_hucrast,guagerast)
  kw75q_hucrast <- overlay(kw75q_stk,fun=f_mask) #fast
  kw75q.pct <- (cellStats(kw75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(kw75q_hucrast,kw75q_stk)
  ## Average % bareground exposed
  bgm_hucrast <- crop(bgm,e, progress="text")
  bgm_stk <- stack(bgm_hucrast,guagerast)
  bgm_hucrast <- overlay(bgm_stk,fun=f_mask) #fast
  bgm.ave <- cellStats(bgm_hucrast, stat='mean') # was wrong (sum) in 3/2 first run: can divide huc area by cell area to get pixel count to create average
  rm(bgm_hucrast,bgm_stk)
  ## Average flowlength (m)
  flen_hucrast <- crop(flen,e, progress="text")
  flen_stk <- stack(flen_hucrast,guagerast)
  flen_hucrast <- overlay(flen_stk,fun=f_mask) #fast
  flen.ave <- cellStats(flen_hucrast, stat='mean') # was wrong (sum) in 3/2 first run: can divide huc area by cell area to get pixel count to create average
  rm(flen_hucrast,flen_stk)
  gc()
  ## Average upstream (flow) area accumulation (m)
  facc_hucrast <- crop(facc,e, progress="text")
  facc_stk <- stack(facc_hucrast,guagerast)
  facc_hucrast <- overlay(facc_stk,fun=f_mask) #fast
  facc.ave <- cellStats(facc_hucrast, stat='mean') # was wrong (sum) in 3/2 first run: can divide huc area by cell area to get pixel count to create average
  rm(facc_hucrast,facc_stk)
  ## % of area with bareground >75th quantile within USGS GAP macrogroup
  bgm75q_hucrast <- crop(bgm75q,e, progress="text")
  bgm75q_stk <- stack(bgm75q_hucrast,guagerast)
  bgm75q_hucrast <- overlay(bgm75q_stk,fun=f_mask) #fast
  bgm75q.pct <- (cellStats(bgm75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(bgm75q_hucrast,bgm75q_stk)
  ## % of area with bareground >90th quantile within USGS GAP macrogroup
  bgm90q_hucrast <- crop(bgm90q,e, progress="text")
  bgm90q_stk <- stack(bgm90q_hucrast,guagerast)
  bgm90q_hucrast <- overlay(bgm90q_stk,fun=f_mask) #fast
  bgm90q.pct <- (cellStats(bgm90q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(bgm90q_hucrast,bgm90q_stk)
  ## % of area with bareground >75th quantile within USGS GAP macrogroup or > 30% absolute
  bgm75q30p_hucrast <- crop(bgm75q30p,e, progress="text")
  bgm75q30p_stk <- stack(bgm75q30p_hucrast,guagerast)
  bgm75q30p_hucrast <- overlay(bgm75q30p_stk,fun=f_mask) #fast
  bgm75q30p.pct <- (cellStats(bgm75q30p_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(bgm75q30p_hucrast,bgm75q30p_stk)
  gc()
  ### Complex risk indexes
  ## % of area with ec >75th quan, Kw > 75th quan
  ec75q_kw75q_hucrast <- crop(ec75q_kw75q,e, progress="text")
  ec75q_kw75q_stk <- stack(ec75q_kw75q_hucrast,guagerast)
  ec75q_kw75q_hucrast <- overlay(ec75q_kw75q_stk,fun=f_mask) #fast
  ec75q_kw75q.pct <- (cellStats(ec75q_kw75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_kw75q_hucrast,ec75q_kw75q_stk)
  ## % of area with ec >75th quan, facc > 75th quan
  ec75q_facc75q_hucrast <- crop(ec75q_facc75q,e, progress="text")
  ec75q_facc75q_stk <- stack(ec75q_facc75q_hucrast,guagerast)
  ec75q_facc75q_hucrast <- overlay(ec75q_facc75q_stk,fun=f_mask) #fast
  ec75q_facc75q.pct <- (cellStats(ec75q_facc75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_facc75q_hucrast,ec75q_facc75q_stk)
  ## % of area with ec >75th quan, flen < 500 meters
  ec75q_f500m_hucrast <- crop(ec75q_f500m,e, progress="text")
  ec75q_f500m_stk <- stack(ec75q_f500m_hucrast,guagerast)
  ec75q_f500m_hucrast <- overlay(ec75q_f500m_stk,fun=f_mask) #fast
  ec75q_f500m.pct <- (cellStats(ec75q_f500m_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_f500m_hucrast,ec75q_f500m_stk)
  ## % of area with ec >75th quan, bg >75th quan within USGS GAP macrogroup or > 30% absolute
  ec75q_bgm75q30p_hucrast <- crop(ec75q_bgm75q30p,e, progress="text")
  ec75q_bgm75q30p_stk <- stack(ec75q_bgm75q30p_hucrast,guagerast)
  ec75q_bgm75q30p_hucrast <- overlay(ec75q_bgm75q30p_stk,fun=f_mask) #fast
  ec75q_bgm75q30p.pct <- (cellStats(ec75q_bgm75q30p_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_bgm75q30p_hucrast,ec75q_bgm75q30p_stk)
  ## % of area with ec >90th quan, bg >90th quan within USGS GAP macrogroup or > 40% absolute, kw > 90th quan, facc > 90q, flen < 500m
  ec90q_bgm90q40p_kw90q_facc90q_f500m_hucrast <- crop(ec90q_bgm90q40p_kw90q_facc90q_f500m,e, progress="text")
  ec90q_bgm90q40p_kw90q_facc90q_f500m_stk <- stack(ec90q_bgm90q40p_kw90q_facc90q_f500m_hucrast,guagerast)
  ec90q_bgm90q40p_kw90q_facc90q_f500m_hucrast <- overlay(ec90q_bgm90q40p_kw90q_facc90q_f500m_stk,fun=f_mask) #fast
  ec90q_bgm90q40p_kw90q_facc90q_f500m.pct <- (cellStats(ec90q_bgm90q40p_kw90q_facc90q_f500m_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec90q_bgm90q40p_kw90q_facc90q_f500m_hucrast,ec90q_bgm90q40p_kw90q_facc90q_f500m_stk)
  ## % of area with ec >75th quan, bg >75th quan within USGS GAP macrogroup or > 30% absolute, kw > 75th quan, facc > 75 quan, flen < 500m
  ec75q_bgm75q30p_kw75q_facc75q_f500m_hucrast <- crop(ec75q_bgm75q30p_kw75q_facc75q_f500m,e, progress="text")
  ec75q_bgm75q30p_kw75q_facc75q_f500m_stk <- stack(ec75q_bgm75q30p_kw75q_facc75q_f500m_hucrast,guagerast)
  ec75q_bgm75q30p_kw75q_facc75q_f500m_hucrast <- overlay(ec75q_bgm75q30p_kw75q_facc75q_f500m_stk,fun=f_mask) #fast
  ec75q_bgm75q30p_kw75q_facc75q_f500m.pct <- (cellStats(ec75q_bgm75q30p_kw75q_facc75q_f500m_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_bgm75q30p_kw75q_facc75q_f500m_hucrast,ec75q_bgm75q30p_kw75q_facc75q_f500m_stk)
  ## % of area with ec >50th quan, bg >75th quan within USGS GAP macrogroup, kw > 75th quan, facc > 75 quan, flen < 500m
  ec50q_bgm75q_kw75q_facc75q_f500m_hucrast <- crop(ec50q_bgm75q_kw75q_facc75q_f500m,e, progress="text")
  ec50q_bgm75q_kw75q_facc75q_f500m_stk <- stack(ec50q_bgm75q_kw75q_facc75q_f500m_hucrast,guagerast)
  ec50q_bgm75q_kw75q_facc75q_f500m_hucrast <- overlay(ec50q_bgm75q_kw75q_facc75q_f500m_stk,fun=f_mask) #fast
  ec50q_bgm75q_kw75q_facc75q_f500m.pct <- (cellStats(ec50q_bgm75q_kw75q_facc75q_f500m_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec50q_bgm75q_kw75q_facc75q_f500m_hucrast,ec50q_bgm75q_kw75q_facc75q_f500m_stk)
  gc()
  ### USPED Predictions
  ## USPED average with erosion and deposition limited to theoretical 648 kg/pixel max (~30 cm loss with EC of 18 dS/m [max])
  edskg648_hucrast <- crop(edskg648,e, progress="text")
  edskg648_stk <- stack(edskg648_hucrast,guagerast)
  edskg648_hucrast <- overlay(edskg648_stk,fun=f_mask) #fast
  edskg648.ave <- cellStats(edskg648_hucrast, stat='mean')
  rm(edskg648_hucrast,edskg648_stk)
  ## USPED average with abs(erosion or deposition) limited to theoretical 648 kg/pixel max (~30 cm loss with EC of 18 dS/m [max])
  edskg648abs_hucrast <- crop(edskg648abs,e, progress="text")
  edskg648abs_stk <- stack(edskg648abs_hucrast,guagerast)
  edskg648abs_hucrast <- overlay(edskg648abs_stk,fun=f_mask) #fast
  edskg648abs.ave <- cellStats(edskg648abs_hucrast, stat='mean')
  rm(edskg648abs_hucrast,edskg648abs_stk)
  gc()
  ## Average probability of bedrock within 2m of soil surface - seems well linked to soil depth
  Brock_hucrast <- crop(Brock,e, progress="text")
  Brock_stk <- stack(Brock_hucrast,guagerast)
  Brock_hucrast <- overlay(Brock_stk,fun=f_mask) #fast
  Brock.ave <- cellStats(Brock_hucrast, stat='mean')
  rm(Brock_hucrast,Brock_stk)
  ## Average sodium absorpion ratio in surface soil horizon
  sar_hucrast <- crop(sar,e, progress="text")
  sar_stk <- stack(sar_hucrast,guagerast)
  sar_hucrast <- overlay(sar_stk,fun=f_mask) #fast
  sar.ave <- cellStats(sar_hucrast, stat='mean')
  rm(sar_hucrast,sar_stk)
  gc()
  ## Average % mass rock content in surface soil horizon
  rock_hucrast <- crop(rock,e, progress="text")
  rock_stk <- stack(rock_hucrast,guagerast)
  rock_hucrast <- overlay(rock_stk,fun=f_mask) #fast
  rock.ave <- cellStats(rock_hucrast, stat='mean')
  rm(rock_hucrast,rock_stk)
  ## Average % mass fine sand + very fine sand content in surface soil horizon
  fs_hucrast <- crop(fs,e, progress="text")
  fs_stk <- stack(fs_hucrast,guagerast)
  fs_hucrast <- overlay(fs_stk,fun=f_mask) #fast
  fs.ave <- cellStats(fs_hucrast, stat='mean')
  rm(fs_hucrast,fs_stk)
  gc()
  ## Average available water content from 15B (wilting point) to air dry (residual moisture) in surface soil horizon
  awc_hucrast <- crop(awc,e, progress="text")
  awc_stk <- stack(awc_hucrast,guagerast)
  awc_hucrast <- overlay(awc_stk,fun=f_mask) #fast
  awc.ave <- cellStats(awc_hucrast, stat='mean')
  rm(awc_hucrast,awc_stk)
  ### DART variables
  ## Average elevation (m)
  elev_hucrast <- crop(elev,e, progress="text")
  elev_stk <- stack(elev_hucrast,guagerast)
  elev_hucrast <- overlay(elev_stk,fun=f_mask) #fast
  elev.ave <- cellStats(elev_hucrast, stat='mean')
  rm(elev_hucrast,elev_stk)
  gc()
  ## Average annual precipitation 
  ppt_hucrast <- crop(ppt,e, progress="text")
  ppt_stk <- stack(ppt_hucrast,guagerast)
  ppt_hucrast <- overlay(ppt_stk,fun=f_mask) #fast
  ppt.ave <- cellStats(ppt_hucrast, stat='mean')
  rm(ppt_hucrast,ppt_stk)
  ## Average annual precipitation ratio (summer/annual precip) summer = June-Sept
  pptratio_hucrast <- crop(pptratio,e, progress="text")
  pptratio_stk <- stack(pptratio_hucrast,guagerast)
  pptratio_hucrast <- overlay(pptratio_stk,fun=f_mask) #fast
  pptratio.ave <- cellStats(pptratio_hucrast, stat='mean')
  rm(pptratio_hucrast,pptratio_stk)
  gc()
  ## Average topographic Protection Index (SAGA GIS)
  protind_hucrast <- crop(protind,e, progress="text")
  protind_stk <- stack(protind_hucrast,guagerast)
  protind_hucrast <- overlay(protind_stk,fun=f_mask) #fast
  protind.ave <- cellStats(protind_hucrast, stat='mean')
  rm(protind_hucrast,protind_stk)
  ## Average slope (deg)
  slp_hucrast <- crop(slp,e, progress="text")
  slp_stk <- stack(slp_hucrast,guagerast)
  slp_hucrast <- overlay(slp_stk,fun=f_mask) #fast
  slp.ave <- cellStats(slp_hucrast, stat='mean')
  rm(slp_hucrast,slp_stk)
  ## Average aspect on south-north axis 1 = south, -1 = north
  sness_hucrast <- crop(sness,e, progress="text")
  sness_stk <- stack(sness_hucrast,guagerast)
  sness_hucrast <- overlay(sness_stk,fun=f_mask) #fast
  sness.ave <- cellStats(sness_hucrast, stat='mean')
  rm(sness_hucrast,sness_stk)
  gc()
  ### Miller et al. (2017) BCMs from Flint and Flint (2007)
  ## ave Excess climatic water
  exc_hucrast <- crop(exc,e, progress="text")
  exc_stk <- stack(exc_hucrast,guagerast)
  exc_hucrast <- overlay(exc_stk,fun=f_mask) #fast
  exc.ave <- cellStats(exc_hucrast, stat='mean')
  rm(exc_hucrast,exc_stk)
  ## Ave climatic water deficit
  cwd_hucrast <- crop(cwd,e, progress="text")
  cwd_stk <- stack(cwd_hucrast,guagerast)
  cwd_hucrast <- overlay(cwd_stk,fun=f_mask) #fast
  cwd.ave <- cellStats(cwd_hucrast, stat='mean')
  rm(cwd_hucrast,cwd_stk)
  gc()
  ## Ave snowmelt
  mlt_hucrast <- crop(mlt,e, progress="text")
  mlt_stk <- stack(mlt_hucrast,guagerast)
  mlt_hucrast <- overlay(mlt_stk,fun=f_mask) #fast
  mlt.ave <- cellStats(mlt_hucrast, stat='mean')
  rm(mlt_hucrast,mlt_stk)
  ## Ave recharge
  rch_hucrast <- crop(rch,e, progress="text")
  rch_stk <- stack(rch_hucrast,guagerast)
  rch_hucrast <- overlay(rch_stk,fun=f_mask) #fast
  rch.ave <- cellStats(rch_hucrast, stat='mean')
  rm(rch_hucrast,rch_stk)
  gc()
  ## Pull all the parameters together
  hucdf <- data.frame(g,cal.load,ag_load,geo_load,sqkm,ec0_10.sqkm,ec10_25.sqkm,ec25_50.sqkm,ec50_75.sqkm,ec75_90.sqkm,ec90_100.sqkm,ec75_100F.sqkm,ec0_75F.sqkm,ec0_75N.sqkm,ec75_100N.sqkm,sprg.load,ec75q.pct,ec50q.pct,ec.ave,ec.max,ec.75q,kw.ave,kw.max,kw.75q,ec90q.pct,kw75q.pct,bgm.ave,flen.ave,facc.ave,bgm75q.pct,bgm90q.pct,bgm75q30p.pct,ec75q_kw75q.pct,ec75q_facc75q.pct,ec75q_f500m.pct,ec75q_bgm75q30p.pct,ec90q_bgm90q40p_kw90q_facc90q_f500m.pct,ec75q_bgm75q30p_kw75q_facc75q_f500m.pct,ec50q_bgm75q_kw75q_facc75q_f500m.pct,edskg648.ave,edskg648abs.ave,Brock.ave,sar.ave,rock.ave,fs.ave,awc.ave,elev.ave,ppt.ave,pptratio.ave,protind.ave,slp.ave,sness.ave,exc.ave,cwd.ave,mlt.ave,rch.ave)
  names(hucdf) <- c("guage","cal_tonsyr","ag_load","geo_load","sqkm","ec0_10.sqkm","ec10_25.sqkm","ec25_50.sqkm","ec50_75.sqkm","ec75_90.sqkm","ec90_100.sqkm","ec75_100F.sqkm","ec0_75F.sqkm","ec0_75N.sqkm","ec75_100N.sqkm","sprg.load","ec75q.pct","ec50q.pct","ec.ave","ec.max","ec.75q","kw.ave","kw.max","kw.75q","ec90q.pct","kw75q.pct","bgm.ave","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec50q_bgm75q_kw75q_facc75q_f500m.pct","edskg648.ave","edskg648abs.ave","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
  gc()
  return(hucdf)
}

## Setup up parallel list apply for guage huc summaries
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)
snowfall::sfInit(parallel=TRUE, cpus=30) 
snowfall::sfExport("guages","hucs", "huc_summary_fn", "guages_df")
snowfall::sfLibrary(plyr)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(maptools)
snowfall::sfLibrary(raster)
snowfall::sfLibrary(ncdf4)
snowfall::sfLibrary(rgeos)
Sys.time()
huc_sum <- snowfall::sfLapply(guages, function(g){try(huc_summary_fn(g, hucs=hucs, guages_df=guages_df))})
Sys.time()
snowfall::sfStop()
huc_sum_df <- huc_sum[[1]] ## must be data.frame
huc_sum_df <- huc_sum_df[FALSE,]
for(i in seq(1:length(huc_sum))){
  newrow <- huc_sum[[i]]
  if(class(newrow)=="data.frame"){
    huc_sum_df <- rbind(huc_sum_df, newrow)
  }
  print(paste("Done with ", i, sep=""))
}
huc_sum_df$guageid <- huc_sum_df$guage
huc_sum_guage_df = merge(huc_sum_df,guages_df, by="guageid")

## Create source percentage variables for direct yield-based modeling
huc_sum_guage_df$ec0_10.pct <- huc_sum_guage_df$ec0_10.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec10_25.pct <- huc_sum_guage_df$ec10_25.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec25_50.pct <- huc_sum_guage_df$ec25_50.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec50_75.pct <- huc_sum_guage_df$ec50_75.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec75_90.pct <- huc_sum_guage_df$ec75_90.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec90_100.pct <- huc_sum_guage_df$ec90_100.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec0_75F.pct <- huc_sum_guage_df$ec0_75F.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec75_100F.pct <- huc_sum_guage_df$ec75_100F.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec0_75N.pct <- huc_sum_guage_df$ec0_75N.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$ec75_100N.pct <- huc_sum_guage_df$ec75_100N.sqkm/huc_sum_guage_df$sqkm
huc_sum_guage_df$sprg_load.persqkm <- huc_sum_guage_df$sprg.load/huc_sum_guage_df$sqkm
huc_sum_guage_df$cal_tonsyrsqkm <- huc_sum_guage_df$cal_tonsyr/huc_sum_guage_df$sqkm

## Files
setwd('/home/tnaum/data/BLM_salinity/DSM_SPARROW')
write.table(huc_sum_guage_df, "UCRB_guages_DSM_SPARROW_oldKw_WithDiversions.txt", sep = "\t", row.names = FALSE)
huc_sum_guage_df <- read.delim("UCRB_guages_DSM_SPARROW_oldKw_WithDiversions.txt")

## Hucs from original network
allhuc_sum_df <- read.delim("UCRB_allHUCs_DSM_SPARROW_oldKw.txt")
## Merge with huc layer
hucs_w_covs = merge(hucs,allhuc_sum_df, by="WATERID")
hucs_w_covs$ec0_10.pct <- hucs_w_covs$ec0_10.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec10_25.pct <- hucs_w_covs$ec10_25.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec25_50.pct <- hucs_w_covs$ec25_50.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec50_75.pct <- hucs_w_covs$ec50_75.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec75_90.pct <- hucs_w_covs$ec75_90.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec90_100.pct <- hucs_w_covs$ec90_100.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec0_75F.pct <- hucs_w_covs$ec0_75F.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec75_100F.pct <- hucs_w_covs$ec75_100F.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec0_75N.pct <- hucs_w_covs$ec0_75N.sqkm/hucs_w_covs$sqkm
hucs_w_covs$ec75_100N.pct <- hucs_w_covs$ec75_100N.sqkm/hucs_w_covs$sqkm
hucs_w_covs$sprg_load.persqkm <- hucs_w_covs$sprg.load/hucs_w_covs$sqkm

## Still have some negative values...
huc_sum_guage_dfc <- subset(huc_sum_guage_df, cal_tonsyr > 0)

#### Random Forest Predictions ########################
library(randomForest)
varlist <- c("ec0_10.sqkm","ec10_25.sqkm","ec25_50.sqkm","ec50_75.sqkm","ec75_90.sqkm","ec90_100.sqkm","ec75_100F.sqkm","ec0_75F.sqkm","ec0_75N.sqkm","ec75_100N.sqkm","sprg.load","ec75q.pct","ec50q.pct","ec.ave","ec.75q","kw.ave","kw.75q","kw75q.pct","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec75q_bgm75q.pct","ec75q_bgm75q_kw75q.pct","ec75q_bgm75q_facc75q.pct","ec75q_bgm75q_f500m.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec75q_bgm75q_kw75q_facc75q_f500m.pct","ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
### RF for load: Log or not?
#varlist <- c("ec0_10.pct","ec10_25.pct","ec25_50.pct","ec50_75.pct","ec75_90.pct","ec90_100.pct","ec75_100F.pct","ec0_75F.pct","ec0_75N.pct","ec75_100N.pct","sprg.load","ec75q.pct","ec50q.pct","ec.ave","ec.max","ec.75q","kw.ave","kw.max","kw.75q","ec90q.pct","kw75q.pct","bgm.ave","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec50q_bgm75q_kw75q_facc75q_f500m.pct","edskg648.ave","edskg648abs.ave","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
formulaStringRF_adj_load <- as.formula(paste('log10(cal_tonsyr) ~', paste(varlist, collapse="+")))# put in dep variable name
adj_load_rf = randomForest(formulaStringRF_adj_load, data = huc_sum_guage_dfc, importance=TRUE, proximity=FALSE, ntree=200, keep.forest=TRUE, nodesize=1) 
adj_load_rf #summary: record for 10x runs then record 
## Check prediction of UCRB total load
allHucloadslog10 <- unname(predict(adj_load_rf, newdata=hucs_w_covs))
allHucloads <- 10^(as.numeric(allHucloadslog10))
UCRBload <- sum(allHucloads, na.rm=T) # Record for approach testing




### RF for yield ###########################
# Full list
varlist_yield <- c("ec0_10.pct","ec10_25.pct","ec25_50.pct","ec50_75.pct","ec75_90.pct","ec90_100.pct","ec75_100F.pct","ec0_75F.pct","ec0_75N.pct","ec75_100N.pct","sprg_load.persqkm","ec75q.pct","ec50q.pct","ec.ave","ec.75q","kw.ave","kw.75q","kw75q.pct","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec75q_bgm75q.pct","ec75q_bgm75q_kw75q.pct","ec75q_bgm75q_facc75q.pct","ec75q_bgm75q_f500m.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec75q_bgm75q_kw75q_facc75q_f500m.pct","ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
# Modified list
#varlist_yield <- c("ec0_10.sqkm","ec10_25.sqkm","ec25_50.sqkm","ec50_75.sqkm","ec75_90.sqkm","ec90_100.sqkm","ec75_100F.sqkm","ec0_75F.sqkm","ec0_75N.sqkm","ec75_100N.sqkm","sprg.load","ec0_10.pct","ec10_25.pct","ec25_50.pct","ec50_75.pct","ec75_90.pct","ec90_100.pct","ec75_100F.pct","ec0_75F.pct","ec0_75N.pct","ec75_100N.pct","sprg.load","ec75q.pct","ec50q.pct","ec.ave","ec.max","ec.75q","kw.ave","kw.max","kw.75q","ec90q.pct","kw75q.pct","bgm.ave","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec50q_bgm75q_kw75q_facc75q_f500m.pct","edskg648.ave","edskg648abs.ave","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
formulaStringRF_adj_yield <- as.formula(paste('log10(cal_tonsyrsqkm) ~', paste(varlist_yield, collapse="+")))# put in dep variable name
adj_yield_rf = randomForest(formulaStringRF_adj_yield, data = huc_sum_guage_dfc, importance=TRUE, proximity=FALSE, ntree=200, keep.forest=TRUE)
adj_yield_rf # Record for 10x runs for approach testing
varImpPlot(adj_yield_rf)
## Check prediction of UCRB total load
hucs_w_covs$adj_yield_rf_log10 <- unname(predict(adj_yield_rf, newdata=hucs_w_covs))
hucs_w_covs$adj_yield_rf_yield <- 10^(hucs_w_covs$adj_yield_rf_log10)
hucs_w_covs$adj_yield_rf_load <- hucs_w_covs$adj_yield_rf_yield*hucs_w_covs$sqkm
UCRBload <- sum(hucs_w_covs$adj_yield_rf_load, na.rm=T) # Record for approach testing

###################################################################
### Script to
### 1) integrate upstream reaches to stream guage pour points for FULL Watershed approach
### 2) Summarize raster layers by those upstream reaches
### 3) Adjust streamgauge loads by diversions
### 4) Correlate the raster summaries to stream gauge metrics
### 5) Create salinity yield and load random forests
### 6) Compute and graph results
##
###################################################################

## Load packages
required.packages <- c("raster", "sp", "rgdal","snow", "snowfall","parallel", "itertools","doParallel", "plyr", "ncdf4","maptools", "rgeos","stats","spdep","randomForest")# maybe need dplyr??
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)

#### Read in data sources
guages_df <- read.delim("/home/tnaum/data/BLM_Salinity/UCRB_streamguage_ec_ttab.txt", stringsAsFactors=F)
hucs <- readOGR("/home/tnaum/data/BLM_Salinity/HydroGDB/SIR20175009_UCRB_HydroNetwork.gdb/SIR20175009_UCRB_HydroNetwork.gdb", "sir20175009_UCRB_SPARROW_catchment_p")
strm.reaches <- readOGR("/home/tnaum/data/BLM_Salinity/HydroGDB/SIR20175009_UCRB_HydroNetwork.gdb/SIR20175009_UCRB_HydroNetwork.gdb", "sir20175009_UCRB_SPARROW_network_l")
hucs$STAID <- as.numeric(as.character(hucs$STAID))
guages <- as.numeric(guages_df$guageid)
guagehuc.path <- "/home/tnaum/data/BLM_salinity/guagehucs/"

#### Function/Loop to create shps of all reaches above each stream guage dissolved together (Full Watershed approach)
huc_fn <- function(g,hucs,guagehuc.path){
  #for(g in guages){
  huclist <- c()
  starthuc <- hucs[hucs$STAID %in% g,]
  huclist <- c(huclist,starthuc$WATERID)
  tnode <- starthuc$FNODE
  while(length(tnode)>0){
    newhucs <- hucs[hucs$TNODE %in% tnode,]
    huclist <- c(huclist, newhucs$WATERID)
    tnode <- newhucs$FNODE
  }
  guagehucs <- hucs[hucs$WATERID %in% huclist,]
  guagehucs$dissid <- 1
  guagehuc <- gUnaryUnion(guagehucs, id = guagehucs$dissid)
  filenmpoly <- paste(guagehuc.path,g,"hucpoly.rds",sep="")
  filenmhuclist <- paste(guagehuc.path,g,"huclist.rds",sep="")
  saveRDS(huclist, filenmhuclist)
  saveRDS(guagehuc, filenmpoly)
  #print(paste("Done with ", g, sep=""))
  gc()
}
## Setup up parallel list apply
snowfall::sfInit(parallel=TRUE, cpus=30) ## Choose number of cpus available
snowfall::sfExport("guages","hucs", "guagehuc.path", "huc_fn")
snowfall::sfLibrary(plyr)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(maptools)
snowfall::sfLibrary(raster)
snowfall::sfLibrary(ncdf4)
snowfall::sfLibrary(rgeos)
Sys.time()
# Worked on 2/13/18 with 2 errors
snowfall::sfLapply(guages, function(g){huc_fn(g, hucs=hucs, guagehuc.path=guagehuc.path)})
Sys.time()
snowfall::sfStop()



### Now compute downstream load contributions for all guages after accounting for diversions
hucs$FRAC <- NULL ## Need updated version
huc_frac_tab <- read.delim("/home/tnaum/data/BLM_Salinity/huc_frac_diversion.txt", stringsAsFactors = F)
huc_frac_tab$FRAC <- ifelse(huc_frac_tab$FRAC == 0, 0.01,huc_frac_tab$FRAC) # Adjust zero value to 99% diversion to avoid calculation errors
hucs <- merge(hucs,huc_frac_tab, by="WATERID") # New diversion fractions provided by Matt Miller 3/6/2018
hucs_div <- subset(hucs, FRAC<1)
#hucs_div$FRAC <- hucs_div$FRAC + 0.01 # Eliminate errors from zeros. ## Old approach before new diversion data provided
hucdivlist <- hucs_div$WATERID
guages_df$upstrmFrac <- 1
## Loop to attribute diversion losses to upstream loads for calibration reach calcs
for(h in hucdivlist){
  starthuc <- hucs[hucs$WATERID %in% h,]
  startreach <- strm.reaches[strm.reaches$WATERID %in% h,]
  frac <- starthuc$FRAC
  ## Key assumption that modeled Q is after diversion and needs to be adjusted up at 
  ## diversion reach and outlet to reflect the diverted flow in creating the load diversion factor.
  divflow <- (startreach$Corrected_Q_cms*(1/frac))-startreach$Corrected_Q_cms # amount of flow lost at diversion
  newguage <- ifelse(starthuc$STAID > 0, starthuc$STAID, NA)
  newguage <- newguage[!is.na(newguage)]
  fracadd <- ((startreach$Corrected_Q_cms+divflow)/startreach$Corrected_Q_cms)-1
  guages_df$upstrmFrac <- ifelse(guages_df$guageid %in% newguage, fracadd+guages_df$upstrmFrac, guages_df$upstrmFrac)
  fnode <- starthuc$TNODE
  while(length(fnode)>0){
    newhuc <- hucs[hucs$FNODE %in% fnode,]
    newguage <- ifelse(newhuc$STAID > 0, newhuc$STAID, NA)
    newguage <- newguage[!is.na(newguage)]
    if (length(newguage)>0){ ## Finding proportion of flow from headwater reaches to correct,
      wid <- newhuc$WATERID
      newreach <- strm.reaches[strm.reaches$WATERID %in% wid,]
      newfracadd <- ((newreach$Corrected_Q_cms+divflow)/newreach$Corrected_Q_cms)-1
      guages_df$upstrmFrac <- ifelse(guages_df$guageid == newguage, newfracadd+guages_df$upstrmFrac, guages_df$upstrmFrac)
    }
    fnode <- newhuc$TNODE
  }
  print(paste("Finished with ", h, sep=""))
  gc()
}
## Save file
setwd('/home/tnaum/data/BLM_Salinity/DSM_SPARROW')
write.table(guages_df, "UCRB_guages_with_DiversionFractions_fullwatersheds_20190320.txt", sep = "\t", row.names = FALSE)


#### Summarize rasters and SPARROW Reach Yields by Guage calibration reaches and append to guages_df
huc_summary_fn <- function(g, hucs){ # guages will be list
  hucpolypath <- paste("/home/tnaum/data/BLM_salinity/guagehucs/", g, "hucpoly.rds", sep="")
  huclistpath <- paste("/home/tnaum/data/BLM_salinity/guagehucs/", g, "huclist.rds", sep="")
  hucpoly <- readRDS(hucpolypath)
  huclist <- readRDS(huclistpath)
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
  bg <- raster("/home/tnaum/data/BLM_salinity/risk_index/bareground_mask.tif")
  flen <- raster("/home/tnaum/data/BLM_salinity/risk_index/flength_mask.tif")
  facc<- raster("/home/tnaum/data/BLM_salinity/risk_index/CAlog10_mask.tif")
  bgm75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/bgmacro75.tif")
  bgm90q <- raster("/home/tnaum/data/BLM_salinity/risk_index/bgmacro90_done.tif")
  bgm75q30p <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_bgm75q30p.tif")
  ## Risk indexes
  ec75q_kw75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_kw75q.tif")
  ec75q_facc75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_facc75q.tif")
  ec75q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_f500m.tif")
  ec75q_bgm75q30p <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_b75m30p.tif")
  ec75q_bgm75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_bgm75q.tif")
  ec75q_bgm75q_kw75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_bgm75q_kw75q.tif")
  ec75q_bgm75q_facc75q <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_bgm75q_facc75q.tif")
  ec75q_bgm75q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_bgm75q_f500m.tif")
  ec90q_bgm90q40p_kw90q_facc90q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec90q_b90m40p_kw90q_facc90q_f500m.tif")
  ec75q_bgm75q30p_kw75q_facc75q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_b75m30p_kw75q_facc75q_f500m.tif")
  ec75q_bgm75q_kw75q_facc75q_f500m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec75q_bgm75q_kw75q_facc75q_f500m.tif")
  ec50q_bgm75q30p_kw50q_facc50q_f1000m <- raster("/home/tnaum/data/BLM_salinity/risk_index/srisk_ec50q_bgm75q30p_kw50q_facc50q_f1000m.tif")
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
  guagerast.path <- (paste("/home/tnaum/data/BLM_salinity/guagerasts/", g, "rast.tif", sep=""))
  #guagerast <- rasterize(hucpoly, ec0_10_hucrast,  field=hucpoly$rastfield, progress="text", datatype='INT1U', filename=guagerast.path) # If guage rasts don't exist
  guagerast <- raster(guagerast.path) # for runs after the guage rasts have been made
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
  bg_hucrast <- crop(bg,e, progress="text")
  bg_stk <- stack(bg_hucrast,guagerast)
  bg_hucrast <- overlay(bg_stk,fun=f_mask) #fast
  bg.ave <- cellStats(bg_hucrast, stat='mean') 
  rm(bg_hucrast,bg_stk)
  ## Average flowlength (m)
  flen_hucrast <- crop(flen,e, progress="text")
  flen_stk <- stack(flen_hucrast,guagerast)
  flen_hucrast <- overlay(flen_stk,fun=f_mask) #fast
  flen.ave <- cellStats(flen_hucrast, stat='mean') 
  rm(flen_hucrast,flen_stk)
  gc()
  ## Average upstream (flow) area accumulation (m)
  facc_hucrast <- crop(facc,e, progress="text")
  facc_stk <- stack(facc_hucrast,guagerast)
  facc_hucrast <- overlay(facc_stk,fun=f_mask) #fast
  facc.ave <- cellStats(facc_hucrast, stat='mean') 
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
  ## % of area with ec >75th quan, bg >75th quan within USGS GAP macrogroup
  ec75q_bgm75q_hucrast <- crop(ec75q_bgm75q,e, progress="text")
  ec75q_bgm75q_stk <- stack(ec75q_bgm75q_hucrast,guagerast)
  ec75q_bgm75q_hucrast <- overlay(ec75q_bgm75q_stk,fun=f_mask) #fast
  ec75q_bgm75q.pct <- (cellStats(ec75q_bgm75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_bgm75q_hucrast,ec75q_bgm75q_stk)
  ## % of area with ec >75th quan, bg >75th quan within USGS GAP macrogroup, and kw >75th percentile
  ec75q_bgm75q_kw75q_hucrast <- crop(ec75q_bgm75q_kw75q,e, progress="text")
  ec75q_bgm75q_kw75q_stk <- stack(ec75q_bgm75q_kw75q_hucrast,guagerast)
  ec75q_bgm75q_kw75q_hucrast <- overlay(ec75q_bgm75q_kw75q_stk,fun=f_mask) #fast
  ec75q_bgm75q_kw75q.pct <- (cellStats(ec75q_bgm75q_kw75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_bgm75q_kw75q_hucrast,ec75q_bgm75q_kw75q_stk)
  ## % of area with ec >75th quan, bg >75th quan within USGS GAP macrogroup, and facc >75th percentile
  ec75q_bgm75q_facc75q_hucrast <- crop(ec75q_bgm75q_facc75q,e, progress="text")
  ec75q_bgm75q_facc75q_stk <- stack(ec75q_bgm75q_facc75q_hucrast,guagerast)
  ec75q_bgm75q_facc75q_hucrast <- overlay(ec75q_bgm75q_facc75q_stk,fun=f_mask) #fast
  ec75q_bgm75q_facc75q.pct <- (cellStats(ec75q_bgm75q_facc75q_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_bgm75q_facc75q_hucrast,ec75q_bgm75q_facc75q_stk)
  ## % of area with ec >75th quan, bg >75th quan within USGS GAP macrogroup, and flen < 500m
  ec75q_bgm75q_f500m_hucrast <- crop(ec75q_bgm75q_f500m,e, progress="text")
  ec75q_bgm75q_f500m_stk <- stack(ec75q_bgm75q_f500m_hucrast,guagerast)
  ec75q_bgm75q_f500m_hucrast <- overlay(ec75q_bgm75q_f500m_stk,fun=f_mask) #fast
  ec75q_bgm75q_f500m.pct <- (cellStats(ec75q_bgm75q_f500m_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_bgm75q_f500m_hucrast,ec75q_bgm75q_f500m_stk)
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
  ## % of area with ec >75th quan, bg >75th quan within USGS GAP macrogroup, kw > 75th quan, facc > 75 quan, flen < 500m
  ec75q_bgm75q_kw75q_facc75q_f500m_hucrast <- crop(ec75q_bgm75q_kw75q_facc75q_f500m,e, progress="text")
  ec75q_bgm75q_kw75q_facc75q_f500m_stk <- stack(ec75q_bgm75q_kw75q_facc75q_f500m_hucrast,guagerast)
  ec75q_bgm75q_kw75q_facc75q_f500m_hucrast <- overlay(ec75q_bgm75q_kw75q_facc75q_f500m_stk,fun=f_mask) #fast
  ec75q_bgm75q_kw75q_facc75q_f500m.pct <- (cellStats(ec75q_bgm75q_kw75q_facc75q_f500m_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec75q_bgm75q_kw75q_facc75q_f500m_hucrast,ec75q_bgm75q_kw75q_facc75q_f500m_stk)
  ## % of area with ec >50th quan, bg >75th quan within USGS GAP macrogroup, kw > 75th quan, facc > 75 quan, flen < 500m
  ec50q_bgm75q30p_kw50q_facc50q_f1000m_hucrast <- crop(ec50q_bgm75q30p_kw50q_facc50q_f1000m,e, progress="text")
  ec50q_bgm75q30p_kw50q_facc50q_f1000m_stk <- stack(ec50q_bgm75q30p_kw50q_facc50q_f1000m_hucrast,guagerast)
  ec50q_bgm75q30p_kw50q_facc50q_f1000m_hucrast <- overlay(ec50q_bgm75q30p_kw50q_facc50q_f1000m_stk,fun=f_mask) #fast
  ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct <- (cellStats(ec50q_bgm75q30p_kw50q_facc50q_f1000m_hucrast, stat='sum')*900)/hucarea_sqm #fast
  rm(ec50q_bgm75q30p_kw50q_facc50q_f1000m_hucrast,ec50q_bgm75q30p_kw50q_facc50q_f1000m_stk)
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
  hucdf <- data.frame(g,ag_load,geo_load,sqkm,ec0_10.sqkm,ec10_25.sqkm,ec25_50.sqkm,ec50_75.sqkm,ec75_90.sqkm,ec90_100.sqkm,ec75_100F.sqkm,ec0_75F.sqkm,ec0_75N.sqkm,ec75_100N.sqkm,sprg.load,ec75q.pct,ec50q.pct,ec.ave,ec.max,ec.75q,kw.ave,kw.max,kw.75q,ec90q.pct,kw75q.pct,bg.ave,flen.ave,facc.ave,bgm75q.pct,bgm90q.pct,bgm75q30p.pct,ec75q_kw75q.pct,ec75q_facc75q.pct,ec75q_f500m.pct,ec75q_bgm75q30p.pct,ec75q_bgm75q.pct,ec75q_bgm75q_kw75q.pct,ec75q_bgm75q_facc75q.pct,ec75q_bgm75q_f500m.pct,ec90q_bgm90q40p_kw90q_facc90q_f500m.pct,ec75q_bgm75q30p_kw75q_facc75q_f500m.pct,ec75q_bgm75q_kw75q_facc75q_f500m.pct,ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct,Brock.ave,sar.ave,rock.ave,fs.ave,awc.ave,elev.ave,ppt.ave,pptratio.ave,protind.ave,slp.ave,sness.ave,exc.ave,cwd.ave,mlt.ave,rch.ave)
  names(hucdf) <- c("guage","ag_load","geo_load","sqkm","ec0_10.sqkm","ec10_25.sqkm","ec25_50.sqkm","ec50_75.sqkm","ec75_90.sqkm","ec90_100.sqkm","ec75_100F.sqkm","ec0_75F.sqkm","ec0_75N.sqkm","ec75_100N.sqkm","sprg.load","ec75q.pct","ec50q.pct","ec.ave","ec.max","ec.75q","kw.ave","kw.max","kw.75q","ec90q.pct","kw75q.pct","bg.ave","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec75q_bgm75q.pct","ec75q_bgm75q_kw75q.pct","ec75q_bgm75q_facc75q.pct","ec75q_bgm75q_f500m.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec75q_bgm75q_kw75q_facc75q_f500m.pct","ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
  gc()
  return(hucdf)
}

## Setup up parallel list apply for guage huc summaries
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)
snowfall::sfInit(parallel=TRUE, cpus=30) # Choose available cpus
snowfall::sfExport("guages","hucs", "huc_summary_fn")
snowfall::sfLibrary(plyr)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(maptools)
snowfall::sfLibrary(raster)
snowfall::sfLibrary(ncdf4)
snowfall::sfLibrary(rgeos)
Sys.time()
huc_sum <- snowfall::sfLapply(guages, function(g){try(huc_summary_fn(g, hucs=hucs))}) # This can take a while, depending on cpus available
Sys.time()
snowfall::sfStop()
huc_sum_df <- huc_sum[[1]] ## must be data.frame
huc_sum_df <- huc_sum_df[FALSE,]
for(i in seq(2:length(huc_sum))){
  newrow <- huc_sum[[i]]
  if(class(newrow)=="data.frame"){
    huc_sum_df <- rbind(huc_sum_df, newrow)
  }
  print(paste("Done with ", i, sep=""))
}
huc_sum_df$guageid <- huc_sum_df$guage
huc_sum_guage_df = merge(huc_sum_df,guages_df, by="guageid")

### Save/reopen new summary data frame
setwd('/home/tnaum/data/BLM_Salinity/DSM_SPARROW')
write.table(huc_sum_guage_df, "UCRB_guages_DSM_SPARROW_oldKw_fullwatersheds_2013.txt", sep = "\t", row.names = FALSE)
huc_sum_guage_df <- read.delim("UCRB_guages_DSM_SPARROW_oldKw_fullwatersheds_2013.txt")

## Create source percentage variables for yield
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
huc_sum_guage_df$adj_mean_ds_tonyrsqkm <- huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm

## Hucs from original network (HUC_integration_allHUCs_DSM_SPARROW_oldKw_2013_datarelease.R)
allhuc_sum_df <- read.delim("UCRB_allHUCs_DSM_SPARROW_oldKw_2013.txt") # Has covariates summed by all incremental reaches
## Merge with huc layer
hucs_w_covs <- merge(hucs,allhuc_sum_df, by="WATERID")
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

## Guages with diversions
guages_div <- read.delim("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/UCRB_guages_with_DiversionFractions_fullwatersheds_20190320.txt", stringsAsFactors=F)
guages_div$div_adj_load_tonsyr <- guages_div$adj_mean_ds_tonyr*guages_div$upstrmFrac
guages_div <- subset(guages_div, select=c("guageid","div_adj_load_tonsyr"))
huc_sum_guage_dfc <- merge(huc_sum_guage_df,guages_div, by="guageid")
huc_sum_guage_dfc <- huc_sum_guage_dfc[!is.na(huc_sum_guage_dfc$div_adj_load_tonsyr),]
huc_sum_guage_dfc$div_adj_load_tonsyrsqkm <- huc_sum_guage_dfc$div_adj_load_tonsyr / huc_sum_guage_dfc$sqkm

#### Random Forest Predictions
varlist <- c("ec0_10.sqkm","ec10_25.sqkm","ec25_50.sqkm","ec50_75.sqkm","ec75_90.sqkm","ec90_100.sqkm","ec75_100F.sqkm","ec0_75F.sqkm","ec0_75N.sqkm","ec75_100N.sqkm","sprg.load","ec75q.pct","ec50q.pct","ec.ave","ec.75q","kw.ave","kw.75q","kw75q.pct","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec75q_bgm75q.pct","ec75q_bgm75q_kw75q.pct","ec75q_bgm75q_facc75q.pct","ec75q_bgm75q_f500m.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec75q_bgm75q_kw75q_facc75q_f500m.pct","ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
### RF for load: Log or not?
formulaStringRF_adj_load <- as.formula(paste('log10(div_adj_load_tonsyr) ~', paste(varlist, collapse="+")))# put in dep variable name
adj_load_rf <- randomForest(formulaStringRF_adj_load, data = huc_sum_guage_dfc, importance=TRUE, proximity=FALSE, ntree=200, keep.forest=TRUE, na.rm=T)
adj_load_rf #summary
varImpPlot(adj_load_rf)
plot(log10(huc_sum_guage_dfc$adj_mean_ds_tonyr) ~ predict(adj_load_rf, newdata=huc_sum_guage_dfc))
x1 <-c(0,100,10000,100000000)
y1 <-c(0,100,10000,100000000)
lines(x1,y1, col = 'red')#1:1 line

huc_sum_guage_dfc$adj_load_rf_log10pre <- unname(predict(adj_load_rf, newdata=huc_sum_guage_dfc))
attach(huc_sum_guage_dfc)
adj_load_rf_lm <- lm(log10(div_adj_load_tonsyr) ~ adj_load_rf_log10pre)
detach(huc_sum_guage_dfc)
summary(adj_load_rf_lm)
plot(log10(huc_sum_guage_dfc$div_adj_load_tonsyr) ~ predict(adj_load_rf_lm, newdata=huc_sum_guage_dfc))
lines(x1,y1, col = 'red')
## Check prediction of UCRB total load with no linear adjustment
allHucloadslog10 <- unname(predict(adj_load_rf, newdata=hucs_w_covs))
allHucloads <- 10^(as.numeric(allHucloadslog10))
UCRBload.adj_load_rf <- sum(allHucloads, na.rm=T)
## Check prediction of UCRB total load with linear adjustment
hucs_w_covs$adj_load_rf_log10pre <- unname(predict(adj_load_rf, newdata=hucs_w_covs))
hucs_w_covs$adj_load_rf_log10 <- predict(adj_load_rf_lm, hucs_w_covs)
hucs_w_covs$adj_load_rf_load <- 10^(hucs_w_covs$adj_load_rf_log10)
UCRBloadadj_load_rf_lm <- sum(allHucloads, na.rm=T)


## OOB Predictions: Load
huc_sum_guage_df$predOOB_log10_adj_load_rf = predict(adj_load_rf)
huc_sum_guage_df$predOOB_adj_load_rf = 10^(huc_sum_guage_df$predOOB_log10_adj_load_rf)
adj_load_rf_log10.OOB.RMSE = sqrt(mean((log10(huc_sum_guage_df$adj_mean_ds_tonyr) - huc_sum_guage_df$predOOB_log10_adj_load_rf)^2, na.rm=TRUE))
adj_load_rf_log10.OOB.Rsquared <- 1-var(log10(huc_sum_guage_df$adj_mean_ds_tonyr) - huc_sum_guage_df$predOOB_log10_adj_load_rf, na.rm=TRUE)/var(log10(huc_sum_guage_df$adj_mean_ds_tonyr), na.rm=TRUE)
adj_load_rf.OOB.RMSE = sqrt(mean((huc_sum_guage_df$adj_mean_ds_tonyr - huc_sum_guage_df$predOOB_adj_load_rf)^2, na.rm=TRUE))
adj_load_rf.OOB.Rsquared = 1-var(huc_sum_guage_df$adj_mean_ds_tonyr - huc_sum_guage_df$predOOB_adj_load_rf, na.rm=TRUE)/var(huc_sum_guage_df$adj_mean_ds_tonyr, na.rm=TRUE)
## OOB Convert to Yield
adj_load_rf_log10_yield.OOB.RMSE = sqrt(mean((log10(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - (log10(huc_sum_guage_df$predOOB_adj_load_rf/huc_sum_guage_df$sqkm)))^2, na.rm=TRUE))
adj_load_rf_log10_yield.OOB.Rsquared = 1-var(log10(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - log10(huc_sum_guage_df$predOOB_adj_load_rf/huc_sum_guage_df$sqkm), na.rm=TRUE)/var(log10(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm), na.rm=TRUE)
adj_load_rf_yield.OOB.RMSE = sqrt(mean(((huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - (huc_sum_guage_df$predOOB_adj_load_rf/huc_sum_guage_df$sqkm))^2, na.rm=TRUE))
adj_load_rf_yield.OOB.Rsquared = 1-var((huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - (huc_sum_guage_df$predOOB_adj_load_rf/huc_sum_guage_df$sqkm), na.rm=TRUE)/var(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm, na.rm=TRUE)

## Predictions: Load
huc_sum_guage_df$pred_log10_adj_load_rf = predict(adj_load_rf, newdata=huc_sum_guage_df)
huc_sum_guage_df$pred_adj_load_rf = 10^(huc_sum_guage_df$pred_log10_adj_load_rf)
adj_load_rf_log10.RMSE = sqrt(mean((log10(huc_sum_guage_df$adj_mean_ds_tonyr) - huc_sum_guage_df$pred_log10_adj_load_rf)^2, na.rm=TRUE))
adj_load_rf_log10.Rsquared = 1-var(log10(huc_sum_guage_df$adj_mean_ds_tonyr) - huc_sum_guage_df$pred_log10_adj_load_rf, na.rm=TRUE)/var(log10(huc_sum_guage_df$adj_mean_ds_tonyr), na.rm=TRUE)
adj_load_rf.RMSE = sqrt(mean((huc_sum_guage_df$adj_mean_ds_tonyr - huc_sum_guage_df$pred_adj_load_rf)^2, na.rm=TRUE))
adj_load_rf.Rsquared = 1-var(huc_sum_guage_df$adj_mean_ds_tonyr - huc_sum_guage_df$pred_adj_load_rf, na.rm=TRUE)/var(huc_sum_guage_df$adj_mean_ds_tonyr, na.rm=TRUE)
## Convert to Yield
adj_load_rf_log10_yield.RMSE = sqrt(mean((log10(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - (log10(huc_sum_guage_df$pred_adj_load_rf/huc_sum_guage_df$sqkm)))^2, na.rm=TRUE))
adj_load_rf_log10_yield.Rsquared = 1-var(log10(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - log10(huc_sum_guage_df$pred_adj_load_rf/huc_sum_guage_df$sqkm), na.rm=TRUE)/var(log10(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm), na.rm=TRUE)
adj_load_rf_yield.RMSE = sqrt(mean(((huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - (huc_sum_guage_df$pred_adj_load_rf/huc_sum_guage_df$sqkm))^2, na.rm=TRUE))
adj_load_rf_yield.Rsquared = 1-var((huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm) - (huc_sum_guage_df$pred_adj_load_rf/huc_sum_guage_df$sqkm), na.rm=TRUE)/var(huc_sum_guage_df$adj_mean_ds_tonyr/huc_sum_guage_df$sqkm, na.rm=TRUE)
plot(log10(huc_sum_guage_df$adj_mean_ds_tonyrsqkm)~log10(huc_sum_guage_df$pred_adj_load_rf/huc_sum_guage_df$sqkm))


##################### RF for yield ####################
# Full list
varlist_yield <- c("ec0_10.pct","ec10_25.pct","ec25_50.pct","ec50_75.pct","ec75_90.pct","ec90_100.pct","ec75_100F.pct","ec0_75F.pct","ec0_75N.pct","ec75_100N.pct","sprg_load.persqkm","ec75q.pct","ec50q.pct","ec.ave","ec.75q","kw.ave","kw.75q","kw75q.pct","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec75q_bgm75q.pct","ec75q_bgm75q_kw75q.pct","ec75q_bgm75q_facc75q.pct","ec75q_bgm75q_f500m.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec75q_bgm75q_kw75q_facc75q_f500m.pct","ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
formulaStringRF_adj_yield <- as.formula(paste('log10(div_adj_load_tonsyrsqkm) ~', paste(varlist_yield, collapse="+")))# put in dep variable name
adj_yield_rf <- randomForest(formulaStringRF_adj_yield, data = huc_sum_guage_dfc, importance=TRUE, proximity=FALSE, ntree=200, keep.forest=TRUE,nodesize=1) # Run 10x for model framework testing
adj_yield_rf ## log OOB 10x for framework testing
# Check UCRB total load prediction
hucs_w_covs$adj_yield_rf_log10test <- unname(predict(adj_yield_rf, newdata=hucs_w_covs))
hucs_w_covs$adj_yield_rf_yield_test <- 10^(hucs_w_covs$adj_yield_rf_log10test)
hucs_w_covs$adj_yield_rf_load_test <- hucs_w_covs$adj_yield_rf_yield_test*hucs_w_covs$sqkm
UCRBload.yield_rf_test <- sum(hucs_w_covs$adj_yield_rf_load_test, na.rm=T) ## Total UCRB load test, 15 (out of 10879) NAs come up in edge reaches
## Pruning process
adj_yield_rf <- randomForest(formulaStringRF_adj_yield, data = huc_sum_guage_dfc, importance=TRUE, proximity=FALSE, ntree=500, keep.forest=TRUE,nodesize=1,mtry=7) # rerun with varying tune parameters until a relatively high OOB error is found, then prune
adj_yield_rf #summary:
varImpPlot(adj_yield_rf)#Choose importance break: 15 was chosen when running for paper
# Prune list visuallly by importance
varlist_yieldc <- names(sort(adj_yield_rf$importance[,1], decreasing=T)[0:15]) # chose top 15 variables
formulaStringRF_adj_yieldc <- as.formula(paste('log10(div_adj_load_tonsyrsqkm) ~', paste(varlist_yieldc, collapse="+")))# put in dep variable name
adj_yield_rfc <- randomForest(formulaStringRF_adj_yieldc, data = huc_sum_guage_dfc, importance=TRUE, proximity=FALSE, ntree=500, keep.forest=TRUE,nodesize=1,mtry=7)
adj_yield_rfc #Run 10x to see if prune increase OOB Rsq by 3,  after 10 runs, keep running until a high OOB Rsq (relative to 1st 10 runs) is acheived and use that model in next prune step
varImpPlot(adj_yield_rfc)
# Prune by one variable at a time until max oob Rsq acheived
varlist_yieldcc <- names(sort(adj_yield_rf$importance[,1], decreasing=T)[0:14]) # chose top 14 variables
formulaStringRF_adj_yieldcc <- as.formula(paste('log10(div_adj_load_tonsyrsqkm) ~', paste(varlist_yieldcc, collapse="+")))# put in dep variable name
adj_yield_rfcc = randomForest(formulaStringRF_adj_yieldcc, data = huc_sum_guage_dfc, importance=TRUE, proximity=FALSE, ntree=500, keep.forest=TRUE,nodesize=1,mtry=7)
adj_yield_rfcc #Run 10x to see if prune increase OOB Rsq at all, if pruning helps, run these 5 steps again to test pruning another variable, and keep pruning until OOB doesn't improve
varImpPlot(adj_yield_rfcc)
adj_yield_rf <- adj_yield_rfcc ## rename pruned model for future steps
varImpPlot(adj_yield_rf)
partialPlot(adj_yield_rf, pred.data = huc_sum_guage_dfc, x.var=awc.ave) ## Check pi plot for most influential variable(s)
## Start plotting fin
huc_sum_guage_dfc$adj_yield_rf_log10pre <- unname(predict(adj_yield_rf, newdata=huc_sum_guage_dfc))
plot(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm) ~ huc_sum_guage_dfc$adj_yield_rf_log10pre)
lines(x1,y1, col = 'red')
attach(huc_sum_guage_dfc)
adj_yield_rf_lm <- lm(log10(div_adj_load_tonsyrsqkm) ~ adj_yield_rf_log10pre)
detach(huc_sum_guage_dfc)
summary(adj_yield_rf_lm)
plot(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm) ~ predict(adj_yield_rf_lm, newdata=huc_sum_guage_dfc))
lines(x1,y1, col = 'red')
## Check prediction of UCRB total load
hucs_w_covs$adj_yield_rf_log10pre <- unname(predict(adj_yield_rf, newdata=hucs_w_covs))
hucs_w_covs$adj_yield_rf_log10lm <- predict(adj_yield_rf_lm, hucs_w_covs) #data.frame(pprobs = x1); data.frame(adj_yield_rf_log10 = hucs_w_covs$adj_yield_rf_log10)
hucs_w_covs$adj_yield_rf_yield <- 10^(hucs_w_covs$adj_yield_rf_log10lm)
hucs_w_covs$adj_yield_rf_load <- hucs_w_covs$adj_yield_rf_yield*hucs_w_covs$sqkm
UCRBload.yield_rf_lm <- sum(hucs_w_covs$adj_yield_rf_load, na.rm=T)
UCRBload.yield_rf_lm_metric <- UCRBload.yield_rf_lm*0.90718474 ## metric tons
UCRBload.yield_rf_lm <- sum(hucs_w_covs$adj_yield_rf_load, na.rm=T)
saveRDS(adj_yield_rf, "/home/tnaum/data/BLM_salinity/DSM_SPARROW/scripts/FullWatershed_adj_yield_RF_2013.rds")
saveRDS(adj_yield_rf_lm, "/home/tnaum/data/BLM_salinity/DSM_SPARROW/scripts/FullWatershed_adj_yield_RFlm_2013.rds")
adj_yield_rf <- readRDS("/home/tnaum/data/BLM_salinity/DSM_SPARROW/scripts/FullWatershed_adj_yield_RF_2013.rds")
adj_yield_rf_lm <- readRDS("/home/tnaum/data/BLM_salinity/DSM_SPARROW/scripts/FullWatershed_adj_yield_RFlm_2013.rds")
## Fit metrics
huc_sum_guage_dfc$adj_yield_rf_yield_log10 <- predict(adj_yield_rf_lm, newdata=huc_sum_guage_dfc)
huc_sum_guage_dfc$adj_yield_rf_pre.nonlog <- 10^(huc_sum_guage_dfc$adj_yield_rf_log10pre)
huc_sum_guage_dfc$adj_yield_rf_pre.nonlog_m <- huc_sum_guage_dfc$adj_yield_rf_pre.nonlog * 0.90718474 ## Convert to metric tons
adj_yield_rf_log10pre.RMSE_m = sqrt(mean((log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm * 0.90718474) - log10(huc_sum_guage_dfc$adj_yield_rf_pre.nonlog_m))^2, na.rm=TRUE))
adj_yield_rf_log10pre.Rsquared_m = 1-var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474) - log10(huc_sum_guage_dfc$adj_yield_rf_pre.nonlog_m), na.rm=TRUE)/var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474), na.rm=TRUE)
huc_sum_guage_dfc$adj_yield_rf_yield.nonlog <- 10^(huc_sum_guage_dfc$adj_yield_rf_yield_log10)
huc_sum_guage_dfc$adj_yield_rf_yield.nonlog_m <- huc_sum_guage_dfc$adj_yield_rf_yield.nonlog * 0.90718474 ## Convert to metric tons
adj_yield_rf_log10.RMSE_m = sqrt(mean((log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474) - log10(huc_sum_guage_dfc$adj_yield_rf_yield.nonlog_m))^2, na.rm=TRUE))
adj_yield_rf_log10.Rsquared_m = 1-var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474) - log10(huc_sum_guage_dfc$adj_yield_rf_yield.nonlog_m), na.rm=TRUE)/var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474), na.rm=TRUE)
huc_sum_guage_dfc$rf_oob_log10yield <- predict(adj_yield_rf)
huc_sum_guage_dfc$rf_oob_log10yield.nonlog <- 10^(huc_sum_guage_dfc$rf_oob_log10yield)
huc_sum_guage_dfc$rf_oob_log10yield.nonlog_m <- huc_sum_guage_dfc$rf_oob_log10yield.nonlog * 0.90718474 ## Convert to metric tons
adj_yield_rf_log10.OOB.RMSE_m = sqrt(mean((log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474) - log10(huc_sum_guage_dfc$rf_oob_log10yield.nonlog_m))^2, na.rm=TRUE))
adj_yield_rf_log10.OOB.Rsquared_m = 1-var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474) - log10(huc_sum_guage_dfc$rf_oob_log10yield.nonlog_m), na.rm=TRUE)/var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474), na.rm=TRUE)


## Cross Validate the adjusted Random Forest and the linear model
nfolds <- 10
huc_sum_guage_dfc$folds <- sample.int(nfolds,size =length(huc_sum_guage_dfc$guageid),replace=T)
huc_sum_guage_dfc$lmcvpredpre <- "NA"
huc_sum_guage_dfc$lmcvpred <- "NA"
for (g in seq(nfolds)){
  traindf <- subset(huc_sum_guage_dfc, huc_sum_guage_dfc$folds != g)
  rf.pcvm <- randomForest(formulaStringRF_adj_yield, data=traindf, importance=FALSE, proximity=FALSE, ntree=500, keep.forest=TRUE,nodesize=1,mtry=7)
  huc_sum_guage_dfc$lmcvpredpre <- ifelse(huc_sum_guage_dfc$folds == g, predict(rf.pcvm, newdata=huc_sum_guage_dfc),huc_sum_guage_dfc$lmcvpredpre)
  traindf$lmcvpredpre <- predict(rf.pcvm, newdata=traindf)
  attach(traindf)
  lm.pcvm <- lm(log10(div_adj_load_tonsyrsqkm) ~ lmcvpredpre)
  detach(traindf)
  huc_sum_guage_dfc$lmcvpred <- ifelse(huc_sum_guage_dfc$folds == g, predict(adj_yield_rf_lm, data.frame(adj_yield_rf_log10pre = as.numeric(huc_sum_guage_dfc$lmcvpredpre))),huc_sum_guage_dfc$lmcvpred)
  print(g)
}
huc_sum_guage_dfc$lmcvpredpre = as.numeric(huc_sum_guage_dfc$lmcvpredpre)
huc_sum_guage_dfc$lmcvpred = as.numeric(huc_sum_guage_dfc$lmcvpred)
huc_sum_guage_dfc$lmcvpred.nonlog <- 10^(huc_sum_guage_dfc$lmcvpred)
huc_sum_guage_dfc$lmcvpred.nonlog_m <- huc_sum_guage_dfc$lmcvpred.nonlog * 0.90718474 # metric conversion
lmcvm.RMSE_m <- sqrt(mean((log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474) - log10(huc_sum_guage_dfc$lmcvpred.nonlog_m))^2, na.rm=TRUE))
lmcvm.Rsquared_m = 1-var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474) - log10(huc_sum_guage_dfc$lmcvpred.nonlog_m), na.rm=TRUE)/var(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm*0.90718474), na.rm=TRUE)

 
## Look at residuals and autocorrelation
huc_sum_guage_dfc$adj_yield_rf_lm_resid <- residuals(adj_yield_rf_lm)
## Save dataframe for future use
saveRDS(huc_sum_guage_dfc,"/home/tnaum/data/BLM_salinity/DSM_SPARROW/scripts/FullWatershed_adj_yield_RF_2013_dataframe_wlmCV.rds")
write.table(huc_sum_guage_dfc, "/home/tnaum/data/BLM_salinity/DSM_SPARROW/guages_yield_rflm_wholewatershed_2013_validated_ttab.txt", sep = "\t", row.names = FALSE)
huc_sum_guage_dfc <- read.delim("/home/tnaum/data/BLM_salinity/DSM_SPARROW/guages_yield_rflm_wholewatershed_2013_validated_ttab.txt")
## Plots for paper: model fit and cross validation
setwd("C:/Users/Travis/Dropbox/USGS/BLM_projects/Utah_BLM_Salinity/DSM_SPARROW_manuscript/figures/metric")
## Convert to metric
huc_sum_guage_dfc$div_adj_load_tonsyrsqkm_m <- huc_sum_guage_dfc$div_adj_load_tonsyrsqkm * 0.90718474
huc_sum_guage_dfc$adj_yield_rf_log10pre_m <- log10(10^(huc_sum_guage_dfc$adj_yield_rf_log10pre)*0.90718474)
huc_sum_guage_dfc$adj_yield_rf_yield_log10_m <- log10(10^(huc_sum_guage_dfc$adj_yield_rf_yield_log10)*0.90718474)
huc_sum_guage_dfc$adj_yield_rf_lm_resid_sign <- ifelse(huc_sum_guage_dfc$adj_yield_rf_lm_resid < 0, -1, 1)
huc_sum_guage_dfc$adj_yield_rf_lm_resid_unsigned <- abs(huc_sum_guage_dfc$adj_yield_rf_lm_resid)
huc_sum_guage_dfc$adj_yield_rf_lm_resid_nonlog <- 10^(huc_sum_guage_dfc$adj_yield_rf_lm_resid_unsigned)
huc_sum_guage_dfc$adj_yield_rf_lm_resid_nonlog_m <- huc_sum_guage_dfc$adj_yield_rf_lm_resid_nonlog * 0.90718474
huc_sum_guage_dfc$adj_yield_rf_lm_resid_m_log10unsigned <- log10(huc_sum_guage_dfc$adj_yield_rf_lm_resid_nonlog_m)
huc_sum_guage_dfc$adj_yield_rf_lm_resid_m <- huc_sum_guage_dfc$adj_yield_rf_lm_resid_m_log10unsigned * huc_sum_guage_dfc$adj_yield_rf_lm_resid_sign
huc_sum_guage_dfc$lmcvpred_m <- log10(10^(huc_sum_guage_dfc$lmcvpred)*0.90718474)
huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_sign <- ifelse(huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff < 0, -1, 1)
huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_unsigned <- abs(huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff)
huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_nonlog <- 10^(huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_unsigned)
huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_nonlog_m <- huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_nonlog * 0.90718474
huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_m_log10unsigned <- log10(huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_nonlog_m)
huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_m <- huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_m_log10unsigned * huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_sign
## Set up high res figure
tiff("FigX_Model_Performance.tif",width = 8, height = 5, units = 'in', res = 600 ) 
par(mfrow=c(2,3),mar=c(5,5,1,2),cex.lab = 1.0)
plot(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm_m) ~ huc_sum_guage_dfc$adj_yield_rf_log10pre_m,ylim=c(0,3),xlim=c(0,3),xlab=expression("Predicted Yield log10"(Mg ~ yr^{-1} ~ km^{-2})),ylab=expression("Observed Yield log10"(Mg ~ yr^{-1} ~ km^{-2})))
lines(x1,y1, col = 'red')
plot(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm_m) ~ huc_sum_guage_dfc$adj_yield_rf_yield_log10_m,ylim=c(0,3),xlim=c(0,3),xlab=expression("Corrected Yield log10"(Mg ~ yr^{-1} ~ km^{-2})),ylab=expression("Observed Yield log10"(Mg ~ yr^{-1} ~ km^{-2})))
lines(x1,y1, col = 'red')
plot(huc_sum_guage_dfc$adj_yield_rf_lm_resid_m~huc_sum_guage_dfc$adj_yield_rf_yield_log10_m,ylim=c(-1.2,1.2),xlim=c(0,3), xlab=expression("Corrected Yield log10"(Mg ~ yr^{-1} ~ km^{-2})),ylab="Residuals (In log10 units)")
x2 <-c(0,100,10000,100000000)
y2 <-c(0,0,0,0)
lines(x2,y2, col = 'red')#1:1 line
plot(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm_m) ~ log10(10^(predict(adj_yield_rf))*0.90718474),ylim=c(0,3),xlim=c(0,3),xlab=expression("RF Out-of-Bag Yield log10"(Mg ~ yr^{-1} ~ km^{-2})),ylab=expression("Observed Yield log10"(Mg ~ yr^{-1} ~ km^{-2})))
lines(x1,y1, col = 'red')
plot(log10(huc_sum_guage_dfc$div_adj_load_tonsyrsqkm_m) ~ huc_sum_guage_dfc$lmcvpred_m,ylim=c(0,3),xlim=c(0,3),xlab=expression("Cross Validation Yield log10"(Mg ~ yr^{-1} ~ km^{-2})),ylab=expression("Observed Yield log10"(Mg ~ yr^{-1} ~ km^{-2})))
lines(x1,y1, col = 'red')
plot(huc_sum_guage_dfc$adj_yield_rf_cv_lm_diff_m~huc_sum_guage_dfc$lmcvpred_m,ylim=c(-1.2,1.2),xlim=c(0,3), xlab=expression("Cross Validation Yield log10"(Mg ~ yr^{-1} ~ km^{-2})),ylab="Residuals (In log10 units)")
lines(x2,y2, col = 'red')#1:1 line
dev.off() ## finishes rendering high res figure

## Look at residuals relative to guage locations
guage.pts <- readOGR("C:/Dropbox/USGS/Utah_BLM_Salinity/SPARROW/HydroGDB/SIR20175009_UCRB_HydroNetwork.gdb/SIR20175009_UCRB_HydroNetwork.gdb", "guage_locations_rf_yield_model_alb")
guage.pts <- merge(guage.pts,huc_sum_guage_dfc, by="guage", all = FALSE)
# Write Residuals to a shapefile
guage.pts.residuals <- guage.pts[,c("WATERID","guage","adj_yield_rf_lm_resid")]
guage.pts.residuals$yldresid <- guage.pts.residuals$adj_yield_rf_lm_resid
guage.pts.residuals$adj_yield_rf_lm_resid <- NULL
#writeOGR(guage.pts.residuals, dsn="/home/tnaum/data/BLM_salinity/DSM_SPARROW/huc_predictions",driver="ESRI Shapefile", layer="guage_rf_lm_yield_resid_2013")
sa.dw1<-dnearneigh(guage.pts, d1=0, d2=1100000)
#sa.knn<-knearneigh(guage.pts, k = 10)
#sa.dw1<-knn2nb(sa.knn)
sa.dw1.dist <- nbdists(sa.dw1, coordinates(guage.pts))
all.dists <- c()
for(i in seq(1:length(sa.dw1.dist))){
  #sqdist <- as.data.frame(sa.dw1.dist[i])*as.data.frame(sa.dw1.dist[i])
  sa.dw1.dist[i] <- 1/as.data.frame(sa.dw1.dist[i])
  all.dists <- append(all.dists,as.vector(sa.dw1.dist[i]))
  #distances <- (as.data.frame(sa.dw1.dist[i]))^(-2)
  #sqdist <- as.data.frame(distances^2)
  #sa.dw1.dist[i] <- 1/sqdist
  print(i)
}
sa.wtd<-nb2listw(sa.dw1, style="B",glist=sa.dw1.dist ,zero.policy=NULL)
geary.test <- geary.test(guage.pts$adj_yield_rf_lm_resid, listw2U(sa.wtd), alternative="two.sided")
geary.test


## Huc predictions
hucs_w_covs$rfyield <- hucs_w_covs$adj_yield_rf_yield * 1
hucs_w_covs$rfload <- hucs_w_covs$adj_yield_rf_load * 1
#hucs_yield_exp <- hucs_w_covs[,c("WATERID","rfyield","rfload","sqkm")]
#writeOGR(hucs_yield_exp, dsn="/home/tnaum/data/BLM_salinity/DSM_SPARROW/huc_predictions",driver="ESRI Shapefile", layer="hucs_rflm_yield_2013")
## Compare Huc predictions to SPARROW
hucs_w_covs$SPARROWyield <- hucs_w_covs$Geologic_Yield_tons_mi2 + hucs_w_covs$Irrigated_Yield_tons_mi2
hucs_w_covs$sparrowyieldsqkm <- hucs_w_covs$SPARROWyield*0.386102
hucs_w_covs$sprwyld <- hucs_w_covs$sparrowyieldsqkm
hucs_w_covs$diff <- hucs_w_covs$rfyield - hucs_w_covs$sparrowyieldsqkm 
plot(log10(hucs_w_covs$sparrowyieldsqkm+0.01) ~ log10(hucs_w_covs$rfyield+0.01))
lines(x1,y1, col = 'red')
cor.test(log10(hucs_w_covs$sparrowyieldsqkm+0.01), log10(hucs_w_covs$rfyield+0.01))
plot(hucs_w_covs$sparrowyieldsqkm ~ hucs_w_covs$rfyield,xlim=c(0,1500))
lines(x1,y1, col = 'red')
hucs_yield_diff <- hucs_w_covs[,c("WATERID","rfyield","sprwyld","sqkm","diff")]
writeOGR(hucs_yield_diff, dsn="C:/Dropbox/USGS/Utah_BLM_Salinity/DSM_SPARROW_manuscript/data/huc_predictions",driver="ESRI Shapefile", layer="hucs_rflm_sparrow_diff_yield_2013")

## ggplot Hex bin plots
viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # Color ramp
# METRIC!!!
hucs_w_covs$sparrow.pltyield <- log10((hucs_w_covs$sparrowyieldsqkm*0.90718474)+0.1)
hucs_w_covs$rf.pltyield <- log10((hucs_w_covs$rfyield*0.90718474)+0.1)
hucs_w_covs.df <- as.data.frame(hucs_w_covs)
my_breaks <- c(50,100,150,200)
gplt.SPvRF <- ggplot(data=hucs_w_covs.df, aes(rf.pltyield,sparrow.pltyield)) +
  stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + xlim(-1.1,3.17) + ylim(-1.1,3.17) + 
  theme(axis.text=element_text(size=12), legend.text=element_text(size=14), axis.title=element_text(size=14),plot.title = element_text(hjust = 0.5)) + 
  xlab("RF log10(yield)") + ylab("Sparrow log10(yield)") + scale_fill_gradientn(name = "Count", colours = rev(viri))+ggtitle("Total Predicted Yields")
gplt.SPvRF
setwd("C:/Users/Travis/Dropbox/USGS/BLM_projects/Utah_BLM_Salinity/DSM_SPARROW_manuscript/figures/metric")
#ggsave('SPvRF_m.tif', plot = gplt.SPvRF, device = "tiff", dpi = 600, limitsize = TRUE)


## Simulations with reduced disturbance 
hucs_w_covs_dist_sim <- hucs_w_covs
hucs_w_covs_dist_sim$ec75q_bgm75q_kw75q.sqkm <- hucs_w_covs_dist_sim$ec75q_bgm75q_kw75q.pct * hucs_w_covs_dist_sim$sqkm + 0.0001
hucs_w_covs_dist_sim$ec75q_bgm75q_kw75q.pct <- hucs_w_covs_dist_sim$ec75q_bgm75q_kw75q.pct * 0
hucs_w_covs_dist_sim$adj_yield_rf_log10pre <- unname(predict(adj_yield_rf, newdata=hucs_w_covs_dist_sim))
hucs_w_covs_dist_sim$adj_yield_rf_log10lm <- predict(adj_yield_rf_lm, hucs_w_covs_dist_sim) #data.frame(pprobs = x1); data.frame(adj_yield_rf_log10 = hucs_w_covs$adj_yield_rf_log10)
hucs_w_covs_dist_sim$adj_yield_rf_yield <- 10^(hucs_w_covs_dist_sim$adj_yield_rf_log10lm)
hucs_w_covs_dist_sim$adj_yield_rf_load <- hucs_w_covs_dist_sim$adj_yield_rf_yield*hucs_w_covs_dist_sim$sqkm
UCRBload.yield_rf_lm_dist_sim <- sum(hucs_w_covs_dist_sim$adj_yield_rf_load, na.rm=T)
hucs_w_covs_dist_sim$distyld <-  hucs_w_covs_dist_sim$rfyield - hucs_w_covs_dist_sim$adj_yield_rf_yield
hucs_w_covs_dist_sim$distload <- hucs_w_covs_dist_sim$distyld * hucs_w_covs_dist_sim$sqkm
hucs_w_covs_dist_sim$distyld_perpix <- hucs_w_covs_dist_sim$distyld / hucs_w_covs_dist_sim$ec75q_bgm75q_kw75q.sqkm
hucs_w_covs_dist_sim$distload_perpix <- hucs_w_covs_dist_sim$distload / (hucs_w_covs_dist_sim$ec75q_bgm75q_kw75q.sqkm/0.0009)
hucs_w_covs_dist_sim_chnged <- subset(hucs_w_covs_dist_sim, hucs_w_covs_dist_sim$ec75q_bgm75q_kw75q.sqkm > 0.0001)
mean(hucs_w_covs_dist_sim_chnged$distload_perpix, na.rm=T) # average load saved by 'fixing' one identified disturbance pixel)
hucs_yield_distsim <- hucs_w_covs_dist_sim[,c("WATERID","rfyield","rfload","sqkm","distyld")]
writeOGR(hucs_yield_distsim, dsn="C:/Dropbox/USGS/Utah_BLM_Salinity/DSM_SPARROW_manuscript/data/huc_predictions",driver="ESRI Shapefile", layer="hucs_rflm_dist_yield_2013")

## Simulations with reduced flooded ag in high salinity areas 
hucs_w_covs_ag_sim <- hucs_w_covs
hucs_w_covs_ag_sim$ec75_100F.pct <- hucs_w_covs_ag_sim$ec75_100F.pct * 0
hucs_w_covs_ag_sim$adj_yield_rf_log10pre <- unname(predict(adj_yield_rf, newdata=hucs_w_covs_ag_sim))
hucs_w_covs_ag_sim$adj_yield_rf_log10lm <- predict(adj_yield_rf_lm, hucs_w_covs_ag_sim) #data.frame(pprobs = x1); data.frame(adj_yield_rf_log10 = hucs_w_covs$adj_yield_rf_log10)
hucs_w_covs_ag_sim$adj_yield_rf_yield <- 10^(hucs_w_covs_ag_sim$adj_yield_rf_log10lm)
hucs_w_covs_ag_sim$adj_yield_rf_load <- hucs_w_covs_ag_sim$adj_yield_rf_yield*hucs_w_covs_ag_sim$sqkm
UCRBload.yield_rf_lm_ag_sim <- sum(hucs_w_covs_ag_sim$adj_yield_rf_load, na.rm=T)
hucs_w_covs_ag_sim$agyld <-  hucs_w_covs_ag_sim$rfyield - hucs_w_covs_ag_sim$adj_yield_rf_yield
hucs_w_covs_ag_sim$agload <- hucs_w_covs_ag_sim$agyld * hucs_w_covs_ag_sim$sqkm
hucs_w_covs_ag_sim$agyld_perpix <- hucs_w_covs_ag_sim$agyld / hucs_w_covs_ag_sim$ec75_100F.sqkm
hucs_w_covs_ag_sim$agload_perpix <- hucs_w_covs_ag_sim$agload / ((hucs_w_covs_ag_sim$ec75_100F.sqkm+0.0001)/0.0009)
hucs_w_covs_ag_sim_chnged <- subset(hucs_w_covs_ag_sim, hucs_w_covs_dist_sim$ec75_100F.sqkm > 0.0001)
mean(hucs_w_covs_ag_sim_chnged$agload_perpix, na.rm=T) # average load saved by removing one flooded agricultural pixel)
hucs_yield_agsim <- hucs_w_covs_ag_sim[,c("WATERID","rfyield","rfload","sqkm","agyld")]
writeOGR(hucs_yield_agsim, dsn="C:/Dropbox/USGS/Utah_BLM_Salinity/DSM_SPARROW_manuscript/data/huc_predictions",driver="ESRI Shapefile", layer="hucs_rflm_ag_yield_2013")


## Simulations with reduced flooded ag in non-saline areas 
hucs_w_covs_ag_lowsalt_sim <- hucs_w_covs
hucs_w_covs_ag_lowsalt_sim$ec0_75F.pct <- hucs_w_covs_ag_lowsalt_sim$ec0_75F.pct * 0.1
hucs_w_covs_ag_lowsalt_sim$adj_yield_rf_log10pre <- unname(predict(adj_yield_rf, newdata=hucs_w_covs_ag_lowsalt_sim))
hucs_w_covs_ag_lowsalt_sim$adj_yield_rf_log10lm <- predict(adj_yield_rf_lm, hucs_w_covs_ag_lowsalt_sim) #data.frame(pprobs = x1); data.frame(adj_yield_rf_log10 = hucs_w_covs$adj_yield_rf_log10)
hucs_w_covs_ag_lowsalt_sim$adj_yield_rf_yield <- 10^(hucs_w_covs_ag_lowsalt_sim$adj_yield_rf_log10lm)
hucs_w_covs_ag_lowsalt_sim$adj_yield_rf_load <- hucs_w_covs_ag_lowsalt_sim$adj_yield_rf_yield*hucs_w_covs_ag_lowsalt_sim$sqkm
UCRBload.yield_rf_lm_ag_lowsalt_sim <- sum(hucs_w_covs_ag_lowsalt_sim$adj_yield_rf_load, na.rm=T)

## Simulations with no flooded ag in all soils 
hucs_w_covs_ag_allflooded_sim <- hucs_w_covs
hucs_w_covs_ag_allflooded_sim$ec0_75F.pct <- hucs_w_covs_ag_allflooded_sim$ec0_75F.pct * 0
hucs_w_covs_ag_allflooded_sim$ec75_100F.pct <- hucs_w_covs_ag_allflooded_sim$ec75_100F.pct * 0
hucs_w_covs_ag_allflooded_sim$adj_yield_rf_log10pre <- unname(predict(adj_yield_rf, newdata=hucs_w_covs_ag_allflooded_sim))
hucs_w_covs_ag_allflooded_sim$adj_yield_rf_log10lm <- predict(adj_yield_rf_lm, hucs_w_covs_ag_allflooded_sim) #data.frame(pprobs = x1); data.frame(adj_yield_rf_log10 = hucs_w_covs$adj_yield_rf_log10)
hucs_w_covs_ag_allflooded_sim$adj_yield_rf_yield <- 10^(hucs_w_covs_ag_allflooded_sim$adj_yield_rf_log10lm)
hucs_w_covs_ag_allflooded_sim$adj_yield_rf_load <- hucs_w_covs_ag_allflooded_sim$adj_yield_rf_yield*hucs_w_covs_ag_allflooded_sim$sqkm
UCRBload.yield_rf_lm_ag_allflooded_sim <- sum(hucs_w_covs_ag_allflooded_sim$adj_yield_rf_load, na.rm=T)
hucs_w_covs_ag_allflooded_sim$floodAg_RFyield <- hucs_w_covs_ag_allflooded_sim$rfyield - hucs_w_covs_ag_allflooded_sim$adj_yield_rf_yield ## 165 negative values, lowest -24.5, most higher than -4
hucs_w_covs_ag_allflooded_sim$floodAg_RFyield <- ifelse(hucs_w_covs_ag_allflooded_sim$floodAg_RFyield<0,0,hucs_w_covs_ag_allflooded_sim$floodAg_RFyield)
hucs_w_covs_ag_allflooded_sim$Flood_Irrigated_Yield_tons_sqkm <- hucs_w_covs_ag_allflooded_sim$Flood_Irrigated_Yield_tons_mi2 * 0.386102
# Plot flooded ag yields from Sparrow and the RF
plot(log10(hucs_w_covs_ag_allflooded_sim$Flood_Irrigated_Yield_tons_sqkm+0.1)~log10(hucs_w_covs_ag_allflooded_sim$floodAg_RFyield+0.1)) # compare flooded ag yields
lines(x1,y1, col = 'red')

# Plot Geo Yield versus the non ag yields produced by this simulation
hucs_w_covs_ag_allflooded_sim$Geologic_Yield_tons_sqkm <- hucs_w_covs_ag_allflooded_sim$Geologic_Yield_tons_mi2 * 0.386102
plot(log10(hucs_w_covs_ag_allflooded_sim$Geologic_Yield_tons_sqkm+0.1)~log10(hucs_w_covs_ag_allflooded_sim$adj_yield_rf_yield+0.1)) # compare flooded ag yields
lines(x1,y1, col = 'red')
## Now graph with ggplothexbin
hucs_w_covs_ag_allflooded_sim.df <- as.data.frame(hucs_w_covs_ag_allflooded_sim)
# Flood Ag
hucs_w_covs_ag_allflooded_sim.df$floodAg_RFyield.plt <- log10((hucs_w_covs_ag_allflooded_sim.df$floodAg_RFyield*0.90718474)+0.1)
hucs_w_covs_ag_allflooded_sim.df$Flood_Irrigated_Yield_tonnes_sqkm.plt <- log10((hucs_w_covs_ag_allflooded_sim$Flood_Irrigated_Yield_tons_sqkm*0.90718474)+0.1)
lim.yieldfa1 <- range(hucs_w_covs_ag_allflooded_sim.df$floodAg_RFyield.plt, na.rm=TRUE)
lim.yieldfa2 <- range(hucs_w_covs_ag_allflooded_sim.df$Flood_Irrigated_Yield_tonnes_sqkm.plt, na.rm=TRUE)
lim.yieldfa <-range(append(lim.yieldfa1,lim.yieldfa2),na.rm=T)
gplt.SPvRF.fa <- ggplot(data=hucs_w_covs_ag_allflooded_sim.df, aes(floodAg_RFyield.plt, Flood_Irrigated_Yield_tonnes_sqkm.plt)) +
  stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + xlim(-1.1,3.17) + ylim(-1.1,3.17) + 
  theme(axis.text=element_text(size=12), legend.text=element_text(size=14), axis.title=element_text(size=14),plot.title = element_text(hjust = 0.5)) + 
  xlab("RF log10(yield)") + ylab("SPARROW log10(yield)") + scale_fill_gradientn(name = "Count", colours = rev(viri)) +
  ggtitle("Yields from Flooded Agriculture")
gplt.SPvRF.fa
# Geologic sources
hucs_w_covs_ag_allflooded_sim.df$adj_yield_rf_yield.plt <- log10((hucs_w_covs_ag_allflooded_sim.df$adj_yield_rf_yield*0.90718474)+0.1)
hucs_w_covs_ag_allflooded_sim.df$Geologic_Yield_tonnes_sqkm.plt <- log10((hucs_w_covs_ag_allflooded_sim.df$Geologic_Yield_tons_sqkm*0.90718474)+0.1)
my_breaks <- c(50,100,150,200,250)
gplt.SPvRF.geo <- ggplot(data=hucs_w_covs_ag_allflooded_sim.df, aes(adj_yield_rf_yield.plt, Geologic_Yield_tonnes_sqkm.plt)) +
  stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + xlim(-1.1,3.17) + ylim(-1.1,3.17) + 
  theme(axis.text=element_text(size=12), legend.text=element_text(size=14), axis.title=element_text(size=14),plot.title = element_text(hjust = 0.5)) + 
  xlab("RF log10(yield)") + ylab("SPARROW log10(yield)") + scale_fill_gradientn(name = "Count", colours = rev(viri), breaks = my_breaks, labels = my_breaks) +
  ggtitle("Yields exluding Flooded Ag. Sources")
gplt.SPvRF.geo
##Trellis
SPvRF.plots <- grid.arrange(gplt.SPvRF,gplt.SPvRF.fa,gplt.SPvRF.geo,ncol=3)
SPvRF.plots
ggsave('FigX_YieldPlots_m.tif', plot = SPvRF.plots, device = "tiff", dpi = 600, limitsize = TRUE, width = 10, height = 2.75, units = 'in')

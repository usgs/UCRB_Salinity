###################################################################
### Script to
### 1) Summarize all incremental reaches by raster covariates to use in salinity predictions
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
hucs <- readOGR("/home/tnaum/data/BLM_Salinity/HydroGDB/SIR20175009_UCRB_HydroNetwork.gdb/SIR20175009_UCRB_HydroNetwork.gdb", "sir20175009_UCRB_SPARROW_catchment_p")
hucs$STAID <- as.numeric(as.character(hucs$STAID))
huclist <- hucs$WATERID

#### Summarize input rasters by incremental reach catchment (~11k)
huc_summary_fn <- function(h, hucs){ # huclist will be list
  hucpoly <- hucs[hucs$WATERID %in% h,]
  hucarea_sqm <- gArea(hucpoly)
  sqkm <- hucarea_sqm/1000000
  hucarea_sqmi <- sqkm*0.386102
  ####  Rasters to summarize
  ### Source rasters
  ## Non irrigated ag
  ec0_10 <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec0_10s.tif")
  ec10_25 <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec10_25s.tif")
  ec25_50 <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec25_50s.tif")
  ec50_75 <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec50_75s.tif")
  ec75_90 <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec75_90s.tif")
  ec90_100 <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec90_100s.tif")
  ## Irrigated Ag
  # Flooded
  ec0_75F <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec0_75sF.tif")
  ec75_100F <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec75_100sF.tif")
  # Non-flooded irrigation type
  ec75_100N <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec75_100sN.tif")
  ec0_75N <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ec0_75sN.tif")
  ## Point Source springs
  sprg <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/pt_sourcerast.tif")
  ### Basin Characterization Rasters (BCMs)
  ## Main BCM variables
  ec50q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec50q.tif")
  ec75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q.tif")
  ec <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ecave_mask.tif")
  kw <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/kw_m.tif")
  ec90q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec90q.tif")
  kw75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_kw75q.tif")
  bg <- raster("/home/tnaum/data/BLM_Salinity/risk_index/bareground_mask.tif")
  flen <- raster("/home/tnaum/data/BLM_Salinity/risk_index/flength_mask.tif")
  facc<- raster("/home/tnaum/data/BLM_Salinity/risk_index/CAlog10_mask.tif")
  bgm75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/bgmacro75.tif")
  bgm90q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/bgmacro90_done.tif")
  bgm75q30p <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_bgm75q30p.tif")
  ## Risk indexes
  ec75q_kw75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_kw75q.tif")
  ec75q_facc75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_facc75q.tif")
  ec75q_f500m <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_f500m.tif")
  ec75q_bgm75q30p <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_b75m30p.tif")
  ec75q_bgm75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_bgm75q.tif")
  ec75q_bgm75q_kw75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_bgm75q_kw75q.tif")
  ec75q_bgm75q_facc75q <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_bgm75q_facc75q.tif")
  ec75q_bgm75q_f500m <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_bgm75q_f500m.tif")
  ec90q_bgm90q40p_kw90q_facc90q_f500m <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec90q_b90m40p_kw90q_facc90q_f500m.tif")
  ec75q_bgm75q30p_kw75q_facc75q_f500m <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_b75m30p_kw75q_facc75q_f500m.tif")
  ec75q_bgm75q_kw75q_facc75q_f500m <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec75q_bgm75q_kw75q_facc75q_f500m.tif")
  ec50q_bgm75q30p_kw50q_facc50q_f1000m <- raster("/home/tnaum/data/BLM_Salinity/risk_index/srisk_ec50q_bgm75q30p_kw50q_facc50q_f1000m.tif")
  ## Other Soils variables
  Brock <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/SG250_Rprob.tif") # Prob of bedrock <2m
  sar <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/sar_m.tif") # Na Absorp ratio
  rock <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/rock_m.tif")# surface rock content
  fs <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/fs_m.tif") # % fine sand + vfs
  awc <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/awc_m.tif")
  ## DART variables
  elev <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/elevm_m.tif")
  ppt <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ppt_m.tif")
  pptratio <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/ppt_ratio_m.tif")
  protind <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/prot_index_m.tif")
  slp <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/slope_m.tif")
  sness <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/sness_m.tif") # raster of south vs north aspect
  ## Miller et al., (2017) BCMs as pulled from Flint and Flint (2007)
  exc <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/excess_m.tif")# excess water
  cwd <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/deficit_m.tif")# Climate water deficit
  mlt <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/melt_m.tif")# snowmelt
  rch <- raster("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/inputs/rech_m.tif")
  #### Guage raster prep
  hucpoly$rastfield <- 1
  e <- extent(hucpoly)
  ec0_10_hucrast <- crop(ec0_10,e)
  guagerast.path <- (paste("/home/tnaum/data/BLM_Salinity/DSM_SPARROW/hucrasts/hucrasts", h, "rast.tif", sep=""))
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
  hucdf <- data.frame(h,sqkm,ec0_10.sqkm,ec10_25.sqkm,ec25_50.sqkm,ec50_75.sqkm,ec75_90.sqkm,ec90_100.sqkm,ec75_100F.sqkm,ec0_75F.sqkm,ec0_75N.sqkm,ec75_100N.sqkm,sprg.load,ec75q.pct,ec50q.pct,ec.ave,ec.max,ec.75q,kw.ave,kw.max,kw.75q,ec90q.pct,kw75q.pct,bg.ave,flen.ave,facc.ave,bgm75q.pct,bgm90q.pct,bgm75q30p.pct,ec75q_kw75q.pct,ec75q_facc75q.pct,ec75q_f500m.pct,ec75q_bgm75q30p.pct,ec75q_bgm75q.pct,ec75q_bgm75q_kw75q.pct,ec75q_bgm75q_facc75q.pct,ec75q_bgm75q_f500m.pct,ec90q_bgm90q40p_kw90q_facc90q_f500m.pct,ec75q_bgm75q30p_kw75q_facc75q_f500m.pct,ec75q_bgm75q_kw75q_facc75q_f500m.pct,ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct,Brock.ave,sar.ave,rock.ave,fs.ave,awc.ave,elev.ave,ppt.ave,pptratio.ave,protind.ave,slp.ave,sness.ave,exc.ave,cwd.ave,mlt.ave,rch.ave)
  names(hucdf) <- c("WATERID","sqkm","ec0_10.sqkm","ec10_25.sqkm","ec25_50.sqkm","ec50_75.sqkm","ec75_90.sqkm","ec90_100.sqkm","ec75_100F.sqkm","ec0_75F.sqkm","ec0_75N.sqkm","ec75_100N.sqkm","sprg.load","ec75q.pct","ec50q.pct","ec.ave","ec.max","ec.75q","kw.ave","kw.max","kw.75q","ec90q.pct","kw75q.pct","bg.ave","flen.ave","facc.ave","bgm75q.pct","bgm90q.pct","bgm75q30p.pct","ec75q_kw75q.pct","ec75q_facc75q.pct","ec75q_f500m.pct","ec75q_bgm75q30p.pct","ec75q_bgm75q.pct","ec75q_bgm75q_kw75q.pct","ec75q_bgm75q_facc75q.pct","ec75q_bgm75q_f500m.pct","ec90q_bgm90q40p_kw90q_facc90q_f500m.pct","ec75q_bgm75q30p_kw75q_facc75q_f500m.pct","ec75q_bgm75q_kw75q_facc75q_f500m.pct","ec50q_bgm75q30p_kw50q_facc50q_f1000m.pct","Brock.ave","sar.ave","rock.ave","fs.ave","awc.ave","elev.ave","ppt.ave","pptratio.ave","protind.ave","slp.ave","sness.ave","exc.ave","cwd.ave","mlt.ave","rch.ave")
  gc()
  return(hucdf)
}

## Set up parallel list apply for huc summaries
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)
snowfall::sfInit(parallel=TRUE, cpus=30) 
snowfall::sfExport("hucs", "huc_summary_fn")
snowfall::sfLibrary(plyr)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(maptools)
snowfall::sfLibrary(raster)
snowfall::sfLibrary(ncdf4)
snowfall::sfLibrary(rgeos)
Sys.time()
huc_sum <- snowfall::sfLapply(huclist, function(h){try(huc_summary_fn(h, hucs=hucs))})
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


setwd('/home/tnaum/data/BLM_Salinity/UCRB_Salinity/SalinityYieldCode')
write.table(huc_sum_df, "UCRB_allHUCs_DSM_SPARROW_2013_NASIS.txt", sep = "\t", row.names = FALSE) 

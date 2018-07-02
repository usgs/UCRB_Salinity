The SoilMapPredictionFolder includes code and models for the new soil maps produced for salinity yield models. The field data used to train these models all comes from the USDA-NRCS national characterization database (NCD) (https://ncsslabdatamart.sc.egov.usda.gov/). The location table and various attribute tables were joined individually and then used for model building in the following scripts. The covariates used as predictors for these models are available via the following link. https://www.dropbox.com/sh/bdbpzok2rm2sybn/AACqvp1Mh4kn0ckFoql3JIMga?dl=0

If there are any issues accessing these data, please contact Travis Nauman at tnauman@usgs.gov

Files:

EC_model_00cm.R - Code for soil electrical conductivity prediction at 0 cm depth.

EC_model_15cm.R - Code for soil electrical conductivity prediction at 15 cm depth.

EC_model_30cm.R - Code for soil electrical conductivity prediction at 30 cm depth.

AWC_model_00cm.R - Code for air dry water content prediction at 00 cm depth.

FS_model_00cm.R - Code for percent fine + very fine sand content prediction at 00 cm depth.

SAR_model_00cm.R - Code for percent sodium adsorption ratio (sar) prediction at 00 cm depth.

Rock_model_00cm.R - Code for percent rock fragment by mass prediction at 00 cm depth.

Kfactor_model_00cm.R - Code for soil erodibility prediction at 00 cm depth.

UPCO_NCSSlabpeds_Kw_covars.txt - NCD points use to train K-factor model with environmental covariates extracted. This table is the training matrix for the randomForest model.

UCRB_ncss17_w15l2_covarsc.txt - NCD points use to train the air-dry water content (AWC) model with environmental covariates extracted. This table is the training matrix for the randomForest model.

cop_ncss17SAR_pct_covarsc.txt - NCD points use to train the sodium adsorption ratio (SAR) model with environmental covariates extracted. This table is the training matrix for the randomForest model.

cop_ncss17_FS_VFS_pct_covarsc.txt - NCD points use to train the fine + very fine sand (FS) model with environmental covariates extracted. This table is the training matrix for the randomForest model.

cop_ncss17_Rock_Frgmt_2D_covarsc.txt - NCD points use to train the rock fragments (Rock) model with environmental covariates extracted. This table is the training matrix for the randomForest model.

cop_ncss17_ec12_2D_covarsc.txt - NCD points use to train the electrical conductivity (ec) model with environmental covariates extracted. This table is the training matrix for the randomForest model.

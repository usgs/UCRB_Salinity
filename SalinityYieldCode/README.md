Folder containing code and models used to predict salinity yields

File list:

BetweenGaugeModels.R - Script documenting the between gauge incremental reach watershed approach for model calibration. Scipt includes incremental calibration reach determination, accounting for stream diversions in calibration, summarization of covariates for calibration reaches, and initial random forest model building for testing calibration reach approaches.

FullWatershedModels.R - Script documenting approach looking at the entire upstream watershed of each gauge for model calibration. Scipt includes calibration reach determination, accounting for stream diversions in calibration, summarization of covariates for calibration reaches, and initial random forest model building for testing calibration reach approaches. This script also includes the pruning process for the final modeling building along with the results statistics and figures made for the manuscript.

Reach_covariate_summarization.R - Script documenting the summarization of covariates for all incremental reaches in the Upper Colorado River Watershed for use in final salinity yield predictions. These reach polygons were produced by Miller et al., (2017) for SPARROW modeling, and we adapted them as the hydrological topology for use in our modeling.

FullWatershed_adj_yield_RF_2013.rds - The final purned randomForest model used for salinity yield predictions

FullWatershed_adj_yield_RFlm_2013.rds - The linear model used to adjust random forest predictions to correct for bias.

FullWatershed_adj_yield_RF_2013_dataframe_wlmCV.rds - The guage data table with cross validation predictions used in evaluating the linear-adjusted randomForest model used for making salinity yield predictions

Folder containing code and models used to predict salinity yields

File list:

BetweenGaugeModels.R - Script documenting the between gauge incremental reach watershed approach for model calibration. Scipt includes incremental calibration reach determination, accounting for stream diversions in calibration, summarization of covariates for calibration reaches, and initial random forest model building for testing calibration reach approaches.

FullWatershedModels.R - Script documenting approach looking at the entire upstream watershed of each gauge for model calibration. Scipt includes calibration reach determination, accounting for stream diversions in calibration, summarization of covariates for calibration reaches, and initial random forest model building for testing calibration reach approaches. This script also includes the pruning process for the final modeling building along with the results statistics and figures made for the manuscript.

Reach_covariate_summarization.R - Script documenting the summarization of covariates for all incremental reaches in the Upper Colorado River Watershed for use in final salinity yield predictions. These reach polygons were produced by Miller et al., (2017) for SPARROW modeling, and we adapted them as the hydrological topology for use in our modeling.

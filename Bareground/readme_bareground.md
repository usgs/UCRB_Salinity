This folder contains scripts and field data used to make predictive bareground maps for the upper Colorado River basin. These maps were then used in making risk indices used for salinity yield modeling.

Files:

Bareground_model.R - R script using landsat 8 growing season data from 2014 to create a model of exposed bareground based on BLM-AIM field plots. This data is then predicted onto 2013 Landsat 8 data to try to best match dates with salinity yield data. Both a random forest and a bi-directional AIC driven stepwise regression were tested - the stepwise regression provided the best result.

Landsat8_2013_GEE - Google Earth Engine javascript code to produce Landsat 8 2013 growing season median composite reflectance bands. This was used for bareground predictions

Landsat8_2014_GEE -  Google Earth Engine javascript code to produce Landsat 8 2014 growing season median composite reflectance bands. This was used for bareground model building.

AIM14_cover_landsat_covars.txt - This is the compiled BLM-AIM data from 2014 used to train the bareground model.

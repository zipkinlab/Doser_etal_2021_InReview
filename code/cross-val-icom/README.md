# Code for cross-validation using the integrated community occupancy model

+ `cross-val.R`: code to obtain predictive results from cross-validation for the White Mountain National Forest case study. 
+ `lpd-fit-neon-sampler.R`: function to obtain predictive performance metric for NEON data when the predictions are generated from a model using some other combination of data but the fitting model is NEON. 
+ `lpd-pred-neon-sampler.R`: function to obtain predictive performance metric for NEON data when the predictions are generated from a model using some other combination fo data and the fitting model does not only use NEON data. 
+ `lpd-sampler.R`: function to obtain predictive performance metrics at HBEF and BBS. 
+ All files beginning with `main` are used to run the ICOM through NIMBLE using the specific data sources listed in the file name.
+ All files beginning with `main` are files to run the ICOM for cross-validation using the data sources listed in the file names. See `README` file in `code/wmnf` for more details. 
+ `nimble-code`: directory containing files to run the ICOM in NIMBLE.
+ `summary.R`: file to summarize cross-validation results for foliage-gleaning bird case study. 

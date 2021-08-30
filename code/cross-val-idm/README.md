# Code for cross-validation using single species integrated distribution models

+ `cross-val.R`: code to obtain predictive results from cross-validation for the White Mountain National Forest case study using single-species IDMs.  
+ `lpd-sampler.R`: function to obtain predictive performance metrics. 
+ `main-HBEF-NEON-BBS.R`: code used to run single species IDMs using all three data sources for purposes of cross-validation. 
+ `main-HBEF-BBS.R`: code used to run single species IDMs using HBEF and BBS data for cross-validation. This is used for species that are not observed at NEON locations. 
+ `nimble-code`: directory containing files to run the IDM in NIMBLE.

Results from this directory are summarized in the file `code/cross-val-icom/summary.R` along with the cross-validation results from the ICOM. 


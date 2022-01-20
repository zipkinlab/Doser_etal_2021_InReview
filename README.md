# [Integrated community occupancy models: A framework to assess occurrence and biodiversity dynamics using multiple data sources](https://arxiv.org/abs/2109.01894)

### In press, Methods in Ecology and Evolution

### Jeffrey W. Doser, Wendy Leuenberger, T. Scott Sillett, Michael T. Hallworth, Elise F. Zipkin

### Code/Data DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5883950.svg)](https://doi.org/10.5281/zenodo.5883950)

### Please contact the first author for questions about the code or data: Jeffrey W. Doser (doserjef@msu.edu)
__________________________________________________________________________________________________________________________________________

## Abstract:
1. The occurrence and distributions of wildlife populations and communities are shifting as a result of global changes. To evaluate whether these shifts are negatively impacting biodiversity processes, it is critical to monitor the status, trends, and effects of environmental variables on entire communities. However, modeling the dynamics of multiple species simultaneously can require large amounts of diverse data, and few modeling approaches exist to simultaneously provide species and community level inferences.
2. We present an "integrated community occupancy model" (ICOM) that unites principles of data integration and hierarchical community modeling in a single framework to provide inferences on species-specific and community occurrence dynamics using multiple data sources. The ICOM combines replicated and nonreplicated detection-nondetection data sources using a hierarchical framework that explicitly accounts for different detection and sampling processes across data sources. We use simulations to compare the ICOM to previously developed hierarchical community occupancy models and single species integrated distribution models. We then apply our model to assess the occurrence and biodiversity dynamics of foliage-gleaning birds in the White Mountain National Forest in the northeastern USA from 2010-2018 using three independent data sources.
3. Simulations reveal that integrating multiple data sources in the ICOM increased precision and accuracy of species and community level inferences compared to single data source models, although benefits of integration were dependent on data source quality (e.g., amount of replication). Compared to single species models, the ICOM yielded more precise species-level estimates. Within our case study, the ICOM had the highest out-of-sample predictive performance compared to single species models and models that used only a subset of the three data sources. 
4. The ICOM provides more precise estimates of occurrence dynamics compared to multi-species models using single data sources or integrated single-species models. We further found that the ICOM had improved predictive performance across a broad region of interest with an empirical case study of forest birds. The ICOM offers an attractive approach to estimate species and biodiversity dynamics, which is additionally valuable to inform management objectives of both individual species and their broader communities.


## Repository Directory

### [Code](./code)

Contains code for the case study and simulations. See each subdirectory for README files describing the file contents. For running the ICOM and adapting it to different data sets, we recommend starting in the [simulations](./code/simulations) directory with [main.R](./code/simulations/main.R). 

## [Data](./data)

+ `final-bird-data.R`: the foliage-gleaning bird data used in the case study.
+ `nimble-data.R`: contains the same data as `final-bird-data.R` but is processed and stored for direct input into NIMBLE. 
+ `cross-val-indices.R`: site indexes used in the cross-validation for the foliage-gleaning bird case study.

## [Results](./results)

+ [cross-val](./results/cross-val): subdirectory containing posterior predictive performance metrics for different models in the foliage-gleaning bird case study. 
+ `icom-HBEF-results.R`: results of model fit using only HBEF data. 
+ `icom-NEON-results.R`: results of model fit using only NEON data. 
+ `icom-BBS-results.R`: results of model fit using only BBS data. 
+ `icom-HBEF-NEON-results.R`: results of model fit using HBEF and NEON data. 
+ `icom-HBEF-BBS-results.R`: results of model fit using HBEF and BBS data. 
+ `icom-HBEF-NEON-BBS-results.R`: results of model fit using HBEF, NEON, and BBS data. 


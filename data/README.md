## [Data](./data)

+ `final-bird-data.R`: the foliage-gleaning bird data used in the case study.
+ `nimble-data.R`: contains the same data as `final-bird-data.R` but is processed and stored for direct input into NIMBLE. Metadata on the data files contained in both `nimble-data.R` and `final-bird-data.R` are contained in the `main-*` files in the [wmnf code directory](https://github.com/zipkinlab/Doser_etal_2022_MEE/tree/main/code/wmnf).
+ `cross-val-indices.R`: site indexes used in the cross-validation for the foliage-gleaning bird case study. These are indices that are generated to ensure the same sites are held out when fitting the different models for cross-validation. 

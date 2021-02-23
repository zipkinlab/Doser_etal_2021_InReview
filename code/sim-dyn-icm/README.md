# Simulations for Integrated Dynamic Community Occupancy Model

1. `icm-jags.txt`: contains jags code for running integrated community model
2. `isdm-jags.txt`: contains jags code for running the integrated model for a single species. Used in model development. 
3. `main-icm.R`: script to run ICM and produce some initial plots assessing accuracy of the model. 
4. `main-isdm.R`: script to run ISDM. 
5. `main-sim-1-icm.R`: script to run multiple simulations of the ICM at once. This script is specifically used to see how estimates from the model change when differing amounts of repeated detection/non-detection data and NEON data are incorporated into the model. 
6. `main-sim-2-icm.R`: script to run multiple simulations of the ICM at once. This script is specifically used to see how varying the number of years in the model and the number of repeat visits for the Hubbard Brook and NEON data 
7. `results`: directory for results of `main-sim-1-icm.R` and `main-sim-2-icm.R`. 
8. `sim-icm-data.R`: function for simulating data for use in testing models. 
9. `figures`: directory for storing any figures. 
10. `summary-sim-1.R`: summarize results and produce figures for `main-sim-1-icm.R`.
11. `summary-sim-2.R`: summarize results and produce figures for `main-sim-2-icm.R`.

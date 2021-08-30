# Simulation Study for Integrated Community Occupancy Model

1. `main.R`: main R script used to run ICOM via NIMBLE through R using one replicated data set and two nonreplicated data sets. *We recommend using this script to learn how to run the ICOM and adapting the ICOM to specific case studies.* 
2. `main-sim-1-icom.R`: provides code used to implement the first simulation study described in Doser et al (2021) for assessing benefits of data integration in the ICOM. 
3. `main-sim-2-icom.R`: provides code used to implement the second simulation study described in Doser et al (2021) for assessing benefits of the community framework in the ICOM. Used to run the ICOM portion of the simulations. 
4. `main-sim-2-idm.R`: provides code used to implement the second simulation study described in Doser et al (2021) for assessing benefits of the community framework in the ICOM. Used to run the IDM portion of the simulations. 
5. `nimble-code`: directory containing BUGS code and associated information required for running different forms of the ICOM via NIMBLE in R. File names denote the specific model and the types of data included for each model. 
6. `sim-icom-data.R`: function used to simulate data from the ICOM. In particular, this function generates data from a single replicated data source and two nonreplicated data soruces. 
7. `summary-sim-1.R`: code used to summarize results from the first simulation study (`main-sim-1-icom.R`) and produce resulting figures and tables. 
8. `summary-sim-2.R`: code used to summarize results from the second simulation study (`main-sim-2-icom.R` and `main-sim-2-idm.R`) and produce resulting figures and tables. 


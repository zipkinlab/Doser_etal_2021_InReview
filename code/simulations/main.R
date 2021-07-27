# main.R: main R script to run ICOM via NIMBLE through R using one replicated
#         data set and two nonreplicated data sets. This script can be used for 
#         initial exploration of the ICOM and for adapting the ICOM towards 
#         individual case studies.
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(nimble)
library(coda)
library(tidyverse)
# Load function to simulate data. 
source("code/simulations/sim-icom-data.R")
# Set seed for constant results
#set.seed(1817)
# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Simulate Data -----------------------------------------------------------
# Data are simulated from a single replicated data source and two 
# nonreplicated data sources. 
# Number sites with replicated data source
J.rep <- 50
# Number sites with first nonreplicated data source
J.nrep.1 <- 50
# Number of sites with second nonreplicated data source
J.nrep.2 <- 50
# Total number of sites
J <- J.rep + J.nrep.1 + J.nrep.2
# Number of years
n.years <- 5
# Community level parameters ----------------------------------------------
# Occupancy intercept (vary by year)
beta.0.mean <- runif(n.years, -1.5, .15)
# Occupancy spatial covariate effect
beta.1.mean <- 0.1
# Community auto-logistic mean
phi.mean <- 1.4
# Replicated detection intercept (vary by year)
alpha.0.mean <- runif(n.years, -1, 1)
# Replicated detection space/time covariate effect
alpha.1.mean <- 0.3
# Nonreplicated 1 detection intercept (vary by year)
gamma.1.0.mean <- runif(n.years, -2.5, 0)
# Nonreplicated 1 detection space/time covariate effect
gamma.1.1.mean <- 0.1
# Nonreplicated 2 detection intercept (vary by year)
gamma.2.0.mean <- runif(n.years, 0, 2.5)
# Nonreplicated 2 detection space/time covariate effect
gamma.2.1.mean <- -0.2
# Number of replicates for replicated data set
K.rep <- 3
# Number of replicates for first nonreplicated data set (can change if 
# desired to have two replicated data sets)
K.nrep.1 <- 1
# Number of species
I <- 8
# Community-level variance parameters
sigma.sq.beta.0 <- runif(n.years, 0.5, 1.5)
sigma.sq.beta.1 <- runif(1, 0.25, 2)
sigma.sq.phi <- runif(1, 0.25, 1)
sigma.sq.rep.0 <- runif(n.years, 1, 3)
sigma.sq.rep.1 <- runif(1, 0.25, 2)
sigma.sq.nrep.1.0 <- runif(n.years, 1, 3)
sigma.sq.nrep.1.1 <- runif(1, 0.25, 2)
sigma.sq.nrep.2.0 <- runif(n.years, 1, 3)
sigma.sq.nrep.2.1 <- runif(1, 0.25, 2)

# Form species-specific covariates ----------------------------------------
beta.0 <- matrix(NA, nrow = I, ncol = n.years)
alpha.0 <- matrix(NA, nrow = I, ncol = n.years)
gamma.1.0 <- matrix(NA, nrow = I, ncol = n.years)
gamma.2.0 <- matrix(NA, nrow = I, ncol = n.years)
for (t in 1:n.years) {
  beta.0[, t] <- rnorm(I, beta.0.mean[t], sqrt(sigma.sq.beta.0[t]))
  alpha.0[, t] <- rnorm(I, alpha.0.mean[t], sqrt(sigma.sq.rep.0[t]))
  gamma.1.0[, t] <- rnorm(I, gamma.1.0.mean[t], sqrt(sigma.sq.nrep.1.0[t]))
  gamma.2.0[, t] <- rnorm(I, gamma.2.0.mean[t], sqrt(sigma.sq.nrep.2.0[t]))
}
beta.1 <- rnorm(I, beta.1.mean, sqrt(sigma.sq.beta.1)) 
phi <- rnorm(I, phi.mean, sqrt(sigma.sq.phi))
gamma.1.1 <- rnorm(I, gamma.1.1.mean, sqrt(sigma.sq.nrep.1.1))
alpha.1 <- rnorm(I, alpha.1.mean, sqrt(sigma.sq.rep.1))
gamma.2.1 <- rnorm(I, gamma.2.1.mean, sqrt(sigma.sq.nrep.2.1))

# Simulate the data using the data simulation function. 
dat <- sim.icom.data(J.rep, J.nrep.1, J.nrep.2, beta.0, beta.1, 
		    phi, gamma.2.0, gamma.2.1,
		    alpha.0, alpha.1, 
		    gamma.1.0, gamma.1.1, K.rep, K.nrep.1, I, 
		    n.years)
# dat contains the following objects: 
# X.psi: design matrix for occurrence
# X.nrep.2: design matrix for detection of nonreplicated 2 data set
# X.nrep.1: design matrix for detection of nonreplicated 1 data set
# X.rep: design matrix for detection of replicated data set
# psi: latent occurrence probabilities
# z: latent precense/absence array
# p: detection probabilities for replicated data set
# y: detection-nondetection data for replicated data set
# pi.1, pi.2: detection probabilities for nonreplicated data sets 1 and 2. 
# v.1, v.2: detection-nondetection data for nonreplicated data sets 1 and 2. 


# Load BUGS code and initial values/constants/data for use in NIMBLE ------ 
source("code/simulations/nimble-code/icom-rep-nrep1-nrep2.R")
start <- Sys.time()
# Create Model ------------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
  		          data = icom.data, inits = icom.inits)
  
# Configure MCMC ----------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('int.beta.mean', 'beta.1.mean', 
						  'phi.mean', 
						  'beta.0', 'beta.1', 'phi'))
# Create an MCMC function -------------------------------------------------
icom.mcmc <- buildMCMC(icom.conf)
# Compile model -----------------------------------------------------------
icom.c.model <- compileNimble(icom.model)
icom.c.mcmc <- compileNimble(icom.mcmc, project = icom.model)
# Number of iterations ----------------------------------------------------
n.iter <- 8000
n.thin <- 2
n.burn <- 5000
n.chain <- 1
# Run the model -----------------------------------------------------------
samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn,
   	       thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)
# Quick summary of model --------------------------------------------------
summary(samples)
plot(mcmc(samples), density = FALSE)

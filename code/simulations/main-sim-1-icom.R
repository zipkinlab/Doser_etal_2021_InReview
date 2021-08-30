# main-sim-1-icom.R: this script provides code to implement the first 
#                   simulation study described in Doser et al (2021). 
#                   The purpose of this simulation study is to assess
#                   the benefits of integrating multiple data sets together,
#                   specifically assessing the benefits of integrating
#                   one replicated data set, one nonreplicated data set 
#                   with low detection probability, and one nonreplicated
#                   data set with high detection probability. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(nimble)
library(coda)
# Load function to simulate data
source("code/simulations/sim-icom-data.R")
# This index is used to run simulations across multiple cores. Each 
# simulation scenario is repeated n.sims times on its own core. Numbers
# supplied to index when running should be between 1 and 7, which is the 
# number of unique combinations of the three data sources. 
# Values for index correspond to the following models: 
# 1: replicated data only
# 2: nonreplicated data with low detection only
# 3: nonreplicated data with high detection only
# 4: replicated + nonreplicated with low detection
# 5: replicated + nonreplicated with high detection
# 6: nonreplicated with low detection + nonreplicated with high detection
# 7: replicated + nonreplicated with low detection + nonreplicated with 
#    high detection
# The second number when running on the command line should correspond
# to the number of the specific chain you want to run for a given model. 
# This allows you to run the model with different starting values for multiple
# chains across different cores instead of running the chains sequentially, 
# as is the default scenario for NIMBLE. 
# For example: to run the model using all three data sources for chain 1, 
#              run the following on the command line from the home directory 
#              for this project.
#    Rscript code/simulations/main-sim-1-icom.R 7 1 & 
index <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(index) != 2) base::stop('Need to tell NIMBLE which model and chain to run')
# Alternatively, can specify which simulation scenario to run manually. This
# is useful for testing, or if you don't have familiarity with the command line. 
# index <- c(6, 1)
# Determine file with NIMBLE code based on user input value
code.files <- c('icom-rep.R', 'icom-nrep1.R', 'icom-nrep2.R', 
		'icom-rep-nrep1.R', 'icom-rep-nrep2.R', 'icom-nrep1-nrep2.R', 
		'icom-rep-nrep1-nrep2.R')
curr.model <- paste('code/simulations/nimble-code/', code.files[index[1]], sep = '')

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# MCMC Settings -----------------------------------------------------------
n.iter <- 20000
n.thin <- 4
n.burn <- 10000
n.chain <- 1
# For testing
# n.iter <- 5000
# n.thin <- 1
# n.burn <- 1000
# n.chain <- 1
# Number of iterations post-burn in and thinning
n.pb <- (n.iter - n.burn) / n.thin

# Simulation setup --------------------------------------------------------
# Number of simulations
n.sims <- 100
# Set seeds for consistency
set.seed(101)
my.seeds <- sample(1:101748, n.sims, replace = FALSE)
# Initialize list to store results
samples.list <- list()

# Specify parameters ------------------------------------------------------
# Number sites with replicated data source
J.rep <- 50
# Number sites with first nonreplicated data source
J.nrep.1 <- 50
# Number of sites with second nonreplicated data source
J.nrep.2 <- 50
# Total number of sites (data are collected at different locations in the
# same region).
J <- J.rep + J.nrep.1 + J.nrep.2
# Number of years
n.years <- 6
# Community level parameters ----------------------------------------------
# Occupancy intercept (vary by year)
beta.0.mean <- runif(n.years, -1.5, 1.5)
# Occupancy spatial covariate effect
beta.1.mean <- runif(1, -0.5, 0.5)
# Community auto-logistic mean
phi.mean <- runif(1, -0.2, 0.8)
# Replicated detection intercept (vary by year)
alpha.0.mean <- runif(n.years, -1, 1)
# Replicated detection space/time covariate effect
alpha.1.mean <- runif(1, -0.75, 0.75)
# Nonreplicated 1 detection intercept (vary by year)
gamma.1.0.mean <- runif(n.years, -2.5, 0)
# Nonreplicated 1 detection space/time covariate effect
gamma.1.1.mean <- runif(1, -0.75, 0.75)
# Nonreplicated 2 detection intercept (vary by year)
gamma.2.0.mean <- runif(n.years, 0, 2.5)
# Nonreplicated 2 detection space/time covariate effect
gamma.2.1.mean <- runif(1, -0.75, 0.75)
# Number of replicates for replicated data set 
K.rep <- 3
# Number of replicates for first nonreplicated data set (can change if 
# desired to have two replicated data sets)
K.nrep.1 <- 1
# Number of species
I <- 25
# Community-level variance parameters
sigma.sq.beta.0 <- runif(n.years, 0.5, 1.5)
sigma.sq.beta.1 <- runif(1, 0.25, 2)
sigma.sq.phi <- runif(1, 0.25, 1)
sigma.sq.det.0 <- runif(n.years, 1, 3)
sigma.sq.det.1 <- runif(1, 0.25, 2)
sigma.sq.nrep.1.0 <- runif(n.years, 1, 3)
sigma.sq.nrep.1.1 <- runif(1, 0.25, 2)
sigma.sq.nrep.2.0 <- runif(n.years, 1, 3)
sigma.sq.nrep.2.1 <- runif(1, 0.25, 2)

# Form species-specific covariates ----------------------------------------
# Parameter names are analogous to community level means. 
beta.0 <- matrix(NA, nrow = I, ncol = n.years)
alpha.0 <- matrix(NA, nrow = I, ncol = n.years)
gamma.1.0 <- matrix(NA, nrow = I, ncol = n.years)
gamma.2.0 <- matrix(NA, nrow = I, ncol = n.years)
for (t in 1:n.years) {
  beta.0[, t] <- rnorm(I, beta.0.mean[t], sqrt(sigma.sq.beta.0[t]))
  alpha.0[, t] <- rnorm(I, alpha.0.mean[t], sqrt(sigma.sq.det.0[t]))
  gamma.1.0[, t] <- rnorm(I, gamma.1.0.mean[t], sqrt(sigma.sq.nrep.1.0[t]))
  gamma.2.0[, t] <- rnorm(I, gamma.2.0.mean[t], sqrt(sigma.sq.nrep.2.0[t]))
}
beta.1 <- rnorm(I, beta.1.mean, sqrt(sigma.sq.beta.1)) 
phi <- rnorm(I, phi.mean, sqrt(sigma.sq.phi))
gamma.1.1 <- rnorm(I, gamma.1.1.mean, sqrt(sigma.sq.nrep.1.1))
alpha.1 <- rnorm(I, alpha.1.mean, sqrt(sigma.sq.det.1))
gamma.2.1 <- rnorm(I, gamma.2.1.mean, sqrt(sigma.sq.nrep.2.1))

# Run simulations ---------------------------------------------------------
# Simulate full data set ------------
# Set seeed. 
set.seed(my.seeds[1])
dat <- sim.icom.data(J.rep, J.nrep.1, J.nrep.2, beta.0, beta.1,
		    phi, gamma.2.0, gamma.2.1, alpha.0, alpha.1, 
		    gamma.1.0, gamma.1.1, K.rep, K.nrep.1, I, n.years)
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
# Load nimble code, constants, data, and initial values
source(curr.model)
# Create Model ------------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
  		          data = icom.data, inits = icom.inits)
  
# Configure MCMC ----------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('int.beta.mean', 'beta.1.mean',
						  'phi.mean', 
						  'beta.0', 'beta.1', 
						  'phi'))
# Create an MCMC function -------------------------------------------------
icom.mcmc <- buildMCMC(icom.conf)
# Compile model -----------------------------------------------------------
icom.c.model <- compileNimble(icom.model)
icom.c.mcmc <- compileNimble(icom.mcmc, project = icom.model)
# Run first data set for current condition ------------------------------
samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn, 
                 thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
		 setSeed = index[2])
# Subtract true values from samples
true.vals <- c(beta.0, beta.1, mean(beta.1),
               logit.inv(apply(beta.0, 2, mean)), phi, mean(phi))
true.vals <- rep(true.vals, each = n.pb)
true.vals <- matrix(true.vals, nrow = n.pb)
# Generate bias values
samples <- samples - true.vals
# Save bias values
samples.list[[1]] <- samples

for (j in 2:n.sims) {
  set.seed(my.seeds[j])
  print(paste('Currently on simulation ', j, ' out of ', n.sims, sep = ''))
  # Draw new parameter values ---------------------------------------------
  beta.0.mean <- runif(n.years, -1.5, 1.5)
  beta.1.mean <- runif(1, -0.5, 0.5)
  phi.mean <- runif(1, -0.2, 0.8)
  alpha.0.mean <- runif(n.years, -1, 1)
  alpha.1.mean <- runif(1, -0.75, 0.75)
  gamma.1.0.mean <- runif(n.years, -2.5, 0)
  gamma.1.1.mean <- runif(1, -0.75, 0.75)
  gamma.2.0.mean <- runif(n.years, 0, 2.5)
  gamma.2.1.mean <- runif(1, -0.75, 0.75)
  sigma.sq.beta.0 <- runif(n.years, 0.5, 1.5)
  sigma.sq.beta.1 <- runif(1, 0.25, 2)
  sigma.sq.phi <- runif(1, 0.25, 1)
  sigma.sq.det.0 <- runif(n.years, 1, 3)
  sigma.sq.det.1 <- runif(1, 0.25, 2)
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
    alpha.0[, t] <- rnorm(I, alpha.0.mean[t], sqrt(sigma.sq.det.0[t]))
    gamma.1.0[, t] <- rnorm(I, gamma.1.0.mean[t], sqrt(sigma.sq.nrep.1.0[t]))
    gamma.2.0[, t] <- rnorm(I, gamma.2.0.mean[t], sqrt(sigma.sq.nrep.2.0[t]))
  }
  beta.1 <- rnorm(I, beta.1.mean, sqrt(sigma.sq.beta.1)) 
  phi <- rnorm(I, phi.mean, sqrt(sigma.sq.phi))
  gamma.1.1 <- rnorm(I, gamma.1.1.mean, sqrt(sigma.sq.nrep.1.1))
  alpha.1 <- rnorm(I, alpha.1.mean, sqrt(sigma.sq.det.1))
  gamma.2.1 <- rnorm(I, gamma.2.1.mean, sqrt(sigma.sq.nrep.2.1))
  # Simulate a new dataset
  dat <- sim.icom.data(J.rep, J.nrep.1, J.nrep.2, beta.0, beta.1,
		      phi, gamma.2.0, gamma.2.1,
		      alpha.0, alpha.1, 
		      gamma.1.0, gamma.1.1, K.rep, K.nrep.1, I, 
		      n.years)
  # Load constants and data again to overwrite with new data values 
  source(curr.model)
  # Reset data and initial values -----
  icom.c.model$z <- array(1, dim = c(I, icom.consts$J, n.years))
  icom.c.model$X.psi <- dat$X.psi[psi.indices, ]
  # Current model contains replicated data
  if (index[1] %in% c(1, 4, 5, 7)) { 
  icom.c.model$y <- dat$y
  icom.c.model$X.rep <- dat$X.rep
  }
  # Current model contains nonreplicated 1 data
  if (index[1] %in% c(2, 4, 6, 7)) {
  icom.c.model$v.1 <- dat$v.1
  icom.c.model$X.nrep.1 <- dat$X.nrep.1
  }
  # Current model contains nonreplicated 2 data
  if (index[1] %in% c(3, 5, 6, 7)) {
    icom.c.model$v.2 <- dat$v.2
    icom.c.model$X.nrep.2 <- dat$X.nrep.2
  }
  # Fit the model -------------------------------------------------------
  samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn, 
    	       thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
	       setSeed = index[2])
  # Subtract true values from samples
  true.vals <- c(beta.0, beta.1, mean(beta.1),
                 logit.inv(apply(beta.0, 2, mean)), phi, mean(phi))
  true.vals <- rep(true.vals, each = n.pb)
  true.vals <- matrix(true.vals, nrow = n.pb)
  # Generate bias values
  samples <- samples - true.vals
  # Save bias values
  samples.list[[j]] <- samples
}
# Save the results
date <- Sys.Date()

file.name <- paste('results/sim-icom-results-1-', n.sims, '-condition-', index[1], '-simulations-', 
		   index[2], '-chain-', date, '.R', sep = '')

save(samples.list, file = file.name)


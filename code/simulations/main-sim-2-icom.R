# main-sim-2-icom.R: this script provides code to implement the second
#                   simulation study described in Doser et al (2021). 
#                   This study assesses how the integrated community 
#                   occupancy model performs compared to separate 
#                   integrated distribution models.  The code compares 
#                   model performance when using two data sources: a 
#                   replicated data source and a nonreplicated data source
#                   with medium detection probability. This script is used
#                   to generate results under the ICOM. See main-sim-2-idm.R
#                   for code that generates results under the IDM. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
set.seed(101)
library(nimble)
library(coda)
# Load function to simulate data
source("code/simulations/sim-icom-data.R")

# This index is used to run simulations across multiple cores. This enables
# running different chains on separate cores instead of running them all
# within the same script. The number supplied on the command line is the 
# seed that will be used. 
# For example: 
#    Rscript code/simulations/main-sim-2-icom.R 1 & 
index <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(index) == 0) base::stop("Need to indicate the seed for the current chain")
# Alternatively, can just specify the seed within the script itself
# index <- 1
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


# Simulation setup --------------------------------------------------------
# Initialize list to store results.
samples.list <- list()
# Number of simulations.
n.sims <- 100
# Set seed using user provided input.
set.seed(590)
my.seeds <- sample(1:101748, n.sims, replace = FALSE)

# Specify fixed parameters ------------------------------------------------
# NOTE: data are generated for three data sources (one replicated and 
#       two nonreplicated), but only data from one replicated and one
#       nonreplicated are used in the simulation process. 
# Number sites with replicated data.
J.rep <- 50
# Number sites with the first nonreplicated data.
J.nrep.1 <- 50
# Number of sites with the second nonreplicated data.
J.nrep.2 <- 50
# Total number of sites.
J <- J.rep + J.nrep.1 + J.nrep.2
# Number of years.
n.years <- 6
# Number of replicates for replicated data set.
K.rep <- 3
# Number of replicates for first nonreplicated data set
K.nrep.1 <- 1
# Number of species. 
I <- 25

# Form species-specific covariates ----------------------------------------
# Occurrence intercept
beta.0 <- matrix(NA, nrow = I, ncol = n.years)
# Replicated data detection intercept
alpha.0 <- matrix(NA, nrow = I, ncol = n.years)
# Nonreplicated 1 data detection intercept
gamma.1.0 <- matrix(NA, nrow = I, n.years)
# Nonreplicated 2 data detection intercept
gamma.2.0 <- matrix(NA, nrow = I, n.years)
# Generate species-specific covariates from uniform distributions
# so that species are not simulated as coming directly from the ICOM 
# framework. 
for (t in 1:n.years) {
  beta.0[, t] <- runif(I, -1.5, 1.5)
  gamma.1.0[, t] <- runif(I, -1, 1)
  alpha.0[, t] <- runif(I, -1, 1) 
  gamma.2.0[, t] <- runif(I, -1, 1)
}
# Occurrence covariate effect
beta.1 <- runif(I, -0.5, 0.5)
# Autologistic effect
phi <- runif(I, -0.5, 2.0)
# Replicated data detection covariate effect
gamma.1.1 <- runif(I, -0.5, 0.5) 
# Nonreplicated 1 data detection covariate effect
alpha.1 <- runif(I, -0.5, 0.5)
# Nonreplicated 2 data detection covariate effect
gamma.2.1 <- runif(I, -0.5, 0.5)

# Run simulations ---------------------------------------------------------
# Simulate full integrated data set with three data sets
dat <- sim.icom.data(J.rep, J.nrep.1, J.nrep.2, beta.0, beta.1,
		    phi, gamma.2.0, gamma.2.1,
		    alpha.0, alpha.1,
		    gamma.1.0, gamma.1.1, K.rep, K.nrep.1, I,
		    n.years)

# Load nimble code, constants, data, and initial values 
# NOTE: Here we load a model that uses only two data sources, a replicated
#       data source and then a single nonreplicated data source. 
source("code/simulations/nimble-code/icom-rep-nrep1.R")
# Create Model ------------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
  		          data = icom.data, inits = icom.inits)
  
# Configure MCMC ----------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('int.beta.mean', 'beta.1.mean',
						  'phi.mean', 'beta.0',
						  'beta.1', 'phi', 
						  'alpha.0.mean', 'alpha.1.mean',
						  'gamma.1.0.mean', 'gamma.1.1.mean', 
						  'alpha.0', 'alpha.1', 
						  'gamma.1.0', 'gamma.1.1'))
# Create an MCMC function -------------------------------------------------
icom.mcmc <- buildMCMC(icom.conf)
# Compile model -----------------------------------------------------------
icom.c.model <- compileNimble(icom.model)
icom.c.mcmc <- compileNimble(icom.mcmc, project = icom.model)
# Run first data set for current condition ------------------------------
samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn, 
                 thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
		 setSeed = index)
samples.list[[1]] <- samples

for (j in 2:n.sims) {
  set.seed(my.seeds[j])
  print(paste('Currently on simulation ', j, ' out of ', n.sims, sep = ''))
  # Simulate a new dataset
  dat <- sim.icom.data(J.rep, J.nrep.1, J.nrep.2, beta.0, beta.1,
		      phi, gamma.2.0, gamma.2.1,
		      alpha.0, alpha.1, 
		      gamma.1.0, gamma.1.1, K.rep, K.nrep.1, I, 
		      n.years)
  # Load nimble code, constants, data, and initial values 
  source("code/simulations/nimble-code/icom-rep-nrep1.R")
  # Reset data and initial values
  icom.c.model$z <- array(1, dim = c(I, icom.consts$J, n.years))
  icom.c.model$y <- dat$y
  icom.c.model$v.1 <- dat$v.1
  icom.c.model$X.psi <- dat$X.psi[psi.indices, ]
  icom.c.model$X.rep <- dat$X.rep
  icom.c.model$X.nrep.1 <- dat$X.nrep.1

  # Fit the model -------------------------------------------------------
  samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn, 
    	       thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
	       setSeed = index)
  samples.list[[j]] <- samples
}
# Save the results
date <- Sys.Date()

file.name <- paste('results/sim-icom-results-2-', n.sims, '-simulations-',
		  index, '-chain-', date, '.R', sep = '')

save(samples.list, beta.0, beta.1, phi,
     alpha.0, alpha.1, gamma.1.0, gamma.1.1, file = file.name)

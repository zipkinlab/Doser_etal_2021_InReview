# main-icm-nimble.R: R code to run ICM through NIMBLE via R. 
rm(list = ls())
library(nimble)
library(coda)
library(tidyverse)
# Load function to simulate data
source("sim-icm-data.R")
# This contains the BUGS code for the NIMBLE model. 
source("icm-nimble.R")
set.seed(1817)
# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Simulate Data -----------------------------------------------------------
# Total number of pixels.
J <- 100
# Number of cells with eBird data. The eBird cells are assumed to span the 
# entire region of the R pixels. 
J.eb <- 20
# Number pixels with DET data
J.det <- 100
# Number pixels with NEON data
J.neon <- 100
# Number routes
n.route <- 2
# Number of pixels occupied by BBS data. 
# Note, a route is assumed to be the horizontal length of pixels
J.bbs <- n.route * sqrt(J)
# Community level parameters (intercept, cov effect)
# Initial occupancy
beta.psi.mean <- c(0.5, 0.1)
# Persistence
beta.phi.mean <- c(1, 0.1)
# Colonization
beta.gamma.mean <- c(0, -0.1)
# eBird detection
alpha.eb.mean <- c(-1, 0.4)
# DET detection
alpha.det.mean <- c(0.4, 0.2)
# NEON detection
alpha.neon.mean <- c(0.7, 0.05)
# BBS detection
alpha.bbs.mean <- c(0.2, -0.1)
# Number of repeat visits for eBird data. This is 8 in real life, but
# with lots of missing data. 
K.eb <- 5
# Number of DET repeat visits
K.det <- 3
# Number of NEON repeat visits. This is 1 in real life. 
K.neon <- 1
# Number of species
I <- 10
# Number of years
n.years <- 10
# Variance parameters for all effects
sigma.sq.psi <- c(2, 0.3)
sigma.sq.phi <- c(1, 0.4)
sigma.sq.gamma <- c(1, 0.8)
sigma.sq.eb <- c(1, 0.2)
sigma.sq.det <- c(0.6, 0.2)
sigma.sq.neon <- c(0.6, 0.2)
sigma.sq.bbs <- c(0.6, 0.2)

# Form species-specific covariates ----------------------------------------
beta.psi.0 <- rnorm(I, beta.psi.mean[1], sqrt(sigma.sq.psi[1]))
beta.psi.1 <- rnorm(I, beta.psi.mean[2], sqrt(sigma.sq.psi[2])) 
beta.phi.0 <- rnorm(I, beta.phi.mean[1], sqrt(sigma.sq.phi[1]))
beta.phi.1 <- rnorm(I, beta.phi.mean[2], sqrt(sigma.sq.phi[2])) 
beta.gamma.0 <- rnorm(I, beta.gamma.mean[1], sqrt(sigma.sq.gamma[1]))
beta.gamma.1 <- rnorm(I, beta.gamma.mean[2], sqrt(sigma.sq.gamma[2])) 
alpha.eb.0 <- rnorm(I, alpha.eb.mean[1], sqrt(sigma.sq.eb[1]))
alpha.eb.1 <- rnorm(I, alpha.eb.mean[2], sqrt(sigma.sq.eb[2]))
alpha.neon.0 <- rnorm(I, alpha.neon.mean[1], sqrt(sigma.sq.neon[1]))
alpha.neon.1 <- rnorm(I, alpha.neon.mean[2], sqrt(sigma.sq.neon[2]))
alpha.det.0 <- rnorm(I, alpha.det.mean[1], sqrt(sigma.sq.det[1]))
alpha.det.1 <- rnorm(I, alpha.det.mean[2], sqrt(sigma.sq.det[2]))
alpha.bbs.0 <- rnorm(I, alpha.bbs.mean[1], sqrt(sigma.sq.bbs[1]))
alpha.bbs.1 <- rnorm(I, alpha.bbs.mean[2], sqrt(sigma.sq.bbs[2]))

# Explore these parameter values ------
logit.inv(beta.psi.0)
logit.inv(beta.phi.0)
logit.inv(beta.gamma.0)
logit.inv(alpha.eb.0)
logit.inv(alpha.det.0)
logit.inv(alpha.neon.0)
logit.inv(alpha.bbs.0)

# Simulate the data using the data simulation function. 
# Variables that come from this are described in sim-icm-data.R
dat <- sim.icm.data(J, J.eb, J.det, J.neon, J.bbs, beta.psi.0, beta.psi.1, 
		    beta.phi.0, beta.phi.1, beta.gamma.0, beta.gamma.1, 
		    alpha.bbs.0, alpha.bbs.1,
		    alpha.eb.0, alpha.eb.1, alpha.det.0, alpha.det.1, 
		    alpha.neon.0, alpha.neon.1, K.eb, K.det, K.neon, I, 
		    n.years, n.route, seed = 10187)

# Constants ---------------------------------------------------------------
icm.consts <- list(I = I, J = J, X.psi = dat$X.psi, n.years = n.years, 
		  X.phi = dat$X.phi, X.gamma = dat$X.gamma, K.det = K.det, 
		  J.det = J.det, X.det = dat$X.det, 
                  pixel.det = dat$pixel.det, K.eb = K.eb, low = dat$low, 
		  high = dat$high, J.eb = J.eb, X.eb = dat$X.eb,  
		  K.neon = K.neon, J.neon = J.neon, X.neon = dat$X.neon, 
		  pixel.neon = dat$pixel.neon,  
		  pixel.bbs = dat$pixel.bbs, X.bbs = dat$X.bbs, J.bbs = J.bbs
		  )

# Data --------------------------------------------------------------------
icm.data <- list(C = dat$C, x = dat$x, y = dat$y, y.bbs = dat$y.bbs)

# Initial values ----------------------------------------------------------
z.init <- array(1, dim = c(I, J, n.years))
#z.init <- dat$z
icm.inits <- list(z = z.init, int.psi.0.mean = runif(1, 0.1, 0.9), 
		  beta.psi.1.mean = rnorm(1), int.phi.0.mean = runif(1, 0.1, 0.9), 
		  beta.phi.1.mean = rnorm(1), int.gamma.0.mean = runif(1, 0.1, 0.9), 
		  beta.gamma.1.mean = rnorm(1), 
		  int.alpha.eb.0.mean = runif(1, 0.1, 0.9), 
		  alpha.eb.1.mean = rnorm(1), 
		  int.alpha.det.0.mean = runif(1, 0.1, 0.9), 
		  alpha.det.1.mean = rnorm(1), 
		  int.alpha.neon.0.mean = runif(1, 0.1, 0.9), 
		  alpha.neon.1.mean = rnorm(1), 
		  int.alpha.bbs.0.mean = runif(1, 0.1, 0.9), 
		  alpha.bbs.1.mean = rnorm(1), tau.beta.psi.0 = runif(1, 0.1, 2), 
		  tau.beta.psi.1 = runif(1, 0.1, 2), tau.beta.phi.0 = runif(1, 0.1, 2), 
		  tau.beta.phi.1 = runif(1, 0.1, 2), 
		  tau.beta.gamma.0 = runif(1, 0.1, 2), 
		  tau.beta.gamma.1 = runif(1, 0.1, 2), 
		  tau.alpha.eb.0 = runif(1, 0.1, 2),
		  tau.alpha.eb.1 = runif(1, 0.1, 2), tau.alpha.det.0 = runif(1, 0.1, 2),
		  tau.alpha.det.1 = runif(1, 0.1, 2), 
		  tau.alpha.neon.0 = runif(1, 0.1, 2), 
		  tau.alpha.neon.1 = runif(1, 0.1, 2), 
		  tau.alpha.bbs.0 = runif(1, 0.1, 2), 
		  tau.alpha.bbs.1 = runif(1, 0.1, 2))
start <- Sys.time()
# Create Model ------------------------------------------------------------
icm.model <- nimbleModel(code = icm.code, name = 'icm', constants = icm.consts,
  		          data = icm.data, inits = icm.inits)
  
# Configure MCMC ----------------------------------------------------------
icm.conf <- configureMCMC(icm.model, monitors = c('int.psi.0.mean', 'beta.psi.1.mean', 
				         'int.phi.0.mean', 'beta.phi.1.mean', 
					 'beta.gamma.1.mean', 'int.gamma.0.mean', 
					 'alpha.bbs.1.mean', 'beta.psi.0', 'beta.psi.1', 
					 'beta.phi.0', 'beta.phi.1', 'beta.gamma.1', 
					 'beta.gamma.0', 'alpha.bbs.1'))
# Create an MCMC function -------------------------------------------------
icm.mcmc <- buildMCMC(icm.conf)
# Compile model -----------------------------------------------------------
icm.c.model <- compileNimble(icm.model)
icm.c.mcmc <- compileNimble(icm.mcmc, project = icm.model)
# Number of iterations ----------------------------------------------------
n.iter <- 1000
n.thin <- 1
n.burn <- 500
n.chain <- 1
# Run the model -----------------------------------------------------------
samples <- runMCMC(icm.c.mcmc, niter = n.iter, nburnin = n.burn,
   	       thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)
end.time <- Sys.time()
run.time <- end.time - start
run.time
summary(samples)
plot(mcmc(samples), density = FALSE)

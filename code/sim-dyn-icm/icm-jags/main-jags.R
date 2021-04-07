# main-jags.R: R code to run ICM in jags via jagsUI. 
# Author: Jeffrey W. Doser

rm(list = ls())
library(jagsUI)
library(coda)
library(tidyverse)
# Load function to simulate data
source("../sim-icm-data.R")
set.seed(1817)
# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Simulate Data -----------------------------------------------------------
# Total number of pixels.
# This should be a perfect square
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

# Bundle data -------------------------------------------------------------
bugs.data <- list(I = I, J = J, X.psi = dat$X.psi, n.years = n.years, 
		  X.phi = dat$X.phi, X.gamma = dat$X.gamma, J.det = J.det, 
		  K.det = K.det, X.det = dat$X.det, C = dat$C, 
                  pixel.det = dat$pixel.det, J.eb = J.eb, low = dat$low, 
		  high = dat$high, K.eb = K.eb, X.eb = dat$X.eb, y = dat$y, 
		  J.neon = J.neon, K.neon = K.neon, X.neon = dat$X.neon, 
		  x = dat$x, pixel.neon = dat$pixel.neon, y.bbs = dat$y.bbs, 
		  pixel.bbs = dat$pixel.bbs, J.bbs = J.bbs, X.bbs = dat$X.bbs
		  )

# Initial values ----------------------------------------------------------
z.init <- array(1, dim = c(I, J, n.years))
#z.init <- dat$z
inits <- function() {
  list(
    z = z.init, 
    beta.psi.1.mean = rnorm(1), 
    beta.phi.1.mean = rnorm(1), 
    beta.gamma.1.mean = rnorm(1), 
    alpha.eb.1.mean = rnorm(1),
    alpha.det.1.mean = rnorm(1),
    alpha.neon.1.mean = rnorm(1)
  )
}

# Parameters monitored ----------------------------------------------------
params <- c('beta.psi.0.mean', 'beta.psi.1.mean', 'alpha.eb.0.mean', 
	    'beta.phi.0.mean', 'beta.phi.1.mean', 'beta.gamma.0.mean', 
	    'beta.gamma.1.mean',
	    'alpha.eb.1.mean', 'alpha.det.0.mean', 'alpha.det.1.mean', 
	    'alpha.neon.0.mean', 'alpha.neon.1.mean', 
	    'alpha.bbs.0.mean', 'alpha.bbs.1.mean', 
	    'beta.psi.0', 'beta.psi.1', 'beta.phi.0', 'beta.phi.1', 
	    'beta.gamma.0', 'beta.gamma.1')

# MCMC settings -----------------------------------------------------------
n.iter <- 500
n.thin <- 1
n.burn <- 100
n.chain <- 1

out <- jags(bugs.data, inits, params, 'icm-jags.txt', 
	    n.iter = n.iter, n.thin = n.thin, 
	    n.burn = n.burn, n.chain = n.chain, parallel = FALSE)

# Assessment of Model Performance -----------------------------------------
# Initial Occupancy Probability -------
# Intercept
beta.psi.0.samples <- out$sims.list$beta.psi.0
low.beta.psi.0 <- apply(beta.psi.0.samples, 2, quantile, prob = 0.025)
med.beta.psi.0 <- apply(beta.psi.0.samples, 2, median)
high.beta.psi.0 <- apply(beta.psi.0.samples, 2, quantile, prob = 0.975)

ggplot(mapping = aes(x = beta.psi.0, y = med.beta.psi.0)) + 
  geom_point() + 
  theme_bw() + 
  geom_segment(aes(x = beta.psi.0, y = low.beta.psi.0, xend = beta.psi.0, 
		   yend = high.beta.psi.0), col = 'grey', lineend = 'round') + 
  geom_line(aes(x = beta.psi.0, y = beta.psi.0), col = 'grey', lty = 2)
# Covariate 
beta.psi.1.samples <- out$sims.list$beta.psi.1
low.beta.psi.1 <- apply(beta.psi.1.samples, 2, quantile, prob = 0.025)
med.beta.psi.1 <- apply(beta.psi.1.samples, 2, median)
high.beta.psi.1 <- apply(beta.psi.1.samples, 2, quantile, prob = 0.975)

ggplot(mapping = aes(x = beta.psi.1, y = med.beta.psi.1)) + 
  geom_point() + 
  theme_bw() + 
  geom_segment(aes(x = beta.psi.1, y = low.beta.psi.1, xend = beta.psi.1, 
		   yend = high.beta.psi.1), col = 'grey', lineend = 'round') + 
  geom_line(aes(x = beta.psi.1, y = beta.psi.1), col = 'grey', lty = 2)
# Persistence -------------------------
beta.phi.0.samples <- out$sims.list$beta.phi.0
low.beta.phi.0 <- apply(beta.phi.0.samples, 2, quantile, prob = 0.025)
med.beta.phi.0 <- apply(beta.phi.0.samples, 2, median)
high.beta.phi.0 <- apply(beta.phi.0.samples, 2, quantile, prob = 0.975)

ggplot(mapping = aes(x = beta.phi.0, y = med.beta.phi.0)) + 
  geom_point() + 
  theme_bw() + 
  geom_segment(aes(x = beta.phi.0, y = low.beta.phi.0, xend = beta.phi.0, 
		   yend = high.beta.phi.0), col = 'grey', lineend = 'round') + 
  geom_line(aes(x = beta.phi.0, y = beta.phi.0), col = 'grey', lty = 2)
# Covariate 
beta.phi.1.samples <- out$sims.list$beta.phi.1
low.beta.phi.1 <- apply(beta.phi.1.samples, 2, quantile, prob = 0.025)
med.beta.phi.1 <- apply(beta.phi.1.samples, 2, median)
high.beta.phi.1 <- apply(beta.phi.1.samples, 2, quantile, prob = 0.975)

ggplot(mapping = aes(x = beta.phi.1, y = med.beta.phi.1)) + 
  geom_point() + 
  theme_bw() + 
  geom_segment(aes(x = beta.phi.1, y = low.beta.phi.1, xend = beta.phi.1, 
		   yend = high.beta.phi.1), col = 'grey', lineend = 'round') + 
  geom_line(aes(x = beta.phi.1, y = beta.phi.1), col = 'grey', lty = 2)

# Colonization -------------------------
beta.gamma.0.samples <- out$sims.list$beta.gamma.0
low.beta.gamma.0 <- apply(beta.gamma.0.samples, 2, quantile, prob = 0.025)
med.beta.gamma.0 <- apply(beta.gamma.0.samples, 2, median)
high.beta.gamma.0 <- apply(beta.gamma.0.samples, 2, quantile, prob = 0.975)

ggplot(mapping = aes(x = beta.gamma.0, y = med.beta.gamma.0)) + 
  geom_point() + 
  theme_bw() + 
  geom_segment(aes(x = beta.gamma.0, y = low.beta.gamma.0, xend = beta.gamma.0, 
		   yend = high.beta.gamma.0), col = 'grey', lineend = 'round') + 
  geom_line(aes(x = beta.gamma.0, y = beta.gamma.0), col = 'grey', lty = 2)
# Covariate 
beta.gamma.1.samples <- out$sims.list$beta.gamma.1
low.beta.gamma.1 <- apply(beta.gamma.1.samples, 2, quantile, prob = 0.025)
med.beta.gamma.1 <- apply(beta.gamma.1.samples, 2, median)
high.beta.gamma.1 <- apply(beta.gamma.1.samples, 2, quantile, prob = 0.975)

ggplot(mapping = aes(x = beta.gamma.1, y = med.beta.gamma.1)) + 
  geom_point() + 
  theme_bw() + 
  geom_segment(aes(x = beta.gamma.1, y = low.beta.gamma.1, xend = beta.gamma.1, 
		   yend = high.beta.gamma.1), col = 'grey', lineend = 'round') + 
  geom_line(aes(x = beta.gamma.1, y = beta.gamma.1), col = 'grey', lty = 2)

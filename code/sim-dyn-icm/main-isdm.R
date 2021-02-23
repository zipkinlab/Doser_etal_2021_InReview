rm(list = ls())
library(jagsUI)
library(coda)
library(tidyverse)
# Load function to simulate data
source("sim-icm-data.R")
#set.seed(10107)
# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Simulate Data -----------------------------------------------------------
# Total number of pixels
R <- 100
# Number of cells with eBird data. The eBird cells are assumed to span the 
# entire region of the R pixels. 
R.eb <- 20
# Number pixels with HB data
R.hb <- 100
# Number pixels with NEON data
R.neon <- 100
# Community level parameters (intercept, cov effect)
# Initial occupancy
beta.psi.mean <- c(0.5, 0.1)
# Persistence
beta.phi.mean <- c(0, 0.1)
# Colonization
beta.gamma.mean <- c(0, -0.1)
# eBird detection
alpha.eb.mean <- c(-1, 0.4)
# HB detection
alpha.hb.mean <- c(0, 0.2)
# NEON detection
alpha.neon.mean <- c(0, 0.05)
# Number of repeat visits for eBird data. This is 8 in real life, but
# with lots of missing data. 
J.eb <- 5
# Number of HB repeat visits
J.hb <- 3
# Number of NEON repeat visits. This is 1 in real life. 
J.neon <- 3
# Number of species
K <- 10
# Number of years
n.years <- 10
# Variance parameters for all effects
sigma.sq.psi <- c(2, 0.3)
sigma.sq.phi <- c(1, 0.4)
sigma.sq.gamma <- c(1, 0.8)
sigma.sq.eb <- c(1, 0.2)
sigma.sq.hb <- c(0.6, 0.2)
sigma.sq.neon <- c(0.6, 0.2)

# Form species-specific covariates ----------------------------------------
beta.psi.0 <- rnorm(K, beta.psi.mean[1], sqrt(sigma.sq.psi[1]))
beta.psi.1 <- rnorm(K, beta.psi.mean[2], sqrt(sigma.sq.psi[2])) 
beta.phi.0 <- rnorm(K, beta.phi.mean[1], sqrt(sigma.sq.phi[1]))
beta.phi.1 <- rnorm(K, beta.phi.mean[2], sqrt(sigma.sq.phi[2])) 
beta.gamma.0 <- rnorm(K, beta.gamma.mean[1], sqrt(sigma.sq.gamma[1]))
beta.gamma.1 <- rnorm(K, beta.gamma.mean[2], sqrt(sigma.sq.gamma[2])) 
alpha.eb.0 <- rnorm(K, alpha.eb.mean[1], sqrt(sigma.sq.eb[1]))
alpha.eb.1 <- rnorm(K, alpha.eb.mean[2], sqrt(sigma.sq.eb[2]))
alpha.neon.0 <- rnorm(K, alpha.neon.mean[1], sqrt(sigma.sq.neon[1]))
alpha.neon.1 <- rnorm(K, alpha.neon.mean[2], sqrt(sigma.sq.neon[2]))
alpha.hb.0 <- rnorm(K, alpha.hb.mean[1], sqrt(sigma.sq.hb[1]))
alpha.hb.1 <- rnorm(K, alpha.hb.mean[2], sqrt(sigma.sq.hb[2]))

# Explore these parameter values ------
logit.inv(beta.psi.0)
logit.inv(beta.phi.0)
logit.inv(beta.gamma.0)
logit.inv(alpha.eb.0)
logit.inv(alpha.hb.0)
logit.inv(alpha.neon.0)

dat <- sim.icm.data(R, R.eb, R.hb, R.neon, beta.psi.0, beta.psi.1, 
		    beta.phi.0, beta.phi.1, beta.gamma.0, beta.gamma.1, 
		    alpha.eb.0, alpha.eb.1, alpha.hb.0, alpha.hb.1, 
		    alpha.neon.0, alpha.neon.1, J.eb, J.hb, J.neon, K, 
		    n.years)
# Filter data to only grab a single species. 
curr.sp <- 5
logit.inv(beta.psi.0[curr.sp])
logit.inv(alpha.eb.0[curr.sp])
logit.inv(alpha.hb.0[curr.sp])
logit.inv(alpha.neon.0[curr.sp])
table(dat$z[, curr.sp, ])
# Impute missing data for testing. 
#prop.missing <- 0
#n.vals <- ceiling(prop.missing * length(dat$y))
#curr.indices <- sample(length(dat$y), n.vals, replace = FALSE)
#dat$y[curr.indices] <- NA
# Bundle data -------------------------------------------------------------
bugs.data <- list(R = R, X.psi = dat$X.psi, n.years = n.years, 
		  X.phi = dat$X.phi, X.gamma = dat$X.gamma, R.hb = R.hb, 
		  J.hb = J.hb, X.hb = dat$X.hb, C = dat$C[, , curr.sp, ], 
                  pixel.hb = dat$pixel.hb, R.eb = R.eb, low = dat$low, 
		  high = dat$high, J.eb = J.eb, X.eb = dat$X.eb, y = dat$y[, , curr.sp, ], 
		  R.neon = R.neon, J.neon = J.neon, X.neon = dat$X.neon, 
		  x = dat$x[, , curr.sp, ], pixel.neon = dat$pixel.neon
		  )

# Initial values ----------------------------------------------------------
z.init <- array(1, dim = c(R, n.years))
inits <- function() {
  list(
    z = z.init, 
    beta.psi.1 = rnorm(1), 
    beta.phi.1 = rnorm(1), 
    beta.gamma.1 = rnorm(1), 
    alpha.eb.1 = rnorm(1),
    alpha.hb.1 = rnorm(1),
    alpha.neon.1 = rnorm(1)
  )
}

# Parameters monitored ----------------------------------------------------
params <- c('beta.psi.0', 'beta.psi.1', 'beta.phi.0', 'beta.phi.1', 
	    'beta.gamma.0', 'beta.gamma.1')

# MCMC settings -----------------------------------------------------------
n.iter <- 3000
n.thin <- 1
n.burn <- 1000
n.chain <- 1

out <- jags(bugs.data, inits, params, 'isdm-jags.txt', 
	    n.iter = n.iter, n.thin = n.thin, 
	    n.burn = n.burn, n.chain = n.chain, parallel = FALSE)

# Look at estimated values compared to truth. 
out
beta.psi.0[curr.sp]
beta.psi.1[curr.sp]
beta.phi.0[curr.sp]
beta.phi.1[curr.sp]
beta.gamma.0[curr.sp]
beta.gamma.1[curr.sp]

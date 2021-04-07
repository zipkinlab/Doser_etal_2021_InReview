# main-sim-3-icm.R: this script provides code to implement the third
#                   simulation study described in Doser et al (2021). 
#                   This study assess how incorporating different amounts
#                   of the four data sets in the study area of interest
#                   changes the underlying accuracy and precision of 
#                   occupancy dynamics estimates. 
rm(list = ls())
set.seed(101)
library(nimble)
library(coda)
# Load function to simulate data
source("sim-icm-data.R")
# This contains the BUGS code for the NIMBLE model
source("icm-nimble.R")
# This index is used to run simulations across multiple cores. Each 
# simulation scenario is repeated n.sims times on its own core. Numbers
# supplied to index when running should be between 1 and 15. 
index <- as.numeric(commandArgs(trailingOnly = TRUE))
# Alternatively, can specify which simulation scenario to run manually.
#index <- 1
# To keep track of time
start.time <- Sys.time()

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# MCMC Settings -----------------------------------------------------------
# n.iter <- 10000
# n.thin <- 1
# n.burn <- 5000
# n.chain <- 1
n.iter <- 500
n.thin <- 1
n.burn <- 1
n.chain <- 1


# Simulation setup --------------------------------------------------------
# Number of simulations
n.sims <- 100
set.seed(380)
my.seeds <- sample(1:101748, n.sims, replace = FALSE)
# Simulated number of sites with NEON data
J.neon.vals <- c(10, 30, 50)
# Simulated number of sites with DET data
J.det.vals <- c(10, 30, 50)
# Simulated number of BBS routes
n.route.vals <- c(1, 3, 5)
# Simulated number of eBird cells
J.eb.vals <- c(25)
n.scenarios <- length(J.neon.vals) * length(J.eb.vals) * length(n.route.vals) * 
  length(J.det.vals)
param.vals <- expand.grid(J.neon.vals, J.eb.vals, n.route.vals, J.det.vals)
names(param.vals) <- c('J.neon', 'J.eb', 'n.route', 'J.det')

samples.list <- list()

# Specify fixed parameters ------------------------------------------------
# See main-icm.R for definitions. 
J <- 100
beta.psi.mean <- c(0.5, 0.1)
beta.phi.mean <- c(1, 0.1)
beta.gamma.mean <- c(0, -0.1)
alpha.eb.mean <- c(-1, 0.4)
alpha.det.mean <- c(0.4, 0.2)
alpha.neon.mean <- c(0.7, 0.05)
alpha.bbs.mean <- c(0.2, -0.1)
K.eb <- 8
K.det <- 3
K.neon <- 1
I <- 9
n.years <- 5
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

# Run simulations ---------------------------------------------------------
# Simulate full data set ------------
curr.vals <- param.vals[index, ]
J.det <- curr.vals$J.det
J.neon <- curr.vals$J.neon
J.eb <- curr.vals$J.eb
n.route <- curr.vals$n.route
J.bbs <- sqrt(J) * n.route
dat <- sim.icm.data(J, J.eb, J.det, J.neon, J.bbs, beta.psi.0, beta.psi.1,
		    beta.phi.0, beta.phi.1, beta.gamma.0, beta.gamma.1,
		    alpha.bbs.0, alpha.bbs.1,
		    alpha.eb.0, alpha.eb.1, alpha.det.0, alpha.det.1,
		    alpha.neon.0, alpha.neon.1, K.eb, K.det, K.neon, I,
		    n.years, n.route)

# Impute missing eBird data to make a more realistic data set. 
prop.missing <- 0.5
# Year/visit/site combos that could be missing. 
combos <- expand.grid(1:K.eb, 1:J.eb, 1:n.years)
n.vals <- ceiling(prop.missing * nrow(combos))
curr.indices <- sample(nrow(combos), n.vals, replace = FALSE)
combos.missing <- combos[curr.indices, ]
for (q in 1:n.vals) {
  dat$y[, combos.missing$Var2[q], combos.missing$Var1[q], combos.missing$Var3[q]] <- NA
}
# Constants -----------------------------------------------------------
icm.consts <- list(I = I, J = J, n.years = n.years, J.det = J.det, 
    	     K.det = K.det, J.eb = J.eb, low = dat$low, high = dat$high, 
    	     K.eb = K.eb, J.neon = J.neon, K.neon = K.neon, 
    	     pixel.det = dat$pixel.det, pixel.neon = dat$pixel.neon, 
    	     pixel.bbs = dat$pixel.bbs, J.bbs = J.bbs)
# Data --------------------------------------------------------------------
icm.data <- list(C = dat$C, x = dat$x, y = dat$y, y.bbs = dat$y.bbs, 
    	   X.phi = dat$X.phi, X.psi = dat$X.psi, X.gamma = dat$X.gamma, 
    	   X.det = dat$X.det, X.eb = dat$X.eb, X.neon = dat$X.neon, 
    	   X.bbs = dat$X.bbs)

# Initial values ----------------------------------------------------------
z.init <- array(1, dim = c(I, J, n.years))
icm.inits <- list(z = z.init, int.psi.0.mean = runif(1, 0.1, 0.9), 
		      beta.psi.1.mean = rnorm(1), int.phi.0.mean = runif(1, 0.1, 0.9), 
		      beta.phi.1.mean = rnorm(1), 
    	              int.gamma.0.mean = runif(1, 0.1, 0.9), 
		      beta.gamma.1.mean = rnorm(1), 
		      int.alpha.eb.0.mean = runif(1, 0.1, 0.9), 
		      alpha.eb.1.mean = rnorm(1), 
		      int.alpha.det.0.mean = runif(1, 0.1, 0.9), 
		      alpha.det.1.mean = rnorm(1), 
		      int.alpha.neon.0.mean = runif(1, 0.1, 0.9), 
		      alpha.neon.1.mean = rnorm(1), 
		      int.alpha.bbs.0.mean = runif(1, 0.1, 0.9), 
		      alpha.bbs.1.mean = rnorm(1), tau.beta.psi.0 = runif(1, 0.1, 2), 
		      tau.beta.psi.1 = runif(1, 0.1, 2), 
    	              tau.beta.phi.0 = runif(1, 0.1, 2), 
		      tau.beta.phi.1 = runif(1, 0.1, 2), 
		      tau.beta.gamma.0 = runif(1, 0.1, 2), 
		      tau.beta.gamma.1 = runif(1, 0.1, 2), 
		      tau.alpha.eb.0 = runif(1, 0.1, 2),
		      tau.alpha.eb.1 = runif(1, 0.1, 2), 
    	              tau.alpha.det.0 = runif(1, 0.1, 2),
		      tau.alpha.det.1 = runif(1, 0.1, 2), 
		      tau.alpha.neon.0 = runif(1, 0.1, 2), 
		      tau.alpha.neon.1 = runif(1, 0.1, 2), 
		      tau.alpha.bbs.0 = runif(1, 0.1, 2), 
		      tau.alpha.bbs.1 = runif(1, 0.1, 2))
# Create Model ------------------------------------------------------------
icm.model <- nimbleModel(code = icm.code, name = 'icm', constants = icm.consts,
  		          data = icm.data, inits = icm.inits)
  
# Configure MCMC ----------------------------------------------------------
icm.conf <- configureMCMC(icm.model, monitors = c('int.psi.0.mean', 'beta.psi.1.mean',
                                         'int.phi.0.mean', 'beta.phi.1.mean',
                                         'beta.gamma.1.mean', 'int.gamma.0.mean',
                                         'beta.psi.0', 'beta.psi.1',
                                         'beta.phi.0', 'beta.phi.1', 'beta.gamma.1',
                                         'beta.gamma.0'))
# Create an MCMC function -------------------------------------------------
icm.mcmc <- buildMCMC(icm.conf)
# Compile model -----------------------------------------------------------
icm.c.model <- compileNimble(icm.model)
icm.c.mcmc <- compileNimble(icm.mcmc, project = icm.model)
# Run first data set for current condition ------------------------------
samples <- runMCMC(icm.c.mcmc, niter = n.iter, nburnin = n.burn, 
                 thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)
samples.list[[1]] <- samples

for (j in 2:n.sims) {
  set.seed(my.seeds[j])
  print(paste('Currently on simulation ', j, ' out of ', n.sims, sep = ''))
  # Simulate a new dataset with same pixel values (but different covariates)
  dat <- sim.icm.data(J, J.eb, J.det, J.neon, J.bbs, beta.psi.0, beta.psi.1,
   	        beta.phi.0, beta.phi.1, beta.gamma.0, beta.gamma.1,
    	        alpha.bbs.0, alpha.bbs.1,
    	        alpha.eb.0, alpha.eb.1, alpha.det.0, alpha.det.1,
    	        alpha.neon.0, alpha.neon.1, K.eb, K.det, K.neon, I,
    	        n.years, n.route, dat$pixel.det, dat$pixel.bbs,
    		dat$pixel.neon)
  for (q in 1:n.vals) {
    dat$y[, combos.missing$Var2[q], combos.missing$Var1[q], combos.missing$Var3[q]] <- NA
  }
  icm.c.model$y <- dat$y
  icm.c.model$C <- dat$C
  icm.c.model$x <- dat$x
  icm.c.model$y.bbs <- dat$y.bbs
  icm.c.model$X.phi <- dat$X.phi
  icm.c.model$X.psi <- dat$X.psi
  icm.c.model$X.gamma <- dat$X.gamma
  icm.c.model$X.det <- dat$X.det
  icm.c.model$X.eb <- dat$X.eb
  icm.c.model$X.neon <- dat$X.neon
  icm.c.model$X.bbs <- dat$X.bbs

  # Fit the model -------------------------------------------------------
  samples <- runMCMC(icm.c.mcmc, niter = n.iter, nburnin = n.burn, 
    	       thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)
  samples.list[[j]] <- samples
}
end.time <- Sys.time()
time.taken <- end.time - start.time
# Save the results
date <- Sys.Date()

file.name <- paste('results/sim-icm-results-3-', n.sims, '-condition-', index, '-simulations-', date, '.R', sep = '')

save(samples.list, time.taken, beta.psi.0, beta.psi.1,
     beta.phi.0, beta.phi.1, beta.gamma.0, beta.gamma.1, file = file.name)


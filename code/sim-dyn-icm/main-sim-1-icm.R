# main-sim-1-icm.R: this script provides code to implement the first 
#                   simulation study described in Doser et al (2021). 
#                   This study assess general model performance of the 
#                   integrated community model across a wide range of 
#                   parameters, and compares model performance to models
#                   using smaller amounts of the four data sources. 

rm(list = ls())
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
n.iter <- 10000
n.thin <- 1
n.burn <- 5000
n.chain <- 1
# n.iter <- 500
# n.thin <- 1
# n.burn <- 1
# n.chain <- 1
n.pb <- (n.iter - n.burn) / n.thin

# Simulation setup --------------------------------------------------------
# Number of simulations
n.sims <- 10
set.seed(101)
my.seeds <- sample(1:101748, n.sims, replace = FALSE)
# Is eBird included in the model? 
eb.in.vals <- c(FALSE, TRUE)
# Is DET included in the model? 
det.in.vals <- c(FALSE, TRUE)
# Is NEON included in the model? 
neon.in.vals <- c(FALSE, TRUE)
# Is BBS included in the model? 
bbs.in.vals <- c(FALSE, TRUE)
# Total number of simulation scenarios
# The minus one eliminates the case when no data sets are available
n.scenarios <- length(eb.in.vals) * length(det.in.vals) * length(neon.in.vals) * 
  length(bbs.in.vals) - 1
# Data frame containing different conditions
param.vals <- expand.grid(eb.in.vals, det.in.vals, neon.in.vals, bbs.in.vals)
names(param.vals) <- c('eb.in', 'det.in', 'neon.in', 'bbs.in')
param.vals <- param.vals[-1, ]

# Initialize list to store results
samples.list <- list()

# Specify parameters ------------------------------------------------------
# See main-icm.R for definitions. 
J <- 100
J.eb <- 25
J.det <- 50
J.neon <- 50
n.route <- 5
J.bbs <- sqrt(J) * n.route
beta.psi.mean <- c(runif(1, -1, 2), runif(1, -0.5, 0.5))
beta.phi.mean <- c(runif(1, 0, 2), runif(1, -0.5, 0.5))
beta.gamma.mean <- c(runif(1, -3, 2), runif(1, -0.5, 0.5))
alpha.eb.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
alpha.det.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
alpha.neon.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
alpha.bbs.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
K.eb <- 8
K.det <- 3
n.years <- 5
K.neon <- 1
I <- 9
sigma.sq.psi <- c(runif(1, 1, 3), runif(1, 0.25, 2))
sigma.sq.phi <- c(runif(1, 1, 3), runif(1, 0.25, 2))
sigma.sq.gamma <- c(runif(1, 1, 3), runif(1, 0.25, 2))
sigma.sq.eb <- c(runif(1, 1, 3), runif(1, 0.25, 2))
sigma.sq.det <- c(runif(1, 1, 3), runif(1, 0.25, 2))
sigma.sq.neon <- c(runif(1, 1, 3), runif(1, 0.25, 2))
sigma.sq.bbs <- c(runif(1, 1, 3), runif(1, 0.25, 2))

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
eb.in <- curr.vals$eb.in
det.in <- curr.vals$det.in
neon.in <- curr.vals$neon.in
bbs.in <- curr.vals$bbs.in
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
# Set values to NA if those data are not included
if (!eb.in) dat$y <- array(NA, dim = dim(dat$y))
if (!det.in) dat$C <- array(NA, dim = dim(dat$C))
if (!neon.in) dat$x <- array(NA, dim = dim(dat$x))
if (!bbs.in) dat$y.bbs <- array(NA, dim = dim(dat$y.bbs))

# Constants -----------------------------------------------------------
icm.consts <- list(I = I, J = J, n.years = n.years, K.det = K.det, 
    	     J.det = J.det, K.eb = K.eb, low = dat$low, high = dat$high, 
    	     J.eb = J.eb, K.neon = K.neon, J.neon = J.neon, 
    	     pixel.det = dat$pixel.det, pixel.neon = dat$pixel.neon, 
    	     pixel.bbs = dat$pixel.bbs, J.bbs = J.bbs)
# Data --------------------------------------------------------------------
icm.data <- list(C = dat$C, x = dat$x, y = dat$y, y.bbs = dat$y.bbs, 
    	   X.phi = dat$X.phi, X.psi = dat$X.psi, X.gamma = dat$X.gamma, 
    	   X.det = dat$X.det, X.eb = dat$X.eb, X.neon = dat$X.neon, 
    	   X.bbs = dat$X.bbs)

# Initial values ----------------------------------------------------------
z.init <- array(1, dim = c(I, J, n.years))
icm.inits <- list(z = z.init, 
		  int.psi.0.mean = logit.inv(beta.psi.mean[1] + runif(1, -0.5, 0.5)), 
		  beta.psi.1.mean = beta.psi.mean[2] + runif(1, -0.2, 0.2), 
		  int.phi.0.mean = logit.inv(beta.phi.mean[1] + runif(1, -0.5, 0.5)), 
		  beta.phi.1.mean = beta.phi.mean[2] + runif(1, -0.2, 0.2), 
		  int.gamma.0.mean = logit.inv(beta.gamma.mean[1] + runif(1, -0.5, 0.5)), 
		  beta.gamma.1.mean = beta.gamma.mean[2] + runif(1, -0.2, 0.2), 
		  int.alpha.eb.0.mean = logit.inv(alpha.eb.mean[1] + runif(1, -0.5, 0.5)), 
		  alpha.eb.1.mean = alpha.eb.mean[2] + runif(1, -0.2, 0.2),
		  int.alpha.det.0.mean = logit.inv(alpha.det.mean[1] + runif(1, -0.5, 0.5)), 
		  alpha.det.1.mean = alpha.det.mean[2] + runif(1, -0.2, 0.2),
		  int.alpha.neon.0.mean = logit.inv(alpha.neon.mean[1] + runif(1, -0.5, 0.5)), 
		  alpha.neon.1.mean = alpha.neon.mean[2] + runif(1, -0.2, 0.2),
		  int.alpha.bbs.0.mean = logit.inv(alpha.bbs.mean[1] + runif(1, -0.5, 0.5)), 
		  alpha.bbs.1.mean = alpha.bbs.mean[2] + runif(1, -0.2, 0.2),
		  tau.beta.psi.0 = 1 / sigma.sq.psi[1] + runif(1, -0.1, 0.1), 
		  tau.beta.psi.1 = 1 / sigma.sq.psi[2] + runif(1, -0.1, 0.1), 
		  tau.beta.phi.0 = 1 / sigma.sq.phi[1] + runif(1, -0.1, 0.1), 
		  tau.beta.phi.1 = 1 / sigma.sq.phi[2] + runif(1, -0.1, 0.1), 
		  tau.beta.gamma.0 = 1 / sigma.sq.gamma[1] + runif(1, -0.1, 0.1), 
		  tau.beta.gamma.1 = 1 / sigma.sq.gamma[2] + runif(1, -0.1, 0.1), 
		  tau.alpha.eb.0 = 1 / sigma.sq.eb[1] + runif(1, -0.1, 0.1), 
		  tau.alpha.eb.1 = 1 / sigma.sq.eb[2] + runif(1, -0.1, 0.1), 
		  tau.alpha.det.0 = 1 / sigma.sq.det[1] + runif(1, -0.1, 0.1), 
		  tau.alpha.det.1 = 1 / sigma.sq.det[2] + runif(1, -0.1, 0.1), 
		  tau.alpha.neon.0 = 1 / sigma.sq.neon[1] + runif(1, -0.1, 0.1), 
		  tau.alpha.neon.1 = 1 / sigma.sq.neon[2] + runif(1, -0.1, 0.1), 
		  tau.alpha.bbs.0 = 1 / sigma.sq.bbs[1] + runif(1, -0.1, 0.1), 
		  tau.alpha.bbs.1 = 1 / sigma.sq.bbs[2] + runif(1, -0.1, 0.1))
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
# Subtract true values from samples
true.vals <- c(beta.gamma.0, beta.gamma.1, beta.gamma.mean[2], 
	       beta.phi.0, beta.phi.1, beta.phi.mean[2], 
	       beta.psi.0, beta.psi.1, beta.psi.mean[2], 
	       logit.inv(beta.gamma.mean[1]), logit.inv(beta.phi.mean[1]), 
	       logit.inv(beta.psi.mean[1]))
true.vals <- rep(true.vals, each = n.pb)
true.vals <- matrix(true.vals, nrow = n.pb)
samples <- samples - true.vals
samples.list[[1]] <- samples

for (j in 2:n.sims) {
  set.seed(my.seeds[j])
  print(paste('Currently on simulation ', j, ' out of ', n.sims, sep = ''))
  # Draw new parameter values ---------------------------------------------
  beta.psi.mean <- c(runif(1, -1, 2), runif(1, -0.5, 0.5))
  beta.phi.mean <- c(runif(1, 0, 2), runif(1, -0.5, 0.5))
  beta.gamma.mean <- c(runif(1, -3, 2), runif(1, -0.5, 0.5))
  alpha.eb.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
  alpha.det.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
  alpha.neon.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
  alpha.bbs.mean <- c(runif(1, -1, 1), runif(1, -0.75, 0.75))
  sigma.sq.psi <- c(runif(1, 1, 3), runif(1, 0.25, 2))
  sigma.sq.phi <- c(runif(1, 1, 3), runif(1, 0.25, 2))
  sigma.sq.gamma <- c(runif(1, 1, 3), runif(1, 0.25, 2))
  sigma.sq.eb <- c(runif(1, 1, 3), runif(1, 0.25, 2))
  sigma.sq.det <- c(runif(1, 1, 3), runif(1, 0.25, 2))
  sigma.sq.neon <- c(runif(1, 1, 3), runif(1, 0.25, 2))
  sigma.sq.bbs <- c(runif(1, 1, 3), runif(1, 0.25, 2))
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
  if (!eb.in) dat$y <- array(NA, dim = dim(dat$y))
  if (!det.in) dat$C <- array(NA, dim = dim(dat$C))
  if (!neon.in) dat$x <- array(NA, dim = dim(dat$x))
  if (!bbs.in) dat$y.bbs <- array(NA, dim = dim(dat$y.bbs))
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
  # Subtract true values from samples
  true.vals <- c(beta.gamma.0, beta.gamma.1, beta.gamma.mean[2], 
  	       beta.phi.0, beta.phi.1, beta.phi.mean[2], 
  	       beta.psi.0, beta.psi.1, beta.psi.mean[2], 
  	       logit.inv(beta.gamma.mean[1]), logit.inv(beta.phi.mean[1]), 
  	       logit.inv(beta.psi.mean[1]))
  true.vals <- rep(true.vals, each = n.pb)
  true.vals <- matrix(true.vals, nrow = n.pb)
  samples <- samples - true.vals
  samples.list[[j]] <- samples
}
end.time <- Sys.time()
time.taken <- end.time - start.time
# Save the results
date <- Sys.Date()

file.name <- paste('results/sim-icm-results-1-', n.sims, '-condition-', index, '-simulations-', date, '.R', sep = '')

save(samples.list, file = file.name)


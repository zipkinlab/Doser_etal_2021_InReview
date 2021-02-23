rm(list = ls())
library(jagsUI)
# Load function to simulate data
source("sim-icm-data.R")
set.seed(101)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# MCMC Settings -----------------------------------------------------------
n.iter <- 5000
n.thin <- 1
n.burn <- 3000
n.chain <- 3

# Simulation setup --------------------------------------------------------
# Number of simulations
n.sims <- 20
seeds <- sample(1:81010, n.sims)
# Vary the number of years you have data
n.years.vals <- c(5, 10)
# Vary the number of repeat visits for HB data
J.hb.vals <- c(1, 3)
# Vary the number of repeat visits for neon data
J.neon.vals <- c(1, 3)
# Total number of simulation scenarios in a factorial design
n.scenarios <- length(n.years.vals) * length(J.hb.vals) * length(J.neon.vals)
param.vals <- expand.grid(n.years.vals, J.hb.vals, J.neon.vals)
names(param.vals) <- c('n.years', 'J.hb', 'J.neon')

out.model <- list()

# Specify fixed parameters ------------------------------------------------
# See main-icm.R for definitions.
R <- 100
R.eb <- 25
R.hb <- 25
R.neon <- 25
beta.psi.mean <- c(0.5, 0.1)
beta.phi.mean <- c(0, 0.1)
beta.gamma.mean <- c(0, -0.1)
alpha.eb.mean <- c(-1, 0.4)
alpha.hb.mean <- c(0, 0.2)
alpha.neon.mean <- c(0, 0.05)
J.eb <- 5
K <- 10
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

# Run simulations ---------------------------------------------------------
for (i in 1:n.sims) {
  print(paste("Current iteration: ", i, " out of ", n.sims, sep = ''))
  set.seed(seeds[i])
  for (j in 1:n.scenarios) {
    print(paste("Current scenario: ", j, ' out of ', n.scenarios, sep = ''))
    n.years <- param.vals$n.years[j]
    J.hb <- param.vals$J.hb[j]
    J.neon <- param.vals$J.neon[j]
    dat <- sim.icm.data(R, R.eb, R.hb, R.neon, beta.psi.0, beta.psi.1, 
  		        beta.phi.0, beta.phi.1, beta.gamma.0, beta.gamma.1, 
  		        alpha.eb.0, alpha.eb.1, alpha.hb.0, alpha.hb.1, 
  		        alpha.neon.0, alpha.neon.1, J.eb, J.hb, J.neon, K, 
  		        n.years)
    bugs.data <- list(K = K, R = R, X.psi = dat$X.psi, n.years = n.years, 
		      X.phi = dat$X.phi, X.gamma = dat$X.gamma, R.hb = R.hb, 
		      J.hb = J.hb, X.hb = dat$X.hb, C = dat$C, 
                      pixel.hb = dat$pixel.hb, R.eb = R.eb, low = dat$low, 
		      high = dat$high, J.eb = J.eb, X.eb = dat$X.eb, y = dat$y, 
		      R.neon = R.neon, J.neon = J.neon, X.neon = dat$X.neon, 
		      x = dat$x, pixel.neon = dat$pixel.neon
		 )

    # Initial values ----------------------------------------------------------
    z.init <- array(1, dim = c(R, K, n.years))
    inits <- function() {
      list(
        z = z.init, 
        beta.psi.0.mean = rnorm(1), 
        beta.psi.1.mean = rnorm(1), 
        beta.phi.0.mean = rnorm(1), 
        beta.phi.1.mean = rnorm(1), 
        beta.gamma.0.mean = rnorm(1), 
        beta.gamma.1.mean = rnorm(1), 
        alpha.eb.0.mean = rnorm(1), 
        alpha.eb.1.mean = rnorm(1),
        alpha.hb.0.mean = rnorm(1), 
        alpha.hb.1.mean = rnorm(1),
        alpha.neon.0.mean = rnorm(1), 
        alpha.neon.1.mean = rnorm(1)
      )
    }

    # Parameters monitored ----------------------------------------------------
    params <- c('beta.psi.0', 'beta.psi.1', 'beta.phi.0', 'beta.phi.1', 
		'beta.gamma.0', 'beta.gamma.1')
    
    # Fit the model -------------------------------------------------------
    out <- jags(bugs.data, inits, params, 'icm-jags.txt',
	        n.iter = n.iter, n.thin = n.thin,
	        n.burn = n.burn, n.chain = n.chain, parallel = TRUE, verbose = FALSE)
    out.model[[(i-1) * n.scenarios + j]] <- out

  } # j
  # Read out smaller versions of the simulations since they take awhile
  if (i %in% c(5, 10)) {
    date <- Sys.Date()
    file.name <- paste('results/sim-icm-results-2-', i, '-simulations-', date, '.R', sep = '')
    save(out.model, beta.psi.0, beta.psi.1, beta.phi.0, 
	 beta.phi.1, beta.gamma.0, beta.gamma.1, param.vals, file = file.name)
  }
} # i

# Save the results
date <- Sys.Date()
file.name <- paste('results/sim-icm-results-2-', i, '-simulations-', date, '.R', sep = '')
save(out.model, beta.psi.0, beta.psi.1, beta.phi.0, 
     beta.phi.1, beta.gamma.0, beta.gamma.1, param.vals, file = file.name)


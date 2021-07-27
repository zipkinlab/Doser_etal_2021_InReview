# icom-nrep1-nrep2.R: BUGS code for running ICOM with two nonreplicated
#                     data sources. BUGS code is for use with NIMBLE. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 
require(nimble)
icom.code <- nimbleCode ({
  # Priors ----------------------------------------------------------------
  # Time-varying community intercepts and intercept precision parameters
  for (t in 1:n.years) {
    # Occurrence ----------------------
    int.beta.mean[t] ~ dunif(0, 1) 
    beta.0.mean[t] <- logit(int.beta.mean[t])
    tau.beta.0[t] ~ dgamma(0.1, 0.1)
    # Detection -----------------------
    int.gamma.1.0.mean[t] ~ dunif(0, 1)
    gamma.1.0.mean[t] <- logit(int.gamma.1.0.mean[t])
    int.gamma.2.0.mean[t] ~ dunif(0, 1)
    gamma.2.0.mean[t] <- logit(int.gamma.2.0.mean[t])
    tau.gamma.2.0[t] ~ dgamma(0.1, 0.1)
    tau.gamma.1.0[t] ~ dgamma(0.1, 0.1)
  }
  # Community level covariate effects -
  beta.1.mean ~ dnorm(0, 0.1)
  phi.mean ~ dnorm(0, 0.1) # auto-logistic parameter
  # Nonreplicated data 1
  gamma.1.1.mean ~ dnorm(0, 0.1)
  # Nonreplicated data 2
  gamma.2.1.mean ~ dnorm(0, 0.1)
  tau.beta.1 ~ dgamma(0.1, 0.1)
  tau.phi ~ dgamma(0.1, 0.1)
  tau.gamma.1.1 ~ dgamma(0.1, 0.1)
  tau.gamma.2.1 ~ dgamma(0.1, 0.1)

  # Species-specific coefficients -----------------------------------------
  for (i in 1:I) {
    # Intercepts ------------------------
    for (t in 1:n.years) {
      beta.0[i, t] ~ dnorm(beta.0.mean[t], tau.beta.0[t])
      gamma.1.0[i, t] ~ dnorm(gamma.1.0.mean[t], tau.gamma.1.0[t])
      gamma.2.0[i, t] ~ dnorm(gamma.2.0.mean[t], tau.gamma.2.0[t])
    }
    # Covariate Effects ---------------
    beta.1[i] ~ dnorm(beta.1.mean, tau.beta.1)
    phi[i] ~ dnorm(phi.mean, tau.phi)
    gamma.1.1[i] ~ dnorm(gamma.1.1.mean, tau.gamma.1.1)
    gamma.2.1[i] ~ dnorm(gamma.2.1.mean, tau.gamma.2.1)
  }

  # Likelihood and Process Models -----------------------------------------
  for (i in 1:I) {
    # Process Model -------------------------------------------------------
    for (j in 1:J) {
      logit(psi[i, j, 1]) <- beta.0[i, 1] +
	  	             beta.1[i] * X.psi[j, 2] 
      z[i, j, 1] ~ dbern(psi[i, j, 1])
      for (t in 2:n.years) {
        logit(psi[i, j, t]) <- beta.0[i, t] + 
		               beta.1[i] * X.psi[j, 2] + 
			       phi[i] * z[i, j, t - 1]
        z[i, j, t] ~ dbern(psi[i, j, t]) 
      } # t
    } # j

    for (t in 1:n.years) {
      # Likelihoods -------------------------------------------------------
      # Nonreplicated 1 data ----------
      for (j in 1:J.nrep.1) {
        for (k in 1:K.nrep.1) {
          logit(pi.1[i, j, k, t]) <- gamma.1.0[i, t] + 
                                         gamma.1.1[i] * X.nrep.1[j, k, t, 2]
          v.1[i, j, k, t] ~ dbern(pi.1[i, j, k, t] * z[i, site.nrep.1[j], t])
        } # k
      } # j

      # Nonreplicated 2 data ----------
      for (j in 1:J.nrep.2) {
        logit(pi.2[i, j, t]) <- gamma.2.0[i, t] + 
                                gamma.2.1[i] * X.nrep.2[j, t, 2]
        v.2[i, j, t] ~ dbern(pi.2[i, j, t] * z[i, site.nrep.2[j], t])
      } # j
    } # t
  } # i
})

# Constants -----------------------------------------------------------
icom.consts <- list(I = I, J = J.nrep.1 + J.nrep.2, n.years = n.years, 
		   K.nrep.1 = K.nrep.1, J.nrep.1 = J.nrep.1,
    	           site.nrep.1 = 1:J.nrep.1,
    	           site.nrep.2 = 1:J.nrep.2 + J.nrep.1, J.nrep.2 = J.nrep.2)

# Data --------------------------------------------------------------------
psi.indices <- c(1:J.nrep.1 + J.rep, 1:J.nrep.2 + J.nrep.1 + J.rep)
icom.data <- list(v.1 = dat$v.1, v.2 = dat$v.2,
    	         X.psi = dat$X.psi[psi.indices, ], 
		 X.nrep.1 = dat$X.nrep.1, X.nrep.2 = dat$X.nrep.2)

# Initial values -----------------------------e-----------------------------
icom.inits <- list(z = array(1, dim = c(I, J.nrep.1 + J.nrep.2, n.years)),
		  int.beta.mean = runif(n.years, 0.1, 0.9),
		  beta.1.mean = rnorm(1), 
		  phi.mean = rnorm(1), 
		  int.gamma.1.0.mean = runif(n.years, 0.1, 0.9),
		  gamma.1.1.mean = rnorm(1),
		  int.gamma.2.0.mean = runif(n.years, 0.1, 0.9),
		  gamma.2.1.mean = rnorm(1), tau.beta.0 = runif(n.years, 0.1, 2),
		  tau.beta.1 = runif(1, 0.1, 2),
		  tau.phi = runif(1, 0.1, 2),
		  tau.gamma.1.0 = runif(n.years, 0.1, 2),
		  tau.gamma.1.1 = runif(1, 0.1, 2),
		  tau.gamma.2.0 = runif(n.years, 0.1, 2),
		  tau.gamma.2.1 = runif(1, 0.1, 2))


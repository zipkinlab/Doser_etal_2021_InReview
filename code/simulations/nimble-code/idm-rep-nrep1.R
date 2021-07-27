# idm-rep-nrep1.R: BUGS code, initial values, constants, and data used to run 
#                  integrated species distribution model using one replicated
#                  data source and one nonreplicated data source. BUGS code is 
#                  for use with NIMBLE
# Author: Jeffrey W. Doser (doserjef@msu.ed)
# Citation: 
require(nimble)
idm.code <- nimbleCode ({
  # Priors ----------------------------------------------------------------
  # Time-varying intercepts
  for (t in 1:n.years) {
    # Occurrence ----------------------
    int.beta[t] ~ dunif(0, 1) 
    beta.0[t] <- logit(int.beta[t])
    # Detection -----------------------
    int.alpha.0[t] ~ dunif(0, 1) 
    alpha.0[t] <- logit(int.alpha.0[t]) 
    int.gamma.1.0[t] ~ dunif(0, 1) 
    gamma.1.0[t] <- logit(int.gamma.1.0[t]) 
  }
  # Covariate effects -----------------
  beta.1 ~ dnorm(0, 0.1)
  phi ~ dnorm(0, 0.1) # auto-logistic parameter
  # Replicated Data
  alpha.1 ~ dnorm(0, 0.1)
  # Nonreplicated Data 1
  gamma.1.1 ~ dnorm(0, 0.1)

  # Likelihood and Process Models -----------------------------------------
  # Process Model -------------------------------------------------------
  for (j in 1:J) {
    logit(psi[j, 1]) <- beta.0[1] +
        	        beta.1 * X.psi[j, 2]
    z[j, 1] ~ dbern(psi[j, 1])
    for (t in 2:n.years) {
      logit(psi[j, t]) <- beta.0[t] + 
	                  beta.1 * X.psi[j, 2] + 
			  phi * z[j, t - 1]
      z[j, t] ~ dbern(psi[j, t]) 
    } # t
  } # j

  for (t in 1:n.years) {
    # Likelihoods -------------------------------------------------------
    # Replicated data -----------------
    for (j in 1:J.rep) {
      for (k in 1:K.rep) {
        logit(p[j, k, t]) <- alpha.0[t] + alpha.1 * X.rep[j, k, t, 2]
        y[j, k, t] ~ dbern(p[j, k, t] * z[site.rep[j], t])
      } # k 
    } # j 

    # Nonreplicated data --------------
    for (j in 1:J.nrep.1) {
      logit(pi.1[j, t]) <- gamma.1.0[t] + 
                           gamma.1.1 * X.nrep.1[j, t, 2]
      v.1[j, t] ~ dbern(pi.1[j, t] * z[site.nrep.1[j], t])
    } # j
  } # t

})

# Constants -----------------------------------------------------------
idm.consts <- list(J = J, n.years = n.years, K.rep = K.rep,
    	     J.rep = J.rep, K.nrep.1 = K.nrep.1, J.nrep.1 = J.nrep.1,
    	     site.rep = 1:J.rep, site.nrep.1 = 1:J.nrep.1 + J.rep)

# Data --------------------------------------------------------------------
idm.data <- list(y = y, v.1 = v.1,
    	         X.psi = X.psi, X.rep = X.rep,
		 X.nrep.1 = X.nrep.1)

# Initial values ----------------------------------------------------------
idm.inits <- list(z = array(1, dim = c(J, n.years)),
		  int.beta = runif(n.years, 0.1, 0.9),
		  beta.1 = rnorm(1), 
		  phi = rnorm(1), 
		  int.alpha.0 = runif(n.years, 0.1, 0.9),
		  alpha.1 = rnorm(1),
		  int.gamma.1.0 = runif(n.years, 0.1, 0.9),
		  gamma.1.1 = rnorm(1))

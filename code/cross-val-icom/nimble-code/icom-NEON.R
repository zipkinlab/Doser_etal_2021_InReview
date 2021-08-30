# icom-NEON.R: BUGS code to run community model in NIMBLE using only data
#              from NEON for the foliage-gleaning bird case study. No Bayesian
#              p-value included, which speeds things up.
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

require(nimble)
icom.code <- nimbleCode({
  # Priors ------------------------------------------------------------------
  # Occurrence ------------------------
  beta.1.mean ~ dnorm(0, 0.368) # linear elevation
  beta.2.mean ~ dnorm(0, 0.368) # quadratic elevation
  beta.3.mean ~ dnorm(0, 0.368) # percent forest
  phi.mean ~ dnorm(0, 0.368) # auto-logistic parameter
  for (t in 1:n.years) {
    int.beta.mean[t] ~ dunif(0, 1)
    beta.0.mean[t] <- logit(int.beta.mean[t])
    tau.beta.0[t] ~ dgamma(0.1, 0.1)
  } # t
  # Detection -------------------------
  int.gamma.1.mean ~ dunif(0, 1) # overall (species and yer) NEON detection
  gamma.1.0.mean <- logit(int.gamma.1.mean)
  gamma.1.1.mean ~ dnorm(0, 0.368) # day
  gamma.1.2.mean ~ dnorm(0, 0.368) # day^2
  gamma.1.3.mean ~ dnorm(0, 0.368) # hour
  # Precision Parameters --------------
  tau.beta.1 ~ dgamma(0.1, 0.1)
  tau.beta.2 ~ dgamma(0.1, 0.1)
  tau.beta.3 ~ dgamma(0.1, 0.1)
  tau.phi ~ dgamma(0.1, 0.1)
  tau.gamma.1.0 ~ dgamma(0.1, 0.1)
  tau.gamma.1.1 ~ dgamma(0.1, 0.1)
  tau.gamma.1.2 ~ dgamma(0.1, 0.1)
  tau.gamma.1.3 ~ dgamma(0.1, 0.1)

  # Species-specific coefficients -----------------------------------------
  for (i in 1:I) {
    beta.1[i] ~ dnorm(beta.1.mean, tau.beta.1)
    beta.2[i] ~ dnorm(beta.2.mean, tau.beta.2)
    beta.3[i] ~ dnorm(beta.3.mean, tau.beta.3)
    phi[i] ~ dnorm(phi.mean, tau.phi)
    for (t in 1:n.years) {
      gamma.1.0[i, t] ~ dnorm(gamma.1.0.mean, tau.gamma.1.0)
      beta.0[i, t] ~ dnorm(beta.0.mean[t], tau.beta.0[t])
      logit(int.beta[i, t]) <- beta.0[i, t]
    } # t
    gamma.1.1[i] ~ dnorm(gamma.1.1.mean, tau.gamma.1.1)
    gamma.1.2[i] ~ dnorm(gamma.1.2.mean, tau.gamma.1.2)
    gamma.1.3[i] ~ dnorm(gamma.1.3.mean, tau.gamma.1.3)
  } # i

  # Process Models and Likelihoods ----------------------------------------
  for (j in 1:J) {
    for (i in 1:I) {
      logit(psi[i, j, 1]) <- beta.0[i, 1] + 
	                     beta.1[i] * ELEV[j] + 
			     beta.2[i] * pow(ELEV[j], 2) + 
			     beta.3[i] * FOREST[j] 
      z.neon[i, j, 1] ~ dbern(psi[i, j, 1])
      for (t in 2:n.years) {	    
        # Process Model ---------------------------------------------------
        logit(psi[i, j, t]) <- beta.0[i, t] + 
                               beta.1[i] * ELEV[j] + 
          		       beta.2[i] * pow(ELEV[j], 2) + 
		               beta.3[i] * FOREST[j] + 
			       phi[i] * z.neon[i, j, t - 1]
        z.neon[i, j, t] ~ dbern(psi[i, j, t])
      } # t
    } # i
  } # j

  # NEON Data -----------------------------------------------------------
  # Data are stacked in a single vector as opposed to a multi-dimensional 
  # array to improve computational performance. 
  for (i in 1:n.vals.neon) {
    logit(pi.1[i]) <- gamma.1.0[sp.indx.neon[i], year.indx.neon[i]] + 
	                gamma.1.1[sp.indx.neon[i]] * DAY.neon[i] + 
			gamma.1.2[sp.indx.neon[i]] * pow(DAY.neon[i], 2) + 
			gamma.1.3[sp.indx.neon[i]] * HOUR.neon[i]
    v.1[i] ~ dbern(pi.1[i] * z.neon[sp.indx.neon[i], site.neon[i], year.indx.neon[i]])
    DAY.neon[i] ~ dnorm(0, 1)
    HOUR.neon[i] ~ dnorm(0, 1)
  } # i
})

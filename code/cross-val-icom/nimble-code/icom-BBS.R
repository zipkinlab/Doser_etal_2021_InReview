# icom-BBS.R: BUGS code to run community model in NIMBLE using only data
#             from BBS for the foliage-gleaning bird case study. 
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
  int.gamma.2.mean ~ dunif(0, 1) # overall (species and year) BBS detection
  gamma.2.0.mean <- logit(int.gamma.2.mean)
  gamma.2.1.mean ~ dnorm(0, 0.368) # day
  gamma.2.2.mean ~ dnorm(0, 0.368) # day^2
  # Precision parameters --------------
  tau.beta.1 ~ dgamma(0.1, 0.1)
  tau.beta.2 ~ dgamma(0.1, 0.1)
  tau.beta.3 ~ dgamma(0.1, 0.1)
  tau.phi ~ dgamma(0.1, 0.1)
  tau.gamma.2.0 ~ dgamma(0.1, 0.1)
  tau.gamma.2.1 ~ dgamma(0.1, 0.1)
  tau.gamma.2.2 ~ dgamma(0.1, 0.1)
  tau.gamma.2.3 ~ dgamma(0.1, 0.1) # Random observer effect
  for (i in 1:n.obsv.bbs) {
    gamma.2.3[i] ~ dnorm(0, tau.gamma.2.3)
  } # i

  # Species-specific coefficients -----------------------------------------
  for (i in 1:I) {
    beta.1[i] ~ dnorm(beta.1.mean, tau.beta.1)
    beta.2[i] ~ dnorm(beta.2.mean, tau.beta.2)
    beta.3[i] ~ dnorm(beta.3.mean, tau.beta.3)
    phi[i] ~ dnorm(phi.mean, tau.phi)
    for (t in 1:n.years) {
      beta.0[i, t] ~ dnorm(beta.0.mean[t], tau.beta.0[t])
      logit(int.beta[i, t]) <- beta.0[i, t]
      gamma.2.0[i, t] ~ dnorm(gamma.2.0.mean, tau.gamma.2.0)
    } # t
    gamma.2.1[i] ~ dnorm(gamma.2.1.mean, tau.gamma.2.1)
    gamma.2.2[i] ~ dnorm(gamma.2.2.mean, tau.gamma.2.2)
  } # i

  # Process Models and Likelihoods ----------------------------------------
  for (j in 1:J) {
    for (i in 1:I) {
      logit(psi[i, j, 1]) <- beta.0[i, 1] + 
	                     beta.1[i] * ELEV[j] + 
			     beta.2[i] * pow(ELEV[j], 2) + 
			     beta.3[i] * FOREST[j] 
      z.bbs[i, j, 1] ~ dbern(psi[i, j, 1])
      for (t in 2:n.years) {	    
        # Process Model ---------------------------------------------------
        logit(psi[i, j, t]) <- beta.0[i, t] + 
                               beta.1[i] * ELEV[j] + 
          		       beta.2[i] * pow(ELEV[j], 2) + 
		               beta.3[i] * FOREST[j] + 
			       phi[i] * z.bbs[i, j, t - 1]
        z.bbs[i, j, t] ~ dbern(psi[i, j, t])
      } # t
    } # i
  } # j

  # BBS Data --------------------------------------------------------------
  for (i in 1:n.vals.bbs) {
    logit(pi.2[i]) <- gamma.2.0[sp.indx.bbs[i], year.indx.bbs[i]] +
	              gamma.2.1[sp.indx.bbs[i]] * DAY.bbs[i] +
		      gamma.2.2[sp.indx.bbs[i]] * pow(DAY.bbs[i], 2) +
		      gamma.2.3[obsv.bbs[i]]
    v.2[i] ~ dbern(pi.2[i] * z.bbs[sp.indx.bbs[i], site.bbs[i], year.indx.bbs[i]])
  } # i
})



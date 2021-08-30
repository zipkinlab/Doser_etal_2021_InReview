# icom-NEON-BBS.R: BUGS code to run ICOm using NEON and BBS data for the 
#                  foliage-gleaning bird case study. 
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
  int.gamma.2.mean ~ dunif(0, 1) # overall (species and year) BBS detection
  gamma.2.0.mean <- logit(int.gamma.2.mean)
  gamma.2.1.mean ~ dnorm(0, 0.368) # day
  gamma.2.2.mean ~ dnorm(0, 0.368) # day^2
  # Precision parameters --------------
  tau.beta.1 ~ dgamma(0.1, 0.1)
  tau.beta.2 ~ dgamma(0.1, 0.1)
  tau.beta.3 ~ dgamma(0.1, 0.1)
  tau.phi ~ dgamma(0.1, 0.1)
  tau.gamma.1.0 ~ dgamma(0.1, 0.1)
  tau.gamma.1.1 ~ dgamma(0.1, 0.1)
  tau.gamma.1.2 ~ dgamma(0.1, 0.1)
  tau.gamma.1.3 ~ dgamma(0.1, 0.1)
  tau.gamma.2.0 ~ dgamma(0.1, 0.1)
  tau.gamma.2.1 ~ dgamma(0.1, 0.1)
  tau.gamma.2.2 ~ dgamma(0.1, 0.1)
  tau.gamma.2.3 ~ dgamma(0.1, 0.1) # Random observer effect
  for (i in 1:n.obsv.bbs) {
    gamma.2.3[i] ~ dnorm(0, tau.gamma.2.3)
  } # i

  # Species-specific coefficients -----------------------------------------
  # Process ---------------------------
  for (i in 1:I) {
    for (t in 1:n.years) {	
      beta.0[i, t] ~ dnorm(beta.0.mean[t], tau.beta.0[t])
      logit(int.beta[i, t]) <- beta.0[i, t]
    } # t
    beta.1[i] ~ dnorm(beta.1.mean, tau.beta.1)
    beta.2[i] ~ dnorm(beta.2.mean, tau.beta.2)
    beta.3[i] ~ dnorm(beta.3.mean, tau.beta.3)
    phi[i] ~ dnorm(phi.mean, tau.phi)
  } # i
  # NEON ------------------------------
  for (i in 1:I.neon) {
    for (t in 1:n.years.neon) {
      gamma.1.0[i, t] ~ dnorm(gamma.1.0.mean, tau.gamma.1.0)
    } # t
    gamma.1.1[i] ~ dnorm(gamma.1.1.mean, tau.gamma.1.1)
    gamma.1.2[i] ~ dnorm(gamma.1.2.mean, tau.gamma.1.2)
    gamma.1.3[i] ~ dnorm(gamma.1.3.mean, tau.gamma.1.3)
  } # i
  # BBS -------------------------------
  for (i in 1:I.bbs) {
    for (t in 1:n.years) {
      gamma.2.0[i, t] ~ dnorm(gamma.2.0.mean, tau.gamma.2.0)
    } # t
    gamma.2.1[i] ~ dnorm(gamma.2.1.mean, tau.gamma.2.1)
    gamma.2.2[i] ~ dnorm(gamma.2.2.mean, tau.gamma.2.2)
  } # i

  # Process Models and Likelihoods ----------------------------------------
  # NEON ------------------------------
  for (j in 1:J.neon) {
    for (i in 1:I.neon) {
      logit(psi.neon[i, j, 1]) <- beta.0[sp.neon[i], years.neon[1]] + 
	                          beta.1[sp.neon[i]] * ELEV.neon[j] + 
			          beta.2[sp.neon[i]] * pow(ELEV.neon[j], 2) + 
			          beta.3[sp.neon[i]] * FOREST.neon[j] 
      z.neon[i, j, 1] ~ dbern(psi.neon[i, j, 1])
      for (t in 2:n.years.neon) {	    
        # Process Model ---------------------------------------------------
        logit(psi.neon[i, j, t]) <- beta.0[sp.neon[i], years.neon[t]] + 
                                    beta.1[sp.neon[i]] * ELEV.neon[j] + 
          		            beta.2[sp.neon[i]] * pow(ELEV.neon[j], 2) + 
		                    beta.3[sp.neon[i]] * FOREST.neon[j] + 
			            phi[sp.neon[i]] * z.neon[i, j, t - 1]
        z.neon[i, j, t] ~ dbern(psi.neon[i, j, t])
      } # t
    } # i
  } # j
  # BBS -------------------------------
  for (j in 1:J.bbs) {
    for (i in 1:I) {
      logit(psi.bbs[i, j, 1]) <- beta.0[i, 1] + 
	                         beta.1[i] * ELEV.bbs[j] + 
			         beta.2[i] * pow(ELEV.bbs[j], 2) + 
			         beta.3[i] * FOREST.bbs[j] 
      z.bbs[i, j, 1] ~ dbern(psi.bbs[i, j, 1])
      for (t in 2:n.years) {	    
        # Process Model ---------------------------------------------------
        logit(psi.bbs[i, j, t]) <- beta.0[i, t] + 
                                   beta.1[i] * ELEV.bbs[j] + 
          		           beta.2[i] * pow(ELEV.bbs[j], 2) + 
		                   beta.3[i] * FOREST.bbs[j] + 
			           phi[i] * z.bbs[i, j, t - 1]
        z.bbs[i, j, t] ~ dbern(psi.bbs[i, j, t])
      } # t
    } # i
  } # j

  # BBS Data --------------------------------------------------------------
  for (i in 1:n.vals.bbs) {
    logit(pi.2[i]) <- gamma.2.0[sp.indx.bbs.p[i], year.indx.bbs[i]] +
	               gamma.2.1[sp.indx.bbs.p[i]] * DAY.bbs[i] +
		       gamma.2.2[sp.indx.bbs.p[i]] * pow(DAY.bbs[i], 2) +
		       gamma.2.3[obsv.bbs[i]]
    v.2[i] ~ dbern(pi.2[i] * z.bbs[sp.indx.bbs[i], site.bbs[i], year.indx.bbs[i]])
  } # i

  # NEON Data -----------------------------------------------------------
  for (i in 1:n.vals.neon) {
    logit(pi.1[i]) <- gamma.1.0[sp.indx.neon.p[i], year.indx.neon[i]] +
	                gamma.1.1[sp.indx.neon.p[i]] * DAY.neon[i] +
			gamma.1.2[sp.indx.neon.p[i]] * pow(DAY.neon[i], 2) +
			gamma.1.3[sp.indx.neon.p[i]] * HOUR.neon[i]
    v.1[i] ~ dbern(pi.1[i] * z.neon[sp.indx.neon.p[i], site.neon[i], year.indx.neon[i]])
    DAY.neon[i] ~ dnorm(0, 1)
    HOUR.neon[i] ~ dnorm(0, 1)
  } # i
})


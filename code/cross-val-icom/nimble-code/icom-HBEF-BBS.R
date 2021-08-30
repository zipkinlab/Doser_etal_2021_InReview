# icom-HBEF-BBS.R: BUGS code to run ICOM using HBEF and BBS data for the 
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
  int.alpha.mean ~ dunif(0, 1) # overall (species and year) HBEF detection
  alpha.0.mean <- logit(int.alpha.mean)
  alpha.1.mean ~ dnorm(0, 0.368) # day
  alpha.2.mean ~ dnorm(0, 0.368) # day^2
  alpha.3.mean ~ dnorm(0, 0.368) # time of day
  int.gamma.2.mean ~ dunif(0, 1) # overall (species and year) BBS detection
  gamma.2.0.mean <- logit(int.gamma.2.mean)
  gamma.2.1.mean ~ dnorm(0, 0.368) # day
  gamma.2.2.mean ~ dnorm(0, 0.368) # day^2
  # Precision parameters --------------
  tau.beta.1 ~ dgamma(0.1, 0.1)
  tau.beta.2 ~ dgamma(0.1, 0.1)
  tau.beta.3 ~ dgamma(0.1, 0.1)
  tau.phi ~ dgamma(0.1, 0.1)
  tau.alpha.0 ~ dgamma(0.1, 0.1)
  tau.alpha.1 ~ dgamma(0.1, 0.1)
  tau.alpha.2 ~ dgamma(0.1, 0.1)
  tau.alpha.3 ~ dgamma(0.1, 0.1)
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
  # HBEF ------------------------------
  for (i in 1:I.hbef) {
    for (t in 1:n.years) {
      alpha.0[i, t] ~ dnorm(alpha.0.mean, tau.alpha.0)
      logit(int.alpha[i, t]) <- alpha.0[i, t]
    } # t
    alpha.1[i] ~ dnorm(alpha.1.mean, tau.alpha.1)
    alpha.2[i] ~ dnorm(alpha.2.mean, tau.alpha.2)
    alpha.3[i] ~ dnorm(alpha.3.mean, tau.alpha.3)
  } # i
  # BBS -------------------------------
  for (i in 1:I.bbs) {
    for (t in 1:n.years) {
      gamma.2.0[i, t] ~ dnorm(gamma.2.0.mean, tau.gamma.2.0)
      logit(int.gamma.2[i, t]) <- gamma.2.0[i, t]
    } # t
    gamma.2.1[i] ~ dnorm(gamma.2.1.mean, tau.gamma.2.1)
    gamma.2.2[i] ~ dnorm(gamma.2.2.mean, tau.gamma.2.2)
  } # i

  # Process Model ---------------------------------------------------------
  # HBEF ------------------------------
  for (j in 1:J.hbef) {
    for (i in 1:I) {
      logit(psi.hbef[i, j, 1]) <- beta.0[i, 1] + 
	                         beta.1[i] * ELEV.hbef[j] + 
			         beta.2[i] * pow(ELEV.hbef[j], 2) + 
			         beta.3[i] * FOREST.hbef[j] 
      z.hbef[i, j, 1] ~ dbern(psi.hbef[i, j, 1])
      for (t in 2:n.years) {	    
        # Process Model ---------------------------------------------------
        logit(psi.hbef[i, j, t]) <- beta.0[i, t] + 
                                   beta.1[i] * ELEV.hbef[j] + 
          		           beta.2[i] * pow(ELEV.hbef[j], 2) + 
		                   beta.3[i] * FOREST.hbef[j] + 
			           phi[i] * z.hbef[i, j, t - 1]
        z.hbef[i, j, t] ~ dbern(psi.hbef[i, j, t])
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


  # Hubbard Brook Data ----------------------------------------------------
  for (i in 1:n.vals.hbef) {
    logit(p[i]) <- alpha.0[sp.indx.hbef.p[i], year.indx.hbef[i]] +
                      alpha.1[sp.indx.hbef.p[i]] * DAY.hbef[i] +
        	      alpha.2[sp.indx.hbef.p[i]] * pow(DAY.hbef[i], 2) +
        	      alpha.3[sp.indx.hbef.p[i]] * TOD.hbef[i]
    y[i] ~ dbern(p[i] * z.hbef[sp.indx.hbef[i], site.hbef[i], year.indx.hbef[i]])
    DAY.hbef[i] ~ dnorm(0, 1)
    TOD.hbef[i] ~ dnorm(0, 1)
  } # i
  # BBS Data --------------------------------------------------------------
  for (i in 1:n.vals.bbs) {
    logit(pi.2[i]) <- gamma.2.0[sp.indx.bbs.p[i], year.indx.bbs[i]] +
	               gamma.2.1[sp.indx.bbs.p[i]] * DAY.bbs[i] +
		       gamma.2.2[sp.indx.bbs.p[i]] * pow(DAY.bbs[i], 2) +
		       gamma.2.3[obsv.bbs[i]]
    v.2[i] ~ dbern(pi.2[i] * z.bbs[sp.indx.bbs[i], site.bbs[i], year.indx.bbs[i]])
  } # i
})


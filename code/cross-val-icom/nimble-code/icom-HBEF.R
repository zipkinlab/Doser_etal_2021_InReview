# icom-HBEF.R: BUGS code to run community model in NIMBLE using only data
#              from Hubbard Brook Experimental Forest for the foliage-gleaning
#              bird case study, which speeds things up. 
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
  # Precision Parameters --------------
  tau.beta.1 ~ dgamma(0.1, 0.1)
  tau.beta.2 ~ dgamma(0.1, 0.1)
  tau.beta.3 ~ dgamma(0.1, 0.1)
  tau.phi ~ dgamma(0.1, 0.1)
  tau.alpha.0 ~ dgamma(0.1, 0.1)
  tau.alpha.1 ~ dgamma(0.1, 0.1)
  tau.alpha.2 ~ dgamma(0.1, 0.1)
  tau.alpha.3 ~ dgamma(0.1, 0.1)

  # Species-specific coefficients -----------------------------------------
  for (i in 1:I) {
    beta.1[i] ~ dnorm(beta.1.mean, tau.beta.1)
    beta.2[i] ~ dnorm(beta.2.mean, tau.beta.2)
    beta.3[i] ~ dnorm(beta.3.mean, tau.beta.3)
    phi[i] ~ dnorm(phi.mean, tau.phi)
    for (t in 1:n.years) {
      beta.0[i, t] ~ dnorm(beta.0.mean[t], tau.beta.0[t])
      logit(int.beta[i, t]) <- beta.0[i, t]
      alpha.0[i, t] ~ dnorm(alpha.0.mean, tau.alpha.0)
      logit(int.alpha[i, t]) <- alpha.0[i, t]
    } # t
    alpha.1[i] ~ dnorm(alpha.1.mean, tau.alpha.1)
    alpha.2[i] ~ dnorm(alpha.2.mean, tau.alpha.2)
    alpha.3[i] ~ dnorm(alpha.3.mean, tau.alpha.3)
  } # i

  # Process Models and Likelihoods ----------------------------------------
  for (j in 1:J) {
    for (i in 1:I) {
      logit(psi[i, j, 1]) <- beta.0[i, 1] + 
	                     beta.1[i] * ELEV[j] + 
			     beta.2[i] * pow(ELEV[j], 2) + 
			     beta.3[i] * FOREST[j] 
      z.hbef[i, j, 1] ~ dbern(psi[i, j, 1])
      for (t in 2:n.years) {	    
        # Process Model ---------------------------------------------------
        logit(psi[i, j, t]) <- beta.0[i, t] + 
                               beta.1[i] * ELEV[j] + 
          		       beta.2[i] * pow(ELEV[j], 2) + 
		               beta.3[i] * FOREST[j] + 
			       phi[i] * z.hbef[i, j, t - 1]
        z.hbef[i, j, t] ~ dbern(psi[i, j, t])
      } # t
    } # i
  } # j

  # Hubbard Brook Data ----------------------------------------------------
  # Data are stacked in a single vector as opposed to a multi-dimensional 
  # array to improve computational performance. 
  for (i in 1:n.vals.hbef) {
    logit(p[i]) <- alpha.0[sp.indx.hbef[i], year.indx.hbef[i]] +
                   alpha.1[sp.indx.hbef[i]] * DAY.hbef[i] +
        	   alpha.2[sp.indx.hbef[i]] * pow(DAY.hbef[i], 2) +
        	   alpha.3[sp.indx.hbef[i]] * TOD.hbef[i]
    y[i] ~ dbern(p[i] * z.hbef[sp.indx.hbef[i], site.hbef[i], year.indx.hbef[i]])
    DAY.hbef[i] ~ dnorm(0, 1)
    TOD.hbef[i] ~ dnorm(0, 1)
  } # i
})

# idm-HBEF-BBS.R: BUGS code to run IDM in NIMBLE using HBEF and BBS
#                 data for the foliage-gleaning bird case study. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 
require(nimble)
idm.code <- nimbleCode({
  # Priors ------------------------------------------------------------------
  beta.1 ~ dnorm(0, 0.368) # linear elevation
  beta.2 ~ dnorm(0, 0.368) # quadratic elevation
  beta.3 ~ dnorm(0, 0.368) # percent forest
  phi ~ dnorm(0, 0.368) # auto-logistic parameter
  int.alpha.mean ~ dunif(0, 1) # overall (species and year) HBEF detection
  alpha.0.mean <- logit(int.alpha.mean)
  alpha.1 ~ dnorm(0, 0.368) # day
  alpha.2 ~ dnorm(0, 0.368) # day^2
  alpha.3 ~ dnorm(0, 0.368) # time of day
  int.gamma.2.mean ~ dunif(0, 1) # overall (species and year) BBS detection
  gamma.2.0.mean <- logit(int.gamma.2.mean)
  gamma.2.1 ~ dnorm(0, 0.368) # day
  gamma.2.2 ~ dnorm(0, 0.368) # day^2
  tau.gamma.2.3 ~ dgamma(0.1, 0.1) # Random observer effect
  for (i in 1:n.obsv.bbs) {
    gamma.2.3[i] ~ dnorm(0, tau.gamma.2.3)
  } # i
  tau.alpha.0 ~ dgamma(0.1, 0.1)
  tau.gamma.2.0 ~ dgamma(0.1, 0.1)
  # Intercept for all years
  for (t in 1:n.years) {
    int.beta[t] ~ dunif(0, 1)
    beta.0[t] <- logit(int.beta[t])
  } # t
  for (t in 1:n.years) {
    alpha.0[t] ~ dnorm(alpha.0.mean, tau.alpha.0)
    gamma.2.0[t] ~ dnorm(gamma.2.0.mean, tau.gamma.2.0)
  }

  # Process Models and Likelihoods ----------------------------------------
  # HBEF ------------------------------
  for (j in 1:J.hbef) {
    logit(psi.hbef[j, 1]) <- beta.0[1] + 
	                     beta.1 * ELEV.hbef[j] + 
		             beta.2 * pow(ELEV.hbef[j], 2) + 
			     beta.3 * FOREST.hbef[j]
    z.hbef[j, 1] ~ dbern(psi.hbef[j, 1])
    for (t in 2:n.years) {
      logit(psi.hbef[j, t]) <- beta.0[t] + 
                            beta.1 * ELEV.hbef[j] + 
      		            beta.2 * pow(ELEV.hbef[j], 2) + 
			    beta.3 * FOREST.hbef[j] + 
			    phi * z.hbef[j, t - 1]
      z.hbef[j, t] ~ dbern(psi.hbef[j, t])
    } # t
  } # j
  # BBS -------------------------------
  for (j in 1:J.bbs) {
    logit(psi.bbs[j, 1]) <- beta.0[1] + 
                            beta.1 * ELEV.bbs[j] + 
        	            beta.2 * pow(ELEV.bbs[j], 2) + 
        		    beta.3 * FOREST.bbs[j]
    z.bbs[j, 1] ~ dbern(psi.bbs[j, 1])
    for (t in 2:n.years) {
    logit(psi.bbs[j, t]) <- beta.0[t] + 
                            beta.1 * ELEV.bbs[j] + 
      		            beta.2 * pow(ELEV.bbs[j], 2) + 
        		    beta.3 * FOREST.bbs[j] + 
        		    phi * z.bbs[j, t - 1]
      z.bbs[j, t] ~ dbern(psi.bbs[j, t])
    } # t
  } # j
  # Hubbard Brook Data ----------------------------------------------------
  for (i in 1:n.vals.hbef) {
    logit(p[i]) <- alpha.0[year.indx.hbef[i]] + 
                   alpha.1 * DAY.hbef[i] + 
        	   alpha.2 * pow(DAY.hbef[i], 2) + 
        	   alpha.3 * TOD.hbef[i]
    y[i] ~ dbern(p[i] * z.hbef[site.hbef[i], year.indx.hbef[i]])
    DAY.hbef[i] ~ dnorm(0, 1)
    TOD.hbef[i] ~ dnorm(0, 1)
  } # i
  # BBS Data --------------------------------------------------------------
  for (i in 1:n.vals.bbs) {
    logit(pi.2[i]) <- gamma.2.0[year.indx.bbs[i]] + 
                      gamma.2.1 * DAY.bbs[i] + 
        	      gamma.2.2 * pow(DAY.bbs[i], 2) + 
        	      gamma.2.3[obsv.bbs[i]]
    v.2[i] ~ dbern(pi.2[i] * z.bbs[site.bbs[i], year.indx.bbs[i]])
  } # i
})


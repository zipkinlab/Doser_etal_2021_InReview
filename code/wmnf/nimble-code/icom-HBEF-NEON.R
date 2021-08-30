# icom-HBEF-NEON.R: BUGS code to run ICOM in NIMBLE using data from HBEF and 
#                   NEON for the foliage-gleaning bird case study. 
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
  tau.alpha.0 ~ dgamma(0.1, 0.1)
  tau.alpha.1 ~ dgamma(0.1, 0.1)
  tau.alpha.2 ~ dgamma(0.1, 0.1)
  tau.alpha.3 ~ dgamma(0.1, 0.1)
  tau.gamma.1.0 ~ dgamma(0.1, 0.1)
  tau.gamma.1.1 ~ dgamma(0.1, 0.1)
  tau.gamma.1.2 ~ dgamma(0.1, 0.1)
  tau.gamma.1.3 ~ dgamma(0.1, 0.1)

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
  # NEON ------------------------------
  for (i in 1:I.neon) {
    for (t in 1:n.years.neon) {
      gamma.1.0[i, t] ~ dnorm(gamma.1.0.mean, tau.gamma.1.0)
    } # t
    gamma.1.1[i] ~ dnorm(gamma.1.1.mean, tau.gamma.1.1)
    gamma.1.2[i] ~ dnorm(gamma.1.2.mean, tau.gamma.1.2)
    gamma.1.3[i] ~ dnorm(gamma.1.3.mean, tau.gamma.1.3)
  } # i

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

  # Hubbard Brook Data ----------------------------------------------------
  for (i in 1:n.vals.hbef) {
    logit(p[i]) <- alpha.0[sp.indx.hbef.p[i], year.indx.hbef[i]] +
                      alpha.1[sp.indx.hbef.p[i]] * DAY.hbef[i] +
        	      alpha.2[sp.indx.hbef.p[i]] * pow(DAY.hbef[i], 2) +
        	      alpha.3[sp.indx.hbef.p[i]] * TOD.hbef[i]
    y[i] ~ dbern(p[i] * z.hbef[sp.indx.hbef[i], site.hbef[i], year.indx.hbef[i]])
    # Replicate of data for Bayesian p-value
    y.rep[i] ~ dbern(p[i] * z.hbef[sp.indx.hbef[i], site.hbef[i], year.indx.hbef[i]])
    E.hbef[i] <- p[i] * z.hbef[sp.indx.hbef[i], site.hbef[i], year.indx.hbef[i]]
    DAY.hbef[i] ~ dnorm(0, 1)
    TOD.hbef[i] ~ dnorm(0, 1)
  } # i

  # NEON Data -----------------------------------------------------------
  for (i in 1:n.vals.neon) {
    logit(pi.1[i]) <- gamma.1.0[sp.indx.neon.p[i], year.indx.neon[i]] +
	                gamma.1.1[sp.indx.neon.p[i]] * DAY.neon[i] +
			gamma.1.2[sp.indx.neon.p[i]] * pow(DAY.neon[i], 2) +
			gamma.1.3[sp.indx.neon.p[i]] * HOUR.neon[i]
    v.1[i] ~ dbern(pi.1[i] * z.neon[sp.indx.neon.p[i], site.neon[i], year.indx.neon[i]])
    v.1.rep[i] ~ dbern(pi.1[i] * z.neon[sp.indx.neon.p[i], site.neon[i], year.indx.neon[i]])
    E.neon[i] <- pi.1[i] * z.neon[sp.indx.neon.p[i], site.neon[i], year.indx.neon[i]]
    DAY.neon[i] ~ dnorm(0, 1)
    HOUR.neon[i] ~ dnorm(0, 1)
  } # i

  # Bayesian P-Value ------------------------------------------------------
  # HBEF ------------------------------
  for (i in 1:I) {
    for (t in 1:n.years) {
      y.grouped[i, t] <- sum(y[low.bp.y[i, t]:high.bp.y[i, t]])
      y.rep.grouped[i, t] <- sum(y.rep[low.bp.y[i, t]:high.bp.y[i, t]])
      E.grouped.y[i, t] <- sum(E.hbef[low.bp.y[i, t]:high.bp.y[i, t]])
      x2.y[i, t] <- pow((y.grouped[i, t] - E.grouped.y[i, t]), 2) / (E.grouped.y[i, t] + e)
      x2.y.rep[i, t] <- pow((y.rep.grouped[i, t] - E.grouped.y[i, t]), 2) / (E.grouped.y[i, t] + e)
    } # t
  } # i
  # Add up Chi-square discrepancies
  chi.2.y <- sum(x2.y[1:I, 1:n.years])
  chi.2.rep.y <- sum(x2.y.rep[1:I, 1:n.years])
  # NEON ------------------------------
  for (i in 1:I.neon) {
    for (t in 1:n.years.neon) {
      v.1.grouped[i, t] <- sum(v.1[low.bp.v.1[i, t]:high.bp.v.1[i, t]])
      v.1.rep.grouped[i, t] <- sum(v.1.rep[low.bp.v.1[i, t]:high.bp.v.1[i, t]])
      E.grouped.v.1[i, t] <- sum(E.neon[low.bp.v.1[i, t]:high.bp.v.1[i, t]])
      x2.v.1[i, t] <- pow((v.1.grouped[i, t] - E.grouped.v.1[i, t]), 2) / (E.grouped.v.1[i, t] + e)
      x2.v.1.rep[i, t] <- pow((v.1.rep.grouped[i, t] - E.grouped.v.1[i, t]), 2) / (E.grouped.v.1[i, t] + e)
    } # t
  } # i
  # Add up Chi-square discrepancies
  chi.2.v.1 <- sum(x2.v.1[1:I.neon, 1:n.years.neon])
  chi.2.rep.v.1 <- sum(x2.v.1.rep[1:I.neon, 1:n.years.neon])

})


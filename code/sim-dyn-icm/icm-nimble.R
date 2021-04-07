# icm-nimble.R: BUGS code to run ICM through NIMBLE in R.
require(nimble)
icm.code <- nimbleCode ({
  # Priors ----------------------------------------------------------------
  # Initial occupancy
  int.psi.0.mean ~ dunif(0, 1) # intercept
  beta.psi.0.mean <- logit(int.psi.0.mean) # intercept
  beta.psi.1.mean ~ dnorm(0, 0.1) # covariate effect
  # Persistence
  int.phi.0.mean ~ dunif(0, 1) # intercept
  beta.phi.0.mean <- logit(int.phi.0.mean) # intercept
  beta.phi.1.mean ~ dnorm(0, 0.1) # covariate effect
  # Colonization
  int.gamma.0.mean ~ dunif(0, 1) # intercept
  beta.gamma.0.mean <- logit(int.gamma.0.mean) # intercept
  beta.gamma.1.mean ~ dnorm(0, 0.1) # covariate effect
  # eBird detection
  int.alpha.eb.0.mean ~ dunif(0, 1) # intercept 
  alpha.eb.0.mean <- logit(int.alpha.eb.0.mean) # intercept
  alpha.eb.1.mean ~ dnorm(0, 0.1) # covariate effect
  # DET data detection
  int.alpha.det.0.mean ~ dunif(0, 1) # intercept
  alpha.det.0.mean <- logit(int.alpha.det.0.mean) # intercept
  alpha.det.1.mean ~ dnorm(0, 0.1) # covariate effect
  # NEON data detection
  int.alpha.neon.0.mean ~ dunif(0, 1) # intercept
  alpha.neon.0.mean <- logit(int.alpha.neon.0.mean) # intercept
  alpha.neon.1.mean ~ dnorm(0, 0.1) # covariate effect
  # BBS data detection
  int.alpha.bbs.0.mean ~ dunif(0, 1) # intercept 
  alpha.bbs.0.mean <- logit(int.alpha.bbs.0.mean) # intercept
  alpha.bbs.1.mean ~ dnorm(0, 0.1) # covariate effect
  # Precision parameters for community level
  tau.beta.psi.0 ~ dgamma(0.1, 0.1) 
  tau.beta.psi.1 ~ dgamma(0.1, 0.1)
  tau.beta.phi.0 ~ dgamma(0.1, 0.1)
  tau.beta.phi.1 ~ dgamma(0.1, 0.1)
  tau.beta.gamma.0 ~ dgamma(0.1, 0.1)
  tau.beta.gamma.1 ~ dgamma(0.1, 0.1)
  tau.alpha.eb.0 ~ dgamma(0.1, 0.1)
  tau.alpha.eb.1 ~ dgamma(0.1, 0.1)
  tau.alpha.det.0 ~ dgamma(0.1, 0.1)
  tau.alpha.det.1 ~ dgamma(0.1, 0.1)
  tau.alpha.neon.0 ~ dgamma(0.1, 0.1)
  tau.alpha.neon.1 ~ dgamma(0.1, 0.1)
  tau.alpha.bbs.0 ~ dgamma(0.1, 0.1)
  tau.alpha.bbs.1 ~ dgamma(0.1, 0.1)

  # Species-specific coefficients -----------------------------------------
  for (i in 1:I) {
    beta.psi.0[i] ~ dnorm(beta.psi.0.mean, tau.beta.psi.0)
    beta.psi.1[i] ~ dnorm(beta.psi.1.mean, tau.beta.psi.1)
    beta.phi.0[i] ~ dnorm(beta.phi.0.mean, tau.beta.phi.0)
    beta.phi.1[i] ~ dnorm(beta.phi.1.mean, tau.beta.phi.1)
    beta.gamma.0[i] ~ dnorm(beta.gamma.0.mean, tau.beta.gamma.0)
    beta.gamma.1[i] ~ dnorm(beta.gamma.1.mean, tau.beta.gamma.1)
    alpha.eb.0[i] ~ dnorm(alpha.eb.0.mean, tau.alpha.eb.0)
    alpha.eb.1[i] ~ dnorm(alpha.eb.1.mean, tau.alpha.eb.1)
    alpha.det.0[i] ~ dnorm(alpha.det.0.mean, tau.alpha.det.0)
    alpha.det.1[i] ~ dnorm(alpha.det.1.mean, tau.alpha.det.1)
    alpha.neon.0[i] ~ dnorm(alpha.neon.0.mean, tau.alpha.neon.0)
    alpha.neon.1[i] ~ dnorm(alpha.neon.1.mean, tau.alpha.neon.1)
    alpha.bbs.0[i] ~ dnorm(alpha.bbs.0.mean, tau.alpha.bbs.0)
    alpha.bbs.1[i] ~ dnorm(alpha.bbs.1.mean, tau.alpha.bbs.1)
  }

  # Process models --------------------------------------------------------
  for (i in 1:I) {
    # Initial occupancy process -------
    for (j in 1:J) {
        logit(psi.0[i, j]) <- beta.psi.0[i] + beta.psi.1[i] * X.psi[j, 2]
        z[i, j, 1] ~ dbern(psi.0[i, j])
    } # j 

    for (t in 1:(n.years - 1)) {
      # Persistence and colonization --
      for (j in 1:J) {
        logit(phi[i, j, t]) <- beta.phi.0[i] + beta.phi.1[i] * X.phi[j, t, 2]
        logit(gamma[i, j, t]) <- beta.gamma.0[i] + 
                                 beta.gamma.1[i] * X.gamma[j, t, 2]
        z[i, j, t + 1] ~ dbern(z[i, j, t] * phi[i, j, t] + (1 - z[i, j, t]) * gamma[i, j, t])
      } # j
    } # t
    
    for (t in 1:n.years) {
      # Likelihoods -------------------------------------------------------
      # DET Data -----------------------
      for (j in 1:J.det) {
        for (k in 1:K.det) {
          logit(p.det[i, j, k, t]) <- alpha.det.0[i] + alpha.det.1[i] * X.det[j, k, t, 2]
          C[i, j, k, t] ~ dbern(p.det[i, j, k, t] * z[i, pixel.det[j], t])
        } # k 
      } # j 

      # eBird Data --------------------
      for (g in 1:J.eb) {
        # Change of support for eBird - 
        z.eb[i, g, t] <- step(sum(z[i, low[g]:high[g], t]) - 1)
        for (k in 1:K.eb) {
          logit(p.eb[i, g, k, t]) <- alpha.eb.0[i] + alpha.eb.1[i] * X.eb[g, k, t, 2]
          y[i, g, k, t] ~ dbern(z.eb[i, g, t] * p.eb[i, g, k, t])
        } # k
      } # g
     
      # NEON Data ---------------------
      for (j in 1:J.neon) {
        for (k in 1:K.neon) {
          logit(p.neon[i, j, k, t]) <- alpha.neon.0[i] + 
                                       alpha.neon.1[i] * X.neon[j, k, t, 2]
          x[i, j, k, t] ~ dbern(p.neon[i, j, k, t] * z[i, pixel.neon[j], t])
        } # k
      } # j

      # BBS Data ---------------------
      for (j in 1:J.bbs) {
        logit(p.bbs[i, j, t]) <- alpha.bbs.0[i] + 
                                    alpha.bbs.1[i] * X.bbs[j, t, 2]
        y.bbs[i, j, t] ~ dbern(p.bbs[i, j, t] * z[i, pixel.bbs[j], t])
      } # j
    } # t
  } # i
})


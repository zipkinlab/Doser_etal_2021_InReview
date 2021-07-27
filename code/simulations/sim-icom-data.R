# sim-icom-data.R: script to simulate data following the ICOM framework. Data are
#                 simulated for three data sets, specifically for one replicated 
#                 data source and two nonreplicated data sources, although the 
#                 number of replicates for the replicated data source and 
#                 the first replicated data source are allowed to vary. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

sim.icom.data <- function(J.rep, J.nrep.1, J.nrep.2, 
			 beta.0, beta.1, phi,
			 gamma.2.0, gamma.2.1,
			 alpha.0, alpha.1, 
			 gamma.1.0, gamma.1.1,
			 K.rep, K.nrep.1, I, n.years) { 
  # J.rep = number of sites with replicated data set.
  # J.nrep.1 = number of sites with first nonreplicated data set.
  # J.nrep.2 = number of sites with second nonreplicated data set. 
  # beta.0 = (I x n.years) matrix containing species and year specific
  #          occurrence intercepts
  # beta.1 = Length I vector containing species-specific spatial 
  #          covariate effects on occurrence
  # phi = length I vector containing autologistic species-specific 
  #       effects on occurrence
  # gamma.2.0, alpha.0, gamma.1.0 = (I x n.years) matrix
  #       containing species and year specific detection intercepts for 
  #       nonreplicated data set 2, replicated data set, and 
  #       nonreplicated data set 1, respectively
  # gamma.2.1, alpha.1, gamma.1.1 = Length I vector 
  #       containing species specific space/time covariate effects 
  #       on detection for nonreplicated data set 2, replicated data set, 
  #       and nonreplicated data set 1, respectively. 
  # K.rep = number of replicates for replicated data set.
  # K.nrep.1 = number of replicates for first nonreplicated data set. If this
  #            is in fact nonreplicated, should be set to 1. If set higher, 
  #            this will be a replicated data set. 
  # I = number of species in community. 
  # n.years = number of years data are collected. 

  # Subroutines -----------------------------------------------------------
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Form spatial occupancy covariate --------------------------------------
  J <- J.rep + J.nrep.1 + J.nrep.2
  # Occupancy covariate varies across space
  X.psi <- array(1, dim = c(J, 2))
  X.psi[, 2] <- rnorm(J)

  # Form detection variable for replicated data set -----------------------
  # Detection is a function of space and time
  n.alpha <- 2
  X.rep <- array(1, dim = c(J.rep, K.rep, n.years, n.alpha))
  if (n.alpha > 1) {
    for (a in 2:n.alpha) {
      for (j in 1:J.rep) {
        for (t in 1:n.years) {
          X.rep[j, , t, a] <- rnorm(K.rep)
	} # t
      } # j
    } # a
  }

  # Form detection variable for nonreplicated data set 1 ------------------
  # Detection is a function of space and time
  n.gamma.1 <- 2
  X.nrep.1 <- array(1, dim = c(J.nrep.1, K.nrep.1, n.years, n.gamma.1))
  if (n.gamma.1 > 1) {
    for (a in 2:n.gamma.1) {
      for (j in 1:J.nrep.1) {
        for (t in 1:n.years) {
          X.nrep.1[j, , t, a] <- rnorm(K.nrep.1)
	} # t
      } # j
    } # a
  }

  # Form detection variable for nonreplicated data set 2 ------------------
  # Detection is a function of space and time
  n.gamma.2 <- 2
  X.nrep.2 <- array(1, dim = c(J.nrep.2, n.years, n.gamma.2))
  if (n.gamma.2 > 1) {
    for (a in 2:n.gamma.2) {
      for (j in 1:J.nrep.2) {
	  X.nrep.2[j, , a] <- rnorm(n.years)
      } # j
    } # a
  }

  # Define occurrence at each site ----------------------------------------
  # Note: sites are stored in the following order: all replicated sites, 
  #       all nonreplicated 1 sites, all nonreplicated 2 sites. 
  # Occurrence probability array
  psi <- array(NA, dim = c(I, J, n.years))
  # Latent occurrence state (presence/absence)
  z <- array(NA, dim = c(I, J, n.years))
  for (t in 1:n.years) {
    for (i in 1:I) {
      for (j in 1:J) {
        if (t == 1) {
          psi[i, j, t] <- logit.inv(beta.0[i, 1] * X.psi[j, 1] + 
          		            beta.1[i] * X.psi[j, 2])
        } else {
            psi[i, j, t] <- logit.inv(beta.0[i, t] * X.psi[j, 1] + 
          			      beta.1[i] * X.psi[j, 2] + 
          			      phi[i] * z[i, j, t - 1])  
          }
        # Occurrence
        z[i, j, t] <- rbinom(1, 1, psi[i, j, t]) 
      } # j
    } # i 
  } # t

  # Replicated data set formation -----------------------------------------
  # Detection probability for replicated data set. 
  p <- array(NA, dim = c(I, J.rep, K.rep, n.years))
  # y is the replicated detection-nondetection data.
  y <- array(NA, dim = c(I, J.rep, K.rep, n.years))
  for (i in 1:I) {
    for (j in 1:J.rep) {
      for (k in 1:K.rep) {
	for (t in 1:n.years) {
	  p[i, j, k, t] <- logit.inv(alpha.0[i, t] * X.rep[j, k, t, 1] + 
				     alpha.1[i] * X.rep[j, k, t, 2])
          y[i, j, k, t] <- rbinom(1, 1, p[i, j, k, t] * z[i, j, t])
	} # t
      } # k
    } # j
  } # i

  # Nonreplicated data set 1 formation ------------------------------------
  # Detection probability for nonreplicated data set 1. 
  pi.1 <- array(NA, c(I, J.nrep.1, K.nrep.1, n.years))
  # v.1 is the nonreplicated detection-nondetection data.  
  v.1 <- array(NA, dim = c(I, J.nrep.1, K.nrep.1, n.years))
  for (j in 1:J.nrep.1) {
    for (i in 1:I) {
      for (t in 1:n.years) {
        for (k in 1:K.nrep.1) {
          pi.1[i, j, k, t] <- logit.inv(gamma.1.0[i, t] * X.nrep.1[j, k, t, 1] + 
				        gamma.1.1[i] * X.nrep.1[j, k, t, 2])
          # Note the "j + J.rep" which ensures linking to the correct site
          v.1[i, j, k, t] <- rbinom(1, 1, z[i, j + J.rep, t] * pi.1[i, j, k, t])
        } # k
      } # t
    } # i
  } # j
  

  # Nonreplicated data set 2 formation ------------------------------------
  # Detection probability for nonreplicated data set 2. 
  pi.2 <- array(NA, c(I, J.nrep.2, n.years))
  # v.2 is the nonreplicated detection-nondetection data.  
  v.2 <- array(NA, dim = c(I, J.nrep.2, n.years))
  for (j in 1:J.nrep.2) {
    for (i in 1:I) {
      for (t in 1:n.years) {
          pi.2[i, j, t] <- logit.inv(gamma.2.0[i, t] * X.nrep.2[j, t, 1] + 
			             gamma.2.1[i] * X.nrep.2[j, t, 2])
          # Note the "j + J.rep + J.rep.1" which ensures linking to the correct site
          v.2[i, j, t] <- rbinom(1, 1, z[i, j + J.rep + J.nrep.1, t] * pi.2[i, j, t])
      } # t
    } # i
  } # j

  # Output all information needed for model. 
  return(list(X.psi = X.psi, X.nrep.2 = X.nrep.2, X.nrep.1 = X.nrep.1, 
	      X.rep = X.rep, psi = psi, z = z, 
	      p = p, y = y, pi.1 = pi.1, pi.2 = pi.2, 
	      v.1 = v.1, v.2 = v.2))
}

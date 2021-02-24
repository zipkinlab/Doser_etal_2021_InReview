# Function to simulate data for dynamic integrated community model. 
sim.icm.data <- function(R, R.eb, R.hb, R.neon, beta.psi.0, beta.psi.1, 
			 beta.phi.0, beta.phi.1, beta.gamma.0, beta.gamma.1,
			 alpha.eb.0, alpha.eb.1, alpha.hb.0, alpha.hb.1, 
			 alpha.neon.0, alpha.neon.1,
			 J.eb, J.hb, J.neon, K, n.years) { 
  # R = number of pixels (i.e., the leve at which ecological processes
  #                       are simulated)
  # R.eb = number of eBird cells (cells are bigger than pixels)
  # R.neon = number of pixels with NEON data
  # R.hb = number of pixels with HB data
  # beta.psi.0, beta.psi.1 = vectors of species-specific intercept 
  # and covariate effects, respectively, for initial occupancy
  # beta.phi.0, beta.phi.1 = vectors of species-specific intercept and 
  # covariate effects, respectively, for persistence
  # beta.gamma.0, beta.gamma.1 = vectors of species-specific intercept and 
  # covariate effects, respectively, for colonization
  # alpha.eb.0, alpha.eb.1, alpha.neon.0, alpha.neon.1, alpha.hb.0, 
  # alpha.hb.1 = intercepts and covariate effects for detection for the 
  #              three different data sources (species-specific)

  # Subroutines -----------------------------------------------------------
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Form initial occupancy covariate --------------------------------------
  # Initial occupancy covariate varies across space
  n.beta.psi <- 2
  X.psi <- array(1, dim = c(R, n.beta.psi))
  if (n.beta.psi > 1) {
    for (i in 2:n.beta.psi) {
      X.psi[, i] <- rnorm(R)
    } # i
  }
  # Form persistence covariate -----------------------------------------------
  # Persistence varies across space and time
  n.beta.phi <- 2
  X.phi <- array(1, dim = c(R, n.years - 1, n.beta.phi))
  if (n.beta.phi > 1) {
    for (i in 2:n.beta.phi) {
      for (t in 1:(n.years - 1)) {
        X.phi[, t, i] <- rnorm(R)
      } # t
    } # i
  }
  # Form colonization covariate -------------------------------------------
  # Colonization varies across space and time
  n.beta.gamma <- 2
  X.gamma <- array(1, dim = c(R, n.years - 1, n.beta.gamma))
  if (n.beta.gamma > 1) {
    for (i in 2:n.beta.gamma) {
      for (t in 1:(n.years - 1)) {
        X.gamma[, t, i] <- rnorm(R)
      } # t
    } # i
  }

  # Form effort/detection variable for eBird ------------------------------
  # Effort is a function of space and time
  n.alpha.eb <- 2
  X.eb <- array(1, dim = c(R.eb, J.eb, n.years, n.alpha.eb))
  if (n.alpha.eb > 1) {
    for (a in 2:n.alpha.eb) {
      for (i in 1:R.eb) {
        for (t in 1:n.years) {
          X.eb[i, , t, a] <- runif(J.eb, 0.01, 2) 
	} # t
      } # i
    } # a
  }
  
  # Form detection variable for HB ----------------------------------------
  # Detection is a function of space and time
  n.alpha.hb <- 2
  X.hb <- array(1, dim = c(R.hb, J.hb, n.years, n.alpha.hb))
  if (n.alpha.hb > 1) {
    for (a in 2:n.alpha.hb) {
      for (i in 1:R.hb) {
        for (t in 1:n.years) {
          X.hb[i, , t, a] <- runif(J.hb, 0.01, 2) 
	} # t
      } # i
    } # a
  }

  # Form detection variable for NEON --------------------------------------
  # Detection is a function of space and time
  n.alpha.neon <- 2
  X.neon <- array(1, dim = c(R.neon, J.neon, n.years, n.alpha.neon))
  if (n.alpha.neon > 1) {
    for (a in 2:n.alpha.neon) {
      for (i in 1:R.neon) {
        for (t in 1:n.years) {
          X.neon[i, , t, a] <- runif(J.neon, 0.01, 2) 
	} # t
      } # i
    } # a
  }

  # Define initial occupancy at pixel level -------------------------------
  psi.0 <- matrix(NA, nrow = R, ncol = K)
  z <- array(NA, dim = c(R, K, n.years))
  for (k in 1:K) {
    psi.0[, k] <- logit.inv(beta.psi.0[k] * X.psi[, 1] + beta.psi.1[k] * X.psi[, 2])
    # Initial occupancy
    z[, k, 1] <- rbinom(R, 1, psi.0[, k]) 
  } # k 

  # Mechanistic processes at pixel level ----------------------------------
  phi <- array(NA, dim = c(R, K, n.years - 1))
  gamma <- array(NA, dim = c(R, K, n.years - 1))
  for (k in 1:K) {
    for (t in 1:(n.years - 1)) {
      # Persistence
      phi[, k, t] <- logit.inv(beta.phi.0[k] * X.phi[, t, 1] + beta.phi.1[k] * X.phi[, t, 2])
      # Colonization
      gamma[, k, t] <- logit.inv(beta.gamma.0[k] * X.gamma[, t, 1] +
				 beta.gamma.1[k] * X.gamma[, t, 2])
      # Latent occupancy
      z[, k, t + 1] <- rbinom(R, 1, z[, k, t] * phi[, k, t] + (1 - z[, k, t]) * gamma[, k, t])
    } # t
  } # k

  # HB Data Formation -----------------------------------------------------
  # Pixel IDs where HB data are obtained
  # Note that counts could take place in the same pixel, which could be 
  # reasonable depending on the scale at which occupancy is modeled. 
  pixel.hb <- sample(1:R, R.hb, replace = TRUE)
  # Latent abundance process at each hb site
  z.hb <- z[pixel.hb, , ]

  # HB detection
  p.hb <- array(NA, dim = c(R.hb, J.hb, K, n.years))
  # C is the HB data. 
  C <- array(NA, dim = c(R.hb, J.hb, K, n.years))
  # Probably could do this more efficiently, but I'm lazy...
  for (k in 1:K) {
    for (i in 1:R.hb) {
      for (j in 1:J.hb) {
	for (t in 1:n.years) {
	  p.hb[i, j, k, t] <- logit.inv(alpha.hb.0[k] * X.hb[i, j, t, 1] + 
					alpha.hb.1[k] * X.hb[i, j, t, 2])
          C[i, j, k, t] <- rbinom(1, 1, p.hb[i, j, k, t] * z.hb[i, k, t])
	} # t
      } # j
    } # i
  } # k

  # eBird Data Formation --------------------------------------------------
  # Number of pixels within each eBird cell
  n.eb.size <- R / R.eb
  # Lowest pixel number for each eBird cell
  low <- seq(from = 1, to = R, by = n.eb.size)
  # Highest pixel number for each eBird cell
  high <- seq(from = n.eb.size, to = R, by = n.eb.size)
  # Define detection/sampling process varyig over time and space
  p.eb <- array(NA, dim = c(J.eb, R.eb, K, n.years))
  # eBird data stored in y
  y <- array(NA, dim = c(J.eb, R.eb, K, n.years))
  z.eb <- array(NA, dim = c(R.eb, K, n.years))
  for (l in 1:R.eb) {
    for (k in 1:K) {
      for (t in 1:n.years) {
	# Change of support
        z.eb[l, k, t] <- ifelse(sum(z[low[l]:high[l], k, t]) > 0, 1, 0)
        for (j in 1:J.eb) {
	  # eBird detection
          p.eb[j, l, k, t] <- logit.inv(alpha.eb.0[k] * X.eb[l, j, t, 1] + 
				        alpha.eb.1[k] * X.eb[l, j, t, 2])
          y[j, l, k, t] <- rbinom(1, 1, z.eb[l, k, t] * p.eb[j, l, k, t])
        } # j
      } # t
    } # k
  } # l

  # NEON Data Formation ---------------------------------------------------
  p.neon <- array(NA, c(R.neon, J.neon, K, n.years))
  # x is the NEON data. Bad name, will change eventually. 
  x <- array(NA, dim = c(R.neon, J.neon, K, n.years))
  # Here I simulate so that NEON data are not available in the same pixels
  # as the HB data. 
  nums <- which(!(1:R %in% pixel.hb))
  # These are the pixels at which we have NEON data. Could be multiple
  # point counts within a given pixel. 
  pixel.neon <- sample(nums, R.neon, replace = TRUE)
  z.neon <- z[pixel.neon, , ]
  for (i in 1:R.neon) {
    for (k in 1:K) {
      for (t in 1:n.years) {
        for (j in 1:J.neon) {
          p.neon[i, j, k, t] <- logit.inv(alpha.neon.0[k] * X.neon[i, j, t, 1] + 
					  alpha.neon.1[k] * X.neon[i, j, t, 2])
          x[i, j, k, t] <- rbinom(1, 1, z.neon[i, k, t] * p.neon[i, j, k, t])
        } # j
      } # t
    } # k
  } # i
  
  # Output all information needed for model. 
  return(list(X.psi = X.psi, X.phi = X.phi, X.gamma = X.gamma, 
	      X.eb = X.eb, X.neon = X.neon, X.hb = X.hb, psi.0 = psi.0, 
	      z = z, phi = phi, gamma = gamma, pixel.hb = pixel.hb, 
	      p.hb = p.hb, C = C, low = low, high = high, p.eb = p.eb, 
	      y = y, pixel.neon = pixel.neon, p.neon = p.neon, x = x
	      ))
}


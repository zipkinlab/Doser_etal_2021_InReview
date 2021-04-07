# Function to simulate data for dynamic integrated community model. 
sim.icm.data <- function(J, J.eb, J.det, J.neon, J.bbs, 
			 beta.psi.0, beta.psi.1, beta.phi.0, beta.phi.1, 
			 beta.gamma.0, beta.gamma.1, alpha.bbs.0, alpha.bbs.1,
			 alpha.eb.0, alpha.eb.1, alpha.det.0, alpha.det.1, 
			 alpha.neon.0, alpha.neon.1,
			 K.eb, K.det, K.neon, I, n.years, n.route, 
			 pixel.det = NULL, pixel.bbs = NULL, pixel.neon = NULL, 
			 seed = NULL) { 
  # J = number of pixels (i.e., the level at which ecological processes
  #                       are simulated)
  # J.eb = number of eBird cells (cells are bigger than pixels)
  # J.neon = number of pixels with NEON data
  # J.det = number of pixels with DET data
  # beta.psi.0, beta.psi.1 = vectors of species-specific intercept 
  # and covariate effects, respectively, for initial occupancy
  # beta.phi.0, beta.phi.1 = vectors of species-specific intercept and 
  # covariate effects, respectively, for persistence
  # beta.gamma.0, beta.gamma.1 = vectors of species-specific intercept and 
  # covariate effects, respectively, for colonization
  # alpha.eb.0, alpha.eb.1, alpha.neon.0, alpha.neon.1, alpha.det.0, 
  # alpha.det.1, alpha.bbs.0, alpha.bbs.1 = intercepts and covariate effects 
  #	for detection for the three different data sources (species-specific)
  # K.eb, K.det, K.neon = number of repeat visits for eBird, DET, and NEON data
  # I = number of species
  # n.years = number of years
  # n.route = number of BBS routes. 
  # pixel.det, pixel.bbs, pixel.neon, seed = variables used to control 
  # randomization process for simulations. 

  # Assume the study area is a square region with J pixels, where the first
  # pixel is defined in the lower left corner of the region, and the Jth pixel
  # is defined in the upper right corner of the region. 

  # Subroutines -----------------------------------------------------------
  if(length(seed) > 0) set.seed(seed)
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Form initial occupancy covariate --------------------------------------
  # Initial occupancy covariate varies across space
  n.beta.psi <- 2
  X.psi <- array(1, dim = c(J, n.beta.psi))
  if (n.beta.psi > 1) {
    for (i in 2:n.beta.psi) {
      X.psi[, i] <- rnorm(J)
    } # i
  }
  # Form persistence covariate -----------------------------------------------
  # Persistence varies across space and time
  n.beta.phi <- 2
  X.phi <- array(1, dim = c(J, n.years - 1, n.beta.phi))
  if (n.beta.phi > 1) {
    for (i in 2:n.beta.phi) {
      for (t in 1:(n.years - 1)) {
        X.phi[, t, i] <- rnorm(J)
      } # t
    } # i
  }
  # Form colonization covariate -------------------------------------------
  # Colonization varies across space and time
  n.beta.gamma <- 2
  X.gamma <- array(1, dim = c(J, n.years - 1, n.beta.gamma))
  if (n.beta.gamma > 1) {
    for (i in 2:n.beta.gamma) {
      for (t in 1:(n.years - 1)) {
        X.gamma[, t, i] <- rnorm(J)
      } # t
    } # i
  }

  # Form effort/detection variable for eBird ------------------------------
  # Effort is a function of space and time
  n.alpha.eb <- 2
  X.eb <- array(1, dim = c(J.eb, K.eb, n.years, n.alpha.eb))
  if (n.alpha.eb > 1) {
    for (a in 2:n.alpha.eb) {
      for (j in 1:J.eb) {
        for (t in 1:n.years) {
          X.eb[j, , t, a] <- runif(K.eb, 0.01, 2) 
	} # t
      } # j
    } # a
  }
  
  # Form detection variable for DET ----------------------------------------
  # Detection is a function of space and time
  n.alpha.det <- 2
  X.det <- array(1, dim = c(J.det, K.det, n.years, n.alpha.det))
  if (n.alpha.det > 1) {
    for (a in 2:n.alpha.det) {
      for (j in 1:J.det) {
        for (t in 1:n.years) {
          X.det[j, , t, a] <- runif(K.det, 0.01, 2) 
	} # t
      } # j
    } # a
  }

  # Form detection variable for NEON --------------------------------------
  # Detection is a function of space and time
  n.alpha.neon <- 2
  X.neon <- array(1, dim = c(J.neon, K.neon, n.years, n.alpha.neon))
  if (n.alpha.neon > 1) {
    for (a in 2:n.alpha.neon) {
      for (j in 1:J.neon) {
        for (t in 1:n.years) {
          X.neon[j, , t, a] <- runif(K.neon, 0.01, 2) 
	} # t
      } # j
    } # a
  }

  # Form detection variable for BBS --------------------------------------
  # Detection is a function of space and time
  n.alpha.bbs <- 2
  X.bbs <- array(1, dim = c(J.bbs, n.years, n.alpha.bbs))
  if (n.alpha.bbs > 1) {
    for (a in 2:n.alpha.bbs) {
      for (j in 1:J.bbs) {
	  X.bbs[j, , a] <- runif(n.years, 0.01, 2)
      } # j
    } # a
  }
  # Define initial occupancy at pixel level -------------------------------
  psi.0 <- matrix(NA, nrow = I, ncol = J)
  z <- array(NA, dim = c(I, J, n.years))
  for (i in 1:I) {
    psi.0[i, ] <- logit.inv(beta.psi.0[i] * X.psi[, 1] + beta.psi.1[i] * X.psi[, 2])
    # Initial occupancy
    z[i, , 1] <- rbinom(J, 1, psi.0[i, ]) 
  } # i 

  # Mechanistic processes at pixel level ----------------------------------
  phi <- array(NA, dim = c(I, J, n.years - 1))
  gamma <- array(NA, dim = c(I, J, n.years - 1))
  for (i in 1:I) {
    for (t in 1:(n.years - 1)) {
      # Persistence
      phi[i, , t] <- logit.inv(beta.phi.0[i] * X.phi[, t, 1] + beta.phi.1[i] * X.phi[, t, 2])
      # Colonization
      gamma[i, , t] <- logit.inv(beta.gamma.0[i] * X.gamma[, t, 1] +
				 beta.gamma.1[i] * X.gamma[, t, 2])
      # Latent occupancy
      z[i, , t + 1] <- rbinom(J, 1, z[i, , t] * phi[i, , t] + (1 - z[i, , t]) * gamma[i, , t])
    } # t
  } # i

  # DET Data Formation -----------------------------------------------------
  if (length(pixel.det) == 0) {
    start.det <- sample(1:J, 1)
    # Have DET data be longer in the vertical direction than in the 
    # horizontal direction
    n.cols.det <- ceiling(J.det / sqrt(J))
    vert.det <- ifelse(J.det > sqrt(J), sqrt(J), J.det)
    vert.indices <- c(start.det, rep(NA, vert.det - 1))
    for (i in 2:vert.det) {
      if (vert.indices[i - 1] - sqrt(J) > 0 & vert.indices[i - 1] <= start.det) {
        vert.indices[i] <- vert.indices[i - 1] - sqrt(J)
      } else if (vert.indices[i - 1] < vert.indices[1]) {
          vert.indices[i] <- vert.indices[1] + sqrt(J)
        } else {
            vert.indices[i] <- vert.indices[i - 1] + sqrt(J)
          }
    }
    horz.det <- J.det - vert.det
    pixel.det <- c(vert.indices, rep(NA, horz.det))
    if (J.det > sqrt(J)) {
      for (i in (vert.det + 1):J.det) {
        if (((pixel.det[i - vert.det] - 1) %% sqrt(J) != 0) & (pixel.det[i - vert.det] <= 
           ifelse((i %% sqrt(J)) == 0, vert.indices[sqrt(J)], vert.indices[i %% sqrt(J)])))	{
          pixel.det[i] <- pixel.det[i - vert.det] - 1 
        } else if (pixel.det[i - vert.det] > (vert.indices[ifelse((i %% sqrt(J)) == 0, sqrt(J), 
            						      i %% sqrt(J))])) {
            pixel.det[i] <- pixel.det[i - vert.det] + 1
          } else {
              pixel.det[i] <- vert.indices[ifelse((i - vert.det) %% sqrt(J) == 0, sqrt(J), 
            				       i %% sqrt(J))] + 1
            }
      }
    }
  }
  # DET detection
  p.det <- array(NA, dim = c(I, J.det, K.det, n.years))
  # C is the DET data. 
  C <- array(NA, dim = c(I, J.det, K.det, n.years))
  # Probably could do this more efficiently, but I'm lazy...
  for (i in 1:I) {
    for (j in 1:J.det) {
      for (k in 1:K.det) {
	for (t in 1:n.years) {
	  p.det[i, j, k, t] <- logit.inv(alpha.det.0[i] * X.det[j, k, t, 1] + 
					alpha.det.1[i] * X.det[j, k, t, 2])
          C[i, j, k, t] <- rbinom(1, 1, p.det[i, j, k, t] * z[i, pixel.det[j], t])
	} # t
      } # k
    } # j
  } # i

  # eBird Data Formation --------------------------------------------------
  # Number of pixels within each eBird cell
  n.eb.size <- J / J.eb
  # Lowest pixel number for each eBird cell
  low <- seq(from = 1, to = J, by = n.eb.size)
  # Highest pixel number for each eBird cell
  high <- seq(from = n.eb.size, to = J, by = n.eb.size)
  # Define detection/sampling process varyig over time and space
  p.eb <- array(NA, dim = c(I, J.eb, K.eb, n.years))
  # eBird data stored in y
  y <- array(NA, dim = c(I, J.eb, K.eb, n.years))
  z.eb <- array(NA, dim = c(I, J.eb, n.years))
  for (g in 1:J.eb) {
    for (i in 1:I) {
      for (t in 1:n.years) {
	# Change of support
        z.eb[i, g, t] <- ifelse(sum(z[i, low[g]:high[g], t]) > 0, 1, 0)
        for (k in 1:K.eb) {
	  # eBird detection
          p.eb[i, g, k, t] <- logit.inv(alpha.eb.0[i] * X.eb[g, k, t, 1] + 
				        alpha.eb.1[i] * X.eb[g, k, t, 2])
          y[i, g, k, t] <- rbinom(1, 1, z.eb[i, g, t] * p.eb[i, g, k, t])
        } # k
      } # t
    } # i
  } # g

  # NEON Data Formation ---------------------------------------------------
  p.neon <- array(NA, c(I, J.neon, K.neon, n.years))
  # x is the NEON data. Bad name, will change eventually. 
  x <- array(NA, dim = c(I, J.neon, K.neon, n.years))
  if (length(pixel.neon) == 0) {
    start.neon <- sample(1:J, 1)
    # Have NEON data be longer in the vertical direction than in the
    # horizontal direction
    n.cols.neon <- ceiling(J.neon / sqrt(J))
    vert.neon <- ifelse(J.neon > sqrt(J), sqrt(J), J.neon)
    vert.indices <- c(start.neon, rep(NA, vert.neon - 1))
    for (i in 2:vert.neon) {
      if (vert.indices[i - 1] - sqrt(J) > 0 & vert.indices[i - 1] <= start.neon) {
        vert.indices[i] <- vert.indices[i - 1] - sqrt(J)
      } else if (vert.indices[i - 1] < vert.indices[1]) {
          vert.indices[i] <- vert.indices[1] + sqrt(J)
        } else {
            vert.indices[i] <- vert.indices[i - 1] + sqrt(J)
          }
    }
    horz.neon <- J.neon - vert.neon
    pixel.neon <- c(vert.indices, rep(NA, horz.neon))
    if (J.neon > sqrt(J)) {
      for (i in (vert.neon + 1):J.neon) {
        if (((pixel.neon[i - vert.neon] - 1) %% sqrt(J) != 0) & (pixel.neon[i - vert.neon] <=
           ifelse((i %% sqrt(J)) == 0, vert.indices[sqrt(J)], vert.indices[i %% sqrt(J)])))	{
          pixel.neon[i] <- pixel.neon[i - vert.neon] - 1
        } else if (pixel.neon[i - vert.neon] > (vert.indices[ifelse((i %% sqrt(J)) == 0, sqrt(J),
            						      i %% sqrt(J))])) {
            pixel.neon[i] <- pixel.neon[i - vert.neon] + 1
          } else {
              pixel.neon[i] <- vert.indices[ifelse((i - vert.neon) %% sqrt(J) == 0, sqrt(J),
            				       i %% sqrt(J))] + 1
            }
      }
    }
  }
  for (j in 1:J.neon) {
    for (i in 1:I) {
      for (t in 1:n.years) {
        for (k in 1:K.neon) {
          p.neon[i, j, k, t] <- logit.inv(alpha.neon.0[i] * X.neon[j, k, t, 1] + 
					  alpha.neon.1[i] * X.neon[j, k, t, 2])
          x[i, j, k, t] <- rbinom(1, 1, z[i, pixel.neon[j], t] * p.neon[i, j, k, t])
        } # k
      } # t
    } # i
  } # j
  

  # BBS Data Formation ---------------------------------------------------
  p.bbs <- array(NA, c(I, J.bbs, n.years))
  # y.bbs is the BBS data
  y.bbs <- array(NA, dim = c(I, J.bbs, n.years))
  if (length(pixel.bbs) == 0) { 
    start.bbs <- sample(1:sqrt(J), n.route, replace = FALSE) * sqrt(J)
    pixel.bbs <- c(sapply(start.bbs, function(a) seq(a, a - sqrt(J) + 1, by = -1)))
  }
  for (j in 1:J.bbs) {
    for (i in 1:I) {
      for (t in 1:n.years) {
        p.bbs[i, j, t] <- logit.inv(alpha.bbs.0[i] * X.bbs[j, t, 1] + 
	      			  alpha.bbs.1[i] * X.bbs[j, t, 2])
        y.bbs[i, j, t] <- rbinom(1, 1, z[i, pixel.bbs[j], t] * p.bbs[i, j, t])
      } # t
    } # i
  } # j

  # Output all information needed for model. 
  return(list(X.psi = X.psi, X.phi = X.phi, X.gamma = X.gamma, X.bbs = X.bbs,
	      X.eb = X.eb, X.neon = X.neon, X.det = X.det, psi.0 = psi.0, 
	      z = z, phi = phi, gamma = gamma, pixel.det = pixel.det, 
	      p.det = p.det, C = C, low = low, high = high, p.eb = p.eb, 
	      y = y, pixel.neon = pixel.neon, p.neon = p.neon, x = x, 
	      y.bbs = y.bbs, pixel.bbs = pixel.bbs
	      ))
}


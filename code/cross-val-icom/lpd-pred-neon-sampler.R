# Predict NEON data using a model from a data set combination other than NEON. 
elpd.neon.pred <- function(samples.fit, samples.pred, elev.pred, for.pred, my.iter, 
		           n.iter, I, I.full, n.years) {

  param.fit.names <- attr(samples.fit, 'dimnames')[[2]]
  curr.indx <- which(substr(param.fit.names, 1, 9) %in% c('int.beta['))
  beta.0.fit <- logit(samples.fit[my.iter, curr.indx])
  beta.0.fit <- array(beta.0.fit, dim = c(n.iter, I.full, n.years))
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.1['))
  beta.1.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.2['))
  beta.2.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.3['))
  beta.3.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 4) %in% c('phi['))
  phi.fit <- samples.fit[my.iter, curr.indx]
  # Predicted values
  year.indices <- 6:9
  n.years.neon <- length(year.indices)
  sp.indices <- c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12)
  I.neon <- length(sp.indices)
  param.pred.names <- attr(samples.pred, 'dimnames')[[2]]
  curr.indx <- which(substr(param.pred.names, 1, 7) %in% c('z.neon['))
  z.pred.samples <- samples.pred[my.iter, curr.indx]
  J.pred <- ncol(z.pred.samples) / I / n.years.neon
  z.pred.samples <- array(z.pred.samples, dim = c(n.iter, I, J.pred, n.years.neon))
  z.fit.samples <- array(NA, dim = c(n.iter, I.neon, J.pred, n.years.neon))
  pred.dens.hbef <- array(NA, dim (z.fit.samples))

  # Composition sampling algorithm to get cross validation metrics ----------
  for (a in 1:n.iter) {
    print(paste("Currently on iteration ", a, " out of ", n.iter, sep = ''))
    for (i in 1:I.neon) {
      for (j in 1:J.pred) {
        for (t in 1:n.years.neon) {
          if (t == 1) {
            psi <- logit.inv(beta.0.fit[a, sp.indices[i], year.indices[t]] + 
  		           beta.1.fit[a, sp.indices[i]] * elev.pred[j] + 
  		           beta.2.fit[a, sp.indices[i]] * elev.pred[j]^2 + 
  		           beta.3.fit[a, sp.indices[i]] * for.pred[j])
            z.fit.samples[a, i, j, t] <- rbinom(1, 1, psi)
            pred.dens.hbef[a, i, j, t] <- dbinom(z.pred.samples[a, i, j, t], 1, psi)
          } else {
            psi <- logit.inv(beta.0.fit[a, sp.indices[i], year.indices[t]] + 
  		           beta.1.fit[a, sp.indices[i]] * elev.pred[j] + 
  		           beta.2.fit[a, sp.indices[i]] * elev.pred[j]^2 + 
  		           beta.3.fit[a, sp.indices[i]] * for.pred[j] + 
  			   phi.fit[a, sp.indices[i]] * z.fit.samples[a, i, j, t - 1])
            z.fit.samples[a, i, j, t] <- rbinom(1, 1, psi)
            pred.dens.hbef[a, i, j, t] <- dbinom(z.pred.samples[a, i, j, t], 1, psi)
  	}
        } # t
      } # j
    } # i
  }
  # Compute ELPD estimate
  return(apply(pred.dens.hbef, c(2, 3, 4), mean))
}

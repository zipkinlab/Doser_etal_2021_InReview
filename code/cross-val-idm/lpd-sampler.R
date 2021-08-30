elpd <- function(samples.fit, samples.pred, elev.pred, for.pred, my.iter, 
	         n.iter, n.years, type) {

  param.fit.names <- attr(samples.fit, 'dimnames')[[2]]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.0['))
  beta.0.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 6) %in% c('beta.1'))
  beta.1.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 6) %in% c('beta.2'))
  beta.2.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 6) %in% c('beta.3'))
  beta.3.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 3) %in% c('phi'))
  phi.fit <- samples.fit[my.iter, curr.indx]
  # Predicted values
  param.pred.names <- attr(samples.pred, 'dimnames')[[2]]
  val <- ifelse(type == 'hbef', 'z.hbef', ifelse(type == 'bbs', 'z.bbs[', 'z.neon'))
  year.indices <- 1:n.years
  if (type == 'neon') {
    year.indices <- 6:9
    n.years <- length(year.indices)
  }
  curr.indx <- which(substr(param.pred.names, 1, 6) %in% val)
  z.pred.samples <- samples.pred[my.iter, curr.indx]
  J.pred <- ncol(z.pred.samples) / n.years
  z.pred.samples <- array(z.pred.samples, dim = c(n.iter, J.pred, n.years))
  z.fit.samples <- array(NA, dim(z.pred.samples))
  pred.dens.hbef <- array(NA, dim (z.pred.samples))
  # Composition sampling algorithm to get cross validation metrics ----------
  for (a in 1:n.iter) {
    print(paste("Currently on iteration ", a, " out of ", n.iter, sep = ''))
    for (j in 1:J.pred) {
      for (t in 1:n.years) {
        if (t == 1) {
          psi <- logit.inv(beta.0.fit[a, year.indices[t]] + 
      	           beta.1.fit[a] * elev.pred[j] + 
      	           beta.2.fit[a] * elev.pred[j]^2 + 
      	           beta.3.fit[a] * for.pred[j])
          z.fit.samples[a, j, t] <- rbinom(1, 1, psi)
          pred.dens.hbef[a, j, t] <- dbinom(z.pred.samples[a, j, t], 1, psi)
        } else {
          psi <- logit.inv(beta.0.fit[a, year.indices[t]] + 
      	           beta.1.fit[a] * elev.pred[j] + 
      	           beta.2.fit[a] * elev.pred[j]^2 + 
      	           beta.3.fit[a] * for.pred[j] + 
      		   phi.fit[a] * z.fit.samples[a, j, t - 1])
          z.fit.samples[a, j, t] <- rbinom(1, 1, psi)
          pred.dens.hbef[a, j, t] <- dbinom(z.pred.samples[a, j, t], 1, psi)
      }
      } # t
    } # j
  }
  # Compute ELPD estimate
  return(apply(pred.dens.hbef, c(2, 3), mean))
}

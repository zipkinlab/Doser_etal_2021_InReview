elpd <- function(samples.fit, samples.pred, elev.pred, for.pred, my.iter, 
	         n.iter, I, n.years, type) {

  param.fit.names <- attr(samples.fit, 'dimnames')[[2]]
  curr.indx <- which(substr(param.fit.names, 1, 9) %in% c('int.beta['))
  beta.0.fit <- logit(samples.fit[my.iter, curr.indx])
  beta.0.fit <- array(beta.0.fit, dim = c(n.iter, I, n.years))
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.1['))
  beta.1.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.2['))
  beta.2.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.3['))
  beta.3.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 4) %in% c('phi['))
  phi.fit <- samples.fit[my.iter, curr.indx]
  # Predicted values
  param.pred.names <- attr(samples.pred, 'dimnames')[[2]]
  val <- ifelse(type == 'hbef', 'z.hbef', ifelse(type == 'bbs', 'z.bbs[', 'z.neon'))
  curr.indx <- which(substr(param.pred.names, 1, 6) %in% val)
  z.pred.samples <- samples.pred[my.iter, curr.indx]
  J.pred <- ncol(z.pred.samples) / I / n.years
  z.pred.samples <- array(z.pred.samples, dim = c(n.iter, I, J.pred, n.years))
  z.fit.samples <- array(NA, dim(z.pred.samples))
  pred.dens.hbef <- array(NA, dim (z.pred.samples))
  # Composition sampling algorithm to get cross validation metrics ----------
  for (a in 1:n.iter) {
    print(paste("Currently on iteration ", a, " out of ", n.iter, sep = ''))
    for (i in 1:I) {
      for (j in 1:J.pred) {
        for (t in 1:n.years) {
          if (t == 1) {
            psi <- logit.inv(beta.0.fit[a, i, t] + 
  		           beta.1.fit[a, i] * elev.pred[j] + 
  		           beta.2.fit[a, i] * elev.pred[j]^2 + 
  		           beta.3.fit[a, i] * for.pred[j])
            z.fit.samples[a, i, j, t] <- rbinom(1, 1, psi)
            pred.dens.hbef[a, i, j, t] <- dbinom(z.pred.samples[a, i, j, t], 1, psi)
          } else {
            psi <- logit.inv(beta.0.fit[a, i, t] + 
  		           beta.1.fit[a, i] * elev.pred[j] + 
  		           beta.2.fit[a, i] * elev.pred[j]^2 + 
  		           beta.3.fit[a, i] * for.pred[j] + 
  			   phi.fit[a, i] * z.fit.samples[a, i, j, t - 1])
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

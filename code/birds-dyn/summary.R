rm(list = ls())
library(jagsUI)
library(coda)
library(dplyr)
library(ggplot2)
library(ggthemes)
# For use in plots
load("../../data/bird-500m-3km-model-data.R")

# Functions ---------------------------------------------------------------
# Logit transformation
logit <- function(theta, a = 0, b = 1) {
  log((theta-a)/(b-theta))
}

# Inverse logit transformation
logit.inv <- function(z, a = 0, b = 1) {
  b-(b-a)/(1+exp(z))
}


# Read in results ---------------------------------------------------------
load("results/icm-small-bird-results-50000-iterations-2021-02-15.R")

sp <- attr(bugs.data$C, 'dimnames')[[3]]

# Species are: (1) American Redstart, (2) Black and white warbler, 
# (3) Blue-headed vireo, (4) Blackburnian warbler, (5) Blackpoll warbler, 
# (6) Black-throated blue, (7) Black throated green, (8) Canadian warbler, 
# (9) Magnolia warbler, (10) Nashville warbler, (11) Ovenbird, (12) red-eyed vireo

# Elevation Preference

# Low: BTBW, OVEN, REVI, BTNW, BAWW, CAWA, BLBW, AMRE, BHVI (although seems to range)
# High: BLPW, MAWA, NAWA


int.psi.samples <- out$sims.list$int.psi
colnames(int.psi.samples) <- sp
beta.psi.0.samples <- logit(int.psi.samples)
colnames(beta.psi.0.samples) <- sp
beta.psi.1.samples <- out$sims.list$beta.psi.1
colnames(beta.psi.1.samples) <- sp
beta.psi.2.samples <- out$sims.list$beta.psi.2
colnames(beta.psi.2.samples) <- sp
int.phi.samples <- out$sims.list$int.phi
colnames(int.phi.samples) <- sp
beta.phi.0.samples <- logit(int.phi.samples)
beta.phi.1.samples <- out$sims.list$beta.phi.1
colnames(beta.phi.1.samples) <- sp
beta.phi.2.samples <- out$sims.list$beta.phi.2
colnames(beta.phi.2.samples) <- sp
beta.phi.3.samples <- out$sims.list$beta.phi.3
colnames(beta.phi.3.samples) <- sp
beta.phi.4.samples <- out$sims.list$beta.phi.4
colnames(beta.phi.4.samples) <- sp
beta.phi.5.samples <- out$sims.list$beta.phi.5
colnames(beta.phi.5.samples) <- sp
int.gamma.samples <- out$sims.list$int.gamma
colnames(int.gamma.samples) <- sp
beta.gamma.0.samples <- logit(int.gamma.samples)
beta.gamma.1.samples <- out$sims.list$beta.gamma.1
colnames(beta.gamma.1.samples) <- sp
beta.gamma.2.samples <- out$sims.list$beta.gamma.2
colnames(beta.gamma.2.samples) <- sp
beta.gamma.3.samples <- out$sims.list$beta.gamma.3
colnames(beta.gamma.3.samples) <- sp
beta.gamma.4.samples <- out$sims.list$beta.gamma.4
colnames(beta.gamma.4.samples) <- sp
beta.gamma.5.samples <- out$sims.list$beta.gamma.5
colnames(beta.gamma.5.samples) <- sp
rhat.vals <- out$Rhat


n.iter <- nrow(beta.psi.1.samples)

# Community-wide summarys
summary(mcmc(out$sims.list$int.psi.mean))$quantiles # Intercept
summary(mcmc(out$sims.list$beta.psi.1.mean))$quantiles # Linear Elevation
mean(out$sims.list$beta.psi.1.mean < 0)
# Strong support for an overall negative quadratic relationship with elevation
summary(mcmc(out$sims.list$beta.psi.2.mean))$quantiles # Squared Elevation
mean(out$sims.list$beta.psi.2.mean < 0)
summary(mcmc(out$sims.list$int.phi.mean))$quantiles # Persistence intercept
# Moderate support for neg TMEAN effect
summary(mcmc(out$sims.list$beta.phi.1.mean))$quantiles # TMEAN persistence
mean(out$sims.list$beta.phi.1.mean < 0)
# Moderate support for neg PPT effect
summary(mcmc(out$sims.list$beta.phi.2.mean))$quantiles # PPT persistence
mean(out$sims.list$beta.phi.2.mean < 0) 
summary(mcmc(out$sims.list$beta.phi.3.mean))$quantiles # ELEV persistence
summary(mcmc(out$sims.list$beta.phi.4.mean))$quantiles # TMEAN * ELEV persistence
# Fairly large support for positive interaction of temperature and elevation
# e.g., higher elevation, effect of TMEAN switches to positive. 
mean(out$sims.list$beta.phi.4.mean > 0)
summary(mcmc(out$sims.list$beta.phi.5.mean))$quantiles # PPT * ELEV persistence
# Moderate support for neg interaction. e.g., effect of PPT becomes more neg at 
# higher elevations
mean(out$sims.list$beta.phi.5.mean < 0)
# Next to no positive effect
summary(mcmc(out$sims.list$beta.gamma.1.mean))$quantiles # TMEAN colonization
mean(out$sims.list$beta.gamma.1.mean > 0)
# Moderate support for positive PPT effect
summary(mcmc(out$sims.list$beta.gamma.2.mean))$quantiles # PPT colonization
mean(out$sims.list$beta.gamma.2.mean > 0)
summary(mcmc(out$sims.list$beta.gamma.3.mean))$quantiles # ELEV colonization
# No interaction between temperature and elevation for colonization
summary(mcmc(out$sims.list$beta.gamma.4.mean))$quantiles # TMEAN * ELEV colonization
mean(out$sims.list$beta.gamma.4.mean > 0)
# No interaction between ppt and elevation for colonization
summary(mcmc(out$sims.list$beta.gamma.5.mean))$quantiles # PPT * ELEV colonization
mean(out$sims.list$beta.gamma.5.mean > 0)

# Persistence
summary(mcmc(beta.phi.1.samples))$quantiles # TMEAN
plot(mcmc(beta.phi.1.samples), density = FALSE)
# Prob of a negative effect of temperature
apply(beta.phi.1.samples, 2, function(a) mean(a < 0))
summary(mcmc(out$sims.list$beta.phi.2))$quantiles # PPT
# Prob of a negative effect of precipitation
apply(out$sims.list$beta.phi.2, 2, function(a) mean(a < 0))
summary(mcmc(out$sims.list$beta.phi.3))$quantiles  # Elevation
summary(mcmc(out$sims.list$beta.phi.4))$quantiles # TMEAN * Elevation
apply(out$sims.list$beta.phi.4, 2, function(a) mean(a > 0))
summary(mcmc(out$sims.list$beta.phi.5))$quantiles # PPT * Elevation

# Colonization
# Very little evidence of temperature or precipitation effects on 
# colonization. 
summary(mcmc(out$sims.list$beta.gamma.1))$quantiles # TMEAN
apply(out$sims.list$beta.gamma.1, 2, function(a) mean(a < 0))
summary(mcmc(out$sims.list$beta.gamma.2))$quantiles # PPT
apply(out$sims.list$beta.gamma.2, 2, function(a) mean(a > 0))
summary(mcmc(out$sims.list$beta.gamma.3))$quantiles # Elevation
summary(mcmc(out$sims.list$beta.gamma.4))$quantiles
summary(mcmc(out$sims.list$beta.gamma.5))$quantiles

# Mean colonization across years
summary(mcmc(out$sims.list$int.gamma[, 1, ]))$quantiles
summary(mcmc(out$sims.list$int.gamma[, 2, ]))$quantiles
summary(mcmc(out$sims.list$int.gamma[, 3, ]))$quantiles


# Detection covariates
# eBird
# All of these seemingly make sense. 
summary(mcmc(out$sims.list$int.alpha.eb))$quantiles # Intercept
summary(mcmc(out$sims.list$alpha.eb.1))$quantiles # Observers
summary(mcmc(out$sims.list$alpha.eb.2))$quantiles # Distance
summary(mcmc(out$sims.list$alpha.eb.3))$quantiles # Day
summary(mcmc(out$sims.list$alpha.eb.4))$quantiles # Day^2
summary(mcmc(out$sims.list$alpha.eb.5))$quantiles # Time
summary(mcmc(out$sims.list$alpha.eb.6))$quantiles # Length

# Hubbard Brook
summary(mcmc(out$sims.list$int.alpha.hb))$quantiles # Intercept
summary(mcmc(out$sims.list$alpha.hb.1))$quantiles # Day
summary(mcmc(out$sims.list$alpha.hb.2))$quantiles # Day^2
summary(mcmc(out$sims.list$alpha.hb.3))$quantiles # TOD

# NEON
summary(mcmc(out$sims.list$int.alpha.neon))$quantiles # Intercept
summary(mcmc(out$sims.list$alpha.neon.1))$quantiles # Day
summary(mcmc(out$sims.list$alpha.neon.2))$quantiles # Day^2
summary(mcmc(out$sims.list$alpha.neon.3))$quantiles # TOD

# Vague look at overdispersion in raw data, will eventually need 
# to do a BPV on some binned data....
var(bugs.data$y, na.rm = TRUE)
mean(bugs.data$y, na.rm = TRUE)
mean(bugs.data$C, na.rm = TRUE)
var(bugs.data$C, na.rm = TRUE)
mean(bugs.data$x, na.rm = TRUE)
var(bugs.data$x, na.rm = TRUE)

# Heat map of regression coefficients -------------------------------------
temp.quants <- as.data.frame(summary(mcmc(beta.phi.1.samples))$quantiles)
temp.quants$`10%` <- apply(mcmc(beta.phi.1.samples), 2, quantile, prob = 0.10)
temp.quants$`90%` <- apply(mcmc(beta.phi.1.samples), 2, quantile, prob = 0.90)
elev.quants <- as.data.frame(summary(mcmc(beta.phi.3.samples))$quantiles)
elev.quants$`10%` <- apply(mcmc(beta.phi.3.samples), 2, quantile, prob = 0.10)
elev.quants$`90%` <- apply(mcmc(beta.phi.3.samples), 2, quantile, prob = 0.90)
elev.temp.quants <- as.data.frame(summary(mcmc(beta.phi.4.samples))$quantiles)
elev.temp.quants$`10%` <- apply(mcmc(beta.phi.4.samples), 2, quantile, prob = 0.10)
elev.temp.quants$`90%` <- apply(mcmc(beta.phi.4.samples), 2, quantile, prob = 0.90)
elev.ppt.quants <- as.data.frame(summary(mcmc(beta.phi.5.samples))$quantiles)
elev.ppt.quants$`10%` <- apply(mcmc(beta.phi.5.samples), 2, quantile, prob = 0.10)
elev.ppt.quants$`90%` <- apply(mcmc(beta.phi.5.samples), 2, quantile, prob = 0.90)
ppt.quants <- as.data.frame(summary(mcmc(beta.phi.2.samples))$quantiles)
ppt.quants$`10%` <- apply(mcmc(beta.phi.2.samples), 2, quantile, prob = 0.10)
ppt.quants$`90%` <- apply(mcmc(beta.phi.2.samples), 2, quantile, prob = 0.90)
temp.col.quants <- as.data.frame(summary(mcmc(beta.gamma.1.samples))$quantiles)
temp.col.quants$`10%` <- apply(mcmc(beta.gamma.1.samples), 2, quantile, prob = 0.10)
temp.col.quants$`90%` <- apply(mcmc(beta.gamma.1.samples), 2, quantile, prob = 0.90)
elev.col.quants <- as.data.frame(summary(mcmc(beta.gamma.3.samples))$quantiles)
elev.col.quants$`10%` <- apply(mcmc(beta.gamma.3.samples), 2, quantile, prob = 0.10)
elev.col.quants$`90%` <- apply(mcmc(beta.gamma.3.samples), 2, quantile, prob = 0.90)
ppt.col.quants <- as.data.frame(summary(mcmc(beta.gamma.2.samples))$quantiles)
ppt.col.quants$`10%` <- apply(mcmc(beta.gamma.2.samples), 2, quantile, prob = 0.10)
ppt.col.quants$`90%` <- apply(mcmc(beta.gamma.2.samples), 2, quantile, prob = 0.90)
elev.col.temp.quants <- as.data.frame(summary(mcmc(beta.gamma.4.samples))$quantiles)
elev.col.temp.quants$`10%` <- apply(mcmc(beta.gamma.4.samples), 2, quantile, prob = 0.10)
elev.col.temp.quants$`90%` <- apply(mcmc(beta.gamma.4.samples), 2, quantile, prob = 0.90)
elev.col.ppt.quants <- as.data.frame(summary(mcmc(beta.gamma.5.samples))$quantiles)
elev.col.ppt.quants$`10%` <- apply(mcmc(beta.gamma.5.samples), 2, quantile, prob = 0.10)
elev.col.ppt.quants$`90%` <- apply(mcmc(beta.gamma.5.samples), 2, quantile, prob = 0.90)
# Probs of negative effect
temp.probs <- apply(beta.phi.1.samples, 2, function(a) mean(a < 0))
ppt.probs <- apply(beta.phi.2.samples, 2, function(a) mean(a < 0))
elev.probs <- apply(beta.phi.3.samples, 2, function(a) mean(a < 0))
elev.temp.probs <- apply(beta.phi.4.samples, 2, function(a) mean(a < 0))
elev.ppt.probs <- apply(beta.phi.5.samples, 2, function(a) mean(a < 0))
temp.col.probs <- apply(beta.gamma.1.samples, 2, function(a) mean(a < 0))
ppt.col.probs <- apply(beta.gamma.2.samples, 2, function(a) mean(a < 0))
elev.col.probs <- apply(beta.gamma.3.samples, 2, function(a) mean(a < 0))
elev.temp.col.probs <- apply(beta.gamma.4.samples, 2, function(a) mean(a < 0))
elev.ppt.col.probs <- apply(beta.gamma.5.samples, 2, function(a) mean(a < 0))

probs <- round(c(temp.probs, ppt.probs, elev.probs, elev.temp.probs, elev.ppt.probs, 
	   temp.col.probs, ppt.col.probs, elev.col.probs, elev.temp.col.probs, 
	   elev.ppt.col.probs), 2)

colnames(temp.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(ppt.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(elev.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(elev.temp.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(elev.ppt.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(temp.col.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(ppt.col.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(elev.col.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(elev.col.temp.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')
colnames(elev.col.ppt.quants) <- c('lowest', 'low', 'med', 'high', 'highest', 'lower', 'higher')

climate.vals <- as.data.frame(rbind(temp.quants, ppt.quants, elev.quants, elev.temp.quants, 
				    elev.ppt.quants, 
				    temp.col.quants, ppt.col.quants, 
				    elev.col.quants, elev.col.temp.quants, 
				    elev.col.ppt.quants))
climate.vals$param <- factor(rep(c('T', 'P', 'E', 'TE', 'PE', 'T1', 'P1', 'E1', 'TE1', 'PE1'), 
				 each = bugs.data$K),
			     levels = c('T', 'P', 'E', 'TE', 'PE', 
					'T1', 'P1', 'E1', 'TE1', 'PE1'))
rownames(climate.vals) <- NULL
climate.vals$sp <- rep(factor(sp), length(levels(climate.vals$param)))
climate.vals <- climate.vals %>%
  mutate(sig = ifelse(((lowest < 0) & (highest < 0)) | ((lowest > 0) & (highest > 0)), '**',
		      ifelse(((lower < 0) & (higher < 0)) | ((lower > 0) & (higher > 0)), '*', '')))
climate.vals$probs <- probs

heat.map <- ggplot(climate.vals, aes(x = param, y = sp, fill = med)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = 0, low = 'red', mid = 'white', high = 'blue',
	               na.value = 'gray', limits = c(-4, 4)) +
  geom_text(aes(label = sig)) +
  geom_vline(xintercept = 5.5, col = 'black') +
  theme_bw(base_size = 18) +
  labs(x = '', y = 'Species', fill = 'Effect Size') +
  annotate('text', x = 3, y = -1.6, label = 'Persistence', size = 5) +
  annotate('text', x = 8, y = -1.6, label = 'Colonization', size = 5) +
  coord_cartesian(ylim = c(1, 12), clip = 'off') +
  scale_x_discrete(labels = c('T' = 'Temp', 'P' = 'Ppt', 'E' = 'Elev', 
			      'T1' = 'Temp', 'P1' = 'Ppt', 'E1' = 'Elev',
			      'TE' = 'Temp x Elev', 'PE' = 'Ppt x Elev', 
			      'TE1' = 'Temp x Elev', 'PE1' = 'PPt x Elev')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# scale_y_discrete(labels = c('AMRE' = 'American Redstart',
#      		      'BAWW' = 'Black and White Warbler',
#      		      'BHVI' = 'Blue-headed Vireo',
#      		      'BLBW' = 'Blackburnian Warbler',
#      		      'BLPW' = 'Blackpoll Warbler',
#      		      'BTBW' = 'Black-throated Blue Warbler',
#      		      'BTNW' = 'Black-throated Green Warbler',
#      		      'CAWA' = 'Canada Warbler',
#      		      'MAWA' = 'Magnolia Warbler',
#      		      'NAWA' = 'Nashville Warbler',
#      		      'OVEN' = 'Ovenbird',
#      		      'REVI' = 'Red-eyed Vireo'))
ggsave('figures/coefficient-heat-map.png', heat.map, width = 11, height = 9, 
       units = 'in')

# Temprature effect on persistence ----------------------------------------
n.vals <- 500
phi.1.med <- array(NA, dim = c(n.vals, bugs.data$K))
phi.1.low <- array(NA, dim = c(n.vals, bugs.data$K))
phi.1.high <- array(NA, dim = c(n.vals, bugs.data$K))
phi.2.med <- array(NA, dim = c(n.vals, bugs.data$K))
phi.2.low <- array(NA, dim = c(n.vals, bugs.data$K))
phi.2.high <- array(NA, dim = c(n.vals, bugs.data$K))
phi.3.med <- array(NA, dim = c(n.vals, bugs.data$K))
phi.3.low <- array(NA, dim = c(n.vals, bugs.data$K))
phi.3.high <- array(NA, dim = c(n.vals, bugs.data$K))
#tmean.true.vals <- seq(from = min(tmean.vals, na.rm = TRUE), 
#		       to = max(tmean.vals, na.rm = TRUE), length.out = n.vals)
#tmean.plot.vals <- (tmean.true.vals - mean(tmean.vals, na.rm = TRUE)) / 
#	sd(tmean.vals, na.rm = TRUE)
tmean.plot.vals <- seq(from = min(bugs.data$TMEAN, na.rm = TRUE),
		 to = max(bugs.data$TMEAN, na.rm = TRUE), length.out = n.vals)
elev.vals <- c(min(bugs.data$ELEV, na.rm = TRUE), mean(bugs.data$ELEV, na.rm = TRUE),
	       max(bugs.data$ELEV, na.rm = TRUE))
for (k in 1:bugs.data$K) {
  print(k) 
  for (a in 1:n.vals) {
    phi.1.med[a, k] <- mean(logit.inv(beta.phi.0.samples[, k, 5] + 
				    beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				    beta.phi.3.samples[, k] * elev.vals[1] + 
				    beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[1]))
    phi.1.low[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				        beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				        beta.phi.3.samples[, k] * elev.vals[1] + 
				        beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[1]), 
			      prob = 0.025)
    phi.1.high[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				         beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				         beta.phi.3.samples[, k] * elev.vals[1] + 
				         beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[1]), 
			      prob = 0.975)
    phi.2.med[a, k] <- mean(logit.inv(beta.phi.0.samples[, k, 5] + 
				    beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				    beta.phi.3.samples[, k] * elev.vals[2] + 
				    beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[2]))
    phi.2.low[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				        beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				        beta.phi.3.samples[, k] * elev.vals[2] + 
				        beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[2]), 
			      prob = 0.025)
    phi.2.high[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				         beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				         beta.phi.3.samples[, k] * elev.vals[2] + 
				         beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[2]), 
			      prob = 0.975)
    phi.3.med[a, k] <- mean(logit.inv(beta.phi.0.samples[, k, 5] + 
				    beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				    beta.phi.3.samples[, k] * elev.vals[3] + 
				    beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[3]))
    phi.3.low[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				        beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				        beta.phi.3.samples[, k] * elev.vals[3] + 
				        beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[3]), 
			      prob = 0.025)
    phi.3.high[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				         beta.phi.1.samples[, k] * tmean.plot.vals[a] + 
				         beta.phi.3.samples[, k] * elev.vals[3] + 
				         beta.phi.4.samples[, k] * tmean.plot.vals[a] * elev.vals[3]), 
			      prob = 0.975)
  } # a
} # k

meds <- c(phi.1.med, phi.2.med, phi.3.med)
lows <- c(phi.1.low, phi.2.low, phi.3.low)
highs <- c(phi.1.high, phi.2.high, phi.3.high)
tmean.plot.dat <- data.frame(phi = meds, 
			    sp = rep(rep(sp, each = n.vals), times = length(elev.vals)), 
			    tmean = rep(tmean.plot.vals, times = bugs.data$K * length(elev.vals)), 
			    low = lows, 
			    high = highs, 
			    elev = factor(rep(c('Low (130m)', 'Mid (587m)', 'High (1856m)'), 
					      each = bugs.data$K * n.vals), 
					  levels = c('Low (130m)', 'Mid (587m)', 'High (1856m)')))

my.colors <- c('Low (130m)' = 'red', 'Mid (587m)' = 'gray', 'High (1856m)' = 'blue')
tmean.plot <- ggplot(tmean.plot.dat, aes(x = tmean, y = phi, col = elev)) + 
  geom_line(size = 1.5) + 
  theme_bw(base_size = 18) + 
  facet_wrap(vars(sp)) + 
  labs(x = 'Temperature Deviation', y = 'Persistence', col = 'Elevation') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_color_manual(values = my.colors)

ggsave('figures/persistence-temp.png', width = 11, height = 7, units = 'in')

# PPT effect on persistence ----------------------------------------
n.vals <- 500
phi.1.med <- array(NA, dim = c(n.vals, bugs.data$K))
phi.1.low <- array(NA, dim = c(n.vals, bugs.data$K))
phi.1.high <- array(NA, dim = c(n.vals, bugs.data$K))
phi.2.med <- array(NA, dim = c(n.vals, bugs.data$K))
phi.2.low <- array(NA, dim = c(n.vals, bugs.data$K))
phi.2.high <- array(NA, dim = c(n.vals, bugs.data$K))
phi.3.med <- array(NA, dim = c(n.vals, bugs.data$K))
phi.3.low <- array(NA, dim = c(n.vals, bugs.data$K))
phi.3.high <- array(NA, dim = c(n.vals, bugs.data$K))
#ppt.true.vals <- seq(from = min(ppt.vals, na.rm = TRUE), 
#		       to = max(ppt.vals, na.rm = TRUE), length.out = n.vals)
#ppt.plot.vals <- (ppt.true.vals - mean(ppt.vals, na.rm = TRUE)) / 
#	sd(ppt.vals, na.rm = TRUE)
ppt.plot.vals <- seq(from = min(bugs.data$PPT, na.rm = TRUE),
		 to = max(bugs.data$PPT, na.rm = TRUE), length.out = n.vals)
elev.vals <- c(min(bugs.data$ELEV, na.rm = TRUE), mean(bugs.data$ELEV, na.rm = TRUE),
	       max(bugs.data$ELEV, na.rm = TRUE))
for (k in 1:bugs.data$K) {
  print(k) 
  for (a in 1:n.vals) {
    phi.1.med[a, k] <- mean(logit.inv(beta.phi.0.samples[, k, 5] + 
				    beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				    beta.phi.3.samples[, k] * elev.vals[1] + 
				    beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[1]))
    phi.1.low[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				        beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				        beta.phi.3.samples[, k] * elev.vals[1] + 
				        beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[1]), 
			      prob = 0.025)
    phi.1.high[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				         beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				         beta.phi.3.samples[, k] * elev.vals[1] + 
				         beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[1]), 
			      prob = 0.975)
    phi.2.med[a, k] <- mean(logit.inv(beta.phi.0.samples[, k, 5] + 
				    beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				    beta.phi.3.samples[, k] * elev.vals[2] + 
				    beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[2]))
    phi.2.low[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				        beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				        beta.phi.3.samples[, k] * elev.vals[2] + 
				        beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[2]), 
			      prob = 0.025)
    phi.2.high[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				         beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				         beta.phi.3.samples[, k] * elev.vals[2] + 
				         beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[2]), 
			      prob = 0.975)
    phi.3.med[a, k] <- mean(logit.inv(beta.phi.0.samples[, k, 5] + 
				    beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				    beta.phi.3.samples[, k] * elev.vals[3] + 
				    beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[3]))
    phi.3.low[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				        beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				        beta.phi.3.samples[, k] * elev.vals[3] + 
				        beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[3]), 
			      prob = 0.025)
    phi.3.high[a, k] <- quantile(logit.inv(beta.phi.0.samples[, k, 5] + 
				         beta.phi.2.samples[, k] * ppt.plot.vals[a] + 
				         beta.phi.3.samples[, k] * elev.vals[3] + 
				         beta.phi.5.samples[, k] * ppt.plot.vals[a] * elev.vals[3]), 
			      prob = 0.975)
  } # a
} # k

meds <- c(phi.1.med, phi.2.med, phi.3.med)
lows <- c(phi.1.low, phi.2.low, phi.3.low)
highs <- c(phi.1.high, phi.2.high, phi.3.high)
ppt.plot.dat <- data.frame(phi = meds, 
			    sp = rep(rep(sp, each = n.vals), times = length(elev.vals)), 
			    ppt = rep(ppt.plot.vals, times = bugs.data$K * length(elev.vals)), 
			    low = lows, 
			    high = highs, 
			    elev = factor(rep(c('Low (130m)', 'Mid (587m)', 'High (1856m)'), 
					      each = bugs.data$K * n.vals), 
					  levels = c('Low (130m)', 'Mid (587m)', 'High (1856m)')))

my.colors <- c('Low (130m)' = 'red', 'Mid (587m)' = 'gray', 'High (1856m)' = 'blue')
ppt.plot <- ggplot(ppt.plot.dat, aes(x = ppt, y = phi, col = elev)) + 
  geom_line(size = 1.5) + 
  theme_bw(base_size = 18) + 
  facet_wrap(vars(sp)) + 
  labs(x = 'Precipitation Deviation', y = 'Persistence', col = 'Elevation') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_color_manual(values = my.colors)

ggsave('figures/persistence-ppt.png', width = 11, height = 7, units = 'in')

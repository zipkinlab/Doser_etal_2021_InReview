# summary.R: file to summarize additional results from the full ICOM 
#            for the foliage-gleaning bird case study. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggthemes)


# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# Read in data ------------------------------------------------------------
load("data/nimble-data.R")

# Read in results ---------------------------------------------------------
load("results/icom-HBEF-NEON-BBS-results.R")
param.names <- attr(samples[[1]], 'dimnames')[[2]]
# Species names and four letter codes
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW',
        'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
sp.names <- c('American Redstart', 'Black-and-white Warbler', 'Blue-headed Vireo', 
	      'Blackburnian Warbler', 'Blackpoll Warbler', 'Black-throated Blue Warbler', 
	      'Black-throated Green Warbler', 'Canada Warbler', 'Magnolia Warbler', 
	      'Nashville Warbler', 'Ovenbird', 'Red-eyed Vireo')

# Summarize Data for Supplemental Table -----------------------------------
# HBEF total observations
y.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count, na.rm = TRUE))
# HBEF site/year combinations
y.df %>%
  group_by(Species, Site, Year) %>%
  summarize(Count = max(Count, na.rm = TRUE)) %>%
  group_by(Species) %>%
  summarize(counts = sum(Count, na.rm = TRUE))
# NEON
v.1.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count, na.rm = TRUE))
# BBS
v.2.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count, na.rm = TRUE)) 

# Parameter summaries reported in Table S3
curr.indx <- which(substr(param.names, 1, 9) == 'int.beta[')
# Occurrence intercept ----------------
int.beta.samples <- samples[, curr.indx]
summary(int.beta.samples)$quantiles
#plot(int.beta.samples, density = FALSE)
beta.0.samples <- mcmc.list(lapply(int.beta.samples, logit))

# Linear elevation --------------------
curr.indx <- which(substr(param.names, 1, 7) == 'beta.1[')
beta.1.samples <- samples[, curr.indx]
summary(beta.1.samples)$quantiles
#plot(beta.1.samples, density = FALSE)
# Quadratic elevation -----------------
curr.indx <- which(substr(param.names, 1, 7) == 'beta.2[')
beta.2.samples <- samples[, curr.indx]
summary(beta.2.samples)
# plot(beta.2.samples, density = FALSE)
# Percent forest ----------------------
curr.indx <- which(substr(param.names, 1, 7) == 'beta.3[')
beta.3.samples <- samples[, curr.indx]
summary(beta.3.samples)
# plot(beta.3.samples, density = FALSE)
# Autologistic parameters -------------
curr.indx <- which(substr(param.names, 1, 4) == 'phi[')
phi.samples <- samples[, curr.indx]
summary(phi.samples)
#plot(phi.samples, density = FALSE)

# Post-hoc trend analysis and figure --------------------------------------
trend.samples <- beta.samples[, 2, ]
trend.comm.samples <- apply(beta.samples[, 2, ], 2, mean)
apply(trend.samples, 1, mean)
apply(trend.samples, 1, quantile, prob = c(0.025, 0.975))
apply(trend.samples, 1, sd)
mean(trend.comm.samples)
sd(trend.comm.samples)
quantile(trend.comm.samples, prob = c(0.025, 0.975))

trend.probs <- apply(trend.samples, 1, function(a) sum(a > 0) / length(a))
trend.means <- apply(trend.samples, 1, mean)
trend.df <- data.frame(prob = trend.probs, 
		       val = trend.means, 
		       sp = sp.names)
trend.df <- trend.df %>%
  mutate(sp = factor(sp, levels = sp[order(prob)]))
# Figure S3
trend.df %>%
  ggplot(aes(x = prob, y = sp, col = val)) + 
  geom_point(size = 4.5) + 
  theme_bw(base_size = 20) + 
  scale_color_gradient2(midpoint = 0, low = 'red', mid = 'gray', high = 'blue', 
			na.value = 'grey') + 
  geom_vline(xintercept = 0.75, lty = 2) + 
  geom_vline(xintercept = 0.25, lty = 2) + 
  scale_x_continuous(limits = c(0, 1)) + 
  labs(col = 'Trend', y = 'Species', x = 'Probability of Increasing Trend')
  ggsave(device = 'pdf', filename = 'figures/FigS3.pdf', 
         height = 8, width = 10, units = 'in')


# Elevation Plots ---------------------------------------------------------
n.vals <- 500
psi.med <- array(NA, dim = c(n.vals, icom.consts$I))
psi.low <- array(NA, dim = c(n.vals, icom.consts$I))
psi.high <- array(NA, dim = c(n.vals, icom.consts$I))
elev.true.vals <- seq(from = min(elev.real, na.rm = TRUE), 
		      to = max(elev.real, na.rm = TRUE), length.out = n.vals)
elev.vals <- (elev.true.vals - mean(elev.real, na.rm = TRUE)) / 
	sd(elev.real, na.rm = TRUE)
beta.0.samples.all <- do.call('rbind', beta.0.samples[, 13:24])
beta.1.samples.all <- do.call('rbind', beta.1.samples)
beta.2.samples.all <- do.call('rbind', beta.2.samples)
for (i in 1:icom.consts$I) {
  print(i) 
  for (a in 1:n.vals) {
    psi.med[a, i] <- median(logit.inv(beta.0.samples.all[, i] + 
				      beta.1.samples.all[, i] * elev.vals[a] + 
				      beta.2.samples.all[, i] * elev.vals[a]^2))
    psi.low[a, i] <- quantile(logit.inv(beta.0.samples.all[, i] + 
					beta.1.samples.all[, i] * elev.vals[a] + 
					beta.2.samples.all[, i] * elev.vals[a]^2), 
			      prob = 0.025)
    psi.high[a, i] <- quantile(logit.inv(beta.0.samples.all[, i] + 
					 beta.1.samples.all[, i] * elev.vals[a] + 
					 beta.2.samples.all[, i] * elev.vals[a]^2), 
			      prob = 0.975)
  } # a
} # k

elev.plot.dat <- data.frame(vals = c(psi.med), 
			    sp = rep(sp.names, each = n.vals), 
			    elev = elev.true.vals, 
			    low = c(psi.low), 
			    high = c(psi.high)) 
# Figure S4
elev.plot <- ggplot(elev.plot.dat, aes(x = elev, y = vals)) + 
  geom_line(size = 1.5) + 
  theme_bw(base_size = 22) + 
  facet_wrap(vars(sp)) + 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.25, col = 'grey') + 
  labs(x = 'Elevation (m)', y = 'Occurrence Probability')
elev.plot
ggsave(elev.plot, device = 'pdf', filename = 'figures/FigS4.pdf', 
       height = 8, width = 16, units = 'in')

# Richness by Elevation Plot ----------------------------------------------
hb.rich.mean <- apply(rich.hbef.samples, 2, mean)
hb.rich.low <- apply(rich.hbef.samples, 2, quantile, 0.025)
hb.rich.high <- apply(rich.hbef.samples, 2, quantile, 0.975)
neon.rich.mean <- apply(rich.neon.samples, 2, mean)
neon.rich.low <- apply(rich.neon.samples, 2, quantile, 0.025)
neon.rich.high <- apply(rich.neon.samples, 2, quantile, 0.975)
bbs.rich.mean <- apply(rich.bbs.samples, 2, mean)
bbs.rich.low <- apply(rich.bbs.samples, 2, quantile, 0.025)
bbs.rich.high <- apply(rich.bbs.samples, 2, quantile, 0.975)
rich.df <- data.frame(med = c(hb.rich.mean, neon.rich.mean, bbs.rich.mean), 
		      low = c(hb.rich.low, neon.rich.low, bbs.rich.low), 
		      high = c(hb.rich.high, neon.rich.high, bbs.rich.high), 
		      dataset = c(rep('HBEF', length(hb.rich.mean)),
				  rep('NEON', length(neon.rich.mean)), 
				  rep('BBS', length(bbs.rich.mean))), 
		      elev = elev.real) 
# Jaccard by Elevation Plot -----------------------------------------------
hb.jaccard.mean <- apply(jaccard.hbef.samples, 2, mean)
hb.jaccard.low <- apply(jaccard.hbef.samples, 2, quantile, 0.025)
hb.jaccard.high <- apply(jaccard.hbef.samples, 2, quantile, 0.975)
neon.jaccard.mean <- apply(jaccard.neon.samples, 2, mean)
neon.jaccard.low <- apply(jaccard.neon.samples, 2, quantile, 0.025)
neon.jaccard.high <- apply(jaccard.neon.samples, 2, quantile, 0.975)
bbs.jaccard.mean <- apply(jaccard.bbs.samples, 2, mean)
bbs.jaccard.low <- apply(jaccard.bbs.samples, 2, quantile, 0.025)
bbs.jaccard.high <- apply(jaccard.bbs.samples, 2, quantile, 0.975)
jaccard.df <- data.frame(med = c(hb.jaccard.mean, neon.jaccard.mean, bbs.jaccard.mean), 
		      low = c(hb.jaccard.low, neon.jaccard.low, bbs.jaccard.low), 
		      high = c(hb.jaccard.high, neon.jaccard.high, bbs.jaccard.high), 
		      dataset = c(rep('HBEF', length(hb.jaccard.mean)),
				  rep('NEON', length(neon.jaccard.mean)), 
				  rep('BBS', length(bbs.jaccard.mean))), 
		      elev = elev.real) 
rich.plot <- rich.df %>%
  ggplot(aes(x = elev, y = med, col = dataset)) + 
  geom_point(size = 2) + 
  theme_bw(base_size = 24) + 
  labs(x = 'Elevation (m)', y = 'Average Richness', col = 'Data Set') + 
  scale_color_colorblind()
jaccard.plot <- jaccard.df %>%
  ggplot(aes(x = elev, y = med, col = dataset)) + 
  geom_point(size = 2) + 
  theme_bw(base_size = 24) + 
  labs(x = 'Elevation (m)', y = 'Average Jaccard Index', col = 'Data Set') + 
  scale_color_colorblind()

# Figure 3
ggarrange(rich.plot, jaccard.plot, nrow = 1, ncol = 2, 
	  common.legend = TRUE, labels = c('A', 'B'), 
	  font.label = list(size = 24), legend = 'bottom')
ggsave(device = 'pdf', filename = 'figures/Fig3.pdf', 
       height = 8, width = 14, units = 'in')

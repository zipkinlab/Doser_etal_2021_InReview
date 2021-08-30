# summary-model-comparison.R: file to compare results of seven models that
#                             use different amounts of HBEF, NEON, and BBS
#                             data for the foliage-gleaning bird case study. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(tidyverse)
library(ggthemes)
library(ggpubr)

# Load in formatted data used for model fitting process. 
load("data/nimble-data.R")
# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# Read in results ---------------------------------------------------------
# Model HBEF --------------------------
load("results/icom-HBEF-results.R")
# Rename and save objects of interest from each file. 
# Parameter samples 
samples.1 <- samples
# Names of parameters
param.names.1 <- attr(samples.1[[1]], 'dimnames')[[2]]
# Jaccard Index
jac.hbef.samples.1 <- jaccard.hbef.samples
# Richness 
rich.hbef.samples.1 <- rich.hbef.samples
# Average occurrence each year
psi.samples.1 <- psi.avg.sp.year
# Post-hoc trend coefficients
trend.samples.1 <- beta.samples[, 2, ]
# R-hat for convergence
r.hat.vals.1 <- gelman.diag(samples.1)
# Bayesian P-value
bpv.rep.y.1 <- unlist(samples[, which(param.names.1 == 'chi.2.rep.y')])
bpv.y.1 <- unlist(samples[, which(param.names.1 == 'chi.2.y')])
mean(bpv.rep.y.1 > bpv.y.1)
# Model NEON --------------------------
load("results/icom-NEON-results.R")
# Files are same as those described above for each separate model, 
# removing comments for clarity. 
samples.2 <- samples
param.names.2 <- attr(samples.2[[1]], 'dimnames')[[2]]
jac.neon.samples.2 <- jaccard.neon.samples
rich.neon.samples.2 <- rich.neon.samples
psi.samples.2 <- psi.avg.sp.year
trend.samples.2 <- beta.samples[, 2, ]
r.hat.vals.2 <- gelman.diag(samples.2)
# Bayesian P-value
bpv.rep.v.1.2 <- unlist(samples[, which(param.names.2 == 'chi.2.rep.v.1')])
bpv.v.1.2 <- unlist(samples[, which(param.names.2 == 'chi.2.v.1')])
mean(bpv.rep.v.1.2 > bpv.v.1.2)
# Model BBS ---------------------------
load("results/icom-BBS-results.R")
samples.3 <- samples
param.names.3 <- attr(samples.3[[1]], 'dimnames')[[2]]
jac.bbs.samples.3 <- jaccard.bbs.samples
rich.bbs.samples.3 <- rich.bbs.samples
psi.samples.3 <- psi.avg.sp.year
trend.samples.3 <- beta.samples[, 2, ]
r.hat.vals.3 <- gelman.diag(samples.3)
# Bayesian P-value
bpv.rep.v.2.3 <- unlist(samples[, which(param.names.3 == 'chi.2.rep.bbs')])
bpv.v.2.3 <- unlist(samples[, which(param.names.3 == 'chi.2.bbs')])
mean(bpv.rep.v.2.3 > bpv.v.2.3)
# Model HBEF + NEON -------------------
load("results/icom-HBEF-NEON-results.R")
samples.4 <- samples
param.names.4 <- attr(samples.4[[1]], 'dimnames')[[2]]
jac.neon.samples.4 <- jaccard.neon.samples
rich.neon.samples.4 <- rich.neon.samples
jac.hbef.samples.4 <- jaccard.hbef.samples
rich.hbef.samples.4 <- rich.hbef.samples
psi.samples.4 <- psi.avg.sp.year
trend.samples.4 <- beta.samples[, 2, ]
r.hat.vals.4 <- gelman.diag(samples.4)
# Bayesian P-value
bpv.rep.y.4 <- unlist(samples[, which(param.names.4 == 'chi.2.rep.y')])
bpv.y.4 <- unlist(samples[, which(param.names.4 == 'chi.2.y')])
mean(bpv.rep.y.4 > bpv.y.4)
bpv.rep.v.1.4 <- unlist(samples[, which(param.names.4 == 'chi.2.rep.v.1')])
bpv.v.1.4 <- unlist(samples[, which(param.names.4 == 'chi.2.v.1')])
mean(bpv.rep.v.1.4 > bpv.v.1.4)
# Model HBEF + BBS --------------------
load("results/icom-HBEF-BBS-results.R")
samples.5 <- samples
param.names.5 <- attr(samples.5[[1]], 'dimnames')[[2]]
jac.bbs.samples.5 <- jaccard.bbs.samples
rich.bbs.samples.5 <- rich.bbs.samples
jac.hbef.samples.5 <- jaccard.hbef.samples
rich.hbef.samples.5 <- rich.hbef.samples
psi.samples.5 <- psi.avg.sp.year
trend.samples.5 <- beta.samples[, 2, ]
r.hat.vals.5 <- gelman.diag(samples.5)
# Bayesian P-value
bpv.rep.y.5 <- unlist(samples[, which(param.names.5 == 'chi.2.rep.y')])
bpv.y.5 <- unlist(samples[, which(param.names.5 == 'chi.2.y')])
mean(bpv.rep.y.5 > bpv.y.5)
bpv.rep.v.2.5 <- unlist(samples[, which(param.names.5 == 'chi.2.rep.bbs')])
bpv.v.2.5 <- unlist(samples[, which(param.names.5 == 'chi.2.bbs')])
mean(bpv.rep.v.2.5 > bpv.v.2.5)
# Model NEON + BBS --------------------
load("results/icom-NEON-BBS-results.R")
samples.6 <- samples
param.names.6 <- attr(samples.6[[1]], 'dimnames')[[2]]
jac.neon.samples.6 <- jaccard.neon.samples
rich.neon.samples.6 <- rich.neon.samples
jac.bbs.samples.6 <- jaccard.bbs.samples
rich.bbs.samples.6 <- rich.bbs.samples
psi.samples.6 <- psi.avg.sp.year
trend.samples.6 <- beta.samples[, 2, ]
r.hat.vals.6 <- gelman.diag(samples.6)
# Bayesian P-value 
bpv.rep.v.1.6 <- unlist(samples[, which(param.names.6 == 'chi.2.rep.v.1')])
bpv.v.1.6 <- unlist(samples[, which(param.names.6 == 'chi.2.v.1')])
mean(bpv.rep.v.1.6 > bpv.v.1.6)
bpv.rep.v.2.6 <- unlist(samples[, which(param.names.6 == 'chi.2.rep.bbs')])
bpv.v.2.6 <- unlist(samples[, which(param.names.6 == 'chi.2.bbs')])
mean(bpv.rep.v.2.6 > bpv.v.2.6)
# Model HBEF + NEON + BBS -------------
load("results/icom-HBEF-NEON-BBS-results.R")
samples.7 <- samples
param.names.7 <- attr(samples.7[[1]], 'dimnames')[[2]]
jac.neon.samples.7 <- jaccard.neon.samples
rich.neon.samples.7 <- rich.neon.samples
jac.hbef.samples.7 <- jaccard.hbef.samples
rich.hbef.samples.7 <- rich.hbef.samples
jac.bbs.samples.7 <- jaccard.bbs.samples
rich.bbs.samples.7 <- rich.bbs.samples
psi.samples.7 <- psi.avg.sp.year
trend.samples.7 <- beta.samples[, 2, ]
r.hat.vals.7 <- gelman.diag(samples.7)
# Bayesian P-value
bpv.rep.y.7 <- unlist(samples[, which(param.names.7 == 'chi.2.rep.y')])
bpv.y.7 <- unlist(samples[, which(param.names.7 == 'chi.2.y')])
mean(bpv.rep.y.7 > bpv.y.7)
bpv.rep.v.1.7 <- unlist(samples[, which(param.names.7 == 'chi.2.rep.v.1')])
bpv.v.1.7 <- unlist(samples[, which(param.names.7 == 'chi.2.v.1')])
mean(bpv.rep.v.1.7 > bpv.v.1.7)
bpv.rep.v.2.7 <- unlist(samples[, which(param.names.7 == 'chi.2.rep.bbs')])
bpv.v.2.7 <- unlist(samples[, which(param.names.7 == 'chi.2.bbs')])
mean(bpv.rep.v.2.7 > bpv.v.2.7)

sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW',
        'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
sp.names <- c('American Redstart', 'Black-and-white Warbler', 'Blue-headed Vireo', 
	      'Blackburnian Warbler', 'Blackpoll Warbler', 'Black-throated Blue Warbler', 
	      'Black-throated Green Warbler', 'Canada Warbler', 'Magnolia Warbler', 
	      'Nashville Warbler', 'Ovenbird', 'Red-eyed Vireo')
n.sp <- length(sp)
sp.hbef <- sp.names[unique(y.df$Species)]
n.sp.hbef <- length(sp.hbef)
sp.neon <- sp.names[unique(v.1.df$Species)]
n.sp.neon <- length(sp.neon)
sp.bbs <- sp.names[unique(v.2.df$Species)]
n.sp.bbs <- length(sp.bbs)

# Create trend comparison plot --------------------------------------------
years <- 2010:2018
n.years <- length(years)
years.neon <- 2015:2018
n.years.neon <- length(years.neon)
# Model 1
psi.med.1 <- c(apply(psi.samples.1, c(2, 3), median))
psi.low.1 <- c(apply(psi.samples.1, c(2, 3), quantile, 0.025))
psi.high.1 <- c(apply(psi.samples.1, c(2, 3), quantile, 0.975))
trend.df.1 <- data.frame(med = psi.med.1, 
			 low = psi.low.1, 
			 high = psi.high.1, 
			 year = rep(years, each = n.sp),
			 sp = rep(sp.hbef, n.years),  
			 model = 'HBEF')
# Model 2
psi.med.2 <- c(apply(psi.samples.2, c(2, 3), median))
psi.low.2 <- c(apply(psi.samples.2, c(2, 3), quantile, 0.025))
psi.high.2 <- c(apply(psi.samples.2, c(2, 3), quantile, 0.975))
trend.df.2 <- data.frame(med = psi.med.2, 
			 low = psi.low.2, 
			 high = psi.high.2, 
			 year = rep(years.neon, each = n.sp.neon),
			 sp = rep(sp.neon, n.years.neon),  
			 model = 'NEON')
# Model 3
psi.med.3 <- c(apply(psi.samples.3, c(2, 3), median))
psi.low.3 <- c(apply(psi.samples.3, c(2, 3), quantile, 0.025))
psi.high.3 <- c(apply(psi.samples.3, c(2, 3), quantile, 0.975))
trend.df.3 <- data.frame(med = psi.med.3, 
			 low = psi.low.3, 
			 high = psi.high.3, 
			 year = rep(years, each = n.sp),
			 sp = rep(sp.bbs, n.years),  
			 model = 'BBS')
# Model 4
psi.med.4 <- c(apply(psi.samples.4, c(2, 3), median))
psi.low.4 <- c(apply(psi.samples.4, c(2, 3), quantile, 0.025))
psi.high.4 <- c(apply(psi.samples.4, c(2, 3), quantile, 0.975))
trend.df.4 <- data.frame(med = psi.med.4, 
			 low = psi.low.4, 
			 high = psi.high.4, 
			 year = rep(years, each = n.sp),
			 sp = rep(sp.hbef, n.years),  
			 model = 'HBEF + NEON')
# Model 5
psi.med.5 <- c(apply(psi.samples.5, c(2, 3), median))
psi.low.5 <- c(apply(psi.samples.5, c(2, 3), quantile, 0.025))
psi.high.5 <- c(apply(psi.samples.5, c(2, 3), quantile, 0.975))
trend.df.5 <- data.frame(med = psi.med.5, 
			 low = psi.low.5, 
			 high = psi.high.5, 
			 year = rep(years, each = n.sp),
			 sp = rep(sp.hbef, n.years),  
			 model = 'HBEF + BBS')
# Model 6
psi.med.6 <- c(apply(psi.samples.6, c(2, 3), median))
psi.low.6 <- c(apply(psi.samples.6, c(2, 3), quantile, 0.025))
psi.high.6 <- c(apply(psi.samples.6, c(2, 3), quantile, 0.975))
trend.df.6 <- data.frame(med = psi.med.6, 
			 low = psi.low.6, 
			 high = psi.high.6, 
			 year = rep(years, each = n.sp),
			 sp = rep(sp.bbs, n.years),  
			 model = 'NEON + BBS')
# Model 7
psi.med.7 <- c(apply(psi.samples.7, c(2, 3), median))
psi.low.7 <- c(apply(psi.samples.7, c(2, 3), quantile, 0.025))
psi.high.7 <- c(apply(psi.samples.7, c(2, 3), quantile, 0.975))
trend.df.7 <- data.frame(med = psi.med.7, 
			 low = psi.low.7, 
			 high = psi.high.7, 
			 year = rep(years, each = n.sp),
			 sp = rep(sp.hbef, n.years),  
			 model = 'HBEF + NEON + BBS')
# Combine all together (except the full model, which will be added 
#                       later for plotting)
trend.df <- rbind(trend.df.1, trend.df.2, trend.df.3, 
		  trend.df.4, trend.df.5, trend.df.6)


# Trend figure ------------------------------------------------------------
# Reorder factors for proper display
trend.df <- trend.df %>%
  mutate(model = factor(model, levels = c('HBEF', 'NEON', 'BBS', 
					  'HBEF + NEON', 'HBEF + BBS', 
					  'NEON + BBS', 'HBEF + NEON + BBS'), 
			ordered = TRUE))
# Colors for plot. 
my.colors <- c('HBEF' = '#CC79A7', 'NEON' = '#0072B2', 
	       'BBS' = '#F0E442', 'HBEF + NEON' = '#009E73', 
	       'HBEF + BBS' = '#56B4E9', 'NEON + BBS' = '#E69F00', 
	       'HBEF + NEON + BBS' = 'black')
# Full model to make this one at the forefront of the figure. 
trend.df.7 <- trend.df.7 %>%
  mutate(model = factor(model, levels = c('HBEF', 'NEON', 'BBS', 
					  'HBEF + NEON', 'HBEF + BBS', 
					  'NEON + BBS', 'HBEF + NEON + BBS'), 
			ordered = TRUE))
# Figure 4 in manuscript. 
trend.df %>% 
  ggplot(aes(x = year, y = med, col = model)) +
  geom_point(size = 3.2) + 
  geom_line(size = 0.7) + 
  geom_point(data = trend.df.7, aes(x = year, y = med, col = model), 
	     size = 3.2) + 
  geom_line(data = trend.df.7, aes(x = year, y = med, col = model), 
	    size = 0.7) + 
  geom_ribbon(data = trend.df.7, aes(x = year, ymin = low, ymax = high), 
	      alpha = 0.4, fill = 'grey', col = NA) + 
  facet_wrap(vars(sp), ncol = 3, nrow = 4) + 
  theme_bw(base_size = 18) + 
  labs(x = 'Year', y = 'Average Occurrence Probability', 
       col = 'Data Sets Used') + 
  scale_color_manual(values = my.colors, limits = c('HBEF', 'NEON', 'BBS', 
						    'HBEF + NEON', 'HBEF + BBS', 
						    'NEON + BBS', 'HBEF + NEON + BBS')) + 
  scale_y_continuous(limits = c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
	legend.position = 'bottom') 
ggsave(device = 'pdf', filename = 'figures/Fig4.pdf', 
       height = 12, width = 13, units = 'in')

# Community parameters for Summary Table ----------------------------------
# These results are used in Table S4. 
# Model 1 -----------------------------
curr.indx <- which(substr(param.names.1, 1, 13) == 'int.beta.mean')
int.samples.1 <- samples.1[, curr.indx]
int.samples.1 <- mcmc(do.call('rbind', int.samples.1))
beta.0.samples.1 <- summary(logit(int.samples.1))$quantiles
beta.0.samples.1[, 5] - beta.0.samples.1[, 1]
samples.1.quants <- summary(samples.1)$quantiles
samples.1.quants[, 5] - samples.1.quants[, 1]
# Model 2 -----------------------------
curr.indx <- which(substr(param.names.2, 1, 13) == 'int.beta.mean')
int.samples.2 <- samples.2[, curr.indx]
int.samples.2 <- mcmc(do.call('rbind', int.samples.2))
beta.0.samples.2 <- summary(logit(int.samples.2))$quantiles
beta.0.samples.2[, 5] - beta.0.samples.2[, 1]
samples.2.quants <- summary(samples.2)$quantiles
samples.2.quants[, 5] - samples.2.quants[, 1]
# Model 3 -----------------------------
curr.indx <- which(substr(param.names.3, 1, 13) == 'int.beta.mean')
int.samples.3 <- samples.3[, curr.indx]
int.samples.3 <- mcmc(do.call('rbind', int.samples.3))
beta.0.samples.3 <- summary(logit(int.samples.3))$quantiles
beta.0.samples.3[, 5] - beta.0.samples.3[, 1]
samples.3.quants <- summary(samples.3)$quantiles
samples.3.quants[, 5] - samples.3.quants[, 1]
# Model 4 -----------------------------
curr.indx <- which(substr(param.names.4, 1, 13) == 'int.beta.mean')
int.samples.4 <- samples.4[, curr.indx]
int.samples.4 <- mcmc(do.call('rbind', int.samples.4))
beta.0.samples.4 <- summary(logit(int.samples.4))$quantiles
beta.0.samples.4[, 5] - beta.0.samples.4[, 1]
samples.4.quants <- summary(samples.4)$quantiles
samples.4.quants[, 5] - samples.4.quants[, 1]
# Model 5 -----------------------------
curr.indx <- which(substr(param.names.5, 1, 13) == 'int.beta.mean')
int.samples.5 <- samples.5[, curr.indx]
int.samples.5 <- mcmc(do.call('rbind', int.samples.5))
beta.0.samples.5 <- summary(logit(int.samples.5))$quantiles
beta.0.samples.5[, 5] - beta.0.samples.5[, 1]
samples.5.quants <- summary(samples.5)$quantiles
samples.5.quants[, 5] - samples.5.quants[, 1]
# Model 6 -----------------------------
curr.indx <- which(substr(param.names.6, 1, 13) == 'int.beta.mean')
int.samples.6 <- samples.6[, curr.indx]
int.samples.6 <- mcmc(do.call('rbind', int.samples.6))
beta.0.samples.6 <- summary(logit(int.samples.6))$quantiles
beta.0.samples.6[, 5] - beta.0.samples.6[, 1]
samples.6.quants <- summary(samples.6)$quantiles
samples.6.quants[, 5] - samples.6.quants[, 1]
# Model 7 -----------------------------
curr.indx <- which(substr(param.names.7, 1, 13) == 'int.beta.mean')
int.samples.7 <- samples.7[, curr.indx]
int.samples.7 <- mcmc(do.call('rbind', int.samples.7))
beta.0.samples.7 <- summary(logit(int.samples.7))$quantiles
beta.0.samples.7[, 5] - beta.0.samples.7[, 1]
samples.7.quants <- summary(samples.7)$quantiles
samples.7.quants[, 5] - samples.7.quants[, 1]

# Summarize species-specific parameter precision --------------------------
# These results are used in Table S5 of the manuscript. 
# Intercept
# Model 1 -----------------------------
curr.indx <- which(substr(param.names.1, 1, 9) == 'int.beta[')
int.samples.1 <- samples.1[, curr.indx]
int.samples.1 <- mcmc(do.call('rbind', int.samples.1))
int.samples.1 <- summary(int.samples.1)$quantiles
int.mean.vals.1 <- apply(int.samples.1, 2, mean)
# Model 2 -----------------------------
curr.indx <- which(substr(param.names.2, 1, 9) == 'int.beta[')
int.samples.2 <- samples.2[, curr.indx]
int.samples.2 <- mcmc(do.call('rbind', int.samples.2))
int.samples.2 <- summary(int.samples.2)$quantiles
int.mean.vals.2 <- apply(int.samples.2, 2, mean)
# Model 3 -----------------------------
curr.indx <- which(substr(param.names.3, 1, 9) == 'int.beta[')
int.samples.3 <- samples.3[, curr.indx]
int.samples.3 <- mcmc(do.call('rbind', int.samples.3))
int.samples.3 <- summary(int.samples.3)$quantiles
int.mean.vals.3 <- apply(int.samples.3, 2, mean)
# Model 4 -----------------------------
curr.indx <- which(substr(param.names.4, 1, 9) == 'int.beta[')
int.samples.4 <- samples.4[, curr.indx]
int.samples.4 <- mcmc(do.call('rbind', int.samples.4))
int.samples.4 <- summary(int.samples.4)$quantiles
int.mean.vals.4 <- apply(int.samples.4, 2, mean)
# Model 5 -----------------------------
curr.indx <- which(substr(param.names.5, 1, 9) == 'int.beta[')
int.samples.5 <- samples.5[, curr.indx]
int.samples.5 <- mcmc(do.call('rbind', int.samples.5))
int.samples.5 <- summary(int.samples.5)$quantiles
int.mean.vals.5 <- apply(int.samples.5, 2, mean)
# Model 6 ------------------------------
curr.indx <- which(substr(param.names.6, 1, 9) == 'int.beta[')
int.samples.6 <- samples.6[, curr.indx]
int.samples.6 <- mcmc(do.call('rbind', int.samples.6))
int.samples.6 <- summary(int.samples.6)$quantiles
int.mean.vals.6 <- apply(int.samples.6, 2, mean)
# Model 7 -----------------------------
curr.indx <- which(substr(param.names.7, 1, 9) == 'int.beta[')
int.samples.7 <- samples.7[, curr.indx]
int.samples.7 <- mcmc(do.call('rbind', int.samples.7))
int.samples.7 <- summary(int.samples.7)$quantiles
int.mean.vals.7 <- apply(int.samples.7, 2, mean)

int.mean.vals.1[5] - int.mean.vals.1[1]
int.mean.vals.2[5] - int.mean.vals.2[1]
int.mean.vals.3[5] - int.mean.vals.3[1]
int.mean.vals.4[5] - int.mean.vals.4[1]
int.mean.vals.5[5] - int.mean.vals.5[1]
int.mean.vals.6[5] - int.mean.vals.6[1]
int.mean.vals.7[5] - int.mean.vals.7[1]
# Linear elevation 
# Model 1 -----------------------------
curr.indx <- which(substr(param.names.1, 1, 7) == 'beta.1[')
beta.1.samples.1 <- samples.1[, curr.indx]
beta.1.samples.1 <- mcmc(do.call('rbind', beta.1.samples.1))
beta.1.samples.1 <- summary(beta.1.samples.1)$quantiles
beta.1.mean.vals.1 <- apply(beta.1.samples.1, 2, mean)
# Model 2 -----------------------------
curr.indx <- which(substr(param.names.2, 1, 7) == 'beta.1[')
beta.1.samples.2 <- samples.2[, curr.indx]
beta.1.samples.2 <- mcmc(do.call('rbind', beta.1.samples.2))
beta.1.samples.2 <- summary(beta.1.samples.2)$quantiles
beta.1.mean.vals.2 <- apply(beta.1.samples.2, 2, mean)
# Model 3 -----------------------------
curr.indx <- which(substr(param.names.3, 1, 7) == 'beta.1[')
beta.1.samples.3 <- samples.3[, curr.indx]
beta.1.samples.3 <- mcmc(do.call('rbind', beta.1.samples.3))
beta.1.samples.3 <- summary(beta.1.samples.3)$quantiles
beta.1.mean.vals.3 <- apply(beta.1.samples.3, 2, mean)
# Model 4 -----------------------------
curr.indx <- which(substr(param.names.4, 1, 7) == 'beta.1[')
beta.1.samples.4 <- samples.4[, curr.indx]
beta.1.samples.4 <- mcmc(do.call('rbind', beta.1.samples.4))
beta.1.samples.4 <- summary(beta.1.samples.4)$quantiles
beta.1.mean.vals.4 <- apply(beta.1.samples.4, 2, mean)
# Model 5 -----------------------------
curr.indx <- which(substr(param.names.5, 1, 7) == 'beta.1[')
beta.1.samples.5 <- samples.5[, curr.indx]
beta.1.samples.5 <- mcmc(do.call('rbind', beta.1.samples.5))
beta.1.samples.5 <- summary(beta.1.samples.5)$quantiles
beta.1.mean.vals.5 <- apply(beta.1.samples.5, 2, mean)
# Model 6 -----------------------------
curr.indx <- which(substr(param.names.6, 1, 7) == 'beta.1[')
beta.1.samples.6 <- samples.6[, curr.indx]
beta.1.samples.6 <- mcmc(do.call('rbind', beta.1.samples.6))
beta.1.samples.6 <- summary(beta.1.samples.6)$quantiles
beta.1.mean.vals.6 <- apply(beta.1.samples.6, 2, mean)
# Model 7 -----------------------------
curr.indx <- which(substr(param.names.7, 1, 7) == 'beta.1[')
beta.1.samples.7 <- samples.7[, curr.indx]
beta.1.samples.7 <- mcmc(do.call('rbind', beta.1.samples.7))
beta.1.samples.7 <- summary(beta.1.samples.7)$quantiles
beta.1.mean.vals.7 <- apply(beta.1.samples.7, 2, mean)

beta.1.mean.vals.1[5] - beta.1.mean.vals.1[1]
beta.1.mean.vals.2[5] - beta.1.mean.vals.2[1]
beta.1.mean.vals.3[5] - beta.1.mean.vals.3[1]
beta.1.mean.vals.4[5] - beta.1.mean.vals.4[1]
beta.1.mean.vals.5[5] - beta.1.mean.vals.5[1]
beta.1.mean.vals.6[5] - beta.1.mean.vals.6[1]
beta.1.mean.vals.7[5] - beta.1.mean.vals.7[1]

# Quadratic elevation 
# Model 1 -----------------------------
curr.indx <- which(substr(param.names.1, 1, 7) == 'beta.2[')
beta.2.samples.1 <- samples.1[, curr.indx]
beta.2.samples.1 <- mcmc(do.call('rbind', beta.2.samples.1))
beta.2.samples.1 <- summary(beta.2.samples.1)$quantiles
beta.2.mean.vals.1 <- apply(beta.2.samples.1, 2, mean)
# Model 2 -----------------------------
curr.indx <- which(substr(param.names.2, 1, 7) == 'beta.2[')
beta.2.samples.2 <- samples.2[, curr.indx]
beta.2.samples.2 <- mcmc(do.call('rbind', beta.2.samples.2))
beta.2.samples.2 <- summary(beta.2.samples.2)$quantiles
beta.2.mean.vals.2 <- apply(beta.2.samples.2, 2, mean)
# Model 3 -----------------------------
curr.indx <- which(substr(param.names.3, 1, 7) == 'beta.2[')
beta.2.samples.3 <- samples.3[, curr.indx]
beta.2.samples.3 <- mcmc(do.call('rbind', beta.2.samples.3))
beta.2.samples.3 <- summary(beta.2.samples.3)$quantiles
beta.2.mean.vals.3 <- apply(beta.2.samples.3, 2, mean)
# Model 4 -----------------------------
curr.indx <- which(substr(param.names.4, 1, 7) == 'beta.2[')
beta.2.samples.4 <- samples.4[, curr.indx]
beta.2.samples.4 <- mcmc(do.call('rbind', beta.2.samples.4))
beta.2.samples.4 <- summary(beta.2.samples.4)$quantiles
beta.2.mean.vals.4 <- apply(beta.2.samples.4, 2, mean)
# Model 5 -----------------------------
curr.indx <- which(substr(param.names.5, 1, 7) == 'beta.2[')
beta.2.samples.5 <- samples.5[, curr.indx]
beta.2.samples.5 <- mcmc(do.call('rbind', beta.2.samples.5))
beta.2.samples.5 <- summary(beta.2.samples.5)$quantiles
beta.2.mean.vals.5 <- apply(beta.2.samples.5, 2, mean)
# Model 6 -----------------------------
curr.indx <- which(substr(param.names.6, 1, 7) == 'beta.2[')
beta.2.samples.6 <- samples.6[, curr.indx]
beta.2.samples.6 <- mcmc(do.call('rbind', beta.2.samples.6))
beta.2.samples.6 <- summary(beta.2.samples.6)$quantiles
beta.2.mean.vals.6 <- apply(beta.2.samples.6, 2, mean)
# Model 7 -----------------------------
curr.indx <- which(substr(param.names.7, 1, 7) == 'beta.2[')
beta.2.samples.7 <- samples.7[, curr.indx]
beta.2.samples.7 <- mcmc(do.call('rbind', beta.2.samples.7))
beta.2.samples.7 <- summary(beta.2.samples.7)$quantiles
beta.2.mean.vals.7 <- apply(beta.2.samples.7, 2, mean)

beta.2.mean.vals.1[5] - beta.2.mean.vals.1[1]
beta.2.mean.vals.2[5] - beta.2.mean.vals.2[1]
beta.2.mean.vals.3[5] - beta.2.mean.vals.3[1]
beta.2.mean.vals.4[5] - beta.2.mean.vals.4[1]
beta.2.mean.vals.5[5] - beta.2.mean.vals.5[1]
beta.2.mean.vals.6[5] - beta.2.mean.vals.6[1]
beta.2.mean.vals.7[5] - beta.2.mean.vals.7[1]
# Percent Forest
# Model 1 -----------------------------
curr.indx <- which(substr(param.names.1, 1, 7) == 'beta.3[')
beta.3.samples.1 <- samples.1[, curr.indx]
beta.3.samples.1 <- mcmc(do.call('rbind', beta.3.samples.1))
beta.3.samples.1 <- summary(beta.3.samples.1)$quantiles
beta.3.mean.vals.1 <- apply(beta.3.samples.1, 2, mean)
# Model 2 -----------------------------
curr.indx <- which(substr(param.names.2, 1, 7) == 'beta.3[')
beta.3.samples.2 <- samples.2[, curr.indx]
beta.3.samples.2 <- mcmc(do.call('rbind', beta.3.samples.2))
beta.3.samples.2 <- summary(beta.3.samples.2)$quantiles
beta.3.mean.vals.2 <- apply(beta.3.samples.2, 2, mean)
# Model 3 -----------------------------
curr.indx <- which(substr(param.names.3, 1, 7) == 'beta.3[')
beta.3.samples.3 <- samples.3[, curr.indx]
beta.3.samples.3 <- mcmc(do.call('rbind', beta.3.samples.3))
beta.3.samples.3 <- summary(beta.3.samples.3)$quantiles
beta.3.mean.vals.3 <- apply(beta.3.samples.3, 2, mean)
# Model 4 -----------------------------
curr.indx <- which(substr(param.names.4, 1, 7) == 'beta.3[')
beta.3.samples.4 <- samples.4[, curr.indx]
beta.3.samples.4 <- mcmc(do.call('rbind', beta.3.samples.4))
beta.3.samples.4 <- summary(beta.3.samples.4)$quantiles
beta.3.mean.vals.4 <- apply(beta.3.samples.4, 2, mean)
# Model 5 -----------------------------
curr.indx <- which(substr(param.names.5, 1, 7) == 'beta.3[')
beta.3.samples.5 <- samples.5[, curr.indx]
beta.3.samples.5 <- mcmc(do.call('rbind', beta.3.samples.5))
beta.3.samples.5 <- summary(beta.3.samples.5)$quantiles
beta.3.mean.vals.5 <- apply(beta.3.samples.5, 2, mean)
# Model 6 -----------------------------
curr.indx <- which(substr(param.names.6, 1, 7) == 'beta.3[')
beta.3.samples.6 <- samples.6[, curr.indx]
beta.3.samples.6 <- mcmc(do.call('rbind', beta.3.samples.6))
beta.3.samples.6 <- summary(beta.3.samples.6)$quantiles
beta.3.mean.vals.6 <- apply(beta.3.samples.6, 2, mean)
# Model 7 -----------------------------
curr.indx <- which(substr(param.names.7, 1, 7) == 'beta.3[')
beta.3.samples.7 <- samples.7[, curr.indx]
beta.3.samples.7 <- mcmc(do.call('rbind', beta.3.samples.7))
beta.3.samples.7 <- summary(beta.3.samples.7)$quantiles
beta.3.mean.vals.7 <- apply(beta.3.samples.7, 2, mean)

beta.3.mean.vals.1[5] - beta.3.mean.vals.1[1]
beta.3.mean.vals.2[5] - beta.3.mean.vals.2[1]
beta.3.mean.vals.3[5] - beta.3.mean.vals.3[1]
beta.3.mean.vals.4[5] - beta.3.mean.vals.4[1]
beta.3.mean.vals.5[5] - beta.3.mean.vals.5[1]
beta.3.mean.vals.6[5] - beta.3.mean.vals.6[1]
beta.3.mean.vals.7[5] - beta.3.mean.vals.7[1]
# Autologistic
# Model 1 -----------------------------
curr.indx <- which(substr(param.names.1, 1, 4) == 'phi[')
phi.samples.1 <- samples.1[, curr.indx]
phi.samples.1 <- mcmc(do.call('rbind', phi.samples.1))
phi.samples.1 <- summary(phi.samples.1)$quantiles
phi.mean.vals.1 <- apply(phi.samples.1, 2, mean)
# Model 2 -----------------------------
curr.indx <- which(substr(param.names.2, 1, 4) == 'phi[')
phi.samples.2 <- samples.2[, curr.indx]
phi.samples.2 <- mcmc(do.call('rbind', phi.samples.2))
phi.samples.2 <- summary(phi.samples.2)$quantiles
phi.mean.vals.2 <- apply(phi.samples.2, 2, mean)
# Model 3 -----------------------------
curr.indx <- which(substr(param.names.3, 1, 4) == 'phi[')
phi.samples.3 <- samples.3[, curr.indx]
phi.samples.3 <- mcmc(do.call('rbind', phi.samples.3))
phi.samples.3 <- summary(phi.samples.3)$quantiles
phi.mean.vals.3 <- apply(phi.samples.3, 2, mean)
# Model 4 -----------------------------
curr.indx <- which(substr(param.names.4, 1, 4) == 'phi[')
phi.samples.4 <- samples.4[, curr.indx]
phi.samples.4 <- mcmc(do.call('rbind', phi.samples.4))
phi.samples.4 <- summary(phi.samples.4)$quantiles
phi.mean.vals.4 <- apply(phi.samples.4, 2, mean)
# Model 5 -----------------------------
curr.indx <- which(substr(param.names.5, 1, 4) == 'phi[')
phi.samples.5 <- samples.5[, curr.indx]
phi.samples.5 <- mcmc(do.call('rbind', phi.samples.5))
phi.samples.5 <- summary(phi.samples.5)$quantiles
phi.mean.vals.5 <- apply(phi.samples.5, 2, mean)
# Model 6 -----------------------------
curr.indx <- which(substr(param.names.6, 1, 4) == 'phi[')
phi.samples.6 <- samples.6[, curr.indx]
phi.samples.6 <- mcmc(do.call('rbind', phi.samples.6))
phi.samples.6 <- summary(phi.samples.6)$quantiles
phi.mean.vals.6 <- apply(phi.samples.6, 2, mean)
# Model 7 -----------------------------
curr.indx <- which(substr(param.names.7, 1, 4) == 'phi[')
phi.samples.7 <- samples.7[, curr.indx]
phi.samples.7 <- mcmc(do.call('rbind', phi.samples.7))
phi.samples.7 <- summary(phi.samples.7)$quantiles
phi.mean.vals.7 <- apply(phi.samples.7, 2, mean)

phi.mean.vals.1[5] - phi.mean.vals.1[1]
phi.mean.vals.2[5] - phi.mean.vals.2[1]
phi.mean.vals.3[5] - phi.mean.vals.3[1]
phi.mean.vals.4[5] - phi.mean.vals.4[1]
phi.mean.vals.5[5] - phi.mean.vals.5[1]
phi.mean.vals.6[5] - phi.mean.vals.6[1]
phi.mean.vals.7[5] - phi.mean.vals.7[1]

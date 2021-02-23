rm(list = ls())
library(coda)
library(jagsUI)
library(dplyr)
library(ggplot2)
library(ggthemes)

# Read in results ---------------------------------------------------------
load("results/sim-icm-results-1-40-simulations-2021-01-29.R")

# Number of species
K <- length(beta.gamma.0)
# Number of simulations
n.sim <- 40
# Number of scenarios
n.scenarios <- length(R.hb.vals)

# Initial occupancy intercept ---------------------------------------------

beta.psi.0.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.psi.0})
beta.psi.0.med <- sapply(out.model, FUN = function(a) {a$q50$beta.psi.0})
beta.psi.0.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.psi.0})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')

beta.psi.0.df <- data.frame(low = c(beta.psi.0.low), 
			med = c(beta.psi.0.med), 
			high = c(beta.psi.0.high), 
			sp = rep(species.names, n.sim * n.scenarios), 
			scenario = factor(rep(rep(R.hb.vals / max(R.hb.vals), each = K), n.sim)), 
			rep.scen = rep(1:10, each = K * n.scenarios), 
			true = rep(beta.psi.0, n.sim * n.scenarios)
			)

# Compuate the average mean and 95% CI for parameters across all scenarios
plot.dat.beta.psi.0 <- beta.psi.0.df %>%
  group_by(scenario, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
 mutate(sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				   'Species 4', 'Species 5', 'Species 6', 
				   'Species 7', 'Species 8', 'Species 9', 
				   'Species 10'))) 

# Some plots
ggplot(data = plot.dat.beta.psi.0, aes(x = scenario, y = med, col = scenario)) + 
  geom_point(size = 3.5) + 
  facet_wrap(vars(sp), ncol = 3) + 
  theme_bw(base_size = 18) + 
  geom_hline(aes(yintercept = true), linetype = 2) + 
  scale_color_colorblind() + 
  geom_segment(aes(xend = scenario, y = low, yend = high), 
	       size = 1, lineend = 'round') + 
  labs(x = 'Proportion of Sites with Fine Res Data', y = 'Intercept', 
       col = 'Proportion of Sites with Fine Res Data') + 
  theme(legend.position = 'top')

# Initial Occupancy Covariate ---------------------------------------------

beta.psi.1.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.psi.1})
beta.psi.1.med <- sapply(out.model, FUN = function(a) {a$q50$beta.psi.1})
beta.psi.1.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.psi.1})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')

beta.psi.1.df <- data.frame(low = c(beta.psi.1.low), 
			med = c(beta.psi.1.med), 
			high = c(beta.psi.1.high), 
			sp = rep(species.names, n.sim * n.scenarios), 
			scenario = factor(rep(rep(R.hb.vals / max(R.hb.vals), each = K), n.sim)), 
			rep.scen = rep(1:10, each = K * n.scenarios), 
			true = rep(beta.psi.1, n.sim * n.scenarios)
			)


plot.dat.beta.psi.1 <- beta.psi.1.df %>%
  group_by(scenario, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
 mutate(sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				   'Species 4', 'Species 5', 'Species 6', 
				   'Species 7', 'Species 8', 'Species 9', 
				   'Species 10'))) 

ggplot(data = plot.dat.beta.psi.1, aes(x = scenario, y = med, col = scenario)) + 
  geom_point(size = 3.5) + 
  facet_wrap(vars(sp), ncol = 3) + 
  theme_bw(base_size = 18) + 
  geom_hline(aes(yintercept = true), linetype = 2) + 
  scale_color_colorblind() + 
  geom_segment(aes(xend = scenario, y = low, yend = high), 
	       size = 1, lineend = 'round') + 
  labs(x = 'Proportion of Sites with Fine Res Data', y = 'Covariate', 
       col = 'Proportion of Sites with Fine Res Data') + 
  theme(legend.position = 'top')


# Persistence intercept ---------------------------------------------------

beta.phi.0.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.phi.0})
beta.phi.0.med <- sapply(out.model, FUN = function(a) {a$q50$beta.phi.0})
beta.phi.0.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.phi.0})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')

beta.phi.0.df <- data.frame(low = c(beta.phi.0.low), 
			med = c(beta.phi.0.med), 
			high = c(beta.phi.0.high), 
			sp = rep(species.names, n.sim * n.scenarios), 
			scenario = factor(rep(rep(R.hb.vals / max(R.hb.vals), each = K), n.sim)), 
			rep.scen = rep(1:10, each = K * n.scenarios), 
			true = rep(beta.phi.0, n.sim * n.scenarios)
			)


plot.dat.beta.phi.0 <- beta.phi.0.df %>%
  group_by(scenario, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
 mutate(sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				   'Species 4', 'Species 5', 'Species 6', 
				   'Species 7', 'Species 8', 'Species 9', 
				   'Species 10'))) 

ggplot(data = plot.dat.beta.phi.0, aes(x = scenario, y = med, col = scenario)) + 
  geom_point(size = 3.5) + 
  facet_wrap(vars(sp), ncol = 3) + 
  theme_bw(base_size = 18) + 
  geom_hline(aes(yintercept = true), linetype = 2) + 
  scale_color_colorblind() + 
  geom_segment(aes(xend = scenario, y = low, yend = high), 
	       size = 1, lineend = 'round') + 
  labs(x = 'Proportion of Sites with Fine Res Data', y = 'Intercept', 
       col = 'Proportion of Sites with Fine Res Data') + 
  theme(legend.position = 'top')

# Persistence Covariate ---------------------------------------------------

beta.phi.1.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.phi.1})
beta.phi.1.med <- sapply(out.model, FUN = function(a) {a$q50$beta.phi.1})
beta.phi.1.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.phi.1})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')

beta.phi.1.df <- data.frame(low = c(beta.phi.1.low), 
			med = c(beta.phi.1.med), 
			high = c(beta.phi.1.high), 
			sp = rep(species.names, n.sim * n.scenarios), 
			scenario = factor(rep(rep(R.hb.vals / max(R.hb.vals), each = K), n.sim)), 
			rep.scen = rep(1:10, each = K * n.scenarios), 
			true = rep(beta.phi.1, n.sim * n.scenarios)
			)


plot.dat.beta.phi.1 <- beta.phi.1.df %>%
  group_by(scenario, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
 mutate(sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				   'Species 4', 'Species 5', 'Species 6', 
				   'Species 7', 'Species 8', 'Species 9', 
				   'Species 10'))) 

ggplot(data = plot.dat.beta.phi.1, aes(x = scenario, y = med, col = scenario)) + 
  geom_point(size = 3.5) + 
  facet_wrap(vars(sp), ncol = 3) + 
  theme_bw(base_size = 18) + 
  geom_hline(aes(yintercept = true), linetype = 2) + 
  scale_color_colorblind() + 
  geom_segment(aes(xend = scenario, y = low, yend = high), 
	       size = 1, lineend = 'round') + 
  labs(x = 'Proportion of Sites with Fine Res Data', y = 'Covariate', 
       col = 'Proportion of Sites with Fine Res Data') + 
  theme(legend.position = 'top')


# Colonization intercept ---------------------------------------------------

beta.gamma.0.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.gamma.0})
beta.gamma.0.med <- sapply(out.model, FUN = function(a) {a$q50$beta.gamma.0})
beta.gamma.0.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.gamma.0})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')

beta.gamma.0.df <- data.frame(low = c(beta.gamma.0.low), 
			med = c(beta.gamma.0.med), 
			high = c(beta.gamma.0.high), 
			sp = rep(species.names, n.sim * n.scenarios), 
			scenario = factor(rep(rep(R.hb.vals / max(R.hb.vals), each = K), n.sim)), 
			rep.scen = rep(1:10, each = K * n.scenarios), 
			true = rep(beta.gamma.0, n.sim * n.scenarios)
			)


plot.dat.beta.gamma.0 <- beta.gamma.0.df %>%
  group_by(scenario, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
 mutate(sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				   'Species 4', 'Species 5', 'Species 6', 
				   'Species 7', 'Species 8', 'Species 9', 
				   'Species 10'))) 

ggplot(data = plot.dat.beta.gamma.0, aes(x = scenario, y = med, col = scenario)) + 
  geom_point(size = 3.5) + 
  facet_wrap(vars(sp), ncol = 3) + 
  theme_bw(base_size = 18) + 
  geom_hline(aes(yintercept = true), linetype = 2) + 
  scale_color_colorblind() + 
  geom_segment(aes(xend = scenario, y = low, yend = high), 
	       size = 1, lineend = 'round') + 
  labs(x = 'Proportion of Sites with Fine Res Data', y = 'Intercept', 
       col = 'Proportion of Sites with Fine Res Data') + 
  theme(legend.position = 'top')

# Colonization Covariate ---------------------------------------------------

beta.gamma.1.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.gamma.1})
beta.gamma.1.med <- sapply(out.model, FUN = function(a) {a$q50$beta.gamma.1})
beta.gamma.1.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.gamma.1})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')

beta.gamma.1.df <- data.frame(low = c(beta.gamma.1.low), 
			med = c(beta.gamma.1.med), 
			high = c(beta.gamma.1.high), 
			sp = rep(species.names, n.sim * n.scenarios), 
			scenario = factor(rep(rep(R.hb.vals / max(R.hb.vals), each = K), n.sim)), 
			rep.scen = rep(1:10, each = K * n.scenarios), 
			true = rep(beta.gamma.1, n.sim * n.scenarios)
			)


plot.dat.beta.gamma.1 <- beta.gamma.1.df %>%
  group_by(scenario, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
 mutate(sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				   'Species 4', 'Species 5', 'Species 6', 
				   'Species 7', 'Species 8', 'Species 9', 
				   'Species 10'))) 

gamma.1.plot <- ggplot(data = plot.dat.beta.gamma.1, aes(x = scenario, y = med, col = scenario)) + 
  geom_point(size = 3.5) + 
  facet_wrap(vars(sp), ncol = 3) + 
  theme_bw(base_size = 18) + 
  geom_hline(aes(yintercept = true), linetype = 2) + 
  scale_color_colorblind() + 
  geom_segment(aes(xend = scenario, y = low, yend = high), 
	       size = 1, lineend = 'round') + 
  labs(x = 'Proportion of Sites with Fine Res Data', y = 'Covariate', 
       col = 'Proportion of Sites with Fine Res Data') + 
  theme(legend.position = 'top')
ggsave('figures/gamma-1-sim-1.png', gamma.1.plot, width = 11, height = 10, 
       units = 'in')


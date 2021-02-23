rm(list = ls())
library(coda)
library(jagsUI)
library(dplyr)
library(ggplot2)
library(ggthemes)

# Read in results ---------------------------------------------------------
load("results/sim-icm-results-2-40-simulations-2021-01-30.R")

# Number of species
K <- length(beta.gamma.0)
# Number of simulations
n.sim <- 40
# Number of simulation scenarios
n.scenarios <- nrow(param.vals)

# Initial occupancy intercept ---------------------------------------------

beta.psi.0.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.psi.0})
beta.psi.0.med <- sapply(out.model, FUN = function(a) {a$q50$beta.psi.0})
beta.psi.0.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.psi.0})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')
n.years.vals <- c(5, 10)
J.hb.vals <- c(1, 3)
J.neon.vals <- c(1, 3)
n.sp <- length(species.names)
# c() goes down each row, then across columns. 
# Get data in data frame for plotting. 
beta.psi.0.df <- data.frame(low = c(beta.psi.0.low), 
			med = c(beta.psi.0.med), 
			high = c(beta.psi.0.high), 
			sp = rep(1:10, n.sim * n.scenarios), 
			n.years = rep(rep(n.years.vals, each = n.sp), 
				      times = n.sim * n.scenarios / length(n.years.vals)),
                        J.hb = rep(rep(J.hb.vals, each = n.sp * length(n.years.vals)), 
				   times = n.sim * n.scenarios / (length(n.years.vals) * 
								  length(J.hb.vals))), 
			J.neon = rep(rep(J.neon.vals, each = n.sp * length(n.years.vals) * 
					 length(J.hb.vals)), 
				     times = n.sim * n.scenarios / (length(n.years.vals) * 
								    length(J.hb.vals) * 
								    length(J.neon.vals))), 
			true = rep(beta.psi.0, n.sim * n.scenarios)
			)

# Get average medians and 95% CIs across all simulations for each scenario. 
plot.dat.beta.psi.0 <- beta.psi.0.df %>%
  group_by(n.years, J.hb, J.neon, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
  mutate(n.years = factor(n.years), 
	 J.hb = factor(J.hb), 
	 J.neon = factor(J.neon)) %>%
  ungroup()

plot.dat.beta.psi.0 %>%
  filter(n.years == 5) %>% 
ggplot(aes(x = sp, y = med)) + 
  geom_point(size = 3.5, col = 'black') + 
  geom_point(aes(x = sp + 0.25, y = true), col = 'gray', size = 3.5) + 
  geom_segment(aes(x = sp, xend = sp, y = low, yend = high), 
	       size = 1, lineend = 'round') +  
  facet_grid(vars(J.hb), vars(J.neon)) + 
  theme_bw() + 
  labs(x = 'Species', y = 'Estimate')

# Initial occupancy covariate ---------------------------------------------

beta.psi.1.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.psi.1})
beta.psi.1.med <- sapply(out.model, FUN = function(a) {a$q50$beta.psi.1})
beta.psi.1.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.psi.1})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4', 
		   'Species 5', 'Species 6', 'Species 7', 'Species 8', 
		   'Species 9', 'Species 10')
n.years.vals <- c(5, 10)
J.hb.vals <- c(1, 3)
J.neon.vals <- c(1, 3)
n.sp <- length(species.names)
# c goes down each row, then across columns. 
beta.psi.1.df <- data.frame(low = c(beta.psi.1.low), 
			med = c(beta.psi.1.med), 
			high = c(beta.psi.1.high), 
			sp = rep(1:10, n.sim * n.scenarios), 
			n.years = rep(rep(n.years.vals, each = n.sp), 
				      times = n.sim * n.scenarios / length(n.years.vals)),
                        J.hb = rep(rep(J.hb.vals, each = n.sp * length(n.years.vals)), 
				   times = n.sim * n.scenarios / (length(n.years.vals) * 
								  length(J.hb.vals))), 
			J.neon = rep(rep(J.neon.vals, each = n.sp * length(n.years.vals) * 
					 length(J.hb.vals)), 
				     times = n.sim * n.scenarios / (length(n.years.vals) * 
								    length(J.hb.vals) * 
								    length(J.neon.vals))), 
			true = rep(beta.psi.1, n.sim * n.scenarios)
			)

plot.dat.beta.psi.1 <- beta.psi.1.df %>%
  group_by(n.years, J.hb, J.neon, sp, true) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high)) %>%
  mutate(n.years = factor(n.years), 
	 J.hb = factor(J.hb), 
	 J.neon = factor(J.neon)) %>%
  ungroup()

plot.dat.beta.psi.1 %>%
  filter(n.years == 5) %>% 
ggplot(aes(x = sp, y = med)) + 
  geom_point(size = 3.5, col = 'black') + 
  geom_point(aes(x = sp + 0.25, y = true), col = 'gray', size = 3.5) + 
  geom_segment(aes(x = sp, xend = sp, y = low, yend = high), 
	       size = 1, lineend = 'round') +  
  facet_grid(vars(J.hb), vars(J.neon)) + 
  theme_bw() + 
  labs(x = 'Species', y = 'Estimate')

# Persistence intercept ---------------------------------------------------

beta.phi.0.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.phi.0})
beta.phi.0.med <- sapply(out.model, FUN = function(a) {a$q50$beta.phi.0})
beta.phi.0.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.phi.0})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4',
		   'Species 5', 'Species 6', 'Species 7', 'Species 8',
		   'Species 9', 'Species 10')
n.years.vals <- c(5, 10)
J.hb.vals <- c(1, 3)
J.neon.vals <- c(1, 3)
n.sp <- length(species.names)
# c goes down each row, then across columns.
beta.phi.0.df <- data.frame(low = c(beta.phi.0.low),
			med = c(beta.phi.0.med),
			high = c(beta.phi.0.high),
			sp = rep(1:10, n.sim * n.scenarios),
			n.years = rep(rep(n.years.vals, each = n.sp),
				      times = n.sim * n.scenarios / length(n.years.vals)),
                        J.hb = rep(rep(J.hb.vals, each = n.sp * length(n.years.vals)),
				   times = n.sim * n.scenarios / (length(n.years.vals) *
								  length(J.hb.vals))),
			J.neon = rep(rep(J.neon.vals, each = n.sp * length(n.years.vals) *
					 length(J.hb.vals)),
				     times = n.sim * n.scenarios / (length(n.years.vals) *
								    length(J.hb.vals) *
								    length(J.neon.vals))),
			true = rep(beta.phi.0, n.sim * n.scenarios)
			)

plot.dat.beta.phi.0 <- beta.phi.0.df %>%
  group_by(n.years, J.hb, J.neon, sp, true) %>%
  summarize(low = mean(low),
	    med = mean(med),
	    high = mean(high)) %>%
  mutate(n.years = factor(n.years),
	 J.hb = factor(J.hb),
	 J.neon = factor(J.neon)) %>%
  ungroup()

phi.int.plot <- plot.dat.beta.phi.0 %>%
  filter(n.years == 5) %>%
ggplot(aes(x = sp, y = med)) +
  geom_point(size = 3.5, col = 'black') +
  geom_point(aes(x = sp + 0.25, y = true), col = 'gray', size = 3.5) +
  geom_segment(aes(x = sp, xend = sp, y = low, yend = high),
	       size = 1, lineend = 'round') +
  facet_grid(vars(J.hb), vars(J.neon)) +
  theme_bw(base_size = 18) +
  labs(x = 'Species', y = 'Estimate')
#ggsave('figures/phi-int-sim-2.png', phi.int.plot)

# Persistence covariate ---------------------------------------------------

beta.phi.1.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.phi.1})
beta.phi.1.med <- sapply(out.model, FUN = function(a) {a$q50$beta.phi.1})
beta.phi.1.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.phi.1})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4',
		   'Species 5', 'Species 6', 'Species 7', 'Species 8',
		   'Species 9', 'Species 10')
n.years.vals <- c(5, 10)
J.hb.vals <- c(1, 3)
J.neon.vals <- c(1, 3)
n.sp <- length(species.names)
# c goes down each row, then across columns.
beta.phi.1.df <- data.frame(low = c(beta.phi.1.low),
			med = c(beta.phi.1.med),
			high = c(beta.phi.1.high),
			sp = rep(1:10, n.sim * n.scenarios),
			n.years = rep(rep(n.years.vals, each = n.sp),
				      times = n.sim * n.scenarios / length(n.years.vals)),
                        J.hb = rep(rep(J.hb.vals, each = n.sp * length(n.years.vals)),
				   times = n.sim * n.scenarios / (length(n.years.vals) *
								  length(J.hb.vals))),
			J.neon = rep(rep(J.neon.vals, each = n.sp * length(n.years.vals) *
					 length(J.hb.vals)),
				     times = n.sim * n.scenarios / (length(n.years.vals) *
								    length(J.hb.vals) *
								    length(J.neon.vals))),
			true = rep(beta.phi.1, n.sim * n.scenarios)
			)

plot.dat.beta.phi.1 <- beta.phi.1.df %>%
  group_by(n.years, J.hb, J.neon, sp, true) %>%
  summarize(low = mean(low),
	    med = mean(med),
	    high = mean(high)) %>%
  mutate(n.years = factor(n.years),
	 J.hb = factor(J.hb),
	 J.neon = factor(J.neon)) %>%
  ungroup()

plot.dat.beta.phi.1 %>%
  filter(n.years == 5) %>%
ggplot(aes(x = sp, y = med)) +
  geom_point(size = 3.5, col = 'black') +
  geom_point(aes(x = sp + 0.25, y = true), col = 'gray', size = 3.5) +
  geom_segment(aes(x = sp, xend = sp, y = low, yend = high),
	       size = 1, lineend = 'round') +
  facet_grid(vars(J.hb), vars(J.neon)) +
  theme_bw() +
  labs(x = 'Species', y = 'Estimate')


# Colonization intercept --------------------------------------------------

beta.gamma.0.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.gamma.0})
beta.gamma.0.med <- sapply(out.model, FUN = function(a) {a$q50$beta.gamma.0})
beta.gamma.0.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.gamma.0})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4',
		   'Species 5', 'Species 6', 'Species 7', 'Species 8',
		   'Species 9', 'Species 10')
n.years.vals <- c(5, 10)
J.hb.vals <- c(1, 3)
J.neon.vals <- c(1, 3)
n.sp <- length(species.names)
# c goes down each row, then across columns.
beta.gamma.0.df <- data.frame(low = c(beta.gamma.0.low),
			med = c(beta.gamma.0.med),
			high = c(beta.gamma.0.high),
			sp = rep(1:10, n.sim * n.scenarios),
			n.years = rep(rep(n.years.vals, each = n.sp),
				      times = n.sim * n.scenarios / length(n.years.vals)),
                        J.hb = rep(rep(J.hb.vals, each = n.sp * length(n.years.vals)),
				   times = n.sim * n.scenarios / (length(n.years.vals) *
								  length(J.hb.vals))),
			J.neon = rep(rep(J.neon.vals, each = n.sp * length(n.years.vals) *
					 length(J.hb.vals)),
				     times = n.sim * n.scenarios / (length(n.years.vals) *
								    length(J.hb.vals) *
								    length(J.neon.vals))),
			true = rep(beta.gamma.0, n.sim * n.scenarios)
			)

plot.dat.beta.gamma.0 <- beta.gamma.0.df %>%
  group_by(n.years, J.hb, J.neon, sp, true) %>%
  summarize(low = mean(low),
	    med = mean(med),
	    high = mean(high)) %>%
  mutate(n.years = factor(n.years),
	 J.hb = factor(J.hb),
	 J.neon = factor(J.neon)) %>%
  ungroup()

plot.dat.beta.gamma.0 %>%
  filter(n.years == 5) %>%
ggplot(aes(x = sp, y = med)) +
  geom_point(size = 3.5, col = 'black') +
  geom_point(aes(x = sp + 0.25, y = true), col = 'gray', size = 3.5) +
  geom_segment(aes(x = sp, xend = sp, y = low, yend = high),
	       size = 1, lineend = 'round') +
  facet_grid(vars(J.hb), vars(J.neon)) +
  theme_bw() +
  labs(x = 'Species', y = 'Estimate')

# Colonization covariate --------------------------------------------------

beta.gamma.1.low <- sapply(out.model, FUN = function(a) {a$q2.5$beta.gamma.1})
beta.gamma.1.med <- sapply(out.model, FUN = function(a) {a$q50$beta.gamma.1})
beta.gamma.1.high <- sapply(out.model, FUN = function(a) {a$q97.5$beta.gamma.1})

species.names <- c('Species 1', 'Species 2', 'Species 3', 'Species 4',
		   'Species 5', 'Species 6', 'Species 7', 'Species 8',
		   'Species 9', 'Species 10')
n.years.vals <- c(5, 10)
J.hb.vals <- c(1, 3)
J.neon.vals <- c(1, 3)
n.sp <- length(species.names)
# c goes down each row, then across columns.
beta.gamma.1.df <- data.frame(low = c(beta.gamma.1.low),
			med = c(beta.gamma.1.med),
			high = c(beta.gamma.1.high),
			sp = rep(1:10, n.sim * n.scenarios),
			n.years = rep(rep(n.years.vals, each = n.sp),
				      times = n.sim * n.scenarios / length(n.years.vals)),
                        J.hb = rep(rep(J.hb.vals, each = n.sp * length(n.years.vals)),
				   times = n.sim * n.scenarios / (length(n.years.vals) *
								  length(J.hb.vals))),
			J.neon = rep(rep(J.neon.vals, each = n.sp * length(n.years.vals) *
					 length(J.hb.vals)),
				     times = n.sim * n.scenarios / (length(n.years.vals) *
								    length(J.hb.vals) *
								    length(J.neon.vals))),
			true = rep(beta.gamma.1, n.sim * n.scenarios)
			)

plot.dat.beta.gamma.1 <- beta.gamma.1.df %>%
  group_by(n.years, J.hb, J.neon, sp, true) %>%
  summarize(low = mean(low),
	    med = mean(med),
	    high = mean(high)) %>%
  mutate(n.years = factor(n.years),
	 J.hb = factor(J.hb),
	 J.neon = factor(J.neon)) %>%
  ungroup()

plot.dat.beta.gamma.1 %>%
  filter(n.years == 5) %>%
ggplot(aes(x = sp, y = med)) +
  geom_point(size = 3.5, col = 'black') +
  geom_point(aes(x = sp + 0.25, y = true), col = 'gray', size = 3.5) +
  geom_segment(aes(x = sp, xend = sp, y = low, yend = high),
	       size = 1, lineend = 'round') +
  facet_grid(vars(J.hb), vars(J.neon)) +
  theme_bw() +
  labs(x = 'Species', y = 'Estimate')


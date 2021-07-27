# summary-sim-2.R: this script contains code to summarize the results
#                  of the secon simulation scenario (main-sim-2-icom.R and 
#                  main-sim-2-idm.R) and produce resulting figures.
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(nimble)
library(tidyverse)
library(ggthemes)
# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Read in community model results ----------------------------------------- 
load('results/sim-icom-results-2-100-simulations-2021-07-17.R')
# Results consist of the following objects: 
#     samples.list: list of coda MCMC objects for each simulation run of the ICOM
#     beta.0: true species and year specific intercepts
#     beta.1: true species covariate effect
#     phi: true species autologistic effect
#     alpha.0: true replicated data set detection intercepts
#     alpha.1: true replicated data set detection covariate effect
#     gamma.1.0: true nonreplicated data set detection intercept
#     gamma.1.1: true nonreplicated data set detection covariate effect
# Compute mean for each parameter and simulation
comm.mean <- sapply(samples.list, FUN = function(a) {apply(a, 2, mean)})
# Compute 2.5% quantile for Bayesian CI
comm.bci.low <- sapply(samples.list, FUN = function(a) {apply(a, 2, quantile, 0.025)})
# Compute 97.5% quantile for Bayesian CI
comm.bci.high <- sapply(samples.list, FUN = function(a) {apply(a, 2, quantile, 0.975)})
# Compute medians of Bayesian Credible Intervals
bci.low.med <- apply(comm.bci.low, 1, median)
bci.high.med <- apply(comm.bci.high, 1, median)
# Get 2.5% quantile of means
comm.low <- apply(comm.mean, 1, quantile, prob = 0.025)
# Get median of means
comm.med <- apply(comm.mean, 1, quantile, prob = 0.5)
# Get 97.5% quantile of means
comm.high <- apply(comm.mean, 1, quantile, prob = 0.975)
# Read in idm results ----------------------------------------------------
n.sp <- 25
n.params.tracked <- 28
params.low <- matrix(0, nrow = n.sp, ncol = n.params.tracked)
params.med <- matrix(0, nrow = n.sp, ncol = n.params.tracked)
params.high <- matrix(0, nrow = n.sp, ncol = n.params.tracked)
for (i in 1:n.sp) {
  print(i)
  curr.file <- paste('results/sim-idm-results-2-100-species-', i, 
		     '-simulations-2021-07-18.R', sep = '')
  load(curr.file)
  # Results consist of the following objects: 
  #     samples.list: list of coda MCMC objects for each simulation run of the IDM
  #     beta.0: true species and year specific intercepts
  #     beta.1: true species covariate effect
  #     phi: true species autologistic effect
  #     alpha.0: true replicated data set detection intercepts
  #     alpha.1: true replicated data set detection covariate effect
  #     gamma.1.0: true nonreplicated data set detection intercept
  #     gamma.1.1: true nonreplicated data set detection covariate effect
  # Compute mean for each parameter and simulation
  curr.mean <- sapply(samples.list, FUN = function(a) {apply(a, 2, mean)})
  # Get 25% quantile of means
  params.low[i, ] <- apply(curr.mean, 1, quantile, prob = 0.025)
  # Get median of means
  params.med[i, ] <- apply(curr.mean, 1, median)
  # Get 95% quantile of means
  params.high[i, ] <- apply(curr.mean, 1, quantile, prob = 0.975)
}

# Assign parameter names to the column names. 
colnames(params.low) <- rownames(curr.mean)
colnames(params.med) <- rownames(curr.mean)
colnames(params.high) <- rownames(curr.mean)

# Summary of species-specific effects -------------------------------------
# Get in long format ------------------------------------------------------
params.low <- data.frame(params.low)
params.low$species <- 1:nrow(params.low)
params.low <- params.low %>%
  pivot_longer(cols = -species, values_to = 'low', names_to = 'params')
params.med <- data.frame(params.med)
params.med$species <- 1:nrow(params.med)
params.med <- params.med %>%
  pivot_longer(cols = -species, values_to = 'med', names_to = 'params')
params.high <- data.frame(params.high)
params.high$species <- 1:nrow(params.high)
params.high <- params.high %>%
  pivot_longer(cols = -species, values_to = 'high', names_to = 'params')
params.full <- params.low
params.full$med <- params.med$med
params.full$high <- params.high$high
# Remove intercept
params.full <- params.full %>%
  filter(substr(params, 1, 8) != 'int.beta')
# Get community data in long format ---------------------------------------
n.params <- n_distinct(params.full$params)
# Lower 2.5% quantile -----------------
# Create data frame 
comm.low.df <- data.frame(low = comm.low, names = names(comm.low))
# Remove community average parameters
comm.low <- comm.low.df %>%
  filter(!str_detect(names, 'mean')) %>%
  mutate(species = rep(1:n.sp, times = n.params))
# Get parameter names
comm.low$params <- rep(c('alpha.0.1', 'alpha.0.2', 'alpha.0.3', 'alpha.0.4', 
		       'alpha.0.5', 'alpha.0.6', 'alpha.1', 
		       'gamma.1.0.1', 'gamma.1.0.2', 'gamma.1.0.3', 'gamma.1.0.4', 
		       'gamma.1.0.5', 'gamma.1.0.6', 'gamma.1.1', 
		       'beta.0.1', 'beta.0.2', 'beta.0.3', 'beta.0.4', 
		       'beta.0.5', 'beta.0.6', 'beta.1', 'phi'), each = n.sp)
# Repeat for all summary statistics
# Medians ----------------------------- 
comm.med.df <- data.frame(med = comm.med, names = names(comm.med))
comm.med <- comm.med.df %>%
  filter(!str_detect(names, 'mean')) %>%
  mutate(species = rep(1:n.sp, times = n.params))
comm.med$params <- rep(c('alpha.0.1', 'alpha.0.2', 'alpha.0.3', 'alpha.0.4', 
		       'alpha.0.5', 'alpha.0.6', 'alpha.1', 
		       'gamma.1.0.1', 'gamma.1.0.2', 'gamma.1.0.3', 'gamma.1.0.4', 
		       'gamma.1.0.5', 'gamma.1.0.6', 'gamma.1.1', 
		       'beta.0.1', 'beta.0.2', 'beta.0.3', 'beta.0.4', 
		       'beta.0.5', 'beta.0.6', 'beta.1', 'phi'), each = n.sp)
# Upper 97.5% quantile ----------------
comm.high.df <- data.frame(high = comm.high, names = names(comm.high))
comm.high <- comm.high.df %>%
  filter(!str_detect(names, 'mean')) %>%
  mutate(species = rep(1:n.sp, times = n.params))
comm.high$params <- rep(c('alpha.0.1', 'alpha.0.2', 'alpha.0.3', 'alpha.0.4', 
		       'alpha.0.5', 'alpha.0.6', 'alpha.1', 
		       'gamma.1.0.1', 'gamma.1.0.2', 'gamma.1.0.3', 'gamma.1.0.4', 
		       'gamma.1.0.5', 'gamma.1.0.6', 'gamma.1.1', 
		       'beta.0.1', 'beta.0.2', 'beta.0.3', 'beta.0.4', 
		       'beta.0.5', 'beta.0.6', 'beta.1', 'phi'), each = n.sp)
# Combine into single data frame
comm.all <- comm.low
comm.all$med <- comm.med$med
comm.all$high <- comm.high$high
names(comm.all) <- c('comm.low', 'names', 'species', 'params', 'comm.med', 'comm.high')
params.full$params <- sub("[.]$", "", params.full$params)
# Join idm and icom together ----------------------------------------------
full.df <- left_join(params.full, comm.all, by = c('species', 'params'))


# Plot for intercept ------------------------------------------------------
my.colors = c('IDM' = '#23708EFF', 'ICOM' = '#55C667FF', 'True' = 'black')
# Create a plot for a single intercept parameter
full.df %>%
  filter(params == 'beta.0.2') %>%
  mutate(species = as.numeric(factor(species, levels = order(beta.0[, 2])))) %>%
ggplot(aes(x = species, y = comm.med)) + 
  geom_point(aes(col = 'ICOM'), size = 3) + 
  theme_bw(base_size = 21) + 
  geom_segment(aes(x = species, y = comm.low, xend = species, yend = comm.high, 
		   col = 'ICOM'), 
	       lineend = 'round', size = 1.2) + 
  geom_point(aes(x = species - 0.25, y = med, col = 'IDM'), size = 3) + 
  geom_segment(aes(x = species - 0.25, y = low, xend = species - 0.25, yend = high, 
		   col = 'IDM'), 
	       lineend = 'round', size = 1.2) + 
  geom_point(aes(x = species + 0.25, y = beta.0[, 2], col = 'True'), size = 3) + 
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25), labels = c(1, 5, 10, 15, 20, 25)) + 
  scale_color_manual(values = my.colors) + 
  labs(x = 'Species', y = 'Occupancy Intercept', color = 'Model') + 
  theme(legend.position = c(0.9, 0.3))  
# Save plot if desired
# Figure S2
ggsave(device = 'pdf', filename = 'figures/sim-4-int-3.pdf', 
       height = 5, width = 10, units = 'in')

# Compare all parameters and summarize in table ---------------------------
full.df <- full.df %>%
  arrange(params, species)	
full.df <- full.df %>%
  mutate(true = c(alpha.0[, 1], alpha.0[, 2], alpha.0[, 3], 
		  alpha.0[, 4], alpha.0[, 5], alpha.0[, 6], 
		  alpha.1, 
		  beta.0[, 1], beta.0[, 2], beta.0[, 3], beta.0[, 4], beta.0[, 5], beta.0[, 6], 
		  beta.1, 
		  gamma.1.0[, 1], gamma.1.0[, 2], gamma.1.0[, 3], 
		  gamma.1.0[, 4], gamma.1.0[, 5], gamma.1.0[, 6], 
		  gamma.1.1, 
		  phi))
# Results displayed in Table 2. 
full.df %>%
  #mutate(species = as.numeric(factor(species, levels = order(beta.0)))) %>%
  group_by(species, params) %>%
  summarize(perc.precision = (1 - (comm.high - comm.low) / (high - low)) * 100, 
	    comm.bias = abs(comm.med - true), 
	    sp.bias = abs(med - true)) %>%
  arrange(params, species) %>%
  group_by(params) %>%
  summarize(avg.perc.precision = mean(perc.precision), 
	    avg.comm.bias = mean(comm.bias), 
	    avg.sp.bias = mean(sp.bias)) %>%
  arrange(desc(avg.perc.precision)) %>%
  print(n = nrow(.))

# Compute Bayesian 95% CI Coverage for all species specific parameters ----
true <- c(alpha.0[, 1], alpha.0[, 2], alpha.0[, 3],
		  alpha.0[, 4], alpha.0[, 5], alpha.0[, 6],
		  alpha.1,
		  gamma.1.0[, 1], gamma.1.0[, 2], gamma.1.0[, 3],
		  gamma.1.0[, 4], gamma.1.0[, 5], gamma.1.0[, 6],
		  gamma.1.1, beta.0[, 1],
		  beta.0[, 2], beta.0[, 3], beta.0[, 4], beta.0[, 5], beta.0[, 6],
		  beta.1, phi)
curr.low <- comm.bci.low[!str_detect(names(bci.low.med), 'mean'), ]
curr.high <- comm.bci.high[!str_detect(names(bci.high.med), 'mean'), ]

coverage.mat <- matrix(NA, nrow = nrow(curr.low), ncol = ncol(curr.low))
for (i in 1:nrow(coverage.mat)) {
  coverage.mat[i, ] <- ifelse((true[i] > curr.low[i, ]) & (true[i] < curr.high[i, ]), 1, 0) 
}
# Average Bayesian 95% CI for species-specific parameters using the ICOM. 
mean(coverage.mat)

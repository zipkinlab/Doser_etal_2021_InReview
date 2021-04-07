# summary-sim-3.R: this script contains code to summarize the results
#                  of the third simulation scenario and produce 
#                  resulting figures

rm(list = ls())
library(coda)
library(nimble)
library(tidyverse)
library(ggthemes)
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Simulation scenarios ----------------------------------------------------
# Simulated number of sites with NEON data
J.neon.vals <- c(10, 30, 50)
# Simulated number of sites with HB data
J.det.vals <- c(10, 30, 50)
# Simulated number of BBS routes
n.route.vals <- c(1, 3, 5)
# Simulated number of eBird cells
J.eb.vals <- c(25)
n.scenarios <- length(J.neon.vals) * length(J.eb.vals) * length(n.route.vals) *
  length(J.det.vals)
param.vals <- expand.grid(J.neon.vals, J.eb.vals, n.route.vals, J.det.vals)
names(param.vals) <- c('J.neon', 'J.eb', 'n.route', 'J.det')

# Read in results one at a time -------------------------------------------
# Number of parameters monitored
n.params.tracked <- 60
params.lowest <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.low <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.med <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.high <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.highest <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
for (i in 1:n.scenarios) {
  print(i)
  # Read in results
  curr.file <- paste('results/sim-icm-results-1-100-condition-', i, 
		     '-simulations-2021-03-20.R', sep = '')
  load(curr.file)
  # Compute mean for each parameter and simulation
  curr.mean <- sapply(samples.list, FUN = function(a) {apply(a, 2, mean)})
  # Get 25% quantile of means
  params.low[i, ] <- apply(curr.mean, 1, quantile, prob = 0.25)
  # Get median of means
  params.med[i, ] <- apply(curr.mean, 1, median)
  # Get 95% quantile of means
  params.high[i, ] <- apply(curr.mean, 1, quantile, prob = 0.75)
}

colnames(params.low) <- rownames(curr.mean)
colnames(params.med) <- rownames(curr.mean)
colnames(params.high) <- rownames(curr.mean)

# Get in long format ------------------------------------------------------
# Lower 25% quantile
params.low <- data.frame(params.low)
params.low$scenario <- 1:nrow(params.low)
# Get in long format
params.low <- params.low %>%
  select(-contains('mean')) %>% # remove community level parameters
  pivot_longer(cols = -scenario, values_to = 'low', names_to = 'param')
# Assign species identity to each row
params.low$sp <- str_split(params.low$param, '\\.') %>%
  map_chr(~ .[4])
params.low$sp <- paste('Species ', params.low$sp, sep = '')
# Assign parameter to each row
params.low$param <- str_split(params.low$param, '\\.') %>%
  map_chr(~ paste0(.[1:3], collapse = '-'))
# Change order
params.low <- params.low %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10')))
# Same process below for other summary measures
# Median
params.med <- data.frame(params.med)
params.med$scenario <- 1:nrow(params.med)
params.med <- params.med %>%
  select(-contains('mean')) %>%
  pivot_longer(cols = -scenario, values_to = 'med', names_to = 'param')
params.med$sp <- str_split(params.med$param, '\\.') %>%
  map_chr(~ .[4])
params.med$sp <- paste('Species ', params.med$sp, sep = '')
params.med$param <- str_split(params.med$param, '\\.') %>%
  map_chr(~ paste0(.[1:3], collapse = '-'))
params.med <- params.med %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10')))
# Upper 75% quantile
params.high <- data.frame(params.high)
params.high$scenario <- 1:nrow(params.high)
params.high <- params.high %>%
  select(-contains('mean')) %>%
  pivot_longer(cols = -scenario, values_to = 'high', names_to = 'param')
params.high$sp <- str_split(params.high$param, '\\.') %>%
  map_chr(~ .[4])
params.high$sp <- paste('Species ', params.high$sp, sep = '')
params.high$param <- str_split(params.high$param, '\\.') %>%
  map_chr(~ paste0(.[1:3], collapse = '-'))
params.high <- params.high %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10')))
# Combine all summary statistics together
params.full <- params.low
params.full$med <- params.med$med
params.full$high <- params.high$high
# Assign info to params.full DF on what scenario each value corresponds to. 
params.full <- params.full %>%
  mutate(J.neon = param.vals$J.neon[scenario], 
	 n.route = param.vals$n.route[scenario], 
	 J.det = param.vals$J.det[scenario])
# Add true values for plotting
true.vals <- c(beta.gamma.0, beta.gamma.1, beta.phi.0, beta.phi.1, beta.psi.0, 
	       beta.psi.1)
params.full <- params.full %>%
  mutate(true = rep(true.vals, times = n.scenarios))

# Subtract true values from estimated values to get estimated bias
params.grouped <- params.full %>%
  group_by(scenario, param, J.neon, n.route, J.det) %>%
  summarize(low = mean(low - true), 
	    med = mean(med - true), 
	    high = mean(high - true)) %>%
  ungroup()
# Compute IQR for boxplots
params.grouped <- params.grouped %>%
  mutate(iqr.lowest = med - (high - low) * 1.5,
	 iqr.highest = med + (high - low) * 1.5)

# Some summary plots ------------------------------------------------------
my.colors = c('10' = 'lightsalmon2', 
	      '30' = 'gray', '50' = 'lightskyblue2')
# Just rearrange the data a bit to help with plot formatting
params.grouped %>%
  mutate(J.det = ifelse(J.neon == 10, J.det - 5, ifelse(J.neon == 30, J.det, J.det + 5))) %>%
  mutate(J.neon = factor(J.neon), 
	 n.route = factor(ifelse(n.route == 1, '10 BBS Stops', ifelse(n.route == 3, '30 BBS Stops', 
								     '50 BBS Stops')), 
			  levels = c('10 BBS Stops', '30 BBS Stops', '50 BBS Stops'))) %>%
  filter(param == 'beta-phi-1') %>% 
  ggplot(aes(x = J.det, y = med, fill = J.neon)) + 
    geom_boxplot(aes(group = J.det, ymin = iqr.lowest,
		     lower = low, middle = med, upper = high,
		     ymax = iqr.highest), stat = 'identity', size = 0.9) +
  facet_wrap(vars(n.route), ncol = 3) + 
  theme_bw(base_size = 18) + 
  scale_fill_manual(values = my.colors) + 
  geom_vline(xintercept = 20, linetype = 2, col = 'grey') +
  geom_vline(xintercept = 40, linetype = 2, col = 'grey')  +
  geom_hline(yintercept = 0, linetype = 2) + 
  scale_x_continuous(breaks = c(10, 30, 50), labels = c(10, 30, 50)) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank()) +  
  labs(x = 'Number of DET sites', 
       y = 'Average bias in persistence covariate effect', 
       col = 'Number of NEON sites') + 
  theme(legend.position = c(0.94, 0.1)) +
  ggsave(device = 'pdf', filename = 'figures/sim-1-persistence-cov.pdf', 
	 height = 8, width = 13, units = 'in')

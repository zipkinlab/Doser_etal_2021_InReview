# summary-sim-2.R: this script contains code to summarize the results
#                  of the first simulation scenario and produce
#                  resulting figures.

rm(list = ls())
library(coda)
library(nimble)
library(tidyverse)
library(ggthemes)
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Simulation scenarios ----------------------------------------------------
# Is eBird included in the model?
eb.in.vals <- c(FALSE, TRUE)
# Is DET included in the model?
det.in.vals <- c(FALSE, TRUE)
# Is NEON included in the model?
neon.in.vals <- c(FALSE, TRUE)
# Is BBS included in the model?
bbs.in.vals <- c(FALSE, TRUE)
# The minus one eliminates the case when no data sets are available
n.scenarios <- length(eb.in.vals) * length(det.in.vals) * length(neon.in.vals) *
  length(bbs.in.vals) - 1
# Data frame containing different conditions
param.vals <- expand.grid(eb.in.vals, det.in.vals, neon.in.vals, bbs.in.vals)
names(param.vals) <- c('eb.in', 'det.in', 'neon.in', 'bbs.in')
param.vals <- param.vals[-1, ]

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
  curr.file <- paste('results/sim-icm-results-3-100-condition-', i, 
		     '-simulations-2021-03-18.R', sep = '')
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

# Get data in long format -------------------------------------------------
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
  mutate(neon.in = param.vals$neon.in[scenario], 
	 bbs.in = param.vals$bbs.in[scenario], 
	 det.in = param.vals$det.in[scenario], 
	 eb.in = param.vals$eb.in[scenario])
# Add true values for plotting
true.vals <- c(beta.gamma.0, beta.gamma.1, beta.phi.0, beta.phi.1, beta.psi.0, 
	       beta.psi.1)
params.full <- params.full %>%
  mutate(true = rep(true.vals, times = n.scenarios))
# Compute IQR for boxplots
params.full <- params.full %>%
  mutate(iqr.lowest = med - (high - low) * 1.5, 
	 iqr.highest = med + (high - low) * 1.5)

# Some summary plots ------------------------------------------------------
# Rectangles for background colors
rects <- data.frame(xstart = c(0, 4.5, 10.5, 14.5), 
		    xend = c(4.5, 10.5, 14.5, 16), 
		    col = letters[1:4])

# Sorted order
# (1) eBird only, (2) DET only, (3) NEON only, (4) BBS only, (5) eb + DET, 
# (6) eBird + NEON, (7) DET + NEON, (8) eBird + BBS, (9) DET + BBS, (10) NEON + BBS, 
# (11) eBird + DET + NEON, (12) eBird + DET + BBS, (13) eBird + NEON + BBS, 
# (14) DET + NEON + BBS, (15) eBird + DET + NEON + BBS
# Transform intercept parameters to be on real scale. 
plotting.dat <- params.full %>%
  mutate(low = ifelse(param %in% c('beta-gamma-0', 'beta-phi-0', 'beta-psi-0'), 
		      logit.inv(low), low), 
	 med = ifelse(param %in% c('beta-gamma-0', 'beta-phi-0', 'beta-psi-0'), 
		      logit.inv(med), med), 
	 high = ifelse(param %in% c('beta-gamma-0', 'beta-phi-0', 'beta-psi-0'), 
		      logit.inv(high), high), 
	 iqr.lowest = ifelse(param %in% c('beta-gamma-0', 'beta-phi-0', 'beta-psi-0'), 
		      logit.inv(iqr.lowest), iqr.lowest), 
	 iqr.highest = ifelse(param %in% c('beta-gamma-0', 'beta-phi-0', 'beta-psi-0'), 
		      logit.inv(iqr.highest), iqr.highest), 
	 true = ifelse(param %in% c('beta-gamma-0', 'beta-phi-0', 'beta-psi-0'), 
		      logit.inv(true), true) 
	) %>%
  # Change this parameter in the filter and order statements to show different 
  # parameters in the plot
  filter(param == 'beta-phi-1') %>% 
  mutate(sp = factor(sp, levels = paste('Species ', order(-beta.phi.1), sep = '')), 
   	 scenario = as.numeric(factor(scenario, levels = order(apply(param.vals, 1, sum)))))
# Create plot
ggplot(data = plotting.dat) +
  geom_rect(data = rects, aes(xmin = xstart, 
  			      xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  	    alpha = 0.2) + 
  scale_fill_viridis_d() + 
  geom_boxplot(data = plotting.dat, aes(x = scenario, group = scenario, 
					ymin = iqr.lowest,
					lower = low, middle = med, upper = high, 
					ymax = iqr.highest), stat = 'identity', fill = 'gray') + 
  facet_wrap(vars(sp), ncol = 3) + 
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 4.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 10.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 14.5, linetype = 2, col = 'grey') + 
  geom_hline(aes(yintercept = true), linetype = 2) + 
  guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = 'Data Sets Used', 
       y = 'Persistence Covariate Effect') + 
  scale_x_continuous(expand = c(0, 0), 
		     breaks = 1:15, labels = c('EBIRD', 'DET', 'NEON', 'BBS', 
					       'EBIRD + DET', 'EBIRD + NEON', 'DET + NEON', 
					       'EBIRD + BBS', 'DET + BBS', 'NEON + BBS', 
					       'EBIRD + DET + NEON', 'EBIRD + DET + BBS', 
					       'EBIRD + NEON + BBS', 'DET + NEON + BBS', 
					       'EBIRD + DET + NEON + BBS')) +
  ggsave(device = 'pdf', filename = 'figures/sim-3-persistence-cov.pdf', 
	 height = 11, width = 14, units = 'in')


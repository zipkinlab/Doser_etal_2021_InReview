# summary-sim-1.R: this script contains code to summarize the results
#                  of the first simulation scenario (main-sim-1-icom.R) 
#                  and produce resulting figures. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(nimble)
library(tidyverse)
library(ggthemes)
library(ggpubr)
# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Read in results one at a time -------------------------------------------
# Number of parameters monitored
n.params.tracked <- 208
# Number of simulation scenarios (i.e., number of different combinations 
#                                       of the three data sources). 
n.scenarios <- 7
# Number of years simulated
n.years <- 6
# Initialize matrices to contain results summaries
params.lowest <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.low <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.med <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.mean <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.sd <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.high <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.highest <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
n.sims <- 100
# NOTE: simulation results files are too large to store on GitHub and 
#       are required to run this file. Results 
#       can be generated using the main-sim-1-icom.R file, or contact
#       the first author for the files (doserjef@msu.edu) for a link to 
#       a Dropbox folder. 
# Iteratively read in simulation results from each of seven models
for (i in 1:n.scenarios) {
  print(i)
  # Read in results
  curr.file <- paste('results/sim-icom-results-1-100-condition-', i, 
		     '-simulations-', 1, '-chain-2021-08-26.R', sep = '')
  load(curr.file)
  samples.a <- samples.list
  curr.file <- paste('results/sim-icom-results-1-100-condition-', i, 
		     '-simulations-', 2, '-chain-2021-08-26.R', sep = '')
  load(curr.file)
  samples.b <- samples.list
  curr.file <- paste('results/sim-icom-results-1-100-condition-', i, 
		     '-simulations-', 3, '-chain-2021-08-26.R', sep = '')
  load(curr.file)
  samples.c <- samples.list
  samples.list <- list()
  for (a in 1:n.sims) {
    samples.list[[a]] <- mcmc.list(samples.a[[a]], samples.b[[a]], samples.c[[a]])
  }
  # For convergence assessment. Takes a while. 
  # r.hat.vals <- sapply(samples.list, gelman.diag)
  # Each results file consists of a list of MCMC objects, called samples.list
  # for each simulated run of the individual model. 
  # Compute median for each parameter and simulation 
  curr.mean <- sapply(samples.list, FUN = function(a) {apply(do.call('rbind', a), 2, median)})
  # Get 2.5% quantile of medians
  params.low[i, ] <- apply(curr.mean, 1, quantile, prob = 0.025)
  # Get median of medians
  params.med[i, ] <- apply(curr.mean, 1, median)
  # Get mean of medians
  params.mean[i, ] <- apply(curr.mean, 1, mean)
  # Compute standard deviation of medians
  params.sd[i, ] <- apply(curr.mean, 1, sd)
  # Compute 97.5% quantile of medians
  params.high[i, ] <- apply(curr.mean, 1, quantile, prob = 0.975)
}
# Assign parameter names to the column names
colnames(params.low) <- rownames(curr.mean)
colnames(params.med) <- rownames(curr.mean)
colnames(params.mean) <- rownames(curr.mean)
colnames(params.sd) <- rownames(curr.mean)
colnames(params.high) <- rownames(curr.mean)
# Create a copy of parameter values for summary of community level and 
# species-level parameters
params.low.comm <- params.low
params.med.comm <- params.med
params.mean.comm <- params.mean
params.sd.comm <- params.sd
params.high.comm <- params.high

# Species-specific parameter summary --------------------------------------
# Get data in long format -------------------------------------------------
# Lower 2.5% quantile ------------------
I <- 25
params.low <- data.frame(params.low)
params.low$scenario <- 1:nrow(params.low)
# Get in long format
params.low <- params.low %>%
  select(-contains('mean')) %>% # remove community level parameters
  pivot_longer(cols = -scenario, values_to = 'low', names_to = 'param')
n.params <- length(unique(params.low$param))
# Assign species identity to each row
params.low$sp <- rep(1:I, times = n.params / I * n.scenarios)
params.low$sp <- paste('Species ', params.low$sp, sep = '')
# Assign parameter to each row
params.low$param <- str_split(params.low$param, '\\.') %>%
  map_chr(~ paste0(.[1:2], collapse = '-'))
params.low$param <- ifelse(substr(params.low$param, 1, 3) == 'phi', 'phi', 
			   params.low$param)
# Get order of species correct
params.low <- params.low %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10', 'Species 11', 
				    'Species 12', 'Species 13', 'Species 14', 
				    'Species 15', 'Species 16', 'Species 17', 
				    'Species 18', 'Species 19', 'Species 20', 
				    'Species 21', 'Species 22', 'Species 23', 
				    'Species 24', 'Species 25')))
# Same process below for other summary measures
# Median  -----------------------------
params.med <- data.frame(params.med)
params.med$scenario <- 1:nrow(params.med)
# Get in long format
params.med <- params.med %>%
  select(-contains('mean')) %>% # remove community level parameters
  pivot_longer(cols = -scenario, values_to = 'med', names_to = 'param')
n.params <- length(unique(params.med$param))
# Assign species identity to each row
params.med$sp <- rep(1:I, times = n.params / I * n.scenarios)
params.med$sp <- paste('Species ', params.med$sp, sep = '')
# Assign parameter to each row
params.med$param <- str_split(params.med$param, '\\.') %>%
  map_chr(~ paste0(.[1:2], collapse = '-'))
params.med$param <- ifelse(substr(params.med$param, 1, 3) == 'phi', 'phi', 
			   params.med$param)
# Change order
params.med <- params.med %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10', 'Species 11', 
				    'Species 12', 'Species 13', 'Species 14', 
				    'Species 15', 'Species 16', 'Species 17', 
				    'Species 18', 'Species 19', 'Species 20', 
				    'Species 21', 'Species 22', 'Species 23', 
				    'Species 24', 'Species 25')))
# Upper 97.5% quantile ------------------
params.high <- data.frame(params.high)
params.high$scenario <- 1:nrow(params.high)
# Get in long format
params.high <- params.high %>%
  select(-contains('mean')) %>% # remove community level parameters
  pivot_longer(cols = -scenario, values_to = 'high', names_to = 'param')
n.params <- length(unique(params.high$param))
# Assign species identity to each row
params.high$sp <- rep(1:I, times = n.params / I * n.scenarios)
params.high$sp <- paste('Species ', params.high$sp, sep = '')
# Assign parameter to each row
params.high$param <- str_split(params.high$param, '\\.') %>%
  map_chr(~ paste0(.[1:2], collapse = '-'))
params.high$param <- ifelse(substr(params.high$param, 1, 3) == 'phi', 'phi', 
			   params.high$param)
# Change order
params.high <- params.high %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10', 'Species 11', 
				    'Species 12', 'Species 13', 'Species 14', 
				    'Species 15', 'Species 16', 'Species 17', 
				    'Species 18', 'Species 19', 'Species 20', 
				    'Species 21', 'Species 22', 'Species 23', 
				    'Species 24', 'Species 25')))
# Mean --------------------------------
params.mean <- data.frame(params.mean)
params.mean$scenario <- 1:nrow(params.mean)
# Get in long format
params.mean <- params.mean %>%
  select(-contains('mean')) %>% # remove community level parameters
  pivot_longer(cols = -scenario, values_to = 'mean', names_to = 'param')
n.params <- length(unique(params.mean$param))
# Assign species identity to each row
params.mean$sp <- rep(1:I, times = n.params / I * n.scenarios)
params.mean$sp <- paste('Species ', params.mean$sp, sep = '')
# Assign parameter to each row
params.mean$param <- str_split(params.mean$param, '\\.') %>%
  map_chr(~ paste0(.[1:2], collapse = '-'))
params.mean$param <- ifelse(substr(params.mean$param, 1, 3) == 'phi', 'phi', 
			   params.mean$param)
# Change order
params.mean <- params.mean %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10', 'Species 11', 
				    'Species 12', 'Species 13', 'Species 14', 
				    'Species 15', 'Species 16', 'Species 17', 
				    'Species 18', 'Species 19', 'Species 20', 
				    'Species 21', 'Species 22', 'Species 23', 
				    'Species 24', 'Species 25')))
# Standard deviation ------------------
params.sd <- data.frame(params.sd)
params.sd$scenario <- 1:nrow(params.sd)
# Get in long format
params.sd <- params.sd %>%
  select(-contains('mean')) %>% # remove community level parameters
  pivot_longer(cols = -scenario, values_to = 'sd', names_to = 'param')
n.params <- length(unique(params.sd$param))
# Assign species identity to each row
params.sd$sp <- rep(1:I, times = n.params / I * n.scenarios)
params.sd$sp <- paste('Species ', params.sd$sp, sep = '')
# Assign parameter to each row
params.sd$param <- str_split(params.sd$param, '\\.') %>%
  map_chr(~ paste0(.[1:2], collapse = '-'))
params.sd$param <- ifelse(substr(params.sd$param, 1, 3) == 'phi', 'phi', 
			   params.sd$param)
# Change order
params.sd <- params.sd %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10', 'Species 11', 
				    'Species 12', 'Species 13', 'Species 14', 
				    'Species 15', 'Species 16', 'Species 17', 
				    'Species 18', 'Species 19', 'Species 20', 
				    'Species 21', 'Species 22', 'Species 23', 
				    'Species 24', 'Species 25')))
# Combine all summary statistics together in a single data frame
params.full <- params.low
params.full$med <- params.med$med
params.full$high <- params.high$high
params.full$mean <- params.mean$mean
params.full$sd <- params.sd$sd

# Take averages of the individual species averages
params.grouped <- params.full %>%
  group_by(scenario, param) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    mean = mean(mean),
	    sd = mean(sd),
	    high = mean(high)) %>%
  ungroup()

# Species-specific plots --------------------------------------------------
# Rectangles for background colors
rects <- data.frame(xstart = c(0, 3.5, 6.5), 
		    xend = c(3.5, 6.5, 8), 
		    col = letters[1:3])

# Intercept ---------------------------
plotting.dat <- params.grouped %>%
  filter(param == 'beta-0')

# Get data in proper format for plot. Compute interquartile range if 
# desired for boxplots
plotting.dat <- plotting.dat %>%
  mutate(iqr.lowest = med - (high - low) * 1.5, # lower integer quartile range
	 iqr.highest = med + (high - low) * 1.5,  # upper inter quartile range
	 # Creates variable based on how many data sets are in a given scenario
	 group = ifelse(scenario < 4, 1, ifelse(scenario < 7, 2, 3)))
# Need to remove the one nonreplicated data set because of unrealistic 
# estimates due to nonconvergence and nonidentifiability. 
plotting.dat <- plotting.dat %>%
  filter(scenario != 2)
# Create plot for intercept
int.sp.plot <- ggplot(data = plotting.dat) +
  geom_rect(data = rects, aes(xmin = xstart, 
  			      xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  	    alpha = 0.2) + 
  scale_fill_viridis_d() + 
  geom_point(aes(x = scenario, y = med), size = 4.2) +
  geom_segment(aes(x = scenario, xend = scenario, y = low,
 		   yend = high), size = 1.2, lineend = 'round') +
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 3.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 6.5, linetype = 2, col = 'grey') + 
  geom_hline(yintercept = 0, linetype = 2) + 
  guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = 'Data Sets Used', 
       y = 'Bias in intercept') + 
  scale_x_continuous(expand = c(0, 0), 
		     breaks = 1:7, labels = c('REP', 'NREP L', 'NREP H', 
					      'REP + NREP L', 
					      'REP + NREP H', 'NREP L + NREP H', 
					      'REP + NREP L + NREP H'))

# Covariate ---------------------------
plotting.dat <- params.grouped %>%
  filter(param == 'beta-1')

plotting.dat <- plotting.dat %>%
  mutate(iqr.lowest = med - (high - low) * 1.5,
	 iqr.highest = med + (high - low) * 1.5, 
	 group = ifelse(scenario < 5, 1, ifelse(scenario < 10, 2, 
						ifelse(scenario < 15, 3, 4))))
cov.sp.plot <- ggplot(data = plotting.dat) +
   geom_rect(data = rects, aes(xmin = xstart, 
  			      xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  	    alpha = 0.2) + 
  scale_fill_viridis_d() + 
  geom_point(aes(x = scenario, y = med), size = 4.2) +
  geom_segment(aes(x = scenario, xend = scenario, y = low,
 		   yend = high), size = 1.2, lineend = 'round') +
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 3.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 6.5, linetype = 2, col = 'grey') + 
  geom_hline(yintercept = 0, linetype = 2) + 
  guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = 'Data Sets Used', 
       y = 'Bias in covariate effect ') + 
  scale_x_continuous(expand = c(0, 0), 
		     breaks = 1:7, labels = c('REP', 'NREP L', 'NREP H', 
					      'REP + NREP L', 
					      'REP + NREP H', 'NREP L + NREP H', 
					      'REP + NREP L + NREP H'))

# Put plot together and save to hard drive if desired.
# Figure 2. 
ggarrange(int.sp.plot, cov.sp.plot, labels = c('A', 'B'), 
	  font.label = list(size = 18)) 
  ggsave(device = 'pdf', filename = 'figures/Fig2.pdf', 
	 width = 12, height = 7, units = 'in')

# Summary statistics of simulation results --------------------------------
params.grouped %>%
  mutate(param = factor(param, levels = c('beta-0', 'beta-1', 
					  'phi'),
			ordered = TRUE)) %>%
  select(scenario, param, mean, sd) %>%
  arrange(scenario, param) %>%
  print(n = nrow(.))

# Community level parameter summary ---------------------------------------
# Get data in long format -------------------------------------------------
# Lower 2.5% quantile ------------------
params.low.comm <- data.frame(params.low.comm)
params.low.comm$scenario <- 1:nrow(params.low.comm)
# Get in long format
params.low.comm <- params.low.comm %>%
  select(contains('mean'), 'scenario') %>% # community level parameters
  pivot_longer(cols = -scenario, values_to = 'low', names_to = 'param')

# Median ------------------------------
params.med.comm <- data.frame(params.med.comm)
params.med.comm$scenario <- 1:nrow(params.med.comm)
params.med.comm <- params.med.comm %>%
  select(contains('mean'), 'scenario') %>% # community level parameters
  pivot_longer(cols = -scenario, values_to = 'med', names_to = 'param')

# Upper 97.5% quantile ------------------
params.high.comm <- data.frame(params.high.comm)
params.high.comm$scenario <- 1:nrow(params.high.comm)
params.high.comm <- params.high.comm %>%
  select(contains('mean'), 'scenario') %>% # community level parameters
  pivot_longer(cols = -scenario, values_to = 'high', names_to = 'param')

# Mean --------------------------------
params.mean.comm <- data.frame(params.mean.comm)
params.mean.comm$scenario <- 1:nrow(params.mean.comm)
params.mean.comm <- params.mean.comm %>%
  select(contains('mean'), 'scenario') %>% # community level parameters
  pivot_longer(cols = -scenario, values_to = 'mean', names_to = 'param')

# Standard deviation ------------------
params.sd.comm <- data.frame(params.sd.comm)
params.sd.comm$scenario <- 1:nrow(params.sd.comm)
params.sd.comm <- params.sd.comm %>%
  select(contains('mean'), 'scenario') %>% # community level parameters
  pivot_longer(cols = -scenario, values_to = 'sd', names_to = 'param')

# Combine all summary statistics together 
params.full.comm <- params.low.comm
params.full.comm$med <- params.med.comm$med
params.full.comm$high <- params.high.comm$high
params.full.comm$mean <- params.mean.comm$mean
params.full.comm$sd <- params.sd.comm$sd

# Average across the many intercepts
params.full.comm$unique.param <- rep(c('beta.1.mean', 
				     rep('int.beta.mean', n.years), 'phi.mean'), 
				     times = n.scenarios)
params.full.comm <- params.full.comm %>%
  group_by(unique.param, scenario) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    high = mean(high), 
	    mean = mean(mean), 
	    sd = mean(sd)) %>%
  arrange(scenario, unique.param)

# Some summary plots ------------------------------------------------------
# Rectangles for background colors
rects <- data.frame(xstart = c(0, 3.5, 6.5), 
		    xend = c(3.5, 6.5, 8), 
		    col = letters[1:3])

# Intercept ---------------------------
plotting.dat <- params.full.comm %>%
  filter(unique.param == 'int.beta.mean')

# Get data in proper format for plot
plotting.dat <- plotting.dat %>%
  mutate(iqr.lowest = med - (high - low) * 1.5, # lower integer quartile range
	 iqr.highest = med + (high - low) * 1.5,  # upper inter quartile range
	 # Creates variable based on how many data sets are in a given scenario
	 group = ifelse(scenario < 4, 1, ifelse(scenario < 7, 2, 3)))
# Create plot for intercept
int.comm.plot <- ggplot(data = plotting.dat) +
  geom_rect(data = rects, aes(xmin = xstart, 
  			      xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  	    alpha = 0.2) + 
  scale_fill_viridis_d() + 
  geom_point(aes(x = scenario, y = med), size = 4.2) +
  geom_segment(aes(x = scenario, xend = scenario, y = low,
 		   yend = high), size = 1.2, lineend = 'round') +
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 3.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 6.5, linetype = 2, col = 'grey') + 
  geom_hline(yintercept = 0, linetype = 2) + 
  guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = 'Data Sets Used', 
       y = 'Bias in intercept') + 
  scale_x_continuous(expand = c(0, 0), 
		     breaks = 1:7, labels = c('REP', 'NREP L', 'NREP H', 
					      'REP + NREP L', 
					      'REP + NREP H', 'NREP L + NREP H', 
					      'REP + NREP L + NREP H'))

# Covariate effect --------------------
plotting.dat <- params.full.comm %>%
  filter(unique.param == 'beta.1.mean')

plotting.dat <- plotting.dat %>%
  mutate(iqr.lowest = med - (high - low) * 1.5,
	 iqr.highest = med + (high - low) * 1.5, 
	 group = ifelse(scenario < 5, 1, ifelse(scenario < 10, 2, 
						ifelse(scenario < 15, 3, 4))))
cov.comm.plot <- ggplot(data = plotting.dat) +
   geom_rect(data = rects, aes(xmin = xstart, 
  			      xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  	    alpha = 0.2) + 
  scale_fill_viridis_d() + 
  geom_point(aes(x = scenario, y = med), size = 4.2) +
  geom_segment(aes(x = scenario, xend = scenario, y = low,
 		   yend = high), size = 1.2, lineend = 'round') +
  theme_bw(base_size = 18) + 
  geom_vline(xintercept = 3.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 6.5, linetype = 2, col = 'grey') + 
  geom_hline(yintercept = 0, linetype = 2) + 
  guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = 'Data Sets Used', 
       y = 'Bias in covariate effect') + 
  scale_x_continuous(expand = c(0, 0), 
		     breaks = 1:7, labels = c('REP', 'NREP L', 'NREP H', 
					      'REP + NREP L', 
					      'REP + NREP H', 'NREP L + NREP H', 
					      'REP + NREP L + NREP H'))
# Put plot together and save to hard drive if desired.
# Figure S1. 
ggarrange(int.comm.plot, cov.comm.plot, labels = c('A', 'B'), 
	  font.label = list(size = 18)) 
ggsave(device = 'pdf', filename = 'figures/FigS1.pdf', 
       width = 12, height = 7, units = 'in')

# Summary statistics for further summary of results -----------------------
params.full.comm %>%
  print(n = nrow(.))

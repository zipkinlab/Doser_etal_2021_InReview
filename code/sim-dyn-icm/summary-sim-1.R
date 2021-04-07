# summary-sim-1.R: this script contains code to summarize the results
#                  of the first simulation scenario and produce 
#                  resulting figures. 
rm(list = ls())
library(coda)
library(nimble)
library(tidyverse)
library(ggthemes)
library(ggpubr)
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
params.mean <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.sd <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.high <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
params.highest <- matrix(NA, nrow = n.scenarios, ncol = n.params.tracked)
for (i in 1:n.scenarios) {
  print(i)
  # Read in results
  curr.file <- paste('results/sim-icm-results-4-100-condition-', i, 
		     '-simulations-2021-03-19.R', sep = '')
  load(curr.file)
  # Compute mean for each parameter and simulation 
  curr.mean <- sapply(samples.list, FUN = function(a) {apply(a, 2, mean)})
  # Get 25% quantile of means
  params.low[i, ] <- apply(curr.mean, 1, quantile, prob = 0.25)
  # Get median of means
  params.med[i, ] <- apply(curr.mean, 1, median)
  # Get mean of means
  params.mean[i, ] <- apply(curr.mean, 1, mean)
  # Compute standard deviation of means
  params.sd[i, ] <- apply(curr.mean, 1, sd)
  # Copmpute 75% quantile of means
  params.high[i, ] <- apply(curr.mean, 1, quantile, prob = 0.75)
}
# Assign parameter names to the column names
colnames(params.low) <- rownames(curr.mean)
colnames(params.med) <- rownames(curr.mean)
colnames(params.mean) <- rownames(curr.mean)
colnames(params.sd) <- rownames(curr.mean)
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
# Mean 
params.mean <- data.frame(params.mean)
params.mean$scenario <- 1:nrow(params.mean)
params.mean <- params.mean %>%
  select(-contains('mean')) %>%
  pivot_longer(cols = -scenario, values_to = 'mean', names_to = 'param')
params.mean$sp <- str_split(params.mean$param, '\\.') %>%
  map_chr(~ .[4])
params.mean$sp <- paste('Species ', params.mean$sp, sep = '')
params.mean$param <- str_split(params.mean$param, '\\.') %>%
  map_chr(~ paste0(.[1:3], collapse = '-'))
params.mean <- params.mean %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10')))
# Standard deviation
params.sd <- data.frame(params.sd)
params.sd$scenario <- 1:nrow(params.sd)
params.sd <- params.sd %>%
  select(-contains('mean')) %>%
  pivot_longer(cols = -scenario, values_to = 'sd', names_to = 'param')
params.sd$sp <- str_split(params.sd$param, '\\.') %>%
  map_chr(~ .[4])
params.sd$sp <- paste('Species ', params.sd$sp, sep = '')
params.sd$param <- str_split(params.sd$param, '\\.') %>%
  map_chr(~ paste0(.[1:3], collapse = '-'))
params.sd <- params.sd %>%
  mutate(param = factor(param), 
	 sp = factor(sp, levels = c('Species 1', 'Species 2', 'Species 3', 
				    'Species 4', 'Species 5', 'Species 6', 
				    'Species 7', 'Species 8', 
				    'Species 9', 'Species 10')))
# Combine all summary statistics together 
params.full <- params.low
params.full$med <- params.med$med
params.full$high <- params.high$high
params.full$mean <- params.mean$mean
params.full$sd <- params.sd$sd
# Assign info to params.full DF on what scenario each value corresponds to
params.full <- params.full %>%
  mutate(neon.in = param.vals$neon.in[scenario], 
	 bbs.in = param.vals$bbs.in[scenario], 
	 det.in = param.vals$det.in[scenario], 
	 eb.in = param.vals$eb.in[scenario])

# Take averages of the individual species averages
params.grouped <- params.full %>%
  group_by(scenario, param, neon.in, bbs.in, det.in, eb.in) %>%
  summarize(low = mean(low), 
	    med = mean(med), 
	    mean = mean(mean),
	    sd = mean(sd),
	    high = mean(high)) %>%
  ungroup()


# Some summary plots ------------------------------------------------------
# Rectangles for background colors
rects <- data.frame(xstart = c(0, 4.5, 10.5, 14.5), 
		    xend = c(4.5, 10.5, 14.5, 16), 
		    col = letters[1:4])

# Sorted order
# (1) eBird only, (2) HB only, (3) NEON only, (4) BBS only, (5) eb + HB, 
# (6) eBird + NEON, (7) HB + NEON, (8) eBird + BBS, (9) HB + BBS, (10) NEON + BBS, 
# (11) eBird + HB + NEON, (12) eBird + HB + BBS, (13) eBird + NEON + BBS, 
# (14) HB + NEON + BBS, (15) eBird + HB + NEON + BBS
# Intercept
# To create plot for a different parameter, change the variable name in filter
plotting.dat <- params.grouped %>%
  filter(param == 'beta-gamma-0') %>% 
  mutate(scenario = as.numeric(factor(scenario, levels = order(apply(param.vals, 1, sum)))))

# Get data in proper format for plot
plotting.dat <- plotting.dat %>%
  mutate(iqr.lowest = med - (high - low) * 1.5, # lower integer quartile range
	 iqr.highest = med + (high - low) * 1.5,  # upper inter quartile range
	 # Creates variable based on how many data sets are in a given scenario
	 group = ifelse(scenario < 5, 1, ifelse(scenario < 10, 2, 
						ifelse(scenario < 15, 3, 4))))
# Create plot for intercept
int.plot <- ggplot(data = plotting.dat) +
  geom_rect(data = rects, aes(xmin = xstart, 
  			      xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  	    alpha = 0.2) + 
  scale_fill_viridis_d() + 
  geom_boxplot(data = plotting.dat, aes(x = scenario, group = scenario, ymin = iqr.lowest, 
					lower = low, middle = med, upper = high, 
					ymax = iqr.highest), stat = 'identity', size = 00.9, fill = 'grey') + 
  theme_bw(base_size = 21) + 
  geom_vline(xintercept = 4.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 10.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 14.5, linetype = 2, col = 'grey') + 
  geom_hline(yintercept = 0, linetype = 2) + 
  guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = 'Data Sets Used', 
       y = 'Bias in species-level colonization intercept') + 
  scale_x_continuous(expand = c(0, 0), 
		     breaks = 1:15, labels = c('EBIRD', 'DET', 'NEON', 'BBS', 
					       'EBIRD + DET', 'EBIRD + NEON', 'DET + NEON', 
					       'EBIRD + BBS', 'DET + BBS', 'NEON + BBS', 
					       'EBIRD + DET + NEON', 'EBIRD + DET + BBS', 
					       'EBIRD + NEON + BBS', 'DET + NEON + BBS', 
					       'EBIRD + DET + NEON + BBS')) 

# Create plot for covariate 
plotting.dat <- params.grouped %>%
  filter(param == 'beta-gamma-1') %>% 
  mutate(scenario = as.numeric(factor(scenario, levels = order(apply(param.vals, 1, sum)))))

plotting.dat <- plotting.dat %>%
  mutate(iqr.lowest = med - (high - low) * 1.5,
	 iqr.highest = med + (high - low) * 1.5, 
	 group = ifelse(scenario < 5, 1, ifelse(scenario < 10, 2, 
						ifelse(scenario < 15, 3, 4))))
cov.plot <- ggplot(data = plotting.dat) +
   geom_rect(data = rects, aes(xmin = xstart, 
  			      xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  	    alpha = 0.2) + 
  scale_fill_viridis_d() + 
  geom_boxplot(data = plotting.dat, aes(x = scenario, group = scenario, ymin = iqr.lowest, 
					lower = low, middle = med, upper = high, 
					ymax = iqr.highest), stat = 'identity', size = 00.9, fill = 'grey') + 
  theme_bw(base_size = 21) + 
  geom_vline(xintercept = 4.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 10.5, linetype = 2, col = 'grey') + 
  geom_vline(xintercept = 14.5, linetype = 2, col = 'grey') + 
  geom_hline(yintercept = 0, linetype = 2) + 
  guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = 'Data Sets Used', 
       y = 'Bias in species-level colonization covariate effect') + 
  scale_x_continuous(expand = c(0, 0), 
		     breaks = 1:15, labels = c('EBIRD', 'DET', 'NEON', 'BBS', 
					       'EBIRD + DET', 'EBIRD + NEON', 'DET + NEON', 
					       'EBIRD + BBS', 'DET + BBS', 'NEON + BBS', 
					       'EBIRD + DET + NEON', 'EBIRD + DET + BBS', 
					       'EBIRD + NEON + BBS', 'DET + NEON + BBS', 
					       'EBIRD + DET + NEON + BBS')) 

# Put plot together and save to hard drive. 
ggarrange(int.plot, cov.plot, labels = c('A', 'B'), 
	  font.label = list(size = 22)) + 
  ggsave(device = 'pdf', filename = 'figures/sim-4-colonization.pdf', 
	 width = 18, height = 10, units = 'in')

# Summary statistics used in tables ---------------------------------------
# Change 
# summary statistics. 
params.grouped %>%
  mutate(scenario = as.numeric(factor(scenario, levels = order(apply(param.vals, 1, sum)))), 
	 param = factor(param, levels = c('beta-psi-0', 'beta-psi-1', 
					  'beta-phi-0', 'beta-phi-1', 
					  'beta-gamma-0', 'beta-gamma-1'), 
			ordered = TRUE)) %>%
  select(scenario, param, mean, sd) %>%
  arrange(scenario, param) %>%
  print(n = nrow(.))

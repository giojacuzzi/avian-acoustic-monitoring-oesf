####################################################################################
# A multi-species static occupancy model
#
# INPUT:
path_community_survey_data = "data/cache/1_derive_community_array/community_survey_data.rds"
path_environmental_data = "data/environment/PAM_PreHarvest_Habitat_results_DD_WD_TM.xlsx"
path_unit_key = "data/unit_key.csv"
####################################################################################

library(ggplot2)
library(tidyverse)
library(jagsUI)
library(MCMCvis)

# Season `t`
t = "2020"
threshold = 0.95 # naive conversion from continuous score prediction to binary detection-nondetection

message("Loading community survey data")
community_survey_data = readRDS(path_community_survey_data)
community_survey_data = community_survey_data[, , t, ]

# Derive detection-nondetection matricies
message("Deriving detection-nondetection matrix for each species")
species = dimnames(community_survey_data)[["common_name"]]
ylist   = setNames(vector("list", length(species)), species)
for (sp in species) {
  species_data = community_survey_data[, , sp] # Get the data for the species
  
  observation_matrix = matrix(
    unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence >= threshold, na.rm = TRUE)) else NA)),
    nrow = dim(community_survey_data)[1],
    ncol = dim(community_survey_data)[2],
    dimnames = dimnames(community_survey_data)[1:2]
  )
  ylist[[sp]] = observation_matrix
}

# Discard units with no survey observations and surveys with no unit observations
unit_survey_counts = lapply(ylist, function(x) { rowSums(!is.na(x))})
surveys_per_unit = as.data.frame(t(do.call(rbind, unit_survey_counts)))
units_not_surveyed = rownames(surveys_per_unit)[rowSums(surveys_per_unit) == 0]
if (length(units_not_surveyed) > 0) {
  message("Discarding ", length(units_not_surveyed), " units with no survey observations")
  ylist = lapply(ylist, function(mat) {
    mat[!(rownames(mat) %in% units_not_surveyed), , drop = FALSE]
  })
}
survey_unit_counts = lapply(ylist, function(x) { colSums(!is.na(x))})
units_per_survey = as.data.frame(t(do.call(rbind, survey_unit_counts)))
surveys_not_conducted = rownames(units_per_survey)[rowSums(units_per_survey) == 0]
if (length(surveys_not_conducted) > 0) {
  message("Discarding ", length(surveys_not_conducted), " surveys with no unit observations")
  ylist = lapply(ylist, function(mat) {
    mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE]
  })
}

units = rownames(surveys_per_unit)[rowSums(surveys_per_unit) != 0]
surveys = rownames(units_per_survey)[rowSums(units_per_survey) != 0]

lapply(ylist, head) # Inspect the detection history data

message("Frequency distribution of number of surveys conducted per unit:")
print(table(apply(!is.na(ylist[[1]]), 1, sum)))

message("Generating summary plots for detection history data")

# Plot detected number of species per unit
species_per_unit = setNames(rep(0, length(units)), units)
for (sp in ylist) {
  detected = rowSums(sp, na.rm = TRUE) > 0
  species_per_unit[detected] = species_per_unit[detected] + 1
}
species_per_unit = data.frame(
  unit = names(species_per_unit),
  species_detected = as.vector(species_per_unit)
)
mean_species_detected = mean(species_per_unit$species_detected)
p = ggplot(species_per_unit, aes(x = species_detected)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = mean_species_detected, color = "blue") +
  labs(title = "Naive species richness per unit", x = "Number of species detected", y = "Number of units") +
  theme_minimal(); print(p)

# Plot species detection frequency distribution
units_detected = sapply(ylist, function(mat) {
  sum(apply(mat, 1, function(x) any(x == 1, na.rm = TRUE)))
})
units_detected = data.frame(species = names(units_detected), units = units_detected) %>%
  arrange(desc(units)) %>% mutate(species = factor(species, levels = rev(species)))
p = ggplot(units_detected, aes(x = units, y = species)) +
  geom_bar(stat = "identity") +
  labs(title = "Species detections across sampling units", x = "Number of units with a detection", y = "") +
  theme_minimal(); print(p)

# TODO: Exclude species that were never detected
message("The following species were never detected and are excluded from the model:")
species_never_detected = units_detected %>% filter(units == 0) %>% pull(species)
print(species_never_detected)
ylist[as.character(species_never_detected)] <- NULL
species = names(ylist)

# Get environmental data
env_data = read.csv(path_unit_key)
env_data = env_data %>% filter(unit %in% units) %>% select(unit, stratum) %>% distinct()
strata = c('STAND INIT', 'COMP EXCL', 'THINNED', 'MATURE')
x_stratum = env_data %>%
  mutate(unit = factor(unit, levels = units), stratum = factor(stratum, levels = strata)) %>%
  arrange(unit) %>% select(unit, stratum) # arrange to match y

##### Following https://jamesepaterson.github.io/jamespatersonblog/2024-12-08_multispeciesoccupancymodels.html

# Restructure detection history data to summarize how many surveys a species was detected at
y = matrix(data = NA, nrow = length(units), ncol = length(species))
rownames(y) = units
colnames(y) = species
for (i in 1:length(units)) {
  for (sp in 1:length(ylist)) {
    y[i,sp] = sum(ylist[[sp]][i,], na.rm = TRUE)
  }
}

# Convert detection history to numeric
z = (y > 0)*1  # *1 converts to numeric, keeps matrix
z[z == 0] = NA

# Pack up data for JAGS
msom_jags_data = list(y = y, # y = the detection history matrix
                       nSobs = ncol(y), # nSobs = number of species observed
                       nSites = nrow(y),  # nSites = number of sites
                       nOcc = length(surveys), # nOcc = number of surveys
                       z = z) # z =  input detection data

# Look at structure
str(msom_jags_data)

# Hierarchical model
# For fitting basic multi-species occupancy models, we can think of the components as
# 1. Ecological model
# 2. Observation model
# 3. Priors that describe species level variation in psi and p
# 4. Hyperpriors that describe parameters for the community 

# Create and save file name "msom_simple.jags"
# Identify filepath of model file
msom_simple <- tempfile()

#Write model to file
writeLines("
model{
  for(k in 1:nSobs){  # Loop through species
    # Likelihood
    for(i in 1:nSites) { # Loop through sites
      # Ecological model
      z[i, k] ~ dbern(psi[k])
      # Observation model
      y[i, k] ~ dbin(p[k] * z[i, k], nOcc)
    }
    
    # Priors for Psi (species level)
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    psi[k] <- ilogit(lpsi[k])
    
    # Priors for p (species level)
    lp[k] ~ dnorm(mu.lp, tau.lp)
    p[k] <- ilogit(lp[k])
  }
  
  # Hyperpriors for Psi (community level)
  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  # Hyperpriors for p (community level)
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
}
", con = msom_simple)

# Vector of monitored parameters
wanted <- c("psi.mean", "p.mean", "psi", "p")
# Could add other parameters in model: "mu.lpsi", "sd.lpsi", "tau.lpsi",
# "mu.lp", "sd.lp", "tau.lp"

# model run takes <3 min on a macbook
msom_simple_out <- jags(msom_jags_data, NULL, wanted, msom_simple,
                        n.chains = 3, n.adapt = 100, n.iter = 10000, 
                        n.burnin = 5000, n.thin = 2, parallel = TRUE, DIC = FALSE)

## Diagnostics, checking chains for mixing and convergence with trace and density plots
# https://m-clark.github.io/bayesian-basics/diagnostics.html#monitoring-convergence

# We want Rhat values that are very close to 1.0. This is a test statistic for testing if the variance within chains is different than the variance between chains. It is meant to test if each chain was sampling from similar distributions. 
(summary = MCMCsummary(msom_simple_out))
message("The following parameters (and corresponding species) may not have converged:")
params_poor_rhat = rownames(summary[summary$Rhat > 1.01, ])
species_poor_rhat = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", params_poor_rhat))
print(summary[params_poor_rhat, ] %>% mutate(species = species[species_poor_rhat]))

# Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value.
MCMCtrace(msom_simple_out$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)

# Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak.
MCMCtrace(msom_simple_out$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)

# Gelman-Rubin convergence diagnostics (Potential Scale Reduction Factor). Values substantially above 1 indicate lack of convergence. If the chains have not converged, Bayesian credible intervals based on the t-distribution are too wide, and have the potential to shrink by this factor if the MCMC run is continued.
coda::gelman.diag(msom_simple_out$samples)

# Summary of results (using head() to just show higher levels)
# Full summary includes psi and p estimates for each species
head(summary(msom_simple_out))

# Put psi estimates in a dataframe
psi_species_estimates <- summary(msom_simple_out) %>%
  as.data.frame(.) %>%
  mutate(parameter = row.names(.)) %>%
  # Filtering to only the estimates of psi for each species 
  # species psi parameters in our model take the form "psi[1]"
  filter(parameter %in% c(paste0("psi","[",c(1:length(species)),"]"))) %>%
  arrange(mean) %>%
  mutate(species = species[as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter))])

# Plot psi estimates
p = ggplot(psi_species_estimates, aes(x = factor(species, levels = species), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  geom_hline(yintercept = summary(msom_simple_out)["psi.mean", "mean"], color = "blue") +
  geom_hline(yintercept = summary(msom_simple_out)["psi.mean", "2.5%"], color = "blue", linetype = "dashed") +
  geom_hline(yintercept = summary(msom_simple_out)["psi.mean", "97.5%"], color = "blue", linetype = "dashed") +
  labs(title = "Occupancy probability", x = "Species", y = "Occupancy estimate") +
  scale_x_discrete(labels = psi_species_estimates$species)+ 
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  coord_flip() +
  theme_classic(); print(p)

# Put p estimates in a dataframe
p_species_estimates <- summary(msom_simple_out) %>%
  as.data.frame(.) %>%
  mutate(parameter = row.names(.)) %>%
  # Filtering to only the estimates of psi for each species 
  # species psi parameters in our model take the form "psi[1]"
  filter(parameter %in% c(paste0("p","[",c(1:length(species)),"]"))) %>%
  arrange(mean) %>%
  mutate(species = species[as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter))])

# Plot p estimates
p = ggplot(p_species_estimates, aes(x = factor(species, levels = species), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  geom_hline(yintercept = summary(msom_simple_out)["p.mean", "mean"], color = "blue") +
  geom_hline(yintercept = summary(msom_simple_out)["p.mean", "2.5%"], color = "blue", linetype = "dashed") +
  geom_hline(yintercept = summary(msom_simple_out)["p.mean", "97.5%"], color = "blue", linetype = "dashed") +
  labs(title = "Detection probability", x = "Species", y = "Detection probability estimate") +
  scale_x_discrete(labels = p_species_estimates$species)+ 
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  coord_flip() +
  theme_classic(); print(p)




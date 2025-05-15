####################################################################################
# A multi-species static occupancy model with covariates
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

# Derive detection-nondetection and covariate matricies
message("Deriving detection-nondetection and covariate matricies for each species")
species = dimnames(community_survey_data)[["common_name"]]
ylist   = setNames(vector("list", length(species)), species)
xlist_yday = setNames(vector("list", length(species)), species)
for (sp in species) {
  species_data = community_survey_data[, , sp] # get the data for the species
  
  observation_matrix = matrix( # detection-nondetection
    unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence >= threshold, na.rm = TRUE)) else NA)),
    nrow = dim(community_survey_data)[1],
    ncol = dim(community_survey_data)[2],
    dimnames = dimnames(community_survey_data)[1:2]
  )
  ylist[[sp]] = observation_matrix
  
  covariate_yday_matrix = matrix( # survey date
    unlist(lapply(species_data, function(x) if (!is.null(x)) yday(x$survey_date) else NA)),
    nrow = dim(community_survey_data)[1],
    ncol = dim(community_survey_data)[2],
    dimnames = dimnames(community_survey_data)[1:2]
  )
  xlist_yday[[sp]] = covariate_yday_matrix
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
  xlist_yday = lapply(xlist_yday, function(mat) {
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
  xlist_yday = lapply(xlist_yday, function(mat) {
    mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE]
  })
}

units = rownames(surveys_per_unit)[rowSums(surveys_per_unit) != 0]
surveys = rownames(units_per_survey)[rowSums(units_per_survey) != 0]

# Inspect the detection history and covariate data
lapply(ylist, head)
lapply(xlist_yday, head)

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

# Exclude species that were never detected
message("The following species were never detected and are excluded from the model:")
species_never_detected = units_detected %>% filter(units == 0) %>% pull(species)
print(species_never_detected)
ylist[as.character(species_never_detected)]      <- NULL
xlist_yday[as.character(species_never_detected)] <- NULL
species = names(ylist)

# Format y for modeling as a 3D array (site × survey × species)
y = array(NA,
          dim = c(length(units), length(surveys), length(species)),
          dimnames = list(site = units, survey = surveys, species = species))
for (sp in seq_along(ylist)) {
  y[, , sp] = as.matrix(ylist[[sp]])
}

x_yday = array(NA,
               dim = c(length(units), length(surveys), length(species)),
               dimnames = list(site = units, survey = surveys, species = species))
for (sp in seq_along(xlist_yday)) {
  x_yday[, , sp] = as.matrix(xlist_yday[[sp]])
}

# Left-align data (moving any missing NA surveys to the right) to allow for direct indexing by number of surveys per site
left_align_row <- function(x) {
  non_na <- x[!is.na(x)]
  len <- length(x)
  c(non_na, rep(NA, len - length(non_na)))
}

y_new = y
x_yday_new = x_yday
for (sp in dimnames(y)[[3]]) {
  y_mat <- y[, , sp]
  y_aligned <- t(apply(y_mat, 1, left_align_row))
  dimnames(y_aligned) <- dimnames(y_mat)
  y_new[, , sp] <- y_aligned
  
  x_yday_mat <- x_yday[, , sp]
  x_yday_aligned <- t(apply(x_yday_mat, 1, left_align_row))
  dimnames(x_yday_aligned) <- dimnames(x_yday_mat)
  x_yday_new[, , sp] <- x_yday_aligned
}

# Number of surveys per site (excluding NAs)
nOcc <- apply(!is.na(y_new[, , 1]), 1, sum)

# # Get environmental data
# env_data = read.csv(path_unit_key)
# env_data = env_data %>% filter(unit %in% units) %>% select(unit, stratum) %>% distinct()
# strata = c('STAND INIT', 'COMP EXCL', 'THINNED', 'MATURE')
# x_stratum = env_data %>%
#   mutate(unit = factor(unit, levels = units), stratum = factor(stratum, levels = strata)) %>%
#   arrange(unit) %>% select(unit, stratum) # arrange to match y

y_flattened = matrix(data = NA, nrow = length(units), ncol = length(species))
rownames(y_flattened) = units
colnames(y_flattened) = species
for (i in 1:length(units)) {
  for (sp in 1:length(species)) {
    y_flattened[i,sp] = sum(y_new[i, , sp], na.rm = TRUE)
  }
}
z = (y_flattened > 0)*1
z[z == 0] = NA

# Pack up data for JAGS
msom_jags_data = list(y        = y_new,           # detection-nondetection matrix
                      x_yday   = x_yday_new,      # day of year detection covariate matrix
                      nSpecies = length(species), # number of species observed
                      nSites   = length(units),   # number of sites
                      surveys  = as.vector(nOcc), # number of surveys per site vector
                      z = z)
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
  for(k in 1:nSpecies) {  # loop through species
    
    # Likelihood
    for (i in 1:nSites) { # loop through sites
      # Ecological model
      z[i,k] ~ dbern(psi[k])
      
      for (j in 1:surveys[i]) { # loop through surveys at site i
      
        logit(p[i,j,k]) <- alpha[k] + beta[k] * x_yday[i,j,k] # logit of detection prob depends on survey-specific covariate
        
        # Observation model
        y[i,j,k] ~ dbern(p[i,j,k] * z[i,k])
      }
    }
    
    # Priors for Psi (species level)
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) # logit-scale occupancy probability for species k
    psi[k] <- ilogit(lpsi[k])          # occupancy probability for species k
    
    # Priors for p (species level)
    # lp[k] ~ dnorm(mu.lp, tau.lp) # logit-scale detection probability for species k
    # p[k] <- ilogit(lp[k])        # detection probability for species k
    
    # Priors for detection intercept and slope (species level)
    alpha[k] ~ dnorm(alpha.mean, tau.alpha) # baseline detectability of species k when the detection covariate
    beta[k] ~ dnorm(beta.mean, tau.beta)    # species-specific slope describing how detection changes with the covariate
  }
  
  # Hyperpriors for Psi (community level)
  psi.mean ~ dbeta(1, 1)     # community mean occupancy probability
  mu.lpsi <- logit(psi.mean) # mean of the species-specific logit-scale occupancy
  sd.lpsi ~ dunif(0, 5)      # standard deviation of logit-scale occpuancy among species
  tau.lpsi <- 1/sd.lpsi^2    # precision (inverse variance) of logit-scale occupancy
  
  # Hyperpriors for p (community level)
  # p.mean ~ dbeta(1, 1)       # community mean detection probability
  # mu.lp <- logit(p.mean)     # mean of species-specific logit-scale detection probability
  # sd.lp ~ dunif(0, 5)        # standard deviation of logit-scale detection among species
  # tau.lp <- 1/sd.lp^2        # precision of logit-scale detection
  
  # Hyperpriors for detection intercept
  alpha.mean ~ dnorm(0, 0.001) # community-level mean of species-specific intercepts (mean baseline detect prob across all species)
  sd.alpha ~ dunif(0, 5)       # community-level standard deviation of species-specific intercepts (how much detectability varies between species)
  tau.alpha <- 1 / pow(sd.alpha, 2)

  # Hyperpriors for detection slope (covariate effect)
  beta.mean ~ dnorm(0, 0.001) # community-level mean of the species-specific slopes describing how detection changes with x_yday
  sd.beta ~ dunif(0, 5)       # community-level standard deviation
  tau.beta <- 1 / pow(sd.beta, 2)
}
", con = msom_simple)

# Vector of monitored parameters
wanted <- c("psi.mean", "alpha.mean", "beta.mean", "psi", "alpha", "beta")

message("Running JAGS (current time ", time_start <- Sys.time(), ")")
msom_simple_out <- jags(msom_jags_data, NULL, wanted, msom_simple,
                        n.chains = 3, n.adapt = 100, n.iter = 10000, 
                        n.burnin = 5000, n.thin = 2, parallel = TRUE, DIC = FALSE)
message("Finished running JAGS (", round(as.numeric(Sys.time() - time_start) / 60, 2), " minutes)")

(summary = MCMCsummary(msom_simple_out))
MCMCtrace(msom_simple_out$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)
MCMCtrace(msom_simple_out$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)

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
  theme_minimal(); print(p)

## Predict mean detection probabilities across the observed range of yday values
yday_seq = seq(min(x_yday_new, na.rm = TRUE), max(x_yday_new, na.rm = TRUE), by = 1)

# For each species, calculate parameter posterior mean (and 95% credible interval) of detection probability
# by transforming each MCMC sample from log-odds to probability and averaging over the samples
alpha_samples = msom_simple_out$sims.list$alpha
beta_samples  = msom_simple_out$sims.list$beta
species_yday_posterior = data.frame()
for (k in 1:ncol(alpha_samples)) {
  probs = sapply(yday_seq, function(x) {
    plogis(alpha_samples[, k] + beta_samples[, k] * x)
  })
  probs_mean_ci = apply(probs, 2, function(x) {
    c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
  })
  species_yday_posterior = rbind(species_yday_posterior, data.frame(yday = yday_seq, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], species = as.character(k)
  ))
}

# For the community, calculate hyperparameter posterior mean (and 95% credible interval) of detection probability
alpha_mean_samples = msom_simple_out$sims.list$alpha.mean
beta_mean_samples  = msom_simple_out$sims.list$beta.mean
probs = sapply(yday_seq, function(x) {
  plogis(alpha_mean_samples + beta_mean_samples * x)
})
probs_mean_ci = apply(probs, 2, function(x) {
  c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
})
community_yday_posterior = data.frame(yday = yday_seq, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], species = "Community mean")

# Visualize relationship between yday and detection probability for the community
ggplot() +
  geom_line(data = species_yday_posterior, aes(x = yday, y = mean, group = species), alpha = 0.2) +
  geom_line(data = community_yday_posterior, aes(x = yday, y = mean), color = "blue", linewidth = 1.5) +
  geom_ribbon(data = community_yday_posterior, aes(x = yday, y = mean, ymin = lower, ymax = upper), alpha = 0.2, fill = NA, color = "blue", linetype = "dashed") +
  labs(x = "Day of year", y = "Detection probability", title = "Detected community") +
  theme_minimal()

# Visualize relationship between yday and detection probability for a specific species
species_name = "Orange-crowned Warbler"
species_idx = as.character(which(species == species_name))
ggplot(data = species_yday_posterior %>% filter(species == species_idx), aes(x = yday, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  lims(y = c(0.0, 1.0)) +
  labs(x = "Day of year", y = "Detection probability", title = paste0(species_name, " (", species_idx,")")) +
  theme_minimal()



####################################################################################
# A multi-species static occupancy model with species-specific detection covariates
#
# INPUT:
path_community_survey_data = "data/cache/1_derive_community_array/community_survey_data.rds"
path_environmental_data = "data/environment/PAM_PreHarvest_Habitat_results_DD_WD_TM.xlsx"
path_site_key = "data/unit_key.csv"
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

# Discard sites with no survey observations and surveys with no site observations
site_survey_counts = lapply(ylist, function(x) { rowSums(!is.na(x))})
surveys_per_site = as.data.frame(t(do.call(rbind, site_survey_counts)))
sites_not_surveyed = rownames(surveys_per_site)[rowSums(surveys_per_site) == 0]
if (length(sites_not_surveyed) > 0) {
  message("Discarding ", length(sites_not_surveyed), " sites with no survey observations")
  ylist = lapply(ylist, function(mat) {
    mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE]
  })
  xlist_yday = lapply(xlist_yday, function(mat) {
    mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE]
  })
}
survey_site_counts = lapply(ylist, function(x) { colSums(!is.na(x))})
sites_per_survey = as.data.frame(t(do.call(rbind, survey_site_counts)))
surveys_not_conducted = rownames(sites_per_survey)[rowSums(sites_per_survey) == 0]
if (length(surveys_not_conducted) > 0) {
  message("Discarding ", length(surveys_not_conducted), " surveys with no site observations")
  ylist = lapply(ylist, function(mat) {
    mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE]
  })
  xlist_yday = lapply(xlist_yday, function(mat) {
    mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE]
  })
}

sites   = rownames(surveys_per_site)[rowSums(surveys_per_site) != 0]
surveys = rownames(sites_per_survey)[rowSums(sites_per_survey) != 0]

# Inspect the detection history and covariate data
lapply(ylist, head)
lapply(xlist_yday, head)

message("Frequency distribution of number of surveys conducted per site:")
print(table(apply(!is.na(ylist[[1]]), 1, sum)))

message("Generating summary plots for detection history data")

# Plot detected number of species per site
species_per_site = setNames(rep(0, length(sites)), sites)
for (sp in ylist) {
  detected = rowSums(sp, na.rm = TRUE) > 0
  species_per_site[detected] = species_per_site[detected] + 1
}
species_per_site = data.frame(
  site = names(species_per_site),
  species_detected = as.vector(species_per_site)
)
mean_species_detected = mean(species_per_site$species_detected)
p = ggplot(species_per_site, aes(x = species_detected)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = mean_species_detected, color = "blue") +
  labs(title = "Naive species richness per site", x = "Number of species detected", y = "Number of sites") +
  theme_minimal(); print(p)

# Plot species detection frequency distribution
sites_detected = sapply(ylist, function(mat) {
  sum(apply(mat, 1, function(x) any(x == 1, na.rm = TRUE)))
})
sites_detected = data.frame(species = names(sites_detected), sites = sites_detected) %>%
  arrange(desc(sites)) %>% mutate(species = factor(species, levels = rev(species)))
p = ggplot(sites_detected, aes(x = sites, y = species)) +
  geom_bar(stat = "identity") +
  labs(title = "Species detections across sampling sites", x = "Number of sites with a detection", y = "") +
  theme_minimal(); print(p)

# Exclude species that were never detected
message("The following species were never detected and are excluded from the model:")
species_never_detected = sites_detected %>% filter(sites == 0) %>% pull(species)
print(species_never_detected)
ylist[as.character(species_never_detected)]      <- NULL
xlist_yday[as.character(species_never_detected)] <- NULL
species = names(ylist)

# Format observation detection-nondetection and covariate data for modeling as 3D arrays (site × survey × species)

# Observed detection-nondetection data
y = array(NA,
          dim = c(length(sites), length(surveys), length(species)),
          dimnames = list(site = sites, survey = surveys, species = species))
for (sp in seq_along(ylist)) {
  y[, , sp] = as.matrix(ylist[[sp]])
}

# Detection covariate data
# TODO: Standardize
x_yday = array(NA,
               dim = c(length(sites), length(surveys), length(species)),
               dimnames = list(site = sites, survey = surveys, species = species))
for (sp in seq_along(xlist_yday)) {
  x_yday[, , sp] = as.matrix(xlist_yday[[sp]])
}

# Left-align data (moving any missing NA surveys to the right) to allow for direct indexing by number of surveys per site
y_unaligned = y
x_yday_unaligned = x_yday
left_align_row = function(x) {
  non_na = x[!is.na(x)]
  c(non_na, rep(NA, length(x) - length(non_na)))
}
for (sp in dimnames(y)[[3]]) {
  sp_y_mat = y[, , sp]
  sp_y_aligned = t(apply(sp_y_mat, 1, left_align_row))
  dimnames(sp_y_aligned) = dimnames(sp_y_mat)
  y[, , sp] = sp_y_aligned
  
  sp_x_yday_mat = x_yday[, , sp]
  sp_x_yday_aligned = t(apply(sp_x_yday_mat, 1, left_align_row))
  dimnames(sp_x_yday_aligned) = dimnames(sp_x_yday_mat)
  x_yday[, , sp] = sp_x_yday_aligned
}
n_surveys_per_site <- apply(!is.na(y[, , 1]), 1, sum)

# # Get environmental data
# env_data = read.csv(path_site_key)
# env_data = env_data %>% filter(unit %in% sites) %>% select(unit, stratum) %>% distinct()
# strata = c('STAND INIT', 'COMP EXCL', 'THINNED', 'MATURE')
# x_stratum = env_data %>%
#   mutate(unit = factor(unit, levels = sites), stratum = factor(stratum, levels = strata)) %>%
#   arrange(unit) %>% select(unit, stratum) # arrange to match y

# Initialize latent occupancy state z[i] as 1 if a detection occurred at unit i, and 0 otherwise
z = matrix(data = NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (i in 1:length(sites)) {
  for (sp in 1:length(species)) { z[i,sp] = sum(y[i, , sp], na.rm = TRUE) }
}
z = (z > 0) * 1

msom_data = list(
  y = y,                                   # detection-nondetection matrix
  x_yday = x_yday,                         # day of year detection covariate matrix
  K = length(species),                     # number of species observed
  I = length(sites),                       # number of sites
  surveys  = as.vector(n_surveys_per_site) # number of surveys per site
)
str(msom_data)

# Specify hierarhical model and write to file
model_file = tempfile()
writeLines("
model{

  for(k in 1:K) { # for each species
    for (i in 1:I) { # for each site
      
      # Ecological process model for latent occurrence z
      z[i,k] ~ dbern(psi[k])
      
      for (j in 1:surveys[i]) { # for each survey at site i

        # Observation model for observed data y
        logit(p[i,j,k]) <- a0[k] + ayday[k] * x_yday[i,j,k] # logit of detection prob depends on survey-specific covariate
        
        y[i,j,k] ~ dbern(p[i,j,k] * z[i,k])
      }
    }
    
    # Priors for psi (species level)
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) # logit-scale occupancy probability for species k
    psi[k] <- ilogit(lpsi[k])          # occupancy probability for species k

    # Priors for detection intercept and slope (species level)
    a0[k] ~ dnorm(mu.a0, tau.a0) # baseline detectability of species k when the detection covariate
    ayday[k] ~ dnorm(mu.yday, tau.yday)    # species-specific slope describing how detection changes with the covariate
  }
  
  # Hyperpriors for psi (community level)
  psi.mean ~ dbeta(1, 1)     # community mean occupancy probability
  mu.lpsi <- logit(psi.mean) # mean of the species-specific logit-scale occupancy
  sd.lpsi ~ dunif(0, 5)      # standard deviation of logit-scale occpuancy among species
  tau.lpsi <- 1/sd.lpsi^2    # precision (inverse variance) of logit-scale occupancy
  
  # Hyperpriors for detection intercept (community level)
  mu.a0 ~ dnorm(0, 0.001)   # community mean of species-specific intercepts (mean baseline detect prob across all species)
  sd.a0 ~ dunif(0, 5)       # community standard deviation of species-specific intercepts (how much detectability varies between species)
  tau.a0 <- 1 / pow(sd.a0, 2)

  # Hyperpriors for yday covariate effect (community level)
  mu.yday ~ dnorm(0, 0.001) # community mean of the species-specific slopes describing how detection changes with x_yday
  sd.yday ~ dunif(0, 5)     # community standard deviation
  tau.yday <- 1 / pow(sd.yday, 2)
  
  # Derived quantities
  for (k in 1:K) {
    Nocc[k]  <- sum(z[,k]) # estimated number of occupied sites per species (among the sampled population of sites)
  }
  for (i in 1:I) {
    Nsite[i] <- sum(z[i,]) # estimated number of species occuring per site (among the species that were detected anywhere)
  }
}
", con = model_file)

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list(z = z) }, # initial values to avoid data/model conflicts
            parameters.to.save = c("psi.mean", "mu.a0", "mu.yday", "psi", "a0", "ayday", "Nsite", "Nocc"), # monitored parameters
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 10000, n.burnin = 5000, n.thin = 2,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

(msom_summary = summary(msom))
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)

# Put psi estimates in a dataframe
psi_species_estimates <- msom_summary %>%
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
  geom_hline(yintercept = msom_summary["psi.mean", "mean"], color = "blue") +
  geom_hline(yintercept = msom_summary["psi.mean", "2.5%"], color = "blue", linetype = "dashed") +
  geom_hline(yintercept = msom_summary["psi.mean", "97.5%"], color = "blue", linetype = "dashed") +
  labs(title = "Occupancy probability", x = "Species", y = "Occupancy estimate") +
  scale_x_discrete(labels = psi_species_estimates$species)+ 
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  coord_flip() +
  theme_minimal(); print(p)

## Predict mean detection probabilities across the observed range of yday values
yday_seq = seq(min(x_yday, na.rm = TRUE), max(x_yday, na.rm = TRUE), by = 1)

# For each species, calculate parameter posterior mean (and 95% credible interval) of detection probability
# by transforming each MCMC sample from log-odds to probability and averaging over the samples
a0_samples = msom$sims.list$a0
ayday_samples  = msom$sims.list$ayday
species_yday_posterior = data.frame()
for (k in 1:ncol(a0_samples)) {
  probs = sapply(yday_seq, function(x) {
    plogis(a0_samples[, k] + ayday_samples[, k] * x)
  })
  probs_mean_ci = apply(probs, 2, function(x) {
    c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
  })
  species_yday_posterior = rbind(species_yday_posterior, data.frame(yday = yday_seq, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], species = as.character(k)
  ))
}

# For the community, calculate hyperparameter posterior mean (and 95% CRI) of detection probability
a0_mean_samples    = msom$sims.list$mu.a0
yday_mean_samples  = msom$sims.list$mu.yday
probs = sapply(yday_seq, function(x) {
  plogis(a0_mean_samples + yday_mean_samples * x)
})
probs_mean_ci = apply(probs, 2, function(x) {
  c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
})
community_yday_posterior = data.frame(yday = yday_seq, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], species = "Community mean")

# Visualize relationship between yday and detection probability for the community
p = ggplot() +
  geom_line(data = species_yday_posterior, aes(x = yday, y = mean, group = species), alpha = 0.2) +
  geom_line(data = community_yday_posterior, aes(x = yday, y = mean), color = "blue", linewidth = 1.5) +
  geom_ribbon(data = community_yday_posterior, aes(x = yday, y = mean, ymin = lower, ymax = upper), fill = NA, color = "blue", linetype = "dashed") +
  labs(x = "Day of year", y = "Detection probability", title = "Detected community") +
  theme_minimal(); print(p)

# Visualize relationship between yday and detection probability for a specific species
species_name = "Orange-crowned Warbler"
species_idx = as.character(which(species == species_name))
p = ggplot(data = species_yday_posterior %>% filter(species == species_idx), aes(x = yday, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  lims(y = c(0.0, 1.0)) +
  labs(x = "Day of year", y = "Detection probability", title = paste0(species_name, " (", species_idx,")")) +
  theme_minimal(); print(p)

# Plot observed versus estimated number of species per site
Nsite_posterior = as.data.frame(msom_summary[grepl("^Nsite", rownames(msom_summary)), ])
Nsite_posterior$site = rownames(Nsite_posterior)
ggplot(Nsite_posterior, aes(x = site, y = mean)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  labs(x = "Site", y = "Estimated richness", title = "Estimated number of species per site") +
  theme_minimal()

####################################################################################
# A multi-species static occupancy model with occupancy and detection covariates
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

# Exclude species that were detected below a minimum number of sites
message("The following species were never detected and are excluded from the model:")
min_sites_detected = 5
species_to_remove = sites_detected %>% filter(sites < min_sites_detected) %>% pull(species)
print(species_to_remove)
ylist[as.character(species_to_remove)]      <- NULL
xlist_yday[as.character(species_to_remove)] <- NULL
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

# Get environmental data
env_data = read.csv(path_site_key)
env_data = env_data %>% filter(unit %in% sites) %>% select(unit, stratum) %>% distinct()
strata = c('STAND INIT', 'COMP EXCL', 'THINNED', 'MATURE')
x_stratum = env_data %>%
  mutate(site = factor(unit, levels = sites), stratum = factor(stratum, levels = strata)) %>%
  arrange(site) %>% select(site, stratum) # arrange to match y
all(x_stratum$site == sites) # check that sites are aligned between observation and covariate data
x_stratum = x_stratum %>% pull(stratum)

# Distribution of strata among sites
table(x_stratum)

# Initialize latent occupancy state z[i] as 1 if a detection occurred at unit i, and 0 otherwise
z = matrix(data = NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (i in 1:length(sites)) {
  for (sp in 1:length(species)) { z[i,sp] = sum(y[i, , sp], na.rm = TRUE) }
}
z = (z > 0) * 1

msom_data = list(
  # Observed data
  y = y,                                    # detection-nondetection matrix
  K = length(species),                      # number of species observed
  I = length(sites),                        # number of sites
  surveys  = as.vector(n_surveys_per_site), # number of surveys per site
  # Occupancy covariates
  x_stratum = as.numeric(x_stratum),        # forest stratum occupancy covariate site vector
  n_strata = length(strata),                # number of categories for stratum variable
  # Detection covariates
  x_yday = x_yday                           # day of year detection covariate site-survey matrix
)
str(msom_data)

# Specify hierarhical model and write to file
model_file = tempfile()
writeLines("
model{

  for(k in 1:K) { # for each species
    for (i in 1:I) { # for each site
      
      # Ecological process model for latent occurrence z
      logit(psi[i, k]) <- b0[k] + bstratum[k, x_stratum[i]]
      z[i,k] ~ dbern(psi[i, k])
      
      for (j in 1:surveys[i]) { # for each survey at site i

        # Observation model for observed data y
        logit(p[i,j,k]) <- a0[k] + ayday[k] * x_yday[i,j,k] # logit of detection prob depends on survey-specific covariate
        
        y[i,j,k] ~ dbern(p[i,j,k] * z[i,k])
      }
    }
    
    # Priors for occupancy slopes (species level)
    b0[k] ~ dnorm(mu.b0, tau.b0)                            # baseline occupancy of species k when the x_yday variable is zero
    for (s in 1:n_strata) {
      bstratum[k, s] ~ dnorm(mu.stratum[s], tau.stratum[s]) # species-specific effect of stratum on occupancy
    }

    # Priors for detection intercepts and slopes (species level)
    a0[k] ~ dnorm(mu.a0, tau.a0)           # baseline detectability of species k when the x_yday variable is zero
    ayday[k] ~ dnorm(mu.yday, tau.yday)    # species-specific effect of yday on detection
  }
  
  # Hyperpriors for occupancy intercept (community level)
  mu.b0 ~ dnorm(0, 0.01)
  sd.b0 ~ dunif(0, 5) # TODO: choose bounds of uniform by trial and error?
  tau.b0 <- pow(sd.b0, -2)
  
  # Hyperpriors for stratum covariate effect (community level)
  for (s in 1:n_strata) {
    mu.stratum[s] ~ dnorm(0, 0.001)       # community mean effect of stratum s
    sd.stratum[s] ~ dunif(0, 5)
    tau.stratum[s] <- pow(sd.stratum[s], -2)
  }
  
  # Hyperpriors for detection intercept (community level)
  mu.a0 ~ dnorm(0, 0.001)   # community mean of species-specific intercepts on the logit scale (mean baseline detect prob across all species)
  sd.a0 ~ dunif(0, 5)       # community standard deviation of species-specific intercepts (how much detectability varies between species)
  tau.a0 <- pow(sd.a0, -2)

  # Hyperpriors for yday covariate effect (community level)
  mu.yday ~ dnorm(0, 0.001) # community mean of the species-specific slopes on the logit scale describing how detection changes with x_yday
  sd.yday ~ dunif(0, 5)     # community standard deviation
  tau.yday <- pow(sd.yday, -2)
  
  # Derived quantities
  for (k in 1:K) {
    Nocc[k]  <- sum(z[ ,k]) # estimated number of occupied sites per species (among the sampled population of sites)
  }
  for (i in 1:I) {
    Nsite[i] <- sum(z[i, ]) # estimated number of species occuring per site (among the species that were detected anywhere)
  }
}
", con = model_file)

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list(z = z) }, # initial values to avoid data/model conflicts
            parameters.to.save = c("mu.b0", "mu.stratum", "b0", "bstratum", "mu.a0", "mu.yday", "a0", "ayday", "Nsite", "Nocc", "psi.mean"), # monitored parameters
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 10000, n.burnin = 5000, n.thin = 2,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

## Diagnostics, checking chains for mixing and convergence with trace and density plots
# https://m-clark.github.io/bayesian-basics/diagnostics.html#monitoring-convergence

# We want Rhat values that are very close to 1.0. This is a test statistic for testing if the variance within chains is different than the variance between chains. It is meant to test if each chain was sampling from similar distributions.
(msom_summary = summary(msom))
mcmc_summary = MCMCsummary(msom)
message("The following parameters may not have converged:")
params_poor_rhat = rownames(mcmc_summary[mcmc_summary$Rhat > 1.01, ])
# species_poor_rhat = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", params_poor_rhat))
# print(mcmc_summary[params_poor_rhat, ] %>% mutate(species = species[species_poor_rhat]))

# Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value.
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)

# Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak.
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)

# TODO: Gelman-Rubin convergence diagnostics (Potential Scale Reduction Factor). Values substantially above 1 indicate lack of convergence. If the chains have not converged, Bayesian credible intervals based on the t-distribution are too wide, and have the potential to shrink by this factor if the MCMC run is continued.
# coda::gelman.diag(msom$samples)

## Predict mean occupancy probabilities across the categories of forest strata

# For each species, calculate parameter posterior means (and 95% credible intervals) of occupancy probability
# by transforming each MCMC sample from log-odds to probability and averaging over the samples per stratum
b0_samples       = msom$sims.list$b0
bstratum_samples = msom$sims.list$bstratum
species_stratum_posterior = data.frame()
for (k in 1:ncol(b0_samples)) {
  for (l in 1:length(strata)) {
    probs = sapply(l, function(l) {
      plogis(b0_samples[, k] + bstratum_samples[, k, l])
    })
    probs_mean_ci = apply(probs, 2, function(x) {
      c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
    })
    species_stratum_posterior = rbind(species_stratum_posterior, data.frame(
      stratum = l, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], sp = species[k]
    ))
  }
}
species_stratum_posterior = species_stratum_posterior %>% arrange(mean) %>% mutate(sp = factor(sp, levels = species))

# For the community, calculate hyperparameter posterior mean (and 95% CRI) of detection probability
b0_mean_samples       = msom$sims.list$mu.b0
stratum_mean_samples  = msom$sims.list$mu.stratum
community_stratum_posterior = data.frame()
for (l in 1:length(strata)) {
  probs = sapply(l, function(x) {
    plogis(b0_mean_samples + stratum_mean_samples[, l])
  })
  probs_mean_ci = apply(probs, 2, function(x) {
    c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
  })
  community_stratum_posterior = rbind(community_stratum_posterior, data.frame(
    stratum = l, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], sp = "Community mean"
  ))
}

# Visualize relationship between stratum and occupancy probability for the community
for (stratum_name in strata) {
  stratum_idx  = which(strata == stratum_name)
  p = ggplot() +
    geom_point(data = species_stratum_posterior %>% filter(stratum == stratum_idx) %>% arrange(mean) %>% mutate(sp = factor(sp, levels = unique(sp))),
               aes(x = mean, y = sp)) +
    geom_errorbar(data = species_stratum_posterior %>% filter(stratum == stratum_idx) %>% arrange(mean) %>% mutate(sp = factor(sp, levels = unique(sp))),
                  aes(x = mean, y = sp, xmin = lower, xmax = upper), width = 0.2, alpha = 0.2) +
    geom_vline(xintercept = community_stratum_posterior %>% filter(stratum == stratum_idx) %>% pull(mean), color = "blue") +
    geom_vline(xintercept = community_stratum_posterior %>% filter(stratum == stratum_idx) %>% pull(lower), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = community_stratum_posterior %>% filter(stratum == stratum_idx) %>% pull(upper), color = "blue", linetype = "dashed") +
    xlim(0.0, 1.0) +
    labs(title = stratum_name, x = "Occupancy probability", y = "") +
    theme_classic(); print(p)
}

# Visualize relationship between stratum and occupancy probability for a specific species
species_name = "Brown Creeper"
species_idx = which(species == species_name)
p = ggplot(species_stratum_posterior %>% filter(sp == species_name) %>% mutate(stratum = factor(strata[stratum], levels = strata)),
           aes(x = stratum, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  ylim(0.0, 1.0) +
  labs(title = paste0(species_name, " (", species_idx,")"), x = "Stratum", y = "Occupancy probability") +
  theme_classic(); print(p)

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
  theme_classic(); print(p)

# Visualize relationship between yday and detection probability for a specific species
species_name = "Orange-crowned Warbler"
species_idx = as.character(which(species == species_name))
p = ggplot(data = species_yday_posterior %>% filter(species == species_idx), aes(x = yday, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  lims(y = c(0.0, 1.0)) +
  labs(x = "Day of year", y = "Detection probability", title = paste0(species_name, " (", species_idx,")")) +
  theme_classic(); print(p)

# Visualize estimated number of species per site
Nsite_posterior = as.data.frame(msom_summary[grepl("^Nsite", rownames(msom_summary)), ])
Nsite_posterior$site = rownames(Nsite_posterior)
ggplot(Nsite_posterior, aes(x = site, y = mean)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  labs(x = "Site", y = "Estimated richness", title = "Estimated number of species per site") +
  theme_classic()

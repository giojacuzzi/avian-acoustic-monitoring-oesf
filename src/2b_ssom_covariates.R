####################################################################################
# A single-season occupancy model with covariates
#
# INPUT:
path_community_array = "data/cache/1_derive_community_array/community_array.rds"
path_community_survey_data = "data/cache/1_derive_community_array/community_survey_data.rds"
path_environmental_data = "data/environment/PAM_PreHarvest_Habitat_results_DD_WD_TM.xlsx"
path_unit_key = "data/unit_key.csv"
####################################################################################

library(MCMCvis)
library(ggplot2)
library(jagsUI)
library(tidyverse)
library(readxl)
library(unmarked)
library(coda)

s = "Orange-crowned Warbler" # species s
t = "2020"          # season t
threshold = 0.95

community_survey_data = readRDS(path_community_survey_data)
y = community_survey_data[, , t, s] # Model observations for a single species from a single season
x_yday = matrix( # Survey day of year
  unlist(lapply(y, function(x) if (!is.null(x)) yday(x$survey_date) else NA)),
  nrow = dim(y)[1],
  ncol = dim(y)[2],
  dimnames = dimnames(y)[1:2]
)
y = matrix( # Threshold confidence scores to obtain binary presence-absence data
  unlist(lapply(y, function(x) if (!is.null(x)) as.integer(any(x$confidence >= threshold, na.rm = TRUE)) else NA)),
  nrow = dim(y)[1],
  ncol = dim(y)[2],
  dimnames = dimnames(y)[1:2]
)

# Clean matricies
unit_survey_counts = rowSums(!is.na(y))
if (sum(rowSums(!is.na(y)) == 0)) {
  warning("Discarding ", sum(rowSums(!is.na(y)) == 0), " units with no survey observations...")
  x_yday = x_yday[rowSums(!is.na(y)) > 0, ]
  y = y[rowSums(!is.na(y)) > 0, ]
}
survey_unit_counts = colSums(!is.na(y))
if (sum(survey_unit_counts == 0) > 0) {
  warning("Discarding ", sum(survey_unit_counts == 0), " surveys with no unit observations...")
  x_yday = x_yday[ , survey_unit_counts > 0]
  y = y[ , survey_unit_counts > 0]
}

nunit = nrow(y)   # unit i
nsurvey = ncol(y) # survey j

# Get environmental data
env_data = read_csv(path_unit_key)
env_data = env_data[env_data$season == t, ]
env_data = env_data[env_data$unit %in% rownames(y), ]
stratum = env_data %>% mutate(unit = factor(unit, levels = rownames(y))) %>% arrange(unit) %>% select(unit, stratum) # arrange to match y

# env_data_snags_deadwood = read_excel(path_environmental_data, sheet = 6, skip = 1) %>% rename(unit = 1, stratum = 2)
# env_data_snags_deadwood = env_data_snags_deadwood %>% mutate(unit = sapply(strsplit(as.character(unit), "_"), `[`, 1))
# env_data_snags_deadwood[,c('unit', 'vol_alldown_m3')]

# unmarked ---------------------------------------------
hab_levels = c('STAND INIT', 'COMP EXCL', 'THINNED', 'MATURE')
hab = factor(stratum$stratum, levels = hab_levels)

summary((
  y_unmarked = unmarkedFrameOccu(
    y = as.matrix(y),
    siteCovs = data.frame(hab = hab), # occupancy (site) covariate
    obsCovs  = list(yday = x_yday)    # detection (observation) covariate
  )
))

ssom_unmarked = occu(formula = 
                     ~yday   # detection probability varies by day of year
                     ~hab-1, # occupancy varies by habitat (no intercept, each habitat has a coefficient)
                     data = y_unmarked)
summary(ssom_unmarked)

# Predict the mean occupancy probability for each habitat
new_data = data.frame(hab = factor(c("STAND INIT", "COMP EXCL", "THINNED", "MATURE"),
                                    levels = c("STAND INIT", "COMP EXCL", "THINNED", "MATURE")))
pred_psi = predict(ssom_unmarked, type = "state", newdata = new_data)
pred_psi$habitat = new_data$hab
ggplot(pred_psi, aes(x = habitat, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  ylim(0.0, 1.0) +
  labs(title = paste(s, t, "(unmarked)"), x = "Habitat", y = "Occupancy probability") +
  theme_minimal()

# Predict the mean detection probability across the observed range of yday values
yday_seq = seq(min(x_yday, na.rm = T), max(x_yday, na.rm = T), length.out = 100)
pred_rho = predict(ssom_unmarked, type = "det", newdata = data.frame(yday = yday_seq))
pred_rho$yday = yday_seq
ggplot(pred_rho, aes(x = yday, y = Predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  xlim(min(x_yday, na.rm = T), max(x_yday, na.rm = T)) +
  ylim(0.0, 1.0) +
  labs(title = paste(s, t, "(unmarked)"), x = "Day of year", y = "Detection probability") +
  theme_minimal()

# jags -------------------------------------------------
niter = 50000
nburnin = 7500
nchains = 2
nthin = 1

# Jags requires all covariates to be completely defined (no NA values), so here we replace NA covariates with a dummy -1 value.
# Missing observations are excluded from the likelihood, so these dummy values will have no affect 
x_yday_clean = x_yday
x_yday_clean[is.na(x_yday_clean)] = 0

# We also should standardize our covariate
x_yday_clean = scale(x_yday_clean)

jags_file = tempfile(fileext = ".txt")
sink(jags_file)
cat("
model {
  # Priors
  alpha ~ dnorm(0, 0.001)  # Effect of day of year on detection probability
  for (h in 1:H) {
    beta[h] ~ dnorm(0, 1)  # Effect of habitat on occupancy
  }
  
  # Likelihood
  for (i in 1:I){ # For each site
    logit(psi[i]) <- beta[habitat[i]] # occupancy probability (no intercept)
    z[i] ~ dbern(psi[i])              # latent occupancy state
    
    for (j in 1:J) { # For each survey
      logit(p[i,j]) <- alpha * yday[i,j] # detection probability
      y[i,j] ~ dbern(z[i] * p[i,j])      # observed detection-nondetection
    }
  }
}
", fill = T)
sink()

init_func = function() {
  list( # Initial values to avoid data/model conflicts
    z = apply(y, 1, max, na.rm = T) # Initialize z[i] as 1 if a detection occurred at unit i, and 0 otherwise
  )
}

ssom_jags <- jags(
  data = list(
    # Observed data
    y = y,
    I = nunit,
    J = nsurvey,
    # Occupancy covariate
    H = length(unique(hab)),
    habitat = as.numeric(hab),
    # Detection covariate
    yday = x_yday_clean
  ),
  inits = init_func,
  parameters.to.save = c("beta", "alpha"),
  model.file = jags_file,
  n.chains = nchains,
  n.thin = nthin,
  n.iter = niter,
  n.burnin = nburnin
)

MCMCsummary(ssom_jags)

## Diagnostics
MCMCplot(ssom_jags, params = c('beta','alpha'),ci=c(50,95))

# Check for convergence with trace and density plots
MCMCtrace(ssom_jags$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)
MCMCtrace(ssom_jags$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)
gelman.diag(ssom_jags$samples) # Gelman-Rubin diagnostics

## Investigate posterior distributions for covariates

# Visualize posterior distributions for occupancy probability in relation to habitat covariates
beta_posterior_occupancy = as.data.frame(plogis(ssom_jags$sims.list$beta))
colnames(beta_posterior_occupancy) = hab_levels
ggplot(pivot_longer(beta_posterior_occupancy, cols = everything(), names_to = "Habitat", values_to = "Value"),
       aes(x = `Value`, fill = `Habitat`)) +
  geom_histogram(binwidth = 0.01, alpha = 0.6, position = "identity") +
  labs(title = "Posterior occupancy probability distributions by habitat", x = "Occupancy probability", y = "Frequency") +
  theme_minimal()

# Predict the mean occupancy probability for each habitat (beta coefficient)
occupancy_summary = data.frame(habitat = factor(hab_levels, levels = hab_levels))
# Calculate the mean probability and 95% credible intervals for each beta coefficient from the posterior
# by transforming each MCMC sample from log-odds to probability and averaging over the samples
beta_posterior_samples = ssom_jags$sims.list$beta
occupancy_summary$mean  = apply(beta_posterior_samples, 2, function(x) mean(plogis(x)))
occupancy_summary$lower = apply(beta_posterior_samples, 2, function(x) quantile(plogis(x), 0.025))
occupancy_summary$upper = apply(beta_posterior_samples, 2, function(x) quantile(plogis(x), 0.975))
ggplot(occupancy_summary, aes(x = habitat, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  ylim(0.0, 1.0) +
  labs(title = paste(s, t, "(jags)"), x = "Habitat", y = "Occupancy probability") +
  theme_minimal()

# Predict the mean detection probability across the observed range of yday values
yday_seq = seq(min(x_yday_clean, na.rm = T), max(x_yday_clean, na.rm = T), length.out = 100)
# For each yday value, calculate the mean detection probability and 95% credible intervals from the
# posterior by transforming each MCMC sample from log-odds to probability and averaging over the samples
alpha_posterior_samples = ssom_jags$sims.list$alpha
# detection_summary = t(sapply(yday_seq, function(yday_val) { # Retain the standardized scale
#   rho_samples = plogis(alpha_posterior_samples * yday_val)
#   c(yday = yday_val, mean = mean(rho_samples), quantile(rho_samples, 0.025), quantile(rho_samples, 0.975))
# }))
detection_summary = t(sapply(yday_seq, function(yday_val) {
  yday_original_scale = yday_val * sd(x_yday, na.rm = T) + mean(x_yday, na.rm = T) # Convert back to the original (unstandardized) scale
  rho_samples = plogis(alpha_posterior_samples * yday_val)
  c(yday = yday_original_scale, mean = mean(rho_samples), quantile(rho_samples, 0.025), quantile(rho_samples, 0.975))
}))
ggplot(detection_summary, aes(x = yday, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  # lims(x = c(min(x_yday_clean, na.rm = T), max(x_yday_clean, na.rm = T)), y = c(0.0, 1.0)) +
  lims(x = c(min(x_yday, na.rm = T), max(x_yday, na.rm = T)), y = c(0.0, 1.0)) +
  labs(title = paste(s, t, "(jags)"), x = "Day of year", y = "Detection probability") +
  theme_minimal()


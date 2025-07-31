####################################################################################
# A multi-species static occupancy model with occupancy and detection covariates
#
# INPUT:
path_community_survey_data = "data/cache/1_derive_community_array/community_survey_data.rds"
path_plot_scale_data = "data/cache/occurrence_covariates/data_plot_scale.rds"
####################################################################################

library(tidyverse)
library(jagsUI)
library(MCMCvis)
library(ggplot2)
theme_set(theme_classic())

# Get local plot scale data
local_plot_data = readRDS(path_plot_scale_data)
local_plot_data = local_plot_data %>% sf::st_drop_geometry() %>% arrange(site) %>% filter(hs == TRUE)

# Season `t`
t = "2020"
threshold = 0.95 # naive conversion from continuous score prediction to binary detection-nondetection

message("Loading community survey data")
community_survey_data = readRDS(path_community_survey_data)
community_survey_data = community_survey_data[, , t, ]
dimnames(community_survey_data)

# Derive putative observation and survey date matricies
message("Deriving detection-nondetection and yday matricies for each species")
species = dimnames(community_survey_data)[["common_name"]]
ylist   = setNames(vector("list", length(species)), species)
xlist_yday = setNames(vector("list", length(species)), species)
for (sp in species) {
  species_data = community_survey_data[, , sp]

  # Putative observation matrix (detection-nondetection via thresholded confidence scores)
  mat_obs = matrix(
    unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence >= threshold, na.rm = TRUE)) else NA)),
    nrow = dim(community_survey_data)[1],
    ncol = dim(community_survey_data)[2],
    dimnames = dimnames(community_survey_data)[1:2]
  )
  ylist[[sp]] = mat_obs
  
  # Survey date matrix (day of year)
  mat_yday = matrix(
    unlist(lapply(species_data, function(x) if (!is.null(x)) yday(x$survey_date) else NA)),
    nrow = dim(community_survey_data)[1],
    ncol = dim(community_survey_data)[2],
    dimnames = dimnames(community_survey_data)[1:2]
  )
  xlist_yday[[sp]] = mat_yday
}

# Discard sites with no environmental data
sites_missing_environmental_data = setdiff(dimnames(community_survey_data)$unit, local_plot_data$site)
if (length(sites_missing_environmental_data) > 0) {
  message("Discarding ", length(sites_missing_environmental_data), " sites with missing environmental data")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
}
# Discard sites with no observations
sites_with_environmental_data_missing_observations = setdiff(local_plot_data$site, dimnames(community_survey_data)$unit)
if (length(sites_with_environmental_data_missing_observations) > 0) {
  message("Discarding ", length(sites_with_environmental_data_missing_observations), " sites with missing observations")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_with_environmental_data_missing_observations), , drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[!(rownames(mat) %in% sites_with_environmental_data_missing_observations), , drop = FALSE] })
}

# Discard sites with no survey observations and surveys with no site observations
site_survey_counts = lapply(ylist, function(x) { rowSums(!is.na(x))})
surveys_per_site = as.data.frame(t(do.call(rbind, site_survey_counts)))
sites_not_surveyed = rownames(surveys_per_site)[rowSums(surveys_per_site) == 0]
if (length(sites_not_surveyed) > 0) {
  message("Discarding ", length(sites_not_surveyed), " sites with no survey observations")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE] })
}
survey_site_counts = lapply(ylist, function(x) { colSums(!is.na(x))})
sites_per_survey = as.data.frame(t(do.call(rbind, survey_site_counts)))
surveys_not_conducted = rownames(sites_per_survey)[rowSums(sites_per_survey) == 0]
if (length(surveys_not_conducted) > 0) {
  message("Discarding ", length(surveys_not_conducted), " surveys with no site observations")
  ylist = lapply(ylist, function(mat) { mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE] })
}

sites   = rownames(surveys_per_site)[rowSums(surveys_per_site) != 0]
surveys = rownames(sites_per_survey)[rowSums(sites_per_survey) != 0]

# Inspect the detection history and covariate data
message("Total number of sites: ", length(sites))
lapply(ylist, head)
lapply(xlist_yday, head)

n_surveys_per_site = apply(!is.na(ylist[[1]]), 1, sum)
message(sum(n_surveys_per_site), " total sampling periods (surveys) conducted across ", length(sites), " sampling units (sites)")
message("Sampling periods (surveys) conducted per site: median ", median(n_surveys_per_site), ", range ", min(n_surveys_per_site), "–", max(n_surveys_per_site))
print(table(n_surveys_per_site))

# Naive species occurrence
naive_occurrence = sapply(ylist, function(mat) { sum(apply(mat, 1, function(x) any(x == 1, na.rm = TRUE))) })
naive_occurrence = data.frame(species = names(naive_occurrence), nsites = naive_occurrence) %>%
  arrange(desc(nsites)) %>% mutate(species = factor(species, levels = rev(species))) %>% mutate(prob = nsites / length(sites))
message(naive_occurrence %>% filter(nsites > 0) %>% nrow(), " species detected")
message("Naive occurrence rate: mean ", round(mean(naive_occurrence$prob),2),
        ", min ", round(min(naive_occurrence$prob),2), ", max ", round(max(naive_occurrence$prob),2))
message("Most commonly detected species:")
cat(naive_occurrence %>% slice_max(prob, n=1) %>% pull(species) %>% as.character(), '\n')
message("Least commonly detected species:")
cat(naive_occurrence %>% slice_min(prob, n=1) %>% pull(species) %>% as.character(), '\n')
p = ggplot(naive_occurrence, aes(x = nsites, y = species)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = mean(naive_occurrence$nsites), color = "blue") +
  labs(title = "Species occurrence across sites", x = "Number of sites detected as occurrence", y = ""); print(p)

# Naive species richness per site
species_per_site = setNames(rep(0, length(sites)), sites)
for (sp in ylist) {
  detected = rowSums(sp, na.rm = TRUE) > 0
  species_per_site[detected] = species_per_site[detected] + 1
}
species_per_site = data.frame(site = names(species_per_site), species_detected = as.vector(species_per_site))
mean_species_detected = mean(species_per_site$species_detected)
message("Naive species richness per site: mean ", round(mean(species_per_site$species_detected),2), ", range ", min(species_per_site$species_detected), "–", max(species_per_site$species_detected))
p = ggplot(species_per_site, aes(x = species_detected)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = mean_species_detected, color = "blue") +
  labs(title = "Naive species richness per site", x = "Number of species detected", y = "Number of sites"); print(p)

# Exclude species that were detected below a minimum number of sites
message("The following species are excluded from the model:")
min_sites_detected = 5
species_to_remove = naive_occurrence %>% filter(nsites < min_sites_detected) %>% pull(species) %>% as.character() %>% sort()
species_to_remove = c(species_to_remove, "Great Horned Owl", "Hermit Warbler", "Pine Grosbeak", "Townsend's Solitaire") # TODO: incorporate these additional species determined to be absent via manual review above into the naive statistics
cat(species_to_remove)
ylist[species_to_remove]      = NULL
xlist_yday[species_to_remove] = NULL
species = names(ylist)

# Format observation detection-nondetection and covariate data for modeling as 3D arrays (site × survey × species)

# Observed detection-nondetection data
y = array(NA, dim = c(length(sites), length(surveys), length(species)),
          dimnames = list(site = sites, survey = surveys, species = species))
for (sp in seq_along(ylist)) {
  y[, , sp] = as.matrix(ylist[[sp]])
}

# Detection covariate data
x_yday = array(NA, dim = c(length(sites), length(surveys), length(species)),
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
n_surveys_per_site = apply(!is.na(y[, , 1]), 1, sum)

# Get scaled detection covariate data
x_yday_scaled = scale(as.vector(x_yday))

# Get scaled observational covariate data
local_plot_data = local_plot_data %>% filter(site %in% dimnames(y)$site) # discard data for irrelevant sites
all(dimnames(y)$site == local_plot_data$site) # check that covariate data are aligned with observation matrix by site
x_canopy_scaled = scale(local_plot_data$plot_canopy_cover_rs)

# Initialize latent occupancy state z[i] as 1 if a detection occurred at unit i, and 0 otherwise
z = matrix(data = NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (i in 1:length(sites)) {
  for (sp in 1:length(species)) { z[i,sp] = sum(y[i, , sp], na.rm = TRUE) }
}
z = (z > 0) * 1

# Prepare all data for the model
msom_data = list(
  # Observed data
  y = y,                                             # detection-nondetection matrix
  I = length(species),                               # number of species observed
  J = length(sites),                                 # number of sites sampled
  sampling_periods  = as.vector(n_surveys_per_site), # number of sampling periods (surveys) per site
  # Occupancy covariates
  x_canopy = x_canopy_scaled[,1],                    # canopy cover covariate site vector
  # Detection covariates
  x_yday = array(x_yday_scaled, dim = dim(x_yday)) # day of year detection covariate site-survey matrix
)
str(msom_data)

# Specify hierarhical model and write to file
# Following:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x
# - https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2293
model_file = tempfile()
writeLines("
model{

  for(i in 1:I) { # for each species
    for (j in 1:J) { # for each site
      
      # Ecological process model for latent occurrence z
      logit(psi[j, i]) <- u[i] + acanopy[i]*x_canopy[j]
      z[j,i] ~ dbern(psi[j, i])
      
      for (k in 1:sampling_periods[j]) { # for each sampling period (survey) at site j

        # Observation model for observed data y
        logit(p[j,k,i]) <- v[i] + byday[i]*x_yday[j,k,i] # logit of detection prob depends on survey-specific covariate
        
        y[j,k,i] ~ dbern(p[j,k,i] * z[j,i])
      }
    }
    
    # Priors for occupancy coefficients (species level)
    u[i] ~ dnorm(mu.u, tau.u)           # baseline occupancy of species i when the x_canopy variable is zero
    acanopy[i] ~ dnorm(mu.canopy, tau.canopy) # species-specific effect of canopy on occupancy

    # Priors for detection coefficients (species level)
    v[i] ~ dnorm(mu.v, tau.v)           # baseline detectability of species i when the x_yday variable is zero
    byday[i] ~ dnorm(mu.yday, tau.yday)    # species-specific effect of yday on detection
  }
  
  # Hyperpriors for occupancy intercept (community level)
  mu.u ~ dnorm(0, 0.01) # TODO: reasoning?
  sd.u ~ dunif(0, 5)    # TODO: choose bounds of uniform by trial and error?
  tau.u <- pow(sd.u,-2)
  
  # Hyperpriors for canopy covariate effect (community level)
  mu.canopy ~ dnorm(0, 0.001)
  sd.canopy ~ dunif(0, 5)
  tau.canopy <- pow(sd.canopy,-2)
  
  # Hyperpriors for detection intercept (community level)
  mu.v ~ dnorm(0, 0.001)   # community mean of species-specific intercepts on the logit scale (mean baseline detect prob across all species)
  sd.v ~ dunif(0, 5)       # community standard deviation of species-specific intercepts (how much detectability varies between species)
  tau.v <- pow(sd.v,-2)

  # Hyperpriors for yday covariate effect (community level)
  mu.yday ~ dnorm(0, 0.001) # community mean of the species-specific slopes on the logit scale describing how detection changes with x_yday
  sd.yday ~ dunif(0, 5)     # community standard deviation
  tau.yday <- pow(sd.yday,-2)
  
  # Derived quantities
  for (i in 1:I) {
    Nocc[i]  <- sum(z[ ,i]) # estimated number of occupied sites per species (among the sampled population of sites)
  }
  for (j in 1:J) {
    Nsite[j] <- sum(z[j, ]) # estimated number of species occuring per site (among the species that were detected anywhere)
  }
}
", con = model_file)

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list(z = z) }, # initial values to avoid data/model conflicts
            parameters.to.save = c( # monitored parameters
              "u", "mu.u", "sd.u",
              "v", "mu.v", "sd.v",
              "acanopy", "mu.canopy", "sd.canopy",
              "byday", "mu.yday", "sd.canopy",
              "Nsite", "Nocc"
            ),
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 10000, n.burnin = 5000, n.thin = 2,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

## Diagnostics, checking chains for mixing and convergence with trace and density plots
# https://m-clark.github.io/bayesian-basics/diagnostics.html#monitoring-convergence

# Gelman-Rubin statistic (Potential Scale Reduction Factor) Rhat values as a convergence diagnostic (Gelman and Rubin 1992). This is a test statistic for testing if the variance within chains is different than the variance between chains, and is meant to test if each chain was sampling from similar distributions -- if all the chains are “the same”, then the between chain variation should be close to zero. Rhat values substantially above 1.0 indicate lack of convergence; 1.2 is sometimes used as a guideline for “approximate convergence” (Brooks and Gelman 1998), but in practice a more stringent rule of Rhat < 1.1 is often used to declare convergence. If the chains have not converged, Bayesian credible intervals based on the t-distribution are too wide, and have the potential to shrink by this factor if the MCMC run is continued.
(msom_summary = summary(msom))
mcmc_summary = MCMCsummary(msom)
mcmc_summary %>% arrange(desc(Rhat))
message("The following parameters may not have converged:")
mcmc_summary %>% filter(Rhat > 1.0)

# Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value.
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)

# Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak.
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)

# If the model indicates convergence issues, we may need to:
# - Increase the burn-in period
# - Increase iterations
# - Use more informative priors
# - Reparametrize the model
msom_summary = summary(msom) %>% as.data.frame() %>% mutate(parameter = row.names(.)) %>%
  mutate(prob = plogis(mean), prob_lower95 = plogis(`2.5%`), prob_upper95 = plogis(`97.5%`))

# Visualize estimated number of species per site
Nsite_posterior = msom_summary %>% filter(stringr::str_starts(parameter, "Nsite")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(site_idx = as.integer(str_extract(parameter, "\\d+"))) %>% mutate(site = sites[site_idx])
Nsite_mean = mean(Nsite_posterior$mean)
message("Mean estimated species richness across all sites: ", round(Nsite_mean,1), " (range ", round(min(Nsite_posterior$mean),1), "–", round(max(Nsite_posterior$mean),1), ")")
ggplot(Nsite_posterior, aes(x = as.factor(plot_order), y = mean)) +
  geom_hline(yintercept = Nsite_mean, linetype = "solid", color = "blue") +
  geom_point() + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  scale_x_discrete(labels = Nsite_posterior$site) + 
  labs(title = "Estimated species richness per site", x = "Site", y = "Estimated species richness") +
  coord_flip()

# As covariates on occupancy were standardized, the inverse-logit (`plogis()`) of u[i] is the baseline occurrence probability for species i at a site with ‘average’ habitat characteristics.
occurrence_prob = msom_summary %>% filter(stringr::str_starts(parameter, "u")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(str_extract(parameter, "\\d+"))) %>% mutate(species_name = species[species_idx])
mean_occurrence_prob = msom_summary %>% filter(parameter == "mu.u")
message("Mean species occurrence probability: ", round(mean_occurrence_prob$prob,2), " (95% BCI ", round(mean_occurrence_prob$prob_lower95,2), "–", round(mean_occurrence_prob$prob_upper95,2), ")")
message("Species occurrence probability range: ", round(min(occurrence_prob$prob),2), "–", round(max(occurrence_prob$prob),2), " (", occurrence_prob %>% slice_min(prob, n=1) %>% pull(species_name), ", ", occurrence_prob %>% slice_max(prob, n=1) %>% pull(species_name), ")")
ggplot(occurrence_prob, aes(x = as.factor(plot_order), y = prob)) +
  geom_hline(yintercept = mean_occurrence_prob$prob, linetype = "solid", color = "blue") +
  geom_hline(yintercept = mean_occurrence_prob$prob_lower95, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mean_occurrence_prob$prob_upper95, linetype = "dashed", color = "blue") +
  geom_point() + geom_errorbar(aes(ymin = `prob_lower95`, ymax = `prob_upper95`), width = 0) +
  scale_x_discrete(labels = occurrence_prob$species_name) + 
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Baseline occurrence probability", x = "Species", y = "Occurrence probability") +
  coord_flip()

# Similarly, the inverse-logit of v[i] is the detection probability for species i under 'average' detection conditions.
detection_prob = msom_summary %>% filter(stringr::str_starts(parameter, "v")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(str_extract(parameter, "\\d+"))) %>% mutate(species_name = species[species_idx])
mean_detection_prob = msom_summary %>% filter(parameter == "mu.v")
message("Mean species detection probability: ", round(mean_detection_prob$prob,2), " (95% BCI ", round(mean_detection_prob$prob_lower95,2), "–", round(mean_detection_prob$prob_upper95,2), ")")
message("Species detection probability range: ", round(min(detection_prob$prob),2), "–", round(max(detection_prob$prob),2), " (", detection_prob %>% slice_min(prob, n=1) %>% pull(species_name), ", ", detection_prob %>% slice_max(prob, n=1) %>% pull(species_name), ")")
ggplot(detection_prob, aes(x = as.factor(plot_order), y = prob)) +
  geom_hline(yintercept = mean_detection_prob$prob, linetype = "solid", color = "blue") +
  geom_hline(yintercept = mean_detection_prob$prob_lower95, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mean_detection_prob$prob_upper95, linetype = "dashed", color = "blue") +
  geom_point() + geom_errorbar(aes(ymin = `prob_lower95`, ymax = `prob_upper95`), width = 0) +
  scale_x_discrete(labels = detection_prob$species_name) + 
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Baseline detection probability", x = "Species", y = "Detection probability") +
  coord_flip()

# Community-level summaries of the hyper-parameters for the detection and occupancy covariates
# mu is the community response (mean across species) to a given covariate and sd is the standard deviation (among species). Thus, the hyper-parameters are simply the mean and variance for each covariate as measured across species (Kéry & Royle 2009)
message("Community-level summaries of hyper-parameters for occurrence and detection covariates:")
msom_summary %>% filter(parameter %in% c(
  'mu.canopy', 'sd.canopy'
)) %>% select(mean, sd, `2.5%`, `97.5%`) %>% round(2)

# TODO: The coefficient acanopy is the effect of canopy cover on occupancy of species i.

# canopy_coefficients = summary(msom) %>% as.data.frame() %>% mutate(parameter = row.names(.)) %>%
#   filter(stringr::str_starts(parameter, "acanopy")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
#   mutate(species_idx = as.integer(str_extract(parameter, "\\d+"))) %>% mutate(species_name = species[species_idx]) %>%
#   mutate(bci_includes_zero = `2.5%` < 0.0 & `97.5%` > 0.0)
# 
# mu_canopy_summary = summary(msom) %>% 
#   as.data.frame() %>%
#   mutate(parameter = row.names(.)) %>%
#   filter(parameter == "mu.canopy") %>%
#   select(mean, `2.5%`, `97.5%`)
# 
# ggplot(canopy_coefficients, aes(x = as.factor(plot_order), y = mean)) +
#   geom_hline(yintercept = 0, color = "gray") +
#   geom_point(aes(color = bci_includes_zero)) +
#   geom_errorbar(aes(ymin = `25%`, ymax = `75%`, color = bci_includes_zero), width = 0, size = 1) +
#   geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = bci_includes_zero), width = 0) +
#   geom_hline(yintercept = mu_canopy_summary$mean, linetype = "solid", color = "blue") +
#   geom_hline(yintercept = mu_canopy_summary$`2.5%`, linetype = "dashed", color = "blue") +
#   geom_hline(yintercept = mu_canopy_summary$`97.5%`, linetype = "dashed", color = "blue") +
#   labs(x = "Species", y = "Canopy coefficient estimate on occupancy") +
#   scale_x_discrete(labels = canopy_coefficients$species_name) + 
#   # scale_y_continuous(limits = c(-0.25, 1.15), breaks = c(-0.25, 0, 0.25, 0.5, 1.0)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "gray")) +
#   coord_flip() +
#   theme(legend.position = "none") +
#   theme_classic()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## Predict mean occupancy probabilities across the categories of forest strata
# 
# # For each species, calculate parameter posterior means (and 95% credible intervals) of occupancy probability
# # by transforming each MCMC sample from log-odds to probability
# b0_samples      = msom$sims.list$b0
# bcanopy_samples = msom$sims.list$bcanopy
# species_posterior = data.frame()
# for (k in 1:ncol(b0_samples)) {
#     probs = plogis(b0_samples[, k] + bcanopy_samples[, k])
#     probs_mean_ci = apply(probs, 2, function(x) {
#       c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
#     })
#     species_posterior = rbind(species_posterior, data.frame(mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], sp = species[k]
#     ))
# }
# species_posterior = species_posterior %>% arrange(mean) %>% mutate(sp = factor(sp, levels = species))
# 
# # For the community, calculate hyperparameter posterior mean (and 95% CRI) of detection probability
# b0_mean_samples       = msom$sims.list$mu.b0
# stratum_mean_samples  = msom$sims.list$mu.stratum
# community_stratum_posterior = data.frame()
# for (l in 1:length(strata)) {
#   probs = sapply(l, function(x) {
#     plogis(b0_mean_samples + stratum_mean_samples[, l])
#   })
#   probs_mean_ci = apply(probs, 2, function(x) {
#     c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
#   })
#   community_stratum_posterior = rbind(community_stratum_posterior, data.frame(
#     stratum = l, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], sp = "Community mean"
#   ))
# }
# 
# # Visualize relationship between stratum and occupancy probability for the community
# for (stratum_name in strata) {
#   stratum_idx  = which(strata == stratum_name)
#   p = ggplot() +
#     geom_point(data = species_stratum_posterior %>% filter(stratum == stratum_idx) %>% arrange(mean) %>% mutate(sp = factor(sp, levels = unique(sp))),
#                aes(x = mean, y = sp)) +
#     geom_errorbar(data = species_stratum_posterior %>% filter(stratum == stratum_idx) %>% arrange(mean) %>% mutate(sp = factor(sp, levels = unique(sp))),
#                   aes(x = mean, y = sp, xmin = lower, xmax = upper), width = 0.2, alpha = 0.2) +
#     geom_vline(xintercept = community_stratum_posterior %>% filter(stratum == stratum_idx) %>% pull(mean), color = "blue") +
#     geom_vline(xintercept = community_stratum_posterior %>% filter(stratum == stratum_idx) %>% pull(lower), color = "blue", linetype = "dashed") +
#     geom_vline(xintercept = community_stratum_posterior %>% filter(stratum == stratum_idx) %>% pull(upper), color = "blue", linetype = "dashed") +
#     xlim(0.0, 1.0) +
#     labs(title = stratum_name, x = "Occupancy probability", y = "") +
#     theme_classic(); print(p)
# }
# 
# # Visualize relationship between stratum and occupancy probability for a specific species
# species_name = "Brown Creeper"
# species_idx = which(species == species_name)
# p = ggplot(species_stratum_posterior %>% filter(sp == species_name) %>% mutate(stratum = factor(strata[stratum], levels = strata)),
#            aes(x = stratum, y = mean)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
#   ylim(0.0, 1.0) +
#   labs(title = paste0(species_name, " (", species_idx,")"), x = "Stratum", y = "Occupancy probability") +
#   theme_classic(); print(p)
# 
# ## Predict mean detection probabilities across the observed range of yday values
# yday_seq = seq(min(x_yday, na.rm = TRUE), max(x_yday, na.rm = TRUE), by = 1)
# 
# # For each species, calculate parameter posterior mean (and 95% credible interval) of detection probability
# # by transforming each MCMC sample from log-odds to probability and averaging over the samples
# a0_samples = msom$sims.list$a0
# byday_samples  = msom$sims.list$byday
# species_yday_posterior = data.frame()
# for (k in 1:ncol(a0_samples)) {
#   probs = sapply(yday_seq, function(x) {
#     plogis(a0_samples[, k] + byday_samples[, k] * x)
#   })
#   probs_mean_ci = apply(probs, 2, function(x) {
#     c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
#   })
#   species_yday_posterior = rbind(species_yday_posterior, data.frame(yday = yday_seq, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], species = as.character(k)
#   ))
# }
# 
# # For the community, calculate hyperparameter posterior mean (and 95% CRI) of detection probability
# a0_mean_samples    = msom$sims.list$mu.a0
# yday_mean_samples  = msom$sims.list$mu.yday
# probs = sapply(yday_seq, function(x) {
#   plogis(a0_mean_samples + yday_mean_samples * x)
# })
# probs_mean_ci = apply(probs, 2, function(x) {
#   c(mean = mean(x), lower = quantile(x, 0.025, names = FALSE), upper = quantile(x, 0.975, names = FALSE))
# })
# community_yday_posterior = data.frame(yday = yday_seq, mean = probs_mean_ci["mean", ], lower = probs_mean_ci["lower", ], upper = probs_mean_ci["upper", ], species = "Community mean")
# 
# # Visualize relationship between yday and detection probability for the community
# p = ggplot() +
#   geom_line(data = species_yday_posterior, aes(x = yday, y = mean, group = species), alpha = 0.2) +
#   geom_line(data = community_yday_posterior, aes(x = yday, y = mean), color = "blue", linewidth = 1.5) +
#   geom_ribbon(data = community_yday_posterior, aes(x = yday, y = mean, ymin = lower, ymax = upper), fill = NA, color = "blue", linetype = "dashed") +
#   labs(x = "Day of year", y = "Detection probability", title = "Detected community") +
#   theme_classic(); print(p)
# 
# # Visualize relationship between yday and detection probability for a specific species
# species_name = "Orange-crowned Warbler"
# species_idx = as.character(which(species == species_name))
# p = ggplot(data = species_yday_posterior %>% filter(species == species_idx), aes(x = yday, y = mean)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
#   lims(y = c(0.0, 1.0)) +
#   labs(x = "Day of year", y = "Detection probability", title = paste0(species_name, " (", species_idx,")")) +
#   theme_classic(); print(p)
# 
# # Visualize estimated number of species per site
# Nsite_posterior = as.data.frame(msom_summary[grepl("^Nsite", rownames(msom_summary)), ])
# Nsite_posterior$site = rownames(Nsite_posterior)
# ggplot(Nsite_posterior, aes(x = site, y = mean)) +
#   geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
#   labs(x = "Site", y = "Estimated richness", title = "Estimated number of species per site") +
#   theme_classic()

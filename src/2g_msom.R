####################################################################################
# A multi-species static occupancy model with multiple occupancy and detection covariates
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
local_plot_data = local_plot_data %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>% filter(hs == TRUE)

# Season `t`
t = "2020"
year = as.numeric(t)
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

newxlist_yday = matrix(
  unlist(lapply(community_survey_data[, , 1], function(x) if (!is.null(x)) yday(x$survey_date) else NA)),
  nrow = dim(community_survey_data)[1],
  ncol = dim(community_survey_data)[2],
  dimnames = dimnames(community_survey_data)[1:2]
)



# Discard sites with no environmental data
sites_missing_environmental_data = setdiff(dimnames(community_survey_data)$unit, local_plot_data$site)
if (length(sites_missing_environmental_data) > 0) {
  message("Discarding ", length(sites_missing_environmental_data), " sites with missing environmental data")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
  
  newxlist_yday = newxlist_yday[!rownames(newxlist_yday) %in% sites_missing_environmental_data, ]
}
# Discard sites with no observations
sites_with_environmental_data_missing_observations = setdiff(local_plot_data$site, dimnames(community_survey_data)$unit)
if (length(sites_with_environmental_data_missing_observations) > 0) {
  message("Discarding ", length(sites_with_environmental_data_missing_observations), " sites with missing observations")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_with_environmental_data_missing_observations), , drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[!(rownames(mat) %in% sites_with_environmental_data_missing_observations), , drop = FALSE] })
  
  newxlist_yday = newxlist_yday[!rownames(newxlist_yday) %in% sites_with_environmental_data_missing_observations, ]
}

# Discard sites with no survey observations and surveys with no site observations
site_survey_counts = lapply(ylist, function(x) { rowSums(!is.na(x))})
surveys_per_site = as.data.frame(t(do.call(rbind, site_survey_counts)))
sites_not_surveyed = rownames(surveys_per_site)[rowSums(surveys_per_site) == 0]
if (length(sites_not_surveyed) > 0) {
  message("Discarding ", length(sites_not_surveyed), " sites with no survey observations")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE] })
  
  newxlist_yday = newxlist_yday[!rownames(newxlist_yday) %in% sites_not_surveyed, ]
}
survey_site_counts = lapply(ylist, function(x) { colSums(!is.na(x))})
sites_per_survey = as.data.frame(t(do.call(rbind, survey_site_counts)))
surveys_not_conducted = rownames(sites_per_survey)[rowSums(sites_per_survey) == 0]
if (length(surveys_not_conducted) > 0) {
  message("Discarding ", length(surveys_not_conducted), " surveys with no site observations")
  ylist = lapply(ylist, function(mat) { mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE] })
  xlist_yday = lapply(xlist_yday, function(mat) { mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE] })
  
  newxlist_yday = newxlist_yday[, !colnames(newxlist_yday) %in% surveys_not_conducted]
}

sites   = rownames(surveys_per_site)[rowSums(surveys_per_site) != 0]
surveys = rownames(sites_per_survey)[rowSums(sites_per_survey) != 0]

# Inspect the detection history and covariate data
message("Total number of sites: ", length(sites))
lapply(ylist, head)
lapply(xlist_yday, head)
head(newxlist_yday)

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

# newx_yday = array(NA, dim = c(length(sites), length(surveys)),
                  # dimnames = list(site = sites, survey = surveys))
# newx_yday = as.matrix(newxlist_yday)
newx_yday = newxlist_yday

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

newx_yday_unaligned = newx_yday
newx_yday = t(apply(newx_yday_unaligned, 1, left_align_row))
dimnames(newx_yday) = dimnames(newx_yday_unaligned)

# Get detection covariate data
detection_data = readRDS("data/cache/detection_covariates/data_detection.rds")
#############################################
# # Initialize empty array of same shape
# x_prcp = x_yday  # same dimensions and dimnames
# x_prcp[] = NA    # replace values with NA first
# for (i in seq_len(dim(x_yday)[1])) {        # sites
#   for (j in seq_len(dim(x_yday)[2])) {      # surveys
#     for (k in seq_len(dim(x_yday)[3])) {    # species
#       # message('site ', i,' survey ', j, ' species ', k)
#       # message('yday', x_yday[i, j, k])
#       site_i = dimnames(x_yday)$site[i]
#       yday_j = x_yday[i, j, k]
#       data = detection_data %>% filter(year == t, site == site_i, yday == yday_j)
# 
#       prcp = data %>% pull(tmax_deg_c)
#       
#       message('site ', site_i,' survey ', j, ' species ', k, ' yday ', yday_j, ' prcp ', prcp)
#       
#       x_prcp[i, j, k] = prcp
#       
#       # if (!is.null(val)) {
#       #   x_prcp[i, j, k] <- val
#       # }
#     }
#   }
# }

#############################################

# Standardize detection covariate data to z-scale (mean 0, standard deviation 1)
x_yday_scaled = scale(as.vector(x_yday))
newx_yday_scaled = scale(as.vector(newx_yday))

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
  x_yday = array(newx_yday_scaled, dim = dim(newx_yday)) # day of year detection covariate site-survey matrix
)
str(msom_data)

# Specify hierarhical model and write to file
# Following:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x
# - https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2293
model_file = tempfile()
writeLines("
model{

  ## Community level hyperpriors

  # Baseline occurrence (intercept)
  mu.u ~ dnorm(0, 0.01) # TODO: reasoning?
  sd.u ~ dunif(0, 5)    # TODO: choose bounds of uniform by trial and error?
  tau.u <- pow(sd.u,-2)
  
  # Covariate effects on occurrence
  mu.canopy ~ dnorm(0, 0.001)
  sd.canopy ~ dunif(0, 5)
  tau.canopy <- pow(sd.canopy,-2)
  
  # Baseline detection (intercept)
  mu.v ~ dnorm(0, 0.001)   # community mean of species-specific intercepts on the logit scale (mean baseline detect prob across all species)
  sd.v ~ dunif(0, 5)       # community standard deviation of species-specific intercepts (how much detectability varies between species)
  tau.v <- pow(sd.v,-2)

  # Covariate effects on detection
  mu.yday ~ dnorm(0, 0.001) # community mean of the species-specific slopes on the logit scale describing how detection changes with x_yday
  sd.yday ~ dunif(0, 5)     # community standard deviation
  tau.yday <- pow(sd.yday,-2)

  for(i in 1:I) { # for each species
    for (j in 1:J) { # for each site
      
      # Ecological process model for latent occurrence z
      logit(psi[j, i]) <- u[i] + acanopy[i]*x_canopy[j]
      z[j,i] ~ dbern(psi[j, i])
      
      for (k in 1:sampling_periods[j]) { # for each sampling period (survey) at site j

        # Observation model for observed data y
        logit(p[j,k,i]) <- v[i] + byday[i]*x_yday[j,k] # logit of detection prob depends on survey-specific covariate
        
        y[j,k,i] ~ dbern(p[j,k,i] * z[j,i])
      }
    }
    
    # Species level priors for occupancy coefficients
    u[i] ~ dnorm(mu.u, tau.u)                 # baseline occupancy of species i under 'average' conditions
    acanopy[i] ~ dnorm(mu.canopy, tau.canopy) # species-specific effect of canopy on occupancy

    # Species level priors for detection coefficients
    v[i] ~ dnorm(mu.v, tau.v)              # baseline detectability of species i under 'average' conditions
    byday[i] ~ dnorm(mu.yday, tau.yday)    # species-specific effect of yday on detection
  }
  
  ## Derived quantities
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
              "byday", "mu.yday", "sd.yday",
              "Nsite", "Nocc"
            ),
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 10000, n.burnin = 5000, n.thin = 2,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

## Diagnostics, checking chains for mixing and convergence with trace and density plots
# https://m-clark.github.io/bayesian-basics/diagnostics.html#monitoring-convergence

# Gelman-Rubin statistic (Potential Scale Reduction Factor) Rhat values as a convergence diagnostic (Gelman and Rubin 1992). This is a test statistic for testing if the variance within chains is different than the variance between chains, and is meant to test if each chain was sampling from similar distributions -- if all the chains are “the same”, then the between chain variation should be close to zero. Rhat values substantially above 1.0 indicate lack of convergence; 1.2 is sometimes used as a guideline for “approximate convergence” (Brooks and Gelman 1998), but in practice a more stringent rule of Rhat < 1.1 is often used to declare convergence. If the chains have not converged, Bayesian credible intervals based on the t-distribution are too wide, and have the potential to shrink by this factor if the MCMC run is continued.
msom_summary = summary(msom)
msom_summary = msom_summary %>% as_tibble() %>%
  mutate(param = rownames(msom_summary), overlap0 = as.factor(overlap0)) %>% relocate(param, .before = 1) %>%
  mutate(prob = plogis(mean), prob_lower95 = plogis(`2.5%`), prob_upper95 = plogis(`97.5%`))
message("The following parameters may not have converged:")
msom_summary %>% filter(Rhat >= 1.1)

# Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value.
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)

# Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak.
MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)

# If the model indicates convergence issues, we may need to:
# - Increase the burn-in period
# - Increase iterations
# - Use more informative priors
# - Reparametrize the model

# Visualize estimated number of species per site
Nsite_posterior = msom_summary %>% filter(stringr::str_starts(param, "Nsite")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(site_idx = as.integer(str_extract(param, "\\d+"))) %>% mutate(site = sites[site_idx])
Nsite_mean = mean(Nsite_posterior$mean)
message("Mean estimated species richness across all sites: ", round(Nsite_mean,1), " (range ", round(min(Nsite_posterior$mean),1), "–", round(max(Nsite_posterior$mean),1), ")")
ggplot(Nsite_posterior, aes(x = as.factor(plot_order), y = mean)) +
  geom_hline(yintercept = Nsite_mean, linetype = "solid", color = "blue") +
  geom_point() + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  scale_x_discrete(labels = Nsite_posterior$site) + 
  labs(title = "Estimated species richness per site", x = "Site", y = "Estimated species richness") +
  coord_flip()
ggplot(Nsite_posterior, aes(x = mean)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = Nsite_mean, color = "blue") +
  labs(title = "Estimated species richness per site", x = "Number of species estimated", y = "Number of sites")

# As covariates on occupancy were standardized, the inverse-logit (`plogis()`) of u[i] is the baseline occurrence probability for species i at a site with ‘average’ habitat characteristics.

occurrence_prob = msom_summary %>% filter(stringr::str_starts(param, "u")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+"))) %>% mutate(species_name = species[species_idx])
mean_occurrence_prob = msom_summary %>% filter(param == "mu.u")
message("Mean species occurrence probability: ", round(mean_occurrence_prob$prob,2), " (95% BCI ", round(mean_occurrence_prob$prob_lower95,2), "–", round(mean_occurrence_prob$prob_upper95,2), ")")
message("Species occurrence probability range: ", round(min(occurrence_prob$prob),2), "–", round(max(occurrence_prob$prob),2), " (", occurrence_prob %>% slice_min(prob, n=1) %>% pull(species_name), ", ", occurrence_prob %>% slice_max(prob, n=1) %>% pull(species_name), ")")
ggplot(occurrence_prob, aes(x = as.factor(plot_order), y = prob)) +
  geom_hline(yintercept = mean_occurrence_prob$prob,         linetype = "solid", color = "blue") +
  geom_hline(yintercept = mean_occurrence_prob$prob_lower95, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mean_occurrence_prob$prob_upper95, linetype = "dashed", color = "blue") +
  geom_point() + geom_errorbar(aes(ymin = `prob_lower95`, ymax = `prob_upper95`), width = 0) +
  scale_x_discrete(labels = occurrence_prob$species_name) + 
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Baseline occurrence probability", x = "Species", y = "Occurrence probability") +
  coord_flip()

# Similarly, the inverse-logit of v[i] is the detection probability for species i under 'average' detection conditions.
detection_prob = msom_summary %>% filter(stringr::str_starts(param, "v")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+"))) %>% mutate(species_name = species[species_idx])
mean_detection_prob = msom_summary %>% filter(param == "mu.v")
message("Mean species detection probability: ", round(mean_detection_prob$prob,2), " (95% BCI ", round(mean_detection_prob$prob_lower95,2), "–", round(mean_detection_prob$prob_upper95,2), ")")
message("Species detection probability range: ", round(min(detection_prob$prob),2), "–", round(max(detection_prob$prob),2), " (", detection_prob %>% slice_min(prob, n=1) %>% pull(species_name), ", ", detection_prob %>% slice_max(prob, n=1) %>% pull(species_name), ")")
ggplot(detection_prob, aes(x = as.factor(plot_order), y = prob)) +
  geom_hline(yintercept = mean_detection_prob$prob,         linetype = "solid", color = "blue") +
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
msom_summary %>% filter(param %in% c(
  'mu.canopy', 'sd.canopy'
)) %>% select(param, mean, sd, `2.5%`, `97.5%`)

# The coefficient acanopy is the effect of canopy cover on occurrence of species i
canopy_coef = msom_summary %>% filter(stringr::str_starts(param, "acanopy")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+"))) %>% mutate(species_name = species[species_idx])
mu_canopy_summary = msom_summary %>% filter(param == "mu.canopy") %>% select(mean, `2.5%`, `97.5%`)
ggplot(canopy_coef, aes(x = as.factor(plot_order), y = mean)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, size = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = overlap0), width = 0) +
  geom_hline(yintercept = mu_canopy_summary$mean,    linetype = "solid",  color = "blue") +
  geom_hline(yintercept = mu_canopy_summary$`2.5%`,  linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mu_canopy_summary$`97.5%`, linetype = "dashed", color = "blue") +
  labs(title = "Effect of canopy cover on occurrence", x = "Species", y = "Canopy coefficient estimate") +
  scale_x_discrete(labels = canopy_coef$species_name) +
  scale_color_manual(values = c("black", "gray")) +
  coord_flip() +
  theme(legend.position = "none")

# Mean marginal probabilities of occurrence for the metacommunity in relation to canopy cover
mean_canopy = attributes(x_canopy_scaled)$`scaled:center` # to transform between z-scale and original scale
sd_canopy   = attributes(x_canopy_scaled)$`scaled:scale`
original_canopy = seq(0, 100, by = 1) # range of possible canopy cover values
standardized_canopy = (original_canopy - mean_canopy) / sd_canopy
intercepts = msom_summary %>% # species-specific occurrence intercepts u[i]
  filter(str_starts(param, "u")) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+")), species_name = species[species_idx]) %>%
  select(species_name, u_i = mean)
intercepts_and_coeffs = canopy_coef %>% # species-specific canopy coefficients
  rename(acanopy_i = mean) %>%
  select(species_name, acanopy_i) %>%
  left_join(intercepts, by = "species_name")
canopy_preds = intercepts_and_coeffs %>% # predict species-specific occurrence probabilities
  rowwise() %>% do({
    i <- .
    tibble(
      species_name = i$species_name,
      canopy = original_canopy,
      psi = plogis(i$u_i + i$acanopy_i * standardized_canopy)
    )
  }) %>% bind_rows()
mu_u_samples      = as.matrix(msom$sims.list$mu.u) # predict meta-community occurrence probabilities (across posterior samples to calculate BCI)
mu_canopy_samples = as.matrix(msom$sims.list$mu.canopy)
meta_preds = map_dfr(1:nrow(mu_u_samples), function(i) {
  tibble(
    canopy = original_canopy,
    psi = plogis(mu_u_samples[i, 1] + mu_canopy_samples[i, 1] * standardized_canopy),
    draw = i
  )
})
meta_summary = meta_preds %>% # calculate means and 95% BCIs
  group_by(canopy) %>% summarise(
    psi_mean = mean(psi),
    psi_lower = quantile(psi, 0.025),
    psi_upper = quantile(psi, 0.975)
  )
ggplot(canopy_preds, aes(x = canopy, y = psi, group = species_name)) +
  geom_line(aes(color = species_name), alpha = 0.4) +
  geom_ribbon(data = meta_summary, aes(x = canopy, ymin = psi_lower, ymax = psi_upper), fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = canopy, y = psi_mean), color = "blue", size = 1.2, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Metacommunity occurrence probability in relation to canopy cover", x = "Canopy cover (%)", y = "Occurrence probability") +
  theme(legend.position = "none")

# Mean marginal probabilities of occurrence for specific species in relation to canopy cover
species_of_interest = c("Orange-crowned Warbler", "Brown Creeper", "Pacific Wren")
canopy_preds_filtered = canopy_preds %>% filter(species_name %in% species_of_interest)
species_idx = match(species_of_interest, species)
u_samples = as.matrix(msom$sims.list$u)[, species_idx] # posterior samples for intercept and covariate slopes for each species
acanopy_samples = as.matrix(msom$sims.list$acanopy)[, species_idx]
n_draws = nrow(u_samples) # posterior draws
n_canopy = length(original_canopy)
posterior_preds = map2_df( # Predict posterior occurrence across draws for each species
  .x = seq_along(species_idx), .y = species_of_interest,
  ~ {
    i <- .x
    name <- .y
    map_dfr(1:n_draws, function(draw) {
      tibble(
        species_name = name,
        canopy = original_canopy,
        psi = plogis(u_samples[draw, i] + acanopy_samples[draw, i] * standardized_canopy),
        draw = draw
      )
    })
  }
)
posterior_summary = posterior_preds %>%
  group_by(species_name, canopy) %>% summarise(
    psi_mean = mean(psi),
    psi_lower = quantile(psi, 0.025),
    psi_upper = quantile(psi, 0.975),
    .groups = "drop"
  )
ggplot(posterior_summary, aes(x = canopy, y = psi_mean, color = species_name, fill = species_name)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = psi_lower, ymax = psi_upper), alpha = 0.2, color = NA) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Occurrence probability in relation to canopy cover", x = "Canopy cover (%)", y = "Occurrence probability", color = "Species", fill = "Species")

# The coefficient byday is the effect of day of year on detection of species i
yday_coef = msom_summary %>% filter(stringr::str_starts(param, "byday")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+"))) %>% mutate(species_name = species[species_idx])
mu_yday_summary = msom_summary %>% filter(param == "mu.yday") %>% select(mean, `2.5%`, `97.5%`)
ggplot(yday_coef, aes(x = as.factor(plot_order), y = mean)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, size = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = overlap0), width = 0) +
  geom_hline(yintercept = mu_yday_summary$mean,    linetype = "solid",  color = "blue") +
  geom_hline(yintercept = mu_yday_summary$`2.5%`,  linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mu_yday_summary$`97.5%`, linetype = "dashed", color = "blue") +
  labs(title = "Effect of day of year on detection", x = "Species", y = "Day of year coefficient estimate") +
  scale_x_discrete(labels = yday_coef$species_name) +
  scale_color_manual(values = c("black", "gray")) +
  coord_flip() +
  theme(legend.position = "none")

# Mean marginal probabilities of detection for the metacommunity in relation to day of year
mean_yday = attributes(x_yday_scaled)$`scaled:center` # to transform between z-scale and original scale
sd_yday   = attributes(x_yday_scaled)$`scaled:scale`
original_yday = seq(min(x_yday, na.rm = TRUE), max(x_yday, na.rm = TRUE), by = 1) # range of observed yday values
standardized_yday = (original_yday - mean_yday) / sd_yday
intercepts = msom_summary %>% # species-specific detection intercepts v[i]
  filter(str_starts(param, "v")) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+")), species_name = species[species_idx]) %>%
  select(species_name, v_i = mean)
intercepts_and_coeffs = yday_coef %>% # species-specific yday coefficients
  rename(byday_i = mean) %>%
  select(species_name, byday_i) %>%
  left_join(intercepts, by = "species_name")
yday_preds = intercepts_and_coeffs %>% # predict species-specific occurrence probabilities
  rowwise() %>% do({
    i <- .
    tibble(
      species_name = i$species_name,
      yday = original_yday,
      psi = plogis(i$v_i + i$byday_i * standardized_yday)
    )
  }) %>% bind_rows()
mu_v_samples      = as.matrix(msom$sims.list$mu.v) # predict meta-community occurrence probabilities (across posterior samples to calculate BCI)
mu_yday_samples = as.matrix(msom$sims.list$mu.yday)
meta_preds = map_dfr(1:nrow(mu_v_samples), function(i) {
  tibble(
    yday = original_yday,
    psi = plogis(mu_v_samples[i, 1] + mu_yday_samples[i, 1] * standardized_yday),
    draw = i
  )
})
meta_summary = meta_preds %>% # calculate means and 95% BCIs
  group_by(yday) %>% summarise(
    psi_mean = mean(psi),
    psi_lower = quantile(psi, 0.025),
    psi_upper = quantile(psi, 0.975)
  )
ggplot(yday_preds, aes(x = yday, y = psi, group = species_name)) +
  geom_line(aes(color = species_name), alpha = 0.4) +
  geom_ribbon(data = meta_summary, aes(x = yday, ymin = psi_lower, ymax = psi_upper), fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = yday, y = psi_mean), color = "blue", size = 1.2, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(min(original_yday), max(original_yday))) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Metacommunity detection probability in relation to day of year", x = "Day of year", y = "Detection probability") +
  theme(legend.position = "none")

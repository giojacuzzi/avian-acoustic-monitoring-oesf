####################################################################################
# A multi-species static occupancy model with multiple occupancy and detection covariates
# Formulated according to:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x#b14
# - https://wildlife.onlinelibrary.wiley.com/doi/10.1002/jwmg.442
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
}

# Survey date matrix (day of year)
x_yday = matrix(
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
  
  x_yday = x_yday[!rownames(x_yday) %in% sites_missing_environmental_data, ]
}
# Discard sites with no observations
sites_with_environmental_data_missing_observations = setdiff(local_plot_data$site, dimnames(community_survey_data)$unit)
if (length(sites_with_environmental_data_missing_observations) > 0) {
  message("Discarding ", length(sites_with_environmental_data_missing_observations), " sites with missing observations")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_with_environmental_data_missing_observations), , drop = FALSE] })
  
  x_yday = x_yday[!rownames(x_yday) %in% sites_with_environmental_data_missing_observations, ]
}

# Discard sites with no survey observations and surveys with no site observations
site_survey_counts = lapply(ylist, function(x) { rowSums(!is.na(x))})
surveys_per_site = as.data.frame(t(do.call(rbind, site_survey_counts)))
sites_not_surveyed = rownames(surveys_per_site)[rowSums(surveys_per_site) == 0]
if (length(sites_not_surveyed) > 0) {
  message("Discarding ", length(sites_not_surveyed), " sites with no survey observations")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE] })
  
  x_yday = x_yday[!rownames(x_yday) %in% sites_not_surveyed, ]
}
survey_site_counts = lapply(ylist, function(x) { colSums(!is.na(x))})
sites_per_survey = as.data.frame(t(do.call(rbind, survey_site_counts)))
surveys_not_conducted = rownames(sites_per_survey)[rowSums(sites_per_survey) == 0]
if (length(surveys_not_conducted) > 0) {
  message("Discarding ", length(surveys_not_conducted), " surveys with no site observations")
  ylist = lapply(ylist, function(mat) { mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE] })
  
  x_yday = x_yday[, !colnames(x_yday) %in% surveys_not_conducted]
}

sites   = rownames(surveys_per_site)[rowSums(surveys_per_site) != 0]
surveys = rownames(sites_per_survey)[rowSums(sites_per_survey) != 0]

# Inspect the detection history and covariate data
message("Total number of sites: ", length(sites))
lapply(ylist, head)
head(x_yday)

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
species = names(ylist)

# Format observation detection-nondetection and covariate data for modeling as 3D arrays (site × survey × species)

# Observed detection-nondetection data
y = array(NA, dim = c(length(sites), length(surveys), length(species)),
          dimnames = list(site = sites, survey = surveys, species = species))
for (sp in seq_along(ylist)) {
  y[, , sp] = as.matrix(ylist[[sp]])
}

# Left-align data (moving any missing NA surveys to the right) to allow for direct indexing by number of surveys per site
y_unaligned = y
left_align_row = function(x) {
  non_na = x[!is.na(x)]
  c(non_na, rep(NA, length(x) - length(non_na)))
}
for (sp in dimnames(y)[[3]]) {
  sp_y_mat = y[, , sp]
  sp_y_aligned = t(apply(sp_y_mat, 1, left_align_row))
  dimnames(sp_y_aligned) = dimnames(sp_y_mat)
  y[, , sp] = sp_y_aligned
}
n_surveys_per_site = apply(!is.na(y[, , 1]), 1, sum)

x_yday_unaligned = x_yday
x_yday = t(apply(x_yday_unaligned, 1, left_align_row))
dimnames(x_yday) = dimnames(x_yday_unaligned)

# Get occurrence and detection covariate data
local_plot_data = local_plot_data %>% filter(site %in% dimnames(y)$site) # discard data for irrelevant sites
all(dimnames(y)$site == local_plot_data$site) # check that covariate data are aligned with observation matrix by site
detection_data = readRDS("data/cache/detection_covariates/data_detection.rds") %>% filter(year == t)
#############################################

x_yday_df = as.data.frame(x_yday)
x_yday_df$site = rownames(x_yday_df)
x_yday_long = tidyr::pivot_longer(
  x_yday_df,
  cols = -site,
  names_to = "survey",
  values_to = "yday"
)

get_var_matrix = function(variable) {
  detection_data_long = x_yday_long %>%
    left_join(detection_data, by = c("site", "yday")) %>%
    select(site, survey, !!sym(variable))
  x = tidyr::pivot_wider(
    detection_data_long,
    names_from = survey,
    values_from = !!sym(variable)
  )
  x = as.data.frame(x)
  rownames(x) = x$site
  x$site = NULL
  return(as.matrix(x))
}

x_tmax = get_var_matrix("tmax_deg_c")
x_prcp = get_var_matrix("prcp_mm_day")

# Standardize detection covariate data to z-scale (mean 0, standard deviation 1)
x_tmax_scaled = scale(as.vector(x_tmax))
x_prcp_scaled = scale(as.vector(x_prcp))
x_yday_scaled = scale(as.vector(x_yday))

# Standardize occurrence covariate data
x_alpha1_scaled = scale(local_plot_data$plot_canopy_cover_rs)
x_alpha2_scaled = scale(local_plot_data$plot_downvol_hs)
x_alpha3_scaled = scale(local_plot_data$plot_ht_cv_hs)

# Initialize latent occupancy state z[i] as 1 if a detection occurred at unit i, and 0 otherwise
z = matrix(data = NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (i in 1:length(sites)) {
  for (sp in 1:length(species)) { z[i,sp] = sum(y[i, , sp], na.rm = TRUE) }
}
z = (z > 0) * 1

# "The two-point finite mixture model proposed by Royle and Link (2006)... first formulated an occupancy model for false positives which applies to the standard occupancy design and data structure. The model has one additional parameter p10, which is the probability of a false detection, or Prðy ¼ 1jz ¼ 0Þ,wherez is the latent occupancy state. In this model, no distinction is made between true or false positives, and instead, the patterns of detections in encounter histories are used to determine the most likely values of the false-negative (p01) and false-positive error rates (p10)." (Kery and Royle 2020)

# Prepare all data for the model
msom_data = list(
  # Observed data
  y = y,                                             # detection-nondetection matrix
  I = length(species),                               # number of species observed
  J = length(sites),                                 # number of sites sampled
  K  = as.vector(n_surveys_per_site),                # number of sampling periods (surveys) per site
  # Occupancy covariates
  x_alpha1  = x_alpha1_scaled[,1],                   # alpha1 cover covariate site vector
  x_alpha2  = x_alpha2_scaled[,1],
  x_alpha3  = x_alpha3_scaled[,1],
  # Detection covariates
  # x_tmax = array(x_tmax_scaled, dim = dim(x_tmax)), # TODO: shared among species?
  # x_prcp = array(x_prcp_scaled, dim = dim(x_prcp)), # TODO: shared among species?
  x_beta1 = array(x_yday_scaled, dim = dim(x_yday)) # day of year detection covariate site-survey matrix
)
str(msom_data)

# Specify hierarhical model and write to file
# Following:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x
# - https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2293
# - https://www.sciencedirect.com/science/article/abs/pii/S0006320709004819
model_file = tempfile()
writeLines("
model{

  ## Community level hyperpriors

  # Baseline occurrence (intercept)
  mean.occ ~ dunif(0,1)
  mu.u <- log(mean.occ) - log(1-mean.occ) # `a`
  tau.u ~ dgamma(0.1,0.1) # `tau1`
  sigma.u <- 1/sqrt(tau.u)
  
  # Covariate effects on occurrence
  mu.alpha1 ~ dnorm(0, 0.001)
  mu.alpha2 ~ dnorm(0, 0.001)
  mu.alpha3 ~ dnorm(0, 0.001)
  tau.alpha1 ~ dgamma(0.1,0.1)
  tau.alpha2 ~ dgamma(0.1,0.1)
  tau.alpha3 ~ dgamma(0.1,0.1)
  
  # Baseline true-positive detection (intercept)
  mean.tp ~ dunif(0,1)
  mu.v <- log(mean.tp) - log(1-mean.tp) # `b` # community mean of species-specific intercepts on the logit scale (mean baseline detect prob across all species)
  tau.v ~ dgamma(0.1,0.1) # `tau2`
  rho ~ dunif(-1,1) # Because high abundance species are likely to be both easier to detect and more prevalent across the landscape, we modelled a correlation `rho` between occurrence and detection in the model
  var.v <- tau.v/(1.-pow(rho,2))
  sigma.v <- 1/sqrt(tau.v)

  # Covariate effects on true-positive detection
  mu.beta1  ~ dnorm(0, 0.001) # community mean of the species-specific slopes on the logit scale describing how detection changes with x_yday
  tau.beta1 ~ dgamma(0.1,0.1)
  
  # Baseline false-positive detection
  mean.fp ~ dunif(0,1)
  mu.fp <- log(mean.fp) - log(1-mean.fp)
  tau.fp ~ dgamma(0.1,0.1)
  

  for (i in 1:I) { # for each species
  
      # Species level priors for occupancy coefficients
      u[i] ~ dnorm(mu.u, tau.u)                # baseline occupancy of species i under 'average' conditions
      alpha1[i] ~ dnorm(mu.alpha1, tau.alpha1) # species-specific effect of alpha1 on occupancy
      alpha2[i] ~ dnorm(mu.alpha2, tau.alpha2)
      alpha3[i] ~ dnorm(mu.alpha3, tau.alpha3)
  
      # Species level priors for true-positive detection coefficients
      v[i] ~ dnorm(mu.v + (rho*sigma.v/sigma.u)*(u[i]-mu.u), var.v)              # baseline detectability of species i under 'average' conditions
      beta1[i] ~ dnorm(mu.beta1, tau.beta1)  # species-specific effect of yday on detection
      
      # Species level priors for false-positive detection
      fp[i] ~ dnorm(mu.fp, tau.fp)
  
      for (j in 1:J) { # for each site
        
          # Ecological process model for latent occurrence z
          logit(psi[j, i]) <- u[i] + alpha1[i]*x_alpha1[j] + alpha2[i]*x_alpha2[j] + alpha3[i]*x_alpha3[j]
          z[j,i] ~ dbern(psi[j, i])
          
          for (k in 1:K[j]) { # for each sampling period (survey) at site j
              
              # True positive detection probability
              logit(p11[j,k,i]) <- v[i] + beta1[i]*x_beta1[j,k] # logit of detection prob depends on survey-specific covariate
              
              # False positive detection probability
              logit(p10[j,k,i]) <- fp[i]  # species-specific constant false positive rate
              
              # Detection probability mixture
              p[j,k,i] <- z[j,i]*p11[j,k,i] + (1 - z[j,i])*p10[j,k,i]
              
              # Observation model for observed data y
              y[j,k,i]    ~ dbern(p[j,k,i])
              ynew[j,k,i] ~ dbern(p[j,k,i])
              
              # Create simulated dataset to calculate the Bayesian p-value
              d[j,k,i] <- abs(y[j,k,i] - p[j,k,i]) 
              dnew[j,k,i] <- abs(ynew[j,k,i] - p[j,k,i]) 
              d2[j,k,i] <- pow(d[j,k,i],2)  
              dnew2[j,k,i] <- pow(dnew[j,k,i],2)
          }
          dsum[j,i] <- sum(d2[j,1:K[j],i]) 
          dnewsum[j,i] <- sum(dnew2[j,1:K[j],i])
      }
  }
  
  ## Derived quantities
  p.fit <- sum(dsum[1:J,1:I]) # discrepancy measure, then defined as the mean(p.fit > p.fitnew)
  p.fitnew <- sum(dnewsum[1:J,1:I])
  for (i in 1:I) {
    Nocc[i] <- sum(z[ ,i]) # estimated number of occupied sites per species (among the sampled population of sites)
  }
  for (j in 1:J) {
    Nsite[j] <- sum(z[j, ]) # estimated number of species occuring per site (among the species that were detected anywhere)
  }
}
", con = model_file)

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

logit <- function(p) {
  log(p / (1 - p))
}

msom = jags(data = msom_data,
            inits = function() { list(
              # initial values to avoid data/model conflicts
              z = z,
              # enforce p11 > p10 for multimodal likelihood
              v = rep(logit(0.70), length(species)),
              fp = rep(logit(0.05), length(species))
            ) },
            parameters.to.save = c( # monitored parameters
              "u", "mu.u", "tau.u", # TODO: also calculate sigma (standard deviation) directly?
              "v", "mu.v", "tau.v",
              "fp", "mu.fp", "tau.fp",
              "alpha1", "mu.alpha1", "tau.alpha1",
              "alpha2", "mu.alpha2", "tau.alpha2",
              "alpha3", "mu.alpha3", "tau.alpha3",
              "beta1", "mu.beta1", "tau.beta1",
              "p.fit", "p.fitnew", "Nsite", "Nocc"
            ),
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 10000, n.burnin = 5000, n.thin = 2,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

## Diagnostics, checking chains for mixing and convergence with trace and density plots
# https://m-clark.github.io/bayesian-basics/diagnostics.html#monitoring-convergence

# Gelman-Rubin statistic (i.e. "potential scale reduction factor") Rhat values serve as a convergence diagnostic (Gelman and Rubin 1992). This is a test statistic for testing if the variance within chains is different than the variance between chains, and is meant to test if each chain was sampling from similar distributions -- if all the chains are “the same”, then the between chain variation should be close to zero. Rhat values substantially above 1.0 indicate lack of convergence; 1.2 is sometimes used as a guideline for “approximate convergence” (Brooks and Gelman 1998), but in practice a more stringent rule of Rhat < 1.1 is often used to declare convergence. If the chains have not converged, Bayesian credible intervals based on the t-distribution are too wide, and have the potential to shrink by this factor if the MCMC run is continued.
msom_summary = summary(msom)
msom_summary = msom_summary %>% as_tibble() %>%
  mutate(param = rownames(msom_summary), overlap0 = as.factor(overlap0)) %>% relocate(param, .before = 1) %>%
  mutate(prob = plogis(mean), prob_lower95 = plogis(`2.5%`), prob_upper95 = plogis(`97.5%`))
message("The following parameters may not have converged:")
msom_summary %>% filter(Rhat >= 1.1)

# Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value.
MCMCtrace(msom$samples, excl = c('Nsite', 'Nocc'), post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE, pdf = F)

# Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak.
MCMCtrace(msom$samples, excl = c('Nsite', 'Nocc'), post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE, pdf = F)

# If the model indicates convergence issues, we may need to:
# - Increase the burn-in period
# - Increase iterations
# - Use more informative priors
# - Reparametrize the model

# "We assessed the adequacy of the model using the approach suggested by Gelman et al. (1996) referred to as a Bayesian p-value. We defined a discrepancy measure for the observations yi and their expected values pi under the model. This discrepancy statistic is computed at each iteration of the MCMC algorithm. A reference distribution is computed by simulating data sets from the posterior distribution and computing the discrepancy measure, Dsim, for the simulated data set. The Bayesian p-value is defined as the probability: Pr(D > Dsim). Extreme values (e.g. less than 0.05 or greater than 0.95) indicate that the model is inadequate." (Zipkin et al. 2009) A value of 0.5 indicates good fit, while a value of 0.25 or 0.75 indicates moderate fit.
# "Our P-value, also known as a posterior predictive check, followed the same process as Carrillo-Rubio et al. (2014), Kroll et al. (2014), and Tobler et al. (2015), but was based on the deviance residuals rather than the Pearson's residuals because of its relationship to information criterion theory (Spiegelhalter et al. 1998). We describe the approach in detail in Appendix S1." (TODO: compare to Broms et al. 2016)
p_fit    = msom$sims.list$p.fit
p_fitnew = msom$sims.list$p.fitnew
p_value  = mean(p_fit > p_fitnew) # proportion of samples for which the model fits the observed data worse than the simulated data
message("Bayesian p-value: ", round(p_value, 3))

# Visualize estimated number of species per site
Nsite_posterior = msom_summary %>% filter(stringr::str_starts(param, "Nsite")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(site_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(site = sites[site_idx])
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
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
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
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
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
(occurrence_coeff_summary = msom_summary %>% filter(param %in% c(
  'mu.alpha1', 'tau.alpha1',
  'mu.alpha2', 'tau.alpha2',
  'mu.alpha3', 'tau.alpha3'
)) %>% select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0))
(detection_coeff_summary = msom_summary %>% filter(param %in% c(
  'mu.beta1', 'tau.beta1'
)) %>% select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0))

# Compare community level effect sizes for occurrence coefficients
ggplot(occurrence_coeff_summary %>% filter(str_starts(param, "mu")), aes(x = mean, y = as.factor(param))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
  scale_color_manual(values = c("black", "gray")) +
  labs(title = "Community level effect sizes for occurrence covariates", x = "Coefficient estimate", y = "Parameter")

# The coefficient alpha1 is the effect of alpha1 cover on occurrence of species i
alpha1_coef = msom_summary %>% filter(stringr::str_starts(param, "alpha1")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
mu_alpha1_summary = msom_summary %>% filter(param == "mu.alpha1") %>% select(mean, `2.5%`, `97.5%`)
ggplot(alpha1_coef, aes(x = as.factor(plot_order), y = mean)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = overlap0), width = 0) +
  geom_hline(yintercept = mu_alpha1_summary$mean,    linetype = "solid",  color = "blue") +
  geom_hline(yintercept = mu_alpha1_summary$`2.5%`,  linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mu_alpha1_summary$`97.5%`, linetype = "dashed", color = "blue") +
  labs(title = "Effect of alpha1 on occurrence", x = "Species", y = "alpha1 coefficient estimate") +
  scale_x_discrete(labels = alpha1_coef$species_name) +
  scale_color_manual(values = c("black", "gray")) +
  coord_flip() +
  theme(legend.position = "none")

# Mean marginal probabilities of occurrence for the metacommunity in relation to alpha1 cover
mean_alpha1 = attributes(x_alpha1_scaled)$`scaled:center` # to transform between z-scale and original scale
sd_alpha1   = attributes(x_alpha1_scaled)$`scaled:scale`
original_alpha1 = seq(0, 100, by = 1) # range of possible alpha1 cover values
standardized_alpha1 = (original_alpha1 - mean_alpha1) / sd_alpha1
intercepts = msom_summary %>% # species-specific occurrence intercepts u[i]
  filter(str_starts(param, "u")) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+")), species_name = species[species_idx]) %>%
  select(species_name, u_i = mean)
intercepts_and_coeffs = alpha1_coef %>% # species-specific alpha1 coefficients
  rename(alpha1_i = mean) %>%
  select(species_name, alpha1_i) %>%
  left_join(intercepts, by = "species_name")
alpha1_preds = intercepts_and_coeffs %>% # predict species-specific occurrence probabilities
  rowwise() %>% do({
    i <- .
    tibble(
      species_name = i$species_name,
      alpha1 = original_alpha1,
      psi = plogis(i$u_i + i$alpha1_i * standardized_alpha1)
    )
  }) %>% bind_rows()
mu_u_samples      = as.matrix(msom$sims.list$mu.u) # predict meta-community occurrence probabilities (across posterior samples to calculate BCI)
mu_alpha1_samples = as.matrix(msom$sims.list$mu.alpha1)
meta_preds = map_dfr(1:nrow(mu_u_samples), function(i) {
  tibble(
    alpha1 = original_alpha1,
    psi = plogis(mu_u_samples[i, 1] + mu_alpha1_samples[i, 1] * standardized_alpha1),
    draw = i
  )
})
meta_summary = meta_preds %>% # calculate means and 95% BCIs
  group_by(alpha1) %>% summarise(
    psi_mean = mean(psi),
    psi_lower = quantile(psi, 0.025),
    psi_upper = quantile(psi, 0.975)
  )
ggplot(alpha1_preds, aes(x = alpha1, y = psi, group = species_name)) +
  geom_line(aes(color = species_name), alpha = 0.4) +
  geom_ribbon(data = meta_summary, aes(x = alpha1, ymin = psi_lower, ymax = psi_upper), fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = alpha1, y = psi_mean), color = "blue", size = 1.2, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Metacommunity occurrence probability in relation to alpha1 cover", x = "alpha1 cover (%)", y = "Occurrence probability") +
  theme(legend.position = "none")

# Mean marginal probabilities of occurrence for specific species in relation to alpha1 cover
species_of_interest = c("Orange-crowned Warbler", "Brown Creeper", "Pacific Wren")
alpha1_preds_filtered = alpha1_preds %>% filter(species_name %in% species_of_interest)
species_idx = match(species_of_interest, species)
u_samples = as.matrix(msom$sims.list$u)[, species_idx] # posterior samples for intercept and covariate slopes for each species
alpha1_samples = as.matrix(msom$sims.list$alpha1)[, species_idx]
n_draws = nrow(u_samples) # posterior draws
n_alpha1 = length(original_alpha1)
posterior_preds = map2_df( # Predict posterior occurrence across draws for each species
  .x = seq_along(species_idx), .y = species_of_interest,
  ~ {
    i <- .x
    name <- .y
    map_dfr(1:n_draws, function(draw) {
      tibble(
        species_name = name,
        alpha1 = original_alpha1,
        psi = plogis(u_samples[draw, i] + alpha1_samples[draw, i] * standardized_alpha1),
        draw = draw
      )
    })
  }
)
posterior_summary = posterior_preds %>%
  group_by(species_name, alpha1) %>% summarise(
    psi_mean = mean(psi),
    psi_lower = quantile(psi, 0.025),
    psi_upper = quantile(psi, 0.975),
    .groups = "drop"
  )
ggplot(posterior_summary, aes(x = alpha1, y = psi_mean, color = species_name, fill = species_name)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = psi_lower, ymax = psi_upper), alpha = 0.2, color = NA) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Occurrence probability in relation to alpha1 cover", x = "alpha1 cover (%)", y = "Occurrence probability", color = "Species", fill = "Species")

# The coefficient beta1 is the effect of day of year on detection of species i
yday_coef = msom_summary %>% filter(stringr::str_starts(param, "beta1")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
mu_yday_summary = msom_summary %>% filter(param == "mu.beta1") %>% select(mean, `2.5%`, `97.5%`)
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
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param)), species_name = species[species_idx]) %>%
  select(species_name, v_i = mean)
intercepts_and_coeffs = yday_coef %>% # species-specific yday coefficients
  rename(beta1_i = mean) %>%
  select(species_name, beta1_i) %>%
  left_join(intercepts, by = "species_name")
yday_preds = intercepts_and_coeffs %>% # predict species-specific occurrence probabilities
  rowwise() %>% do({
    i <- .
    tibble(
      species_name = i$species_name,
      yday = original_yday,
      psi = plogis(i$v_i + i$beta1_i * standardized_yday)
    )
  }) %>% bind_rows()
mu_v_samples      = as.matrix(msom$sims.list$b) # predict meta-community occurrence probabilities (across posterior samples to calculate BCI)
mu_yday_samples = as.matrix(msom$sims.list$mu.beta1)
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

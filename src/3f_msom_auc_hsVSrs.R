####################################################################################
# Compare sets of predictor variables among a shared number of sites:
# - Local plot variables (measured by habitat survey)
# - Local plot variables (measured by remote sensing)
# - Homerange (i.e. neighborhood/landscape) variables
# QUESTION: Is community assembly at a given OESF forest location influenced more by local structural characteristics of the stand, or by the surrounding composition and configuration of the landscape?
#
# A multi-species static occupancy model with predictive performance measured via ROC AUC for model selection
#
# Formulated according to:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x#b14
# - https://wildlife.onlinelibrary.wiley.com/doi/10.1002/jwmg.442
#
# INPUT:
model_name = "RSHR" # "HS" habitat survey local or "RS" remote sensing local or "HR" homerange
# Naive conversion from continuous score prediction to binary detection-nondetection
naive_threshold = NA # set to NA to use probabilistic threshold from 1c_calculate_species_thresholds.R
generate_diagnostic_plots = FALSE
####################################################################################

path_community_survey_data = "data/cache/1_derive_community_survey_data/community_survey_data_2025-08-15.rds"
path_species_thresholds    = "data/cache/1_calculate_species_thresholds/species_thresholds_manual_selection.csv"
path_plot_scale_data       = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_homerange_scale_data  = "data/cache/occurrence_covariates/data_homerange_scale.rds"
path_out = paste0("data/cache/models/msom_", model_name, "_", format(Sys.Date(), "%Y-%m-%d"), ".rds")

library(progress)
library(car)
library(tidyverse)
library(jagsUI)
library(MCMCvis)
library(glue)
library(ggplot2)
library(ggrepel)
theme_set(theme_classic())

# Get occurrence data (local plot and homerange scales)
occ_data_plot_shared = readRDS(path_plot_scale_data) %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, hs, elev, dist_watercourse_major
)

occ_data_rs_shared = readRDS(path_homerange_scale_data)[['plot']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, homerange_canopy_cover_mean, homerange_canopy_cover_cv
)

occ_data_plot_hs = readRDS(path_plot_scale_data) %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, plot_downvol_hs, plot_ht_cv_hs, plot_ht_hs, plot_lcr_hs,
  plot_qmd_all_hs, plot_snagden_hs, plot_tree_all_diversity,
  plot_treeden_all_hs, plot_treeden_gt10cmDbh_hs, plot_treeden_lt10cmDbh_hs, plot_understory_vol
)

occ_data_plot_rs = readRDS(path_homerange_scale_data)[['plot']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, homerange_downvol_mean, homerange_htmax_cv, homerange_qmd_mean, homerange_treeden_all_mean, homerange_treeden_gt4in_dbh_mean, homerange_ba_mean
)

occ_data_homerange = readRDS(path_homerange_scale_data)[['median']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, cover_forest_diversity, density_edge_cw, density_roads_paved, density_streams_major,
  focalpatch_area_homeange_pcnt, prop_abund_standinit, prop_abund_lsog, prop_abund_comthin, shape_idx
)

occ_data_hs = occ_data_plot_shared %>%
  left_join(occ_data_rs_shared, by = "site") %>%
  left_join(occ_data_plot_hs, by = "site")

occ_data_rs = occ_data_plot_shared %>%
  left_join(occ_data_rs_shared, by = "site") %>%
  left_join(occ_data_plot_rs, by = "site")

occ_data_hr = occ_data_plot_shared %>%
  left_join(occ_data_rs_shared, by = "site") %>%
  left_join(occ_data_homerange, by = "site")

if (model_name == "HS") {
  occ_data = occ_data_hs %>% select(
    -homerange_canopy_cover_mean, -homerange_canopy_cover_cv
  ) %>% filter(hs == TRUE)
} else if (model_name == "RS") {
  occ_data = occ_data_rs %>% filter(hs == TRUE)
} else if (model_name == "HR") {
  occ_data = occ_data_hr %>% filter(hs == TRUE)
} else if (model_name == "RSHR") {
  occ_data = full_join(occ_data_rs, occ_data_hr)
}
occ_data = occ_data %>% mutate(across(where(~ !is.character(.) & !is.logical(.)), as.numeric))

# Check for multicollinearity via pairwise correlation and VIF
cor_matrix = cor(occ_data %>% select(where(is.numeric)), use = "pairwise.complete.obs", method = "pearson")
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= 0.8)
if (nrow(collinearity_candidates) > 0) {
  message(crayon::yellow("WARNING:", nrow(collinearity_candidates), "covariates with high collinearity"))
  print(collinearity_candidates)
}

# Remove correlated variable(s)
if (model_name == "HS") {
  occ_data = occ_data %>% select(
    -plot_treeden_all_hs
  )
} else if (model_name == "RS") {
  occ_data = occ_data %>% select(
    -homerange_canopy_cover_mean, -homerange_treeden_gt4in_dbh_mean
  )
} else if (model_name == "HR") {
  occ_data = occ_data %>% select(
    -homerange_canopy_cover_mean
  )
} else if (model_name == "RSHR") {
  occ_data = occ_data %>% select(
    -homerange_canopy_cover_mean, # ~ homerange_htmax_cv, homerange_qmd_mean, homerange_treeden_gt4in_dbh_mean, density_edge_cw
    -homerange_canopy_cover_cv, # TODO: swap for homerange_htmax_cv?
    -homerange_treeden_gt4in_dbh_mean, # ~ homerange_htmax_cv
    -prop_abund_standinit, # ~ homerange_ba_mean, homerange_qmd_mean
    -homerange_ba_mean,
    -dist_watercourse_major
  )
}
vif_model = lm(rep(1, nrow(occ_data)) ~ ., data = occ_data %>% select(where(is.numeric)))
vif_results = sort(vif(vif_model))
if (max(vif_results) > 10) {
  message(crayon::yellow("WARNING: covariates with high multi-collinearity"))
}
print(vif_results)

# Season `t`
t = "2020"
year = as.numeric(t)

message("Loading community survey data")
community_survey_data = readRDS(path_community_survey_data)
community_survey_data = community_survey_data[, , t, ]
dimnames(community_survey_data)

# Load species-specific thresholds
message("Loading species-specific thresholds")
species_thresholds = read.csv(path_species_thresholds)

# Derive putative observation and survey date matricies
message("Deriving detection-nondetection and yday matricies for each species")
species = dimnames(community_survey_data)[["common_name"]]
ylist   = setNames(vector("list", length(species)), species)
xlist_yday = setNames(vector("list", length(species)), species)

species_discrepancies = sort(c(setdiff(species, species_thresholds$species), setdiff(species_thresholds$species, species)))
if (length(species_discrepancies) > 0) {
  message(crayon::yellow("WARNING:", length(species_discrepancies), "species discrepancies"))
  message(crayon::yellow(paste(species_discrepancies, collapse = ", ")))
}

for (sp in species) {
  # Populate putative observation matrix (detection-nondetection via thresholded confidence scores)
  n_row = dim(community_survey_data)[1]
  n_col = dim(community_survey_data)[2]
  dim_names = dimnames(community_survey_data)[1:2]
  species_data = community_survey_data[, , sp]
  
  if (is.na(naive_threshold)) {
    # Use probablistic thresholding
    if (sp %in% species_thresholds$species) {
      sp_threshdata = species_thresholds %>% filter(species == sp)
      model     = sp_threshdata %>% pull(model)
      threshold = sp_threshdata %>% pull(threshold)
      if (model == "source") {
        mat_obs = matrix(
          unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence_source >= threshold, na.rm = TRUE)) else NA)),
          nrow = n_row, ncol = n_col, dimnames = dim_names)
      } else if (model == "target") {
        mat_obs = matrix(
          unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence_target >= threshold, na.rm = TRUE)) else NA)),
          nrow = n_row, ncol = n_col, dimnames = dim_names)
      }
    } else {
      # There is no threshold for this species
      threshold = NA
      mat_obs = matrix(
        unlist(lapply(species_data, function(x) if (!is.null(x)) 0.0 else NA)),
        nrow = n_row, ncol = n_col, dimnames = dim_names)
    }
  } else {
    # Use naive arbitrary thresholding with both models (i.e. "confidence" column)
    mat_obs = matrix(
      unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence >= naive_threshold, na.rm = TRUE)) else NA)),
      nrow = n_row, ncol = n_col, dimnames = dim_names)
  }
  ylist[[sp]] = mat_obs # Store the resulting putative observations
}

# Survey date matrix (day of year)
x_yday = matrix(
  unlist(lapply(community_survey_data[, , 1], function(x) if (!is.null(x)) yday(x$survey_date) else NA)),
  nrow = dim(community_survey_data)[1],
  ncol = dim(community_survey_data)[2],
  dimnames = dimnames(community_survey_data)[1:2]
)

# Discard sites with no environmental data
sites_missing_environmental_data = setdiff(dimnames(community_survey_data)$site, occ_data$site)
if (length(sites_missing_environmental_data) > 0) {
  message("Discarding ", length(sites_missing_environmental_data), " sites with missing environmental data")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
  
  x_yday = x_yday[!rownames(x_yday) %in% sites_missing_environmental_data, ]
}
# Discard sites with no observations
sites_with_environmental_data_missing_observations = setdiff(occ_data$site, dimnames(community_survey_data)$site)
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

# Naive species detections
naive_occurrence = sapply(ylist, function(mat) { sum(apply(mat, 1, function(x) any(x == 1, na.rm = TRUE))) })
naive_occurrence = data.frame(species = names(naive_occurrence), nsites = naive_occurrence) %>%
  arrange(desc(nsites)) %>% mutate(species = factor(species, levels = rev(species))) %>% mutate(prob = nsites / length(sites))
message(naive_occurrence %>% filter(nsites > 0) %>% nrow(), " species detected")
message("Most commonly detected species:")
print(naive_occurrence %>% slice_max(prob, n=1) %>% pull(species) %>% as.character())
message("Least commonly detected species:")
print(naive_occurrence %>% slice_min(prob, n=1) %>% pull(species) %>% as.character())
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

min_sites_detected = 1
message("Excluding species that were detected at fewer than the minimum number of sites (", min_sites_detected, "):")
species_to_remove = naive_occurrence %>% filter(nsites < min_sites_detected) %>% pull(species) %>% as.character() %>% sort()
# TODO: incorporate these additional species determined to be absent via manual review above into the naive statistics
print(species_to_remove)
ylist[species_to_remove] = NULL
species = names(ylist)

# Naive species occurrence
naive_occurrence = naive_occurrence %>% filter(species %in% names(ylist))
message("Naive occurrence rate: mean ", round(mean(naive_occurrence$prob),2),
        ", min ", round(min(naive_occurrence$prob),2), ", max ", round(max(naive_occurrence$prob),2))
p = ggplot(naive_occurrence, aes(x = prob, y = species)) +
  geom_point(stat = "identity") +
  geom_vline(xintercept = mean(naive_occurrence$prob), color = "blue") +
  labs(title = "Naive species occurrence", x = "Proportion of sites detected as occurrence", y = ""); print(p)

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

# Get detection covariate data
detection_data = readRDS("data/cache/detection_covariates/data_detection.rds") %>% filter(year == t)

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

occ_data = occ_data %>% filter(site %in% dimnames(y)$site) # discard data for irrelevant sites
stopifnot(dimnames(y)$site == occ_data$site) # check that covariate data are aligned with observation matrix by site

# Assemble occurrence covariate data
if (model_name == "HS") {
  param_alpha_names = c(
    "elev",
    "plot_downvol_hs",
    "plot_ht_cv_hs",
    "plot_lcr_hs",
    "plot_qmd_all_hs",
    "plot_snagden_hs",
    "plot_tree_all_diversity",
    "plot_treeden_gt10cmDbh_hs",
    "plot_treeden_lt10cmDbh_hs",
    "plot_understory_vol"
  )
} else if (model_name == "RS") {
  param_alpha_names = c(
    "elev",
    "homerange_downvol_mean",
    "homerange_htmax_cv",
    "homerange_treeden_all_mean",
    "homerange_qmd_mean"
  )
} else if (model_name == "HR") {
  param_alpha_names = c(
    "elev",
    "cover_forest_diversity",
    "density_edge_cw",
    "density_roads_paved",
    "density_streams_major",
    "focalpatch_area_homeange_pcnt",
    "prop_abund_comthin",
    "prop_abund_lsog",
    "prop_abund_standinit",
    "shape_idx"
  )
} else if (model_name == "RSHR") {
  param_alpha_names = c(
    "elev",
    "cover_forest_diversity",
    "density_edge_cw",
    "density_roads_paved",
    "density_streams_major",
    "focalpatch_area_homeange_pcnt",
    "homerange_htmax_cv",
    "homerange_qmd_mean",
    "homerange_treeden_all_mean",
    "prop_abund_comthin",
    "prop_abund_lsog",
    "shape_idx"
  )
}
# Store alpha parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_alpha_data = tibble(param = paste0("alpha", 1:length(param_alpha_names)), name  = param_alpha_names)
param_alpha_data = param_alpha_data %>% rowwise() %>% mutate(scaled = list(scale(occ_data[[name]]))) %>% ungroup()
n_alpha_params = nrow(param_alpha_data)

# Assemble detection covariate data
detect_data = list(
  yday        = x_yday,
  prcp_mm_day = x_prcp,
  tmax_deg_c  = x_tmax
)
# Store beta parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_beta_data = tibble(param = paste0("beta", seq_along(detect_data)), name = names(detect_data))
param_beta_data = param_beta_data %>% rowwise() %>% mutate(scaled = list(scale(as.vector(detect_data[[name]])))) %>% ungroup()
n_beta_params = nrow(param_beta_data)

# Initialize latent occupancy state z[i] as 1 if a detection occurred at site i, and 0 otherwise
z = matrix(data = NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (i in 1:length(sites)) {
  for (sp in 1:length(species)) { z[i,sp] = sum(y[i, , sp], na.rm = TRUE) }
}
z = (z > 0) * 1

# Prepare all data for the model
I = length(species)
J = length(sites)
K = as.vector(n_surveys_per_site)
msom_data = list(
  y = y, # observed (detection-nondetection) data matrix
  I = I, # number of species observed
  J = J, # number of sites sampled
  K = K  # number of sampling periods (surveys) per site
)
for (a in seq_len(n_alpha_params)) { # Add occupancy covariates
  msom_data[[paste0("x_", param_alpha_data$param[a])]] <- as.vector(param_alpha_data$scaled[[a]])
}
for (b in seq_len(n_beta_params)) { # Add detection covariates
  mat_dim <- dim(detect_data[[param_beta_data$name[b]]])
  msom_data[[paste0("x_", param_beta_data$param[b])]] <- array(param_beta_data$scaled[[b]], dim = mat_dim)
}
str(msom_data)

# Specify hierarchical model and write to file
# Following:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x
# - https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2293
# - https://www.sciencedirect.com/science/article/abs/pii/S0006320709004819
model_template = "
model{

  ## Community level hyperpriors

  # Community mean occurrence
  psi.mean ~ dunif(0,1)                   # probability scale
  mu.u <- log(psi.mean) - log(1-psi.mean) # logit scale
  sigma.u ~ dunif(0,5)                    # standard deviation
  tau.u <- pow(sigma.u,-2)                # precision
  
  # Covariate effects on occurrence
  __OCCURRENCE_HYPERPRIORS__
  
  # Community mean detection
  p.mean  ~ dunif(0,1)                 # probability scale
  mu.v  <- log(p.mean) - log(1-p.mean) # logit scale
  sigma.v ~ dunif(0,5)                 # standard deviation
  tau.v <- pow(sigma.v,-2)             # precision

  # Covariate effects on detection
  mu.beta1    ~ dnorm(0,0.01)
  sigma.beta1 ~ dunif(0,5)
  tau.beta1  <- pow(sigma.beta1,-2)
  mu.beta2    ~ dnorm(0,0.01)
  sigma.beta2 ~ dunif(0,5)
  tau.beta2  <- pow(sigma.beta2,-2)
  mu.beta3    ~ dnorm(0,0.01)
  sigma.beta3 ~ dunif(0,5)
  tau.beta3  <- pow(sigma.beta3,-2)

  for (i in 1:I) { # for each species
  
      # Species level priors for occupancy coefficients (note that dnorm in JAGS is parametrized with precision [tau], not sd [sigma])
      u[i] ~ dnorm(mu.u, tau.u)
      __SPECIES_OCCURRENCE_PRIORS__
  
      # Species level priors for detection coefficients
      v[i] ~ dnorm(mu.v, tau.v)
      beta1[i] ~ dnorm(mu.beta1,tau.beta1)
      beta2[i] ~ dnorm(mu.beta2,tau.beta2)
      beta3[i] ~ dnorm(mu.beta3,tau.beta3)
  
      for (j in 1:J) { # for each site
        
          # Ecological process model for latent occurrence z
          logit(psi[j,i]) <- u[i] __OCCURRENCE_COVARIATE_EQ__
          z[j,i] ~ dbern(psi[j,i])
          
          for (k in 1:K[j]) { # for each sampling period (survey) at site j
      
              # Observation model for observed data y
              logit(p[j,k,i]) <- v[i] + beta1[i]*x_beta1[j,k] + beta2[i]*x_beta2[j,k] + beta3[i]*x_beta3[j,k]
              p.obs[j,k,i] <- p[j,k,i]*z[j,i]
              y[j,k,i] ~ dbern(p.obs[j,k,i])
              
              # Create simulated dataset and calculate discrepancies to inform bayesian p-value
              y.sim[j,k,i] ~ dbern(p.obs[j,k,i])
              d.obs[j,k,i] <- pow(abs(y[j,k,i] - p.obs[j,k,i]), 2) # TODO: consider alternative deviance residuals (e.g. Broms et al. 2016)
              d.sim[j,k,i] <- pow(abs(y.sim[j,k,i] - p.obs[j,k,i]), 2)
          }
          d.obs.sum[j,i] <- sum(d.obs[j,1:K[j],i]) 
          d.sim.sum[j,i] <- sum(d.sim[j,1:K[j],i])
      }
  }
  
  ## Derived quantities
  
  # Discrepancy measure between observed and simulated data is defined as mean(D.obs > D.sim)
  p.dobs <- sum(d.obs.sum[1:J,1:I])
  p.dsim <- sum(d.sim.sum[1:J,1:I])
  
  # Estimated number of occuring sites per species (among the sampled population of sites)
  for (i in 1:I) {
    Nocc[i] <- sum(z[ ,i])
  }
  
  # Estimated number of occuring species per site (among the species that were detected anywhere)
  for (j in 1:J) {
    Nsite[j] <- sum(z[j, ])
  }
}
"
# Dynamically generate model specs for alpha parameters
p_alpha = param_alpha_data$param
occurrence_hyperpriors = paste(paste0(
  "mu.",    p_alpha, " ~ dnorm(0,0.01)\n",
  "sigma.", p_alpha, " ~ dunif(0,5)\n",
  "tau.",   p_alpha, " <- pow(sigma.", p_alpha, ",-2)\n"
), collapse = "")
species_occurrence_priors = paste(paste0(
  p_alpha, "[i] ~ dnorm(mu.", p_alpha, ",tau.", p_alpha, ")\n"
), collapse = "")
occurrence_covariate_eq = paste(paste0(
  " + ", p_alpha, "[i]*x_", p_alpha, "[j]"
), collapse = "")

model_spec = model_template
model_spec = gsub("__OCCURRENCE_HYPERPRIORS__", occurrence_hyperpriors, model_spec)
model_spec = gsub("__SPECIES_OCCURRENCE_PRIORS__", species_occurrence_priors, model_spec)
model_spec = gsub("__OCCURRENCE_COVARIATE_EQ__", occurrence_covariate_eq, model_spec)

cat(strsplit(model_spec, "\n")[[1]], sep = "\n") # Print model specification to console

model_file = tempfile()
writeLines(model_spec, con = model_file)

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list(z = z) }, # initial values to avoid data/model conflicts
            parameters.to.save = c( # monitored parameters
              "psi", "z",
              "mu.u", "sigma.u", "u",
              "mu.v", "sigma.v", "v",
              "p.dobs", "p.dsim",
              "Nsite", "Nocc",
              paste0("mu.alpha", 1:n_alpha_params), paste0("sigma.alpha", 1:n_alpha_params), paste0("alpha", 1:n_alpha_params),
              paste0("mu.beta",  1:n_beta_params),  paste0("sigma.beta",  1:n_beta_params),  paste0("beta",  1:n_beta_params)
            ),
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 3000, n.burnin = 1000, n.thin = 1,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

## Convergence and goodness-of-fit diagnostics, checking chains for mixing and convergence with trace and density plots
# https://m-clark.github.io/bayesian-basics/diagnostics.html#monitoring-convergence

# Gelman-Rubin statistic (i.e. "potential scale reduction factor") Rhat values serve as a convergence diagnostic (Gelman and Rubin 1992). This is a test statistic for testing if the variance within chains is different than the variance between chains, and is meant to test if each chain was sampling from similar distributions -- if all the chains are “the same”, then the between chain variation should be close to zero. Rhat values substantially above 1.0 indicate lack of convergence; 1.2 is sometimes used as a guideline for “approximate convergence” (Brooks and Gelman 1998), but in practice a more stringent rule of Rhat < 1.1 is often used to declare convergence. If the chains have not converged, Bayesian credible intervals based on the t-distribution are too wide, and have the potential to shrink by this factor if the MCMC run is continued.
msom_summary = summary(msom)
msom_summary = msom_summary %>% as_tibble() %>%
  mutate(param = rownames(summary(msom)), overlap0 = as.factor(overlap0)) %>% relocate(param, .before = 1) %>%
  mutate(prob = plogis(mean), prob_lower95 = plogis(`2.5%`), prob_upper95 = plogis(`97.5%`))
rhat_threshold = 1.1
suspected_nonconvergence = msom_summary %>% filter(Rhat >= rhat_threshold) %>% filter(!str_starts(param, "z\\[") & !str_starts(param, "psi\\["))
if (nrow(suspected_nonconvergence) > 1) {
  message("The following ", nrow(suspected_nonconvergence), " parameters may not have converged:")
  print(suspected_nonconvergence)
} else {
  message("All parameters appear to have converged (rhat < ", rhat_threshold, ")")
}

if (generate_diagnostic_plots) {
  # Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value.
  MCMCtrace(msom$samples, excl = c('p.dobs', 'p.dsim', 'Nsite', 'Nocc'), post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE, pdf = F)
  
  # Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak.
  MCMCtrace(msom$samples, excl = c('p.dobs', 'p.dsim', 'Nsite', 'Nocc'), post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE, pdf = F)
}

# If the model indicates convergence issues, we may need to:
# - Increase the burn-in period
# - Increase iterations
# - Use more informative priors
# - Reparametrize the model

# "We assessed the adequacy of the model using the approach suggested by Gelman et al. (1996) referred to as a Bayesian p-value. We defined a discrepancy measure for the observations yi and their expected values pi under the model. This discrepancy statistic is computed at each iteration of the MCMC algorithm. A reference distribution is computed by simulating data sets from the posterior distribution and computing the discrepancy measure, Dsim, for the simulated data set. The Bayesian p-value is defined as the probability: Pr(D > Dsim). Extreme values (e.g. less than 0.05 or greater than 0.95) indicate that the model is inadequate." (Zipkin et al. 2009) A value of 0.5 indicates good fit, while a value of 0.25 or 0.75 indicates moderate fit.
# "Our P-value, also known as a posterior predictive check, followed the same process as Carrillo-Rubio et al. (2014), Kroll et al. (2014), and Tobler et al. (2015), but was based on the deviance residuals rather than the Pearson's residuals because of its relationship to information criterion theory (Spiegelhalter et al. 1998). We describe the approach in detail in Appendix S1." (TODO: compare to Broms et al. 2016)
# "If the observed data set is consistent with the model in question, then the Bayesian p-value should be close to 0.50. In practice, a p-value close to 0 or 1 indicates that the model is inadequate in some way -- close to 0 suggests a lack of fit and close to 1 suggests that the model over-fits the data, which may occur when it is too complex." (MacKenzie et al. 2018)
p_dobs = msom$sims.list$p.dobs
p_dsim = msom$sims.list$p.dsim
p_val  = mean(p_dobs > p_dsim) # proportion of samples for which the model fits the observed data worse than the simulated data
message("Bayesian p-value: ", round(p_val, 3))

## Model predictive performance (in-sample ROC AUC) for model selection

# "We evaluated model performance by computing the area under the curve of the receiver operating characteristic (AUC). When applied to the dataset used during model construction, AUC measures a model’s goodness-of-fit by estimating the probability that a randomly chosen occupied sampling point (where zij=1) has a higher probability of occupancy than a randomly chosen unoccupied sampling point (where zij=0). If a model fits well, then it consistently predicts a higher probability of occupancy for occupied sites yielding an AUC closer to 1.0. Conversely, if a model fits poorly, it will perform the same as chance yielding an AUC closer to 0.5. We utilized AUC for evaluating our model in two ways..." (Mattsson et al. 2013)
# "The AUC (ranging from 0–1) measures the discriminatory ability of a model, which in this case corresponds to the ability to correctly project which areas are occupied. A value of 0.5 indicates that the model performs no better than random. Values greater than 0.5 indicate progressively better discriminatory capabilities (Hosmer and Lemeshow 2000)." ()

# Get posterior samples (i.e. MCMC simulated draws) for estimated occurrence psi and latent state z
psi_samples = msom$sims.list$psi
z_samples   = msom$sims.list$z
# These should be of dimension: samples x sites x species
n_samples = dim(psi_samples)[1]
stopifnot(identical(J, dim(psi_samples)[2]))
stopifnot(identical(I, dim(psi_samples)[3]))
stopifnot(identical(J, dim(z_samples)[2]))
stopifnot(identical(I, dim(z_samples)[3]))

# "First, we calculated mean and 95% Bayesian credibility interval (BCI) AUC values reflecting goodness-of-fit for the model based on the vector of AUC values across MCMC iterations (henceforth, consolidated AUC values) for all species combined."
# This effectively estimates the probability that a randomly chosen occupied species-site combination has a higher probability of occurrence than a random unoccupied one. As such, it reflects an overal measure of model fit. Note that the number of site-species pairs is dominated by common species, so these can disproportionately influence the combined AUC value. This is why it's important to also quantify species-specific AUC.

auc_combined = rep(NA, n_samples)
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = n_samples, clear = FALSE)
for (s in 1:n_samples) {
  psi_all = as.vector(psi_samples[s, , ]) # combine all species into a single vector
  z_all   = as.vector(z_samples[s, , ])
  
  if (length(unique(z_all)) > 1) {
    auc_combined[s] = as.numeric(pROC::auc(pROC::roc(z_all, psi_all, quiet=TRUE)))
  } else {
    stop("Cannot calculate ROC") # latent z states are identical for this sample, something is wrong
  }
  pb$tick()
}
auc_combined_mean = round(mean(auc_combined, na.rm=TRUE), 3)
auc_combined_bci  = round(quantile(auc_combined, probs=c(0.025, 0.975), na.rm=TRUE), 3)
auc_combined = data.frame(
  mean_auc = auc_combined_mean,
  lower95_auc = auc_combined_bci[[1]],
  upper95_auc = auc_combined_bci[[2]]
)

message("Combined mean AUC: ", auc_combined_mean, " (95% BCI ", auc_combined_bci[1], "–", auc_combined_bci[2], ")")

# "Second, we calculated AUC values reflecting model goodness-of-fit for each species under each model rendering a mean and 95% BCI AUC value for each species-model combination (henceforth, species-specific AUC values). Calculating AUCs for each species in each model is made possible by examining the species-specific binary occupancy predictions along with predicted occupancy probabilities for each sampling point across the respective vectors of MCMC iterations."
# Check for outlier species with poor fit. Note that AUC is less stable when there are fewer positive occurrence examples for a species. Interpret rare species with caution.

auc_species = matrix(NA, nrow=n_samples, ncol=I)
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = n_samples, clear = FALSE)
for (s in 1:n_samples) {
  for (i in 1:I) {
    psi_vec = psi_samples[s, , i]
    z_vec   = z_samples[s, , i]

    if (length(unique(z_vec)) > 1) {
      auc_species[s, i] = as.numeric(pROC::auc(pROC::roc(z_vec, psi_vec, quiet=TRUE)))
    } else {
      # species[i] latent z state is identical at all sites for this sample, cannot calculate ROC
    }
  }
  pb$tick()
}
auc_species_mean = round(apply(auc_species, 2, mean, na.rm=TRUE), 3)
auc_species_bci  = round(apply(auc_species, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE), 3)
nocc_samples = msom$sims.list$Nocc # posterior mean number of occupied sites (z=1) per species (i.e. effective positive sample size)
nocc_mean = round(apply(nocc_samples, 2, mean),1)
nocc_bci = apply(nocc_samples, 2, quantile, probs = c(0.025, 0.975))
auc_species = data.frame(
  species = species[1:I],
  mean_auc = auc_species_mean,
  lower95_auc = auc_species_bci["2.5%", ],
  upper95_auc = auc_species_bci["97.5%", ],
  mean_nocc = nocc_mean,
  lower95_nocc = nocc_bci["2.5%", ],
  upper95_nocc = nocc_bci["97.5%", ]
)

message("Species-specific mean AUC: ", round(mean(auc_species$mean_auc, na.rm = TRUE),3),
        " (range ",  round(min(auc_species$mean_auc, na.rm = TRUE),3), "–", round(max(auc_species$mean_auc, na.rm = TRUE),3), ")")
print(auc_species)
ggplot(auc_species, aes(x = mean_nocc, y = mean_auc)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin = lower95_auc, ymax = upper95_auc), width = 0, color = "gray") +
  geom_errorbarh(aes(xmin = lower95_nocc, xmax = upper95_nocc), height = 0, color = "gray") +
  geom_point() +
  scale_y_continuous(breaks = seq(0.5, 1.0, by = 0.1)) +
  geom_text_repel(aes(label = species), size = 2) +
  labs(
    x = "Mean estimated number of occurring sites", y = "Mean ROC AUC",
    title = "Species-specific ROC AUC and effective positive sample size"
  )

# "We concluded a statistically significant difference between posterior distributions when the 95% BCI (2.5th to 97.5th percentile of the posterior distribution) for one posterior excluded the BCI of the opposing posterior."

# Write results to cache
msom_results = list(
  msom_summary = msom_summary,
  p_val        = p_val,
  auc_combined = auc_combined,
  auc_species  = auc_species,
  param_alpha_data = param_alpha_data,
  param_beta_data = param_beta_data,
  sites = sites,
  species = species
)
path_out = paste0("data/cache/models/msom_", model_name, "_", format(Sys.Date(), "%Y-%m-%d"), ".rds")
if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
saveRDS(msom_results, file = path_out)
message(crayon::green("Cached model and results to ", path_out))

## Explore model parameters and covariate effects

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
occurrence_coeff_summary = msom_summary %>%
    filter(str_detect(param, "^mu\\.alpha|^sigma\\.alpha")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)
detection_coeff_summary = msom_summary %>%
    filter(str_detect(param, "^mu\\.beta|^sigma\\.beta")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)
print(occurrence_coeff_summary)
print(detection_coeff_summary)

# Compare community level effect sizes for occurrence coefficients
occurrence_effect_sizes = full_join(occurrence_coeff_summary %>% filter(str_starts(param, "mu")), param_alpha_data %>% mutate(param = paste0("mu.", param)), by='param')
plt = ggplot(occurrence_effect_sizes, aes(x = mean, y = as.factor(name))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
  scale_color_manual(values = c("black", "gray")) +
  labs(title = "Community level effect sizes for occurrence covariates", x = "Coefficient estimate", y = "Parameter"); print(plt)

# Compare species level effects of each covariate on occurrence
for (alpha_param in param_alpha_data$param) {
  alpha_name = param_alpha_data %>% filter(param == alpha_param) %>% pull(name)
  alpha_coef = msom_summary %>% filter(str_detect(param, paste0("^", alpha_param, "(?!\\d)", "\\["))) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
    mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
  mu_alpha_summary = msom_summary %>% filter(param == paste0("mu.", alpha_param)) %>% select(mean, `2.5%`, `97.5%`)
  plt = ggplot(alpha_coef, aes(x = as.factor(plot_order), y = mean)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_point(aes(color = overlap0)) +
    geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = overlap0), width = 0) +
    geom_hline(yintercept = mu_alpha_summary$mean,    linetype = "solid",  color = "blue") +
    geom_hline(yintercept = mu_alpha_summary$`2.5%`,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = mu_alpha_summary$`97.5%`, linetype = "dashed", color = "blue") +
    labs(title = paste("Effect of", alpha_name, "on occurrence"), x = "Species", y = paste(alpha_param, "coefficient estimate")) +
    scale_x_discrete(labels = alpha_coef$species_name) +
    scale_color_manual(values = c("0" = "black", "1" = "gray")) +
    coord_flip() +
    theme(legend.position = "none"); print(plt)
}

# Compare community level effect sizes for detection coefficients
detection_effect_sizes = full_join(detection_coeff_summary %>% filter(str_starts(param, "mu")), param_beta_data %>% mutate(param = paste0("mu.", param)), by='param')
plt = ggplot(detection_effect_sizes, aes(x = mean, y = as.factor(name))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
  scale_color_manual(values = c("black", "gray")) +
  labs(title = "Community level effect sizes for detection covariates", x = "Coefficient estimate", y = "Parameter"); print(plt)

# Compare species level effects of each covariate on detection
for (beta_param in param_beta_data$param) {
  beta_name = param_beta_data %>% filter(param == beta_param) %>% pull(name)
  beta_coef = msom_summary %>% filter(str_detect(param, paste0("^", beta_param, "(?!\\d)", "\\["))) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
    mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
  mu_beta_summary = msom_summary %>% filter(param == paste0("mu.", beta_param)) %>% select(mean, `2.5%`, `97.5%`)
  plt = ggplot(beta_coef, aes(x = as.factor(plot_order), y = mean)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_point(aes(color = overlap0)) +
    geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = overlap0), width = 0) +
    geom_hline(yintercept = mu_beta_summary$mean,    linetype = "solid",  color = "blue") +
    geom_hline(yintercept = mu_beta_summary$`2.5%`,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = mu_beta_summary$`97.5%`, linetype = "dashed", color = "blue") +
    labs(title = paste("Effect of", beta_name, "on detection"), x = "Species", y = paste(beta_param, "coefficient estimate")) +
    scale_x_discrete(labels = beta_coef$species_name) +
    scale_color_manual(values = c("black", "gray")) +
    coord_flip() +
    theme(legend.position = "none"); print(plt)
}

# TODO: Save results
# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from=plots.png.paths, to="data/temp")

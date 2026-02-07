####################################################################################
# Multi-species occupancy model using multiple detection method described by Miller et al. 2011

model_name = "MDM"
naive_threshold = NA # set to NA to use probabilistic threshold from 1c_calculate_species_thresholds.R
generate_diagnostic_plots = FALSE
####################################################################################

path_community_survey_data = "data/cache/1_derive_community_survey_data/community_survey_data_2025-08-15.rds"
path_species_thresholds    = "data/cache/1_calculate_species_thresholds/species_thresholds.csv"
path_plot_scale_data       = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_homerange_scale_data  = "data/cache/occurrence_covariates/data_homerange_scale.rds"
path_out = paste0("data/cache/models/", model_name, "_", format(Sys.Date(), "%yU-%m-%d"), ".rds")

library(progress)
library(car)
library(tidyverse)
library(jagsUI)
library(MCMCvis)
library(glue)
library(ggplot2)
library(ggrepel)
theme_set(theme_classic())

####################################################################################
# Load occurrence covariate and survey data

# Get occurrence data (local plot and homerange scales)
occ_data_plot_shared = readRDS(path_plot_scale_data) %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, hs, elev, dist_watercourse_major
)
occ_data_rs_shared = readRDS(path_homerange_scale_data)[['plot']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, homerange_canopy_cover_mean, homerange_canopy_cover_cv
)
occ_data_plot_rs = readRDS(path_homerange_scale_data)[['plot']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, homerange_downvol_mean, homerange_htmax_cv, homerange_qmd_mean, homerange_treeden_all_mean, homerange_treeden_gt4in_dbh_mean, homerange_ba_mean, homerange_snagden_gt15dbh_mean
)
occ_data_homerange = readRDS(path_homerange_scale_data)[['median']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, cover_forest_diversity, density_edge_cw, density_roads_paved, density_streams_major,
  focalpatch_area_homeange_pcnt, prop_abund_standinit, prop_abund_lsog, prop_abund_comthin, shape_idx
)
occ_data_rs = occ_data_plot_shared %>% left_join(occ_data_rs_shared, by = "site") %>% left_join(occ_data_plot_rs, by = "site")
occ_data_hr = occ_data_plot_shared %>% left_join(occ_data_rs_shared, by = "site") %>% left_join(occ_data_homerange, by = "site")
occ_data = full_join(occ_data_rs, occ_data_hr)
occ_data = occ_data %>% mutate(across(where(~ !is.character(.) & !is.logical(.)), as.numeric))

# Select occupancy covariates
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
occ_data = occ_data %>% select(site, all_of(param_alpha_names))

# Test for multicollinearity
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

####################################################################################
# Load species-specific thresholds
message("Loading species-specific thresholds")
species_thresholds = read.csv(path_species_thresholds)
species_thresholds_source = species_thresholds %>% filter(model == 'source')
species_thresholds_target = species_thresholds %>% filter(model == 'target')

# Identify species only supported by the source model
species_only_in_source = species_thresholds_source %>% filter(!species %in% species_thresholds_target$species)

# Manually choose model for specific species
target_species = species_thresholds_target %>% filter(species %in% c(
  "marbled murrelet",
  "western screech-owl",
  "sooty grouse",
  "northern pygmy-owl",
  "western wood-pewee",
  "red-breasted nuthatch",
  "northern saw-whet owl",
  "white-crowned sparrow",
  "townsend's warbler",
  "dark-eyed junco",
  "hermit thrush",
  "golden-crowned kinglet",
  "song sparrow",
  "band-tailed pigeon",
  "pileated woodpecker",
  "rufous hummingbird",
  "red crossbill"
))
source_species = species_thresholds_source %>% filter(model == 'source') %>% filter(species %in% c(
  "ruby-crowned kinglet",
  "violet-green swallow",
  "american robin",
  "wilson's warbler",
  "spotted towhee",
  "purple finch",
  "olive-sided flycatcher",
  "western tanager",
  "hutton's vireo",
  "black-throated gray warbler",
  "varied thrush",
  "pacific-slope flycatcher",
  "pacific wren",
  "swainson's thrush",
  "barred owl",
  "belted kingfisher",
  "hairy woodpecker",
  "northern flicker",
  "hammond's flycatcher",
  "common raven"
))
species_thresholds_manual_selection = rbind(species_only_in_source, source_species)
species_thresholds_manual_selection = rbind(species_thresholds_manual_selection, target_species)
stopifnot(nrow(species_thresholds_manual_selection) == length(species_thresholds_source$species))

# Set minimum threshold to 0.5
species_thresholds_manual_selection = species_thresholds_manual_selection %>%
  mutate(
    threshold = ifelse(n_pos == 0, NA, ifelse(t_conf_tp >= 0.5, t_conf_tp, 0.5)),
    precision = ifelse(t_conf_tp >= 0.5, precision_tp, precision_0.5),
    recall = ifelse(t_conf_tp >= 0.5, recall_tp, recall_0.5)
  ) %>% select(species, model, threshold, precision, recall, auc_pr, auc_roc, n_pos, n_neg) %>% arrange(species)
# print(species_thresholds_manual_selection)
species_thresholds = species_thresholds_manual_selection

####################################################################################
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
  print(sites_missing_environmental_data)
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
  
  x_yday = x_yday[!rownames(x_yday) %in% sites_missing_environmental_data, ]
}
# Discard sites with no observations
sites_with_environmental_data_missing_observations = setdiff(occ_data$site, dimnames(community_survey_data)$site)
if (length(sites_with_environmental_data_missing_observations) > 0) {
  message("Discarding ", length(sites_with_environmental_data_missing_observations), " sites with missing observations")
  print(sites_with_environmental_data_missing_observations)
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
# p = ggplot(naive_occurrence, aes(x = nsites, y = species)) +
#   geom_bar(stat = "identity") +
#   geom_vline(xintercept = mean(naive_occurrence$nsites), color = "blue") +
#   labs(title = "Species occurrence across sites", x = "Number of sites detected as occurrence", y = ""); print(p)

min_sites_detected = 1
message("Excluding species that were detected at fewer than the minimum number of sites (", min_sites_detected, "):")
species_to_remove = naive_occurrence %>% filter(nsites < min_sites_detected) %>% pull(species) %>% as.character() %>% sort()
# TODO: incorporate these additional species determined to be absent via manual review above into the naive statistics
print(species_to_remove)
ylist[species_to_remove] = NULL
species = names(ylist)

# Format observation detection-nondetection and covariate data for modeling as 3D arrays (site × survey × species)

# Unconfirmed observation data
yU = array(NA, dim = c(length(sites), length(surveys), length(species)),
          dimnames = list(site = sites, survey = surveys, species = species))
for (i in seq_along(ylist)) {
  yU[, , i] = as.matrix(ylist[[i]])
}

# Left-align data (moving any missing NA surveys to the right) to allow for direct indexing by number of surveys per site
y_unaligned = yU
left_align_row = function(x) {
  non_na = x[!is.na(x)]
  c(non_na, rep(NA, length(x) - length(non_na)))
}
for (i in dimnames(yU)[[3]]) {
  sp_y_mat = yU[, , i]
  sp_y_aligned = t(apply(sp_y_mat, 1, left_align_row))
  dimnames(sp_y_aligned) = dimnames(sp_y_mat)
  yU[, , i] = sp_y_aligned
}
n_surveys_per_site = apply(!is.na(yU[, , 1]), 1, sum)

x_yday_unaligned = x_yday
x_yday = t(apply(x_yday_unaligned, 1, left_align_row))
dimnames(x_yday) = dimnames(x_yday_unaligned)

# Populate confirmed observation data using yday
site_key = read.csv("data/sites/site_key_long.csv") %>% mutate(site = tolower(site), site_agg = tolower(site_agg))
site_key$date = as.Date(site_key$date, format = "%m/%d/%y")
yC = array(NA, dim = dim(yU), dimnames = dimnames(yU))

confirmed_detections = arrow::read_parquet("data/debug/predictions_source.parquet") %>% arrange(file, label_predicted)
for (i in dimnames(yC)[[3]]) {
  print(i)
  species_cd = confirmed_detections %>% filter(label_truth == i)
  species_cd = species_cd %>% mutate(
    parsed = str_match(file, "_([A-Za-z0-9]+)_(\\d{8})_(\\d{6})$"),
    serialno = parsed[,2],
    date = parsed[,3] %>% ymd(),
    time = parsed[,4]
  )
  species_cd$yday = yday(species_cd$date)
  # Get associated site
  species_cd = left_join(species_cd, site_key %>% filter(year == t), by = c('serialno', 'date'))
  species_yC = x_yday
  lookup <- with(species_cd, paste(site, yday))
  species_yC_new <- species_yC
  for (r in seq_len(nrow(species_yC))) {
    for (c in seq_len(ncol(species_yC))) {
      site <- rownames(species_yC)[r]
      yday <- species_yC[r, c]
      if (!is.na(yday) && paste(site, yday) %in% lookup) {
        species_yC_new[r, c] <- 1
      } else {
        species_yC_new[r, c] <- NA
      }
    }
  }
  yC[, , i] = species_yC_new
}

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

# Assemble occurrence covariate data
occ_data = occ_data %>% filter(site %in% dimnames(yU)$site) # discard data for irrelevant sites
stopifnot(dimnames(yU)$site == occ_data$site) # check that covariate data are aligned with observation matrix by site

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

# Initialize latent occupancy state z[j,i] as 1 if a detection OR CONFIRMATION occurred at site i, and 0 otherwise
z = matrix(data = NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (j in 1:length(sites)) {
  for (i in 1:length(species)) { z[j,i] = sum(yU[j, , i], na.rm = TRUE) + sum(yC[j, , i], na.rm = TRUE) }
}
z = (z > 0) * 1

# Prepare all data for the model
I = length(species)
J = length(sites)
K = as.vector(n_surveys_per_site)
msom_data = list(
  yU = yU, # unconfirmed (classifier) observation data matrix
  yC = yC, # confirmed (human) observation data matrix
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

####################################################################################
# Specify hierarchical model and write to file

# TODO: Use validation data on classifier false-positive rates to set informative priors for mu.eta0/sigma.eta0 or for species eta0[i]
model_template = "
model{

  ## Community level hyperpriors

  # Community mean occurrence
  psi.mean ~ dunif(0,1)                   # probability scale
  mu.u <- log(psi.mean) - log(1-psi.mean) # logit scale
  sigma.u ~ dunif(0,5)                    # standard deviation
  tau.u <- pow(sigma.u,-2)                # precision
  
  # Covariate effects on occurrence
  mu.alpha1 ~ dnorm(0,0.01)
  sigma.alpha1 ~ dunif(0,5)
  tau.alpha1 <- pow(sigma.alpha1,-2)
  mu.alpha2 ~ dnorm(0,0.01)
  sigma.alpha2 ~ dunif(0,5)
  tau.alpha2 <- pow(sigma.alpha2,-2)
  mu.alpha3 ~ dnorm(0,0.01)
  sigma.alpha3 ~ dunif(0,5)
  tau.alpha3 <- pow(sigma.alpha3,-2)
  mu.alpha4 ~ dnorm(0,0.01)
  sigma.alpha4 ~ dunif(0,5)
  tau.alpha4 <- pow(sigma.alpha4,-2)
  mu.alpha5 ~ dnorm(0,0.01)
  sigma.alpha5 ~ dunif(0,5)
  tau.alpha5 <- pow(sigma.alpha5,-2)
  mu.alpha6 ~ dnorm(0,0.01)
  sigma.alpha6 ~ dunif(0,5)
  tau.alpha6 <- pow(sigma.alpha6,-2)
  mu.alpha7 ~ dnorm(0,0.01)
  sigma.alpha7 ~ dunif(0,5)
  tau.alpha7 <- pow(sigma.alpha7,-2)
  mu.alpha8 ~ dnorm(0,0.01)
  sigma.alpha8 ~ dunif(0,5)
  tau.alpha8 <- pow(sigma.alpha8,-2)
  mu.alpha9 ~ dnorm(0,0.01)
  sigma.alpha9 ~ dunif(0,5)
  tau.alpha9 <- pow(sigma.alpha9,-2)
  mu.alpha10 ~ dnorm(0,0.01)
  sigma.alpha10 ~ dunif(0,5)
  tau.alpha10 <- pow(sigma.alpha10,-2)
  
  # Community hyperpriors for unconfirmed detection method (classifier)
  pU.mean  ~ dunif(0,1)                 # probability scale
  mu.vU  <- log(pU.mean) - log(1-pU.mean) # logit scale
  sigma.vU ~ dunif(0,5)                 # standard deviation
  tau.vU <- pow(sigma.vU,-2)             # precision

  # Covariate effects on detection (unconfirmed classifier)
  mu.beta1    ~ dnorm(0,0.01)
  sigma.beta1 ~ dunif(0,5)
  tau.beta1  <- pow(sigma.beta1,-2)
  mu.beta2    ~ dnorm(0,0.01)
  sigma.beta2 ~ dunif(0,5)
  tau.beta2  <- pow(sigma.beta2,-2)
  mu.beta3    ~ dnorm(0,0.01)
  sigma.beta3 ~ dunif(0,5)
  tau.beta3  <- pow(sigma.beta3,-2)
  
  # Community hyperpriors for confirmed detection method (human verification)
  pC.mean ~ dunif(0,1)
  mu.vC  <- log(pC.mean) - log(1-pC.mean)
  sigma.vC ~ dunif(0,5)
  tau.vC <- pow(sigma.vC,-2)
  
  # Hyperpriors for false-positive unconfirmed detection method (classifier)
  mu.eta0 ~ dnorm(0,0.01)
  sigma.eta0 ~ dunif(0,5)
  tau.eta0 <- pow(sigma.eta0,-2)

  for (i in 1:I) { # for each species
  
      # Species level priors for occupancy coefficients (note that dnorm in JAGS is parametrized with precision [tau], not sd [sigma])
      u[i] ~ dnorm(mu.u, tau.u)
      alpha1[i] ~ dnorm(mu.alpha1,tau.alpha1)
      alpha2[i] ~ dnorm(mu.alpha2,tau.alpha2)
      alpha3[i] ~ dnorm(mu.alpha3,tau.alpha3)
      alpha4[i] ~ dnorm(mu.alpha4,tau.alpha4)
      alpha5[i] ~ dnorm(mu.alpha5,tau.alpha5)
      alpha6[i] ~ dnorm(mu.alpha6,tau.alpha6)
      alpha7[i] ~ dnorm(mu.alpha7,tau.alpha7)
      alpha8[i] ~ dnorm(mu.alpha8,tau.alpha8)
      alpha9[i] ~ dnorm(mu.alpha9,tau.alpha9)
      alpha10[i] ~ dnorm(mu.alpha10,tau.alpha10)
  
      # Species level priors for unconfirmed detection covariates (classifier)
      vU[i] ~ dnorm(mu.vU, tau.vU)
      beta1[i] ~ dnorm(mu.beta1,tau.beta1)
      beta2[i] ~ dnorm(mu.beta2,tau.beta2)
      beta3[i] ~ dnorm(mu.beta3,tau.beta3)
      
      # Species level priors for confirmed detection (human)
      vC[i] ~ dnorm(mu.vC, tau.vC)
      
      # Species level priors for unconfirmed detection false-positive covariates (classifier)
      eta0[i] ~ dnorm(mu.eta0,tau.eta0)
  
      for (j in 1:J) { # for each site
        
          # Ecological process model for latent occurrence z
          logit(psi[j,i]) <- u[i] + alpha1[i]*x_alpha1[j] + alpha2[i]*x_alpha2[j] + alpha3[i]*x_alpha3[j] + alpha4[i]*x_alpha4[j] + alpha5[i]*x_alpha5[j] + alpha6[i]*x_alpha6[j] + alpha7[i]*x_alpha7[j] + alpha8[i]*x_alpha8[j] + alpha9[i]*x_alpha9[j] + alpha10[i]*x_alpha10[j]
          z[j,i] ~ dbern(psi[j,i])
          
          for (k in 1:K[j]) { # for each sampling period (survey) at site j
          
              # Note that if yU or yU contains NA, JAGS should ignore it
      
              ## Observation mixture model for unconfirmed data yU
              
              # p11 is true-positive detection probability given z=1
              logit(p11[j,k,i]) <- vU[i] + beta1[i]*x_beta1[j,k] + beta2[i]*x_beta2[j,k] + beta3[i]*x_beta3[j,k]
              
              # p10 is false-positive detection probability given z=0
              logit(p10[j,k,i]) <- eta0[i] # TODO: extend with relevant covariates that are site or survey-specific
            
              muU[j,k,i] <- z[j,i] * p11[j,k,i] + (1 - z[j,i]) * p10[j,k,i]
              yU[j,k,i] ~ dbern(muU[j,k,i])
              
              ## Observation model for confirmed data yU
              
              logit(pC[j,k,i]) <- vC[i]
              muC[j,k,i] <- z[j,i] * pC[j,k,i]
              yC[j,k,i] ~ dbern(muC[j,k,i])
          }
      }
  }
}
"
model_spec = model_template
cat(strsplit(model_spec, "\n")[[1]], sep = "\n") # print model specification to console
model_file = tempfile()
writeLines(model_spec, con = model_file)

####################################################################################
# Run JAGS

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

# TODO: for eta0, vU, vC, reasonable initialize near 0 (logit(0.5))?
msom = jags(data = msom_data,
            inits = function() { list(z = z) }, # initial values to avoid data/model conflicts
            parameters.to.save = c( # monitored parameters
              "psi", "z",
              "mu.u",  "sigma.u",  "u",
              "mu.vU", "sigma.vU", "vU",
              "mu.vC", "sigma.vC", "vC",
              "mu.eta0", "sigma.eta0", "eta0",
              paste0("mu.alpha", 1:n_alpha_params), paste0("sigma.alpha", 1:n_alpha_params), paste0("alpha", 1:n_alpha_params),
              paste0("mu.beta",  1:n_beta_params),  paste0("sigma.beta",  1:n_beta_params),  paste0("beta",  1:n_beta_params)
            ),
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 3000, n.burnin = 1000, n.thin = 1,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

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

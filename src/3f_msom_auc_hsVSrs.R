####################################################################################
# A multi-species static occupancy model with predictive performance measured via ROC AUC for model selection
# Comparing a set of related strucutral variables calculated at the habitat survey plot and 100-m remote sensing scales
#
# Formulated according to:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x#b14
# - https://wildlife.onlinelibrary.wiley.com/doi/10.1002/jwmg.442
#
# INPUT:
HS_OR_RS = "HS" # "HS" or "RS"
#
path_community_survey_data = "data/cache/1_derive_community_array/community_survey_data.rds"
path_plot_scale_data = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_homerange_scale_data = "data/cache/occurrence_covariates/data_homerange_scale.rds"
####################################################################################

generate_diagnostic_plots = FALSE

library(progress)
library(car)
library(tidyverse)
library(jagsUI)
library(MCMCvis)
library(ggplot2)
library(ggrepel)
theme_set(theme_classic())

# Get occurrence data (local plot and homerange scales)
occ_data_plot_shared = readRDS(path_plot_scale_data) %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, hs, elev
)

occ_data_rs_shared = readRDS(path_homerange_scale_data)[['min']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site #, homerange_canopy_cover_mean
)

occ_data_plot_hs = readRDS(path_plot_scale_data) %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, plot_downvol_hs, plot_ht_cv_hs, plot_ht_hs, plot_qmd_all_hs, plot_snagden_hs, plot_treeden_all_hs, plot_understory_vol, plot_tree_all_diversity
)

occ_data_plot_rs = readRDS(path_homerange_scale_data)[['min']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, homerange_downvol_mean, homerange_htmax_cv, homerange_qmd_mean, homerange_treeden_all_mean
)

occ_data_homerange = readRDS(path_homerange_scale_data)[['median']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, prop_abund_standinit, prop_abund_undstryreinit, prop_abund_oldgrowth, prop_abund_comthin,
  focalpatch_area_homeange_pcnt, focalpatch_isolation,
  density_edge_cw, shape_idx,
  density_roads_paved, density_streams_major
)

occ_data_hs = occ_data_plot_shared %>%
  left_join(occ_data_rs_shared, by = "site") %>%
  left_join(occ_data_plot_hs, by = "site") %>%
  left_join(occ_data_homerange, by = "site") %>% filter(hs == TRUE)

occ_data_rs = occ_data_plot_shared %>%
  left_join(occ_data_rs_shared, by = "site") %>%
  left_join(occ_data_plot_rs, by = "site") %>%
  left_join(occ_data_homerange, by = "site")

if (HS_OR_RS == "HS") {
  occ_data = occ_data_hs
} else {
  occ_data = occ_data_rs
}
occ_data = occ_data %>% mutate(across(where(~ !is.character(.) & !is.logical(.)), as.numeric))

# Check for multicollinearity via pairwise correlation and VIF
cor_matrix = cor(occ_data %>% select(where(is.numeric)), use = "pairwise.complete.obs", method = "pearson")
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= 0.8))

if (HS_OR_RS == "HS") {
  occ_data = occ_data %>% select(
    -prop_abund_standinit,   # plot_ht_hs
  )
} else {
  occ_data = occ_data %>% select(
    -prop_abund_standinit,   # homerange_qmd_mean @ 100m
  )
}
vif_model = lm(rep(1, nrow(occ_data)) ~ ., data = occ_data %>% select(where(is.numeric)))
sort(vif(vif_model))


# Season `t`
t = "2020"
year = as.numeric(t)
threshold = 0.75 # naive conversion from continuous score prediction to binary detection-nondetection

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
sites_missing_environmental_data = setdiff(dimnames(community_survey_data)$unit, occ_data$site)
if (length(sites_missing_environmental_data) > 0) {
  message("Discarding ", length(sites_missing_environmental_data), " sites with missing environmental data")
  ylist = lapply(ylist, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
  
  x_yday = x_yday[!rownames(x_yday) %in% sites_missing_environmental_data, ]
}
# Discard sites with no observations
sites_with_environmental_data_missing_observations = setdiff(occ_data$site, dimnames(community_survey_data)$unit)
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
species_to_remove = c(species_to_remove, "Great Horned Owl", "Hermit Warbler", "Pine Grosbeak", "Townsend's Solitaire", "Cooper's Hawk", "Bald Eagle", "Merlin", "Cassin's Vireo", "Sharp-shinned Hawk", "Clark's Nutcracker", "Fox Sparrow", "Spotted Owl", "American Crow", "Osprey", "Black-capped Chickadee", "House Finch", "American Dipper", "Anna's Hummingbird", "Common Merganser", "Northern Rough-winged Swallow", "Barn Swallow", "American Kestrel", "Chipping Sparrow", "Tree Swallow", "Lincoln's Sparrow", "Willow Flycatcher", "White-throated Sparrow", "Western Bluebird", "Bewick's Wren", "Horned Lark", "Common Yellowthroat", "Savannah Sparrow", "Barn Owl", "Peregrine Falcon", "Eurasian Collared-Dove", "Marbled Murrelet") # TODO: incorporate these additional species determined to be absent via manual review above into the naive statistics
cat(species_to_remove, "\n")
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

# Standardize occurrence covariate data
params_alpha_names = data.frame()
if (HS_OR_RS == "HS") {
  x_alpha1_scaled = scale(occ_data$plot_elev_rs)
  x_alpha2_scaled = scale(occ_data$plot_downvol_hs)
  x_alpha3_scaled = scale(occ_data$plot_ht_cv_hs)
  x_alpha4_scaled = scale(occ_data$plot_ht_hs)
  x_alpha5_scaled = scale(occ_data$plot_qmd_all_hs)
  x_alpha6_scaled = scale(occ_data$plot_snagden_hs)
  x_alpha7_scaled = scale(occ_data$plot_treeden_all_hs)
  x_alpha8_scaled = scale(occ_data$plot_understory_vol)
  x_alpha9_scaled = scale(occ_data$tree_all_diversity)
  x_alpha10_scaled = scale(occ_data$focalpatch_pcnt)
  x_alpha11_scaled = scale(occ_data$cw_edge_density)
  x_alpha12_scaled = scale(occ_data$prop_abund_5)
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha1", name = "plot_elev_rs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha2", name = "plot_downvol_hs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha3", name = "plot_ht_cv_hs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha4", name = "plot_ht_hs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha5", name = "plot_qmd_all_hs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha6", name = "plot_snagden_hs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha7", name = "plot_treeden_all_hs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha8", name = "plot_understory_vol"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha9", name = "tree_all_diversity"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha10", name = "focalpatch_pcnt"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha11", name = "cw_edge_density"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha12", name = "prop_abund_5"))
} else {
  x_alpha1_scaled = scale(occ_data$plot_elev_rs)
  x_alpha2_scaled = scale(occ_data$homerange_downvol_mean)
  x_alpha3_scaled = scale(occ_data$homerange_htmax_cv)
  x_alpha4_scaled = scale(occ_data$shape_idx)
  x_alpha5_scaled = scale(occ_data$homerange_qmd_mean)
  x_alpha6_scaled = scale(occ_data$prop_abund_4)
  x_alpha7_scaled = scale(occ_data$homerange_treeden_all_mean)
  x_alpha8_scaled = scale(occ_data$density_roads_paved)
  x_alpha9_scaled = scale(occ_data$density_streams_major)
  x_alpha10_scaled = scale(occ_data$focalpatch_pcnt)
  x_alpha11_scaled = scale(occ_data$cw_edge_density)
  x_alpha12_scaled = scale(occ_data$prop_abund_5)
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha1", name = "plot_elev_rs"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha2", name = "homerange_downvol_mean"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha3", name = "homerange_htmax_cv"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha4", name = "shape_idx"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha5", name = "homerange_qmd_mean"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha6", name = "prop_abund_4"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha7", name = "homerange_treeden_all_mean"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha8", name = "density_roads_paved"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha9", name = "density_streams_major"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha10", name = "focalpatch_pcnt"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha11", name = "cw_edge_density"))
  params_alpha_names = rbind(params_alpha_names, data.frame(param = "alpha12", name = "prop_abund_5"))
}

# Standardize detection covariate data to z-scale (mean 0, standard deviation 1)
params_beta_names = data.frame()
x_beta1_scaled = scale(as.vector(x_yday))
x_beta2_scaled = scale(as.vector(x_prcp))
x_beta3_scaled = scale(as.vector(x_tmax))
params_beta_names = rbind(params_beta_names, data.frame(param = "beta1", name = "yday"))
params_beta_names = rbind(params_beta_names, data.frame(param = "beta2", name = "prcp_mm_day"))
params_beta_names = rbind(params_beta_names, data.frame(param = "beta3", name = "tmax_deg_c"))

# Initialize latent occupancy state z[i] as 1 if a detection occurred at unit i, and 0 otherwise
z = matrix(data = NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (i in 1:length(sites)) {
  for (sp in 1:length(species)) { z[i,sp] = sum(y[i, , sp], na.rm = TRUE) }
}
z = (z > 0) * 1

# Prepare all data for the model
n_alpha_params = 12
n_beta_params = 3
I = length(species)
J = length(sites)
msom_data = list(
  # Observed data
  y = y,                             # detection-nondetection matrix
  I = I,                             # number of species observed
  J = J,                             # number of sites sampled
  K = as.vector(n_surveys_per_site), # number of sampling periods (surveys) per site
  # Occupancy covariates
  x_alpha1   = x_alpha1_scaled[,1],
  x_alpha2   = x_alpha2_scaled[,1],
  x_alpha3   = x_alpha3_scaled[,1],
  x_alpha4   = x_alpha4_scaled[,1],
  x_alpha5   = x_alpha5_scaled[,1],
  x_alpha6   = x_alpha6_scaled[,1],
  x_alpha7   = x_alpha7_scaled[,1],
  x_alpha8   = x_alpha8_scaled[,1],
  x_alpha9   = x_alpha9_scaled[,1],
  x_alpha10   = x_alpha10_scaled[,1],
  x_alpha11   = x_alpha11_scaled[,1],
  x_alpha12   = x_alpha12_scaled[,1],
  # Detection covariates
  x_beta1 = array(x_beta1_scaled, dim = dim(x_yday)),
  x_beta2 = array(x_beta2_scaled, dim = dim(x_prcp)),
  x_beta3 = array(x_beta3_scaled, dim = dim(x_tmax))
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

  # Community mean occurrence
  psi.mean ~ dunif(0,1)                   # probability scale
  mu.u <- log(psi.mean) - log(1-psi.mean) # logit scale
  sigma.u ~ dunif(0,5)                    # standard deviation
  tau.u <- pow(sigma.u,-2)                # precision
  
  # Covariate effects on occurrence
  mu.alpha1     ~ dnorm(0,0.01)
  sigma.alpha1  ~ dunif(0,5)
  tau.alpha1   <- pow(sigma.alpha1,-2)
  mu.alpha2     ~ dnorm(0,0.01)
  sigma.alpha2  ~ dunif(0,5)
  tau.alpha2   <- pow(sigma.alpha2,-2)
  mu.alpha3     ~ dnorm(0,0.01)
  sigma.alpha3  ~ dunif(0,5)
  tau.alpha3   <- pow(sigma.alpha3,-2)
  mu.alpha4     ~ dnorm(0,0.01)
  sigma.alpha4  ~ dunif(0,5)
  tau.alpha4   <- pow(sigma.alpha4,-2)
  mu.alpha5     ~ dnorm(0,0.01)
  sigma.alpha5  ~ dunif(0,5)
  tau.alpha5   <- pow(sigma.alpha5,-2)
  mu.alpha6     ~ dnorm(0,0.01)
  sigma.alpha6  ~ dunif(0,5)
  tau.alpha6   <- pow(sigma.alpha6,-2)
  mu.alpha7     ~ dnorm(0,0.01)
  sigma.alpha7  ~ dunif(0,5)
  tau.alpha7   <- pow(sigma.alpha7,-2)
  mu.alpha8     ~ dnorm(0,0.01)
  sigma.alpha8  ~ dunif(0,5)
  tau.alpha8   <- pow(sigma.alpha8,-2)
  mu.alpha9     ~ dnorm(0,0.01)
  sigma.alpha9  ~ dunif(0,5)
  tau.alpha9   <- pow(sigma.alpha9,-2)
  mu.alpha10     ~ dnorm(0,0.01)
  sigma.alpha10  ~ dunif(0,5)
  tau.alpha10   <- pow(sigma.alpha10,-2)
  mu.alpha11     ~ dnorm(0,0.01)
  sigma.alpha11  ~ dunif(0,5)
  tau.alpha11   <- pow(sigma.alpha11,-2)
  mu.alpha12     ~ dnorm(0,0.01)
  sigma.alpha12  ~ dunif(0,5)
  tau.alpha12   <- pow(sigma.alpha12,-2)
  
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
      alpha11[i] ~ dnorm(mu.alpha11,tau.alpha11)
      alpha12[i] ~ dnorm(mu.alpha12,tau.alpha12)
  
      # Species level priors for detection coefficients
      v[i] ~ dnorm(mu.v, tau.v)
      beta1[i] ~ dnorm(mu.beta1,tau.beta1)
      beta2[i] ~ dnorm(mu.beta2,tau.beta2)
      beta3[i] ~ dnorm(mu.beta3,tau.beta3)
  
      for (j in 1:J) { # for each site
        
          # Ecological process model for latent occurrence z
          logit(psi[j,i]) <- u[i] + alpha1[i]*x_alpha1[j] + alpha2[i]*x_alpha2[j] + alpha3[i]*x_alpha3[j] + alpha4[i]*x_alpha4[j] + alpha5[i]*x_alpha5[j] + alpha6[i]*x_alpha6[j] + alpha7[i]*x_alpha7[j] + alpha8[i]*x_alpha8[j] + alpha9[i]*x_alpha9[j] + alpha10[i]*x_alpha10[j] + alpha11[i]*x_alpha11[j] + alpha12[i]*x_alpha12[j]
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
", con = model_file)

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
              paste0("mu.beta", 1:n_beta_params), paste0("sigma.beta", 1:n_beta_params), paste0("beta", 1:n_beta_params)
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
suspected_nonconvergence = msom_summary %>% filter(Rhat >= rhat_threshold) %>% filter(!str_starts(param, "z\\["))
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

message("Species-specific mean AUC: ", round(mean(auc_species$mean_auc),3),
        " (range ",  round(min(auc_species$mean_auc),3), "–", round(max(auc_species$mean_auc),3), ")")
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
(occurrence_coeff_summary = msom_summary %>%
    filter(str_detect(param, "^mu\\.alpha|^sigma\\.alpha")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0))
(detection_coeff_summary = msom_summary %>%
    filter(str_detect(param, "^mu\\.beta|^sigma\\.beta")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0))

# Compare community level effect sizes for occurrence coefficients
occurrence_effect_sizes = full_join(occurrence_coeff_summary %>% filter(str_starts(param, "mu")), params_alpha_names %>% mutate(param = paste0("mu.", param)), by='param')
plt = ggplot(occurrence_effect_sizes, aes(x = mean, y = as.factor(name))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
  scale_color_manual(values = c("black", "gray")) +
  labs(title = "Community level effect sizes for occurrence covariates", x = "Coefficient estimate", y = "Parameter"); print(plt)

# Compare species level effects of each covariate on occurrence
for (alpha_param in params_alpha_names$param) {
  alpha_name = params_alpha_names %>% filter(param == alpha_param) %>% pull(name)
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
detection_effect_sizes = full_join(detection_coeff_summary %>% filter(str_starts(param, "mu")), params_beta_names %>% mutate(param = paste0("mu.", param)), by='param')
plt = ggplot(detection_effect_sizes, aes(x = mean, y = as.factor(name))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
  scale_color_manual(values = c("black", "gray")) +
  labs(title = "Community level effect sizes for detection covariates", x = "Coefficient estimate", y = "Parameter"); print(plt)

# Compare species level effects of each covariate on detection
for (beta_param in params_beta_names$param) {
  beta_name = params_beta_names %>% filter(param == beta_param) %>% pull(name)
  beta_coef = msom_summary %>% filter(str_detect(param, paste0("^", beta_param, "(?!\\d)", "\\["))) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
    mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
  mu_beta_summary = msom_summary %>% filter(param == paste0("mu.", beta_param)) %>% select(mean, `2.5%`, `97.5%`)
  plt = ggplot(beta_coef, aes(x = as.factor(plot_order), y = mean)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_point(aes(color = overlap0)) +
    geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, size = 1) +
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

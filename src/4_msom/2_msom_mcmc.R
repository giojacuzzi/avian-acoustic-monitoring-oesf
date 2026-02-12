## 2_msom_mcmc.R #########################################################################################
# Fit a multispecies occupancy model using MCMC
#
## CONFIG:
grouping = "all" # Species grouping ("all", "diet"...)
param_alpha_names = c(
  "elev",
  # "homerange_htmax_cv",
  "dist_road_paved",
  "dist_watercourse_major",
  "homerange_qmd_mean",
  "homerange_treeden_all_mean"
)
param_delta_names = c(
  # "cover_forest_diversity",
  "shape_idx",
  # "density_roads_paved", # TODO: density_roads_all?
  # "density_streams_major",
  "focalpatch_area_homeange_pcnt",
  "prop_abund_comthin",
  "prop_abund_lsog",
  "prop_abund_standinit"
)
param_beta_names = c(
  "yday",
  "prcp_mm_day",
  "tmax_deg_c"
)
OVERRIDE_SPECIES_SCALE_WITH_MEDIAN = TRUE # DEBUG
#
## OUTPUT:
dir_out = "data/cache/models"
#
## INPUT:
model_file                               = "src/4_msom/msom.txt"
path_y                                   = "data/cache/4_msom/1_assemble_msom_data/y.rds"
path_occurrence_predictor_plot_data      = "data/cache/4_msom/1_assemble_msom_data/occurrence_predictor_plot_data.rds"
path_occurrence_predictor_homerange_data = "data/cache/4_msom/1_assemble_msom_data/occurrence_predictor_homerange_data.rds"
path_detection_predictor_data            = "data/cache/4_msom/1_assemble_msom_data/detection_predictor_data.rds"
path_trait_data                          = "data/cache/trait_data/trait_data.csv"
###########################################################################################################

source("src/global.R")

model_name = tools::file_path_sans_ext(basename(model_file))
path_out = paste0(dir_out, "/", model_name, "_", grouping, "_", format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), ".rds")
if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)

# Load dependencies ----------------------------------------------------------------------------

message("Loading occurrence predictor plot scale data from ", path_occurrence_predictor_plot_data)
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)

message("Loading occurrence predictor homerange scale data from ", path_occurrence_predictor_homerange_data)
occurrence_predictor_homerange_data = readRDS(path_occurrence_predictor_homerange_data)

message("Loading detection predictor data from ", path_detection_predictor_data)
detection_predictor_data = readRDS(path_detection_predictor_data)

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)

message("Loading species detection histories 'y' from ", path_y)
y = readRDS(path_y)

species = dimnames(y)$species
seasons = dimnames(y)$season
surveys = dimnames(y)$survey
sites   = dimnames(y)$site

# Prepare occurrence plot scale covariate data ----------------------------------------------------------------------

# Store alpha parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_alpha_data = tibble(param = paste0("alpha", 1:length(param_alpha_names)), name  = param_alpha_names)
param_alpha_data = param_alpha_data %>% rowwise() %>% mutate(scaled = list(scale(occurrence_predictor_plot_data[[name]]))) %>% ungroup()
n_alpha_params = nrow(param_alpha_data)

# Prepare occurrence homerange scale covariate data ----------------------------------------------------------------------

# Standardize species names of occ homerange data
names(occurrence_predictor_homerange_data)[names(occurrence_predictor_homerange_data) == "western flycatcher"] = "pacific-slope flycatcher"

# Store minimum (i.e. floor) plot homerange data for those species who have smaller estimated home range sizes
if (!OVERRIDE_SPECIES_SCALE_WITH_MEDIAN) {
  occ_data_homerange_floor = occurrence_predictor_homerange_data[['plot']] %>% filter(site %in% dimnames(y)$site)
  stopifnot(dimnames(y)$site == occ_data_homerange_floor$site)
}

# Store median plot homerange data for debugging
occ_data_homerange_median = occurrence_predictor_homerange_data[['median']] %>% filter(site %in% dimnames(y)$site)
stopifnot(dimnames(y)$site == occ_data_homerange_median$site)

# Store delta parameter ID, variable name, and standardized data in species-specific matricies
delta_data = setNames(vector('list', length(param_delta_names)), param_delta_names)
for (param in param_delta_names) {
  message(param)
  
  species_delta_data = setNames(vector('list', length(species)), species)
  
  for (i in 1:length(species)) {
    species_name = species[i]
    # message(i, " ", species_name)
    if (OVERRIDE_SPECIES_SCALE_WITH_MEDIAN) {
      param_data = occ_data_homerange_median %>% select(site, all_of(param))
    } else {
      # discard data for irrelevant sites
      species_occ_data = occ_data_homerange[[species_name]] %>% filter(site %in% dimnames(y)$site)
      # check that data are aligned with observation matrix by site
      stopifnot(dimnames(y)$site == species_occ_data$site)
      if (unique(species_occ_data %>% pull(buffer_radius_m)) >= 100) {
        param_data = species_occ_data %>% select(site, all_of(param))
      } else {
        param_data = occ_data_homerange_floor %>% select(site, all_of(param))
      }
    }
    species_delta_data[[species_name]] = scale(param_data[[param]])
  }
  delta_data[[param]] = species_delta_data
}
param_delta_data = tibble(param = paste0("delta", 1:length(param_delta_names)), name  = param_delta_names)
param_delta_data = param_delta_data %>% rowwise() %>% mutate(data = list(delta_data[[name]])) %>% ungroup()
n_delta_params = nrow(param_delta_data)

# Prepare detection covariate data ----------------------------------------------------------------------------

# Store beta parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_beta_data = tibble(param = paste0("beta", seq_along(detection_predictor_data)), name = param_beta_names)
param_beta_data = param_beta_data %>% rowwise() %>%
  mutate(scaled = list(scale(unlist(detection_predictor_data[[name]], use.names = FALSE)))) %>% ungroup()
n_beta_params = nrow(param_beta_data)

param_season_data = tibble(param = "season", name = "season", scaled = list(scale(1:length(seasons))))

# Assign species group membership ------------------------------------------------------------------------------

group_col = paste0("group_", grouping)
st = species_traits %>% filter(common_name %in% species)
groups = st %>%
  mutate(
    group       = st[[group_col]],
    group_idx   = as.integer(factor(group)),
    grouping    = grouping
  ) %>%
  select(common_name, group, group_idx, grouping)

# Order to match species
groups = groups %>% arrange(match(common_name, species))

# Package all model data -------------------------------------------------------------------------------

message("Packaging MSOM data")

n_surveys_per_site = list() # number of non-NA surveys per site (for Kmax)
for (t in dimnames(y)[["season"]]) {
  n_surveys_per_site[[t]] = apply(!is.na(y[, , t, 1]), 1, sum)
}

# Model data constants and covariates
J = length(sites)
K = as.matrix(as.data.frame(lapply(n_surveys_per_site, as.vector)))
Kmax = max(K)
T = length(seasons)
I = length(species)
G = length(unique(groups$group_idx))

# "JAGS dcat() requires integer categories 1..L. Miller's paper uses y = 0 (no detect), 1 (uncertain), 2 (certain). You must recode your y input to JAGS so the values are 1=no detect, 2=uncertain, 3=certain. I will assume that mapping in the code below. If you prefer to keep 0/1/2 in R, transform before sending data to JAGS: y_jags <- y + 1."
y = y + 1

msom_data = list(
  J = J,                    # number of sites sampled
  K = K,                    # number of secondary sampling periods (surveys) per site per season (site x season)
  Kmax = Kmax,              # maximum number of surveys across all sites and seasons
  T = T,                    # number of primary sampling periods (seasons)
  I = I,                    # number of species observed
  G = G,                    # number of species groups
  group = groups$group_idx, # species group membership
  y = y                     # observed (detection-nondetection) data matrix
)
# Add alpha parameters
for (a in seq_len(n_alpha_params)) {
  msom_data[[paste0("x_", param_alpha_data$param[a])]] <- as.vector(param_alpha_data$scaled[[a]])
}
# Add delta parameters
for (d in seq_len(n_delta_params)) {
  param_data_new = do.call(cbind, param_delta_data %>% filter(name == param_delta_data$name[d]) %>% pull(data) %>% .[[1]])
  msom_data[[paste0("x_", param_delta_data$param[d])]] = param_data_new
}
# Add season parameter
msom_data[["x_season"]] = as.vector(param_season_data$scaled[[1]])

# Add beta parameters
for (b in seq_len(n_beta_params)) {
  arr = detection_predictor_data[[param_beta_data$name[1]]]
  vec = as.vector(param_beta_data$scaled[[b]]) # flatten to 1D
  
  stopifnot(length(vec) == prod(dim(arr)))  # safety check
  
  msom_data[[paste0("x_", param_beta_data$param[b])]] = array(
    vec,
    dim = dim(arr),
    dimnames = dimnames(arr)
  )
}

str(msom_data)

# Initialize latent occupancy state z[i] as 1 if a detection occurred at site i, and 0 otherwise
z = array(NA, dim = c(length(sites), length(seasons), length(species)), dimnames = list(sites, seasons, species))
for (j in seq_along(sites)) {
  for (t in seq_along(seasons)) {
    for (i in seq_along(species)) {
      z[j, t, i] = (sum(y[j, , t, i], na.rm = TRUE) > 0) * 1
    }
  }
}

# Run JAGS ----------------------------------------------------------------------------------------------------

message("\n", "System CPU: "); print(as.data.frame(t(benchmarkme::get_cpu())))
message("System RAM: "); print(benchmarkme::get_ram())

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list( # initial values to avoid data/model conflicts
              z = z
            ) },
            parameters.to.save = c( # monitored parameters
              "mu.u", "sigma.u", "u",
              "mu.v", "sigma.v", "v",
              "mu.w", "sigma.w", "w",
              "mu.b", "sigma.b", "b",
              paste0("mu.alpha", 1:n_alpha_params), paste0("sigma.alpha", 1:n_alpha_params), paste0("alpha", 1:n_alpha_params),
              paste0("mu.delta", 1:n_delta_params), paste0("sigma.delta", 1:n_delta_params), paste0("delta", 1:n_delta_params),
              paste0("mu.beta",  1:n_beta_params),  paste0("sigma.beta",  1:n_beta_params),  paste0("beta",  1:n_beta_params),
              "mu.alphaseason", "sigma.alphaseason", "alphaseason",
              "mu.betaseason", "sigma.betaseason", "betaseason",
              "D_obs", "D_sim", "z"
            ),
            model.file = model_file,
            n.chains = 3, n.adapt = 200, n.iter = 2000, n.burnin = 200, n.thin = 1, parallel = FALSE, # ETA: 11 hr
            # n.chains = 3, n.adapt = 5000, n.iter = 30000, n.burnin = 10000, n.thin = 3, parallel = TRUE, # ETA: TODO
            DIC = FALSE, verbose=TRUE)

message("Finished running JAGS (", round(msom$mcmc.info$elapsed.mins / 60, 2), " hr)")

message("MCMC information:")
print(data.frame(
  n.chains = msom$mcmc.info$n.chains,
  n.adapt  = msom$mcmc.info$n.adapt[1],
  n.iter   = msom$mcmc.info$n.iter,
  n.burnin = msom$mcmc.info$n.burnin,
  n.thin   = msom$mcmc.info$n.thin,
  samples  = msom$mcmc.info$n.samples
))

message("Raw model size: ", format(utils::object.size(msom), units = "auto"))

# Retrieve summary data and investigate goodness-of-fit ----------------------------------------------------

msom_summary = summary(msom)
msom_summary = msom_summary %>% as_tibble() %>%
  mutate(param = rownames(summary(msom)), overlap0 = as.factor(overlap0)) %>% relocate(param, .before = 1) %>%
  mutate(prob = plogis(mean), prob_lower95 = plogis(`2.5%`), prob_upper95 = plogis(`97.5%`))
rhat_threshold = 1.1
suspected_nonconvergence = msom_summary %>% filter(Rhat >= rhat_threshold) %>% filter(!str_starts(param, "z\\[") & !str_starts(param, "psi\\["))
suspected_nonconvergence = suspected_nonconvergence %>% mutate(
  index = str_extract(param, "(?<=\\[)\\d+(?=\\])"),
  index = as.integer(index),
  species = ifelse(!is.na(index), species[index], NA)
)
if (nrow(suspected_nonconvergence) > 1) {
  message("The following ", nrow(suspected_nonconvergence), " parameters may not have converged:")
  print(suspected_nonconvergence)
} else {
  message("All parameters appear to have converged (rhat < ", rhat_threshold, ")")
}

## Posterior predictive check - Bernoulli deviance contribution (Broms et al. 2016)
# "If the observed data are consistent with the model in question, then the Bayesian p-value should be close to 0.5. In practice, a p-value close to 0 or 1 indicates that the model is inadequate in some way -- close to 0 suggests a lack of fit and close to 1 suggests that the model over-fits the data, which may occur when it is too complex... A Bayesian p-value is calculated as the proportion of times the selected summary statistic calculated for the generated data is greater than the value calculated from the observed data" (MacKenzie et al. 2018)
# Is the overall likelihood of the observed detection histories under the fitted model about the same as the likelihood of new data generated from that model?
# This Bayesian p-value is the probability (proportion of iterations) that the simulated deviance is greater than the observed deviance.
p_val = mean(msom$sims.list$D_sim > msom$sims.list$D_obs)
message("Baysian p-value (deviance): ", round(p_val,3))

if (FALSE) {
  # "Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value."
  MCMCvis::MCMCtrace(msom$samples, ISB = FALSE, pdf = TRUE, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)
  
  # "Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak."
  MCMCvis::MCMCtrace(msom$samples, ISB = FALSE, pdf = TRUE, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)
}

# Inspect the mean and 95% BCI of hyperparameter estimates
whiskerplot(msom, c(paste0('mu.', param_alpha_data$param)))
whiskerplot(msom, c(paste0('mu.', param_delta_data$param)))
whiskerplot(msom, c(paste0('mu.', param_beta_data$param)))

# Write results to cache
msom_results = list(
  msom              = msom,
  msom_summary      = msom_summary,
  p_val             = p_val,
  param_alpha_data  = param_alpha_data,
  param_beta_data   = param_beta_data,
  sites             = sites,
  species           = species,
  groups            = groups
)
saveRDS(msom_results, file = path_out)
message(crayon::green("Cached model and results to", path_out, "-", format(structure(file.info(path_out)$size, class = "object_size"), units = "auto")))


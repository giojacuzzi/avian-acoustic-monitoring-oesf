## 2_msom_mcmc.R #########################################################################################
# Fit a multispecies occupancy model using MCMC with stage as a categorical variable
#
## CONFIG:
grouping = "nest" # Species grouping ("all", "diet", "nest" ...)
model_type = "nofp" # nofp, fp
param_alpha_stage = "stage_3"
param_alpha_season = "season"
param_alpha_point_names = c(
  "elevation",
  "dist_road_paved",
  "dist_watercourse_major"
)
param_alpha_plot_names = c( # stage-specific effects
  "homerange_ba_mean"
)
param_alpha_homerange_names = c(
  "aggregation_idx",
  "density_edge_cw"
)
# param_alpha_names = c(
#   "elev",
#   "dist_road_paved",
#   "dist_watercourse_major",
#   "homerange_qmd_mean",
#   "homerange_treeden_all_mean"
# )
# param_delta_names = c(
#   # "cover_forest_diversity",
#   "shape_idx",
#   # "density_roads_paved", # TODO: density_roads_all?
#   # "density_streams_major",
#   "focalpatch_area_homeange_pcnt",
#   "prop_abund_comthin",
#   "prop_abund_lsog",
#   "prop_abund_standinit"
# )
param_beta_point_names = c(
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
path_y                                   = "data/cache/4_msom/1_assemble_msom_data/y.rds"
path_occurrence_predictor_plot_data      = "data/cache/4_msom/1_assemble_msom_data/occurrence_predictor_plot_data.rds"
path_occurrence_predictor_homerange_data = "data/cache/4_msom/1_assemble_msom_data/occurrence_predictor_homerange_data.rds"
path_detection_predictor_data            = "data/cache/4_msom/1_assemble_msom_data/detection_predictor_data.rds"
path_trait_data                          = "data/cache/trait_data/trait_data.csv"
###########################################################################################################

source("src/global.R")

if (model_type == "nofp") {
  model_file = "src/4_msom/msom_stage_cat.txt"
# } else if (model_type == "fp") {
#   model_file = "src/4_msom/msom_fp.txt"
} else {
  stop("ERROR: Unsupported model type")
}

model_name = tools::file_path_sans_ext(basename(model_file))
path_out = paste0(dir_out, "/", model_name, "_", model_type, "_", grouping, ".rds")
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

# Stage categorical
stages = occurrence_predictor_plot_data %>% select(stage_3) %>% mutate(stage_idx = as.integer(stage_3))

# Season
season_data_list = list((scale(1:length(seasons))))
param_season_data = tibble(
  param = c("alpha_season"),
  name = "season",
  data = c(season_data_list)
)

# Store alpha parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_alpha_point_data = tibble(param = paste0("alpha_point", 1:length(param_alpha_point_names)), name  = param_alpha_point_names)
param_alpha_point_data = param_alpha_point_data %>% rowwise() %>% mutate(data = list(scale(occurrence_predictor_plot_data[[name]]))) %>% ungroup()
n_alpha_point_params = nrow(param_alpha_point_data)

param_alpha_plot_data = tibble(param = paste0("alpha_plot", 1:length(param_alpha_plot_names)), name  = param_alpha_plot_names)
param_alpha_plot_data = param_alpha_plot_data %>% rowwise() %>% mutate(data = list(scale(occurrence_predictor_plot_data[[name]]))) %>% ungroup()
n_alpha_plot_params = nrow(param_alpha_plot_data)

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
alpha_homerange_data = setNames(vector('list', length(param_alpha_homerange_names)), param_alpha_homerange_names)
for (param in param_alpha_homerange_names) {
  message(param)
  
  species_alpha_data = setNames(vector('list', length(species)), species)
  
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
    species_alpha_data[[species_name]] = scale(param_data[[param]])
  }
  alpha_homerange_data[[param]] = species_alpha_data
}
param_alpha_homerange_data = tibble(param = paste0("alpha_homerange", 1:length(param_alpha_homerange_names)), name  = param_alpha_homerange_names)
param_alpha_homerange_data = param_alpha_homerange_data %>% rowwise() %>% mutate(data = list(alpha_homerange_data[[name]])) %>% ungroup()
n_alpha_homerange_params = nrow(param_alpha_homerange_data)

# Prepare detection covariate data ----------------------------------------------------------------------------

# Store beta parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_beta_point_data = tibble(param = paste0("beta_point", seq_along(detection_predictor_data)), name = param_beta_point_names)
param_beta_point_data = param_beta_point_data %>% rowwise() %>%
  mutate(data = list(scale(unlist(detection_predictor_data[[name]], use.names = FALSE)))) %>% ungroup()
n_beta_point_params = nrow(param_beta_point_data)

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
Tseason = length(seasons)
I = length(species)
G = length(unique(groups$group_idx))
S = length(unique(stages$stage_idx))

# "JAGS dcat() requires integer categories 1..L. Miller's paper uses y = 0 (no detect), 1 (uncertain), 2 (certain). You must recode your y input to JAGS so the values are 1=no detect, 2=uncertain, 3=certain. I will assume that mapping in the code below. If you prefer to keep 0/1/2 in R, transform before sending data to JAGS: y_jags <- y + 1."
if (model_type == "fp") {
  y = y + 1
} else {
  y[y == 2] <- 1
}

msom_data = list(
  J = J,                    # number of sites sampled
  K = K,                    # number of secondary sampling periods (surveys) per site per season (site x season)
  Kmax = Kmax,              # maximum number of surveys across all sites and seasons
  Tseason = Tseason,        # number of primary sampling periods (seasons)
  I = I,                    # number of species observed
  G = G,                    # number of species groups
  group = groups$group_idx, # species group membership
  S = S,                    # number of site stages
  stage = stages$stage_idx, # site stage membership
  y = y                     # observed (detection-nondetection) data matrix
)
# Add alpha parameters
for (a in seq_len(n_alpha_point_params)) {
  msom_data[[paste0("x_", param_alpha_point_data$param[a])]] <- as.vector(param_alpha_point_data$data[[a]])
}
for (a in seq_len(n_alpha_plot_params)) {
  msom_data[[paste0("x_", param_alpha_plot_data$param[a])]] <- as.vector(param_alpha_plot_data$data[[a]])
}
for (a in seq_len(n_alpha_homerange_params)) {
  param_data_new = do.call(cbind, param_alpha_homerange_data %>% filter(name == param_alpha_homerange_data$name[a]) %>% pull(data) %>% .[[1]])
  msom_data[[paste0("x_", param_alpha_homerange_data$param[a])]] = param_data_new
}
# Season
msom_data[["x_season"]] = as.vector(param_season_data$data[[1]])

# Add beta parameters
for (b in seq_len(n_beta_point_params)) {
  arr = detection_predictor_data[[param_beta_point_data$name[1]]]
  vec = as.vector(param_beta_point_data$data[[b]]) # flatten to 1D
  
  stopifnot(length(vec) == prod(dim(arr)))  # safety check
  
  msom_data[[paste0("x_", param_beta_point_data$param[b])]] = array(
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

params_to_monitor = c(
  "mu.u", "sigma.u", "u",
  "mu.alpha_stage", "sigma.alpha_stage", "alpha_stage",
  "mu.alpha_season", "sigma.alpha_season", "alpha_season",
  paste0("mu.alpha_point",     1:n_alpha_point_params),     paste0("sigma.alpha_point", 1:n_alpha_point_params),        paste0("alpha_point", 1:n_alpha_point_params),
  # paste0("mu.alpha_plot",      1:n_alpha_plot_params),      paste0("sigma.alpha_plot", 1:n_alpha_plot_params),          paste0("alpha_plot", 1:n_alpha_plot_params),
  paste0("mu.alpha_homerange", 1:n_alpha_homerange_params), paste0("sigma.alpha_homerange", 1:n_alpha_homerange_params), paste0("alpha_homerange", 1:n_alpha_homerange_params),
  "mu.v", "sigma.v", "v",
  paste0("mu.beta_point",      1:n_beta_point_params),      paste0("sigma.beta_point",  1:n_beta_point_params),         paste0("beta_point",  1:n_beta_point_params),
  "D_obs", "D_sim", "z"
)

if (model_type == "fp") {
  params_to_monitor = c(params_to_monitor,
                        "mu.w", "sigma.w", "w",
                        "mu.b", "sigma.b", "b")
}

# Run JAGS ----------------------------------------------------------------------------------------------------

message("Model file:", model_file)

message("\n", "System CPU: "); print(as.data.frame(t(benchmarkme::get_cpu())))
message("System RAM: "); print(benchmarkme::get_ram())

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list( # initial values to avoid data/model conflicts
              z = z
              # TODO
              # "This model has multimodal likelihood, resulting in identical support for different parameter values. We addressed this issue by constraining the parameters so that p11 > p10, an assumption that is supported by the validation  data at the chosen threshold17. We imposed this constraint by providing arbitrarily chosen starting values that  align with this condition (0.7 for p11 and 0.1 for p10)23." 
              # v = rep(logit(0.70), length(species)), # informative priors necessary to avoid invalid PPC log(0) values
              # w = rep(logit(0.05), length(species))
            ) },
            parameters.to.save = params_to_monitor,
            model.file = model_file,
            # Test run, ETA: 11/25 hr (fp) 6 hr (parallel fp), 3 hr (nofp)
            n.chains = 1, n.adapt = 200, n.iter = 2000, n.burnin = 200, n.thin = 1, parallel = FALSE,
            # Formal run, ETA TODO
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

# Write results to cache
msom_results = list(
  msom              = msom,
  msom_summary      = msom_summary,
  p_val             = p_val,
  param_alpha_data  = list(
    param_alpha_stage = param_alpha_stage,
    param_alpha_season = param_alpha_season,
    param_alpha_point_data = param_alpha_point_data,
    param_alpha_plot_data = param_alpha_plot_data,
    param_alpha_homerange_data = param_alpha_homerange_data
  ),
  param_beta_data = list(
    param_beta_point_data = param_beta_point_data
  ),
  sites   = sites,
  seasons = seasons,
  species = species,
  groups  = groups,
  stages  = stages
)
saveRDS(msom_results, file = path_out)
message(crayon::green("Cached model and results to", path_out, "-", format(structure(file.info(path_out)$size, class = "object_size"), units = "auto")))

# INSPECT MSOM ===============================================================================

stop("READY FOR INSPECTION")

model_data = read_rds(path_out)
msom = model_data$msom
stages = model_data$stages
groups = model_data$groups
species = model_data$species
param_alpha_stage = model_data$param_alpha_data$param_alpha_stage
param_alpha_season = model_data$param_alpha_data$param_alpha_season
param_alpha_point_data = model_data$param_alpha_data$param_alpha_point_data
param_alpha_plot_data = model_data$param_alpha_data$param_alpha_plot_data
param_alpha_homerange_data = model_data$param_alpha_data$param_alpha_homerange_data

if (FALSE) {
  # "Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value."
  MCMCvis::MCMCtrace(msom$samples, ISB = FALSE, pdf = TRUE, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)
  
  # "Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak."
  MCMCvis::MCMCtrace(msom$samples, ISB = FALSE, pdf = TRUE, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)
}

match_s = stages %>% distinct() %>% arrange(stage_idx)
match_g = groups %>% select(group, group_idx) %>% distinct() %>% arrange(group_idx)
match_i = tibble(species = species, species_idx = 1:length(species))

# Inspect the mean and 95% BCI of hyperparameter estimates
whiskerplot(msom, 'mu.u')
whiskerplot(msom, 'mu.alpha_stage')
# whiskerplot(msom, c(paste0('mu.', param_alpha_plot_data$param)))
whiskerplot(msom, c(paste0('mu.', param_alpha_point_data$param)))
whiskerplot(msom, c(paste0('mu.', param_alpha_homerange_data$param)))
whiskerplot(msom, 'mu.alpha_season')

whiskerplot(msom, 'mu.v')
whiskerplot(msom, c(paste0('mu.', param_beta_point_data$param)))

{
  mu_alpha_stage = summary(msom) %>% as.data.frame() %>% rownames_to_column("param") %>%
    filter(str_starts(param, "mu.alpha_stage\\[")) %>%
    mutate(
      s = as.integer(str_match(param, "mu\\.alpha_stage\\[(\\d+),(\\d+)\\]")[,2]),
      g = as.integer(str_match(param, "mu\\.alpha_stage\\[(\\d+),(\\d+)\\]")[,3])
    ) %>% left_join(match_s, by = c("s" = "stage_idx")) %>% left_join(match_g, by = c("g" = "group_idx"))
  ggplot(mu_alpha_stage, aes(x = mean, y = stage_3, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5))
  
  alpha_stage = summary(msom) %>% as.data.frame() %>% rownames_to_column("param") %>%
    filter(str_starts(param, "alpha_stage\\[")) %>%
    mutate(
      s = as.integer(str_match(param, "alpha_stage\\[(\\d+),(\\d+)\\]")[,2]),
      i = as.integer(str_match(param, "alpha_stage\\[(\\d+),(\\d+)\\]")[,3])
    ) %>% left_join(match_s, by = c("s" = "stage_idx")) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  ggplot(alpha_stage %>% filter(stage_3 == "standinit") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "standinit") + theme(legend.position = "bottom")

  ggplot(alpha_stage %>% filter(stage_3 == "compex") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "compex") + theme(legend.position = "bottom")
  
  ggplot(alpha_stage %>% filter(stage_3 == "mature") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "mature")
  
  mu_alpha_point = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_point")) %>%
    mutate(param = str_extract(parameter, "alpha_point\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "mu\\.alpha_point(\\d+)\\[(\\d+)\\]")[,2]),
      g = as.integer(str_match(parameter, "mu\\.alpha_point(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_alpha_point_data %>% select(param, name), by = c("param")) %>%
    left_join(match_g, by = c("g" = "group_idx"))
  ggplot(mu_alpha_point, aes(x = mean, y = name, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5))
  
  mu_alpha_homerange = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_homerange")) %>%
    mutate(param = str_extract(parameter, "alpha_homerange\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "mu\\.alpha_homerange(\\d+)\\[(\\d+)\\]")[,2]),
      g = as.integer(str_match(parameter, "mu\\.alpha_homerange(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_alpha_homerange_data %>% select(param, name), by = c("param")) %>%
    left_join(match_g, by = c("g" = "group_idx"))
  ggplot(mu_alpha_homerange, aes(x = mean, y = name, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5))
  
  mu_alpha_season = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_season")) %>%
    mutate(
      g = as.integer(str_match(parameter, "mu\\.alpha_season\\[(\\d+)\\]")[,2])
    ) %>% left_join(match_g, by = c("g" = "group_idx"))
  ggplot(mu_alpha_season, aes(x = mean, y = parameter, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5))
}


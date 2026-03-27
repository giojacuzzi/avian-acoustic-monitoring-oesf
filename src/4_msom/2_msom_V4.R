## 2_msom_mcmc.R #########################################################################################
# Fit a multispecies occupancy model using MCMC with stage as a continuous % variable, using compex as the baseline reference
#
## CONFIG:
grouping = "all" # Species grouping ("all", "diet", "nest", "nest_ps" ...)
model_type = "nofp" # nofp, fp
param_alpha_stage = "stratum_4"
param_alpha_season = "season"
param_alpha_point_names = c(
  "elevation",
  "dist_road_paved",
  "dist_watercourse_major"
)
param_alpha_plot_names = c( # stage-specific effects
  # "ba_mean",
  "tree_acre_6_mean",
  "qmd_6_mean",
  "bap_hwd_mean"
  # "ht_t100_mean"           # ~ Across most stages: qmd_6_mean, and to a slightly lesser extent tree_acre_6_mean
  # "ht_t100_cv",   # ~
  # "sdi_sum_mean", # ~
  # "canopy_cover_mean"    # ~ Across most stages: tree_acre_6_mean, and to a slightly lesser extent qmd_6_mean
  # "canopy_cover_cv",   # ~ 
  # "canopy_layers_mean"     # ~ Across most stages: tree_acre_6_mean, and to a slightly lesser extent qmd_6_mean
  # "cfvol_ddwm_mean"    # ~ qmd_6_mean
)
param_alpha_homerange_names = c(
  "pcnt_standinit",
  "pcnt_thin",
  "pcnt_mature"
  # TODO: patch area?
  # "density_edge_cw", # ~ pcnt_standinit
  # "aggregation_idx"  # ~ various pcnt covers depending on scale (very small [western flycatcher] or large [bald eagle])
  # "shape_idx"        # ~ various pcnt covers depending on scale
)
param_beta_point_names = c(
  "prcp_mm_day",
  "tmax_deg_c",
  "yday"
)
#
## OUTPUT:
dir_out = "data/cache/models"
#
## INPUT:
path_y                                   = "data/cache/4_msom/1_assemble_msom_data/y.rds"
path_occurrence_predictor_plot_data      = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_plot_data.rds"
path_occurrence_predictor_homerange_data = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_homerange_data.rds"
path_detection_predictor_data            = "data/cache/4_msom/1_assemble_msom_data/V3_detection_predictor_data.rds"
path_trait_data                          = "data/cache/2_traits/1_agg_traits/trait_data.csv"
###########################################################################################################

source("src/global.R")

if (model_type == "nofp") {
  model_file = "src/4_msom/msom_V4_nofp.txt"
} else if (model_type == "fp") {
  model_file = "src/4_msom/msom_V4_fp.txt"
} else {
  stop("ERROR: Unsupported model type")
}

model_name = tools::file_path_sans_ext(basename(model_file))
path_out = paste0(dir_out, "/V4_", model_name, "_", model_type, "_", grouping, ".rds")
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
stages  = occurrence_predictor_plot_data[[1]] %>% select(site, stratum_4) %>% filter(site %in% sites)

# Naive occupancy ---------------------------------
naive_occ <- apply(y, c(3, 4), function(mat) {
  site_max <- apply(mat, 1, function(x) { # mat is [sites x visits] for one season x species combo
    if (all(is.na(x))) NA else max(x, na.rm = TRUE)
  })
  detected  <- sum(site_max >= 1, na.rm = TRUE)
  surveyed  <- sum(!is.na(site_max))
  if (surveyed == 0) NA else detected / surveyed
})
print(round(naive_occ, 3))
(avg_occ = data.frame(
  species    = colnames(naive_occ),
  avg_naive_occ = round(colMeans(naive_occ, na.rm = TRUE),3),
  row.names  = NULL
))

# Naive occupancy by stage -----------------------
site_detected <- apply(y, c(1, 3, 4), function(x) {
  if (all(is.na(x))) NA else as.integer(max(x, na.rm = TRUE) >= 1)
})
library(reshape2)
det_df <- melt(site_detected, varnames = c("site", "season", "species"), value.name = "detected")
det_df <- merge(det_df, stages, by = "site")
occ_by_stage_species <- det_df |>
  filter(!is.na(detected)) |>
  group_by(stratum_4, species, season) |>
  summarise(naive_occ = mean(detected), .groups = "drop") |>
  group_by(stratum_4, species) |>
  summarise(avg_naive_occ = mean(naive_occ), .groups = "drop") |>
  split(~species)
# Access by name, e.g.:
occ_by_stage_species[["brown creeper"]]
occ_by_stage_species[["rufous hummingbird"]]

# Total number of detections per species:
n_detections = apply(y, 4, function(x) sum(x > 1, na.rm = TRUE))

# Total number of site detections per species:
(n_detected_sites_per_year = apply(y, c(3,4), function(mat) {
  sum(rowSums(mat > 0, na.rm = TRUE) > 0)
}))
n_detected_sites = colSums(n_detected_sites_per_year)

# TODO: Exclude extremely rare species?
species_to_include = species[n_detected_sites > 2] # species[n_detected_sites >= round(0.05 * length(sites))]
species_to_exclude = setdiff(species, species_to_include)
message("Excluding ", length(species_to_exclude), " species with insufficient observations: ")
print(species_to_exclude)
species = species_to_include
species_idx = match(species, dimnames(y)$species)
y = y[,,,species_idx, drop = FALSE]
dimnames(y)$species = species
str(y)

# Ensure plot sites match
stopifnot(all(sapply(occurrence_predictor_plot_data, function(df) {
  identical(sites, df$site)
})))
# Ensure homerange sites match
stopifnot(all(sapply(occurrence_predictor_homerange_data, function(season_list) {
  sapply(season_list, function(df) {
    identical(sites, df$site)
  })
})))
# Ensure detection sites match
stopifnot(all(unlist(
  lapply(detection_predictor_data, function(arr) {
    sapply(dimnames(arr)[[3]], function(season) {
      identical(sites, rownames(arr[, , season]))
    })
  })
)))

# Transpose plot data
vars_plot <- setdiff(names(occurrence_predictor_plot_data[[1]]), "site")
x_plot <- lapply(vars_plot, function(v) {
  mat <- sapply(seasons, function(y) {
    occurrence_predictor_plot_data[[y]][[v]]
  })
  rownames(mat) <- sites
  colnames(mat) <- seasons
  mat
})
names(x_plot) <- vars_plot
stopifnot(all(sapply(occurrence_predictor_plot_data, function(df) identical(df$site, sites)))) # check

# Impute any missing plot data from 2020
# x_plot[["tree_acre_6_mean"]] # before
for (var in names(x_plot)) {
  mat <- x_plot[[var]]  # extract the matrix
  if ("2020" %in% colnames(mat)) {
    for (year in setdiff(colnames(mat), "2020")) {
      na_rows <- is.nan(mat[, year])
      if (any(na_rows)) {
        mat[na_rows, year] <- mat[na_rows, "2020"]
      }
    }
  }
  x_plot[[var]] <- mat  # put it back in the list
}
# x_plot[["tree_acre_6_mean"]] # after

# Ensure no missing values for variables we are modeling
stopifnot(!any(sapply(param_alpha_plot_names, function(var) {
  any(is.na(x_plot[[var]]) | is.nan(x_plot[[var]]))
})))

# # Keep only homerange scales for scales being modeled
# occurrence_predictor_homerange_data = lapply(
#   occurrence_predictor_homerange_data,
#   function(t) {
#     t[names(t) %in% c(species, 'plot')]
#   }
# )

# Transpose homerange data
scales <- names(occurrence_predictor_homerange_data[[1]])  # e.g., "plot", "median", species
n_sites <- length(sites)
n_scales <- length(scales)
n_seasons <- length(seasons)
vars_homerange <- setdiff(names(occurrence_predictor_homerange_data[[1]][[2]]), c("site", "scale", "buffer_radius_m"))

# Initialize list of 3D arrays: x_homerange[["var"]][site, scale, season]
x_homerange <- lapply(vars_homerange, function(v) {
  array(NA_real_, dim = c(n_sites, n_scales, n_seasons), dimnames = list(sites, scales, seasons))
})
names(x_homerange) <- vars_homerange

# Fill arrays
for (season in seasons) {
  # season <- seasons[s]
  for (scale in scales) {
    df <- occurrence_predictor_homerange_data[[season]][[scale]]
    for (v in vars_homerange) {
      if (v %in% names(df)) {
        x_homerange[[v]][, scale, season] <- df[[v]]
      } else {
        warning("Variable ", v, " not found in scale ", scale, " for season ", season)
      }
    }
  }
}
# x_homerange[["pcnt_standinit"]]["aa014i", "american crow", "2021"] # access example
# Ensure no missing values remain
stopifnot(!any(sapply(x_homerange, function(x) any(is.na(x)))))

# Prepare occurrence plot scale covariate data ----------------------------------------------------------------------

# Categorical stage/stratum
stages = occurrence_predictor_plot_data[[1]] %>% mutate(stratum_4 = factor(stratum_4)) %>% select(stratum_4) %>% mutate(stage_idx = as.integer(stratum_4))
# NOTE: Redefine 'compex' as idx 1, i.e. baseline stratum
stages$stage_idx = as.numeric(
  factor(stages$stratum_4, levels = c("compex", "standinit", "mature", "thin"))
)

# Season
season_data_list = list((scale(1:length(seasons))))
param_season_data = tibble(
  param = c("alpha_season"),
  name = "season",
  data = c(season_data_list)
)

# Store alpha parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_alpha_point_data = tibble(param = paste0("alpha_point", 1:length(param_alpha_point_names)), name = param_alpha_point_names)
param_alpha_point_data = param_alpha_point_data %>% rowwise() %>%
  mutate(
    data = list(scale(as.numeric(unlist(x_plot[[name]])))) # flatten across all seasons and standardize [j,t]
  ) %>% ungroup()
n_alpha_point_params = nrow(param_alpha_point_data)

param_alpha_plot_data = tibble(param = paste0("alpha_plot", 1:length(param_alpha_plot_names)), name = param_alpha_plot_names)
param_alpha_plot_data = param_alpha_plot_data %>% rowwise() %>%
  mutate(
    data = list(scale(as.numeric(unlist(x_plot[[name]])))) # flatten across all seasons and standardize [j,t]
  ) %>% ungroup()
n_alpha_plot_params = nrow(param_alpha_plot_data)

# Prepare occurrence homerange scale covariate data ----------------------------------------------------------------------

# Enforce minimum homerange scale to plot scale
for (i in species) {
  buffer_radius_m = occurrence_predictor_homerange_data[[1]][[i]]$buffer_radius_m[1]
  if (buffer_radius_m < 100) {
    message("Enforcing minimum homerange scale for species ", i)
    for (v in vars_homerange) {
      for (t in seasons) {
        x_homerange[[v]][ , i, t] = x_homerange[[v]][ , 'plot', t]
      }
    }
  }
}

# Store delta parameter ID, variable name, and standardized data in species-specific matricies
alpha_homerange_data = setNames(vector('list', length(param_alpha_homerange_names)), param_alpha_homerange_names)
for (param in param_alpha_homerange_names) {
  message(param)
  
  species_alpha_data = setNames(vector('list', length(species)), species)
  
  for (i in species) {
    species_alpha_data[[i]] = scale(as.numeric(unlist(x_homerange[[param]][ , i, ])))
  }
  alpha_homerange_data[[param]] = species_alpha_data
}
param_alpha_homerange_data = tibble(param = paste0("alpha_homerange", 1:length(param_alpha_homerange_names)), name  = param_alpha_homerange_names)
param_alpha_homerange_data = param_alpha_homerange_data %>% rowwise() %>% mutate(data = list(alpha_homerange_data[[name]])) %>% ungroup()
n_alpha_homerange_params = nrow(param_alpha_homerange_data)

# TODO: Is the above correct?

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
  y[y == 2] = 1
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
  msom_data[[paste0("x_", param_alpha_point_data$param[a])]] <- matrix(param_alpha_point_data$data[[a]], nrow = J, ncol = Tseason)
}
for (a in seq_len(n_alpha_plot_params)) {
  msom_data[[paste0("x_", param_alpha_plot_data$param[a])]] <- matrix(param_alpha_plot_data$data[[a]], nrow = J, ncol = Tseason)
}
for (a in seq_len(n_alpha_homerange_params)) {
  tmp = do.call(cbind, lapply(param_alpha_homerange_data$data[[a]], as.vector))
  msom_data[[paste0("x_", param_alpha_homerange_data$param[a])]] <- array(
    tmp,
    dim = c(J, Tseason, I)
  )
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

# Check predictor correlation -------------------------------------------------------

# Do not retain any predictors with correlation coefficients > 0.8
threshold = 0.8
x_alpha_names = grep("^x_alpha", names(msom_data), value = TRUE)
x_alpha_names = setdiff(x_alpha_names, grep("^x_alpha_homerange", x_alpha_names, value = TRUE))
flat = lapply(msom_data[x_alpha_names], as.vector)
x_alpha = flat |> as.data.frame() |> as_tibble()
x_alpha$stage = rep(msom_data$stage, Tseason)

# Check for correlation among point- and plot-scale predictors
pairwise_collinearity(x_alpha, threshold)

# Check for correlation among homerange predictors across all species spatial scales
x_alpha_homerange_names = grep("^x_alpha_homerange", names(msom_data), value = TRUE)
for (i in 1:length(species)) {
  species_name = species[i]
  md_i = lapply(msom_data[x_alpha_homerange_names], function(x) x[, , i])
  flat = lapply(md_i, as.vector)
  x_alpha_homerange = flat |> as.data.frame() |> as_tibble()
  x_alpha_all = cbind(x_alpha, x_alpha_homerange)
  pc = pairwise_collinearity_by_group(x_alpha_all, "stage", threshold)
  if (nrow(pc)) {
    message(yellow("Correlation >", threshold, "for", species_name, "radius =", unique(occurrence_predictor_homerange_data[[1]][[species_name]]$buffer_radius_m)))
    print(pc)
  }
}

# Check among all predictors, using median homerange values
for (s in c('plot', 'median')) {
  for (t in 1:Tseason) {
    message("Multicollinearity at '", s, "' homerange scale season ", t)
    x_alpha_homerange = occurrence_predictor_homerange_data[[t]][[s]][param_alpha_homerange_names]
    x_alpha_all = cbind(x_alpha, x_alpha_homerange)
    
    # message("Pairwise collinearity for all predictors across all stages:")
    # print(pairwise_collinearity(x_alpha_all |> select(-stage), threshold))
    message("Pairwise collinearity by stage for all predictors:")
    print(pairwise_collinearity_by_group(x_alpha_all, "stage", threshold))
    
    # # VIF analysis for multicollinearity (consider dropping variable(s) with high VIF values (> 10))
    # message("VIF analysis for multicollinearity (all):")
    # model = lm(rep(1, nrow(x_alpha_all |> select(-stage))) ~ ., data = x_alpha_all |> select(-stage))
    # print(sort(vif(model)))
  }
}

# VIF at median homerange values PER STAGE flattening all seasons
scale = 'median'
mango = map2_dfr( # flatten data across all seasons
  occurrence_predictor_homerange_data, names(occurrence_predictor_homerange_data), 
  ~ .x[[scale]] %>% mutate(season = .y)
) %>% select(season, all_of(param_alpha_homerange_names))
x_alpha_all = cbind(x_alpha, mango)

for (s in 1:S) {
  message(s)
  x_alpha_all_stage = x_alpha_all %>% filter(stage == s) %>% as.data.frame() %>% select(-stage, -season)
  model <- lm(rep(1, nrow(x_alpha_all_stage)) ~ ., data = x_alpha_all_stage)
  print(vif(model))
}

# Initialize JAGS ----------------------------------------------------------------------------

# Initialize latent occupancy state z[i] as 1 if a detection occurred at site i, and 0 otherwise
z = array(NA, dim = c(length(sites), length(seasons), length(species)), dimnames = list(sites, seasons, species))
for (j in seq_along(sites)) {
  for (t in seq_along(seasons)) {
    for (i in seq_along(species)) {
      z[j, t, i] = (sum(y[j, , t, i], na.rm = TRUE) > 0) * 1
    }
  }
}

# Initial values to avoid data/model conflicts
init_vals = list(z = z)
# "This model has multimodal likelihood, resulting in identical support for different parameter values. We addressed this issue by constraining the parameters so that p11 > p10, an assumption that is supported by the validation  data at the chosen threshold17. We imposed this constraint by providing arbitrarily chosen starting values that  align with this condition (0.7 for p11 and 0.1 for p10)23." 
if (model_type == "fp") {
  # informative priors necessary to avoid invalid PPC log(0) values
  init_vals[["v"]] = matrix(qlogis(0.70), nrow = S, ncol = length(species))
  init_vals[["w"]] = rep(logit(0.05), length(species)) # NOTE: Ubiquitous species may have identifiability problems (wide `w` BCI) because there is very little z=0 data to inform w.
}

# Monitor parameter values
params_to_monitor = c(
  # Occupancy
  "mu.u", "sigma.u", "u",
  paste0("mu.alpha_point",     1:n_alpha_point_params),     paste0("sigma.alpha_point", 1:n_alpha_point_params),        paste0("alpha_point", 1:n_alpha_point_params),
  paste0("mu.alpha_plot",      1:n_alpha_plot_params),      paste0("sigma.alpha_plot", 1:n_alpha_plot_params),          paste0("alpha_plot", 1:n_alpha_plot_params),
  paste0("mu.alpha_homerange", 1:n_alpha_homerange_params), paste0("sigma.alpha_homerange", 1:n_alpha_homerange_params), paste0("alpha_homerange", 1:n_alpha_homerange_params),
  "mu.alpha_season", "sigma.alpha_season", "alpha_season",
  # True positive detection
  "mu.v", "sigma.v", "v",
  paste0("mu.beta_point",      1:n_beta_point_params),      paste0("sigma.beta_point",  1:n_beta_point_params),         paste0("beta_point",  1:n_beta_point_params),
  "mu.beta_point3_sq", "sigma.beta_point3_sq", "beta_point3_sq",
  # Others
  "D_obs", "D_sim", "z"
)
if (model_type == "fp") {
  params_to_monitor = c(params_to_monitor,
                        "epsilon", "sigma.epsilon", # False positive detection
                        "mu.w", "sigma.w", "w",
                        "mu.b", "sigma.b", "b")
}

message("Model file:", model_file)

message("\n", "System CPU: "); print(as.data.frame(t(benchmarkme::get_cpu())))
message("System RAM: "); print(benchmarkme::get_ram())

# Run JAGS ----------------------------------------------------------------------------------------------------

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { init_vals },
            parameters.to.save = params_to_monitor,
            model.file = model_file,
            # Longer test for lab:
            # n.chains = 3, n.adapt  = 2000, n.iter   = 5000, n.burnin = 1000, n.thin   = 1, parallel = TRUE,
            # Shorter test for home:
            n.chains = 3, n.adapt = 500, n.iter = 2500, n.burnin = 500, n.thin = 1, parallel = TRUE, # 12 hr for fp; 3 hr for nofp
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
  print(suspected_nonconvergence, n = Inf)
} else {
  message("All parameters appear to have converged (rhat < ", rhat_threshold, ")")
}

## Posterior predictive check - Bernoulli deviance contribution (Broms et al. 2016)

# Sanity check: confirm no Inf or NA values in observed or simulated data
summary(msom$sims.list$D_obs)
summary(msom$sims.list$D_sim)

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

model_data = read_rds(path_out) # DEBUG: "data/cache/models/Mar8_mscom_pcnt_nofp_all.rds"
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
# whiskerplot(msom, 'mu.alpha_stage')
whiskerplot(msom, c(paste0('mu.', param_alpha_plot_data$param)))
whiskerplot(msom, c(paste0('mu.', param_alpha_point_data$param)))
whiskerplot(msom, c(paste0('mu.', param_alpha_homerange_data$param)))
whiskerplot(msom, 'mu.alpha_season')

whiskerplot(msom, 'mu.v')
whiskerplot(msom, c(paste0('mu.', param_beta_point_data$param)))

{
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
  
  alpha_homerange = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "alpha_homerange")) %>%
    mutate(param = str_extract(parameter, "alpha_homerange\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "alpha_homerange(\\d+)\\[(\\d+)\\]")[,2]),
      i = as.integer(str_match(parameter, "alpha_homerange(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_alpha_homerange_data %>% select(param, name), by = c("param")) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  ggplot(alpha_homerange %>% filter(name == "pcnt_standinit") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "pcnt_standinit") + theme(legend.position = "bottom")
  ggplot(alpha_homerange %>% filter(name == "pcnt_mature") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "pcnt_mature") + theme(legend.position = "bottom")
  ggplot(alpha_homerange %>% filter(name == "pcnt_thin") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "pcnt_thin") + theme(legend.position = "bottom")
  
  mu_alpha_plot = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_plot")) %>%
    mutate(param = str_extract(parameter, "alpha_plot\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "mu\\.alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,2]),
      s = as.integer(str_match(parameter, "mu\\.alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,3]),
      g = as.integer(str_match(parameter, "mu\\.alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,4])
    ) %>% left_join(param_alpha_plot_data %>% select(param, name), by = c("param")) %>%
    left_join(match_s, by = c("s" = "stage_idx")) %>% left_join(match_g, by = c("g" = "group_idx"))
  ggplot(mu_alpha_plot, aes(x = mean, y = name, color = group, shape = stratum_4)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5))
  
  alpha_plot = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "alpha_plot")) %>%
    mutate(param = str_extract(parameter, "alpha_plot\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,2]),
      s = as.integer(str_match(parameter, "alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,3]),
      i = as.integer(str_match(parameter, "alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,4])
    ) %>% left_join(param_alpha_plot_data %>% select(param, name), by = c("param")) %>%
    left_join(match_s, by = c("s" = "stage_idx")) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  ggplot(alpha_plot %>% filter(name == "tree_acre_6_mean", stratum_4 == "mature") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group, shape = stratum_4)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "alpha_plot x stratum_4") + theme(legend.position = "bottom")
  
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
  
  mu_alpha_season = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_season")) %>%
    mutate(
      g = as.integer(str_match(parameter, "mu\\.alpha_season\\[(\\d+)\\]")[,2])
    ) %>% left_join(match_g, by = c("g" = "group_idx"))
  ggplot(mu_alpha_season, aes(x = mean, y = parameter, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5))
  
  alpha_season = summary(msom) %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "alpha_season")) %>%
    mutate(
      i = as.integer(str_match(parameter, "alpha_season\\[(\\d+)\\]")[,2])
    ) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  ggplot(alpha_season %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
    geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
    labs(title = "season") + theme(legend.position = "bottom")
}

path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"

source("src/global.R")
library(reshape2)

message("Loading site key from ", path_site_key)
site_key = read_csv(path_site_key, show_col_types = FALSE) %>% mutate(site = str_to_lower(site))

path_y = "data/cache/4_msom/1_assemble_msom_data/y.rds"
message("Loading species detection histories 'y' from ", path_y)
y = readRDS(path_y)


# INSPECT MSOM ===============================================================================

stop("READY FOR INSPECTION")

model_data = read_rds(path_msom) # DEBUG: "data/cache/models/Mar8_mscom_pcnt_nofp_all.rds"
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

# Compare naive site occupancy to estimated occupancy from MSOM ==================================

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
stages = data.frame(
  site = dimnames(y)$site,
  stage = site_key$stratum[ match(dimnames(y)$site, site_key$site) ]
)
det_df <- merge(det_df, stages, by = "site")
occ_by_stage_species <- det_df |>
  filter(!is.na(detected)) |>
  group_by(stage, species, season) |>
  summarise(naive_occ = mean(detected), .groups = "drop") |>
  group_by(stage, species) |>
  summarise(avg_naive_occ = mean(naive_occ), .groups = "drop") |>
  split(~species)
# Access by name, e.g.:
occ_by_stage_species[["brown creeper"]]
occ_by_stage_species[["rufous hummingbird"]]


# Load MSOM data

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
strata = as.factor(site_key$stratum[ match(sites, site_key$site) ])

# Posterior mean occupancy per site x season x species
z_mean <- apply(msom$sims.list$z, c(2, 3, 4), mean)
# Result: [224 sites x 4 seasons x 67 species]
# Assign dimnames to match y
dimnames(z_mean) <- list(
  site    = sites,
  season  = seasons,
  species = species
)

# Only compare species present in both
shared_species <- intersect(dimnames(y)[[4]], species)

# Subset z to shared species by name (not position)
z_mean_aligned <- z_mean[, , shared_species]

z_df <- melt(z_mean_aligned, varnames = c("site", "season", "species"), value.name = "z_mean")
z_df <- merge(z_df, stages, by = "site")

# Average by stratum x species (same averaging structure as before)
z_by_stage <- z_df |>
  group_by(stage, species, season) |>
  summarise(z_occ = mean(z_mean), .groups = "drop") |>
  group_by(stage, species) |>
  summarise(avg_z_occ = mean(z_occ), .groups = "drop") |>
  split(~species)

# Now compare directly
z_by_stage[["rufous hummingbird"]]
occ_by_stage_species[["rufous hummingbird"]]


z_by_stage[["barred owl"]]
occ_by_stage_species[["barred owl"]]

z_by_stage[["marbled murrelet"]]
occ_by_stage_species[["marbled murrelet"]]


z_by_stage[["pileated woodpecker"]]
occ_by_stage_species[["pileated woodpecker"]]


# =======================

# Pull mean z per species per stage from standard model
z_stage_species <- z_df |>
  group_by(stage, species) |>
  summarise(avg_z = mean(z_mean), .groups = "drop")

# Which species are most enriched in STAND INIT vs other stages?
z_wide <- z_stage_species |>
  pivot_wider(names_from = stage, values_from = avg_z) |>
  mutate(standinit_enrichment = `STAND INIT` - rowMeans(across(c(`COMP EXCL`, MATURE, THINNED)))) |>
  arrange(desc(standinit_enrichment))

head(z_wide, 15)

# Mature vs others?
z_wide %>%
  mutate(mature_enrichment = MATURE - rowMeans(across(c(`COMP EXCL`, `STAND INIT`, THINNED)))) %>%
  arrange(desc(mature_enrichment)) %>%
  head(15)

path_occurrence_predictor_plot_data      = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_plot_data.rds"
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)

# Look at differences by age
df_ordered <- occurrence_predictor_plot_data[[1]] %>%
  filter(site %in% sites) %>%
  select(site, age_mean) %>%
  slice(match(sites, site))
stopifnot(all(df_ordered$site == sites))
df_ordered = df_ordered %>% left_join(stages, by = c("site"))

library(tidyverse)
library(mgcv)

# --- 1. Average z_mean across seasons ---
# z_mean is [sites × seasons × species], average over seasons (dim 2)
z_avg <- apply(z_mean, c(1, 3), mean, na.rm = TRUE)  # now [224 × 67]

# --- 2. Convert to long format and join site metadata ---
z_df <- as.data.frame(z_avg) |>
  rownames_to_column("site") |>
  pivot_longer(-site, names_to = "species", values_to = "occupancy") |>
  left_join(df_ordered, by = "site")

# --- 3. Filter to MATURE only ---
mature_df <- z_df |>
  filter(stage == "MATURE")

# --- 4. Quick look at the age distribution within mature ---
mature_df |> 
  distinct(site, age_mean) |> 
  summary()

focal_spp <- c("brown creeper", "pileated woodpecker", "vaux's swift",
               "red-breasted nuthatch", "hairy woodpecker",
               "chestnut-backed chickadee", "barred owl")

mature_df |>
  filter(species %in% focal_spp) |>
  ggplot(aes(x = age_mean, y = occupancy)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),  
              # k=4 conservative given likely small n of mature sites
              color = "steelblue", se = TRUE) +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Mean stand age", y = "Occupancy probability",
       title = "Old-forest species occupancy within mature stands ~ age") +
  theme_minimal()





# --- The real question: at what age does the community 
#     transition away from comp excl character? ---

# Use all forest types, all focal species
# This is where your age gradient actually has variance to work with

transition_df <- z_df |>
  filter(species %in% focal_spp) |>
  left_join(df_ordered |> select(site, age_mean, stage), by = "site")

# Per-species GAM across all stages using continuous age
transition_models <- z_df |>
  group_by(species) |>
  nest() |>
  mutate(
    gam_fit = map(data, ~ gam(occupancy ~ s(age_mean, k = 5),
                              data = .x)),
    gam_tidy = map(gam_fit, ~ {
      s <- summary(.x)
      tibble(edf = s$edf, p.value = s$s.pv, r.sq = s$r.sq)
    })
  )

transition_models |>
  unnest(gam_tidy) |>
  select(species, edf, r.sq, p.value) |>
  arrange(desc(r.sq))

# Visualize transition across full age gradient
z_df |>
  filter(species %in% focal_spp) |>
  ggplot(aes(x = age_mean, y = occupancy, color = stage)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(aes(group = 1), method = "gam", 
              formula = y ~ s(x, k = 5),
              color = "black", se = TRUE) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Mean stand age", y = "Occupancy probability",
       title = "Occupancy across full age gradient — all forest types",
       color = "Stage") +
  theme_minimal()

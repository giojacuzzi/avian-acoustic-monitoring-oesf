## Calculate landscape matrix variables and test for their effects on richness, FDis, and old growth specialist presence
#
source("src/global.R")

t = 2020 # year: 2020, 2021, 2022, 2023
cover_classification = "clean_strata_4"

path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"

# Get response variables ====================================================

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)
message("Loading site key from ", path_site_key)
site_key = read_csv(path_site_key, show_col_types = FALSE) %>% mutate(site = str_to_lower(site))
message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
stages = model_data$stages

strata = as.factor(site_key$stratum[ match(sites, site_key$site) ])

species_traits = species_traits %>% filter(common_name %in% species)

# Get the expected occurrence for year t
z = msom$sims.list$z
z = z[ , , match(t , c(2020, 2021, 2022, 2023)), ]
z_mean <- apply(z, c(2, 3), mean)
rownames(z_mean) = sites
colnames(z_mean) = species

# Sum across species (dim 3) for each posterior draw and site
# Result: 6000 draws × 224 sites
richness_posterior <- apply(z, c(1, 2), sum)

# Summarize across draws (dim 1)
richness_mean  <- apply(richness_posterior, 2, mean)
richness_sd    <- apply(richness_posterior, 2, sd)
richness_lower <- apply(richness_posterior, 2, quantile, 0.025)
richness_upper <- apply(richness_posterior, 2, quantile, 0.975)

# Combine into a data frame
richness_df <- data.frame(
  site         = sites,
  mean         = richness_mean,
  sd           = richness_sd,
  lower_95     = richness_lower,
  upper_95     = richness_upper
)
stopifnot(richness_df$site == sites)
richness_df$strata = stages$stratum_4
SR = richness_df %>% select(mean, site, strata) %>% rename(SR = mean)

# Get FDis
# Create species-by-traits matrix
trait_matrix = species_traits %>% as.data.frame() %>%
  select(group_nest_ps,
         group_forage_substrate,
         group_diet,
         group_migrant,
         mass)
rownames(trait_matrix) = species_traits$common_name
# TODO: Subset species to those with at least one occurrence?

stopifnot(rownames(trait_matrix) == colnames(z_mean))  # must be TRUE

FDis = fdisp(gowdis(trait_matrix), z_mean)
FDis = stack(FDis$FDis)
colnames(FDis) = c("FDis", "site")
stopifnot(all(FDis$site == sites))
FDis$strata = stages$stratum_4

ggplot(richness_df, aes(x = strata, y = mean, fill = strata)) + geom_boxplot() +
  ggplot(FDis, aes(x = strata, y = FDis, fill = strata)) + geom_boxplot()

# Calculate landscape predictors ===========================================

path_rast_cover = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_", t, "_", cover_classification, ".tif")
path_site_cover_class = paste0("data/cache/3_gis/2_gen_cover_rasters/site_cover_class_", t, "_sf.rds")

sites_sf = read_rds(path_site_cover_class)

message('Loading raster cover data from cache ', path_rast_cover)
rast_cover = rast(path_rast_cover)

mapview(rast_cover) + mapview(sites_sf)

# For each site buffer:
# Fragment size — how much mature forest is nearby?
# Isolation — how far to nearest mature patch?
# Matrix proportion — what % is COMP EXCL?
# Edge density — interface between mature and matrix

library(landscapemetrics)
metrics <- sample_lsm(
  rast_cover,
  y = sites_sf,
  plot_id = sites_sf$site,
  shape = "circle",
  size = 3200,
  what = c(
    "lsm_c_pland",
    "lsm_c_ca",      # total class area (mature forest area in buffer)
    "lsm_c_np",      # number of patches
    "lsm_c_ed",      # edge density
    # "lsm_l_shdi",    # Shannon diversity index of landscape
    "lsm_c_enn_mn"   # mean Euclidean nearest neighbor distance (isolation)
    #     "lsm_p_area",    # area of focal patch (patch the point falls in)
    # "lsm_c_lpi",     # largest patch index (% buffer = largest patch)
    # "lsm_c_clumpy"   # clumpiness - how aggregated vs dispersed is mature forest
  )
)
levels(rast_cover)

# Create class lookup from your raster categories
class_lookup <- data.frame(
  class = levels(rast_cover)[[1]][,1],  # adjust based on levels() output
  stage = levels(rast_cover)[[1]][,2]  # adjust order
)

# Pivot to wide format: one row per site, one column per metric x stage
metrics_wide <- metrics %>%
  left_join(class_lookup, by = "class") %>%
  filter(stage != "other") %>%  # optionally drop non-forest
  pivot_wider(
    id_cols    = plot_id,
    names_from = c(metric, stage),
    values_from = value,
    names_sep  = "_"
  ) %>%
  rename(site = plot_id)

head(metrics_wide)

metrics_filtered <- metrics_wide %>%
  filter(site %in% sites) %>%
  slice(match(sites, site))
stopifnot(all(metrics_filtered$site == sites))


# --- 2. Distance to nearest mature patch ---
# landscapemetrics can extract mature patches as polygons directly

# Categorical raster — need to match on label
# Use numeric level 4 which corresponds to "mature"
mature_class_value <- 4  # numeric ID for mature from levels()

# classify works on the underlying numeric codes
mature_mask <- classify(
  rast_cover,
  rcl = matrix(c(4, 1), ncol = 2),
  others = NA
)

plot(mature_mask)  # sanity check — should show only mature patches

# Vectorize
mature_patches <- mature_mask |>
  as.polygons(dissolve = TRUE) |>
  st_as_sf() |>
  st_transform(st_crs(sites_sf))  # ensure CRS matches your sites

# Compute distance
dist_matrix <- st_distance(sites_sf, mature_patches)

dist_to_mature <- sites_sf |>
  st_drop_geometry() |>
  select(site) |>
  mutate(
    dist_to_mature_m = as.numeric(dist_matrix[, 1])
    # [,1] because st_union gives single geometry
    # if you skipped st_union use apply(dist_matrix, 1, min)
  )
dist_to_mature_filtered <- dist_to_mature %>%
  filter(site %in% sites) %>%
  slice(match(sites, site))
stopifnot(all(dist_to_mature_filtered$site == sites))

# stages is in the same row order as metrics_filtered
# bind site column directly
stages_with_sites <- stages |>
  mutate(site = metrics_filtered$site) |>
  rename(stage = stratum_4) |>
  select(site, stage, stage_idx)

# Quick check — do stages align with dist_to_mature?
# aa019i has dist=0 so should be mature
stages_with_sites |> filter(site %in% c("aa019i", "aa014i", "aa016i"))

# Build full landscape metrics
landscape_metrics <- stages_with_sites |>
  left_join(metrics_filtered, by = "site") |>
  left_join(dist_to_mature_filtered, by = "site") |>
  mutate(dist_to_mature_km = dist_to_mature_m / 1000)

# Sanity check
landscape_metrics |>
  group_by(stage) |>
  summarise(
    mean_prop_mature     = mean(pland_mature, na.rm = TRUE),
    mean_dist_mature_km  = mean(dist_to_mature_km, na.rm = TRUE),
    median_dist_mature_km = median(dist_to_mature_km, na.rm = TRUE),
    n = n()
  )

# Confirm by looking at the actual distributions
landscape_metrics |>
  filter(stage == "compex") |>
  summarise(
    dist_sd     = sd(dist_to_mature_km),
    dist_min    = min(dist_to_mature_km),
    dist_max    = max(dist_to_mature_km),
    pland_sd    = sd(pland_mature),
    pland_min   = min(pland_mature),
    pland_max   = max(pland_mature)
  )

# Visualize distributions
landscape_metrics |>
  filter(stage == "compex") |>
  pivot_longer(cols = c(dist_to_mature_km, pland_mature),
               names_to = "metric", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of landscape predictors in comp excl sites")


# Convert to long format
z_df <- z_mean |>
  as.data.frame() |>
  rownames_to_column("site") |>
  pivot_longer(-site, names_to = "species", values_to = "occupancy")

# Then join stage info from landscape_metrics
z_df <- z_df |>
  left_join(landscape_metrics, by = "site")

# Verify
glimpse(z_df)
z_df |> count(stage)  # should match your n per stage from earlier


focal_spp = species
# Join occupancy to landscape metrics, filter to comp excl
source_sink_df <- z_df |>
  filter(species %in% focal_spp) |>
  filter(stage == "compex")

# Check n and predictor ranges before modeling
source_sink_df |>
  filter(species == "brown creeper") |>  # just one species to check
  summarise(
    n               = n(),
    dist_range      = paste(round(min(dist_to_mature_km), 2), 
                            round(max(dist_to_mature_km), 2), sep = "–"),
    pland_range     = paste(round(min(pland_mature), 1),
                            round(max(pland_mature), 1), sep = "–"),
    occ_range       = paste(round(min(occupancy), 3),
                            round(max(occupancy), 3), sep = "–")
  )

# Models
source_sink_models <- source_sink_df |>
  group_by(species) |>
  nest() |>
  mutate(
    fit_dist  = map(data, ~ lm(occupancy ~ dist_to_mature_km, data = .x)),
    fit_prop  = map(data, ~ lm(occupancy ~ pland_mature, data = .x)),
    fit_both  = map(data, ~ lm(occupancy ~ dist_to_mature_km + pland_mature,
                               data = .x)),
    tidy_dist  = map(fit_dist, broom::tidy),
    tidy_prop  = map(fit_prop, broom::tidy),
    tidy_both  = map(fit_both, broom::tidy),
    glance_dist = map(fit_dist, broom::glance),
    glance_prop = map(fit_prop, broom::glance),
    glance_both = map(fit_both, broom::glance)
  )

# Results table — all predictors together for easy comparison
results_dist <- source_sink_models |>
  unnest(tidy_dist) |>
  filter(term == "dist_to_mature_km") |>
  select(species, estimate_dist = estimate, p_dist = p.value)

results_prop <- source_sink_models |>
  unnest(tidy_prop) |>
  filter(term == "pland_mature") |>
  select(species, estimate_prop = estimate, p_prop = p.value)

results_r2 <- source_sink_models |>
  mutate(
    r2_dist = map_dbl(glance_dist, ~ .x$r.squared),
    r2_prop = map_dbl(glance_prop, ~ .x$r.squared),
    r2_both = map_dbl(glance_both, ~ .x$r.squared)
  ) |>
  select(species, r2_dist, r2_prop, r2_both)

source_sink_results <- results_dist |>
  left_join(results_prop, by = "species") |>
  left_join(results_r2, by = "species") |>
  arrange(p_dist)

source_sink_results %>% print(n = Inf)

# What patches do your mature survey sites actually fall in?
# Get the patch-level area for each mature survey site

sites_with_patch <- sites_sf |>
  filter(site %in% (stages_with_sites |>
                      filter(stage == "mature") |>
                      pull(site))) |>
  st_join(
    mature_patches |> select(patch_id, patch_area_ha),
    join = st_within
  ) |>
  st_drop_geometry() |>
  select(site, patch_id, patch_area_ha)

# Some sites may not fall strictly within a polygon due to 
# raster/vector edge misalignment — use nearest patch for those
sites_missing_patch <- sites_with_patch |>
  filter(is.na(patch_id)) |>
  pull(site)

if (length(sites_missing_patch) > 0) {
  cat("Sites not within any patch:", length(sites_missing_patch), "\n")
  # Use nearest feature instead
  nearest_patch <- sites_sf |>
    filter(site %in% sites_missing_patch) |>
    st_join(mature_patches |> select(patch_id, patch_area_ha),
            join = st_nearest_feature) |>
    st_drop_geometry() |>
    select(site, patch_id, patch_area_ha)
  
  sites_with_patch <- sites_with_patch |>
    filter(!is.na(patch_id)) |>
    bind_rows(nearest_patch)
}

# Distribution of patch sizes your mature sites actually fall in
sites_with_patch |>
  mutate(size_class = case_when(
    patch_area_ha < 1    ~ "< 1 ha",
    patch_area_ha < 10   ~ "1–10 ha",
    patch_area_ha < 100  ~ "10–100 ha",
    patch_area_ha < 1000 ~ "100–1000 ha",
    TRUE                 ~ "> 1000 ha"
  )) |>
  count(size_class) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  arrange(size_class)

# Summary stats
sites_with_patch |>
  summarise(
    mean_patch_ha   = round(mean(patch_area_ha, na.rm = TRUE), 1),
    median_patch_ha = round(median(patch_area_ha, na.rm = TRUE), 1),
    min_patch_ha    = round(min(patch_area_ha, na.rm = TRUE), 1),
    max_patch_ha    = round(max(patch_area_ha, na.rm = TRUE), 1)
  )

# Revisit core area specifically for patches your sites are in
sampled_patch_ids <- sites_with_patch |> 
  pull(patch_id) |> 
  unique()

sampled_patches <- mature_patches |>
  filter(patch_id %in% sampled_patch_ids)

# Core area at depths appropriate for low-contrast matrix
edge_depths_lowcontrast <- c(30, 50, 100)

core_sampled <- map_dfr(edge_depths_lowcontrast, function(depth) {
  cores <- sampled_patches |>
    st_buffer(-depth) |>
    mutate(
      core_area_ha = pmax(as.numeric(st_area(geometry)) / 10000, 0),
      has_core     = core_area_ha > 0
    ) |>
    st_drop_geometry()
  
  tibble(
    edge_depth_m          = depth,
    n_patches             = nrow(cores),
    n_with_core           = sum(cores$has_core),
    pct_with_core         = round(mean(cores$has_core) * 100, 1),
    total_core_ha         = round(sum(cores$core_area_ha), 1),
    median_core_ha        = round(median(cores$core_area_ha[cores$has_core]), 1),
    pct_area_retained     = round(sum(cores$core_area_ha) / 
                                    sum(sampled_patches$patch_area_ha) * 100, 1)
  )
})

cat("=== CORE AREA FOR SAMPLED PATCHES ONLY ===\n")
print(core_sampled)

# Per-site: how much core area does each survey point's patch have?
# And is the survey point itself in the core or edge zone?
site_core_context <- sites_with_patch |>
  left_join(
    map_dfr(c(50, 100), function(depth) {
      sampled_patches |>
        st_buffer(-depth) |>
        mutate(
          core_area_ha = pmax(as.numeric(st_area(geometry)) / 10000, 0),
          edge_depth_m = depth
        ) |>
        st_drop_geometry() |>
        select(patch_id, edge_depth_m, core_area_ha)
    }),
    by = "patch_id"
  ) |>
  left_join(
    landscape_metrics |> select(site, dist_to_mature_edge_m),
    by = "site"
  ) |>
  mutate(
    in_core_50m  = dist_to_mature_edge_m > 50,
    in_core_100m = dist_to_mature_edge_m > 100
  )

# Summary: what fraction of mature survey points are 
# in genuine interior by different depth thresholds?
site_core_context |>
  filter(edge_depth_m == 50) |>
  summarise(
    pct_sites_in_core_50m  = round(mean(in_core_50m,  na.rm = TRUE) * 100, 1),
    pct_sites_in_core_100m = round(mean(in_core_100m, na.rm = TRUE) * 100, 1),
    median_dist_to_edge    = round(median(dist_to_mature_edge_m, na.rm = TRUE), 1),
    mean_dist_to_edge      = round(mean(dist_to_mature_edge_m, na.rm = TRUE), 1)
  )

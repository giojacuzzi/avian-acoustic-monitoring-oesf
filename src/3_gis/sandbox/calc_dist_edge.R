## Calculate distance to edge of survey sites
#
source("src/global.R")

t = 2020 # year: 2020, 2021, 2022, 2023
cover_classification = "clean_strata_4"

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

path_rast_cover = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_", t, "_", cover_classification, ".tif")
path_site_cover_class = paste0("data/cache/3_gis/2_gen_cover_rasters/site_cover_class_", t, "_sf.rds")

sites_sf = read_rds(path_site_cover_class)

message('Loading raster cover data from cache ', path_rast_cover)
rast_cover = rast(path_rast_cover)

# =============================================================================
# EDGE PENETRATION ANALYSIS
# Oregon Coast Range — mature forest fragments in comp excl matrix
# =============================================================================

library(tidyverse)
library(terra)
library(sf)
library(mgcv)
library(broom)

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

# =============================================================================
# 1. REBUILD CORE DATA OBJECTS
# Assumes in environment: rast_cover, sites_sf, stages, 
#                         metrics_filtered, dist_to_mature_filtered, z_mean
# =============================================================================

# --- Site metadata ---
stages_with_sites <- stages |>
  mutate(site = metrics_filtered$site) |>
  rename(stage = stratum_4) |>
  select(site, stage, stage_idx)

# --- Mature patches (individual polygons) ---
mature_patches <- rast_cover |>
  classify(rcl = matrix(c(4, 1), ncol = 2), others = NA) |>
  as.polygons(dissolve = TRUE) |>
  st_as_sf() |>
  st_transform(st_crs(sites_sf)) |>
  st_cast("POLYGON") |>
  mutate(
    patch_id      = row_number(),
    patch_area_ha = as.numeric(st_area(geometry)) / 10000
  )

cat("Mature patches:", nrow(mature_patches), "\n")
cat("Total mature area:", round(sum(mature_patches$patch_area_ha)), "ha\n")

mapview(mature_patches, zcol = "patch_area_ha") +
  mapview(sites_sf) + mapview(st_buffer(sites_sf, 250))

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

# --- Distance from each site to nearest mature forest edge ---
mature_boundary <- mature_patches |>
  st_union() |>
  st_boundary()

dist_to_edge <- sites_sf |>
  mutate(
    dist_to_mature_edge_m = as.numeric(
      st_distance(geometry, mature_boundary)
    )
  ) |>
  st_drop_geometry() |>
  select(site, dist_to_mature_edge_m)

# --- Check site count alignment ---
cat("\nSites in stages_with_sites:", nrow(stages_with_sites), "\n")
cat("Sites in dist_to_edge:", nrow(dist_to_edge), "\n")

extra_in_edge <- dist_to_edge |> 
  anti_join(stages_with_sites, by = "site")
cat("Sites in dist_to_edge but not stages:", nrow(extra_in_edge), "\n")

missing_edge <- stages_with_sites |> 
  anti_join(dist_to_edge, by = "site")
cat("Sites in stages but missing edge distance:", nrow(missing_edge), "\n")

# --- Build landscape_metrics ---
landscape_metrics <- stages_with_sites |>
  left_join(metrics_filtered, by = "site") |>
  left_join(dist_to_mature_filtered, by = "site") |>
  mutate(dist_to_mature_km = dist_to_mature_m / 1000) |>
  left_join(dist_to_edge, by = "site")

# --- Rebuild z_df (long format occupancy) ---
z_df <- z_mean |>
  as.data.frame() |>
  rownames_to_column("site") |>
  pivot_longer(-site, names_to = "species", values_to = "occupancy") |>
  left_join(landscape_metrics, by = "site")

cat("\nz_df rows:", nrow(z_df), "\n")
cat("NA in dist_to_mature_edge_m:", sum(is.na(z_df$dist_to_mature_edge_m)), "\n")

# =============================================================================
# 2. PATCH SIZE DISTRIBUTION
# =============================================================================

patch_size_summary <- mature_patches |>
  st_drop_geometry() |>
  mutate(size_class = case_when(
    patch_area_ha < 1    ~ "< 1 ha",
    patch_area_ha < 10   ~ "1–10 ha",
    patch_area_ha < 100  ~ "10–100 ha",
    patch_area_ha < 1000 ~ "100–1000 ha",
    TRUE                 ~ "> 1000 ha"
  )) |>
  group_by(size_class) |>
  summarise(
    n_patches    = n(),
    pct_patches  = round(n() / nrow(mature_patches) * 100, 1),
    total_ha     = round(sum(patch_area_ha), 1),
    pct_area     = round(sum(patch_area_ha) / sum(mature_patches$patch_area_ha) * 100, 1)
  ) |>
  arrange(size_class)

cat("\n=== PATCH SIZE DISTRIBUTION ===\n")
print(patch_size_summary)

ggplot(mature_patches, aes(x = patch_area_ha)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Patch area (ha, log scale)", y = "Count",
       title = "Mature patch size distribution",
       subtitle = paste0("n = ", nrow(mature_patches), " patches, ",
                         round(sum(mature_patches$patch_area_ha)), " ha total")) +
  theme_minimal()

# =============================================================================
# 3. CORE AREA AT MULTIPLE EDGE DEPTHS
# =============================================================================

edge_depths <- c(50, 100, 200)

core_summary <- map_dfr(edge_depths, function(depth) {
  cores <- mature_patches |>
    st_buffer(-depth) |>
    mutate(
      core_area_ha = pmax(as.numeric(st_area(geometry)) / 10000, 0),
      has_core     = core_area_ha > 0
    ) |>
    st_drop_geometry()
  
  tibble(
    edge_depth_m          = depth,
    n_patches_total       = nrow(cores),
    n_patches_with_core   = sum(cores$has_core),
    pct_patches_with_core = round(mean(cores$has_core) * 100, 1),
    total_core_ha         = round(sum(cores$core_area_ha), 1),
    pct_area_remaining    = round(sum(cores$core_area_ha) / 
                                    sum(mature_patches$patch_area_ha) * 100, 1),
    median_core_ha        = round(median(cores$core_area_ha[cores$has_core]), 2)
  )
})

cat("\n=== CORE AREA BY EDGE DEPTH ===\n")
print(core_summary)

# =============================================================================
# 4. SURVEY SITE DISTANCE TO EDGE — BY FOREST TYPE
# =============================================================================

edge_by_stage <- landscape_metrics |>
  group_by(stage) |>
  summarise(
    mean_dist_m   = round(mean(dist_to_mature_edge_m, na.rm = TRUE), 1),
    median_dist_m = round(median(dist_to_mature_edge_m, na.rm = TRUE), 1),
    min_dist_m    = round(min(dist_to_mature_edge_m, na.rm = TRUE), 1),
    max_dist_m    = round(max(dist_to_mature_edge_m, na.rm = TRUE), 1),
    pct_over_50m  = round(mean(dist_to_mature_edge_m > 50,  na.rm = TRUE) * 100, 1),
    pct_over_100m = round(mean(dist_to_mature_edge_m > 100, na.rm = TRUE) * 100, 1),
    pct_over_200m = round(mean(dist_to_mature_edge_m > 200, na.rm = TRUE) * 100, 1),
    n = n()
  )

cat("\n=== SURVEY SITE DISTANCE TO MATURE EDGE BY STAGE ===\n")
print(edge_by_stage)

# Distribution plot
landscape_metrics |>
  ggplot(aes(x = dist_to_mature_edge_m, fill = stage)) +
  geom_histogram(bins = 30, color = "white") +
  geom_vline(xintercept = c(50, 100, 200), 
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  annotate("text", x = c(50, 100, 200), y = Inf,
           label = c("50m", "100m", "200m"),
           vjust = 1.5, hjust = -0.1, size = 3, color = "grey30") +
  facet_wrap(~ stage, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Distance to nearest mature forest edge (m)", y = "Count",
       title = "Survey site distance to mature forest edge",
       subtitle = "Dashed lines = common edge depth thresholds") +
  theme_minimal() +
  theme(legend.position = "none")

# =============================================================================
# 5. OCCUPANCY ~ DISTANCE TO EDGE WITHIN MATURE STANDS
# Tests Hypothesis 1: edge effects inflate comp excl species in mature patches
# Prediction: old-forest species occupancy INCREASES with distance from edge
#             generalist species occupancy DECREASES with distance from edge
# =============================================================================

# Species to test — mix of old-forest associated and generalists
old_forest_spp <- c("brown creeper", "pileated woodpecker", 
                    "red-breasted nuthatch", "barred owl",
                    "vaux's swift", "hairy woodpecker")

generalist_spp <- c("american robin", "steller's jay", 
                    "song sparrow", "wilson's warbler")

focal_edge_spp <- c(old_forest_spp, generalist_spp)

# Models within mature only
edge_models_mature <- z_df |>
  filter(stage == "mature",
         species %in% focal_edge_spp) |>
  group_by(species) |>
  nest() |>
  mutate(
    lm_fit    = map(data, ~ lm(occupancy ~ dist_to_mature_edge_m, data = .x)),
    lm_tidy   = map(lm_fit, broom::tidy),
    lm_glance = map(lm_fit, broom::glance)
  )

edge_results_mature <- edge_models_mature |>
  unnest(lm_tidy) |>
  filter(term == "dist_to_mature_edge_m") |>
  left_join(
    edge_models_mature |>
      unnest(lm_glance) |>
      select(species, r.squared),
    by = "species"
  ) |>
  select(species, estimate, std.error, statistic, p.value, r.squared) |>
  mutate(
    guild = ifelse(species %in% old_forest_spp, 
                   "old-forest associated", "generalist"),
    direction = ifelse(estimate > 0, "increases with depth", 
                       "decreases with depth")
  ) |>
  arrange(p.value)

cat("\n=== OCCUPANCY ~ DISTANCE TO EDGE WITHIN MATURE STANDS ===\n")
print(edge_results_mature)

# Visualization
z_df |>
  filter(stage == "mature", species %in% focal_edge_spp) |>
  mutate(guild = ifelse(species %in% old_forest_spp,
                        "Old-forest associated", "Generalist")) |>
  ggplot(aes(x = dist_to_mature_edge_m, y = occupancy, color = guild)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  geom_vline(xintercept = 100, linetype = "dashed", 
             color = "grey50", linewidth = 0.5) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Distance to nearest mature forest edge (m)",
       y = "Occupancy probability",
       color = "Guild",
       title = "Edge penetration within mature stands",
       subtitle = "Dashed line = 100m edge depth threshold") +
  theme_minimal() +
  theme(legend.position = "bottom")

# =============================================================================
# 6. RUN ACROSS ALL 67 SPECIES IN MATURE
# Identify which species show significant edge effects
# =============================================================================

edge_models_all <- z_df |>
  filter(stage == "mature") |>
  group_by(species) |>
  nest() |>
  mutate(
    # check variance first — ceiling species will be uninformative
    occ_sd  = map_dbl(data, ~ sd(.x$occupancy, na.rm = TRUE)),
    lm_fit  = map(data, ~ lm(occupancy ~ dist_to_mature_edge_m, data = .x)),
    lm_tidy = map(lm_fit, broom::tidy),
    r2      = map_dbl(lm_fit, ~ summary(.x)$r.squared)
  )

# Fix for the forest plot — occ_sd join issue
edge_plot_data <- edge_models_all |>
  unnest(lm_tidy) |>
  filter(term == "dist_to_mature_edge_m") |>
  select(species, estimate, std.error, p.value, r2, occ_sd) |>
  filter(occ_sd >= 0.02) |>
  mutate(sig = p.value < 0.05)

edge_plot_data |>
  ggplot(aes(x = estimate, y = reorder(species, estimate),
             color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error
  ), width = 0.3, orientation = "y") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("grey60", "firebrick"),
                     labels = c("p ≥ 0.05", "p < 0.05")) +
  labs(x = "Effect of distance to edge on occupancy (per meter)",
       y = NULL, color = NULL,
       title = "Edge penetration effects within mature stands",
       subtitle = "Positive = occupancy increases toward interior | Ceiling species excluded") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 7))


# Check: how many of your mature survey sites fall within 
# the 4 large patches vs all other patches?

large_patches <- mature_patches |>
  filter(patch_area_ha > 1000)

cat("Large patches (>1000 ha):", nrow(large_patches), "\n")
cat("Total area:", round(sum(large_patches$patch_area_ha)), "ha\n")

# Which mature survey sites are in large patches?
sites_in_large <- sites_sf |>
  filter(site %in% (stages_with_sites |> 
                      filter(stage == "mature") |> 
                      pull(site))) |>
  st_join(large_patches |> select(patch_id, patch_area_ha),
          join = st_within) |>
  mutate(in_large_patch = !is.na(patch_id)) |>
  st_drop_geometry()

sites_in_large |>
  count(in_large_patch) |>
  mutate(pct = round(n / sum(n) * 100, 1))

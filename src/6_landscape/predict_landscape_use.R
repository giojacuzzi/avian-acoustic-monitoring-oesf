####################################################################################
# Predict habitat use across the landscape
#
## CONFIG:
sp_name  = "olive-sided flycatcher"
season_t = "2020"
units_to_predict = c("Upper Clearwater", "Willy - Huel", "Kalaloch", "Queets", "Copper Mine")
mask_riparian_zones = TRUE # Exclude "Riparian Management Zones" interior-core buffers for type 1-4 streams
#
## INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
####################################################################################

source("src/global.R")

# Load MSOM ---------------------------------------------------------------------------
message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)
msom_summary = model_data$msom_summary
msom         = model_data$msom
groups       = model_data$groups %>% arrange(common_name)
sites        = model_data$sites
species      = model_data$species
seasons      = model_data$seasons
stages       = model_data$stages

z         = msom$sims.list$z
n_iter    = dim(z)[1]
n_sites   = dim(z)[2]
n_seasons = dim(z)[3]
n_species = dim(z)[4]

# Load plot and homerange scale data for all points -----------------------------------
path_rast_cover = "data/cache/3_gis/2_gen_cover_rasters/rast_cover_2020_clean_strata_4.tif"
message("Loading raster cover data from cache ", path_rast_cover)
rast_cover = rast(path_rast_cover)

path_homerange = "data/cache/7_landscape/OPT_landscape_data_homerange_scale_2020_clean_strata_4_landscape.rds"
path_point     = "data/cache/7_landscape/OPT_landscape_data_plot_scale_2020_clean_strata_4_landscape.rds"

data_homerange = read_rds(path_homerange)
points_point   = read_rds(path_point)

points_plot = st_as_sf(
  data_homerange[["plot"]] %>% select("x", "y", "bap_hwd_mean", "qmd_6_mean", "tree_acre_6_mean"),
  coords = c("x", "y"), crs = crs_m
)

points_homerange = st_as_sf(
  data_homerange[[sp_name]] %>% select("x", "y", "pcnt_standinit", "pcnt_compex", "pcnt_thin", "pcnt_mature"),
  coords = c("x", "y"), crs = crs_m
)

points = points_homerange %>%
  st_join(points_point, join = st_equals) %>%
  st_join(points_plot, join = st_equals)

# TODO: Impute NAs from nearby values, or leave NAs as they are?
if (anyNA(points)) {
  message(crayon::yellow("WARNING: Imputing NA values from layer means"))
  print(colSums(is.na(st_drop_geometry(points))))
  points = points %>%
    mutate(across(
      where(~ is.numeric(.) && anyNA(.)),
      ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
    ))
}

# Get stage from raster cover
stage_vals   = extract(rast_cover, vect(points))
points$stage = stage_vals$stage
points = points %>%
  left_join(
    stages %>%
      distinct() %>%
      mutate(stratum_4 = as.character(stratum_4)) %>%
      select(stratum_4, stage_idx),
    by = c("stage" = "stratum_4")
  )

# Build raster stack ------------------------------------------------------------------
point_data = cbind(st_coordinates(points), st_drop_geometry(points))

vars = c(
  "pcnt_standinit", "pcnt_compex", "pcnt_mature", "pcnt_thin",  # homerange
  "elevation", "dist_road_paved", "dist_watercourse_major",      # point
  "tree_acre_6_mean", "qmd_6_mean", "bap_hwd_mean",              # plot
  "stage_idx"
)

rasters = lapply(vars, function(v) {
  message(v)
  rast(point_data[, c("X", "Y", v)], type = "xyz", crs = crs_m_rast)
})
names(rasters) = vars
rast_stack = rast(rasters)

# Mask by management unit and forest-only upland extent -------------------------------------------------------
mask_poly = landscape_planning_units %>%
  filter(unit %in% units_to_predict) %>%
  st_transform(crs = crs(rast_stack)) %>%
  st_collection_extract("POLYGON") %>%  # extract only polygon parts
  st_union() %>%                         # dissolve into single geometry
  st_sf()                                # back to sf object
mask_poly = st_transform(mask_poly, crs = crs(rast_stack))

# Mask and crop
rast_stack = rast_stack %>% crop(vect(mask_poly)) %>% mask(vect(mask_poly))

if (mask_riparian_zones) {
  source("src/3_gis/1_preprocess_gis_data.R")
  
  # Interior-core buffers (per WADNR FEIS)
  riparian_buffers = watercourses %>%
    filter(sl_wtrty_c %in% c(1, 2, 3, 4)) %>%
    mutate(buffer_m = case_when(
      # 150 feet for type 1 and 2 streams
      sl_wtrty_c %in% c(1, 2) ~ 150 * conv_ft_to_m,
      # 100 feet for type 3 and 4 streams
      sl_wtrty_c %in% c(3, 4) ~ 100 * conv_ft_to_m
    )) %>%
    st_buffer(dist = .$buffer_m) %>%
    st_union() %>%  # dissolve overlapping buffers into single polygon
    st_sf()
  
  # mapview(watercourses %>% filter(sl_wtrty_c %in% c(1, 2, 3, 4))) + mapview(stream_buffers)
  rast_stack = rast_stack %>% mask(vect(stream_buffers), inverse = TRUE)
  names(rast_stack) = vars
}

# Predict occupancy probability -------------------------------------------------------

# Helper: apply training scale parameters to new data
sc = function(x, scaled_vec) {
  (x - attr(scaled_vec, "scaled:center")) / attr(scaled_vec, "scaled:scale")
}

sp_idx = which(species == sp_name)
t_idx  = which(seasons == season_t)

# Unpack param data
pad    = model_data$param_alpha_data
apoint = pad$param_alpha_point_data     # alpha_point1/2/3 → elevation/road/water
aplot  = pad$param_alpha_plot_data      # alpha_plot1/2/3  → tree_acre/qmd/bap_hwd
ahome  = pad$param_alpha_homerange_data # alpha_homerange1/2/3 → standinit/thin/mature

# Build prediction data frame and scale predictors
df_pred = as.data.frame(rast_stack, xy = TRUE) %>% na.omit()

# Point-scale (shared scaling across all species)
df_pred$elev_sc  = sc(df_pred$elevation,             apoint$data[[1]])
df_pred$road_sc  = sc(df_pred$dist_road_paved,       apoint$data[[2]])
df_pred$water_sc = sc(df_pred$dist_watercourse_major, apoint$data[[3]])

# Plot-scale (shared scaling)
df_pred$tree_acre_sc = sc(df_pred$tree_acre_6_mean, aplot$data[[1]])
df_pred$qmd_sc       = sc(df_pred$qmd_6_mean,       aplot$data[[2]])
df_pred$bap_hwd_sc   = sc(df_pred$bap_hwd_mean,     aplot$data[[3]])

# Homerange-scale (species-specific scaling)
df_pred$standinit_sc = sc(df_pred$pcnt_standinit, ahome$data[[1]][[sp_name]])
df_pred$thin_sc      = sc(df_pred$pcnt_thin,       ahome$data[[2]][[sp_name]])
df_pred$mature_sc    = sc(df_pred$pcnt_mature,     ahome$data[[3]][[sp_name]])

# Season scalar (reproduce training scaling exactly)
x_season_t = as.vector(scale(seq_along(seasons)))[t_idx]

# Extract posterior samples for this species
sl     = msom$sims.list
n_iter  = dim(sl$u)[1]
n_cells = nrow(df_pred)
sv      = df_pred$stage_idx  # stage index vector [n_cells]

u_i   = sl$u[, sp_idx]
ap1   = sl$alpha_plot1[, , sp_idx]  # [n_iter, S]
ap2   = sl$alpha_plot2[, , sp_idx]
ap3   = sl$alpha_plot3[, , sp_idx]
apt1  = sl$alpha_point1[, sp_idx]
apt2  = sl$alpha_point2[, sp_idx]
apt3  = sl$alpha_point3[, sp_idx]
ahr1  = sl$alpha_homerange1[, sp_idx]
ahr2  = sl$alpha_homerange2[, sp_idx]
ahr3  = sl$alpha_homerange3[, sp_idx]
aseas = sl$alpha_season[, sp_idx]

# Chunk loop: compute psi summaries + accumulate area samples in one pass
cell_area_ha = prod(res(rast_stack)) / 10000
stage_levels = c("compex", "standinit", "mature", "thin")  # compex=1, standinit=2, mature=3, thin=4

chunk_size = 1000
n_chunks   = ceiling(n_cells / chunk_size)

psi_mean = numeric(n_cells)
psi_sd   = numeric(n_cells)
psi_q025 = numeric(n_cells)
psi_q975 = numeric(n_cells)

expected_ha_samples  = numeric(n_iter)
expected_ha_by_stage = matrix(0, nrow = n_iter, ncol = length(stage_levels),
                               dimnames = list(NULL, stage_levels))

for (ch in seq_len(n_chunks)) {
  idx   = ((ch - 1) * chunk_size + 1):min(ch * chunk_size, n_cells)
  sv_ch = sv[idx]
  
  logit_psi_ch =
    outer(u_i, rep(1, length(idx))) +
    ap1[, sv_ch, drop = FALSE] * outer(rep(1, n_iter), df_pred$tree_acre_sc[idx]) +
    ap2[, sv_ch, drop = FALSE] * outer(rep(1, n_iter), df_pred$qmd_sc[idx]) +
    ap3[, sv_ch, drop = FALSE] * outer(rep(1, n_iter), df_pred$bap_hwd_sc[idx]) +
    outer(apt1, df_pred$elev_sc[idx]) +
    outer(apt2, df_pred$road_sc[idx]) +
    outer(apt3, df_pred$water_sc[idx]) +
    outer(ahr1, df_pred$standinit_sc[idx]) +
    outer(ahr2, df_pred$thin_sc[idx]) +
    outer(ahr3, df_pred$mature_sc[idx]) +
    outer(aseas, rep(x_season_t, length(idx)))
  
  psi_ch = plogis(logit_psi_ch)  # [n_iter x chunk_size]
  
  psi_mean[idx] = colMeans(psi_ch)
  psi_sd[idx]   = apply(psi_ch, 2, sd)
  psi_q025[idx] = apply(psi_ch, 2, quantile, 0.025)
  psi_q975[idx] = apply(psi_ch, 2, quantile, 0.975)
  
  # Accumulate area across chunks — rowSums gives [n_iter] per-iteration totals
  expected_ha_samples = expected_ha_samples + rowSums(psi_ch) * cell_area_ha
  
  for (s in seq_along(stage_levels)) {
    mask = sv_ch == s
    if (any(mask)) {
      expected_ha_by_stage[, s] = expected_ha_by_stage[, s] +
        rowSums(psi_ch[, mask, drop = FALSE]) * cell_area_ha
    }
  }
  
  if (ch %% 10 == 0) message("Chunk ", ch, " / ", n_chunks)
}

df_pred$psi_mean = psi_mean
df_pred$psi_sd   = psi_sd
df_pred$psi_q025 = psi_q025
df_pred$psi_q975 = psi_q975

# Back to raster
rast_psi = rast(
  df_pred[, c("x", "y", "psi_mean", "psi_sd", "psi_q025", "psi_q975")],
  type = "xyz", crs = crs(rast_stack)
)

# Visualize species distribution maps
stop("READY FOR INSPECTION")
# mapview(rast_psi[["psi_sd"]]) + mapview(rast_psi[["psi_mean"]])

# Predicted distribution map expressed in terms of model estimated probability of the species using an X by X unit, where X is the resolution of a cell
p_mean = ggplot() +
  geom_spatraster(data = rast_psi[["psi_mean"]]) +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1), na.value = NA) +
  labs(title = sp_name, fill = "Ψ") +
  theme(legend.position = "right")

# Corresponding map of standard deviation
p_sd = ggplot() +
  geom_spatraster(data = rast_psi[["psi_sd"]]) +
  scale_fill_viridis_c(option = "magma", na.value = NA) +
  labs(title = "", fill = "SD") +
  theme(legend.position = "right",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()
        )

p_mean + p_sd

# Posterior mean use probability distribution
ggplot(df_pred, aes(x = psi_mean, fill = after_stat(x))) +
  geom_histogram(bins = 100, color = NA) +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Posterior mean use probability (Ψ)", y = "Number of cells", title = sp_name) +
  theme(legend.position = "none")

# Posterior credible interval on total expected area --------------------
total_landscape_ha = nrow(df_pred) * cell_area_ha

area_summary = tibble(
  mean_ha            = mean(expected_ha_samples),
  sd_ha              = sd(expected_ha_samples),
  q025_ha            = quantile(expected_ha_samples, 0.025),
  median_ha          = quantile(expected_ha_samples, 0.500),
  q975_ha            = quantile(expected_ha_samples, 0.975),
  total_landscape_ha = total_landscape_ha,
  pct_landscape      = mean(expected_ha_samples) / total_landscape_ha * 100
)

message(sp_name, " — expected area used:")
print(area_summary)

# Posterior distribution of expected area used
ggplot(tibble(ha = expected_ha_samples), aes(x = ha)) +
  geom_histogram(bins = 60, fill = "gray", color = NA, alpha = 0.8) +
  geom_vline(xintercept = quantile(expected_ha_samples, c(0.025, 0.975)), linetype = "dashed") +
  geom_vline(xintercept = mean(expected_ha_samples), linetype = "solid") +
  labs(x = "Expected area used (ha)", y = "Posterior iterations",
    subtitle = paste0(
      "Mean: ", round(mean(expected_ha_samples)), " ha",
      " (", round(quantile(expected_ha_samples, 0.025)), "–",
      round(quantile(expected_ha_samples, 0.975)), " 95% BCI)"
    )
  )

# Breakdown by stage ----------------------------------------------------
stage_total_ha = tibble(
  stratum   = stage_levels,
  stage_idx = seq_along(stage_levels),
  total_ha  = as.numeric(table(factor(sv, levels = seq_along(stage_levels)))) * cell_area_ha
)

stage_summary = as.data.frame(expected_ha_by_stage) %>%
  pivot_longer(everything(), names_to = "stratum", values_to = "ha") %>%
  group_by(stratum) %>%
  summarise(
    mean_ha   = mean(ha),
    q025_ha   = quantile(ha, 0.025),
    median_ha = quantile(ha, 0.500),
    q975_ha   = quantile(ha, 0.975),
    .groups   = "drop"
  ) %>%
  left_join(stage_total_ha, by = "stratum") %>%
  mutate(
    pct_of_stage = mean_ha / total_ha * 100,
    stratum      = factor(stratum, levels = c("standinit", "compex", "thin", "mature"))
  )

message(sp_name, " — expected area occupied by stage:")
print(stage_summary)

# Gray is the total stage area
# Color is the expected used area with 95% BCI errorbars
ggplot(stage_summary, aes(x = stratum)) +
  geom_col(aes(y = total_ha), fill = "gray80") +
  geom_col(aes(y = mean_ha, fill = stratum), alpha = 0.9) +
  geom_errorbar(aes(ymin = q025_ha, ymax = q975_ha), width = 0.2, color = "gray10") +
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = unname(stage_colors)) +
  labs(x = NULL, y = "Area (ha)") +
  theme(legend.position = "none")

####################################################################################
### Predict habitat use across the landscape — ALL SPECIES
### CONFIG:
season_t            = "2020"
units_to_predict    = c("Upper Clearwater", "Willy - Huel", "Kalaloch", "Queets", "Copper Mine")
mask_riparian_zones = TRUE
richness_thin_k     = 10   # keep every k-th iteration for richness SD/quantiles
# memory ~ n_cells × (n_iter/k) × 8 bytes
# e.g. 1M cells × 400 thinned iters = ~3.2 GB
path_out_species_results = "data/cache/6_landscape/species_results.rds"
path_out_richness_results = "data/cache/6_landscape/richness_results.rds"

### INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
####################################################################################

source("src/global.R")

# Load MSOM ---------------------------------------------------------------------------
message("Loading data for multi-species occupancy model ", path_msom)
model_data   = readRDS(path_msom)
msom_summary = model_data$msom_summary
msom         = model_data$msom
groups       = model_data$groups %>% arrange(common_name)
sites        = model_data$sites
species      = model_data$species
seasons      = model_data$seasons
stages       = model_data$stages

sl     = msom$sims.list
n_iter = dim(sl$z)[1]

# DEBUG: Only predict conservation priority species for now
sp_to_predict = species
sp_to_predict = intersect(species, conpri_species)
# DEBUG

# Load raster cover -------------------------------------------------------------------
path_rast_cover = "data/cache/3_gis/2_gen_cover_rasters/rast_cover_2020_clean_strata_4.tif"
message("Loading raster cover data from cache ", path_rast_cover)
rast_cover = rast(path_rast_cover)

# Load trait data
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% mutate(species = common_name)

# Load landscape data -----------------------------------------------------------------
path_homerange = "data/cache/7_landscape/OPT_landscape_data_homerange_scale_2020_clean_strata_4_landscape.rds"
path_point     = "data/cache/7_landscape/OPT_landscape_data_plot_scale_2020_clean_strata_4_landscape.rds"
data_homerange = read_rds(path_homerange)
points_point   = read_rds(path_point)

# Coerce point-scale data to sf once (carries elevation / road / water / age_point) ---
points_point_sf = st_as_sf(
  points_point,
  coords = c("x", "y"),
  crs    = crs_m
)

# Plot-scale sf (same for all species) ------------------------------------------------
points_plot_sf = st_as_sf(
  data_homerange[["plot"]] %>% select("x", "y", "bap_hwd_mean", "qmd_6_mean", "tree_acre_6_mean"),
  coords = c("x", "y"),
  crs    = crs_m
)

# Shared scaling parameters -----------------------------------------------------------
pad    = model_data$param_alpha_data
apoint = pad$param_alpha_point_data
aplot  = pad$param_alpha_plot_data
ahome  = pad$param_alpha_homerange_data

# Season scalar -----------------------------------------------------------------------
x_season_t = as.vector(scale(seq_along(seasons)))[which(seasons == season_t)]

# Helper: apply training scale to new data --------------------------------------------
sc = function(x, scaled_vec) {
  (x - attr(scaled_vec, "scaled:center")) / attr(scaled_vec, "scaled:scale")
}

# Stage levels ------------------------------------------------------------------------
stage_levels = c("compex", "standinit", "mature", "thin")

# Build shared mask polygon (once) ----------------------------------------------------
mask_poly_raw = landscape_planning_units %>%
  filter(unit %in% units_to_predict) %>%
  st_collection_extract("POLYGON") %>%
  st_union() %>%
  st_sf()

# Pre-build riparian buffers (once — expensive) ---------------------------------------
if (mask_riparian_zones) {
  source("src/3_gis/1_preprocess_gis_data.R")
  riparian_buffers = watercourses %>%
    filter(sl_wtrty_c %in% c(1, 2, 3, 4)) %>%
    mutate(buffer_m = case_when(
      sl_wtrty_c %in% c(1, 2) ~ 150 * conv_ft_to_m,
      sl_wtrty_c %in% c(3, 4) ~ 100 * conv_ft_to_m
    )) %>%
    st_buffer(dist = .$buffer_m) %>%
    st_union() %>%
    st_sf()
}

vars = c(
  "pcnt_standinit", "pcnt_compex", "pcnt_mature", "pcnt_thin",  # homerange
  "elevation", "dist_road_paved", "dist_watercourse_major",      # point
  "tree_acre_6_mean", "qmd_6_mean", "bap_hwd_mean",              # plot
  "stage_idx"
)

# Thinned iteration indices for richness uncertainty ----------------------------------
thin_idx = seq(1, n_iter, by = richness_thin_k)
n_thin   = length(thin_idx)
message("Richness uncertainty will use ", n_thin, " thinned draws (every ", richness_thin_k, "th of ", n_iter, ")")

# =====================================================================================
# MAIN LOOP
# =====================================================================================

species_results = vector("list", length(sp_to_predict))
names(species_results) = sp_to_predict

# Richness accumulators — initialized after first species when n_cells is known
richness_sum      = NULL   # [n_cells] exact: ∑_{iter,sp} psi(iter,cell), for posterior mean
richness_thin_mat = NULL   # [n_thin × n_cells] thinned: ∑_sp psi_sp for thinned iters, for SD/quantiles
richness_xy       = NULL

i = 1
for (sp_name in sp_to_predict) {
  
  sp_num = which(species == sp_name)
  message("\n========== ", sp_name, " (", i, "/", length(sp_to_predict), ") ==========")
  
  if (species_traits %>% filter(common_name == sp_name) %>% pull(home_range_radius_m) < 100) {
    sp_scale = "plot"
  } else {
    sp_scale = sp_name
  }
  
  # --- Build points for this species -------------------------------------------------
  points_homerange = st_as_sf(
    data_homerange[[sp_scale]] %>%
      select("x", "y", "pcnt_standinit", "pcnt_compex", "pcnt_thin", "pcnt_mature"),
    coords = c("x", "y"),
    crs    = crs_m
  )
  
  points = points_homerange %>%
    st_join(points_point_sf, join = st_equals) %>%   # elevation / road / water / age_point
    st_join(points_plot_sf,  join = st_equals)        # tree_acre / qmd / bap_hwd
  
  if (anyNA(st_drop_geometry(points))) {
    message(crayon::yellow("WARNING: Imputing NA values from layer means for", sp_name))
    print(colSums(is.na(st_drop_geometry(points))))
    points = points %>%
      mutate(across(
        where(~ is.numeric(.) && anyNA(.)),
        ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
      ))
  }
  
  # Attach forest stage
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
  
  # --- Build raster stack ------------------------------------------------------------
  point_data = cbind(st_coordinates(points), st_drop_geometry(points))
  
  missing_vars = setdiff(vars, colnames(point_data))
  if (length(missing_vars) > 0) {
    stop("Missing columns after join for '", sp_name, "': ",
         paste(missing_vars, collapse = ", "))
  }
  
  rast_stack = rast(lapply(setNames(vars, vars), function(v) {
    rast(point_data[, c("X", "Y", v)], type = "xyz", crs = crs_m_rast)
  }))
  
  mask_poly  = st_transform(mask_poly_raw, crs = crs(rast_stack))
  rast_stack = rast_stack %>% crop(vect(mask_poly)) %>% mask(vect(mask_poly))
  
  if (mask_riparian_zones) {
    rast_stack   = rast_stack %>%
      mask(vect(st_transform(riparian_buffers, crs = crs(rast_stack))), inverse = TRUE)
    names(rast_stack) = vars
  }
  
  # --- Scale predictors --------------------------------------------------------------
  df_pred = as.data.frame(rast_stack, xy = TRUE) %>% na.omit()
  
  df_pred$elev_sc      = sc(df_pred$elevation,              apoint$data[[1]])
  df_pred$road_sc      = sc(df_pred$dist_road_paved,        apoint$data[[2]])
  df_pred$water_sc     = sc(df_pred$dist_watercourse_major, apoint$data[[3]])
  df_pred$tree_acre_sc = sc(df_pred$tree_acre_6_mean,       aplot$data[[1]])
  df_pred$qmd_sc       = sc(df_pred$qmd_6_mean,             aplot$data[[2]])
  df_pred$bap_hwd_sc   = sc(df_pred$bap_hwd_mean,          aplot$data[[3]])
  df_pred$standinit_sc = sc(df_pred$pcnt_standinit, ahome$data[[1]][[sp_name]])
  df_pred$thin_sc      = sc(df_pred$pcnt_thin,       ahome$data[[2]][[sp_name]])
  df_pred$mature_sc    = sc(df_pred$pcnt_mature,     ahome$data[[3]][[sp_name]])
  
  # --- Extract posterior parameters --------------------------------------------------
  sp_idx = which(species == sp_name)
  
  u_i  = sl$u[, sp_idx]
  ap1  = sl$alpha_plot1[, , sp_idx]
  ap2  = sl$alpha_plot2[, , sp_idx]
  ap3  = sl$alpha_plot3[, , sp_idx]
  apt1 = sl$alpha_point1[, sp_idx]
  apt2 = sl$alpha_point2[, sp_idx]
  apt3 = sl$alpha_point3[, sp_idx]
  ahr1 = sl$alpha_homerange1[, sp_idx]
  ahr2 = sl$alpha_homerange2[, sp_idx]
  ahr3 = sl$alpha_homerange3[, sp_idx]
  aseas = sl$alpha_season[, sp_idx]
  
  # --- Chunk loop --------------------------------------------------------------------
  n_cells      = nrow(df_pred)
  chunk_size   = 1000
  n_chunks     = ceiling(n_cells / chunk_size)
  cell_area_ha = prod(res(rast_stack)) / 10000
  sv           = df_pred$stage_idx
  
  psi_mean = numeric(n_cells)
  psi_sd   = numeric(n_cells)
  psi_q025 = numeric(n_cells)
  psi_q975 = numeric(n_cells)
  
  expected_ha_samples  = numeric(n_iter)
  expected_ha_by_stage = matrix(0, nrow = n_iter, ncol = length(stage_levels),
                                dimnames = list(NULL, stage_levels))
  
  # Initialize richness accumulators after first species (n_cells now known)
  if (is.null(richness_sum)) {
    message("  Initializing richness accumulators: ",
            n_cells, " cells × ", n_thin, " thinned draws (",
            round(n_cells * n_thin * 8 / 1e9, 2), " GB)")
    richness_sum      = numeric(n_cells)
    richness_thin_mat = matrix(0, nrow = n_thin, ncol = n_cells)
    richness_xy       = df_pred[, c("x", "y")]
  }
  
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
    
    psi_ch = plogis(logit_psi_ch)  # [n_iter × chunk_size]
    
    psi_mean[idx] = colMeans(psi_ch)
    psi_sd[idx]   = apply(psi_ch, 2, sd)
    psi_q025[idx] = apply(psi_ch, 2, quantile, 0.025)
    psi_q975[idx] = apply(psi_ch, 2, quantile, 0.975)
    
    expected_ha_samples = expected_ha_samples + rowSums(psi_ch) * cell_area_ha
    
    for (s in seq_along(stage_levels)) {
      mask_s = sv_ch == s
      if (any(mask_s)) {
        expected_ha_by_stage[, s] = expected_ha_by_stage[, s] +
          rowSums(psi_ch[, mask_s, drop = FALSE]) * cell_area_ha
      }
    }
    
    # Richness: exact sum over all iterations (for posterior mean)
    richness_sum[idx] = richness_sum[idx] + colSums(psi_ch)
    
    # Richness: thinned draws accumulated for SD / quantiles
    richness_thin_mat[, idx] = richness_thin_mat[, idx] + psi_ch[thin_idx, , drop = FALSE]
    
    if (ch %% 10 == 0) message("  Chunk ", ch, " / ", n_chunks)
  }
  
  df_pred$psi_mean = psi_mean
  df_pred$psi_sd   = psi_sd
  df_pred$psi_q025 = psi_q025
  df_pred$psi_q975 = psi_q975
  
  rast_psi = rast(
    df_pred[, c("x", "y", "psi_mean", "psi_sd", "psi_q025", "psi_q975")],
    type = "xyz", crs = crs(rast_stack)
  )
  
  # --- Summaries ---------------------------------------------------------------------
  total_landscape_ha = n_cells * cell_area_ha
  
  area_summary = tibble(
    species            = sp_name,
    mean_ha            = mean(expected_ha_samples),
    sd_ha              = sd(expected_ha_samples),
    q025_ha            = quantile(expected_ha_samples, 0.025),
    median_ha          = quantile(expected_ha_samples, 0.500),
    q975_ha            = quantile(expected_ha_samples, 0.975),
    total_landscape_ha = total_landscape_ha,
    pct_landscape      = mean(expected_ha_samples) / total_landscape_ha * 100
  )
  
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
  
  # --- Figures -----------------------------------------------------------------------
  fig_map = ggplot() +
    geom_spatraster(data = rast_psi[["psi_mean"]]) +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 1), na.value = NA) +
    labs(title = sp_name, fill = "Ψ") +
    theme(legend.position = "right")
  
  fig_map_sd = ggplot() +
    geom_spatraster(data = rast_psi[["psi_sd"]]) +
    scale_fill_viridis_c(option = "magma", na.value = NA) +
    labs(title = "", fill = "SD") +
    theme(legend.position = "right",
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
  
  fig_hist = ggplot(df_pred, aes(x = psi_mean, fill = after_stat(x))) +
    geom_histogram(bins = 100, color = NA) +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Posterior mean use probability (Ψ)", y = "Number of cells", title = sp_name) +
    theme(legend.position = "none")
  
  fig_stage = ggplot(stage_summary, aes(x = stratum)) +
    geom_col(aes(y = total_ha), fill = "gray80") +
    geom_col(aes(y = mean_ha, fill = stratum), alpha = 0.9) +
    geom_errorbar(aes(ymin = q025_ha, ymax = q975_ha), width = 0.2, color = "gray10") +
    scale_fill_manual(values = unname(stage_colors)) +
    labs(x = NULL, y = "Area (ha)", title = sp_name) +
    theme(legend.position = "none")
  
  species_results[[sp_name]] = list(
    area_summary  = area_summary,
    stage_summary = stage_summary,
    rast_psi      = rast_psi,
    df_pred       = df_pred[, c("x", "y", "psi_mean", "psi_sd", "psi_q025", "psi_q975")],
    figures       = list(
      map       = fig_map,
      map_sd    = fig_map_sd,
      histogram = fig_hist,
      by_stage  = fig_stage
    )
  )
  
  message(sp_name, " — expected area used:")
  print(area_summary)
  message(sp_name, " — expected area used by stage:")
  print(stage_summary)
  i = i + 1
}  # end species loop

# =====================================================================================
# SPECIES RICHNESS
# =====================================================================================

# Posterior mean: exact (sum of all iterations / n_iter, accumulated across species)
# SD and quantiles: from thinned draws accumulated across species
# richness_thin_mat[i, c] = ∑_sp psi_sp(thin_iter_i, cell_c) — i.e., one richness draw per thinned iter

richness_df = richness_xy %>%
  mutate(
    richness_mean = richness_sum / n_iter,                              # exact
    richness_sd   = apply(richness_thin_mat, 2, sd),                   # from thinned draws
    richness_q025 = apply(richness_thin_mat, 2, quantile, 0.025),      # from thinned draws
    richness_q975 = apply(richness_thin_mat, 2, quantile, 0.975)       # from thinned draws
  )

rast_richness = rast(
  richness_df[, c("x", "y", "richness_mean", "richness_sd", "richness_q025", "richness_q975")],
  type = "xyz", crs = crs(rast_stack)
)

# Landscape-level posterior from thinned draws: rowSums gives total species×ha per iteration
richness_landscape_iter    = rowSums(richness_thin_mat) * cell_area_ha
richness_landscape_summary = tibble(
  mean_species_x_ha   = sum(richness_sum) / n_iter * cell_area_ha,   # exact mean
  sd_species_x_ha     = sd(richness_landscape_iter),                  # from thinned
  q025_species_x_ha   = quantile(richness_landscape_iter, 0.025),
  median_species_x_ha = quantile(richness_landscape_iter, 0.500),
  q975_species_x_ha   = quantile(richness_landscape_iter, 0.975)
)

message("Landscape-level expected total species × ha:")
print(richness_landscape_summary)

# --- Richness figures ----------------------------------------------------------------
fig_richness_map = ggplot() +
  geom_spatraster(data = rast_richness[["richness_mean"]]) +
  scale_fill_viridis_c(option = "turbo", na.value = NA) +
  labs(title = "Expected species richness E[S]", fill = "E[S]") +
  theme(legend.position = "right")

fig_richness_sd = ggplot() +
  geom_spatraster(data = rast_richness[["richness_sd"]]) +
  scale_fill_viridis_c(option = "magma", na.value = NA) +
  labs(title = "Richness uncertainty (posterior SD)", fill = "SD") +
  theme(legend.position = "right",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

fig_richness_ci_width = ggplot() +
  geom_spatraster(data = rast_richness[["richness_q975"]] - rast_richness[["richness_q025"]]) +
  scale_fill_viridis_c(option = "magma", na.value = NA) +
  labs(title = "Richness 95% BCI width", fill = "Width") +
  theme(legend.position = "right",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

fig_richness_hist = ggplot(richness_df, aes(x = richness_mean, fill = after_stat(x))) +
  geom_histogram(bins = 60, color = NA) +
  scale_fill_viridis_c(option = "turbo") +
  labs(x = "Expected species richness E[S]", y = "Number of cells",
       title = "Distribution of expected richness across cells") +
  theme(legend.position = "none")

# --- Bundle all outputs --------------------------------------------------------------
all_area_summaries = bind_rows(lapply(species_results, `[[`, "area_summary"))

richness_results = list(
  richness_df                = richness_df,
  rast_richness              = rast_richness,
  richness_landscape_summary = richness_landscape_summary,
  all_area_summaries         = all_area_summaries,
  figures = list(
    richness_map      = fig_richness_map,
    richness_sd       = fig_richness_sd,
    richness_ci_width = fig_richness_ci_width,
    richness_hist     = fig_richness_hist
  )
)

dir.create(dirname(path_out_species_results), recursive = TRUE, showWarnings = FALSE)
# saveRDS(species_results, path_out_species_results)
# 
dir.create(dirname(path_out_richness_results), recursive = TRUE, showWarnings = FALSE)
# saveRDS(richness_results, path_out_richness_results)

# ---- saving (end of prediction script) ---------------------------------------------

# Wrap all SpatRasters before caching — converts external pointer to serializable object
richness_results_cache = richness_results
richness_results_cache$rast_richness = wrap(richness_results$rast_richness)
richness_results_cache$figures       = NULL   # never cache ggplot objects with rasters

saveRDS(richness_results_cache, path_out_richness_results)

# Same for species_results
species_results_cache = lapply(species_results, function(x) {
  x$rast_psi = wrap(x$rast_psi)
  x$figures  = NULL
  x
})
saveRDS(species_results_cache, path_out_species_results)

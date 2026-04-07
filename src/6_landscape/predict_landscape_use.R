#####################################################################################
# Predict habitat use across the landscape
#
# INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"

source("src/global.R")

# Load MSOM ---------------------------------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
stages = model_data$stages
strata = factor(stages$stratum_4, levels = c("standinit", "compex", "thin", "mature"))

z = msom$sims.list$z
n_iter    = dim(z)[1]
n_sites   = dim(z)[2]
n_seasons = dim(z)[3]
n_species = dim(z)[4]

# Load plot and homerange scale data for all points -----------------------------------
path_rast_cover = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_2020_clean_strata_4.tif")

message('Loading raster cover data from cache ', path_rast_cover)
rast_cover = rast(path_rast_cover)

path_homerange = "data/cache/7_landscape/OPT_landscape_data_homerange_scale_2020_clean_strata_4_landscape.rds"
path_point      = "data/cache/7_landscape/OPT_landscape_data_plot_scale_2020_clean_strata_4_landscape.rds"

data_homerange = read_rds(path_homerange)
names(data_homerange)

points_point = read_rds(path_point)
str(points_point)

points_plot = st_as_sf(data_homerange[["plot"]] %>% select(
  "x", "y", "bap_hwd_mean", "qmd_6_mean", "tree_acre_6_mean"
), coords = c("x", "y"), crs = crs_m)

sp_name = "brown creeper"
str(data_homerange[[sp_name]])

points_homerange = st_as_sf(data_homerange[[sp_name]] %>% select(
  "x", "y", "pcnt_standinit", "pcnt_compex", "pcnt_thin", "pcnt_mature"
), coords = c("x", "y"), crs = crs_m)

points = points_homerange %>% st_join(points_point, join = st_equals)
points = points %>% st_join(points_plot, join = st_equals)

# TODO: Impute NA layers smarter
if (anyNA(points)) {
  message(crayon::yellow("WARNING: Imputing NA values from layer means for now"))
  colSums(is.na(st_drop_geometry(points)))
  
  points = points %>%
    mutate(across(
      where(~ is.numeric(.) && anyNA(.)),  # only numeric columns with NAs
      ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
    ))
  colSums(is.na(st_drop_geometry(points)))
}

# Get stage from rast cover
stage_vals = extract(rast_cover, vect(points))
points$stage = stage_vals$stage
points = points %>% left_join(stages %>% distinct() %>% select(stratum_4, stage_idx), by = c("stage" = "stratum_4"))

# mapview(points)

point_data = cbind(st_coordinates(points), st_drop_geometry(points))
str(point_data)

vars = c(
  # param_alpha_homerange_data
  "pcnt_standinit", "pcnt_compex", "pcnt_mature", "pcnt_thin",
  # param_alpha_point_data
   "elevation", "dist_road_paved", "dist_watercourse_major",
  # param_alpha_plot_data
  "tree_acre_6_mean", "qmd_6_mean", "bap_hwd_mean",
  # stage
  "stage_idx"
)

rasters = lapply(vars, function(v) {
  message(v)
  rast(point_data[, c("X", "Y", v)], type = "xyz", crs = crs_m_rast)
})
names(rasters) = vars

mapview(rasters[["pcnt_mature"]])
plot(rasters[["stage_idx"]])

rast_stack = rast(rasters)
plot(rast_stack)

# Predict use probability -----------------------------------

# ── Helper ───────────────────────────────────────────────────────────────────
sc <- function(x, scaled_vec) {
  (x - attr(scaled_vec, "scaled:center")) / attr(scaled_vec, "scaled:scale")
}

# ── CONFIG ────────────────────────────────────────────────────────────────────
season_t <- "2020"

sp_idx <- which(species == sp_name)
t_idx  <- which(seasons == season_t)

# Unpack param data for readability
pad <- model_data$param_alpha_data
apoint <- pad$param_alpha_point_data    # alpha_point1/2/3 → elevation/road/water
aplot  <- pad$param_alpha_plot_data     # alpha_plot1/2/3  → tree_acre/qmd/bap_hwd
ahome  <- pad$param_alpha_homerange_data # alpha_homerange1/2/3 → standinit/thin/mature

# ── Scale predictor raster cells using training attributes ────────────────────
df_pred <- as.data.frame(rast_stack, xy = TRUE) |> na.omit()

# Point-scale (shared scaling across all species)
df_pred$elev_sc  <- sc(df_pred$elevation,             apoint$data[[1]])
df_pred$road_sc  <- sc(df_pred$dist_road_paved,       apoint$data[[2]])
df_pred$water_sc <- sc(df_pred$dist_watercourse_major, apoint$data[[3]])

# Plot-scale (shared scaling)
df_pred$tree_acre_sc <- sc(df_pred$tree_acre_6_mean, aplot$data[[1]])
df_pred$qmd_sc       <- sc(df_pred$qmd_6_mean,       aplot$data[[2]])
df_pred$bap_hwd_sc   <- sc(df_pred$bap_hwd_mean,     aplot$data[[3]])

# Homerange-scale (species-specific scaling — index into the species sublist)
df_pred$standinit_sc <- sc(df_pred$pcnt_standinit, ahome$data[[1]][[sp_name]])
df_pred$thin_sc      <- sc(df_pred$pcnt_thin,       ahome$data[[2]][[sp_name]])
df_pred$mature_sc    <- sc(df_pred$pcnt_mature,      ahome$data[[3]][[sp_name]])

# ── Season scalar ─────────────────────────────────────────────────────────────
x_season_t <- as.vector(pad$param_alpha_season |> 
                          (\(.) scale(1:length(seasons)))())[t_idx]
# Or more simply, since you have it already:
x_season_t <- as.vector(scale(1:length(seasons)))[t_idx]

# ── Posterior samples for this species ───────────────────────────────────────
sl  <- msom$sims.list
n_iter  <- dim(sl$u)[1]
n_cells <- nrow(df_pred)
sv      <- df_pred$stage_idx  # stage vector [n_cells]

# Pre-extract all iterations for this species [n_iter x S] or [n_iter]
u_i   <- sl$u[, sp_idx]
ap1   <- sl$alpha_plot1[, , sp_idx]   # [n_iter, S]
ap2   <- sl$alpha_plot2[, , sp_idx]
ap3   <- sl$alpha_plot3[, , sp_idx]
apt1  <- sl$alpha_point1[, sp_idx]
apt2  <- sl$alpha_point2[, sp_idx]
apt3  <- sl$alpha_point3[, sp_idx]
ahr1  <- sl$alpha_homerange1[, sp_idx]
ahr2  <- sl$alpha_homerange2[, sp_idx]
ahr3  <- sl$alpha_homerange3[, sp_idx]
aseas <- sl$alpha_season[, sp_idx]

# ── Vectorised posterior summary (avoids R loop over iterations) ──────────────
# Build design matrices: [n_iter x n_cells], broadcasting stage-indexed coefficients
# ap1[, sv] selects the correct stage column for each cell → [n_iter x n_cells]

# logit_psi_mat <-
#   outer(u_i,   rep(1, n_cells)) +
#   ap1[, sv]  * outer(rep(1, n_iter), df_pred$tree_acre_sc) +
#   ap2[, sv]  * outer(rep(1, n_iter), df_pred$qmd_sc) +
#   ap3[, sv]  * outer(rep(1, n_iter), df_pred$bap_hwd_sc) +
#   outer(apt1, df_pred$elev_sc) +
#   outer(apt2, df_pred$road_sc) +
#   outer(apt3, df_pred$water_sc) +
#   outer(ahr1, df_pred$standinit_sc) +
#   outer(ahr2, df_pred$thin_sc) +
#   outer(ahr3, df_pred$mature_sc) +
#   outer(aseas, rep(x_season_t, n_cells))
# 
# psi_mat <- plogis(logit_psi_mat)  # [n_iter x n_cells]
# 
# df_pred$psi_mean <- colMeans(psi_mat)
# df_pred$psi_sd   <- apply(psi_mat, 2, sd)
# df_pred$psi_q025 <- apply(psi_mat, 2, quantile, 0.025)
# df_pred$psi_q975 <- apply(psi_mat, 2, quantile, 0.975)

chunk_size <- 1000  # tune this down if still hitting memory
n_chunks   <- ceiling(n_cells / chunk_size)

psi_mean <- numeric(n_cells)
psi_sd   <- numeric(n_cells)
psi_q025 <- numeric(n_cells)
psi_q975 <- numeric(n_cells)

for (ch in seq_len(n_chunks)) {
  idx <- ((ch - 1) * chunk_size + 1):min(ch * chunk_size, n_cells)
  sv_ch <- sv[idx]
  
  logit_psi_ch <-
    outer(u_i, rep(1, length(idx))) +
    ap1[, sv_ch] * outer(rep(1, n_iter), df_pred$tree_acre_sc[idx]) +
    ap2[, sv_ch] * outer(rep(1, n_iter), df_pred$qmd_sc[idx]) +
    ap3[, sv_ch] * outer(rep(1, n_iter), df_pred$bap_hwd_sc[idx]) +
    outer(apt1, df_pred$elev_sc[idx]) +
    outer(apt2, df_pred$road_sc[idx]) +
    outer(apt3, df_pred$water_sc[idx]) +
    outer(ahr1, df_pred$standinit_sc[idx]) +
    outer(ahr2, df_pred$thin_sc[idx]) +
    outer(ahr3, df_pred$mature_sc[idx]) +
    outer(aseas, rep(x_season_t, length(idx)))
  
  psi_ch <- plogis(logit_psi_ch)  # [n_iter x chunk_size]
  
  psi_mean[idx] <- colMeans(psi_ch)
  psi_sd[idx]   <- apply(psi_ch, 2, sd)
  psi_q025[idx] <- apply(psi_ch, 2, quantile, 0.025)
  psi_q975[idx] <- apply(psi_ch, 2, quantile, 0.975)
  
  if (ch %% 10 == 0) message("Chunk ", ch, " / ", n_chunks)
}

df_pred$psi_mean <- psi_mean
df_pred$psi_sd   <- psi_sd
df_pred$psi_q025 <- psi_q025
df_pred$psi_q975 <- psi_q975

# ── Back to raster ────────────────────────────────────────────────────────────
rast_psi <- rast(
  df_pred[, c("x", "y", "psi_mean", "psi_sd", "psi_q025", "psi_q975")],
  type = "xyz", crs = crs(rast_stack)
)
mapview(rast_psi[["psi_mean"]]) + mapview(rast_psi[["psi_sd"]])

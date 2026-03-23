####################################################################################
# Marginal plots of varying homerange stage composition
#
# CONFIG:
path_msom = "data/cache/models/V3_msom_pcnt_fp_fp_all.rds"
#
####################################################################################

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

(msom_summary = model_data$msom_summary)
(msom = model_data$msom)
(groups = model_data$groups %>% arrange(common_name))
(sites = model_data$sites)
(species = model_data$species)
(stages = model_data$stages)

sp_name <- "hairy woodpecker"

# ── 1. SPECIES SETUP ──────────────────────────────────────────────────────────
sp_idx  <- which(model_data$species == sp_name)  # 39

sims    <- model_data$msom$sims.list
n_draws <- dim(sims$u)[1]

# ── 2. PREDICTION GRID ────────────────────────────────────────────────────────
# Raw proportions: only compex <-> mature gradient, SI and thin fixed at 0
p_mature_raw <- seq(0, 1, by = 0.01)   # 101 grid points
p_compex_raw <- 1 - p_mature_raw        # implicit, never enters model
p_si_raw     <- 0
p_thin_raw   <- 0
n_grid       <- length(p_mature_raw)

# Helper to extract species-specific scaling attributes
get_hr_scaling <- function(hr_param_name, species_name) {
  x <- model_data$param_alpha_data$param_alpha_homerange_data |>
    filter(name == hr_param_name) |>
    pull(data) |> _[[1]]   # named list of 67
  list(
    center = attr(x[[species_name]], "scaled:center"),
    scale  = attr(x[[species_name]], "scaled:scale")
  )
}

# Scale raw proportions using species-specific attributes
# For SI and thin: raw value is 0 throughout, but still needs scaling
hr3_sc <- get_hr_scaling("pcnt_mature",    sp_name)
hr1_sc <- get_hr_scaling("pcnt_standinit", sp_name)
hr2_sc <- get_hr_scaling("pcnt_thin",      sp_name)

p_mature_scaled <- (p_mature_raw - hr3_sc$center) / hr3_sc$scale
p_si_scaled     <- (p_si_raw    - hr1_sc$center)  / hr1_sc$scale   # scalar
p_thin_scaled   <- (p_thin_raw  - hr2_sc$center)  / hr2_sc$scale   # scalar

# ── 3. STAGE-CONDITIONAL PLOT COVARIATE MEANS ─────────────────────────────────
# Stage coding: 1=compex, 2=standinit, 3=mature, 4=thin
# Stage vector is site-level (224 sites); expand to match 896-row covariate arrays
stage_vec      <- model_data$stages$stage_idx
n_sites        <- length(stage_vec)
n_rows         <- nrow(model_data$param_alpha_data$param_alpha_plot_data$data[[1]])
n_reps         <- n_rows / n_sites   # 4 seasons
stage_expanded <- rep(stage_vec, times = n_reps)  # length 896

# Extract globally-scaled plot covariate vectors
plot_data <- model_data$param_alpha_data$param_alpha_plot_data
x_plot1   <- as.vector(plot_data$data[[1]])  # tree_acre_6_mean
x_plot2   <- as.vector(plot_data$data[[2]])  # qmd_6_mean
x_plot3   <- as.vector(plot_data$data[[3]])  # bap_hwd_mean

# Stage-conditional means on globally-scaled values
stage_plot_means <- function(stage_code) {
  idx <- stage_expanded == stage_code
  c(
    mean(x_plot1[idx], na.rm = TRUE),
    mean(x_plot2[idx], na.rm = TRUE),
    mean(x_plot3[idx], na.rm = TRUE)
  )
}

mean_plot_compex <- stage_plot_means(1)
mean_plot_mature <- stage_plot_means(3)

# Sanity check — mature should be above-zero on QMD relative to global mean
data.frame(
  covariate   = c("tree_acre_6_mean", "qmd_6_mean", "bap_hwd_mean"),
  mean_compex = round(mean_plot_compex, 3),
  mean_mature = round(mean_plot_mature, 3)
)

# ── 4. EXTRACT POSTERIOR SLICES FOR SPECIES ───────────────────────
# alpha_plot: [n_draws, n_stages, n_species]  — verify with dim(sims$alpha_plot1)
# alpha_season: [n_draws, n_species]          — verify with dim(sims$alpha_season)
u_sp    <- sims$u[, sp_idx]

# Plot coefficients for compex (stage 1) and mature (stage 3)
ap1_ce  <- sims$alpha_plot1[, 1, sp_idx]
ap2_ce  <- sims$alpha_plot2[, 1, sp_idx]
ap3_ce  <- sims$alpha_plot3[, 1, sp_idx]

ap1_mat <- sims$alpha_plot1[, 3, sp_idx]
ap2_mat <- sims$alpha_plot2[, 3, sp_idx]
ap3_mat <- sims$alpha_plot3[, 3, sp_idx]

# Homerange coefficients
ahr1    <- sims$alpha_homerange1[, sp_idx]   # pcnt_standinit (fixed at 0 raw)
ahr2    <- sims$alpha_homerange2[, sp_idx]   # pcnt_thin      (fixed at 0 raw)
ahr3    <- sims$alpha_homerange3[, sp_idx]   # pcnt_mature    (varies)

# ── 5. BUILD LINEAR PREDICTORS ────────────────────────────────────────────────
# Point covariates and season: global scaling → grand mean = 0 → drop out
# Homerange: SI and thin are fixed scalars; mature varies across grid

# Fixed homerange contribution from SI and thin [n_draws scalar]
hr_fixed <- ahr1 * p_si_scaled + ahr2 * p_thin_scaled

# Varying homerange contribution from mature [n_draws x n_grid]
hr_mature <- outer(ahr3, p_mature_scaled)   # [n_draws x n_grid]

# Full homerange contribution [n_draws x n_grid]
hr_contrib <- hr_mature + hr_fixed   # hr_fixed recycles across columns

# Plot contributions [n_draws scalar each]
plot_contrib_compex <- ap1_ce  * mean_plot_compex[1] +
  ap2_ce  * mean_plot_compex[2] +
  ap3_ce  * mean_plot_compex[3]

plot_contrib_mature <- ap1_mat * mean_plot_mature[1] +
  ap2_mat * mean_plot_mature[2] +
  ap3_mat * mean_plot_mature[3]

# Full linear predictors [n_draws x n_grid]
# u_sp, plot_contrib_* are [n_draws] vectors — recycle across columns (column-major)
eta_compex <- hr_contrib + (u_sp + plot_contrib_compex)
eta_mature <- hr_contrib + (u_sp + plot_contrib_mature)

# ── 6. OPTION B: WEIGHTED AVERAGE ACROSS STAGE PROBABILITIES ─────────────────
# psi[d,g] = p_compex[g] * plogis(eta_compex[d,g])
#           + p_mature[g]  * plogis(eta_mature[d,g])
# SI and thin weights are 0 throughout so drop out
p_compex_mat <- matrix(p_compex_raw, nrow = n_draws, ncol = n_grid, byrow = TRUE)
p_mature_mat <- matrix(p_mature_raw, nrow = n_draws, ncol = n_grid, byrow = TRUE)

psi_draws <- p_compex_mat * plogis(eta_compex) +
  p_mature_mat  * plogis(eta_mature)   # [n_draws x n_grid]

# ── 7. POSTERIOR SUMMARY ──────────────────────────────────────────────────────
pred_df <- tibble(
  p_mature  = p_mature_raw,
  p_compex  = p_compex_raw,
  psi_mean  = apply(psi_draws, 2, mean),
  psi_lower = apply(psi_draws, 2, quantile, 0.025),
  psi_upper = apply(psi_draws, 2, quantile, 0.975)
)

# ── 8. PLOT ───────────────────────────────────────────────────────────────────
ggplot(pred_df, aes(x = p_mature)) +
  geom_ribbon(aes(ymin = psi_lower, ymax = psi_upper),
              alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = psi_mean), colour = "steelblue", linewidth = 1) +
  scale_x_continuous(
    labels   = scales::percent,
    sec.axis = sec_axis(
      ~ 1 - .,
      name   = "% Competitive exclusion",
      labels = scales::percent
    )
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x        = "% Mature",
    y        = "Predicted occupancy probability",
    title    = sp_name,
    subtitle = paste(
      "Marginal occupancy across competitive exclusion \u2013 mature gradient",
      "p(SI) = 0, p(thin) = 0 | point covariates at grand mean | plot covariates at stage-conditional means",
      sep = "\n"
    )
  )

### find the range of covariate values that were actually observed in the data -------------------

# Extract the raw observed values for species's pcnt_mature
hr3_data <- model_data$param_alpha_data$param_alpha_homerange_data |>
  filter(name == "pcnt_mature") |>
  pull(data) |> _[[1]]   # named list of 67 species

pw_vals_scaled <- as.vector(hr3_data[[sp_name]])

# Back-transform to raw proportions
pw_vals_raw <- pw_vals_scaled * hr3_sc$scale + hr3_sc$center

summary(pw_vals_raw)
quantile(pw_vals_raw, c(0.05, 0.25, 0.5, 0.75, 0.95))

obs_range <- quantile(pw_vals_raw, c(0.05, 0.95))

ggplot(pred_df, aes(x = p_mature)) +
  # Extrapolation zones
  annotate("rect",
           xmin = 0, xmax = obs_range[1],
           ymin = 0, ymax = 1,
           fill = "grey80", alpha = 0.4) +
  annotate("rect",
           xmin = obs_range[2], xmax = 1,
           ymin = 0, ymax = 1,
           fill = "grey80", alpha = 0.4) +
  # Credible interval and mean
  geom_ribbon(aes(ymin = psi_lower, ymax = psi_upper),
              alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = psi_mean), colour = "steelblue", linewidth = 1) +
  scale_x_continuous(
    labels   = scales::percent,
    sec.axis = sec_axis(~ 1 - ., name = "% Competitive exclusion",
                        labels = scales::percent)
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x        = "% Mature",
    y        = "Predicted occupancy probability",
    title    = sp_name,
    subtitle = paste(
      "Marginal occupancy across competitive exclusion \u2013 mature gradient",
      "Shaded regions indicate extrapolation beyond 5th\u201395th percentile of observed data",
      sep = "\n"
    )
  )

# check the observed ranges for all three homerange predictors at once to get a full picture of where the data actually are
hr_names <- c("pcnt_standinit", "pcnt_thin", "pcnt_mature")
hr_attrs <- list(hr1_sc, hr2_sc, hr3_sc)

map2_dfr(hr_names, hr_attrs, function(nm, sc) {
  vals_scaled <- model_data$param_alpha_data$param_alpha_homerange_data |>
    filter(name == nm) |>
    pull(data) |> _[[1]] |> _[[sp_name]] |> as.vector()
  vals_raw <- vals_scaled * sc$scale + sc$center
  tibble(
    covariate = nm,
    min       = min(vals_raw),
    p05       = quantile(vals_raw, 0.05),
    p25       = quantile(vals_raw, 0.25),
    median    = median(vals_raw),
    p75       = quantile(vals_raw, 0.75),
    p95       = quantile(vals_raw, 0.95),
    max       = max(vals_raw)
  )
})

## Clip bounds to the observed max, and add a rug for observations

ggplot(pred_df, aes(x = p_mature)) +
  annotate("rect",
           xmin = obs_range[2], xmax = 1,
           ymin = 0, ymax = 1,
           fill = "grey80", alpha = 0.4) +
  geom_ribbon(aes(ymin = psi_lower, ymax = psi_upper),
              alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = psi_mean), colour = "steelblue", linewidth = 1) +
  geom_rug(
    data = tibble(x = pw_vals_raw),
    aes(x = x), inherit.aes = FALSE,
    alpha = 0.3, length = unit(0.03, "npc")
  ) +
  scale_x_continuous(
    limits   = c(0, max(pw_vals_raw)),   # clip to observed max
    labels   = scales::percent,
    sec.axis = sec_axis(~ 1 - ., name = "% Competitive exclusion",
                        labels = scales::percent)
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x        = "% Mature",
    y        = "Predicted occupancy probability",
    title    = sp_name,
    subtitle = paste(
      "Marginal occupancy across competitive exclusion \u2013 mature gradient",
      "Shaded region: extrapolation beyond 95th percentile | Rug: observed homerange values",
      sep = "\n"
    )
  )

# DEBUGGING

# Baseline occupancy
msom_summary %>% filter(param == paste0("u[", sp_idx, "]"))

# 1. How far does 100% mature extrapolate in standardised units for the species?
hr3_sc_bc <- get_hr_scaling("pcnt_mature", sp_name)
(scaled_100pct <- (1 - hr3_sc_bc$center) / hr3_sc_bc$scale)

# 2. What is the species' homerange3 coefficient?
msom_summary |> filter(param == "alpha_homerange3[11]")

# 3. What is the resulting linear predictor contribution at 100% mature?
mean(sims$alpha_homerange3[, sp_idx]) * scaled_100pct

# 4. What was actually observed?
hr3_data <- model_data$param_alpha_data$param_alpha_homerange_data |>
  filter(name == "pcnt_mature") |>
  pull(data) |> _[[1]]

bc_vals_raw <- as.vector(hr3_data[[sp_name]]) * hr3_sc_bc$scale + hr3_sc_bc$center
quantile(bc_vals_raw, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

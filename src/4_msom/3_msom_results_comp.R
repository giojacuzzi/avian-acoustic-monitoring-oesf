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

library(tidyverse)
library(patchwork)

library(tidyverse)
library(patchwork)

library(tidyverse)

# ── SETUP (run once) ──────────────────────────────────────────────────────────
sims    <- model_data$msom$sims.list
n_draws <- dim(sims$u)[1]

# Stage name -> index mapping (1=compex, 2=standinit, 3=mature, 4=thin)
stage_idx_map <- c(compex = 1, standinit = 2, mature = 3, thin = 4)

# Stage-conditional plot covariate means (globally scaled, computed once)
stage_vec      <- model_data$stages$stage_idx
n_rows         <- nrow(model_data$param_alpha_data$param_alpha_plot_data$data[[1]])
stage_expanded <- rep(stage_vec, times = n_rows / length(stage_vec))

plot_data <- model_data$param_alpha_data$param_alpha_plot_data
x_plot    <- lapply(plot_data$data, as.vector)  # list of 3 globally-scaled vectors

stage_plot_means <- function(stage_code) {
  idx <- stage_expanded == stage_code
  sapply(x_plot, \(x) mean(x[idx], na.rm = TRUE))
}

mean_plot <- lapply(stage_idx_map, stage_plot_means)
names(mean_plot) <- names(stage_idx_map)

# Helpers
get_hr_scaling <- function(hr_param_name, species_name) {
  x <- model_data$param_alpha_data$param_alpha_homerange_data |>
    filter(name == hr_param_name) |>
    pull(data) |> _[[1]]
  list(center = attr(x[[species_name]], "scaled:center"),
       scale  = attr(x[[species_name]], "scaled:scale"))
}

get_hr_raw <- function(hr_param_name, species_name, sc) {
  model_data$param_alpha_data$param_alpha_homerange_data |>
    filter(name == hr_param_name) |>
    pull(data) |> _[[1]] |> _[[species_name]] |> as.vector() |>
    (\(x) x * sc$scale + sc$center)()
}

scale_proportion <- function(p, sc) (p - sc$center) / sc$scale

# ── PREDICTION FUNCTION ───────────────────────────────────────────────────────
predict_hr_gradient <- function(species_name,
                                from = "compex",
                                to   = "mature") {
  
  stopifnot(from %in% names(stage_idx_map), to %in% names(stage_idx_map))
  stopifnot(from != to)
  
  sp_idx     <- which(model_data$species == species_name)
  n_grid     <- 101
  p_to_raw   <- seq(0, 1, length.out = n_grid)  # "to" stage goes 0 -> 1
  p_from_raw <- 1 - p_to_raw                     # "from" stage goes 1 -> 0
  
  # Homerange scaling for all three non-compex stages
  hr_sc <- list(
    standinit = get_hr_scaling("pcnt_standinit", species_name),
    thin      = get_hr_scaling("pcnt_thin",      species_name),
    mature    = get_hr_scaling("pcnt_mature",     species_name)
  )
  
  # Posterior slices
  u_sp   <- sims$u[, sp_idx]
  get_ap <- function(stage, k) sims[[paste0("alpha_plot", k)]][, stage_idx_map[stage], sp_idx]
  ap_from <- lapply(1:3, \(k) get_ap(from, k))
  ap_to   <- lapply(1:3, \(k) get_ap(to,   k))
  
  ahr <- list(
    standinit = sims$alpha_homerange1[, sp_idx],
    thin      = sims$alpha_homerange2[, sp_idx],
    mature    = sims$alpha_homerange3[, sp_idx]
  )
  
  # ── Homerange contributions ──────────────────────────────────────────────────
  # For each non-compex stage:
  #   "to" or "from" stage: p varies across grid -> outer() -> [n_draws x n_grid]
  #   fixed stage (held at 0%): p is scalar -> [n_draws] vector, recycles across columns
  hr_contrib <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (s in c("standinit", "thin", "mature")) {
    p_s_scaled <- if (s == to) {
      scale_proportion(p_to_raw,   hr_sc[[s]])   # vector, length n_grid
    } else if (s == from) {
      scale_proportion(p_from_raw, hr_sc[[s]])   # vector, length n_grid
    } else {
      scale_proportion(0, hr_sc[[s]])             # scalar
    }
    
    contribution <- if (length(p_s_scaled) > 1) {
      outer(ahr[[s]], p_s_scaled)   # [n_draws x n_grid]
    } else {
      ahr[[s]] * p_s_scaled         # [n_draws] vector, recycles across columns
    }
    
    hr_contrib <- hr_contrib + contribution
  }
  
  # ── Plot contributions [n_draws vector] ──────────────────────────────────────
  plot_contrib <- function(stage, ap_list) {
    ap_list[[1]] * mean_plot[[stage]][1] +
      ap_list[[2]] * mean_plot[[stage]][2] +
      ap_list[[3]] * mean_plot[[stage]][3]
  }
  pc_from <- plot_contrib(from, ap_from)
  pc_to   <- plot_contrib(to,   ap_to)
  
  # ── Full linear predictors [n_draws x n_grid] ────────────────────────────────
  # u_sp and pc_* are [n_draws] vectors — recycle across columns (column-major)
  eta_from <- hr_contrib + (u_sp + pc_from)
  eta_to   <- hr_contrib + (u_sp + pc_to)
  
  # ── Option B: weighted average of stage probabilities ────────────────────────
  p_from_mat <- matrix(p_from_raw, nrow = n_draws, ncol = n_grid, byrow = TRUE)
  p_to_mat   <- matrix(p_to_raw,   nrow = n_draws, ncol = n_grid, byrow = TRUE)
  
  psi_draws <- p_from_mat * plogis(eta_from) +
    p_to_mat   * plogis(eta_to)
  
  # ── Posterior summary ─────────────────────────────────────────────────────────
  pred_df <- tibble(
    p_to      = p_to_raw,
    p_from    = p_from_raw,
    psi_mean  = apply(psi_draws, 2, mean),
    psi_lower = apply(psi_draws, 2, quantile, 0.025),
    psi_upper = apply(psi_draws, 2, quantile, 0.975)
  )
  
  # ── Observed range for the "to" stage ────────────────────────────────────────
  # Compex is implicit (no homerange predictor); approximate as 1 - sum of others
  obs_raw <- if (to == "compex") {
    1 - get_hr_raw("pcnt_standinit", species_name, hr_sc$standinit) -
      get_hr_raw("pcnt_thin",      species_name, hr_sc$thin)      -
      get_hr_raw("pcnt_mature",    species_name, hr_sc$mature)
  } else {
    get_hr_raw(paste0("pcnt_", to), species_name, hr_sc[[to]])
  }
  obs_range <- quantile(obs_raw, c(0.05, 0.95))
  
  list(pred_df   = pred_df,
       obs_raw   = obs_raw,
       obs_range = obs_range,
       from      = from,
       to        = to,
       sp_name   = species_name)
}

# ── PLOT FUNCTION ─────────────────────────────────────────────────────────────
plot_hr_gradient <- function(result) {
  
  pred_df   <- result$pred_df
  obs_raw   <- result$obs_raw
  obs_range <- result$obs_range
  from      <- result$from
  to        <- result$to
  sp_name   <- result$sp_name
  
  stage_labels <- c(compex    = "Competitive exclusion",
                    standinit = "Stand initiation",
                    mature    = "Mature",
                    thin      = "Thinned")
  
  x_max      <- max(obs_raw)
  to_label   <- stage_labels[to]
  from_label <- stage_labels[from]
  
  ggplot(pred_df, aes(x = p_to)) +
    annotate("rect",
             xmin = obs_range["95%"], xmax = x_max,
             ymin = 0, ymax = 1, fill = "grey80", alpha = 0.4) +
    geom_ribbon(aes(ymin = psi_lower, ymax = psi_upper),
                alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = psi_mean), colour = "steelblue", linewidth = 1) +
    geom_rug(data = tibble(x = obs_raw), aes(x = x),
             inherit.aes = FALSE, alpha = 0.3, length = unit(0.03, "npc")) +
    scale_x_continuous(
      limits   = c(0, x_max),
      labels   = scales::percent,
      sec.axis = sec_axis(~ 1 - ., name = paste0("% ", from_label),
                          labels = scales::percent)
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      x        = paste0("% ", to_label),
      y        = "Predicted occupancy probability",
      title    = sp_name,
      subtitle = paste(
        paste0("Marginal occupancy: ", from_label, " \u2013 ", to_label, " gradient"),
        "Other stages fixed at 0% | Point covariates at grand mean | Plot covariates at stage-conditional means",
        "Shaded: extrapolation beyond 95th percentile | Rug: observed homerange values",
        sep = "\n"
      )
    ) +
    theme_classic()
}

# ── USAGE ─────────────────────────────────────────────────────────────────────
sp <- "brown creeper"

plot_hr_gradient(predict_hr_gradient(sp, from = "compex",    to = "mature"))
plot_hr_gradient(predict_hr_gradient(sp, from = "compex",    to = "standinit"))
plot_hr_gradient(predict_hr_gradient(sp, from = "compex",    to = "thin"))
plot_hr_gradient(predict_hr_gradient(sp, from = "standinit", to = "mature"))
plot_hr_gradient(predict_hr_gradient(sp, from = "thin",      to = "standinit"))

##########################################################################################

library(tidyverse)
library(ggtern)

# ── SETUP (run once) ──────────────────────────────────────────────────────────
sims    <- model_data$msom$sims.list
n_draws <- dim(sims$u)[1]

stage_idx_map <- c(compex = 1, standinit = 2, mature = 3, thin = 4)

# Homerange predictor slot for each non-compex stage
# (compex is reference — no homerange predictor)
stage_hr_map <- c(standinit = "alpha_homerange1",
                  thin      = "alpha_homerange2",
                  mature    = "alpha_homerange3")

# Stage-conditional plot covariate means (globally scaled, computed once)
stage_vec      <- model_data$stages$stage_idx
n_rows         <- nrow(model_data$param_alpha_data$param_alpha_plot_data$data[[1]])
stage_expanded <- rep(stage_vec, times = n_rows / length(stage_vec))

plot_data <- model_data$param_alpha_data$param_alpha_plot_data
x_plot    <- lapply(plot_data$data, as.vector)

stage_plot_means <- function(stage_code) {
  idx <- stage_expanded == stage_code
  sapply(x_plot, \(x) mean(x[idx], na.rm = TRUE))
}

mean_plot <- lapply(stage_idx_map, stage_plot_means)
names(mean_plot) <- names(stage_idx_map)

# Helpers
get_hr_scaling <- function(hr_param_name, species_name) {
  x <- model_data$param_alpha_data$param_alpha_homerange_data |>
    filter(name == hr_param_name) |>
    pull(data) |> _[[1]]
  list(center = attr(x[[species_name]], "scaled:center"),
       scale  = attr(x[[species_name]], "scaled:scale"))
}

scale_proportion <- function(p, sc) (p - sc$center) / sc$scale

# ── TERNARY PREDICTION FUNCTION ───────────────────────────────────────────────
# axes:  character vector of exactly 3 stage names to vary on the ternary surface
#        the 4th stage is automatically identified and fixed at 0
# step:  grid resolution (e.g. 0.05 = 5% steps -> ~231 points)
predict_ternary <- function(species_name,
                            axes = c("compex", "standinit", "mature"),
                            step = 0.02) {
  
  stopifnot(length(axes) == 3)
  stopifnot(all(axes %in% names(stage_idx_map)))
  stopifnot(length(unique(axes)) == 3)
  
  sp_idx      <- which(model_data$species == species_name)
  fixed_stage <- setdiff(names(stage_idx_map), axes)  # the one stage held at 0
  
  # Homerange scaling for all non-compex stages
  hr_sc <- list(
    standinit = get_hr_scaling("pcnt_standinit", species_name),
    thin      = get_hr_scaling("pcnt_thin",      species_name),
    mature    = get_hr_scaling("pcnt_mature",     species_name)
  )
  
  # ── Ternary grid ─────────────────────────────────────────────────────────────
  # Build over axes[1] and axes[2]; axes[3] = 1 - axes[1] - axes[2]
  grid_raw <- expand.grid(
    p1 = seq(0, 1, by = step),
    p2 = seq(0, 1, by = step)
  ) |>
    filter(p1 + p2 <= 1 + step * 0.01) |>
    mutate(
      p1 = round(p1, 10),
      p2 = round(p2, 10),
      p3 = round(1 - p1 - p2, 10)
    ) |>
    as_tibble()
  
  # Name columns by actual stage names
  names(grid_raw) <- c(paste0("p_", axes[1]),
                       paste0("p_", axes[2]),
                       paste0("p_", axes[3]))
  
  # Add fixed stage column (always 0)
  grid_raw[[paste0("p_", fixed_stage)]] <- 0
  n_grid <- nrow(grid_raw)
  
  message(species_name, ": ", n_grid, " grid points | axes: ",
          paste(axes, collapse = ", "), " | fixed: ", fixed_stage, " = 0")
  
  # ── Posterior slices ─────────────────────────────────────────────────────────
  u_sp   <- sims$u[, sp_idx]
  get_ap <- function(stage, k) sims[[paste0("alpha_plot", k)]][, stage_idx_map[stage], sp_idx]
  
  ap <- setNames(
    lapply(names(stage_idx_map), \(s) lapply(1:3, \(k) get_ap(s, k))),
    names(stage_idx_map)
  )
  
  ahr <- list(
    standinit = sims$alpha_homerange1[, sp_idx],
    thin      = sims$alpha_homerange2[, sp_idx],
    mature    = sims$alpha_homerange3[, sp_idx]
  )
  
  # ── Homerange contribution [n_draws x n_grid] ────────────────────────────────
  # For each non-compex stage:
  #   - if on an axis: p varies across grid -> outer() -> [n_draws x n_grid]
  #   - if fixed at 0: scalar -> [n_draws] vector, recycles across columns
  hr_contrib <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (s in names(stage_hr_map)) {  # standinit, thin, mature
    p_raw <- grid_raw[[paste0("p_", s)]]
    
    p_s_scaled <- if (length(unique(p_raw)) > 1) {
      scale_proportion(p_raw, hr_sc[[s]])   # vector [n_grid]
    } else {
      scale_proportion(0, hr_sc[[s]])        # scalar
    }
    
    contribution <- if (length(p_s_scaled) > 1) {
      outer(ahr[[s]], p_s_scaled)    # [n_draws x n_grid]
    } else {
      ahr[[s]] * p_s_scaled          # [n_draws] vector, recycles across columns
    }
    
    hr_contrib <- hr_contrib + contribution
  }
  
  # ── Plot contributions [n_draws vectors] ─────────────────────────────────────
  plot_contrib <- function(stage) {
    ap[[stage]][[1]] * mean_plot[[stage]][1] +
      ap[[stage]][[2]] * mean_plot[[stage]][2] +
      ap[[stage]][[3]] * mean_plot[[stage]][3]
  }
  
  # ── Option B: weighted average across active stages ───────────────────────────
  # Fixed stage has weight 0 everywhere so is skipped automatically
  psi_draws <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (s in names(stage_idx_map)) {
    p_s <- grid_raw[[paste0("p_", s)]]
    if (all(p_s == 0)) next   # fixed stage — skip
    
    eta_s <- hr_contrib + (u_sp + plot_contrib(s))
    w_s   <- matrix(p_s, nrow = n_draws, ncol = n_grid, byrow = TRUE)
    psi_draws <- psi_draws + w_s * plogis(eta_s)
  }
  
  # ── Posterior summary ─────────────────────────────────────────────────────────
  grid_raw |>
    mutate(
      psi_mean  = apply(psi_draws, 2, mean),
      psi_lower = apply(psi_draws, 2, quantile, 0.025),
      psi_upper = apply(psi_draws, 2, quantile, 0.975),
      ci_width  = psi_upper - psi_lower,
      sp_name   = species_name
    )
}

# ── TERNARY PLOT FUNCTION ─────────────────────────────────────────────────────
plot_ternary <- function(pred_df,
                         axes  = c("compex", "standinit", "mature"),
                         value = "psi_mean") {
  
  sp_name     <- unique(pred_df$sp_name)
  fixed_stage <- setdiff(names(stage_idx_map), axes)
  
  stage_labels <- c(compex    = "Competitive exclusion",
                    standinit = "Stand initiation",
                    mature    = "Mature",
                    thin      = "Thinned")
  
  value_label <- switch(value,
                        psi_mean = "Mean\noccupancy\nprobability",
                        ci_width = "95% CI\nwidth"
  )
  
  ggtern(pred_df,
         aes(x = .data[[paste0("p_", axes[1])]],
             y = .data[[paste0("p_", axes[2])]],
             z = .data[[paste0("p_", axes[3])]])) +
    geom_point(aes(color = .data[[value]]), shape = 16, size = 1.8) +
    scale_color_viridis_c(
      option = "viridis",
      limits = if (value == "psi_mean") c(0, 1) else NULL,
      labels = scales::percent,
      name   = value_label
    ) +
    labs(
      title    = sp_name,
      subtitle = paste0(
        "p(", stage_labels[fixed_stage], ") = 0 | ",
        "Point covariates at grand mean | Plot covariates at stage-conditional means"
      ),
      x = stage_labels[axes[1]],
      y = stage_labels[axes[2]],
      z = stage_labels[axes[3]]
    ) +
    theme_bw() +
    theme_showarrows()
}

# ── USAGE ─────────────────────────────────────────────────────────────────────
sp <- "golden-crowned kinglet"

# Default: compex, standinit, mature (thin fixed at 0)
pred <- predict_ternary(sp, axes = c("standinit", "compex", "mature"), step = 0.02)
plot_ternary(pred, axes = c("standinit", "compex", "mature"), value = "psi_mean")
plot_ternary(pred, axes = c("standinit", "compex", "mature"), value = "ci_width")

# Alternative: standinit, thin, mature (compex fixed at 0)
pred2 <- predict_ternary(sp, axes = c("standinit", "thin", "mature"), step = 0.02)
plot_ternary(pred2, axes = c("standinit", "thin", "mature"), value = "psi_mean")

# All species, default axes
all_preds <- lapply(model_data$species, predict_ternary,
                    axes = c("compex", "standinit", "mature"), step = 0.05)
plots <- lapply(all_preds, plot_ternary,
                axes = c("compex", "standinit", "mature"), value = "psi_mean")

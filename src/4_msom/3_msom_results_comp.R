####################################################################################
# Marginal plots of varying homerange stage composition
#
# CONFIG:
path_msom = "data/cache/models/V3_msom_pcnt_fp_fp_all.rds"
#
# INPUT:
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
####################################################################################

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

(msom_summary = model_data$msom_summary)
(msom = model_data$msom)
(groups = model_data$groups %>% arrange(common_name))
(sites = model_data$sites)
(species = model_data$species)
(stages = model_data$stages)

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

# ── PREDICTION FUNCTION (returns draws + summary) ────────────────────────────
predict_hr_gradient <- function(species_name,
                                from = "compex",
                                to   = "thin") {
  
  stopifnot(from %in% names(stage_idx_map), to %in% names(stage_idx_map))
  stopifnot(from != to)
  
  sp_idx     <- which(model_data$species == species_name)
  n_grid     <- 101
  p_to_raw   <- seq(0, 1, length.out = n_grid)
  p_from_raw <- 1 - p_to_raw
  
  hr_sc <- list(
    standinit = get_hr_scaling("pcnt_standinit", species_name),
    thin      = get_hr_scaling("pcnt_thin",      species_name),
    mature    = get_hr_scaling("pcnt_mature",     species_name)
  )
  
  u_sp   <- sims$u[, sp_idx]
  get_ap <- function(stage, k) sims[[paste0("alpha_plot", k)]][, stage_idx_map[stage], sp_idx]
  ap_from <- lapply(1:3, \(k) get_ap(from, k))
  ap_to   <- lapply(1:3, \(k) get_ap(to,   k))
  
  ahr <- list(
    standinit = sims$alpha_homerange1[, sp_idx],
    thin      = sims$alpha_homerange2[, sp_idx],
    mature    = sims$alpha_homerange3[, sp_idx]
  )
  
  hr_contrib <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (s in c("standinit", "thin", "mature")) {
    p_s_scaled <- if (s == to) {
      scale_proportion(p_to_raw,   hr_sc[[s]])
    } else if (s == from) {
      scale_proportion(p_from_raw, hr_sc[[s]])
    } else {
      scale_proportion(0, hr_sc[[s]])
    }
    
    contribution <- if (length(p_s_scaled) > 1) {
      outer(ahr[[s]], p_s_scaled)
    } else {
      ahr[[s]] * p_s_scaled
    }
    
    hr_contrib <- hr_contrib + contribution
  }
  
  plot_contrib <- function(stage, ap_list) {
    ap_list[[1]] * mean_plot[[stage]][1] +
      ap_list[[2]] * mean_plot[[stage]][2] +
      ap_list[[3]] * mean_plot[[stage]][3]
  }
  pc_from <- plot_contrib(from, ap_from)
  pc_to   <- plot_contrib(to,   ap_to)
  
  eta_from <- hr_contrib + (u_sp + pc_from)
  eta_to   <- hr_contrib + (u_sp + pc_to)
  
  p_from_mat <- matrix(p_from_raw, nrow = n_draws, ncol = n_grid, byrow = TRUE)
  p_to_mat   <- matrix(p_to_raw,   nrow = n_draws, ncol = n_grid, byrow = TRUE)
  
  psi_draws <- p_from_mat * plogis(eta_from) +
    p_to_mat   * plogis(eta_to)
  
  pred_df <- tibble(
    p_to      = p_to_raw,
    p_from    = p_from_raw,
    psi_mean  = apply(psi_draws, 2, mean),
    psi_lower = apply(psi_draws, 2, quantile, 0.025),
    psi_upper = apply(psi_draws, 2, quantile, 0.975)
  )
  
  obs_raw <- if (to == "compex") {
    1 - get_hr_raw("pcnt_standinit", species_name, hr_sc$standinit) -
      get_hr_raw("pcnt_thin",      species_name, hr_sc$thin)      -
      get_hr_raw("pcnt_mature",    species_name, hr_sc$mature)
  } else {
    get_hr_raw(paste0("pcnt_", to), species_name, hr_sc[[to]])
  }
  obs_range <- quantile(obs_raw, c(0.05, 0.95))
  
  list(pred_df    = pred_df,
       psi_draws  = psi_draws,   # [n_draws x n_grid] — for community averaging
       p_to_raw   = p_to_raw,
       obs_raw    = obs_raw,
       obs_range  = obs_range,
       from       = from,
       to         = to,
       sp_name    = species_name)
}

# ── SINGLE-SPECIES PLOT FUNCTION ─────────────────────────────────────────────
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

# ── COMMUNITY GRADIENT PREDICTION ────────────────────────────────────────────
# Returns both per-species summaries (for ghost lines) and a community mean
# (for the bold line). Averaging is done on draws before summarising so that
# uncertainty in the community mean is propagated correctly.
#
# species_subset: character vector; NULL = all species
predict_hr_gradient_community <- function(from           = "compex",
                                          to             = "thin",
                                          species_subset = NULL) {
  
  sp_list <- if (is.null(species_subset)) model_data$species else species_subset
  n_sp    <- length(sp_list)
  n_grid  <- 101
  
  message("Community gradient: ", n_sp, " species | ", from, " -> ", to)
  
  # Collect per-species draws and summaries
  psi_sum      <- matrix(0, nrow = n_draws, ncol = n_grid)
  species_dfs  <- vector("list", n_sp)
  
  for (k in seq_along(sp_list)) {
    sp <- sp_list[k]
    message("  ", sp)
    res <- predict_hr_gradient(sp, from = from, to = to)
    psi_sum        <- psi_sum + res$psi_draws
    species_dfs[[k]] <- res$pred_df |> mutate(sp_name = sp)
  }
  
  # Community mean draws -> summary
  psi_comm_draws <- psi_sum / n_sp
  p_to_raw       <- seq(0, 1, length.out = n_grid)
  
  comm_df <- tibble(
    p_to      = p_to_raw,
    p_from    = 1 - p_to_raw,
    psi_mean  = apply(psi_comm_draws, 2, mean),
    psi_lower = apply(psi_comm_draws, 2, quantile, 0.025),
    psi_upper = apply(psi_comm_draws, 2, quantile, 0.975),
    sp_name   = paste0("Community mean (", n_sp, " species)")
  )
  
  list(
    species_df = bind_rows(species_dfs),
    comm_df    = comm_df,
    from       = from,
    to         = to,
    n_sp       = n_sp
  )
}

# ── MULTI-SPECIES GRADIENT PLOT ───────────────────────────────────────────────
# Plots per-species lines as translucent ghosts and the community mean as a
# bold line with a credible interval ribbon.
#
# result:          output of predict_hr_gradient_community()
# species_subset:  optionally highlight a named subset in a different color
# show_ci:         if TRUE draws the 95% CI ribbon on the community mean line
plot_hr_gradient_community <- function(result,
                                       show_ci = TRUE) {
  
  species_df <- result$species_df
  comm_df    <- result$comm_df
  from       <- result$from
  to         <- result$to
  
  stage_labels <- c(compex    = "Competitive exclusion",
                    standinit = "Stand initiation",
                    mature    = "Mature",
                    thin      = "Thinned")
  
  to_label   <- stage_labels[to]
  from_label <- stage_labels[from]
  
  p <- ggplot() +
    # Ghost lines — one per species
    geom_line(
      data    = species_df,
      mapping = aes(x = p_to, y = psi_mean, group = sp_name),
      color   = "steelblue", alpha = 0.25, linewidth = 0.4
    )
  
  # Optional CI ribbon on community mean
  if (show_ci) {
    p <- p +
      geom_ribbon(
        data    = comm_df,
        mapping = aes(x = p_to, ymin = psi_lower, ymax = psi_upper),
        fill    = "steelblue", alpha = 0.2
      )
  }
  
  p +
    # Community mean — bold foreground line
    geom_line(
      data    = comm_df,
      mapping = aes(x = p_to, y = psi_mean),
      color   = "steelblue", linewidth = 1.2
    ) +
    scale_x_continuous(
      labels   = scales::percent,
      sec.axis = sec_axis(~ 1 - ., name = paste0("% ", from_label),
                          labels = scales::percent)
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      x        = paste0("% ", to_label),
      y        = "Predicted occupancy probability",
      title    = comm_df$sp_name,
      subtitle = paste(
        paste0("Marginal occupancy: ", from_label, " \u2013 ", to_label, " gradient"),
        "Ghost lines: individual species | Bold line: community mean \u00B1 95% CI",
        "Other stages fixed at 0% | Point covariates at grand mean | Plot covariates at stage-conditional means",
        sep = "\n"
      )
    ) +
    theme_classic()
}

# ── USAGE ─────────────────────────────────────────────────────────────────────

# Single species gradient (unchanged)
plot_hr_gradient(predict_hr_gradient("hairy woodpecker", from = "compex", to = "thin"))

# Community gradient — all species
comm_result <- predict_hr_gradient_community(from = "compex", to = "thin")
plot_hr_gradient_community(comm_result)

# Community gradient — focal subset
focal_species <- c("pileated woodpecker", "hairy woodpecker", "brown creeper",
                   "red-breasted nuthatch", "pacific wren", "varied thrush")
comm_focal <- predict_hr_gradient_community(from           = "compex",
                                            to             = "thin",
                                            species_subset = focal_species)
plot_hr_gradient_community(comm_focal)

# Without CI ribbon
plot_hr_gradient_community(comm_result, show_ci = FALSE)

for (sp in species) {
  p = plot_hr_gradient(predict_hr_gradient(sp, from = "compex", to = "thin"))
  print(p)
}

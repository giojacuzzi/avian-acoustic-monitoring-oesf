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

library(ggtern)

# ── SETUP (run once) ──────────────────────────────────────────────────────────
sims    <- model_data$msom$sims.list
n_draws <- dim(sims$u)[1]

stage_idx_map <- c(compex = 1, standinit = 2, mature = 3, thin = 4)

stage_hr_map <- c(standinit = "alpha_homerange1",
                  thin      = "alpha_homerange2",
                  mature    = "alpha_homerange3")

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

get_hr_scaling <- function(hr_param_name, species_name) {
  x <- model_data$param_alpha_data$param_alpha_homerange_data |>
    filter(name == hr_param_name) |>
    pull(data) |> _[[1]]
  list(center = attr(x[[species_name]], "scaled:center"),
       scale  = attr(x[[species_name]], "scaled:scale"))
}

scale_proportion <- function(p, sc) (p - sc$center) / sc$scale

# ── GRID BUILDER ──────────────────────────────────────────────────────────────
make_ternary_grid <- function(axes, fixed_stage, step) {
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
  
  names(grid_raw) <- c(paste0("p_", axes[1]),
                       paste0("p_", axes[2]),
                       paste0("p_", axes[3]))
  grid_raw[[paste0("p_", fixed_stage)]] <- 0
  grid_raw
}

# ── SINGLE-SPECIES PREDICTION (returns raw draws) ─────────────────────────────
predict_species_draws <- function(species_name, axes, grid_raw) {
  
  sp_idx <- which(model_data$species == species_name)
  n_grid <- nrow(grid_raw)
  
  hr_sc <- list(
    standinit = get_hr_scaling("pcnt_standinit", species_name),
    thin      = get_hr_scaling("pcnt_thin",      species_name),
    mature    = get_hr_scaling("pcnt_mature",     species_name)
  )
  
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
  
  hr_contrib <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (s in names(stage_hr_map)) {
    p_raw <- grid_raw[[paste0("p_", s)]]
    
    p_s_scaled <- if (length(unique(p_raw)) > 1) {
      scale_proportion(p_raw, hr_sc[[s]])
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
  
  plot_contrib <- function(stage) {
    ap[[stage]][[1]] * mean_plot[[stage]][1] +
      ap[[stage]][[2]] * mean_plot[[stage]][2] +
      ap[[stage]][[3]] * mean_plot[[stage]][3]
  }
  
  psi_draws <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (s in names(stage_idx_map)) {
    p_s <- grid_raw[[paste0("p_", s)]]
    if (all(p_s == 0)) next
    eta_s <- hr_contrib + (u_sp + plot_contrib(s))
    w_s   <- matrix(p_s, nrow = n_draws, ncol = n_grid, byrow = TRUE)
    psi_draws <- psi_draws + w_s * plogis(eta_s)
  }
  
  psi_draws
}

# ── SINGLE-SPECIES TERNARY PREDICTION ─────────────────────────────────────────
predict_ternary <- function(species_name,
                            axes = c("compex", "standinit", "mature"),
                            step = 0.02) {
  
  stopifnot(length(axes) == 3)
  stopifnot(all(axes %in% names(stage_idx_map)))
  stopifnot(length(unique(axes)) == 3)
  
  fixed_stage <- setdiff(names(stage_idx_map), axes)
  grid_raw    <- make_ternary_grid(axes, fixed_stage, step)
  message(species_name, ": ", nrow(grid_raw), " grid points")
  
  psi_draws <- predict_species_draws(species_name, axes, grid_raw)
  
  grid_raw |>
    mutate(
      psi_mean  = apply(psi_draws, 2, mean),
      psi_lower = apply(psi_draws, 2, quantile, 0.025),
      psi_upper = apply(psi_draws, 2, quantile, 0.975),
      ci_width  = psi_upper - psi_lower,
      sp_name   = species_name
    )
}

# ── COMMUNITY TERNARY PREDICTION ──────────────────────────────────────────────
predict_ternary_community <- function(axes           = c("compex", "standinit", "mature"),
                                      step           = 0.02,
                                      species_subset = NULL) {
  
  stopifnot(length(axes) == 3)
  stopifnot(all(axes %in% names(stage_idx_map)))
  stopifnot(length(unique(axes)) == 3)
  
  fixed_stage <- setdiff(names(stage_idx_map), axes)
  grid_raw    <- make_ternary_grid(axes, fixed_stage, step)
  n_grid      <- nrow(grid_raw)
  
  sp_list <- if (is.null(species_subset)) model_data$species else species_subset
  n_sp    <- length(sp_list)
  
  message("Community prediction: ", n_sp, " species | ", n_grid, " grid points")
  
  psi_sum <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (sp in sp_list) {
    message("  ", sp)
    psi_sum <- psi_sum + predict_species_draws(sp, axes, grid_raw)
  }
  
  psi_mean_draws <- psi_sum / n_sp
  
  grid_raw |>
    mutate(
      psi_mean  = apply(psi_mean_draws, 2, mean),
      psi_lower = apply(psi_mean_draws, 2, quantile, 0.025),
      psi_upper = apply(psi_mean_draws, 2, quantile, 0.975),
      ci_width  = psi_upper - psi_lower,
      sp_name   = paste0("Community mean (", n_sp, " species)")
    )
}

# ── SPECIES MAXIMA DENSITY PREDICTION ─────────────────────────────────────────
# For each species, finds the grid point with the highest posterior mean
# occupancy. Then counts how many species share each grid point as their
# maximum. Returns a grid data frame where n_species is the count at each point
# (0 for grid points that are no species' maximum).
#
# species_subset: character vector of species to include; NULL = all species
predict_ternary_maxima <- function(axes           = c("compex", "standinit", "mature"),
                                   step           = 0.05,
                                   species_subset = NULL) {
  
  stopifnot(length(axes) == 3)
  stopifnot(all(axes %in% names(stage_idx_map)))
  stopifnot(length(unique(axes)) == 3)
  
  fixed_stage <- setdiff(names(stage_idx_map), axes)
  grid_raw    <- make_ternary_grid(axes, fixed_stage, step)
  n_grid      <- nrow(grid_raw)
  
  sp_list <- if (is.null(species_subset)) model_data$species else species_subset
  n_sp    <- length(sp_list)
  
  message("Species maxima: ", n_sp, " species | ", n_grid, " grid points")
  
  # For each species, record which grid row has the highest posterior mean
  # occupancy. Near-ties (within 0.5% of max) are collapsed to their centroid
  # row — the closest grid point to that centroid is used.
  maxima_rows <- integer(n_sp)
  
  for (k in seq_along(sp_list)) {
    sp <- sp_list[k]
    message("  ", sp)
    
    psi_draws  <- predict_species_draws(sp, axes, grid_raw)
    psi_mean   <- colMeans(psi_draws)
    
    # Find centroid of near-tied maxima, snap to nearest grid point
    tied_idx   <- which(psi_mean >= max(psi_mean) - 0.005)
    centroid   <- colMeans(grid_raw[tied_idx, paste0("p_", axes), drop = FALSE])
    
    # Nearest grid point to centroid (Euclidean distance in ternary space)
    dists         <- rowSums(sweep(as.matrix(grid_raw[, paste0("p_", axes)]),
                                   2, centroid, "-")^2)
    maxima_rows[k] <- which.min(dists)
  }
  
  # Build per-species maxima table (useful for downstream labelling)
  maxima_df <- grid_raw[maxima_rows, ] |>
    mutate(sp_name = sp_list)
  
  # Count how many species share each grid point as their maximum
  counts <- tabulate(maxima_rows, nbins = n_grid)
  
  grid_raw |>
    mutate(n_species = counts)
}

# ── CONTOUR HELPER ────────────────────────────────────────────────────────────
# Computes contour lines directly from the prediction grid using base R
# contourLines(), guaranteeing alignment with the color fill values.
compute_ternary_contours <- function(pred_df, axes, breaks) {
  
  x_col <- paste0("p_", axes[1])
  y_col <- paste0("p_", axes[2])
  z_col <- paste0("p_", axes[3])
  
  x_vals <- sort(unique(round(pred_df[[x_col]], 10)))
  y_vals <- sort(unique(round(pred_df[[y_col]], 10)))
  
  z_mat <- matrix(NA_real_, nrow = length(x_vals), ncol = length(y_vals),
                  dimnames = list(x_vals, y_vals))
  
  for (i in seq_len(nrow(pred_df))) {
    xi <- as.character(round(pred_df[[x_col]][i], 10))
    yi <- as.character(round(pred_df[[y_col]][i], 10))
    z_mat[xi, yi] <- pred_df$psi_mean[i]
  }
  
  breaks <- breaks[breaks > min(pred_df$psi_mean, na.rm = TRUE) &
                     breaks < max(pred_df$psi_mean, na.rm = TRUE)]
  if (length(breaks) == 0) return(NULL)
  
  cl <- contourLines(x = x_vals, y = y_vals, z = z_mat, levels = breaks)
  if (length(cl) == 0) return(NULL)
  
  map_dfr(seq_along(cl), function(k) {
    tibble(
      !!x_col    := cl[[k]]$x,
      !!y_col    := cl[[k]]$y,
      !!z_col    := pmax(0, 1 - cl[[k]]$x - cl[[k]]$y),
      level      = cl[[k]]$level,
      contour_id = k
    )
  })
}

# ── TERNARY PLOT FUNCTION ──────────────────────────────────────────────────────
plot_ternary <- function(pred_df,
                         axes           = c("compex", "standinit", "mature"),
                         value          = "psi_mean",
                         scale_limits   = c("full", "species"),
                         annotate_max   = TRUE,
                         contour_breaks = NULL,
                         no_mask        = TRUE) {
  
  scale_limits <- match.arg(scale_limits)
  sp_name      <- unique(pred_df$sp_name)
  fixed_stage  <- setdiff(names(stage_idx_map), axes)
  
  stage_labels <- c(compex    = "Competitive exclusion",
                    standinit = "Stand initiation",
                    mature    = "Mature",
                    thin      = "Thinned")
  
  stage_abbrev <- c(compex    = "CE",
                    standinit = "SI",
                    mature    = "Mat",
                    thin      = "Thin")
  
  value_label <- switch(value,
                        psi_mean = "Mean\noccupancy\nprobability",
                        ci_width = "95% CI\nwidth"
  )
  
  col_limits <- if (value == "psi_mean") {
    switch(scale_limits,
           full    = c(0, 1),
           species = c(0, max(pred_df$psi_mean))
    )
  } else {
    NULL
  }
  
  max_row <- pred_df |>
    filter(psi_mean >= max(psi_mean) - 0.005) |>
    summarise(
      across(starts_with("p_"), mean),
      psi_mean  = max(psi_mean),
      psi_lower = mean(psi_lower),
      psi_upper = mean(psi_upper)
    )
  
  max_point    <- max_row
  caption_text <- if (annotate_max && value == "psi_mean") {
    paste0(
      "Max \u03C8: ",
      scales::percent(max_point$psi_mean, accuracy = 0.1), " at ",
      scales::percent(max_point[[paste0("p_", axes[1])]], accuracy = 1),
      " ", stage_abbrev[axes[1]], " / ",
      scales::percent(max_point[[paste0("p_", axes[2])]], accuracy = 1),
      " ", stage_abbrev[axes[2]], " / ",
      scales::percent(max_point[[paste0("p_", axes[3])]], accuracy = 1),
      " ", stage_abbrev[axes[3]]
    )
  } else {
    NULL
  }
  
  p <- ggtern(pred_df,
              aes(x = .data[[paste0("p_", axes[1])]],
                  y = .data[[paste0("p_", axes[2])]],
                  z = .data[[paste0("p_", axes[3])]])) +
    geom_point(aes(color = .data[[value]]), shape = 16, size = 1.8) +
    scale_color_viridis_c(
      option = "viridis",
      limits = col_limits,
      labels = scales::percent,
      name   = value_label
    ) +
    labs(
      title    = sp_name,
      subtitle = paste0(
        "p(", stage_labels[fixed_stage], ") = 0 | ",
        "Point covariates at grand mean | Plot covariates at stage-conditional means"
      ),
      x       = stage_labels[axes[1]],
      y       = stage_labels[axes[2]],
      z       = stage_labels[axes[3]],
      caption = caption_text
    ) +
    theme_bw() +
    theme_showarrows()
  
  if (no_mask) p <- p + theme_nomask()
  
  # ── Contour lines ─────────────────────────────────────────────────────────────
  if (!is.null(contour_breaks) && value == "psi_mean") {
    contour_df <- compute_ternary_contours(pred_df, axes, contour_breaks)
    
    if (!is.null(contour_df)) {
      contour_labels_df <- contour_df |>
        group_by(contour_id, level) |>
        slice(floor(n() / 2)) |>
        ungroup() |>
        mutate(label = scales::percent(level, accuracy = 1))
      
      p <- p +
        geom_path(
          data      = contour_df,
          mapping   = aes(x     = .data[[paste0("p_", axes[1])]],
                          y     = .data[[paste0("p_", axes[2])]],
                          z     = .data[[paste0("p_", axes[3])]],
                          group = contour_id),
          color     = "white",
          linewidth = 0.6
        ) +
        geom_text(
          data     = contour_labels_df,
          mapping  = aes(x     = .data[[paste0("p_", axes[1])]],
                         y     = .data[[paste0("p_", axes[2])]],
                         z     = .data[[paste0("p_", axes[3])]],
                         label = label),
          color    = "white",
          size     = 2.8,
          fontface = "bold"
        )
    }
  }
  
  # ── Maximum annotation ────────────────────────────────────────────────────────
  if (annotate_max && value == "psi_mean") {
    p <- p +
      geom_point(data = max_point, shape = 21, size = 5,
                 fill = "white", color = "black", stroke = 1.5) +
      geom_point(data = max_point, shape = 16, size = 2, color = "black")
  }
  
  p
}

# ── SPECIES MAXIMA DENSITY PLOT ───────────────────────────────────────────────
# Plots the number of species whose occupancy maximum falls at each grid point.
# Grid points with zero species are shown in a neutral grey; occupied points
# use the viridis scale. Point size scales with n_species to reinforce the
# count signal visually.
plot_ternary_maxima <- function(maxima_df,
                                axes           = c("compex", "standinit", "mature"),
                                species_subset = NULL,
                                no_mask        = TRUE) {
  
  fixed_stage  <- setdiff(names(stage_idx_map), axes)
  sp_list      <- if (is.null(species_subset)) model_data$species else species_subset
  n_sp         <- length(sp_list)
  
  stage_labels <- c(compex    = "Competitive exclusion",
                    standinit = "Stand initiation",
                    mature    = "Mature",
                    thin      = "Thinned")
  
  # Split into zero and non-zero counts for separate layers
  df_zero    <- maxima_df |> filter(n_species == 0)
  df_nonzero <- maxima_df |> filter(n_species  > 0)
  
  ggtern(maxima_df,
         aes(x = .data[[paste0("p_", axes[1])]],
             y = .data[[paste0("p_", axes[2])]],
             z = .data[[paste0("p_", axes[3])]])) +
    # Empty grid points — neutral backdrop
    geom_point(data  = df_zero,
               shape = 16, size = 1.5, color = "grey88") +
    # Occupied grid points — colored and sized by count
    geom_point(data    = df_nonzero,
               mapping = aes(color = n_species, size = n_species),
               shape   = 16) +
    scale_color_viridis_c(
      option = "viridis",
      breaks = scales::breaks_pretty(),
      name   = "Number\nof species"
    ) +
    scale_size_continuous(
      range  = c(2, 8),
      breaks = scales::breaks_pretty(),
      name   = "Number\nof species",
      guide  = "legend"
    ) +
    # Merge color and size legends into one
    guides(color = guide_legend(), size = guide_legend()) +
    labs(
      title    = paste0("Species occupancy maxima (n = ", n_sp, ")"),
      subtitle = paste0(
        "p(", stage_labels[fixed_stage], ") = 0 | ",
        "Each point = number of species whose predicted occupancy peaks at that composition"
      ),
      x = stage_labels[axes[1]],
      y = stage_labels[axes[2]],
      z = stage_labels[axes[3]]
    ) +
    theme_bw() +
    theme_showarrows() +
    { if (no_mask) theme_nomask() else NULL }
}

# ── USAGE ─────────────────────────────────────────────────────────────────────
axes <- c("standinit", "compex", "mature")

# Single species
pred_sp <- predict_ternary("pileated woodpecker", axes = axes, step = 0.02)
plot_ternary(pred_sp, axes = axes, scale_limits = "species",
             annotate_max = TRUE, contour_breaks = NULL, no_mask = TRUE)
plot_ternary(pred_sp, axes = axes, scale_limits = "species",
             annotate_max = TRUE, contour_breaks = seq(0,1,0.1), no_mask = FALSE)

# Community mean with contours
pred_comm <- predict_ternary_community(axes = axes, step = 0.05)
plot_ternary(pred_comm, axes = axes, scale_limits = "species",
             annotate_max = TRUE, contour_breaks = c(0.25, 0.50, 0.75), no_mask = TRUE)

# Group mean with contours
focal_species = species_traits %>% filter(group_nest_ps == "tree") %>% pull(common_name)
# focal_species = species_traits %>% filter(group_habitat == "early seral") %>% pull(common_name)
focal_species = intersect(focal_species, species)
pred_focal = predict_ternary_community(axes = axes, step = 0.02, species_subset = focal_species)
plot_ternary(pred_focal, axes = axes, scale_limits = "species",
             annotate_max = TRUE, contour_breaks = c(0.25, 0.50, 0.75), no_mask = FALSE)

# Species maxima density — all species
maxima_all <- predict_ternary_maxima(axes = axes, step = 0.05)
plot_ternary_maxima(maxima_all, axes = axes, no_mask = TRUE)

# Species maxima density — focal subset
focal_species <- c("pileated woodpecker", "hairy woodpecker", "brown creeper",
                   "red-breasted nuthatch", "pacific wren", "varied thrush")
maxima_focal <- predict_ternary_maxima(axes = axes, step = 0.05,
                                       species_subset = focal_species)
plot_ternary_maxima(maxima_focal, axes = axes, species_subset = focal_species,
                    no_mask = TRUE)


# for (sp in species) {
#   p = plot_ternary(predict_ternary(sp, axes = axes, step = 0.02), axes = axes, scale_limits = "species", annotate_max = TRUE)
#   print(p)
# }

# compex max:
# steller's jay
# pacific-slope flycatcher
# hutton's vireo <
# black-throated gray warbler <
# band-tailed pigeon
# american robin (generalist)

# thinned max:


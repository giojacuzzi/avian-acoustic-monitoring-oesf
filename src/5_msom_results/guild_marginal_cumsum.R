####################################################################################
# Marginal plots of varying homerange stage composition
#
# CONFIG:
path_msom = "data/cache/models/V4_msom_V4_fp_all.rds"
#
# INPUT:
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
####################################################################################

source("src/global.R")

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

stage_idx_map <- c(compex = 1, standinit = 2, mature = 3, thin = 4)

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

# Stage colors
stage_colors <- c(
  standinit = "#E8A838",
  compex    = "#4CAF6B",
  thin      = "#8B5DB5",
  mature    = "#8B5E3C"
)

# ── HELPERS ───────────────────────────────────────────────────────────────────
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

interpolate_colors <- function(color_from, color_to, t) {
  ramp     <- colorRamp(c(color_from, color_to), space = "rgb")
  rgb_vals <- ramp(t)
  rgb(rgb_vals[, 1], rgb_vals[, 2], rgb_vals[, 3], maxColorValue = 255)
}

# ── PAIRWISE GRADIENT PREDICTION ──────────────────────────────────────────────
predict_hr_gradient <- function(species_name,
                                from   = "compex",
                                to     = "mature",
                                n_grid = 101) {
  
  stopifnot(from %in% names(stage_idx_map), to %in% names(stage_idx_map))
  stopifnot(from != to)
  
  sp_idx     <- which(model_data$species == species_name)
  p_to_raw   <- seq(0, 1, length.out = n_grid)
  p_from_raw <- 1 - p_to_raw
  
  hr_sc <- list(
    standinit = get_hr_scaling("pcnt_standinit", species_name),
    thin      = get_hr_scaling("pcnt_thin",      species_name),
    mature    = get_hr_scaling("pcnt_mature",     species_name)
  )
  
  u_sp    <- sims$u[, sp_idx]
  get_ap  <- function(stage, k) sims[[paste0("alpha_plot", k)]][, stage_idx_map[stage], sp_idx]
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
  
  list(pred_df   = pred_df,
       psi_draws = psi_draws,
       p_to_raw  = p_to_raw,
       obs_raw   = obs_raw,
       obs_range = obs_range,
       from      = from,
       to        = to,
       sp_name   = species_name)
}

# ── CONCATENATED GRADIENT PREDICTION ──────────────────────────────────────────
predict_gradient_concat <- function(species_name, n_grid = 501) {
  
  segments <- list(
    list(from = "standinit", to = "compex", label = "SI \u2192 Compex"),
    list(from = "compex",    to = "thin",   label = "Compex \u2192 Thin"),
    list(from = "thin",      to = "mature", label = "Thin \u2192 Mature")
  )
  
  all_dfs  <- vector("list", length(segments))
  seg_meta <- vector("list", length(segments))
  x_offset <- 0
  
  for (k in seq_along(segments)) {
    seg <- segments[[k]]
    res <- predict_hr_gradient(species_name, from = seg$from, to = seg$to,
                               n_grid = n_grid)
    
    x_vals <- seq(x_offset, x_offset + 1, length.out = n_grid)
    
    hex_colors <- interpolate_colors(
      color_from = stage_colors[seg$from],
      color_to   = stage_colors[seg$to],
      t          = res$pred_df$p_to
    )
    
    df <- res$pred_df |>
      mutate(
        x         = x_vals,
        segment   = k,
        seg_label = seg$label,
        from      = seg$from,
        to        = seg$to,
        p_to_raw  = p_to,
        hex_color = hex_colors
      )
    
    if (k < length(segments)) df <- df[-nrow(df), ]
    
    all_dfs[[k]] <- df
    
    seg_meta[[k]] <- list(
      seg_label = seg$label,
      from      = seg$from,
      to        = seg$to,
      segment   = k,
      x_min     = x_offset,
      x_max     = x_offset + 1,
      obs_lower = x_offset + as.numeric(res$obs_range["5%"]),
      obs_upper = x_offset + as.numeric(res$obs_range["95%"]),
      obs_max   = x_offset + min(max(res$obs_raw), 1)
    )
    
    x_offset <- x_offset + 1
  }
  
  # ── Rug: observed combination of from/to proportions ──────────────────────
  # For each segment, get raw observed proportions of both stages.
  # Position: p_to / (p_from + p_to) — relative balance of the two stages,
  #   normalized to [0,1] within the segment.
  # Alpha weight: p_from + p_to — how much of the homerange is accounted for
  #   by these two stages. Sites where both are near 0 (dominated by other
  #   stages) get faint ticks; sites where they jointly dominate get dark ticks.
  get_obs_raw <- function(stage, sp) {
    if (stage == "compex") {
      hr_sc_si  <- get_hr_scaling("pcnt_standinit", sp)
      hr_sc_th  <- get_hr_scaling("pcnt_thin",      sp)
      hr_sc_mat <- get_hr_scaling("pcnt_mature",     sp)
      1 - get_hr_raw("pcnt_standinit", sp, hr_sc_si) -
        get_hr_raw("pcnt_thin",      sp, hr_sc_th) -
        get_hr_raw("pcnt_mature",    sp, hr_sc_mat)
    } else {
      sc <- get_hr_scaling(paste0("pcnt_", stage), sp)
      get_hr_raw(paste0("pcnt_", stage), sp, sc)
    }
  }
  
  rug_df <- map_dfr(seq_along(seg_meta), \(k) {
    seg        <- seg_meta[[k]]
    p_from_obs <- get_obs_raw(seg$from, species_name)
    p_to_obs   <- get_obs_raw(seg$to,   species_name)
    combined   <- p_from_obs + p_to_obs
    
    tibble(
      x     = seg$x_min + if_else(combined > 0, p_to_obs / combined, NA_real_),
      alpha = pmin(combined, 1)
    ) |>
      filter(!is.na(x))
  })
  
  list(
    pred_df = bind_rows(all_dfs),
    meta_df = map_dfr(seg_meta, as_tibble),
    rug_df  = rug_df,
    sp_name = species_name
  )
}

# ── SINGLE-SPECIES CONCATENATED GRADIENT PLOT ─────────────────────────────────
plot_gradient_concat <- function(result) {
  
  pred_df <- result$pred_df
  meta_df <- result$meta_df
  rug_df  <- result$rug_df
  sp_name <- result$sp_name
  
  stage_labels <- c(
    standinit = "Stand initiation",
    compex    = "Competitive exclusion",
    thin      = "Thinned",
    mature    = "Mature"
  )
  
  # ── Line segments ─────────────────────────────────────────────────────────
  line_segs <- pred_df |>
    mutate(
      x_end      = lead(x),
      y_end      = lead(psi_mean),
      color_next = lead(hex_color),
      color_mid  = map2_chr(hex_color, color_next,
                            \(a, b) if (is.na(b)) NA_character_
                            else interpolate_colors(a, b, t = 0.5))
    ) |>
    filter(!is.na(x_end), !is.na(color_mid))
  
  # ── Ribbon strips ─────────────────────────────────────────────────────────
  ribbon_strips <- pred_df |>
    mutate(
      x_end      = lead(x),
      lower_end  = lead(psi_lower),
      upper_end  = lead(psi_upper),
      color_next = lead(hex_color),
      color_mid  = map2_chr(hex_color, color_next,
                            \(a, b) if (is.na(b)) NA_character_
                            else interpolate_colors(a, b, t = 0.5)),
      y_lower    = (psi_lower + lower_end) / 2,
      y_upper    = (psi_upper + upper_end) / 2
    ) |>
    filter(!is.na(x_end), !is.na(color_mid))
  
  # ── Extrapolation zones ───────────────────────────────────────────────────
  # Shade from observed maximum of the "to" stage to the segment end (x_max).
  # obs_max is already capped at x_offset + 1 so xmax never exceeds segment end.
  extrap_df <- meta_df |>
    mutate(xmin = obs_max, xmax = x_max) |>
    filter(xmax > xmin)
  
  # ── Junction labels ───────────────────────────────────────────────────────
  junction_labels <- tibble(
    x     = c(0,           1,        2,         3),
    label = c("Stand init", "Compex", "Thinned", "Mature"),
    color = stage_colors[c("standinit", "compex", "thin", "mature")]
  )
  
  ggplot() +
    
    # Ribbon
    geom_rect(
      data    = ribbon_strips,
      mapping = aes(xmin = x, xmax = x_end,
                    ymin = y_lower, ymax = y_upper,
                    fill = color_mid),
      alpha       = 0.25,
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    
    # Extrapolation shading — beyond observed maximum
    geom_rect(
      data    = extrap_df,
      mapping = aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1),
      fill    = "grey", alpha = 0.25,
      inherit.aes = FALSE
    ) +
    
    # Segment junction lines
    geom_vline(
      xintercept = c(1, 2),
      linewidth  = 0.4,
      color      = "grey80"
    ) +
    
    # Mean line
    geom_segment(
      data    = line_segs,
      mapping = aes(x = x, xend = x_end,
                    y = psi_mean, yend = y_end,
                    color = color_mid),
      linewidth   = 1.1,
      lineend     = "round",
      inherit.aes = FALSE
    ) +
    
    # Rug — one tick per site per segment
    # Position: p_to / (p_from + p_to) — relative balance of the two stages
    # Alpha:    p_from + p_to — how much of the homerange these two stages account for
    #           Dark ticks = site dominated by this pair; faint = diluted by other stages
    geom_rug(
      data        = rug_df,
      mapping     = aes(x = x, alpha = alpha),
      sides       = "b",
      color       = "black",
      length      = unit(0.03, "npc"),
      inherit.aes = FALSE
    ) +
    scale_alpha_continuous(range = c(0.02, 0.4), guide = "none") +
    
    # Junction labels — centered, colored by stage, clip = "off" allows
    # leftmost and rightmost to spill beyond panel edges
    geom_text(
      data        = junction_labels,
      mapping     = aes(x = x, y = 1.0, label = label, color = color),
      hjust       = 0.5,
      vjust       = -0.4,
      size        = 2.8,
      fontface    = "bold",
      inherit.aes = FALSE
    ) +
    
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c("0:100%", "50%", "0:100%", "50%", "0:100%", "50%", "0:100%"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = scales::percent,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x     = "Pairwise homerange composition",
      y     = "Predicted occupancy probability",
      title = sp_name
    ) +
    theme_classic() +
    theme(
      plot.margin        = margin(20, 20, 20, 20),
      legend.position    = "none",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ── MULTI-SPECIES CONCATENATED GRADIENT PREDICTION ────────────────────────────
# Runs predict_gradient_concat() for each species and returns:
#   species_df:  stacked per-species pred_df with sp_name column
#   cumsum_df:   cumulative occupancy summed across species (draws-based)
predict_gradient_concat_multi <- function(species_names, n_grid = 501) {
  
  message("Multi-species gradient: ", length(species_names), " species")
  n_sp      <- length(species_names)
  psi_sum   <- NULL   # will be [n_draws x n_grid_total] accumulated sum
  sp_dfs    <- vector("list", n_sp)
  ref_x     <- NULL   # x and hex_color template from first species
  
  for (k in seq_along(species_names)) {
    sp  <- species_names[k]
    message("  ", sp)
    res <- predict_gradient_concat(sp, n_grid = n_grid)
    
    # Accumulate raw draws for cumulative sum
    # predict_gradient_concat doesn't return psi_draws directly, so we
    # recompute the summary columns from psi_draws via predict_hr_gradient
    # for each segment and stack them. Instead, use the simpler approach:
    # re-run predict_gradient_concat and extract psi_draws per segment.
    sp_dfs[[k]] <- res$pred_df |> mutate(sp_name = sp)
    
    if (is.null(ref_x)) {
      ref_x <- res$pred_df |> select(x, hex_color, segment, from, to, p_to_raw)
    }
  }
  
  # Cumulative sum on posterior mean values
  # For proper uncertainty propagation this would need raw draws, but since
  # predict_gradient_concat does not return them, we sum the posterior means
  # and approximate the CI as sqrt(sum of squared CI half-widths) — valid
  # under approximate independence across species
  species_df <- bind_rows(sp_dfs)
  
  cumsum_df <- species_df |>
    group_by(x) |>
    summarise(
      psi_sum    = sum(psi_mean),
      psi_se     = sqrt(sum(((psi_upper - psi_lower) / (2 * 1.96))^2)),
      psi_lower  = psi_sum - 1.96 * psi_se,
      psi_upper  = psi_sum + 1.96 * psi_se,
      .groups    = "drop"
    ) |>
    left_join(ref_x, by = "x") |>
    rename(psi_mean = psi_sum)
  
  list(
    species_df = species_df,
    cumsum_df  = cumsum_df,
    n_sp       = n_sp
  )
}

# ── MULTI-SPECIES CONCATENATED GRADIENT PLOT ──────────────────────────────────
# One mean line per species. No ribbons or rugs.
# Stage color gradient is blended with a per-species color so lines remain
# distinguishable while the stage gradient stays readable.
#
# pred_df:   output of predict_gradient_concat_multi()
# sp_colors: optional named character vector of hex colors, one per species.
#            If NULL, uses a discrete palette automatically.
plot_gradient_concat_multi <- function(multi_result,
                                       sp_colors = NULL) {
  
  pred_df  <- multi_result$species_df
  sp_names <- unique(pred_df$sp_name)
  n_sp     <- length(sp_names)
  
  if (is.null(sp_colors)) {
    sp_colors <- setNames(scales::hue_pal()(n_sp), sp_names)
  }
  
  # Build line segments — colored by species only, no stage gradient blending
  line_segs <- pred_df |>
    group_by(sp_name) |>
    mutate(
      x_end = lead(x),
      y_end = lead(psi_mean)
    ) |>
    ungroup() |>
    filter(!is.na(x_end))
  
  junction_labels <- tibble(
    x     = c(0,            1,        2,         3),
    label = c("Stand init", "Compex", "Thinned", "Mature")
  )
  
  ggplot() +
    
    geom_vline(xintercept = c(1, 2), linewidth = 0.4, color = "grey80") +
    
    geom_segment(
      data    = line_segs,
      mapping = aes(x = x, xend = x_end,
                    y = psi_mean, yend = y_end,
                    color = sp_name,
                    group = sp_name),
      linewidth   = 1.0,
      lineend     = "round",
      inherit.aes = FALSE
    ) +
    
    scale_color_manual(
      values = sp_colors,
      name   = NULL
    ) +
    
    geom_text(
      data        = junction_labels,
      mapping     = aes(x = x, y = 1.0, label = label),
      color       = "grey30",
      hjust       = 0.5, vjust = -0.4,
      size        = 2.8, fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c("0:100%", "50%", "0:100%", "50%", "0:100%", "50%", "0:100%"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = scales::percent,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Pairwise homerange composition",
      y = "Predicted occupancy probability"
    ) +
    theme_classic() +
    theme(
      plot.margin        = margin(20, 20, 20, 20),
      legend.position    = "right",
      legend.text        = element_text(size = 8),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ── FACETED MULTI-SPECIES PLOT ─────────────────────────────────────────────────
# Vertically stacked panels, one per species, sharing the x axis.
# Species are ordered by which segment's peak occupancy is highest, with
# stage priority standinit < compex < thin < mature — so species that peak
# in SI appear at the top, mature-peaking species at the bottom.
#
# pred_df: output of predict_gradient_concat_multi()
plot_gradient_concat_facet <- function(multi_result) {
  
  pred_df <- multi_result$species_df
  
  # ── Species ordering ────────────────────────────────────────────────────────
  sp_order <- pred_df |>
    group_by(sp_name) |>
    slice_max(psi_mean, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(x) |>
    pull(sp_name)
  
  pred_df <- pred_df |>
    mutate(sp_name = factor(sp_name, levels = sp_order))
  
  # Line segments colored by stage gradient as in single-species plot
  line_segs <- pred_df |>
    group_by(sp_name) |>
    mutate(
      x_end      = lead(x),
      y_end      = lead(psi_mean),
      color_next = lead(hex_color),
      color_mid  = map2_chr(hex_color, color_next,
                            \(a, b) if (is.na(b)) NA_character_
                            else interpolate_colors(a, b, t = 0.5))
    ) |>
    ungroup() |>
    filter(!is.na(x_end), !is.na(color_mid))
  
  # Ribbon strips
  ribbon_strips <- pred_df |>
    group_by(sp_name) |>
    mutate(
      x_end      = lead(x),
      lower_end  = lead(psi_lower),
      upper_end  = lead(psi_upper),
      color_next = lead(hex_color),
      color_mid  = map2_chr(hex_color, color_next,
                            \(a, b) if (is.na(b)) NA_character_
                            else interpolate_colors(a, b, t = 0.5)),
      y_lower    = (psi_lower + lower_end) / 2,
      y_upper    = (psi_upper + upper_end) / 2
    ) |>
    ungroup() |>
    filter(!is.na(x_end), !is.na(color_mid))
  
  # Junction labels — drawn once, replicated across facets via geom_text
  junction_labels <- tibble(
    x     = c(0,            1,        2,         3),
    label = c("Stand init", "Compex", "Thinned", "Mature"),
    color = stage_colors[c("standinit", "compex", "thin", "mature")]
  )
  
  ggplot() +
    
    scale_fill_identity() +
    scale_color_identity() +
    
    # Ribbon
    geom_rect(
      data    = ribbon_strips,
      mapping = aes(xmin = x, xmax = x_end,
                    ymin = y_lower, ymax = y_upper,
                    fill = color_mid,
                    group = sp_name),
      alpha       = 0.25,
      inherit.aes = FALSE
    ) +
    
    # Junction lines
    geom_vline(xintercept = c(1, 2), linewidth = 0.4, color = "grey80") +
    
    # Mean line
    geom_segment(
      data    = line_segs,
      mapping = aes(x = x, xend = x_end,
                    y = psi_mean, yend = y_end,
                    color = color_mid,
                    group = sp_name),
      linewidth   = 0.9,
      lineend     = "round",
      inherit.aes = FALSE
    ) +
    
    # Junction labels — top of each panel
    geom_text(
      data        = junction_labels,
      mapping     = aes(x = x, y = 1.0, label = label),
      color       = junction_labels$color,
      hjust       = 0.5, vjust = -0.3,
      size        = 2.5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    facet_wrap(
      ~ sp_name,
      ncol   = 1,
      scales = "free_y",   # y free so each species uses its own range
      strip.position = "right"
    ) +
    
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c("0:100%", "50%", "0:100%", "50%", "0:100%", "50%", "0:100%"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = scales::percent,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Pairwise homerange composition",
      y = "Predicted occupancy probability"
    ) +
    theme_classic() +
    theme(
      plot.margin        = margin(20, 20, 20, 20),
      strip.background   = element_blank(),
      strip.text.y.right = element_text(angle = 0, hjust = 0, size = 8,
                                        face = "italic"),
      panel.spacing      = unit(0.3, "lines"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      axis.text.x        = element_text(size = 7),
      # Only show x axis labels on the bottom panel
      axis.title.x       = element_text(size = 9)
    )
}

# ── CUMULATIVE OCCUPANCY CONCATENATED GRADIENT PLOT ───────────────────────────
# Plots the sum of posterior mean occupancy across all species in multi_result,
# with an approximate 95% CI. The y-axis represents cumulative expected number
# of species occupying a site, ranging from 0 to n_species.
#
# CI approximation: sqrt(sum of squared species-level CI half-widths),
# valid under approximate independence across species.
plot_gradient_concat_cumsum <- function(multi_result, show_ci = TRUE) {
  
  cumsum_df <- multi_result$cumsum_df
  n_sp      <- multi_result$n_sp
  
  junction_labels <- tibble(
    x     = c(0,            1,        2,         3),
    label = c("Stand init", "Compex", "Thinned", "Mature"),
    color = stage_colors[c("standinit", "compex", "thin", "mature")]
  )
  
  line_segs <- cumsum_df |>
    mutate(
      x_end      = lead(x),
      y_end      = lead(psi_mean),
      color_next = lead(hex_color),
      color_mid  = map2_chr(hex_color, color_next,
                            \(a, b) if (is.na(b)) NA_character_
                            else interpolate_colors(a, b, t = 0.5))
    ) |>
    filter(!is.na(x_end), !is.na(color_mid))
  
  ribbon_strips <- if (show_ci) {
    cumsum_df |>
      mutate(
        x_end      = lead(x),
        lower_end  = lead(psi_lower),
        upper_end  = lead(psi_upper),
        color_next = lead(hex_color),
        color_mid  = map2_chr(hex_color, color_next,
                              \(a, b) if (is.na(b)) NA_character_
                              else interpolate_colors(a, b, t = 0.5)),
        y_lower    = (psi_lower + lower_end) / 2,
        y_upper    = (psi_upper + upper_end) / 2
      ) |>
      filter(!is.na(x_end), !is.na(color_mid))
  } else NULL
  
  p <- ggplot() +
    scale_fill_identity() +
    scale_color_identity()
  
  if (show_ci) {
    p <- p +
      geom_rect(
        data    = ribbon_strips,
        mapping = aes(xmin = x, xmax = x_end,
                      ymin = y_lower, ymax = y_upper,
                      fill = color_mid),
        alpha       = 0.25,
        inherit.aes = FALSE
      )
  }
  
  p +
    geom_vline(xintercept = c(1, 2), linewidth = 0.4, color = "grey80") +
    
    geom_segment(
      data    = line_segs,
      mapping = aes(x = x, xend = x_end,
                    y = psi_mean, yend = y_end,
                    color = color_mid),
      linewidth   = 1.1,
      lineend     = "round",
      inherit.aes = FALSE
    ) +
    
    geom_text(
      data        = junction_labels,
      mapping     = aes(x = x, y = n_sp, label = label, color = color),
      hjust       = 0.5, vjust = -0.4,
      size        = 2.8, fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c("0:100%", "50%", "0:100%", "50%", "0:100%", "50%", "0:100%"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, n_sp),
      breaks = pretty(c(0, n_sp), n = 5),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x        = "Pairwise homerange composition",
      y        = "Cumulative predicted occupancy (no. species)",
      title    = paste0("Cumulative occupancy (", n_sp, " species)"),
      subtitle = "Sum of posterior mean occupancy probabilities across all species"
    ) +
    theme_classic() +
    theme(
      plot.margin        = margin(20, 20, 20, 20),
      legend.position    = "none",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}
# Predicts occupancy for a "typical species" using community hyperparameters
# (mu.u, mu.alpha_plot*, mu.alpha_homerange*) in place of species-specific
# parameters. Group index g = 1 throughout (single "all" group).
#
# Raw proportions (0-1) are passed directly to mu.alpha_homerange* without
# species-specific scaling. The x-axis therefore represents raw proportions,
# and the coefficients operate on that scale. Point covariates and season are
# fixed at 0 (grand mean on scaled units).
#
# Option B marginalization still applies: the survey point is weighted by
# the from/to stage proportions at each grid point.
predict_hr_gradient_community <- function(from   = "compex",
                                          to     = "mature",
                                          n_grid = 501) {
  
  stopifnot(from %in% names(stage_idx_map), to %in% names(stage_idx_map))
  stopifnot(from != to)
  
  g          <- 1   # group index
  p_to_raw   <- seq(0, 1, length.out = n_grid)
  p_from_raw <- 1 - p_to_raw
  
  # Community hyperparameter draws — [n_draws] or [n_draws, G] sliced to g
  mu_u   <- sims$mu.u                         # [n_draws] if G=1, else [n_draws,G]
  if (is.matrix(mu_u)) mu_u <- mu_u[, g]
  
  mu_ahr <- list(
    standinit = if (is.matrix(sims$mu.alpha_homerange1)) sims$mu.alpha_homerange1[, g]
    else sims$mu.alpha_homerange1,
    thin      = if (is.matrix(sims$mu.alpha_homerange2)) sims$mu.alpha_homerange2[, g]
    else sims$mu.alpha_homerange2,
    mature    = if (is.matrix(sims$mu.alpha_homerange3)) sims$mu.alpha_homerange3[, g]
    else sims$mu.alpha_homerange3
  )
  
  # mu.alpha_plot*[s, g] — [n_draws, n_stages, n_groups] or [n_draws, n_stages]
  get_mu_ap <- function(param, stage) {
    arr <- sims[[param]]
    if (length(dim(arr)) == 3) arr[, stage_idx_map[stage], g]
    else                       arr[, stage_idx_map[stage]]
  }
  
  mu_ap_from <- lapply(paste0("mu.alpha_plot", 1:3), get_mu_ap, stage = from)
  mu_ap_to   <- lapply(paste0("mu.alpha_plot", 1:3), get_mu_ap, stage = to)
  
  # ── Homerange contribution [n_draws x n_grid] ──────────────────────────────
  # Raw proportions passed directly — no species-specific scaling
  hr_contrib <- matrix(0, nrow = n_draws, ncol = n_grid)
  
  for (s in c("standinit", "thin", "mature")) {
    p_s <- if (s == to)   p_to_raw   else
      if (s == from) p_from_raw else 0   # scalar for fixed stages
    
    contribution <- if (length(p_s) > 1) outer(mu_ahr[[s]], p_s)
    else                 mu_ahr[[s]] * p_s
    
    hr_contrib <- hr_contrib + contribution
  }
  
  # ── Plot contributions [n_draws] ───────────────────────────────────────────
  plot_contrib <- function(ap_list, stage) {
    ap_list[[1]] * mean_plot[[stage]][1] +
      ap_list[[2]] * mean_plot[[stage]][2] +
      ap_list[[3]] * mean_plot[[stage]][3]
  }
  pc_from <- plot_contrib(mu_ap_from, from)
  pc_to   <- plot_contrib(mu_ap_to,   to)
  
  # ── Full linear predictors [n_draws x n_grid] ──────────────────────────────
  eta_from <- hr_contrib + (mu_u + pc_from)
  eta_to   <- hr_contrib + (mu_u + pc_to)
  
  # ── Option B: weighted average ─────────────────────────────────────────────
  p_from_mat <- matrix(p_from_raw, nrow = n_draws, ncol = n_grid, byrow = TRUE)
  p_to_mat   <- matrix(p_to_raw,   nrow = n_draws, ncol = n_grid, byrow = TRUE)
  
  psi_draws <- p_from_mat * plogis(eta_from) +
    p_to_mat   * plogis(eta_to)
  
  tibble(
    p_to      = p_to_raw,
    p_from    = p_from_raw,
    psi_mean  = apply(psi_draws, 2, mean),
    psi_lower = apply(psi_draws, 2, quantile, 0.025),
    psi_upper = apply(psi_draws, 2, quantile, 0.975),
    from      = from,
    to        = to
  )
}

# ── COMMUNITY MEAN CONCATENATED GRADIENT PREDICTION ───────────────────────────
predict_gradient_concat_community <- function(n_grid = 501) {
  
  segments <- list(
    list(from = "standinit", to = "compex", label = "SI \u2192 Compex"),
    list(from = "compex",    to = "thin",   label = "Compex \u2192 Thin"),
    list(from = "thin",      to = "mature", label = "Thin \u2192 Mature")
  )
  
  all_dfs  <- vector("list", length(segments))
  seg_meta <- vector("list", length(segments))
  x_offset <- 0
  
  for (k in seq_along(segments)) {
    seg <- segments[[k]]
    df  <- predict_hr_gradient_community(from   = seg$from,
                                         to     = seg$to,
                                         n_grid = n_grid)
    
    x_vals     <- seq(x_offset, x_offset + 1, length.out = n_grid)
    hex_colors <- interpolate_colors(stage_colors[seg$from], stage_colors[seg$to],
                                     t = df$p_to)
    
    df <- df |>
      mutate(x = x_vals, segment = k, seg_label = seg$label, hex_color = hex_colors)
    
    if (k < length(segments)) df <- df[-nrow(df), ]
    all_dfs[[k]] <- df
    
    seg_meta[[k]] <- tibble(
      seg_label = seg$label, from = seg$from, to = seg$to,
      segment = k, x_min = x_offset, x_max = x_offset + 1
    )
    
    x_offset <- x_offset + 1
  }
  
  list(
    pred_df = bind_rows(all_dfs),
    meta_df = bind_rows(seg_meta)
  )
}

# ── COMMUNITY MEAN CONCATENATED GRADIENT PLOT ─────────────────────────────────
# Identical visual style to plot_gradient_concat() but for the community mean.
# Includes the stage color gradient on the line and optional BCI ribbon.
plot_gradient_concat_community <- function(result, show_ci = TRUE) {
  
  pred_df <- result$pred_df
  
  junction_labels <- tibble(
    x     = c(0,            1,        2,         3),
    label = c("Stand init", "Compex", "Thinned", "Mature"),
    color = stage_colors[c("standinit", "compex", "thin", "mature")]
  )
  
  line_segs <- pred_df |>
    mutate(
      x_end      = lead(x),
      y_end      = lead(psi_mean),
      color_next = lead(hex_color),
      color_mid  = map2_chr(hex_color, color_next,
                            \(a, b) if (is.na(b)) NA_character_
                            else interpolate_colors(a, b, t = 0.5))
    ) |>
    filter(!is.na(x_end), !is.na(color_mid))
  
  ribbon_strips <- if (show_ci) {
    pred_df |>
      mutate(
        x_end      = lead(x),
        lower_end  = lead(psi_lower),
        upper_end  = lead(psi_upper),
        color_next = lead(hex_color),
        color_mid  = map2_chr(hex_color, color_next,
                              \(a, b) if (is.na(b)) NA_character_
                              else interpolate_colors(a, b, t = 0.5)),
        y_lower    = (psi_lower + lower_end) / 2,
        y_upper    = (psi_upper + upper_end) / 2
      ) |>
      filter(!is.na(x_end), !is.na(color_mid))
  } else NULL
  
  p <- ggplot() +
    scale_fill_identity() +
    scale_color_identity()
  
  if (show_ci) {
    p <- p + geom_rect(
      data    = ribbon_strips,
      mapping = aes(xmin = x, xmax = x_end, ymin = y_lower, ymax = y_upper,
                    fill = color_mid),
      alpha = 0.25, inherit.aes = FALSE
    )
  }
  
  p +
    geom_vline(xintercept = c(1, 2), linewidth = 0.4, color = "grey80") +
    geom_segment(
      data    = line_segs,
      mapping = aes(x = x, xend = x_end, y = psi_mean, yend = y_end,
                    color = color_mid),
      linewidth = 1.1, lineend = "round", inherit.aes = FALSE
    ) +
    geom_text(
      data        = junction_labels,
      mapping     = aes(x = x, y = 1.0, label = label),
      color       = junction_labels$color,
      hjust       = 0.5, vjust = -0.4,
      size        = 2.8, fontface = "bold", inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c("0:100%", "50%", "0:100%", "50%", "0:100%", "50%", "0:100%"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = c(0, 0.5, 1),
      labels = scales::percent, expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x        = "Pairwise homerange composition",
      y        = "Predicted occupancy probability",
      title    = "Community mean (typical species)",
      subtitle = paste(
        "Hyperparameter predictions using mu.u, mu.alpha_homerange*, mu.alpha_plot*",
        "Raw proportions passed directly | Point covariates at grand mean",
        sep = "\n"
      )
    ) +
    theme_classic() +
    theme(
      plot.margin        = margin(20, 20, 20, 20),
      legend.position    = "none",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ── COMMUNITY MEAN CONCATENATED GRADIENT PREDICTION ───────────────────────────
# Predicts occupancy for a hypothetical average species using community
# hyperparameters (mu.u, mu.alpha_homerange*, mu.alpha_plot*).
#
# Raw proportions (0-1) are passed directly to mu.alpha_homerange — no
# species-specific scaling is applied. The x-axis therefore represents raw
# proportions; the coefficients operate on that scale.
#
# Option B marginalization is applied exactly as for species-specific
# predictions: occupancy at each grid point is a weighted average of the
# per-stage predictions, weighted by the homerange proportions.
#
# g: group index (default 1 = "all")
predict_gradient_concat_community <- function(n_grid = 501, g = 1) {
  
  segments <- list(
    list(from = "standinit", to = "compex", label = "SI \u2192 Compex"),
    list(from = "compex",    to = "thin",   label = "Compex \u2192 Thin"),
    list(from = "thin",      to = "mature", label = "Thin \u2192 Mature")
  )
  
  # Community hyperparameter draws [n_draws] — index group g
  # mu.u and mu.alpha_homerange* may be [n_draws] (G=1) or [n_draws, G]
  # Use drop=FALSE and then [,g] to handle both cases safely
  extract_mu <- function(param) {
    x <- sims[[param]]
    if (is.matrix(x)) x[, g] else x
  }
  
  mu_u   <- extract_mu("mu.u")
  mu_hr1 <- extract_mu("mu.alpha_homerange1")  # pcnt_standinit
  mu_hr2 <- extract_mu("mu.alpha_homerange2")  # pcnt_thin
  mu_hr3 <- extract_mu("mu.alpha_homerange3")  # pcnt_mature
  
  # mu.alpha_plot[stage, group] — [n_draws, n_stages, n_groups]
  # or [n_draws, n_stages] if G=1
  extract_mu_plot <- function(param, stage_idx) {
    x <- sims[[param]]
    if (length(dim(x)) == 3) x[, stage_idx, g] else x[, stage_idx]
  }
  
  all_dfs  <- vector("list", length(segments))
  seg_meta <- vector("list", length(segments))
  x_offset <- 0
  
  for (k in seq_along(segments)) {
    seg        <- segments[[k]]
    from       <- seg$from
    to         <- seg$to
    s_from     <- stage_idx_map[from]
    s_to       <- stage_idx_map[to]
    n_grid_seg <- n_grid
    
    p_to_raw   <- seq(0, 1, length.out = n_grid_seg)
    p_from_raw <- 1 - p_to_raw
    
    # Plot covariate contributions [n_draws] using stage-conditional means
    pc <- function(stage) {
      s  <- stage_idx_map[stage]
      extract_mu_plot("mu.alpha_plot1", s) * mean_plot[[stage]][1] +
        extract_mu_plot("mu.alpha_plot2", s) * mean_plot[[stage]][2] +
        extract_mu_plot("mu.alpha_plot3", s) * mean_plot[[stage]][3]
    }
    pc_from <- pc(from)
    pc_to   <- pc(to)
    
    # Homerange contributions — raw proportions, no species-specific scaling
    # For each non-compex stage: varies (outer) or fixed at 0 (scalar)
    hr_contrib <- matrix(0, nrow = n_draws, ncol = n_grid_seg)
    
    for (s_name in c("standinit", "thin", "mature")) {
      mu_hr <- switch(s_name,
                      standinit = mu_hr1,
                      thin      = mu_hr2,
                      mature    = mu_hr3
      )
      p_s <- if (s_name == to)   p_to_raw   else
        if (s_name == from) p_from_raw else
          0
      
      contribution <- if (length(p_s) > 1) outer(mu_hr, p_s) else mu_hr * p_s
      hr_contrib   <- hr_contrib + contribution
    }
    
    # Full linear predictors [n_draws x n_grid]
    eta_from <- hr_contrib + (mu_u + pc_from)
    eta_to   <- hr_contrib + (mu_u + pc_to)
    
    # Option B: weighted average
    p_from_mat <- matrix(p_from_raw, nrow = n_draws, ncol = n_grid_seg, byrow = TRUE)
    p_to_mat   <- matrix(p_to_raw,   nrow = n_draws, ncol = n_grid_seg, byrow = TRUE)
    
    psi_draws <- p_from_mat * plogis(eta_from) +
      p_to_mat   * plogis(eta_to)
    
    x_vals     <- seq(x_offset, x_offset + 1, length.out = n_grid_seg)
    hex_colors <- interpolate_colors(
      color_from = stage_colors[from],
      color_to   = stage_colors[to],
      t          = p_to_raw
    )
    
    df <- tibble(
      x         = x_vals,
      segment   = k,
      seg_label = seg$label,
      from      = from,
      to        = to,
      p_to_raw  = p_to_raw,
      hex_color = hex_colors,
      psi_mean  = apply(psi_draws, 2, mean),
      psi_lower = apply(psi_draws, 2, quantile, 0.025),
      psi_upper = apply(psi_draws, 2, quantile, 0.975)
    )
    
    if (k < length(segments)) df <- df[-nrow(df), ]
    all_dfs[[k]] <- df
    
    seg_meta[[k]] <- list(
      seg_label = seg$label,
      from      = from,
      to        = to,
      segment   = k,
      x_min     = x_offset,
      x_max     = x_offset + 1
    )
    
    x_offset <- x_offset + 1
  }
  
  list(
    pred_df = bind_rows(all_dfs),
    meta_df = map_dfr(seg_meta, as_tibble)
  )
}

# ── COMMUNITY MEAN CONCATENATED GRADIENT PLOT ─────────────────────────────────
# Plots the community mean occupancy probability across the concatenated stage
# gradient. Includes optional 95% CI ribbon.
plot_gradient_concat_community <- function(result, show_ci = TRUE) {
  
  pred_df <- result$pred_df
  
  junction_labels <- tibble(
    x     = c(0,            1,        2,         3),
    label = c("Stand init", "Compex", "Thinned", "Mature"),
    color = stage_colors[c("standinit", "compex", "thin", "mature")]
  )
  
  line_segs <- pred_df |>
    mutate(
      x_end      = lead(x),
      y_end      = lead(psi_mean),
      color_next = lead(hex_color),
      color_mid  = map2_chr(hex_color, color_next,
                            \(a, b) if (is.na(b)) NA_character_
                            else interpolate_colors(a, b, t = 0.5))
    ) |>
    filter(!is.na(x_end), !is.na(color_mid))
  
  ribbon_strips <- if (show_ci) {
    pred_df |>
      mutate(
        x_end      = lead(x),
        lower_end  = lead(psi_lower),
        upper_end  = lead(psi_upper),
        color_next = lead(hex_color),
        color_mid  = map2_chr(hex_color, color_next,
                              \(a, b) if (is.na(b)) NA_character_
                              else interpolate_colors(a, b, t = 0.5)),
        y_lower    = (psi_lower + lower_end) / 2,
        y_upper    = (psi_upper + upper_end) / 2
      ) |>
      filter(!is.na(x_end), !is.na(color_mid))
  } else NULL
  
  p <- ggplot() +
    scale_fill_identity() +
    scale_color_identity()
  
  if (show_ci) {
    p <- p +
      geom_rect(
        data    = ribbon_strips,
        mapping = aes(xmin = x, xmax = x_end,
                      ymin = y_lower, ymax = y_upper,
                      fill = color_mid),
        alpha       = 0.25,
        inherit.aes = FALSE
      )
  }
  
  p +
    geom_vline(xintercept = c(1, 2), linewidth = 0.4, color = "grey80") +
    
    geom_segment(
      data    = line_segs,
      mapping = aes(x = x, xend = x_end,
                    y = psi_mean, yend = y_end,
                    color = color_mid),
      linewidth   = 1.1,
      lineend     = "round",
      inherit.aes = FALSE
    ) +
    
    geom_text(
      data        = junction_labels,
      mapping     = aes(x = x, y = 1.0, label = label, color = color),
      hjust       = 0.5, vjust = -0.4,
      size        = 2.8, fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      labels = c("0:100%", "50%", "0:100%", "50%", "0:100%", "50%", "0:100%"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = scales::percent,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x        = "Pairwise homerange composition",
      y        = "Predicted occupancy probability",
      title    = "Community mean",
      subtitle = paste(
        "Predicted for a hypothetical average species (hyperparameter means)",
        "Raw proportions applied directly to community-mean coefficients",
        sep = "\n"
      )
    ) +
    theme_classic() +
    theme(
      plot.margin        = margin(20, 20, 20, 20),
      legend.position    = "none",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ── USAGE ─────────────────────────────────────────────────────────────────────
stop()

# Single species
result <- predict_gradient_concat("pileated woodpecker", n_grid = 101)
plot_gradient_concat(result)

result <- predict_gradient_concat("vaux's swift", n_grid = 101)
plot_gradient_concat(result)

result <- predict_gradient_concat("golden-crowned kinglet", n_grid = 501)
plot_gradient_concat(result)

# Multiple species — overlaid
focal_species <- c("pileated woodpecker", "hairy woodpecker", "red-breasted nuthatch")
multi_pred <- predict_gradient_concat_multi(focal_species, n_grid = 501)
plot_gradient_concat_multi(multi_pred)
plot_gradient_concat_facet(multi_pred)

# Larger set
focal_species <- intersect(species, species_traits %>% filter(group_nest_ps == "ground") %>% pull(common_name))
multi_pred <- predict_gradient_concat_multi(focal_species, n_grid = 101)

plot_gradient_concat_cumsum(multi_pred)
plot_gradient_concat_cumsum(multi_pred, show_ci = FALSE)

# All species cumulative
all_pred <- predict_gradient_concat_multi(model_data$species, n_grid = 101)
plot_gradient_concat_cumsum(all_pred)

# Community mean (hyperparameter-based)
comm_result <- predict_gradient_concat_community(n_grid = 501)
plot_gradient_concat_community(comm_result)
plot_gradient_concat_community(comm_result, show_ci = FALSE)
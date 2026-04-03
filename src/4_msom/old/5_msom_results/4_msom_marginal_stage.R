####################################################################################
# Marginal plots of varying homerange stage composition
#
# CONFIG:
path_msom = "data/cache/models/V3_msom_pcnt_fp_fp_all.rds"
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

# MANGO

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
    
    # # Extrapolation shading — beyond observed maximum
    # geom_rect(
    #   data    = extrap_df,
    #   mapping = aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1),
    #   fill    = "grey", alpha = 0.25,
    #   inherit.aes = FALSE
    # ) +
    
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
      labels = c("100:0%", "50%", "100:0%", "50%", "100:0%", "50%", "100:0%"),
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
      title = str_to_title(sp_name)
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
result <- predict_gradient_concat("pileated woodpecker", n_grid = 101)
plot_gradient_concat(result)

result <- predict_gradient_concat("vaux's swift", n_grid = 101)
plot_gradient_concat(result)

result <- predict_gradient_concat("rufous hummingbird", n_grid = 101)
plot_gradient_concat(result)

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

# Stage colors — user-specified
stage_colors <- c(
  standinit = "#E8A838",   # orange
  compex    = "#4CAF6B",   # green
  thin      = "#8B5DB5",   # purple
  mature    = "#8B5E3C"    # brown
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

# Interpolate between two hex colors using a vector of weights t in [0, 1]
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
# Segments:
#   1: SI   (100%) -> Compex (100%)   x in [0, 1]
#   2: Compex (100%) -> Thin (100%)   x in [1, 2]
#   3: Thin (100%) -> Mature (100%)   x in [2, 3]
#
# Each row in pred_df has a hex_color computed by interpolating between the
# from-stage and to-stage colors using p_to as the blend weight.
predict_gradient_concat <- function(species_name, n_grid = 501) {
  
  segments <- list(
    list(from = "standinit", to = "compex", label = "SI \u2192 Compex"),
    list(from = "compex",    to = "thin",   label = "Compex \u2192 Thin"),
    list(from = "thin",      to = "mature", label = "Thin \u2192 Mature")
  )
  
  all_dfs   <- vector("list", length(segments))
  seg_meta  <- vector("list", length(segments))
  x_offset  <- 0
  
  for (k in seq_along(segments)) {
    seg <- segments[[k]]
    res <- predict_hr_gradient(species_name, from = seg$from, to = seg$to,
                               n_grid = n_grid)
    
    x_vals <- seq(x_offset, x_offset + 1, length.out = n_grid)
    
    # Interpolate color at each grid point: p_to = 0 -> from color, 1 -> to color
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
  
  list(
    pred_df = bind_rows(all_dfs),
    meta_df = map_dfr(seg_meta, as_tibble),
    sp_name = species_name
  )
}

# ── SINGLE-SPECIES CONCATENATED GRADIENT PLOT ─────────────────────────────────
# The continuous color gradient is achieved by:
#   - Line:   geom_segment() between consecutive grid points, colored at midpoint
#   - Ribbon: geom_rect() strips one per grid step, filled at midpoint color
# Both use scale_color_identity() / scale_fill_identity() so hex values are
# used directly without any scale mapping.
plot_gradient_concat <- function(result) {
  
  pred_df <- result$pred_df
  meta_df <- result$meta_df
  sp_name <- result$sp_name
  
  stage_labels <- c(
    standinit = "Stand initiation",
    compex    = "Competitive exclusion",
    thin      = "Thinned",
    mature    = "Mature"
  )
  
  # ── Segment pairs for geom_segment (line) ────────────────────────────────
  # Each row represents one tiny line segment between grid point k and k+1.
  # color_mid is computed row-wise via map2_chr so each segment gets its own
  # pairwise interpolation rather than a ramp through all colors at once.
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
  
  # ── Ribbon strips for geom_rect ───────────────────────────────────────────
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
  extrap_df <- meta_df |>
    mutate(xmin = obs_upper, xmax = obs_max) |>
    filter(xmax > xmin)
  
  # ── Stage labels at segment boundaries ───────────────────────────────────
  # All centered (hjust = 0.5) over their 100% point; clip = "off" allows
  # the leftmost and rightmost labels to spill beyond the panel edges
  junction_labels <- tibble(
    x     = c(0,          1,       2,            3),
    label = c("Stand init", "Compex", "Thinned", "Mature"),
    color = stage_colors[c("standinit", "compex", "thin", "mature")]
  )
  
  # ── Legend data: one point per stage for manual legend ───────────────────
  legend_df <- tibble(
    stage = names(stage_colors),
    label = stage_labels[names(stage_colors)],
    color = stage_colors
  )
  
  ggplot() +
    
    # Ribbon as semi-transparent colored strips
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
    
    # Extrapolation shading (on top of ribbon)
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
      color      = "grey80",
      linetype   = "solid"
    ) +
    
    # Mean line as colored segments — drawn after ribbon so it is always visible
    geom_segment(
      data    = line_segs,
      mapping = aes(x = x, xend = x_end,
                    y = psi_mean, yend = y_end,
                    color = color_mid),
      linewidth   = 1.1,
      lineend     = "round",
      inherit.aes = FALSE
    ) +
    
    # # Segment transition labels
    # geom_text(
    #   data        = meta_df,
    #   mapping     = aes(x = x_min + 0.5, y = 0.97, label = seg_label),
    #   size        = 3,
    #   color       = "grey25",
    #   fontface    = "italic",
    #   inherit.aes = FALSE
    # ) +
    
    # Boundary labels at top — centered over their 100% point, colored by stage
    geom_text(
      data        = junction_labels,
      mapping     = aes(x = x, y = 1.0, label = label, color = color),
      hjust       = 0.5,
      vjust       = -0.4,
      size        = 2.8,
      fontface    = "bold",
      inherit.aes = FALSE
    ) +
    
    # # Manual color legend using dummy points off-plot
    # geom_point(
    #   data    = legend_df,
    #   mapping = aes(x = -Inf, y = -Inf, color = color),
    #   size    = 3,
    #   inherit.aes = FALSE
    # ) +
    # 
    # # Stage color legend via manual guide
    # guides(
    #   color = guide_legend(
    #     title    = "Stage",
    #     override.aes = list(
    #       color = stage_colors,
    #       size  = 3
    #     )
    #   )
    # ) +
    # Inject named legend labels by adding a dummy invisible scale
    # (scale_color_identity doesn't support labels — use annotate approach instead)
    
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
      title    = sp_name
    ) +
    # Stage color swatches as annotation in top-right corner
    # annotate("point",
    #          x     = c(2.55, 2.70, 2.85, 3.00),
    #          y     = rep(0.06, 4),
    #          color = stage_colors[c("standinit", "compex", "thin", "mature")],
    #          size  = 3
    # ) +
    # annotate("text",
    #          x      = c(2.55, 2.70, 2.85, 3.00),
    #          y      = rep(0.02, 4),
    #          label  = c("SI", "CE", "Thin", "Mat"),
    #          size   = 2.5,
    #          color  = "grey30",
    #          hjust  = 0.5
    # ) +
    theme_classic() +
    theme(
      plot.margin        = margin(20, 20, 20, 20),
      legend.position    = "none",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ── USAGE ─────────────────────────────────────────────────────────────────────
# n_grid controls ribbon smoothness — higher = smoother but slower
# 101 = fast/coarse, 501 = good balance, 1001 = very smooth
result <- predict_gradient_concat("pileated woodpecker", n_grid = 501)
plot_gradient_concat(result)

result <- predict_gradient_concat("vaux's swift", n_grid = 101)
plot_gradient_concat(result)

result <- predict_gradient_concat("golden-crowned kinglet", n_grid = 501)
plot_gradient_concat(result)
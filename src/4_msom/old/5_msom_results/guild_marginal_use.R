####################################################################################
#
# CONFIG:
path_msom = "data/cache/models/V4_msom_V4_fp_all.rds"
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

msom_summary %>% filter(str_starts(param, "mu.alpha_homerange"))

####################################################################################
#
# CONFIG: path_msom = "data/cache/models/V4_msom_V4_fp_all.rds"
#
# INPUT: path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
####################################################################################

library(tidyverse)

# ── Load data ──────────────────────────────────────────────────────────────────
message("Loading species trait data from ", path_trait_data)
species_traits <- read_csv(path_trait_data, show_col_types = FALSE)

message("Loading data for multi-species occupancy model ", path_msom)
model_data <- readRDS(path_msom)

msom_summary <- model_data$msom_summary
msom         <- model_data$msom
groups       <- model_data$groups %>% arrange(common_name)
sites        <- model_data$sites
species      <- model_data$species
stages       <- model_data$stages

# ── Extract posterior samples from jagsUI object ──────────────────────────────
post <- as.data.frame(do.call(rbind, msom$samples))

# ── Check homerange hyperparameters ───────────────────────────────────────────
msom_summary %>% filter(str_starts(param, "mu.alpha_homerange"))

####################################################################################
# FUNCTION: plot_homerange_marginal
#
# Produces a marginal occupancy probability plot for a given homerange parameter.
#
# Args:
#   param_num       - Integer 1, 2, or 3 (selects alpha_homerange1/2/3)
#   model_data      - The loaded RDS list
#   post            - Posterior samples data.frame from do.call(rbind, msom$samples)
#   species_filter  - Optional character vector of common_name values to show as
#                     spaghetti lines. NULL (default) shows all species in grey.
#                     When a filter is supplied, lines are coloured by species
#                     with a legend. The community mean ribbon always uses all
#                     species.
#
# Example:
#   plot_homerange_marginal(3, model_data, post)
#   plot_homerange_marginal(3, model_data, post,
#     species_filter = species_traits %>%
#       filter(group_nest_ps == "cavity_p") %>%
#       pull(common_name))
####################################################################################

plot_homerange_marginal <- function(param_num, model_data, post,
                                    species_filter = NULL) {
  
  # ── 1. Pull metadata ──────────────────────────────────────────────────────────
  hr_table   <- model_data$param_alpha_data$param_alpha_homerange_data
  param_row  <- hr_table %>% filter(param == paste0("alpha_homerange", param_num))
  
  stopifnot("param_num must be 1, 2, or 3" = nrow(param_row) == 1)
  
  param_name <- param_row$name
  data_list  <- param_row$data[[1]]
  sp_names   <- names(data_list)
  n_species  <- length(sp_names)
  
  # ── 2. Validate species_filter ────────────────────────────────────────────────
  if (!is.null(species_filter)) {
    unrecognised <- setdiff(species_filter, sp_names)
    if (length(unrecognised) > 0)
      warning("These species were not found and will be ignored:\n  ",
              paste(unrecognised, collapse = ", "))
    sp_show <- intersect(species_filter, sp_names)
    if (length(sp_show) == 0)
      stop("No valid species remaining after filtering.")
  } else {
    sp_show <- sp_names
  }
  
  use_colors <- !is.null(species_filter)
  
  # ── 3. Per-species scale attributes ──────────────────────────────────────────
  sp_attrs <- map_dfr(sp_names, function(sp) {
    v <- data_list[[sp]]
    tibble(
      species = sp,
      center  = attr(v, "scaled:center"),
      scale   = attr(v, "scaled:scale")
    )
  })
  
  # ── 4. True raw range ─────────────────────────────────────────────────────────
  all_raw_vals <- map2(data_list, seq_len(n_species), function(v, i) {
    as.vector(v) * sp_attrs$scale[i] + sp_attrs$center[i]
  }) %>% unlist()
  
  x_seq_raw <- seq(
    max(0, min(all_raw_vals, na.rm = TRUE)),
    min(1, max(all_raw_vals, na.rm = TRUE)),
    length.out = 200
  )
  
  # ── 5. Community-mean marginal predictions ────────────────────────────────────
  mean_center       <- mean(sp_attrs$center)
  mean_scale        <- mean(sp_attrs$scale)
  x_seq_scaled_mean <- (x_seq_raw - mean_center) / mean_scale
  
  mu_u_samp     <- post[["mu.u"]]
  mu_alpha_samp <- post[[paste0("mu.alpha_homerange", param_num)]]
  
  pred_df <- map_dfr(seq_along(x_seq_raw), function(k) {
    psi_samp <- plogis(mu_u_samp + mu_alpha_samp * x_seq_scaled_mean[k])
    tibble(
      x_pct  = x_seq_raw[k] * 100,
      mean   = mean(psi_samp),
      median = median(psi_samp),
      lo95   = quantile(psi_samp, 0.025),
      hi95   = quantile(psi_samp, 0.975),
      lo50   = quantile(psi_samp, 0.250),
      hi50   = quantile(psi_samp, 0.750)
    )
  })
  
  # ── 6. Per-species spaghetti (filtered) ───────────────────────────────────────
  sp_show_idx <- which(sp_names %in% sp_show)
  
  species_pred_df <- map_dfr(sp_show_idx, function(i) {
    alpha_samp     <- post[[paste0("alpha_homerange", param_num, "[", i, "]")]]
    u_samp         <- post[[paste0("u[", i, "]")]]
    x_seq_scaled_i <- (x_seq_raw - sp_attrs$center[i]) / sp_attrs$scale[i]
    
    map_dfr(seq_along(x_seq_raw), function(k) {
      tibble(
        species = sp_names[i],
        x_pct   = x_seq_raw[k] * 100,
        mean    = mean(plogis(u_samp + alpha_samp * x_seq_scaled_i[k]))
      )
    })
  })
  
  # ── 7. Subtitle note ──────────────────────────────────────────────────────────
  sp_label <- if (!use_colors) {
    paste0("all ", n_species, " species")
  } else {
    paste0(length(sp_show), " selected species")
  }
  
  # ── 8. Plot ───────────────────────────────────────────────────────────────────
  p <- ggplot(pred_df, aes(x = x_pct)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#E67E22", alpha = 0.25) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#E67E22", alpha = 0.45) +
    scale_x_continuous(labels = scales::label_number(suffix = "%"),
                       limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(1)) +
    labs(
      x        = paste0(param_name, " in home range (%)"),
      y        = "Occupancy probability",
      title    = paste0("Marginal effect of ", param_name, " on occupancy"),
      subtitle = paste0("Community mean ± 50/95% CI  |  ", sp_label, " shown"),
      color    = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor  = element_blank(),
      legend.position   = if (use_colors) "right" else "none",
      legend.text       = element_text(size = 8),
      legend.key.height = unit(0.6, "lines")
    )
  
  if (use_colors) {
    # Coloured lines + legend, sorted alphabetically for tidy legend
    n_show <- length(sp_show)
    p <- p +
      geom_line(
        data = species_pred_df %>%
          mutate(species = factor(species, levels = sort(unique(species)))),
        aes(y = mean, group = species, color = species),
        alpha = 0.8, linewidth = 0.6
      ) +
      scale_color_manual(
        values = setNames(
          colorRampPalette(RColorBrewer::brewer.pal(min(n_show, 8), "Dark2"))(n_show),
          sort(sp_show)
        )
      )
  } else {
    # Plain grey spaghetti, no legend
    p <- p +
      geom_line(
        data  = species_pred_df,
        aes(y = mean, group = species),
        color = "steelblue", alpha = 0.4, linewidth = 0.5
      )
  }
  
  # Community mean on top
  p <- p +
    geom_line(aes(y = median), color = "#E67E22", linewidth = 1.1)
  
  p
}

# ── Example calls ─────────────────────────────────────────────────────────────
stop()

# All species (grey spaghetti, no legend)
plot_homerange_marginal(3, model_data, post)
plot_homerange_marginal(1, model_data, post)
plot_homerange_marginal(2, model_data, post)

# Cavity nesters only (coloured by species, legend shown)
cavity_spp <- species_traits %>%
  filter(group_nest_ps == "cavity_p") %>%
  pull(common_name)

plot_homerange_marginal(3, model_data, post, species_filter = cavity_spp)
plot_homerange_marginal(1, model_data, post, species_filter = cavity_spp)
plot_homerange_marginal(2, model_data, post, species_filter = cavity_spp)

# Ad-hoc selection
plot_homerange_marginal(3, model_data, post,
                        species_filter = c("pileated woodpecker", "hairy woodpecker", "downy woodpecker"))
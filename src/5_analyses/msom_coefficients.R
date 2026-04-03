####################################################################################
# MSOM coefficient estimates
#
# CONFIG:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
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

# ── EXTRACT HOMERANGE COEFFICIENTS ───────────────────────────────────────────
hr_param_labels <- c(
  alpha_homerange1 = "Stand initiation",
  alpha_homerange2 = "Thinned",
  alpha_homerange3 = "Mature"
)

coef_df <- model_data$msom_summary |>
  filter(str_detect(param, "^alpha_homerange[123]\\[\\d+\\]$")) |>
  mutate(
    param_base = str_extract(param, "alpha_homerange\\d+"),
    sp_idx     = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    stage      = hr_param_labels[param_base],
    sp_name    = model_data$species[sp_idx]
  ) |>
  # Join nest strategy from species_traits
  left_join(
    species_traits |> select(common_name, group_nest_ps),
    by = c("sp_name" = "common_name")
  ) |>
  # overlap0: "0" = CI excludes zero (certain), "1" = CI overlaps zero (uncertain)
  mutate(pt_alpha = if_else(overlap0 == "1", 0.35, 1.0)) |>
  select(sp_name, stage, group_nest_ps, pt_alpha,
         mean, `2.5%`, `25%`, `75%`, `97.5%`, overlap0, f)

# ── ORDER SPECIES ─────────────────────────────────────────────────────────────
coef_df$sp_name = str_to_sentence(coef_df$sp_name)
sp_order <- coef_df |>
  filter(stage == "Stand initiation") |>
  arrange(mean) |>
  pull(sp_name)

coef_df <- coef_df |>
  mutate(
    sp_name       = factor(sp_name, levels = sp_order),
    stage         = factor(stage, levels = c("Stand initiation", "Thinned", "Mature")),
    group_nest_ps = factor(group_nest_ps)
  )

# ── PLOT ──────────────────────────────────────────────────────────────────────
ggplot(coef_df, aes(y = sp_name, color = group_nest_ps)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  # 95% CI — thin, alpha tied to overlap0
  geom_errorbarh(
    aes(xmin = `2.5%`, xmax = `97.5%`, alpha = pt_alpha),
    height    = 0,
    linewidth = 0.4
  ) +
  # 50% CI — thick, alpha tied to overlap0
  # geom_errorbarh(
  #   aes(xmin = `25%`, xmax = `75%`, alpha = pt_alpha),
  #   height    = 0,
  #   linewidth = 1.2
  # ) +
  # Posterior mean point, alpha tied to overlap0
  geom_point(
    aes(x = mean, alpha = pt_alpha),
    shape = 16,
    size  = 2
  ) +
  scale_alpha_identity() +
  # scale_color_primer() +
  # scale_color_npg() +
  scale_color_brewer(palette = "Dark2", name = "Nest strategy") +
  facet_wrap(~ stage, nrow = 1, scales = "free_x") +
  labs(
    x        = "Standardized coefficient estimate",
    y        = NULL,
    title    = "Homerange stage composition effects on occupancy",
  ) +
  theme_classic() +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 8),
    legend.position    = "bottom",
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

# ── EXTRACT POINT COEFFICIENTS ───────────────────────────────────────────────
point_param_labels <- c(
  alpha_point1 = "Point covariate 1",  # replace with your actual covariate names
  alpha_point2 = "Point covariate 2",
  alpha_point3 = "Point covariate 3"
)

coef_df <- model_data$msom_summary |>
  filter(str_detect(param, "^alpha_point[123]\\[\\d+\\]$")) |>
  mutate(
    param_base = str_extract(param, "alpha_point\\d+"),
    sp_idx     = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    stage      = point_param_labels[param_base],
    sp_name    = model_data$species[sp_idx]
  ) |>
  left_join(
    species_traits |> select(common_name, group_nest_ps),
    by = c("sp_name" = "common_name")
  ) |>
  mutate(pt_alpha = if_else(overlap0 == "1", 0.35, 1.0)) |>
  select(sp_name, stage, group_nest_ps, pt_alpha,
         mean, `2.5%`, `25%`, `75%`, `97.5%`, overlap0, f)

# ── EXTRACT POINT AND SEASON COEFFICIENTS ────────────────────────────────────
point_param_labels <- c(
  alpha_point1  = "Elevation",
  alpha_point2  = "Distance to paved road",
  alpha_point3  = "Distance to major watercourse",
  alpha_season  = "Season"
)

coef_df <- model_data$msom_summary |>
  filter(str_detect(param, "^alpha_(point[123]|season)\\[\\d+\\]$")) |>
  mutate(
    param_base = str_extract(param, "alpha_(point\\d+|season)"),
    sp_idx     = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    covariate  = point_param_labels[param_base],
    sp_name    = model_data$species[sp_idx]
  ) |>
  left_join(
    species_traits |> select(common_name, group_nest_ps),
    by = c("sp_name" = "common_name")
  ) |>
  mutate(pt_alpha = if_else(overlap0 == "1", 0.35, 1.0)) |>
  select(sp_name, covariate, group_nest_ps, pt_alpha,
         mean, `2.5%`, `25%`, `75%`, `97.5%`, overlap0, f)

# ── ORDER SPECIES ─────────────────────────────────────────────────────────────
coef_df$sp_name <- str_to_sentence(coef_df$sp_name)

sp_order <- coef_df |>
  filter(covariate == "Season") |>
  arrange(mean) |>
  pull(sp_name)

coef_df <- coef_df |>
  mutate(
    sp_name       = factor(sp_name, levels = sp_order),
    covariate     = factor(covariate, levels = c("Season",
                                                 "Elevation",
                                                 "Distance to paved road",
                                                 "Distance to major watercourse")),
    group_nest_ps = factor(group_nest_ps)
  )

# ── PLOT ──────────────────────────────────────────────────────────────────────
ggplot(coef_df, aes(y = sp_name, color = group_nest_ps)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_errorbarh(
    aes(xmin = `2.5%`, xmax = `97.5%`, alpha = pt_alpha),
    height    = 0,
    linewidth = 0.4
  ) +
  geom_point(
    aes(x = mean, alpha = pt_alpha),
    shape = 16,
    size  = 2
  ) +
  scale_alpha_identity() +
  scale_color_brewer(palette = "Dark2", name = "Nest strategy") +
  facet_wrap(~ covariate, nrow = 1, scales = "free_x") +
  labs(
    x     = "Standardized coefficient estimate",
    y     = NULL,
    title = "Point-scale and season covariate effects on occupancy"
  ) +
  theme_classic() +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 8),
    legend.position    = "bottom",
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

# ── PARAMETER LABELS ─────────────────────────────────────────────────────────
param_labels <- c(
  mu.alpha_homerange1 = "Stand initiation",
  mu.alpha_homerange2 = "Thinned",
  mu.alpha_homerange3 = "Mature",
  mu.alpha_point1     = "Elevation",
  mu.alpha_point2     = "Distance to paved road",
  mu.alpha_point3     = "Distance to major watercourse",
  mu.alpha_season     = "Season"
)

param_groups <- c(
  mu.alpha_homerange1 = "Homerange composition",
  mu.alpha_homerange2 = "Homerange composition",
  mu.alpha_homerange3 = "Homerange composition",
  mu.alpha_point1     = "Point",
  mu.alpha_point2     = "Point",
  mu.alpha_point3     = "Point",
  mu.alpha_season     = "Season"
)

# ── EXTRACT COEFFICIENTS ──────────────────────────────────────────────────────
coef_df <- model_data$msom_summary |>
  filter(str_detect(param, "^mu\\.alpha_(homerange[123]|point[123]|season)(\\[\\d+\\])?$")) |>
  mutate(
    param_base = str_extract(param, "mu\\.alpha_(homerange\\d+|point\\d+|season)"),
    group_idx  = if_else(
      str_detect(param, "\\[\\d+\\]"),
      as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
      1L
    ),
    covariate     = param_labels[param_base],
    covariate_grp = param_groups[param_base],
    pt_alpha      = if_else(overlap0 == "1", 0.35, 1.0)
  ) |>
  left_join(
    model_data$groups |> distinct(group_idx, group),
    by = "group_idx"
  ) |>
  mutate(
    covariate = factor(covariate, levels = rev(c(
      # Homerange — early to late seral
      "Stand initiation", "Thinned", "Mature",
      # Point
      "Elevation", "Distance to paved road", "Distance to major watercourse",
      # Season
      "Season"
    ))),
    covariate_grp = factor(covariate_grp, levels = c(
      "Homerange composition", "Point", "Season"
    ))
  )

# ── PLOT ──────────────────────────────────────────────────────────────────────
ggplot(coef_df, aes(x = mean, y = covariate, color = group)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_errorbarh(
    aes(xmin = `2.5%`, xmax = `97.5%`, alpha = pt_alpha),
    height    = 0.15,
    linewidth = 0.4,
    position  = position_dodge(width = 0.5)
  ) +
  geom_point(
    aes(alpha = pt_alpha),
    shape    = 16,
    size     = 3,
    position = position_dodge(width = 0.5)
  ) +
  scale_alpha_identity() +
  scale_color_brewer(palette = "Dark2", name = "Group") +
  facet_grid(covariate_grp ~ ., scales = "free_y", space = "free_y") +
  labs(
    x     = "Standardized coefficient estimate",
    y     = NULL,
    title = "Community mean covariate effects on occupancy"
  ) +
  theme_classic() +
  theme(
    strip.background   = element_blank(),
    strip.text.y       = element_text(face = "bold", size = 10, angle = 0, hjust = 0),
    axis.text.y        = element_text(size = 9),
    legend.position    = "bottom",
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.spacing      = unit(0.6, "lines")
  )

# ── PARAMETER LABELS ─────────────────────────────────────────────────────────
plot_param_labels <- c(
  mu.alpha_plot1 = "Trees per acre",
  mu.alpha_plot2 = "Quadratic mean diameter",
  mu.alpha_plot3 = "Hardwood basal area"
)

stage_labels <- c(
  "1" = "Complex",
  "2" = "Stand initiation",
  "3" = "Mature",
  "4" = "Thinned"
)

# ── EXTRACT COEFFICIENTS ──────────────────────────────────────────────────────
coef_df <- model_data$msom_summary |>
  filter(str_detect(param, "^mu\\.alpha_plot[123]\\[\\d+,?\\d*\\]$")) |>
  mutate(
    param_base = str_extract(param, "mu\\.alpha_plot\\d+"),
    stage_idx  = str_extract(param, "(?<=\\[)\\d+(?=,|\\])"),
    group_idx  = if_else(
      str_detect(param, ","),
      as.integer(str_extract(param, "(?<=,)\\d+(?=\\])")),
      1L
    ),
    covariate  = plot_param_labels[param_base],
    stage      = stage_labels[stage_idx],
    pt_alpha   = if_else(overlap0 == "1", 0.35, 1.0)
  ) |>
  left_join(
    model_data$groups |> distinct(group_idx, group),
    by = "group_idx"
  ) |>
  mutate(
    covariate = factor(covariate, levels = rev(c(
      "Trees per acre",
      "Quadratic mean diameter",
      "Hardwood basal area"
    ))),
    stage = factor(stage, levels = c(
      "Complex",
      "Stand initiation",
      "Mature",
      "Thinned"
    ))
  )

# ── PLOT ──────────────────────────────────────────────────────────────────────
ggplot(coef_df, aes(x = mean, y = covariate, color = stage, shape = stage)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_errorbarh(
    aes(xmin = `2.5%`, xmax = `97.5%`, alpha = pt_alpha),
    height    = 0.15,
    linewidth = 0.4,
    position  = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(alpha = pt_alpha),
    size     = 3,
    position = position_dodge(width = 0.6)
  ) +
  scale_alpha_identity() +
  scale_shape_manual(
    name   = "Stage",
    values = c(
      "Complex"          = 15,  # filled circle 16
      "Stand initiation" = 16,  # filled triangle 17
      "Mature"           = 17,  # filled square 15
      "Thinned"          = 18   # filled diamond 18
    )
  ) +
  scale_color_manual(
    name   = "Stage",
    values = c(
      "Complex"          = "forestgreen",  # filled circle 16
      "Stand initiation" = "orange",  # filled triangle 17
      "Mature"           = "tan4",  # filled square 15
      "Thinned"          = "purple"   # filled diamond 18
    )
  ) +
  labs(
    x     = "Standardized coefficient estimate",
    y     = NULL,
    title = "Community mean plot-scale covariate effects on occupancy by stage"
  ) +
  theme_classic() +
  theme(
    axis.text.y        = element_text(size = 10),
    legend.position    = "bottom",
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

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
  filter(stage == "Mature") |>
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

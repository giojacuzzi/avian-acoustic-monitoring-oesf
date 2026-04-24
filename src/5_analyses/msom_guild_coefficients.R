source("src/global.R")

path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% mutate(species = common_name)

models = data.frame(
  model = c(
    "All species",
    "Nesting guild",
    "Foraging guild"
  ),
  path = c(
    "data/cache/models/prefinal_msom_jags_nofp_all.rds",
    "data/cache/models/prefinal_msom_jags_nofp_nest_ps.rds",
    "data/cache/models/prefinal_msom_jags_nofp_forage_substrate.rds"
  )
)

results = setNames(
  lapply(models$model, function(x) list()),
  models$model
)
for (m in seq_along(models$model)) {
  
  model = models[m, "model"]
  path = models[m, "path"]
  
  message("Retrieving data for model '", model, "' from ", path)
  
  model_data = read_rds(path)
  msom_summary = summary(model_data$msom)
  stages = model_data$stages
  groups = model_data$groups
  species = model_data$species
  param_alpha_stage = model_data$param_alpha_data$param_alpha_stage
  param_alpha_season = model_data$param_alpha_data$param_alpha_season
  param_alpha_point_data = model_data$param_alpha_data$param_alpha_point_data
  param_alpha_plot_data = model_data$param_alpha_data$param_alpha_plot_data
  param_alpha_homerange_data = model_data$param_alpha_data$param_alpha_homerange_data
  param_beta_point_data = model_data$param_beta_data$param_beta_point_data
  rm(model_data)
  
  match_s = stages %>% distinct() %>% arrange(stage_idx)
  match_g = groups %>% select(group, group_idx) %>% distinct() %>% arrange(group_idx)
  match_i = tibble(species = species, species_idx = 1:length(species))
  
  mu_alpha_homerange = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_homerange")) %>%
    mutate(param = str_extract(parameter, "alpha_homerange\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "mu\\.alpha_homerange(\\d+)\\[(\\d+)\\]")[,2]),
      g = as.integer(str_match(parameter, "mu\\.alpha_homerange(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_alpha_homerange_data %>% select(param, name), by = c("param")) %>%
    left_join(match_g, by = c("g" = "group_idx"))
  results[[model]][["mu_alpha_homerange"]] = mu_alpha_homerange
  
  alpha_homerange = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "alpha_homerange")) %>%
    mutate(param = str_extract(parameter, "alpha_homerange\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "alpha_homerange(\\d+)\\[(\\d+)\\]")[,2]),
      i = as.integer(str_match(parameter, "alpha_homerange(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_alpha_homerange_data %>% select(param, name), by = c("param")) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  results[[model]][["alpha_homerange"]] = alpha_homerange  
  
  # ggplot(alpha_homerange %>% filter(name == "pcnt_standinit") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
  #   geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  #   geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  #   labs(title = "pcnt_standinit") + theme(legend.position = "bottom")
  # ggplot(alpha_homerange %>% filter(name == "pcnt_mature") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
  #   geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  #   geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  #   labs(title = "pcnt_mature") + theme(legend.position = "bottom")
  # ggplot(alpha_homerange %>% filter(name == "pcnt_thin") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
  #   geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  #   geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  #   labs(title = "pcnt_thin") + theme(legend.position = "bottom")
  
  alpha_plot = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "alpha_plot")) %>%
    mutate(param = str_extract(parameter, "alpha_plot\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,2]),
      s = as.integer(str_match(parameter, "alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,3]),
      i = as.integer(str_match(parameter, "alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,4])
    ) %>% left_join(param_alpha_plot_data %>% select(param, name), by = c("param")) %>%
    left_join(match_s, by = c("s" = "stage_idx")) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  results[[model]][["alpha_plot"]] = alpha_plot  
  
  # ggplot(alpha_plot %>% filter(name == "tree_acre_6_mean", stratum_4 == "mature") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group, shape = stratum_4)) +
  #   geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  #   geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  #   labs(title = "alpha_plot x stratum_4") + theme(legend.position = "bottom")
  
  mu_alpha_point = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_point")) %>%
    mutate(param = str_extract(parameter, "alpha_point\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "mu\\.alpha_point(\\d+)\\[(\\d+)\\]")[,2]),
      g = as.integer(str_match(parameter, "mu\\.alpha_point(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_alpha_point_data %>% select(param, name), by = c("param")) %>%
    left_join(match_g, by = c("g" = "group_idx"))
  results[[model]][["mu_alpha_point"]] = mu_alpha_point  
  
  mu_alpha_season = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_season")) %>%
    mutate(param = str_extract(parameter, "alpha_season"), name = "season") %>%
    mutate(
      p = NA,
      g = as.integer(str_match(parameter, "mu\\.alpha_season\\[(\\d+)\\]")[,2])
    ) %>% left_join(match_g, by = c("g" = "group_idx"))
  results[[model]][["mu_alpha_season"]] = mu_alpha_season  
  
  alpha_season = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "alpha_season")) %>%
    mutate(
      i = as.integer(str_match(parameter, "alpha_season\\[(\\d+)\\]")[,2])
    ) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  results[[model]][["alpha_season"]] = alpha_season
  
  alpha_point = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "alpha_point")) %>%
    mutate(param = str_extract(parameter, "alpha_point\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "alpha_point(\\d+)\\[(\\d+)\\]")[,2]),
      i = as.integer(str_match(parameter, "alpha_point(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_alpha_point_data %>% select(param, name), by = c("param")) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))
  results[[model]][["alpha_point"]] = alpha_point
  
  # ggplot(alpha_season %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group)) +
  #   geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  #   geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  #   labs(title = "season") + theme(legend.position = "bottom")
  
  mu_alpha_plot = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.alpha_plot")) %>%
    mutate(param = str_extract(parameter, "alpha_plot\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "mu\\.alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,2]),
      s = as.integer(str_match(parameter, "mu\\.alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,3]),
      g = as.integer(str_match(parameter, "mu\\.alpha_plot(\\d+)\\[(\\d+),(\\d+)\\]")[,4])
    ) %>% left_join(param_alpha_plot_data %>% select(param, name), by = c("param")) %>%
    left_join(match_s, by = c("s" = "stage_idx")) %>% left_join(match_g, by = c("g" = "group_idx"))
  results[[model]][["mu_alpha_plot"]] = mu_alpha_plot
  
  mu_beta_point = msom_summary %>% as.data.frame() %>% rownames_to_column("parameter") %>%
    filter(str_starts(parameter, "mu.beta_point")) %>%
    mutate(param = str_extract(parameter, "beta_point\\d+")) %>%
    mutate(
      p = as.integer(str_match(parameter, "mu\\.beta_point(\\d+)\\[(\\d+)\\]")[,2]),
      g = as.integer(str_match(parameter, "mu\\.beta_point(\\d+)\\[(\\d+)\\]")[,3])
    ) %>% left_join(param_beta_point_data %>% select(param, name), by = c("param")) %>% left_join(match_g, by = c("g" = "group_idx"))
  results[[model]][["mu_beta_point"]] = mu_beta_point
}

# Guild coefficients ─────────────────────────────────────────────────────────────

mu_alpha_site_all = rbind(
  results[["All species"]][["mu_alpha_homerange"]],
  results[["All species"]][["mu_alpha_point"]],
  results[["All species"]][["mu_alpha_season"]]
)
mu_alpha_site_all$group = "All species"
mu_alpha_site_all$model = "All species"

mu_alpha_site_nesting = rbind(
  results[["Nesting guild"]][["mu_alpha_homerange"]],
  results[["Nesting guild"]][["mu_alpha_point"]],
  results[["Nesting guild"]][["mu_alpha_season"]]
)
mu_alpha_site_nesting$group = str_to_sentence(mu_alpha_site_nesting$group)
mu_alpha_site_nesting$model = "Nesting guild"

mu_alpha_site_foraging = rbind(
  results[["Foraging guild"]][["mu_alpha_homerange"]],
  results[["Foraging guild"]][["mu_alpha_point"]],
  results[["Foraging guild"]][["mu_alpha_season"]]
)
mu_alpha_site_foraging$group = str_to_sentence(mu_alpha_site_foraging$group)
mu_alpha_site_foraging$model = "Foraging guild"

mu_alpha_site = bind_rows(mu_alpha_site_all, mu_alpha_site_nesting, mu_alpha_site_foraging)
mu_alpha_site$model = factor(mu_alpha_site$model, levels = c("Foraging guild", "Nesting guild", "All species"))
mu_alpha_site$name = factor(mu_alpha_site$name,
                             levels = c("pcnt_standinit", "pcnt_thin", "pcnt_mature", 
                                        "season", "elevation", "dist_road_paved", "dist_watercourse_major"),
                             labels = c("Stand Initiation", "Thinned", "Mature", 
                                        "Season", "Elevation", "Distance Road", "Distance Watercourse")
)
mu_alpha_site$group = factor(mu_alpha_site$group,
                             levels = c("All species",
                                        "Cavity_p", "Cavity_s", "Tree", "Shrub", "Ground",
                                        "Aerial forager", "Bark forager", "Foliage gleaner", "Ground forager"),
                             labels = c("All species",
                                        "Cavity, primary", "Cavity, secondary", "Tree", "Shrub", "Ground",
                                        "Aerial forager", "Bark forager", "Foliage gleaner", "Ground forager"))
mu_alpha_site = mu_alpha_site %>%
  mutate(pt_alpha = ifelse(`25%` < 0 & `75%` > 0, 0.3, 1))

fig_guilds = ggplot(mu_alpha_site %>% filter(name %in% c("Stand Initiation", "Thinned", "Mature")) %>% filter(!is.na(group)),
       aes(x = mean, y = (model), color = group, alpha = pt_alpha)) +
  facet_wrap(~name) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_point(position = position_dodge(width = 0.5, reverse = TRUE), size = 1.25) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, linewidth = 0.4, position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbarh(aes(xmin = `25%`, xmax = `75%`), width = 0, linewidth = 0.8, position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_alpha_identity(guide = "none") +
  labs(x = "Standardized coefficient estimate", y = "", color = "Guild") +
  theme(
    panel.border = element_blank(),
    axis.line.x  = element_line(color = "grey70"),
    axis.line.y  = element_line(color = "grey70"),
    legend.position = "left",
    panel.grid.major.x = element_line(color = "gray96")
  ); print(fig_guilds)

ggplot(mu_alpha_site %>% filter(name %in% c("Season", "Elevation", "Distance Road", "Distance Watercourse")),
       aes(x = mean, y = (model), color = group)) +
  facet_wrap(~name) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_point(position = position_dodge(width = 0.5, reverse = TRUE), size = 1.25) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0, linewidth = 0.4, position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbarh(aes(xmin = `25%`, xmax = `75%`), height = 0, linewidth = 0.8, position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  labs(x = "Standardized coefficient estimate", y = "", color = "Guild") +
  theme(
    legend.position = "left",
    panel.grid.major.x = element_line(color = "gray96")
  )

# Species coefficients ─────────────────────────────────────────────────────────────

alpha_homerange = results[["All species"]][["alpha_homerange"]]
alpha_homerange$name = factor(alpha_homerange$name,
                            levels = c("pcnt_standinit", "pcnt_thin", "pcnt_mature"),
                            labels = c("Stand Initiation", "Thinned", "Mature")
)

alpha_homerange$sp_name = str_to_sentence(alpha_homerange$species)
sp_order <- alpha_homerange |>
  filter(name == "Stand Initiation") |>
  arrange(mean) |>
  pull(sp_name)
alpha_homerange = alpha_homerange %>%  mutate(
  sp_name = factor(sp_name, levels = sp_order),
  pt_alpha = if_else(overlap0 == "1", 0.35, 1.0))

fig_species = ggplot(alpha_homerange, aes(y = sp_name, color = name)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, alpha = pt_alpha), width = 0, linewidth = 0.4) +
  geom_point(aes(x = mean, alpha = pt_alpha), shape = 16, size  = 1.2) +
  scale_x_continuous(expand = expansion(mult = 0.01), breaks = c(-1, 0, 1)) +
  scale_y_discrete(expand = expansion(add = c(0.8, 0.5))) +
  scale_alpha_identity() +
  scale_color_manual(values = stage_colors) +
  facet_wrap(~ name, nrow = 1) +
  labs(x = "Standardized coefficient estimate", y = NULL) +
  theme(
    panel.border = element_blank(),
    axis.line.x  = element_line(color = "grey70"),
    axis.line.y  = element_line(color = "grey70"),
    axis.text.y = element_text(size = 7),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "none"
  ); print(fig_species)

# Season effect

ggplot(results[["All species"]][["alpha_season"]] %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(title = "season") + theme(legend.position = "bottom")

ggplot(results[["All species"]][["alpha_point"]], aes(x = mean, y = species)) +
  facet_wrap(~name) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(title = "point") + theme(legend.position = "bottom")

# Plot-level effects

mu_alpha_plot = results[["All species"]][["mu_alpha_plot"]]
mu_alpha_plot$group = "All species"
mu_alpha_plot$stage = factor(mu_alpha_plot$stratum_4, levels = c("standinit", "compex", "thin", "mature"),
                             labels = c("Stand Initiation", "Competitive Exclusion", "Thinned", "Mature"))
mu_alpha_plot$name = factor(mu_alpha_plot$name, levels = c("bap_hwd_mean", "qmd_6_mean", "tree_acre_6_mean"),
                            labels = c("Proportion hardwood", "Quadratic mean diameter", "Tree density"))
mu_alpha_plot = mu_alpha_plot %>%
  mutate(pt_alpha = ifelse(`25%` < 0 & `75%` > 0, 0.3, 1))

ggplot(mu_alpha_plot, aes(x = mean, y = fct_rev(name), color = stage, alpha = pt_alpha)) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(position = position_dodge(width = 0.5, reverse = TRUE), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, linewidth = 0.4, position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbarh(aes(xmin = `25%`, xmax = `75%`), width = 0, linewidth = 0.8, position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_color_manual(values = stage_colors) +
  scale_alpha_identity(guide = "none") +
  labs(x = "Standardized coefficient estimate", y = "", color = "Management stage") +
  theme(
    panel.border = element_blank(),
    axis.line.x  = element_line(color = "grey70"),
    axis.line.y  = element_line(color = "grey70"),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
  )

ggplot(results[["Nesting guild"]][["mu_alpha_plot"]], aes(x = mean, y = name, color = group, shape = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))

ggplot(results[["Foraging guild"]][["mu_alpha_plot"]], aes(x = mean, y = name, color = group, shape = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))

# Species-specific within-stage habitat effects
ggplot(results[["All species"]][["alpha_plot"]] %>%
         filter(stratum_4 == "mature", name == "bap_hwd_mean") %>% mutate(species = fct_reorder(species, mean)),
       aes(x = mean, y = species, color = group, shape = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))
ggplot(results[["All species"]][["alpha_plot"]] %>%
         filter(stratum_4 == "standinit", name == "qmd_6_mean") %>% mutate(species = fct_reorder(species, mean)),
       aes(x = mean, y = species, color = group, shape = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))
ggplot(results[["All species"]][["alpha_plot"]] %>%
         filter(stratum_4 == "mature", name == "tree_acre_6_mean") %>% mutate(species = fct_reorder(species, mean)),
       aes(x = mean, y = species, color = group, shape = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))


ggplot(results[["All species"]][["alpha_point"]] %>%
         filter(name == "dist_road_paved") %>% mutate(species = fct_reorder(species, mean)),
       aes(x = mean, y = species, color = group)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))


# Detection effects

ggplot(results[["All species"]][["mu_beta_point"]],
       aes(x = mean, y = parameter, color = group)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))

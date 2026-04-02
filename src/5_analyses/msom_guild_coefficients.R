source("src/global.R")

path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% mutate(species = common_name)

models = data.frame(
  model = c(
    "All species",
    "Nesting strategy"
  ),
  path = c(
    "data/cache/models/V4_msom_V4_nofp_nofp_all.rds",
    "data/cache/models/V4_msom_V4_nofp_nofp_nest_ps.rds"
  )
)

results = setNames(
  lapply(models$model, function(x) list()),
  models$model
)
for (m in seq_along(models)) {
  
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
}

stop("DEBUGGY")

### Site and homerange-level effects

mu_alpha_site_all = rbind(
  results[["All species"]][["mu_alpha_homerange"]],
  results[["All species"]][["mu_alpha_point"]],
  results[["All species"]][["mu_alpha_season"]]
)
mu_alpha_site_all$group = "All species"

mu_alpha_site_nesting = rbind(
  results[["Nesting strategy"]][["mu_alpha_homerange"]],
  results[["Nesting strategy"]][["mu_alpha_point"]],
  results[["Nesting strategy"]][["mu_alpha_season"]]
)
mu_alpha_site_nesting$group = str_to_sentence(mu_alpha_site_nesting$group)

mu_alpha_site = rbind(mu_alpha_site_all, mu_alpha_site_nesting)
mu_alpha_site$group = factor(
  mu_alpha_site$group,
  levels = c("All species", "Cavity_p", "Cavity_s", "Tree", "Shrub", "Ground", "Burrow", "Cliff")
)
mu_alpha_site$name <- factor(mu_alpha_site$name, levels = rev(unique(mu_alpha_site$name)))

ggplot(mu_alpha_site %>% filter(name %in% c("pcnt_standinit", "pcnt_thin", "pcnt_mature")), aes(x = mean, y = name, color = fct_rev(group))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 1) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  scale_color_brewer(palette = "Dark2", name = "Guild")

ggplot(mu_alpha_site %>% filter(name %in% c("elevation", "dist_road_paved", "dist_watercourse_major", "season")), aes(x = mean, y = name, color = fct_rev(group))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 1) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  scale_color_brewer(palette = "Dark2", name = "Guild")

alpha_homerange = left_join(results[["All species"]][["alpha_homerange"]], species_traits, by = "species")
ggplot(alpha_homerange %>% filter(name == "pcnt_standinit") %>%
         mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group_nest_ps)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(title = "Homerange % Stand Initiation") + theme(legend.position = "bottom")
ggplot(alpha_homerange %>% filter(name == "pcnt_thin") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group_nest_ps)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(title = "Homerange % Thinned") + theme(legend.position = "bottom")
ggplot(alpha_homerange %>% filter(name == "pcnt_mature") %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group_nest_ps)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(title = "Homerange % Mature") + theme(legend.position = "bottom")

alpha_season = left_join(results[["All species"]][["alpha_season"]], species_traits, by = "species")
ggplot(alpha_season %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group_nest_ps)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(title = "season") + theme(legend.position = "bottom")

# alpha_point = left_join(results[["All species"]][["alpha_point"]], species_traits, by = "species")
# ggplot(alpha_season %>% mutate(species = fct_reorder(species, mean)), aes(x = mean, y = species, color = group_nest_ps)) +
#   geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
#   geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
#   labs(title = "season") + theme(legend.position = "bottom")

### Plot-level effects

mu_alpha_plot = results[["All species"]][["mu_alpha_plot"]]
mu_alpha_plot$group = "All species"
ggplot(mu_alpha_plot, aes(x = mean, y = name, color = stratum_4, shape = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("forestgreen", "tan4", "orange", "purple"))

ggplot(results[["Nesting strategy"]][["mu_alpha_plot"]], aes(x = mean, y = name, color = group, shape = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.1, position = position_dodge(width = 0.5))

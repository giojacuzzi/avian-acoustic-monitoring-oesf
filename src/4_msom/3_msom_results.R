# 3_msom_results_multimodel.R ####################################################################################
# Inspect MSOM results across model variants
#
# CONFIG:
path_msom = "data/cache/models/msom_fp_fp_diet_2026-02-13_18:14:56.rds"
# "data/cache/models/msom_groups_multiseason_fp_Miller_habitat_association_2025-11-30_13:01:50.rds"
# "data/cache/models/msom_all_2026-02-11_19:12:08.rds"
#
# OUTPUT:
out_cache_dir  = "data/cache/4_msom/3_msom_results"
#
# INPUT:
path_y          = "data/cache/4_msom/1_assemble_msom_data/y.rds"
path_trait_data = "data/cache/trait_data/trait_data.csv"
##################################################################################################################

source("src/global.R")

# if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)

message("Loading species detection histories 'y' from ", path_y)
y = readRDS(path_y)

# Inspect group membership ---------------------------------------------------------------------

ggplot(groups %>% left_join(species_traits, by = "common_name"), aes(x = 1, fill = group_nest)) + geom_bar(position = "fill")
ggplot(groups %>% left_join(species_traits, by = "common_name"), aes(x = 1, fill = group_nest)) + geom_bar(position = "fill")
ggplot(groups %>% left_join(species_traits, by = "common_name"), aes(x = 1, fill = group_status)) + geom_bar(position = "fill")

# Visualize results ---------------------------------------------------------------------

## Naive occupancy

y[ , , , "american crow"]

occurrence_per_season <- apply(y, c("site", "season", "species"), function(x) any(x >= 1, na.rm = TRUE))
occurrence_per_season_sum <- apply(occurrence_per_season, c("season", "species"), sum)
total_occurrences <- apply(occurrence_per_season_sum, "species", sum)
total_occurrences_tibble <- tibble(
  species = names(total_occurrences),
  sites_present = as.numeric(total_occurrences)
)
print(total_occurrences_tibble, n = Inf)

species_per_site_season <- apply(occurrence_per_season, c("site", "season"), sum)

## Retrieve baseline responses

community_baselines <- model_data$msom_summary %>%
  mutate(param_base = str_remove(param, "\\[\\d+\\]")) %>%   # strip index first
  filter(param_base %in% c("mu.u", "sigma.u", "mu.v", "sigma.v",
                           "mu.w", "sigma.w", "mu.b", "sigma.b")) %>%
  mutate(group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer()) %>%
  select(param, group_idx, prob, prob_lower95, prob_upper95)

community_baselines = community_baselines %>% left_join(groups %>% distinct(group_idx, group), by = "group_idx")

message("Baseline occurrence probability:")
print(subset(community_baselines, grepl("^mu.u", param)))

message("Baseline unconfirmed detection probability:")
print(subset(community_baselines, grepl("^mu.v", param)))

message("Baseline unconfirmed false positive detection probability:")
print(subset(community_baselines, grepl("^mu.w", param)))

message("Baseline confirmed true positive probability:")
print(subset(community_baselines, grepl("^mu.b", param)))

ggplot(community_baselines %>% filter(startsWith(param, "mu")), aes(x = prob, y = param, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = prob_lower95, xmax = prob_upper95), width = 0.2, position = position_dodge(width = 0.5)) +
  lims(x = c(0, 1)) +
  labs(x = "Posterior probability", y = "Parameter",
       title = "Hyperparameter estimates across models (mean and 95% BCI)",
       subtitle = "
u - occupancy probability across sites
v - true positive detection probability given presence
w - false positive detection probability
b - confirmed true positive probability
") + theme(legend.position = "bottom")

species_baselines = model_data$msom_summary %>%
  filter(grepl("^[uvwb]\\[", param)) %>%
  select(param, prob, prob_lower95, prob_upper95) %>%
  mutate(species_idx = str_extract(param, "(?<=\\[)[0-9]+") %>% as.integer()) %>%
  mutate(common_name = species[species_idx]) %>%
  left_join(groups, by = "common_name")

message("Species occurrence probability range:")
species_baselines %>% filter(startsWith(param, "u")) %>% summarise(
  prob_min = min(prob),
  species_min = common_name[which.min(prob)],
  prob_max = max(prob),
  species_max = common_name[which.max(prob)],
  .groups = "drop"
)

ggplot(species_baselines %>% filter(startsWith(param, "u")), aes(x = prob, y = reorder(common_name, prob), color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index",
       title = "Species-specific occurrence probability `u` across models (mean and 95% BCI)"
  ) + theme(legend.position = "bottom")

# Compare occupancy parameters --------------------------------------------------

message("Summarizing occupancy parameters")

# Get occupancy hyperparameters
param_families <- c("alpha", "delta")
combined_occ_effects <- map_dfr(c(path_msom), function(m) {
  
  groups_m <- groups %>% arrange(common_name)
  
  map_dfr(param_families, function(pf) {
    
    param_data <- model_data[[paste0("param_", pf, "_data")]]
    
    param_pattern <- c(paste0("mu.", pf), paste0("sigma.", pf))
    
    coef_summary <- msom_summary %>%
      filter(Reduce(`|`, lapply(
        param_pattern,
        \(p) str_starts(param, p)
      ))) %>%
      mutate(
        group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
        param = str_remove(param, "\\[\\d+\\]"),
        param_family = pf,
      ) %>%
      select(param_family, param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
      left_join(groups_m %>% distinct(group_idx, group), by = "group_idx")
    
    effect_sizes <- full_join(
      coef_summary, # %>% filter(str_starts(param, "mu")),
      param_data %>%
        mutate(param = paste0("mu.", param),
               param_family = pf,
               model = m),
      by = c("param", "param_family")
    )
    
    effect_sizes
  })
})

# Get means
combined_occ_effects %>% filter(str_starts(param, "mu"))

ggplot(combined_occ_effects %>% filter(str_starts(param, "mu")),
       aes(x = mean, y = name, group = fct_rev(group), color = factor(group))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`), width = 0, linewidth = 1, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, position = position_dodge(width = 0.5))

param_occ_data = bind_rows(
  model_data$param_alpha_data,
  model_data$param_delta_data,
  model_data$param_season_data
)

species_effects = param_occ_data %>%
  mutate(coef = map(param, ~ model_data$msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
  unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(common_name = species[species_idx])

species_effects = species_effects %>% left_join(species_traits, by ="common_name") %>%
  left_join(groups, by = "common_name")


# Inspect a specific parameter
p = ggplot(species_effects %>% filter(name == "prop_abund_lsog"),
           aes(x = coef_mean, y = reorder(common_name, coef_mean), color = group)) +
  geom_vline(xintercept = 0, color = "gray80") + geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) +
  geom_point(); print(p)

p = ggplot(species_effects %>% filter(name == "prop_abund_standinit"),
           aes(x = coef_mean, y = reorder(common_name, coef_mean), color = group)) +
  geom_vline(xintercept = 0, color = "gray80") + geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) +
  geom_point(); print(p)

p = ggplot(species_effects %>% filter(name == "prop_abund_comthin"),
           aes(x = coef_mean, y = reorder(common_name, coef_mean), color = group)) +
  geom_vline(xintercept = 0, color = "gray80") + geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) +
  geom_point(); print(p)

ggplot(species_effects, aes(x = coef_mean, y = name, group = group, color = group)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0.5, color = "gray90", linetype = "dashed") +
  geom_vline(xintercept = -0.5, color = "gray90", linetype = "dashed") +
  geom_beeswarm(aes(shape = coef_overlap0), dodge.width = 0.5, cex = 0.2, priority = "density", size = 0.5, alpha = 0.6) +
  geom_errorbar(data = combined_occ_effects, aes(x = mean, y = as.factor(name), color = as.factor(group), xmin = `2.5%`, xmax = `97.5%`), linewidth = 0.5, width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = combined_occ_effects, aes(x = mean, y = as.factor(name), color = as.factor(group), xmin = `25%`, xmax = `75%`), linewidth = 1, width = 0, position = position_dodge(width = 0.5)) +
  geom_point(data = combined_occ_effects, aes(x = mean, y = as.factor(name), color = as.factor(group), group = as.factor(group)), size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 1)) +
  labs(x = "Occupancy coefficient mean") +
  theme_sleek() + guides(shape = "none")

# Compare detection parameters --------------------------------------------------

param_detect_data = bind_rows(
  model_data$param_beta_data,
  model_data$param_season_data
)

detect_coef_summary = model_data$msom_summary %>%
  filter(Reduce(`|`, lapply(
    c("mu.beta", "sigma.beta"),
    \(p) str_starts(param, p)
  ))) %>%
  mutate(
    group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param = str_remove(param, "\\[\\d+\\]")
  ) %>%
  select(param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
  left_join(groups %>% distinct(group_idx, group), by = "group_idx")

detect_effect_sizes = full_join(detect_coef_summary %>% filter(str_starts(param, "mu")), param_detect_data %>% mutate(param = paste0("mu.", param)), by='param')

detect_species_effects = param_detect_data %>%
  mutate(coef = map(param, ~ model_data$msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
  unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(common_name = species[species_idx])
detect_species_effects = detect_species_effects %>% left_join(species_traits, by ="common_name") %>%
  left_join(groups, by = "common_name")

binwidth = 0.04
detect_species_effects$coef_bin = round(detect_species_effects$coef_mean / binwidth) * binwidth
p_detect = ggplot(detect_species_effects, aes(x = coef_bin, y = name, group = group)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_errorbar(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`, color = as.factor(group)), size = 0.5, width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `25%`, xmax = `75%`, color = as.factor(group)), size = 1.0, width = 0, position = position_dodge(width = 0.5)) +
  geom_point(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group), group = as.factor(group)), size = 3, position = position_dodge(width = 0.5)) +
  geom_beeswarm(aes(color = group, alpha = coef_f),
                shape = 16, dodge.width = 0.5, cex = 1, priority = "density", size = 1.1, side=1L) +
  # geom_text_repel(data = detect_species_effects, aes(x = coef_mean, y = name, label = common_name, color = group), size = 1, direction = "y", hjust = 0.05, max.overlaps = 20, position = position_dodge(0.5)) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  labs(x = "Detection coefficient mean", y = "Predictor", color = "Diet group", alpha = "coef_f") +
  theme_sleek() + guides(shape = "none") + theme(
    legend.position = "none"
  ); print(p_detect)

p_yday = ggplot(detect_species_effects %>% filter(name == "yday"), aes(x = coef_mean, y = reorder(common_name, coef_mean), color = group)) +
  geom_vline(xintercept = 0, color = "gray80") + geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) + geom_point(); print(p_yday)

# Marginal responses --------------------------------------------------------------------

# Single dimension "alpha" covariate
effect_name = "homerange_qmd_mean"

message("Marginal responses for ", effect_name)

param_data = param_occ_data %>% filter(name == effect_name)
param_name = param_data %>% pull(param)

# Mean marginal probabilities of occurrence for the metacommunity in relation to alpha
coefs = model_data$msom_summary %>% filter(stringr::str_starts(param, param_name)) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(common_name = species[species_idx])
param_scaled_data = param_data %>% pull(data)
p_mu = attr(param_scaled_data[[1]], "scaled:center") # to transform between z-scale and pred_range_original scale
p_sd = attr(param_scaled_data[[1]], "scaled:scale")
bound_low  = min(param_scaled_data[[1]]) * p_sd + p_mu
bound_high = max(param_scaled_data[[1]]) * p_sd + p_mu
pred_range_original = seq(bound_low, bound_high, length.out = 100) # range of possible alpha values
pred_range_standardized = (pred_range_original - p_mu) / p_sd

mu_u_samples      = as.matrix(model_data$msom$sims.list[["mu.u"]])
mu_param_samples  = as.matrix(model_data$msom$sims.list[[paste0("mu.", param_name)]])
n_groups = ncol(mu_u_samples)

meta_summary <- map_dfr(seq_len(n_groups), function(g) {
  
  eta <- outer(
    mu_param_samples[, g],
    pred_range_standardized
  ) + mu_u_samples[, g]
  
  psi <- plogis(eta)
  
  tibble(
    group_idx = g,
    idx = pred_range_original,
    psi_mean  = colMeans(psi),
    psi_lower = matrixStats::colQuantiles(psi, probs = 0.025),
    psi_upper = matrixStats::colQuantiles(psi, probs = 0.975)
  )
})

# Predict species-specific occurrence probabilities
# species-specific occurrence intercepts u[i]
intercepts = model_data$msom_summary %>% filter(str_starts(param, "u")) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+")), common_name = species[species_idx]) %>%
  select(common_name, u_i = mean)
# species-specific coefficients
intercepts_and_coeffs = coefs %>% rename(alpha6_i = mean) %>% select(common_name, alpha6_i) %>%
  left_join(intercepts, by = "common_name")
predictions = intercepts_and_coeffs %>%
  rowwise() %>% do({
    i <- .
    tibble(
      common_name = i$common_name,
      idx = pred_range_original,
      psi = plogis(i$u_i + i$alpha6_i * pred_range_standardized)
    )
  }) %>% bind_rows()
predictions = predictions %>% left_join(groups %>% select(common_name, group), by = "common_name")

meta_summary = meta_summary %>%
  left_join(
    groups %>% select(group_idx, group) %>% distinct(),
    by = "group_idx"
  )

ggplot() +
  geom_line(data = predictions, aes(x = idx, y = psi, group = common_name, color = group), alpha = 0.2) +
  geom_ribbon(data = meta_summary, aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group, group = group), alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
  # scale_x_continuous(limits = c(bound_low, bound_high)) +
  facet_wrap(~ group) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = param_data$name, y = "Occurrence probability")


# Multiple dimension "delta" covariate
effect_name = "prop_abund_standinit"

message("Marginal responses for ", effect_name)
s = "american robin"

param_data = param_occ_data %>% filter(name == effect_name)
param_name = param_data %>% pull(param)

# Mean marginal probabilities of occurrence for the metacommunity in relation to alpha
coefs = model_data$msom_summary %>% filter(stringr::str_starts(param, param_name)) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(common_name = species[species_idx])
param_scaled_data = param_data %>% pull(data)
p_mu = attr(param_scaled_data[[1]][[s]], "scaled:center") # to transform between z-scale and pred_range_original scale
p_sd = attr(param_scaled_data[[1]][[s]], "scaled:scale")
bound_low  = min(param_scaled_data[[1]][[s]]) * p_sd + p_mu
bound_high = max(param_scaled_data[[1]][[s]]) * p_sd + p_mu
pred_range_original = seq(bound_low, bound_high, length.out = 100) # range of possible alpha values
pred_range_standardized = (pred_range_original - p_mu) / p_sd

mu_u_samples      = as.matrix(model_data$msom$sims.list[["mu.u"]])
mu_param_samples  = as.matrix(model_data$msom$sims.list[[paste0("mu.", param_name)]])
n_groups = ncol(mu_u_samples)

meta_summary <- map_dfr(seq_len(n_groups), function(g) {
  
  eta <- outer(
    mu_param_samples[, g],
    pred_range_standardized
  ) + mu_u_samples[, g]
  
  psi <- plogis(eta)
  
  tibble(
    group_idx = g,
    idx = pred_range_original,
    psi_mean  = colMeans(psi),
    psi_lower = matrixStats::colQuantiles(psi, probs = 0.025),
    psi_upper = matrixStats::colQuantiles(psi, probs = 0.975)
  )
})

# Predict species-specific occurrence probabilities
# species-specific occurrence intercepts u[i]
intercepts = model_data$msom_summary %>% filter(str_starts(param, "u")) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+")), common_name = species[species_idx]) %>%
  select(common_name, u_i = mean)
# species-specific coefficients
intercepts_and_coeffs = coefs %>% rename(alpha6_i = mean) %>% select(common_name, alpha6_i) %>%
  left_join(intercepts, by = "common_name")
predictions = intercepts_and_coeffs %>%
  rowwise() %>% do({
    i <- .
    tibble(
      common_name = i$common_name,
      idx = pred_range_original,
      psi = plogis(i$u_i + i$alpha6_i * pred_range_standardized)
    )
  }) %>% bind_rows()
predictions = predictions %>% left_join(groups %>% select(common_name, group), by = "common_name")

meta_summary = meta_summary %>%
  left_join(
    groups %>% select(group_idx, group) %>% distinct(),
    by = "group_idx"
  )

ggplot() +
  geom_line(data = predictions, aes(x = idx, y = psi, group = common_name, color = group), alpha = 0.2) +
  geom_ribbon(data = meta_summary, aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group, group = group), alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
  # scale_x_continuous(limits = c(bound_low, bound_high)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = param_data$name, y = "Occurrence probability")

ggplot() +
  geom_line(data = predictions, aes(x = idx, y = psi, group = common_name, color = group), alpha = 0.2) +
  geom_ribbon(data = meta_summary, aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group, group = group), alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
  # scale_x_continuous(limits = c(bound_low, bound_high)) +
  facet_wrap(~ group) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = param_data$name, y = "Occurrence probability")

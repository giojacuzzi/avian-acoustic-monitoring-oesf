# "To determine support for the presence of false positives in the data sets we compared the full parameterization to one where false positive detections were assumed not to occur by making the constraint that p10 = 0."

library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(theme_bw())

####################################################################################################################
# Load data for each model

species_metadata = read.csv("data/traits/species_metadata(included).csv", nrows = 107) %>% select(common_name, scientific_name, home_range_radius_m, residency, habitat_association, habitat_ebird) %>% mutate(species_name = tolower(common_name))

model_paths = c(
  "data/cache/models/nofp_2025-09-22.rds",
  # "data/cache/models/fp_RoyleLink_2025-09-22.rds",
  # "data/cache/models/fp_Miller_2025-09-19.rds"
  "data/cache/models/multiseason_2025-09-23.rds"
)

model_data = list()
model_fits = data.frame()
for (path in model_paths) {
  message("Loading data for model ", path)
  msom_results = readRDS(path)
  
  model_data[[path]] = list(
    msom_summary     = transform(msom_results$msom_summary, model = path),
    param_alpha_data = msom_results$param_alpha_data,
    param_delta_data = msom_results$param_delta_data,
    param_beta_data  = msom_results$param_beta_data,
    sites            = msom_results$sites,
    species          = msom_results$species
  )
  
  model_fits = rbind(model_fits, data.frame(
    model = tools::file_path_sans_ext(basename(path)),
    n_sites = length(msom_results$sites),
    p_se = msom_results$p.se,
    p_dev = msom_results$p.dev,
    combined_mean_auc = msom_results$auc_combined['mean_auc'],
    combined_mean_auc_lower95 = msom_results$auc_combined['lower95_auc'],
    combined_mean_auc_upper95 = msom_results$auc_combined['upper95_auc'],
    species_mean_auc = mean(msom_results$auc_species$mean_auc, na.rm = TRUE),
    species_min_auc  = min(msom_results$auc_species$mean_auc, na.rm = TRUE),
    species_max_auc  = max(msom_results$auc_species$mean_auc, na.rm = TRUE)
  ))
}
model_fits = model_fits %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

####################################################################################################################
# Compare model fits and predictive performance

message("Comparing model fits and predictive performance")

plt = ggplot(model_fits, aes(x = model, y = p_se)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, color = "gray50") +
  geom_hline(yintercept = c(0.25, 0.75), color = "red", linetype = "dashed") +
  ylim(0.0, 1) +
  labs(title = "Bayesian p-value (squared error)"); print(plt + plot_annotation(title = ""))

plt = ggplot(model_fits, aes(x = model, y = p_dev)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, color = "gray50") +
  geom_hline(yintercept = c(0.25, 0.75), color = "red", linetype = "dashed") +
  ylim(0.0, 1) +
  labs(title = "Bayesian p-value (deviance)"); print(plt + plot_annotation(title = ""))

plt = ggplot(model_fits, aes(x = model, y = mean_auc)) +
  geom_errorbar(aes(ymin = lower95_auc, ymax = upper95_auc), width = 0.1) +
  geom_point(size = 2) +
  ylim(0.8, 1) +
  labs(title = "Combined ROC AUC (95% BCI)"); print(plt + plot_annotation(title = ""))

plt = ggplot(model_fits, aes(x = model, y = species_mean_auc)) +
  geom_errorbar(aes(ymin = species_min_auc, ymax = species_max_auc), width = 0.1) +
  geom_point(size = 2) +
  ylim(0.5, 1) +
  labs(title = "Species-specific ROC AUC (range)"); print(plt + plot_annotation(title = ""))

####################################################################################################################
# Compare latent occurrence estimates

# TODO

####################################################################################################################
# Compare baseline community hyperparameter estimates for occupancy and TP/FP detection probabilities

message("Comparing model fits and predictive performance")

model_summaries = do.call(rbind, lapply(model_data, `[[`, "msom_summary"))

community_baselines = model_summaries %>%
  filter(param %in% c("mu.u", "sigma.u",
                      "mu.v", "sigma.v",
                      "mu.w", "sigma.w",
                      "mu.b", "sigma.b")) %>%
  select(model, param, prob, prob_lower95, prob_upper95)
rownames(community_baselines) = NULL
print(community_baselines)

message("Baseline occurrence probability across models:")
subset(community_baselines, param == "mu.u")

message("Baseline TP detection probability across models:")
subset(community_baselines, param == "mu.v")

message("Baseline FP detection probability across models:")
subset(community_baselines, param == "mu.w")

message("Baseline confirmation detection probability across models:")
subset(community_baselines, param == "mu.b")

ggplot(community_baselines, aes(x = prob, y = param, shape = model, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = prob_lower95, xmax = prob_upper95), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Parameter",
    shape = "Model", color = "Model",
    title = "Hyperparameter estimates across models (mean and 95% BCI)",
    subtitle = "NOTE: Ignore `w` for no FP model
u - occupancy probability across sites
v - true positive detection probability given presence
w - false positive detection probability given absence
b - certain confirmation probability given detection
")

species = model_data[[1]]$species

species_baselines = model_summaries %>%
  filter(grepl("^[uvwb]\\[", param)) %>%
  select(model, param, prob, prob_lower95, prob_upper95) %>%
  mutate(species_idx = str_extract(param, "(?<=\\[)[0-9]+") %>% as.integer()) %>%
  mutate(species_name = species[species_idx])

message("Species occurrence probability range:")
species_baselines %>% filter(startsWith(param, "u")) %>% group_by(model) %>% summarise(
  prob_min = min(prob),
  species_min = species_name[which.min(prob)],
  prob_max = max(prob),
  species_max = species_name[which.max(prob)],
  .groups = "drop"
)

# # Aggregate min and max per species and model
# range_prob <- aggregate(prob ~ species + model, data = u_params, 
#                         FUN = function(x) c(min = min(x), max = max(x)))
# 
# # Fix column names (aggregate creates a matrix in the column)
# range_prob <- data.frame(
#   species = range_prob$species,
#   model = range_prob$model,
#   prob_min = range_prob$prob[, "min"],
#   prob_max = range_prob$prob[, "max"]
# )

ggplot(species_baselines %>% filter(startsWith(param, "u")), aes(x = prob, y = species_name, shape = model, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index", shape = "Model", color = "Model",
    title = "Species-specific occupancy probability `u` across models (mean and 95% BCI)"
  )

ggplot(species_baselines %>% filter(startsWith(param, "v")), aes(x = prob, y = species_name, shape = model, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index", shape = "Model", color = "Model",
       title = "Species-specific TP detection probability given presence `v` across models (mean and 95% BCI)"
  )

ggplot(species_baselines %>% filter(startsWith(param, "w")), aes(x = prob, y = species_name, shape = model, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index", shape = "Model", color = "Model",
       title = "Species-specific FP detection probability given absence `w` across models (mean and 95% BCI)"
  )

ggplot(species_baselines %>% filter(startsWith(param, "b")), aes(x = prob, y = species_name, shape = model, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index", shape = "Model", color = "Model",
       title = "Species-specific certain confirmation probability given detection `b` across models (mean and 95% BCI)"
  )

message("Species false positive probability range:")
species_baselines %>% filter(startsWith(param, "w")) %>% group_by(model) %>% summarise(
  prob_min = min(prob),
  species_min = species_name[which.min(prob)],
  prob_max = max(prob),
  species_max = species_name[which.max(prob)],
  .groups = "drop"
)


####################################################################################################################
# Compare estimated species richness per site

for (path in model_paths) {
  message(path)
  msom_summary = model_data[[path]]$msom_summary
  sites = model_data[[path]]$sites
  Nsite_posterior = msom_summary %>% filter(stringr::str_starts(param, "Nsite")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
    mutate(site_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(site = sites[site_idx])
  Nsite_mean = mean(Nsite_posterior$mean)
  message("Mean estimated species richness across all sites: ", round(Nsite_mean,1), " (range ", round(min(Nsite_posterior$mean),1), "–", round(max(Nsite_posterior$mean),1), ")")
  ggplot(Nsite_posterior, aes(x = as.factor(plot_order), y = mean)) +
    geom_hline(yintercept = Nsite_mean, linetype = "solid", color = "blue") +
    geom_point() + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
    scale_x_discrete(labels = Nsite_posterior$site) + 
    labs(title = "Estimated species richness per site", x = "Site", y = "Estimated species richness") +
    coord_flip()
  ggplot(Nsite_posterior, aes(x = mean)) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = Nsite_mean, color = "blue") +
    labs(title = "Estimated species richness per site", x = "Number of species estimated", y = "Number of sites")
}

####################################################################################################################
# Compare coefficients on occupancy

message("Comparing model coefficients on occupancy")

for (path in model_paths) {
  message("Getting occupancy coefs for model ", path)
  
  msom_summary = model_data[[path]]$msom_summary
  occ_coef_summary = msom_summary %>%
    filter(str_detect(param, "^(mu|sigma)\\.(alpha|delta)")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)
  
  detect_coef_summary = msom_summary %>%
    filter(str_detect(param, "^(mu|sigma)\\.beta")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)
  
  param_alpha_data = model_data[[path]]$param_alpha_data
  param_delta_data = model_data[[path]]$param_delta_data
  
  # Compare community level effects of each covariate on occurrence
  param_occ_data = param_alpha_data
  if (!is.null(param_delta_data)) {
    param_occ_data = rbind(param_alpha_data %>% rename(data = scaled), param_delta_data)
  }
  
  occ_effect_sizes = full_join(occ_coef_summary %>% filter(str_starts(param, "mu")), param_occ_data %>% mutate(param = paste0("mu.", param)), by='param')
  
  # plt = ggplot(occ_effect_sizes, aes(x = mean, y = as.factor(name))) +
  #   geom_vline(xintercept = 0, color = "gray") +
  #   geom_point(aes(color = overlap0)) +
  #   geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  #   geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
  #   scale_color_manual(values = c("black", "gray")) +
  #   xlim(c(-1, 1)) +
  #   labs(title = "Community level effect sizes for occurrence covariates",
  #        subtitle = tools::file_path_sans_ext(basename(path)),
  #        x = "Coefficient estimate", y = "Parameter")
  # print(plt)

  # Compare species level effects of each covariate on occurrence
  species_effects = param_occ_data %>%
    mutate(coef = map(param, ~ msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
    unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(species_name = species[species_idx])
  species_effects = species_effects %>% mutate(species_name = recode(species_name, "pacific-slope flycatcher" = "western flycatcher"))
  species_effects = full_join(species_effects, species_metadata, by = c('species_name'))
  species_effects = species_effects %>% filter(!is.na(name) & name != "")
  
  plt = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = species_effects, aes(x = coef_mean, y = name, color = habitat_association, shape = coef_overlap0), position = position_jitter(height = 0.2), size = 3, alpha = 0.95) +
    # geom_dotplot(data = species_effects, aes(x = coef_mean, y = name, fill = name, alpha = coef_overlap0), color = "transparent", binwidth = 0.03, stackdir = "center") +
    scale_color_manual(values = c("darkgray", "#d8c18a", "#b2675e", "#3c8273")) + #c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', '#b2675e', 'darkgray', '#6495ed')
    geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`), width = 0) +
    geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), shape = overlap0), size = 3.5) +
      xlim(c(-2.5, 2.5)) +
      labs(title = "Effect sizes for occurrence covariates at community and species levels",
           subtitle = tools::file_path_sans_ext(basename(path)),
           x = "Coefficient estimate", y = "Parameter") +
    geom_text_repel(
      data = species_effects %>% filter(coef_overlap0 == 0),
      aes(x = coef_mean, y = name, label = species_name, color = habitat_association),
      size = 3, nudge_x = 0.05, direction = "y", hjust = 0.05
    ); print(plt)
}

####################################################################################################################
# Compare coefficients on TP detection

for (path in model_paths) {
  message("Getting TP detection coefs for model ", path)
  
  msom_summary = model_data[[path]]$msom_summary
  detect_coef_summary = msom_summary %>%
    filter(str_detect(param, "^(mu|sigma)\\.(beta)")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)
  
  param_beta_data = model_data[[path]]$param_beta_data
  
  # Compare community level effects of each covariate on TP detection
  param_detect_data = param_beta_data
  
  detect_effect_sizes = full_join(detect_coef_summary %>% filter(str_starts(param, "mu")), param_detect_data %>% mutate(param = paste0("mu.", param)), by='param')
  
  # Compare species level effects of each covariate on TP detection
  species_effects = param_detect_data %>%
    mutate(coef = map(param, ~ msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
    unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(species_name = species[species_idx])
  species_effects = species_effects %>% mutate(species_name = recode(species_name, "pacific-slope flycatcher" = "western flycatcher"))
  species_effects = full_join(species_effects, species_metadata, by = c('species_name'))
  species_effects = species_effects %>% filter(!is.na(name) & name != "")
  
  plt = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = species_effects, aes(x = coef_mean, y = name, color = habitat_association, shape = coef_overlap0), position = position_jitter(height = 0.2), size = 3, alpha = 0.95) +
    # geom_dotplot(data = species_effects, aes(x = coef_mean, y = name, fill = name, alpha = coef_overlap0), color = "transparent", binwidth = 0.03, stackdir = "center") +
    scale_color_manual(values = c("darkgray", "#d8c18a", "#b2675e", "#3c8273")) + #c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', '#b2675e', 'darkgray', '#6495ed')
    geom_errorbar(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`), width = 0) +
    geom_point(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), shape = overlap0), size = 3.5) +
    xlim(c(-2.5, 2.5)) +
    labs(title = "Effect sizes for TP detection covariates at community and species levels",
         subtitle = tools::file_path_sans_ext(basename(path)),
         x = "Coefficient estimate", y = "Parameter") +
    geom_text_repel(
      data = species_effects %>% filter(coef_overlap0 == 0),
      aes(x = coef_mean, y = name, label = species_name, color = habitat_association),
      size = 3, nudge_x = 0.05, direction = "y", hjust = 0.05
    ); print(plt)
}

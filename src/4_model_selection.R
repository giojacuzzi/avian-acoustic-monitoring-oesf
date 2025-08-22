library(dplyr)
library(ggplot2)
theme_set(theme_classic())

model_paths = c(
  "data/cache/models/msom_HS_2025-08-21.rds",
  "data/cache/models/msom_RS_2025-08-21.rds",
  "data/cache/models/msom_RSalt_2025-08-21.rds",
  "data/cache/models/msom_HR_2025-08-21.rds",
  "data/cache/models/msom_RSHR_2025-08-21.rds"
)

model_fits = data.frame()
for (path in model_paths) {
  msom_results = readRDS(path)
  msom_summary = msom_results$msom_summary
  param_alpha_data = msom_results$param_alpha_data
  param_beta_data = msom_results$param_beta_data
  sites = msom_results$sites
  species = msom_results$species
  
  model_fits = rbind(model_fits, data.frame(
    model = tools::file_path_sans_ext(basename(path)),
    n_sites = length(msom_results$sites),
    p_val = msom_results$p_val,
    combined_mean_auc = msom_results$auc_combined['mean_auc'],
    combined_mean_auc_lower95 = msom_results$auc_combined['lower95_auc'],
    combined_mean_auc_upper95 = msom_results$auc_combined['upper95_auc'],
    species_mean_auc = mean(msom_results$auc_species$mean_auc, na.rm = TRUE),
    species_min_auc  = min(msom_results$auc_species$mean_auc, na.rm = TRUE),
    species_max_auc  = max(msom_results$auc_species$mean_auc, na.rm = TRUE)
  ))
  
  # Community-level summaries of the hyper-parameters for the detection and occupancy covariates
  # mu is the community response (mean across species) to a given covariate and sd is the standard deviation (among species). Thus, the hyper-parameters are simply the mean and variance for each covariate as measured across species (Kéry & Royle 2009)
  message("Community-level summaries of hyper-parameters for occurrence and detection covariates:")
  occurrence_coeff_summary = msom_summary %>%
    filter(str_detect(param, "^mu\\.alpha|^sigma\\.alpha")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)
  detection_coeff_summary = msom_summary %>%
    filter(str_detect(param, "^mu\\.beta|^sigma\\.beta")) %>% arrange(param) %>%
    select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)
  print(occurrence_coeff_summary)
  print(detection_coeff_summary)
  
  # Compare community level effect sizes for occurrence coefficients
  occurrence_effect_sizes = full_join(occurrence_coeff_summary %>% filter(str_starts(param, "mu")), param_alpha_data %>% mutate(param = paste0("mu.", param)), by='param')
  plt = ggplot(occurrence_effect_sizes, aes(x = mean, y = as.factor(name))) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(aes(color = overlap0)) +
    geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
    scale_color_manual(values = c("black", "gray")) +
    labs(title = "Community level effect sizes for occurrence covariates",
         subtitle = tools::file_path_sans_ext(basename(path)),
         x = "Coefficient estimate", y = "Parameter"); print(plt)
}
model_fits = model_fits %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

theme_set(theme_bw())

plt = ggplot(model_fits, aes(x = model, y = p_val)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, color = "gray50") +
  geom_hline(yintercept = c(0.25, 0.75), color = "red", linetype = "dashed") +
  ylim(0.0, 1) +
  labs(title = "Bayesian p-value"); print(plt)

plt = ggplot(model_fits, aes(x = model, y = mean_auc)) +
  geom_errorbar(aes(ymin = lower95_auc, ymax = upper95_auc), width = 0.1) +
  geom_point(size = 2) +
  ylim(0.9, 1) +
  labs(title = "Combined ROC AUC (95% BCI)"); print(plt)

plt = ggplot(model_fits, aes(x = model, y = species_mean_auc)) +
  geom_errorbar(aes(ymin = species_min_auc, ymax = species_max_auc), width = 0.1) +
  geom_point(size = 2) +
  ylim(0.5, 1) +
  labs(title = "Species-specific ROC AUC (range)"); print(plt)

stop("Read interpretation below")

# >>> Remote sensing data, despite exhibiting measurement error and a restricted set of candidate variables for modeling (due to collinearity among variables), can nevertheless serve as a practical surrogate for field-measured data for our purposes, as they capture complementary signals for the strongest effects of local-scale forest structure on bird occurrence and have virtually equivalent model fit and predictive performance.

# >>> The context of the surrounding landscape provides comparable, if not more, explanatory power than local-scale forest structure.

# >>> Local fine-scale habitat features are important in influencing patterns of species occurrence, but accounting for the context of the surrounding landscape reveals that community-level habitat use is more strongly structured by landscape configuration and composition –– namely patch area, edge density, and forest heterogeneity. The effect of fine-scale structure may be subsumed by landscape context, at least in driving habitat use at the community level. This is consistent with the idea that occurrence depends on the multi-scale resource requirements of birds across their home range, rather than only at the local point of sampling.

# >>> Structural and compositional heterogeneity drives metacommunity assembly. Strong positive influences of disturbance and edge-related metrics at the community level, and no discernible effect of late-successional or old-growth abundance, may reflect how forest management practices foster a metacommunity biased towards generalist and early seral adapted species, disadvantaging interior and mature-forest specialists. As such, these specialists are underrepresented in the aggregate effects (and MSOM groupings by guild should be explored). This highlights a potential trade-off between managing for high occurrence (i.e. occupancy) across the community and conserving sufficient habitat for mature forest specialists.

# NOTES:
# It could be the case that some mature-forest obligates are missing from LSOG patches because of a combination of limited patch area and high landscape fragmentation, reducing connectivity and increasing edge.
# Early seral habitat is ephemeral, while mature forest is not. How does this play into the conservation implications of our results?
# Conventional management strategy may foster high taxonomic diversity and overall bird occurrence, but at the expense of mature-forest specialists.

############################################################################################################

# Visualize estimated number of species per site
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

# As covariates on occupancy were standardized, the inverse-logit (`plogis()`) of u[i] is the baseline occurrence probability for species i at a site with ‘average’ habitat characteristics.

occurrence_prob = msom_summary %>% filter(stringr::str_starts(param, "u")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
mean_occurrence_prob = msom_summary %>% filter(param == "mu.u")
message("Mean species occurrence probability: ", round(mean_occurrence_prob$prob,2), " (95% BCI ", round(mean_occurrence_prob$prob_lower95,2), "–", round(mean_occurrence_prob$prob_upper95,2), ")")
message("Species occurrence probability range: ", round(min(occurrence_prob$prob),2), "–", round(max(occurrence_prob$prob),2), " (", occurrence_prob %>% slice_min(prob, n=1) %>% pull(species_name), ", ", occurrence_prob %>% slice_max(prob, n=1) %>% pull(species_name), ")")
ggplot(occurrence_prob, aes(x = as.factor(plot_order), y = prob)) +
  geom_hline(yintercept = mean_occurrence_prob$prob,         linetype = "solid", color = "blue") +
  geom_hline(yintercept = mean_occurrence_prob$prob_lower95, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mean_occurrence_prob$prob_upper95, linetype = "dashed", color = "blue") +
  geom_point() + geom_errorbar(aes(ymin = `prob_lower95`, ymax = `prob_upper95`), width = 0) +
  scale_x_discrete(labels = occurrence_prob$species_name) + 
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Baseline occurrence probability", x = "Species", y = "Occurrence probability") +
  coord_flip()

# Similarly, the inverse-logit of v[i] is the detection probability for species i under 'average' detection conditions.
detection_prob = msom_summary %>% filter(stringr::str_starts(param, "v")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
mean_detection_prob = msom_summary %>% filter(param == "mu.v")
message("Mean species detection probability: ", round(mean_detection_prob$prob,2), " (95% BCI ", round(mean_detection_prob$prob_lower95,2), "–", round(mean_detection_prob$prob_upper95,2), ")")
message("Species detection probability range: ", round(min(detection_prob$prob),2), "–", round(max(detection_prob$prob),2), " (", detection_prob %>% slice_min(prob, n=1) %>% pull(species_name), ", ", detection_prob %>% slice_max(prob, n=1) %>% pull(species_name), ")")
ggplot(detection_prob, aes(x = as.factor(plot_order), y = prob)) +
  geom_hline(yintercept = mean_detection_prob$prob,         linetype = "solid", color = "blue") +
  geom_hline(yintercept = mean_detection_prob$prob_lower95, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mean_detection_prob$prob_upper95, linetype = "dashed", color = "blue") +
  geom_point() + geom_errorbar(aes(ymin = `prob_lower95`, ymax = `prob_upper95`), width = 0) +
  scale_x_discrete(labels = detection_prob$species_name) + 
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Baseline detection probability", x = "Species", y = "Detection probability") +
  coord_flip()

# Compare species level effects of each covariate on occurrence
for (alpha_param in param_alpha_data$param) {
  alpha_name = param_alpha_data %>% filter(param == alpha_param) %>% pull(name)
  alpha_coef = msom_summary %>% filter(str_detect(param, paste0("^", alpha_param, "(?!\\d)", "\\["))) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
    mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
  mu_alpha_summary = msom_summary %>% filter(param == paste0("mu.", alpha_param)) %>% select(mean, `2.5%`, `97.5%`)
  plt = ggplot(alpha_coef, aes(x = as.factor(plot_order), y = mean)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_point(aes(color = overlap0)) +
    geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = overlap0), width = 0) +
    geom_hline(yintercept = mu_alpha_summary$mean,    linetype = "solid",  color = "blue") +
    geom_hline(yintercept = mu_alpha_summary$`2.5%`,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = mu_alpha_summary$`97.5%`, linetype = "dashed", color = "blue") +
    labs(title = paste("Effect of", alpha_name, "on occurrence"), x = "Species", y = paste(alpha_param, "coefficient estimate")) +
    scale_x_discrete(labels = alpha_coef$species_name) +
    scale_color_manual(values = c("0" = "black", "1" = "gray")) +
    coord_flip() +
    theme(legend.position = "none"); print(plt)
}

# Compare community level effect sizes for detection coefficients
detection_effect_sizes = full_join(detection_coeff_summary %>% filter(str_starts(param, "mu")), param_beta_data %>% mutate(param = paste0("mu.", param)), by='param')
plt = ggplot(detection_effect_sizes, aes(x = mean, y = as.factor(name))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = overlap0)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
  scale_color_manual(values = c("black", "gray")) +
  labs(title = "Community level effect sizes for detection covariates", x = "Coefficient estimate", y = "Parameter"); print(plt)

# Compare species level effects of each covariate on detection
for (beta_param in param_beta_data$param) {
  beta_name = param_beta_data %>% filter(param == beta_param) %>% pull(name)
  beta_coef = msom_summary %>% filter(str_detect(param, paste0("^", beta_param, "(?!\\d)", "\\["))) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
    mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(species_name = species[species_idx])
  mu_beta_summary = msom_summary %>% filter(param == paste0("mu.", beta_param)) %>% select(mean, `2.5%`, `97.5%`)
  plt = ggplot(beta_coef, aes(x = as.factor(plot_order), y = mean)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_point(aes(color = overlap0)) +
    geom_errorbar(aes(ymin = `25%`,  ymax = `75%`,   color = overlap0), width = 0, size = 1) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = overlap0), width = 0) +
    geom_hline(yintercept = mu_beta_summary$mean,    linetype = "solid",  color = "blue") +
    geom_hline(yintercept = mu_beta_summary$`2.5%`,  linetype = "dashed", color = "blue") +
    geom_hline(yintercept = mu_beta_summary$`97.5%`, linetype = "dashed", color = "blue") +
    labs(title = paste("Effect of", beta_name, "on detection"), x = "Species", y = paste(beta_param, "coefficient estimate")) +
    scale_x_discrete(labels = beta_coef$species_name) +
    scale_color_manual(values = c("black", "gray")) +
    coord_flip() +
    theme(legend.position = "none"); print(plt)
}


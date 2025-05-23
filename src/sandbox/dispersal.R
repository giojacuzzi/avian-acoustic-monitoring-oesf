# Predict home range size from body mass and hand-wing index from generalized linear models
# fit home range size data for a subset of species

library(tidyverse)
library(ggplot2)
library(ggrepel)
theme_set(theme_minimal())

path_species_list = "models/ensemble/ensemble_species_list.txt"
species_list = data.frame(
  scientific_name = sort(sapply(strsplit(readLines(path_species_list), "_"), `[`, 1)),
  common_name = sort(sapply(strsplit(readLines(path_species_list), "_"), `[`, 2))
)

# https://app.box.com/integrations/officeonline/openOfficeOnline?fileId=1728867025786&sharedAccessCode=
metadata = read.csv('/Users/giojacuzzi/Downloads/species_metadata(included).csv') %>% select(common_name, scientific_name, home_range_radius_m, residency)

avonet = readxl::read_xlsx('data/traits/AVONET Supplementary dataset 1.xlsx', sheet = "AVONET2_eBird") %>%
  janitor::clean_names() %>% rename(scientific_name = species2, family = family2, order = order2) %>% filter(scientific_name %in% species_list$scientific_name)

trait_data = left_join(avonet, metadata, by = "scientific_name")
# Manually resolve common name discrepancies
trait_data = trait_data %>% mutate(common_name = if_else(scientific_name == "Accipiter gentilis", "Northern Goshawk", common_name))
trait_data = trait_data %>% mutate(common_name = if_else(scientific_name == "Glaucidium gnoma", "Northern Pygmy-Owl", common_name))
trait_data = trait_data %>% mutate(common_name = if_else(scientific_name == "Tyto alba", "Barn Owl", common_name))

unknown_home_ranges = trait_data %>% filter(is.na(home_range_radius_m))

# Function to predict values with 95% CI with a given GLM
predict_with_ci = function(model, newdata) {
  predicted = predict(model, newdata = newdata, type = "response")
  link_pred = predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  critical_value = qnorm(0.975)
  return(data.frame(
    predicted = predicted,
    lower = exp(link_pred$fit - critical_value * link_pred$se.fit),
    upper = exp(link_pred$fit + critical_value * link_pred$se.fit)
  ))
}

###################################################################################

# Fit a generalized linear model relating natural log of body mass to home range radius
glm_mass = glm(home_range_radius_m ~ log(mass), data = trait_data, family = Gamma(link = "log"))
glm_mass_pseudoR2 = round(1 - (glm_mass$deviance / glm_mass$null.deviance), 2)

# Calculate regression line
prediction_data = data.frame(mass = exp(seq(log(min(trait_data$mass)), log(max(trait_data$mass)), length.out = 100)))
mass_predictions = predict_with_ci(glm_mass, prediction_data) %>% mutate(mass = prediction_data$mass)

ggplot() +
  geom_point(data = trait_data, aes(x = (mass), y = (home_range_radius_m))) +
  geom_text_repel(data = trait_data, aes(x = (mass), y = (home_range_radius_m), label = common_name), size = 3) +
  geom_line(data = mass_predictions, aes(x = (mass), y = (predicted)), color = "blue") +
  geom_ribbon(data = mass_predictions, aes(x = (mass), ymin = (lower), ymax = (upper)), alpha = 0.2, fill = "blue") +
  scale_x_log10() + scale_y_log10() +
  labs(
    x = "Body mass (g)", y = "Home range radius (m)",
    title = paste0("Body mass and approximate home range radius (pseudo R2=", glm_mass_pseudoR2, ")")
  )

# Predict home range radii from mass
prediction_data = data.frame(mass = unknown_home_ranges$mass)
predicted_home_range_by_mass = predict_with_ci(glm_mass, prediction_data) %>% mutate(
  mass = prediction_data$mass,
  common_name = unknown_home_ranges$common_name
)
ggplot(predicted_home_range_by_mass, aes(x = mass, y = predicted)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point(size = 2) +
  scale_x_log10() + scale_y_continuous() +
  labs(
    x = "Body mass (g)", y = "Predicted home range radius (m)",
    title = paste0("Home range predictions from body mass")
  )

###################################################################################

# Fit a generalized linear model relating natural log of hand-wing index to home range radius
glm_hwi = glm(home_range_radius_m ~ log(hand_wing_index), data = trait_data, family = Gamma(link = "log"))
glm_hwi_pseudoR2 = round(1 - (glm_hwi$deviance / glm_hwi$null.deviance), 2)

# Calculate regression line
prediction_data = data.frame(hand_wing_index = exp(seq(log(min(trait_data$hand_wing_index)), log(max(trait_data$hand_wing_index)), length.out = 100)))
hwi_predictions = predict_with_ci(glm_hwi, prediction_data) %>% mutate(hand_wing_index = prediction_data$hand_wing_index)

ggplot() +
  geom_point(data = trait_data, aes(x = (hand_wing_index), y = (home_range_radius_m))) +
  geom_text_repel(data = trait_data, aes(x = (hand_wing_index), y = (home_range_radius_m), label = common_name), size = 3) +
  geom_line(data = hwi_predictions, aes(x = (hand_wing_index), y = (predicted)), color = "blue") +
  geom_ribbon(data = hwi_predictions, aes(x = (hand_wing_index), ymin = (lower), ymax = (upper)), alpha = 0.2, fill = "blue") +
  scale_x_log10() + scale_y_log10() +
  labs(
    x = "Hand-wing index", y = "Home range radius (m)",
    title = paste0("Hand-wing index and approximate home range radius (pseudo R2=", glm_hwi_pseudoR2, ")")
  )

# Predict home range radii from hand-wing-index
prediction_data = data.frame(hand_wing_index = unknown_home_ranges$hand_wing_index)
predicted_home_range_by_hwi = predict_with_ci(glm_hwi, prediction_data) %>% mutate(
  hand_wing_index = prediction_data$hand_wing_index,
  common_name = unknown_home_ranges$common_name
)
ggplot(predicted_home_range_by_hwi, aes(x = hand_wing_index, y = predicted)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point(size = 2) +
  scale_x_log10() + scale_y_continuous() +
  labs(
    x = "Hand-wing index", y = "Predicted home range radius (m)",
    title = paste0("Home range predictions from hand-wing index")
  )

###################################################################################

home_range_sizes = full_join(predicted_home_range_by_mass %>% mutate(home_range_radius_m = predicted), trait_data %>% filter(!is.na(home_range_radius_m)) %>% select(common_name, home_range_radius_m, mass)) %>% arrange(home_range_radius_m)

home_range_sizes$range_bin <- cut(home_range_sizes$home_range_radius_m,
                    breaks = quantile(home_range_sizes$home_range_radius_m, probs = seq(0, 1, length.out = 5), na.rm = TRUE),
                    include.lowest = TRUE,
                    labels = c("Very Small", "Small", "Medium", "Large"))

median_hr = round(median(home_range_sizes$home_range_radius_m), 0)
mean_hr = round(mean(home_range_sizes$home_range_radius_m), 0)
max_hr = round(max(home_range_sizes$home_range_radius_m), 0)
min_hr = round(min(home_range_sizes$home_range_radius_m), 0)

ggplot(home_range_sizes, aes(x = "",y = home_range_radius_m)) +
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black") +
  scale_y_log10() +
  labs(title = "Distribution of home range radii (m)",
       subtitle = paste0("Median ", median_hr, "\nMean ", mean_hr, "\nMax ", max_hr, "\nMin ", min_hr),
       x = "", y = "Estimated home range radius (m)")

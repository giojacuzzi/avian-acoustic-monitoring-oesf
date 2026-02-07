## 1_agg_traits.R #########################################################################################
# Aggregate species functional traits into a single dataframe
#
## OUTPUT:
path_out = "data/cache/trait_data/trait_data.csv"
#
## INPUT:
path_avonet      = "data/traits/AVONET Supplementary dataset 1.xlsx"
path_eltontraits = "data/traits/EltonTraits 1.0/BirdFuncDat.csv"
path_diets       = "data/traits/species_diets.csv"
path_guilds      = "data/traits/species_guilds.csv"
path_metadata    = "data/traits/species_metadata.csv"
###########################################################################################################

source("src/global.R")

# Combine species trait databases -------------------------------------------------------------------------

# Note that the following subspecies traits are inferred from their superspecies in AVONET:
# - American Goshawk "Astur atricapillus" (i.e Eurasian/Northern Goshawk "Accipiter gentilis")
# - Northern Pygmy-owl "Glaucidium californicum" (i.e. Northern/Mountain Pygmy-owl "Glaucidium gnoma")
# - American Barn Owl "Tyto furcata" (i.e. Western Barn Owl "Tyto alba")
message("Loading AVONET from ", path_avonet)
avonet = readxl::read_xlsx(path_avonet, sheet = "AVONET2_eBird") %>%
  janitor::clean_names() %>%
  rename(scientific_name = species2, family = family2, order = order2) %>%
  mutate(scientific_name = tolower(scientific_name)) %>%
  filter(scientific_name %in% class_labels$scientific_name) %>%
  select(scientific_name, family, order, trophic_level, trophic_niche, primary_lifestyle, habitat_density, migration, mass, hand_wing_index)

message("Loading EltonTraits 1.0 from ", path_eltontraits)
eltontraits = read_csv(path_eltontraits, show_col_types = FALSE) %>% clean_names() %>%
  rename(scientific_name = scientific, common_name = english) %>%
  mutate(
    common_name = str_to_lower(common_name), scientific_name = str_to_lower(scientific_name)
  ) %>% mutate(
    common_name = case_when(
      common_name == "grey jay"                    ~ "canada jay",
      common_name == "american black swift"        ~ "black swift",
      common_name == "common starling"             ~ "european starling",
      common_name == "american treecreeper"        ~ "brown creeper",
      common_name == "winter wren"                 ~ "pacific wren",
      common_name == "black-throated grey warbler" ~ "black-throated gray warbler",
      common_name == "common teal"                 ~ "green-winged teal",
      TRUE ~ common_name   # retain all others
    )
  ) %>% filter(common_name %in% class_labels$common_name) %>%
  select(common_name, diet_inv, diet_5cat, body_mass_value) %>% rename(mass_et = body_mass_value)
# Check for any unmatched names
# setdiff(class_labels$common_name, eltontraits %>% filter(common_name %in% class_labels$common_name) %>% pull(common_name))

message("Loading guilds from ", path_guilds)
guilds = read_csv(path_guilds, show_col_types = FALSE) %>% clean_names() %>%
  mutate(common_name = tolower(common_name)) %>% select(common_name, foraging_guild_cornell)

message("Loading metadata from ", path_metadata)
metadata = read.csv(path_metadata, nrows = 107) %>% clean_names() %>%
  mutate(common_name = str_to_lower(common_name)) %>% select(common_name, home_range_radius_m)

message("Combining species trait databases")
trait_data = left_join(class_labels,   avonet,  by = "scientific_name")
trait_data = left_join(trait_data, eltontraits, by = "common_name")
trait_data = left_join(trait_data, guilds,      by = "common_name")
trait_data = left_join(trait_data, metadata,    by = "common_name")

trait_data = trait_data %>% filter(!str_starts(label, "abiotic|biotic")) # retain only avian pecies, no other classes

# All grouping (full community):
trait_data = trait_data %>% mutate(
  group_all = "all"
)

# TODO: Nesting guild
trait_data = trait_data %>% mutate(
  group_nest = "TODO"
)

# TODO: Dietary guild
# "Following Pigot et al. (2020), we assigned all species to nine trophic niches (frugivore; granivore; nectarivore; terrestrial herbivore; aquatic herbivore; invertivore; vertivore; aquatic predator; scavenger) encompassing major resource types utilised by birds."
trait_data = trait_data %>% mutate(
  group_diet = str_to_lower(trophic_niche)
)

# TODO: Foraging behavior
# "Next, we classified each species into five lifestyles (or domains) according to their predominant locomotory niche while foraging: aerial, insessorial, terrestrial, aquatic and generalist. This is a separate dimension to diet inasmuch as species eating fish may be aquatic (e.g. pelican), aerial (e.g. tern), terrestrial (e.g. heron) or insessorial (e.g. kingfisher). Insessorial denotes a perching lifestyle, including arboreal species, but also any species habitually perching on other substrates, including cliffs or manmade structures."
trait_data = trait_data %>% mutate(
  group_forage_behavior = str_to_lower(primary_lifestyle)
)

# TODO: Foraging habitat / substrate
trait_data = trait_data %>% mutate(
  group_forage_substrate = case_when(
    (foraging_guild_cornell %in% c("aerial forager", "soaring", "flycatching", "hovering", "aerial dive")) ~ "aerial",
    (foraging_guild_cornell %in% c("foliage gleaner"))                                                     ~ "gleaner",
    TRUE                                                                                                   ~ "TODO"
  )
)

# TODO: Habitat association (e.g. early seral vs old-forest)
trait_data = trait_data %>% mutate(
  group_habitat = "TODO"
)

# TODO: WADNR species of concern
trait_data = trait_data %>% mutate(
  group_wadnr_concern = "TODO"
)

# TODO: Conservation status
trait_data = trait_data %>% mutate(
  group_status = "TODO"
)

# Predict home range size from traits ------------------------------------------------------------------

message("Predicting home range from traits")

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

unknown_home_ranges = trait_data %>% filter(is.na(home_range_radius_m))

## Fit a generalized linear model relating natural log of body mass to home range radius
glm_mass = glm(home_range_radius_m ~ log(mass), data = trait_data, family = Gamma(link = "log"))
glm_mass_pseudoR2 = round(1 - (glm_mass$deviance / glm_mass$null.deviance), 2)

# Calculate regression line
prediction_data = data.frame(mass = exp(seq(log(min(trait_data$mass)), log(max(trait_data$mass)), length.out = 100)))
mass_predictions = predict_with_ci(glm_mass, prediction_data) %>% mutate(mass = prediction_data$mass)

p_mass = ggplot() +
  geom_point(data = trait_data, aes(x = (mass), y = (home_range_radius_m))) +
  geom_text_repel(data = trait_data, aes(x = (mass), y = (home_range_radius_m), label = common_name), size = 3) +
  geom_line(data = mass_predictions, aes(x = (mass), y = (predicted)), color = "blue") +
  geom_ribbon(data = mass_predictions, aes(x = (mass), ymin = (lower), ymax = (upper)), alpha = 0.2, fill = "blue") +
  scale_x_log10() + scale_y_log10() +
  labs(
    x = "Body mass (g)", y = "Home range radius (m)",
    title = paste0("Body mass and approximate home range radius (pseudo R2=", glm_mass_pseudoR2, ")")
  ); print(p_mass)

# Predict home range radii from mass
prediction_data = data.frame(mass = unknown_home_ranges$mass)
predicted_home_range_by_mass = predict_with_ci(glm_mass, prediction_data) %>% mutate(
  mass = prediction_data$mass,
  common_name = unknown_home_ranges$common_name
)
p_mass_pred = ggplot(predicted_home_range_by_mass, aes(x = mass, y = predicted)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point(size = 2) +
  scale_x_log10() + scale_y_continuous() +
  labs(
    x = "Body mass (g)", y = "Predicted home range radius (m)",
    title = paste0("Home range predictions from body mass")
  ); print(p_mass_pred)

## Fit a generalized linear model relating hand-wing index to home range radius
glm_hwi = glm(home_range_radius_m ~ hand_wing_index, data = trait_data, family = Gamma(link = "log"))
glm_hwi_pseudoR2 = round(1 - (glm_hwi$deviance / glm_hwi$null.deviance), 2)

# Calculate regression line
prediction_data = data.frame(hand_wing_index = seq(min(trait_data$hand_wing_index), max(trait_data$hand_wing_index), length.out = 100))
hwi_predictions = predict_with_ci(glm_hwi, prediction_data) %>% mutate(hand_wing_index = prediction_data$hand_wing_index)

p_hwi = ggplot() +
  geom_point(data = trait_data, aes(x = (hand_wing_index), y = (home_range_radius_m))) +
  geom_text_repel(data = trait_data, aes(x = (hand_wing_index), y = (home_range_radius_m), label = common_name), size = 3) +
  geom_line(data = hwi_predictions, aes(x = (hand_wing_index), y = (predicted)), color = "blue") +
  geom_ribbon(data = hwi_predictions, aes(x = (hand_wing_index), ymin = (lower), ymax = (upper)), alpha = 0.2, fill = "blue") +
  scale_x_continuous() + scale_y_log10() +
  labs(
    x = "Hand-wing index", y = "Home range radius (m)",
    title = paste0("Hand-wing index and approximate home range radius (pseudo R2=", glm_hwi_pseudoR2, ")")
  ); print(p_hwi)

# Predict home range radii from hand-wing-index
prediction_data = data.frame(hand_wing_index = unknown_home_ranges$hand_wing_index)
predicted_home_range_by_hwi = predict_with_ci(glm_hwi, prediction_data) %>% mutate(
  hand_wing_index = prediction_data$hand_wing_index,
  common_name = unknown_home_ranges$common_name
)
p_hwi_pred = ggplot(predicted_home_range_by_hwi, aes(x = hand_wing_index, y = predicted)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point(size = 2) +
  scale_x_continuous() + scale_y_continuous() +
  labs(
    x = "Hand-wing index", y = "Predicted home range radius (m)",
    title = paste0("Home range predictions from hand-wing index")
  ); print(p_hwi_pred)

## Fit a generalized linear model relating natural log of body mass and hand-wing index to home range radius
glm_combined = glm(home_range_radius_m ~ log(mass) + hand_wing_index, data = trait_data, family = Gamma(link = "log"))
glm_combined_pseudoR2 = round(1 - (glm_combined$deviance / glm_combined$null.deviance), 2)

summary(glm_combined)

# Predict home range radii from mass and hwi
prediction_data = data.frame(mass = unknown_home_ranges$mass, hand_wing_index = unknown_home_ranges$hand_wing_index)
predicted_home_range_by_both = predict_with_ci(glm_combined, prediction_data) %>% mutate(
  mass = prediction_data$mass,
  hand_wing_index = prediction_data$hand_wing_index,
  common_name = unknown_home_ranges$common_name
)

# Model selection ---------------------------------------------------------------------

message("Model selection results:")
print(AIC(glm_mass, glm_hwi, glm_combined))
summary(glm_mass)$deviance
summary(glm_hwi)$deviance
summary(glm_combined)$deviance
glm_mass_pseudoR2
glm_hwi_pseudoR2
glm_combined_pseudoR2
summary(predicted_home_range_by_mass$predicted)
summary(predicted_home_range_by_hwi$predicted)
summary(predicted_home_range_by_both$predicted)

# Cache results ---------------------------------------------------------------------

message("Imputing home range radii from mass model")
trait_data = trait_data %>%
  left_join(predicted_home_range_by_mass %>% select(common_name, predicted, lower, upper), by = c("common_name")) %>%
  mutate(home_range_radius_m = ifelse(is.na(home_range_radius_m), predicted, home_range_radius_m))

median_hr = round(median(trait_data$home_range_radius_m), 0)
mean_hr   = round(mean(trait_data$home_range_radius_m), 0)
max_hr    = round(max(trait_data$home_range_radius_m), 0)
min_hr    = round(min(trait_data$home_range_radius_m), 0)

p_hrr = ggplot(trait_data, aes(x = "",y = home_range_radius_m)) +
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black") +
  scale_y_log10() +
  labs(title = "Distribution of home range radii (m)",
       subtitle = paste0("Median ", median_hr, "\nMean ", mean_hr, "\nMax ", max_hr, "\nMin ", min_hr),
       x = "", y = "Estimated home range radius (m)"); print(p_hrr)

message("Creating a categorical body size group")
trait_data$group_size = cut(trait_data$home_range_radius_m,
                           breaks = quantile(trait_data$home_range_radius_m, probs = seq(0, 1, length.out = 5), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c("very small", "small", "medium", "large"))

# Cache results ---------------------------------------------------------------------

message("Trait data aggregation complete:")
print(trait_data, n = Inf)

if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
write.csv(trait_data, path_out)
message(crayon::green("Cached trait data:", path_out))

####################################################################################
# MSOM coefficient estimates
#
# CONFIG:
path_msom = "data/cache/models/prefinal_msom_jags_nofp_all.rds"
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

param_alpha_point_data = model_data$param_alpha_data$param_alpha_point_data
param_alpha_plot_data = model_data$param_alpha_data$param_alpha_plot_data
param_alpha_homerange_data = model_data$param_alpha_data$param_alpha_homerange_data

match_s = stages %>% distinct() %>% arrange(stage_idx)
match_g = groups %>% select(group, group_idx) %>% distinct() %>% arrange(group_idx)
match_i = tibble(species = species, species_idx = 1:length(species))

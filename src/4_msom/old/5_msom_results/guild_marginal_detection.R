####################################################################################
#
# CONFIG:
path_msom = "data/cache/models/V4_msom_V4_fp_all.rds"
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

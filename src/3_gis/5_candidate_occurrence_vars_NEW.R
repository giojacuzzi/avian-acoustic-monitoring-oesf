## 4_candidate_occurrence_vars.R ########################################################################
# Finalize candidate set of variables for occurrence
#
## INPUT:
path_plot_scale_data       = "data/cache/3_gis/3_calc_occurrence_vars/V2_data_plot_scale_2020_clean_strata_4.rds"
path_homerange_scale_data  = "data/cache/3_gis/3_calc_occurrence_vars/V2_data_homerange_scale_2020_clean_strata_4.rds"
###########################################################################################################

source("src/global.R")

message('Loading field measurement data from cache ', path_plot_scale_data)
data_fieldmeasurements = readRDS(path_plot_scale_data)

message('Loading homerange scale data from cache ', path_homerange_scale_data)
data_homerange_scale = readRDS(path_homerange_scale_data)

###########################################################################################
### Collinearity analysis - local plot scale variables from remote sensing
message("Assessing collinearity among variables at local plot scale (remote sensing)")

data_plotscale_rs = data_homerange_scale[['plot']] %>% select(-buffer_radius_m, -scale)
data_plotscale_rs$stage_3 = data_fieldmeasurements$stage_3

# A priori reduce list of candidate plot scale variables (remove irrelevant and poorly performing variables)
var_candidates_plotscale_rs = data_plotscale_rs %>% select(stage_3,
  # "We excluded remote sensing layers having poor agreement (Pearson’s r < 0.6) with field measurements from previous validation efforts across state forest lands (Rickefs, [unpublished report])"
  ba_mean, ba_4_mean, ba_6_mean, ba_t100_mean,
  # ba_cv, ba_4_cv, ba_6_cv, ba_t100_cv,
  bap_hwd_mean,
  canopy_cover_mean,
  # canopy_cover_cv,
  qmd_6_mean, qmd_t100_mean,
  # qmd_6_cv, qmd_t100_cv,
  tree_acre_4_mean, tree_acre_6_mean,
  # tree_acre_4_cv, tree_acre_6_cv,
  htmax_mean,
  # htmax_cv, ht_t100_cv
)
(names(var_candidates_plotscale_rs))
pairwise_collinearity(var_candidates_plotscale_rs)
pairwise_collinearity_by_group(var_candidates_plotscale_rs, "stage_3")

# "As measurements between size classes of basal area, quadratic mean diameter, and tree density were strongly correlated, we constrained size class across these metrics to trees of at least 15 cm diameter at breast height, representing the presence of merchantable trees that are the focus of silvicultural prescriptions."
var_candidates_plotscale_rs = data_plotscale_rs %>% select(stage_3,
  ba_6_mean,
  bap_hwd_mean,
  canopy_cover_mean,
  qmd_6_mean,
  tree_acre_6_mean,
  htmax_mean
)
(names(var_candidates_plotscale_rs))
pairwise_collinearity(var_candidates_plotscale_rs)
pairwise_collinearity_by_group(var_candidates_plotscale_rs, "stage_3")

# "As the number, size, and species composition of trees are fundamental aspects of stand structure manipulated by management (Hansen et al., 1995), we retained tree density, quadratic mean diameter, and hardwood proportion as predictors at the plot scale."
var_candidates_plotscale_rs = data_plotscale_rs %>% select(stage_3,
  bap_hwd_mean,
  qmd_6_mean,
  tree_acre_6_mean
)
pairwise_collinearity(var_candidates_plotscale_rs)
pairwise_collinearity_by_group(var_candidates_plotscale_rs, "stage_3")

###########################################################################################
### Collinearity analysis - homerange scale variables (composition, configuration
message("Assessing collinearity among variables at homerange scale (composition, configuration)")

scale = 'median' # e.g. 'min', 'median', 'mean', 'max'

data_homerange_rs = data_homerange_scale[[scale]] %>% select(-buffer_radius_m, -scale)
data_homerange_rs$stage_3 = data_fieldmeasurements$stage_3

# A priori reduce list of candidate plot scale variables (remove irrelevant variables)
var_candidates_homerangescale = data_homerange_rs %>% select(-all_of(c(
  'density_streams', 'density_streams_major',
  'density_roads', 'density_roads_paved',
  'pcnt_water', 'pcnt_road_paved', 'pcnt_old', 'pcnt_underdev'
)))
(names(var_candidates_homerangescale))

pairwise_collinearity(var_candidates_homerangescale)
pairwise_collinearity_by_group(var_candidates_homerangescale, "stage_3")

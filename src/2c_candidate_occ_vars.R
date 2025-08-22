##############################################################################
# Finalize candidate set of variables for occurrence
##############################################################################
library(progress)
library(car)
library(mapview)
library(terra)
library(sf)
library(dplyr)
library(ggrepel)
theme_set(theme_minimal())

crs_m = 32610 # coordinate reference system EPSG:32610, UTM Zone 10N, (meters)

path_data_plot_scale = "data/cache/occurrence_covariates/data_plot_scale.rds"
message('Loading plot scale data from cache ', path_data_plot_scale)
data_plot_scale = readRDS(path_data_plot_scale)

path_data_homerange_scale = "data/cache/occurrence_covariates/data_homerange_scale.rds"
message('Loading homerange scale data from cache ', path_data_homerange_scale)
data_homerange_scale = readRDS(path_data_homerange_scale)

# path_rast_cover_clean = "data/cache/occurrence_covariates/rast_cover_clean.tif"
# message('Loading raster cover data from cache ', path_rast_cover_clean)
# rast_cover_clean = rast(path_rast_cover_clean)
# 
# stratum_colors = c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', '#b2675e', 'darkgray', '#6495ed')
# mapview(rast_cover_clean,
#         alpha.regions = 1.0,
#         col.regions = stratum_colors) +
#   mapview(data_plot_scale, zcol = 'stratum', col.regions = stratum_colors) +
#   mapview(layer.name = 'plot', st_buffer(data_plot_scale, unique(data_homerange_scale[['plot']]$buffer_radius_m)),
#           alpha.regions = 0.0, lwd = 2) +
#   mapview(layer.name = 'median', st_buffer(data_plot_scale, unique(data_homerange_scale[['median']]$buffer_radius_m)),
#           alpha.regions = 0.0, lwd = 2) +
#   mapview(layer.name = 'mean', st_buffer(data_plot_scale, unique(data_homerange_scale[['mean']]$buffer_radius_m)),
#           alpha.regions = 0.0, lwd = 2) +
#   mapview(layer.name = 'max', st_buffer(data_plot_scale, unique(data_homerange_scale[['max']]$buffer_radius_m)),
#           alpha.regions = 0.0, lwd = 2)

data_plot_scale = data_plot_scale %>% st_drop_geometry()

pairwise_collinearity = function(vars, threshold = 0.8) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

###########################################################################################
### Collinearity analysis - local plot scale variables from habitat survey
message("Assessing collinearity among variables at local plot scale (habitat survey)")

var_candidates_plotscale_hs = data_plot_scale %>% select(where(is.numeric), -stratum)

# A priori reduce list of candidate plot scale variables (remove irrelevant variables)
var_candidates_plotscale_hs = var_candidates_plotscale_hs %>% select(-all_of(c(
  # Age is not a variable of direct interest
  'age_point', 'age_mean', 'age_cv',
  # Individual tree species composition is better summarized with taxonomic diversity metrics
  'plot_tree_gte10cm_density_psme', 'plot_tree_gte10cm_density_thpl', 'plot_tree_gte10cm_density_abam', 'plot_tree_gte10cm_density_tshe', 'plot_tree_gte10cm_density_alru', 'plot_tree_gte10cm_density_pisi',
  'plot_tree_all_density_thpl', 'plot_tree_all_density_abam', 'plot_tree_all_density_tshe', 'plot_tree_all_density_alru', 'plot_tree_all_density_pisi', 'plot_tree_all_density_psme',
  # Drop distance from intermittent streams; favor only consistent watercourses indicative of different riparian habitat.
  'dist_watercourse_all'
)))
(names(var_candidates_plotscale_hs))
pairwise_collinearity(var_candidates_plotscale_hs)

# Reduce list of candidate variables (highly correlated, less preferable, etc.)
var_candidates_plotscale_hs = var_candidates_plotscale_hs %>% select(-all_of(c(
  # plot_treeden_lt10cmDbh_hs is highly correlated with that of all trees; favor breakout of tree sizes for interpretability.
  'plot_treeden_all_hs',
  # Quadratic mean diameter is highly correlated with that of all trees; favor all trees QMD for interpretability.
  'plot_qmd_gt10cmDbh_hs', 'plot_qmd_lt10cmDbh_hs',
  # Choose between understory cover and volume. We keep volume here to better represent vertical structure.
  'plot_understory_cover',
  # Most taxonomic diversity metrics are highly correlated between tree size classes; favor all sizes for interpretability.
  'plot_tree_lt10cm_richness',  'plot_tree_gte10cm_richness',
  'plot_tree_lt10cm_evenness',  'plot_tree_gte10cm_evenness',
  'plot_tree_lt10cm_diversity', 'plot_tree_gte10cm_diversity',
  # Furthermore, diversity and evenness are highly correlated; favor shannon diversity for interpretability.
  'plot_tree_all_evenness', 'plot_tree_all_richness',
  # Tree height, height-to-live-crown, and length-of-live-crown are highly correlated
  'plot_hlc_hs', 'plot_llc_hs'
)))
(names(var_candidates_plotscale_hs))
pairwise_collinearity(var_candidates_plotscale_hs)

# VIF analysis for multicollinearity (consider dropping variable(s) with high VIF values (> 10))
model = lm(rep(1, nrow(var_candidates_plotscale_hs)) ~ ., data = var_candidates_plotscale_hs)
sort(vif(model))
model = lm(rep(1, nrow(var_candidates_plotscale_hs)) ~ . -plot_ba_hs, data = var_candidates_plotscale_hs)
sort(vif(model))

# RESULT: Do not include plot_ba_hs alongside plot_qmd_all_hs
(var_candidates_plotscale_hs = sort(names(sort(vif(model)))))

###########################################################################################
### Collinearity analysis - local plot scale variables from remote sensing
message("Assessing collinearity among variables at local plot scale (remote sensing)")

var_candidates_plotscale_rs = data_homerange_scale[['plot']] %>% select(where(is.numeric), -buffer_radius_m) # NOTE: 'plot' == 100m radius
var_candidates_plotscale_rs = var_candidates_plotscale_rs %>% select(starts_with("homerange_"))

# A priori reduce list of candidate plot scale variables (remove irrelevant variables)
var_candidates_plotscale_rs = var_candidates_plotscale_rs %>% select(-all_of(c(
  # Age is not a variable of direct interest
  'homerange_age_mean', 'homerange_age_cv',
  # The following measures of structural variation are not of direct interest
  'homerange_ba_cv', 'homerange_treeden_all_cv', 'homerange_treeden_gt4in_dbh_cv', 'homerange_qmd_cv',
  'homerange_canopy_layers_cv', 'homerange_snagden_gt15dbh_cv', 'homerange_downvol_cv', 'homerange_canopy_cover_cv',
  # Canopy cover is more applicable than closure
  'homerange_canopy_closure_mean', 'homerange_canopy_closure_cv'
)))
(names(var_candidates_plotscale_rs))
pairwise_collinearity(var_candidates_plotscale_rs)

# Reduce list of candidate variables (highly correlated, less preferable, etc.)
var_candidates_plotscale_rs = var_candidates_plotscale_rs %>% select(-all_of(c(
  # homerange_ba_mean ~ (homerange_qmd_mean, homerange_htmax_mean, homerange_canopy_cover_mean, homerange_snagden_gt15dbh_mean)
  'homerange_ba_mean', 'homerange_htmax_mean',
  # homerange_htmax_cv ~ homerange_treeden_gt4in_dbh_mean
  'homerange_treeden_gt4in_dbh_mean',
  # Drop canopy cover for comparison with analogous habitat survey (hs) variables
  'homerange_canopy_cover_mean',
  'homerange_canopy_layers_mean',
  # homerange_qmd_mean ~ homerange_snagden_gt15dbh_mean
  'homerange_snagden_gt15dbh_mean'
)))
(names(var_candidates_plotscale_rs))
pairwise_collinearity(var_candidates_plotscale_rs)

# VIF analysis for multicollinearity (consider dropping variable(s) with high VIF values (> 10))
model = lm(rep(1, nrow(var_candidates_plotscale_rs)) ~ ., data = var_candidates_plotscale_rs)
sort(vif(model))

# RESULT:
# At the local plot scale (100m radius), there is a higher degree of multicollinearity among structural variables measured via remote sensing than via habitat surveys
# Do not include homerange_ba_mean alongside homerange_qmd_mean or homerange_snagden_gt15dbh_mean
(var_candidates_plotscale_rs = sort(names(sort(vif(model)))))

###########################################################################################
### Collinearity analysis - homerange scale variables (composition, configuration
message("Assessing collinearity among variables at homerange scale (composition, configuration)")

scale = 'median' # e.g. 'min', 'median', 'mean', 'max'
var_candidates_homerangescale = data_homerange_scale[[scale]] %>% select(where(is.numeric), -buffer_radius_m)
# only look at composition/configuration variables
var_candidates_homerangescale = var_candidates_homerangescale[,1:(which(colnames(var_candidates_homerangescale) == "focalpatch_age_mean")-1)]

# A priori reduce list of candidate plot scale variables (remove irrelevant variables)
var_candidates_homerangescale = var_candidates_homerangescale %>% select(-all_of(c(
  # Drop cover diversity metrics in favor of forest cover diversity metrics
  'cover_diversity', 'cover_richness', 'cover_evenness',
  # Drop density streams in favor of density major streams (which are perennial and bordered by riparian habitat)
  'density_streams'
)))
(names(var_candidates_homerangescale))
pairwise_collinearity(var_candidates_homerangescale)

# TODO: re-run pairwise collinearity >= 0.8 across scales 
#
# min: 
# focalpatch_pcnt ~ focalpatch_core_pcnt, aggregation_idx
# focalpatch_core_pcnt ~ cover_forest_evenness, cover_forest_diversity, aggregation_idx, shape_idx
# cover_forest_diversity ~ cover_forest_richness, cover_forest_evenness
# prop_abund_6 ~ density_roads_paved
# Consider dropping: focalpatch_core_pcnt, aggregation_idx, cover_forest_richness, cover_forest_evenness, prop_abund_6
#
# median:
# focalpatch_pcnt ~ focalpatch_core_pcnt, aggregation_idx
# focalpatch_core_pcnt ~ cover_forest_diversity, aggregation_idx
# cover_forest_diversity ~ cover_forest_richness, cover_forest_evenness
# prop_abund_6 ~ density_roads_paved
# prop_abund_1 ~ cw_edge_density
# Consider dropping: focalpatch_core_pcnt, aggregation_idx, cover_forest_richness, cover_forest_evenness, prop_abund_6, prop_abund_1
#
# mean:
# focalpatch_pcnt ~ focalpatch_core_pcnt
# cover_forest_diversity ~ cover_forest_evenness, prop_abund_2
# aggregation_idx ~ cw_edge_density
# prop_abund_6 ~ density_roads_paved
# Consider dropping: focalpatch_core_pcnt, aggregation_idx, cover_forest_evenness, prop_abund_6, prop_abund_2
#
# max:
# focalpatch_pcnt ~ focalpatch_core_pcnt
# cover_forest_richness (no variation)
# cover_forest_diversity ~ cover_forest_evenness
# aggregation_idx ~ cw_edge_density
# prop_abund_6 ~ density_roads_paved
# prop_abund_2 ~ density_streams
# Consider dropping: focalpatch_core_pcnt, aggregation_idx, cover_forest_richness, cover_forest_evenness, prop_abund_6, prop_abund_2

# Reduce list of candidate variables (highly correlated, less preferable, etc.)
var_candidates_homerangescale = var_candidates_homerangescale %>% select(-all_of(c(
  # The following variables cause issues of collinearity across multiple homerange scales
  'focalpatch_core_area_homeange_pcnt', 'aggregation_idx',
  'cover_forest_richness', 'cover_forest_evenness', # ~ cover_forest_diversity
  'prop_abund_undstryreinit', 'prop_abund_oldgrowth', # prop_abund_lsog = prop_abund_undstryreinit + prop_abund_oldgrowth
  'prop_abund_standinit',  # ~ density_edge_cw
  # 'prop_abund_2',        # competitive exclusion
  'prop_abund_roads',      # ~ density_roads_paved
  'prop_abund_water'       # ~ density_streams_major
)))
(names(var_candidates_homerangescale))
pairwise_collinearity(var_candidates_homerangescale)

# VIF analysis for multicollinearity (consider dropping variable(s) with high VIF values (> 10))
model = lm(rep(1, nrow(var_candidates_homerangescale)) ~ ., data = var_candidates_homerangescale)
sort(vif(model))
model = lm(rep(1, nrow(var_candidates_homerangescale)) ~ . -density_roads -prop_abund_stemexcl, data = var_candidates_homerangescale)
sort(vif(model))

# RESULT:
# prop_abund_stemexcl should be dropped as baseline state
# prop_abund_standinit ~ density_edge_cw
# density_roads causes multicollinearity at large spatial scales
(var_candidates_homerangescale = sort(names(sort(vif(model)))))

###########################################################################################
### Collinearity analysis - plot scale (habitat survey) plus homerange scale
message("Assessing collinearity among remaining variables at both the plot (habitat survey) and homerange (composition, configuration) scales")

scale = 'median' # e.g. 'min', 'median', 'mean', 'max'
var_candidates_hs_combined = full_join(
  data_plot_scale %>% select(all_of(c("site", var_candidates_plotscale_hs))),
  data_homerange_scale[[scale]] %>% select(all_of(c("site", var_candidates_homerangescale))),
  by = 'site'
) %>% select(-site)
var_candidates_hs_combined = data.frame(lapply(var_candidates_hs_combined, as.numeric))

pairwise_collinearity(var_candidates_hs_combined, threshold = 0.8)

model = lm(rep(1, nrow(var_candidates_hs_combined)) ~ ., data = var_candidates_hs_combined)
sort(vif(model))

# RESULT: No multicollinearity among variables at min, median, or mean scales. At max scale, density_streams removed.
var_candidates_hs_combined = sort(names(sort(vif(model))))
message("Final list of ", ncol(var_candidates_hs_combined)," uncorrelated variables at both the plot (habitat survey) and homerange (composition, configuration) scales:")
print(var_candidates_hs_combined)

###########################################################################################
### Collinearity analysis - plot scale (remote sensing) plus homerange scale
message("Assessing collinearity among remaining variables at both the plot (remote sensing) and homerange (composition, configuration) scales")

scale = 'median' # e.g. 'plot', 'median', 'mean', 'max'
var_candidates_rs_combined = full_join(
  data_homerange_scale[['plot']] %>% select(all_of(c("site", var_candidates_plotscale_rs))), # 100m local plot scale
  data_homerange_scale[[scale]] %>% select(all_of(c("site", var_candidates_homerangescale))),
  by = 'site'
) %>% select(-site)
var_candidates_rs_combined = data.frame(lapply(var_candidates_rs_combined, as.numeric))

pairwise_collinearity(var_candidates_rs_combined, threshold = 0.8)

model = lm(rep(1, nrow(var_candidates_rs_combined)) ~ ., data = var_candidates_rs_combined)
sort(vif(model))

# RESULT: No multicollinearity among variables at any scale, however, reduced variable set compared to habitat survey set
var_candidates_rs_combined = sort(names(sort(vif(model))))
message("Final list of ", ncol(var_candidates_hs_combined)," uncorrelated variables at both the plot (remote sensing) and homerange (composition, configuration) scales:")
print(var_candidates_rs_combined)

# ### Principal component analysis ##############################################################
# 
# var_candidates_scaled = scale(na.omit(var_candidates))
# pca = prcomp(var_candidates_scaled, center = TRUE, scale. = TRUE)
# summary(pca)
# 
# # Loadings for first 3 PCs
# round(pca$rotation[, 1:3], 2)
# 
# factoextra::fviz_eig(pca)
# factoextra::fviz_pca_ind(pca, geom.ind = "point", col.ind = "cos2", repel = TRUE)
# factoextra::fviz_pca_ind(pca, geom.ind = "point", addEllipses = TRUE, legend.title = "Group", repel = TRUE)
# factoextra::fviz_pca_var(pca, col.var = "contrib", repel = TRUE)

###########################################################################################
### Collinearity analysis - plot scale (remote sensing) plus species-specific homerange scales

# TODO: pairwise r and VIF for each species, keep track of which variables are correlated

# Start with collinearity among species-specific homerange scales
names(data_homerange_scale) = tolower(names(data_homerange_scale))

plot_data_homerange_scale = data_homerange_scale[['plot']] # default minimum at 100m

# Remove irrelevant species that are not included in the analysis
species_to_keep = c(
  "american crow",
  "american goldfinch",
  "american goshawk",
  "american kestrel",
  "american robin",
  "bald eagle",
  "band-tailed pigeon",
  "barred owl",
  "belted kingfisher",
  "black-headed grosbeak",
  "black-throated gray warbler",
  "brown creeper",
  "brown-headed cowbird",
  "canada jay",
  "cedar waxwing",
  "chestnut-backed chickadee",
  "common nighthawk",
  "common raven",
  "dark-eyed junco",
  "downy woodpecker",
  "eurasian collared-dove",
  "evening grosbeak",
  "golden-crowned kinglet",
  "hairy woodpecker",
  "hammond's flycatcher",
  "hermit thrush",
  "hutton's vireo",
  "macgillivray's warbler",
  "marbled murrelet",
  "northern flicker",
  "northern pygmy-owl",
  "northern saw-whet owl",
  "olive-sided flycatcher",
  "orange-crowned warbler",
  "pacific wren",
  "pileated woodpecker",
  "pine siskin",
  "purple finch",
  "red crossbill",
  "red-breasted nuthatch",
  "red-tailed hawk",
  "ruby-crowned kinglet",
  "ruffed grouse",
  "rufous hummingbird",
  "sharp-shinned hawk",
  "song sparrow",
  "sooty grouse",
  "spotted towhee",
  "steller's jay",
  "swainson's thrush",
  "townsend's warbler",
  "varied thrush",
  "vaux's swift",
  "violet-green swallow",
  "warbling vireo",
  "western screech-owl",
  "western tanager",
  "western wood-pewee",
  "white-crowned sparrow",
  "wilson's warbler",
  "yellow warbler",
  "yellow-rumped warbler"
)
data_homerange_scale = data_homerange_scale[species_to_keep]
n_scales = length(names(data_homerange_scale))

# A priori remove some vars
candidates = c(
  "aggregation_idx",
  # "cover_diversity", # Drop cover diversity metrics for forest cover diversity metrics
  # "cover_evenness",
  # "cover_richness",
  "cover_forest_diversity",
  "cover_forest_evenness",
  # "cover_forest_richness", # No variance at highest spatial scales           
  "density_edge_cw",
  "density_roads",                     
  "density_roads_paved",
  # "density_streams", # Drop density streams for major streams (which are perennial and bordered by riparian habitat)
  "density_streams_major",             
  "focalpatch_area_homeange_pcnt",  
  "focalpatch_core_area_homeange_pcnt",
  # "focalpatch_age_cv", # Only look at landscape variables for now
  # "focalpatch_age_mean",
  # "focalpatch_ba_cv",
  # "focalpatch_ba_mean",
  # "focalpatch_canopy_closure_cv",      
  # "focalpatch_canopy_closure_mean",
  # "focalpatch_canopy_cover_cv",
  # "focalpatch_canopy_cover_mean",      
  # "focalpatch_canopy_layers_cv",
  # "focalpatch_canopy_layers_mean",
  # "focalpatch_downvol_cv",
  # "focalpatch_downvol_mean",
  # "focalpatch_htmax_cv",               
  # "focalpatch_htmax_mean",
  # "focalpatch_isolation",
  # "focalpatch_qmd_cv",                 
  # "focalpatch_qmd_mean",
  # "focalpatch_snagden_gt15dbh_cv",
  # "focalpatch_snagden_gt15dbh_mean",   
  # "focalpatch_treeden_all_cv",
  # "focalpatch_treeden_all_mean",
  # "focalpatch_treeden_gt4in_dbh_cv",   
  # "focalpatch_treeden_gt4in_dbh_mean",
  # "homerange_age_cv",
  # "homerange_age_mean",                
  # "homerange_ba_cv",
  # "homerange_ba_mean",
  # "homerange_canopy_closure_cv",       
  # "homerange_canopy_closure_mean",
  # "homerange_canopy_cover_cv",
  # "homerange_canopy_cover_mean",       
  # "homerange_canopy_layers_cv",
  # "homerange_canopy_layers_mean",
  # "homerange_downvol_cv",              
  # "homerange_downvol_mean",
  # "homerange_htmax_cv",
  # "homerange_htmax_mean",              
  # "homerange_qmd_cv",
  # "homerange_qmd_mean",
  # "homerange_snagden_gt15dbh_cv",      
  # "homerange_snagden_gt15dbh_mean",
  # "homerange_treeden_all_cv",
  # "homerange_treeden_all_mean",        
  # "homerange_treeden_gt4in_dbh_cv",
  # "homerange_treeden_gt4in_dbh_mean",
  "prop_abund_comthin",                
  "prop_abund_lsog",
  "prop_abund_oldgrowth",
  "prop_abund_roads",             
  "prop_abund_standinit",
  "prop_abund_stemexcl",
  "prop_abund_undstryreinit",     
  "prop_abund_water",
  "shape_idx"
)

all_species_results = data.frame()
for (i in 1:n_scales) {
  scale = names(data_homerange_scale)[i]
  print(scale)
  species_homerange_data = data_homerange_scale[[scale]]
  buffer_radius_m = unique(species_homerange_data$buffer_radius_m)
  # if buffer_radius_m < 100m, set to 100m and pull data from "plot" scale
  if (buffer_radius_m < 100) {
    species_homerange_data = plot_data_homerange_scale
  }
  species_var_candidates_homerange = species_homerange_data %>% select(where(is.numeric), -buffer_radius_m)
  
  # A priori reduce list of candidate plot scale variables (remove irrelevant variables)
  species_var_candidates_homerange = species_var_candidates_homerange %>% select(all_of(candidates))
  (names(species_var_candidates_homerange))
  
  sd0 = names(species_var_candidates_homerange)[apply(species_var_candidates_homerange, 2, sd, na.rm = TRUE) == 0]
  if (length(sd0) > 0) {
    message("Species ", scale, " has zero-variance variables: ", paste(vars_sd0, collapse = ", "))
  }
  
  species_results = pairwise_collinearity(species_var_candidates_homerange) %>% mutate(scale = scale)
  all_species_results = rbind(all_species_results, species_results)
}

# Make a consistent pair ID with a unique separator
all_species_results$pair = apply(
  all_species_results[,c("Var1","Var2")], 1,
  function(x) paste(sort(x), collapse = "#")
)
# Count number of scales where each pair is flagged
pair_summary = aggregate(scale ~ pair, data = all_species_results, FUN = function(x) length(unique(x)))
pair_summary$var1 = sub("#.*", "", pair_summary$pair)
pair_summary$var2 = sub(".*#", "", pair_summary$pair)
pair_summary = pair_summary[,c("var1","var2","scale")]
pair_summary$pcnt_species = round(pair_summary$scale / length(unique(all_species_results$scale)), 2)
pair_summary = pair_summary %>% arrange(desc(pcnt_species)) %>% rename(freq = scale)
print(pair_summary)

# Reduce list of candidate variables (highly correlated, less preferable, etc.)
candidates = candidates[!candidates %in% c(
  # The following variables cause issues of collinearity across many homerange scales
  'cover_forest_evenness', # ~ cover_forest_diversity
  'prop_abund_roads',                               # ~ density_roads_paved
  'prop_abund_water',                               # ~ density_streams_major
  'focalpatch_core_area_homeange_pcnt', # ~ focalpatch_area_homeange_pcnt, aggregation_idx, cover_forest_diversity
  'aggregation_idx', # ~ focalpatch_area_homeange_pcnt, focalpatch_core_area_homeange_pcnt, density_edge_cw
  'prop_abund_undstryreinit', 'prop_abund_oldgrowth', # prop_abund_lsog = prop_abund_undstryreinit + prop_abund_oldgrowth
  'prop_abund_standinit',  # ~ density_edge_cw
  'prop_abund_stemexcl' # ~ cover_forest_diversity
)]


# VIF analysis
all_species_vif_results = data.frame()
for (i in 1:n_scales) {
  scale = names(data_homerange_scale)[i]
  print(scale)
  species_homerange_data = data_homerange_scale[[scale]]
  buffer_radius_m = unique(species_homerange_data$buffer_radius_m)
  # if buffer_radius_m < 100m, set to 100m and pull data from "plot" scale
  if (buffer_radius_m < 100) {
    species_homerange_data = plot_data_homerange_scale
  }
  # Discard non-candidate variables
  species_homerange_data = species_homerange_data[, names(species_homerange_data) %in% candidates, drop = FALSE]

  model = lm(rep(1, nrow(species_homerange_data)) ~ ., data = species_homerange_data)
  results = as.data.frame(t(vif(model)))
  results = results %>% mutate (scale = scale)
  all_species_vif_results = rbind(all_species_vif_results, results)
}
mean_vif = apply(all_species_vif_results[, -ncol(all_species_vif_results)], 2, mean)
print(sort(mean_vif, decreasing = TRUE))

# Count frequency of VIF >= 10 for each var across species scales
vif_values = all_species_vif_results %>% select(-scale)
freq = colSums(vif_values >= 10)
vif_summary = data.frame(
  var  = names(freq),
  freq = freq,
  pcnt_species = round( freq / nrow(vif_values), 2),
  row.names = NULL
)
vif_summary = vif_summary %>% arrange(desc(pcnt_species))
print(vif_summary)


# TODO: Introduce plot scale (remote sensing) variables


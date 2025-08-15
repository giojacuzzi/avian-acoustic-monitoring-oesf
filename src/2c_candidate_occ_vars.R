##############################################################################
# Finalize candidate set of variables for occurrence
##############################################################################
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

path_rast_cover_clean = "data/cache/occurrence_covariates/rast_cover_clean.tif"
message('Loading raster cover data from cache ', path_rast_cover_clean)
rast_cover_clean = rast(path_rast_cover_clean)

mapview(rast_cover_clean,
        alpha.regions = 1.0,
        col.regions = c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', 'darkgray', '#6495ED')) +
  mapview(data_plot_scale) +
  mapview(layer.name = 'min', st_buffer(data_plot_scale, unique(data_homerange_scale[['min']]$buffer_radius_m)),
          alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'median', st_buffer(data_plot_scale, unique(data_homerange_scale[['median']]$buffer_radius_m)),
          alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'mean', st_buffer(data_plot_scale, unique(data_homerange_scale[['mean']]$buffer_radius_m)),
          alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'max', st_buffer(data_plot_scale, unique(data_homerange_scale[['max']]$buffer_radius_m)),
          alpha.regions = 0.0, lwd = 2)

data_plot_scale = data_plot_scale %>% st_drop_geometry()

pairwise_collinearity = function(vars, threshold = 0.8) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

### Collinearity analysis - local plot scale variables from habitat survey ##############################################################
message("Assessing collinearity among variables at local plot scale (habitat survey)")

var_candidates_plotscale_hs = data_plot_scale %>% select(where(is.numeric), -stratum)

# A priori reduce list of candidate plot scale variables (remove irrelevant variables)
var_candidates_plotscale_hs = var_candidates_plotscale_hs %>% select(-all_of(c(
  # Age is not a variable of direct interest
  'age_point', 'age_mean', 'age_cv',
  # Individual tree species composition is better summarized with taxonomic diversity metrics
  'tree_gte10cm_density_psme', 'tree_gte10cm_density_thpl', 'tree_gte10cm_density_abam', 'tree_gte10cm_density_tshe', 'tree_gte10cm_density_alru', 'tree_gte10cm_density_pisi',
  'tree_all_density_thpl', 'tree_all_density_abam', 'tree_all_density_tshe', 'tree_all_density_alru', 'tree_all_density_pisi', 'tree_all_density_psme',
  # Drop distance from intermittent streams; favor only consistent watercourses indicative of different riparian habitat.
  'dist_watercourses_all'
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
  'tree_lt10cm_richness', 'tree_gte10cm_richness',
  'tree_lt10cm_evenness', 'tree_gte10cm_evenness',
  'tree_lt10cm_diversity', 'tree_gte10cm_diversity',
  # Furthermore, diversity and evenness are highly correlated; favor shannon diversity for interpretability.
  'tree_all_evenness', 'tree_all_richness',
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
var_candidates_plotscale_hs = sort(names(sort(vif(model))))


### Collinearity analysis - local plot scale variables from remote sensing ##############################################################
message("Assessing collinearity among variables at local plot scale (remote sensing)")

var_candidates_plotscale_rs = data_homerange_scale[['min']] %>% select(where(is.numeric), -buffer_radius_m) # NOTE: 'min' == 100m radius
var_candidates_plotscale_rs = var_candidates_plotscale_rs %>% select(starts_with("homerange_"))

# A priori reduce list of candidate plot scale variables (remove irrelevant variables)
var_candidates_plotscale_rs = var_candidates_plotscale_rs %>% select(-all_of(c(
  # Age is not a variable of direct interest
  'homerange_age_mean', 'homerange_age_cv'
)))
(names(var_candidates_plotscale_rs))
pairwise_collinearity(var_candidates_plotscale_rs)

# Reduce list of candidate variables (highly correlated, less preferable, etc.)
var_candidates_plotscale_rs = var_candidates_plotscale_rs %>% select(-all_of(c(
  # Canopy cover is more applicable than closure
  'homerange_canopy_closure_mean', 'homerange_canopy_closure_cv',
  # The following measures of structural variation are not of direct interest
  'homerange_ba_cv', 'homerange_treeden_all_cv', 'homerange_treeden_gt4in_dbh_cv', 'homerange_qmd_cv',
  'homerange_canopy_layers_cv', 'homerange_snagden_gt15dbh_cv', 'homerange_downvol_cv', 'homerange_canopy_cover_cv',
  # homerange_ba_mean ~ (homerange_qmd_mean, homerange_htmax_mean, homerange_canopy_cover_mean, homerange_snagden_gt15dbh_mean)
  'homerange_ba_mean', 'homerange_htmax_mean',
  # homerange_htmax_cv ~ homerange_treeden_gt4in_dbh_mean
  'homerange_treeden_gt4in_dbh_mean',
  # Drop canopy cover for comparison with analogous habitat survey (hs) variables
  'homerange_canopy_cover_mean', 'homerange_canopy_layers_mean',
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
var_candidates_plotscale_rs = sort(names(sort(vif(model))))


### Collinearity analysis - homerange scale variables (composition, configuration) ############################################
message("Assessing collinearity among variables at homerange scale (composition, configuration)")

scale = 'max' # e.g. 'min', 'median', 'mean', 'max'
var_candidates_homerangescale = data_homerange_scale[[scale]] %>% select(where(is.numeric), -buffer_radius_m)
var_candidates_homerangescale = var_candidates_homerangescale[,1:22] # only look at composition/configuration variables

# A priori reduce list of candidate plot scale variables (remove irrelevant variables)
var_candidates_homerangescale = var_candidates_homerangescale %>% select(-all_of(c(
  # Drop cover diversity metrics in favor of forest cover diversity metrics
  'cover_diversity', 'cover_richness', 'cover_evenness'
)))
(names(var_candidates_homerangescale))
pairwise_collinearity(var_candidates_homerangescale)

# Pairwise collinearity >= 0.8 across scales 
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
  'focalpatch_core_pcnt', 'aggregation_idx',
  'cover_forest_richness', 'cover_forest_evenness',
  'prop_abund_1', # stand initiation
  'prop_abund_2', # competitive exclusion
  'prop_abund_6'  # roads
)))
(names(var_candidates_homerangescale))
pairwise_collinearity(var_candidates_homerangescale)

# VIF analysis for multicollinearity (consider dropping variable(s) with high VIF values (> 10))
model = lm(rep(1, nrow(var_candidates_homerangescale)) ~ ., data = var_candidates_homerangescale)
sort(vif(model))
model = lm(rep(1, nrow(var_candidates_homerangescale)) ~ . -cover_forest_diversity -density_roads, data = var_candidates_homerangescale)
sort(vif(model))

# RESULT:
# density_roads causes multicollinearity at large spatial scales
var_candidates_homerangescale = sort(names(sort(vif(model))))


### Collinearity analysis - plot scale (habitat survey) plus homerange scale ##############################################################
message("Assessing collinearity among remaining variables at both the plot (habitat survey) and homerange (composition, configuration) scales")

scale = 'max' # e.g. 'min', 'median', 'mean', 'max'
var_candidates_hs_combined = full_join(
  data_plot_scale %>% select(all_of(c("site", var_candidates_plotscale_hs))),
  data_homerange_scale[[scale]] %>% select(all_of(c("site", var_candidates_homerangescale))),
  by = 'site'
) %>% select(-site)
var_candidates_hs_combined = data.frame(lapply(var_candidates_hs_combined, as.numeric))

pairwise_collinearity(var_candidates_hs_combined, threshold = 0.8)

model = lm(rep(1, nrow(var_candidates_hs_combined)) ~ ., data = var_candidates_hs_combined)
sort(vif(model))
model = lm(rep(1, nrow(var_candidates_hs_combined)) ~ . -density_streams, data = var_candidates_hs_combined)
sort(vif(model))

# RESULT: No multicollinearity among variables at min, median, or mean scales. At max scale, density_streams removed.
var_candidates_hs_combined = sort(names(sort(vif(model))))
message("Final list of ", ncol(var_candidates_hs_combined)," uncorrelated variables at both the plot (habitat survey) and homerange (composition, configuration) scales:")
print(var_candidates_hs_combined)

### Collinearity analysis - plot scale (remote sensing) plus homerange scale ##############################################################
message("Assessing collinearity among remaining variables at both the plot (remote sensing) and homerange (composition, configuration) scales")

scale = 'max' # e.g. 'min', 'median', 'mean', 'max'
var_candidates_rs_combined = full_join(
  data_homerange_scale[['min']] %>% select(all_of(c("site", var_candidates_plotscale_rs))), # 100m local plot scale
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


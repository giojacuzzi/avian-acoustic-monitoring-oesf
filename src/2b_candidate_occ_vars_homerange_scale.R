##############################################################################
# Finalize candidate set of variables for occurrence at the home range scale
#
##############################################################################
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
data_plot_scale = data_plot_scale %>% filter(hs == TRUE)

path_data_homerange_scale = "data/cache/occurrence_covariates/data_homerange_scale.rds"
message('Loading homerange scale data from cache ', path_data_homerange_scale)
data_homerange_scale = readRDS(path_data_homerange_scale)

path_rast_cover_clean = "data/cache/occurrence_covariates/rast_cover_clean.tif"
message('Loading raster cover data from cache ', path_rast_cover_clean)
rast_cover_clean = rast(path_rast_cover_clean)

watersheds = st_read('data/environment/GIS Data/watersheds/watersheds.shp') %>% janitor::clean_names() %>% st_transform(crs = crs_m)

mapview(rast_cover_clean,
        alpha.regions = 1.0,
        col.regions = c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', 'darkgray', '#6495ED')) +
  mapview(data_plot_scale) +
  mapview(layer.name = 'plot', st_buffer(data_plot_scale, 100), alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'median', st_buffer(data_plot_scale, 204), alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'mean', st_buffer(data_plot_scale, 468), alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'nearmax', st_buffer(data_plot_scale, 1000), alpha.regions = 0.0, lwd = 2)

# NOTE: For now, just look at a few representative scales

# Start with all candidate variables
data_homerange_scale_focus = as.data.frame(data_homerange_scale[['mean']])
buffer_size = unique(data_homerange_scale_focus$buffer)

data_homerange_scale_candidates = data_homerange_scale_focus %>% st_drop_geometry() %>% select(where(is.numeric), -buffer)
length(names(data_homerange_scale_candidates))
data_homerange_scale_candidates = as.data.frame(data_homerange_scale_candidates)
data_homerange_scale_candidates = as.data.frame(sapply(data_homerange_scale_candidates, as.numeric))
colnames(data_homerange_scale_candidates)

### Principal component analysis ##############################################################

# PCA of habitat survey variables
dpsc_hs_scaled = scale(data_homerange_scale_candidates)
pca = prcomp(dpsc_hs_scaled, center = TRUE, scale. = TRUE)
summary(pca)

# Loadings for first 3 PCs
round(pca$rotation[, 1:3], 2)

factoextra::fviz_eig(pca)
factoextra::fviz_pca_ind(pca, geom.ind = "point", col.ind = "cos2", repel = TRUE)
factoextra::fviz_pca_ind(pca, geom.ind = "point", addEllipses = TRUE, legend.title = "Group", repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", repel = TRUE)


### Collinearity analysis ##############################################################

cor_matrix = cor(data_homerange_scale_candidates, use = "pairwise.complete.obs", method = "pearson")
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= 0.8))

# Reduce list of candidate variables (highly correlated, less preferable, irrelevant, etc.)
vars_to_drop = c(
  # focal_patch_pcnt ~ focal_patch_core_pcnt
  # Drop focal_patch_core_pcnt in favor of uncorrelated metrics of edge density
  'focal_patch_core_pcnt',
  # Drop cover diversity metrics in favor of forest cover diversity metrics
  'cover_diversity', 'cover_richness', 'cover_evenness',
  # cover_forest_diversity ~ cover_forest_richness ~ cover_forest_evenness
  'cover_forest_richness', 'cover_forest_evenness',
  # prop_abund_6 ~ density_roads_paved
  'prop_abund_6',
  # NOTE: prop_abund_1 ~ cw_edge_density
  'prop_abund_1',
  # NOTE: focal_patch_pcnt ~ aggregation_idx
  'aggregation_idx',
  # NOTE: focal_patch_isolation has no variation for species will small home range sizes
  # TODO: calculate and drop prop_abund_7 in favor of density_streams
  
)
data_homerange_scale_candidates = data_homerange_scale_candidates %>% select(-all_of(vars_to_drop))
cor_matrix = cor(data_homerange_scale_candidates, use = "pairwise.complete.obs", method = "pearson")
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= 0.8))
names(data_homerange_scale_candidates)
length(names(data_homerange_scale_candidates))

# VIF analysis for multicollinearity
library(car)
model = lm(rep(1, nrow(data_homerange_scale_candidates)) ~ ., data = data_homerange_scale_candidates)
sort(vif(model))
# Consider dropping variable(s) with high VIF values (> 10)
model = lm(rep(1, nrow(data_homerange_scale_candidates)) ~ . -prop_abund_2, data = data_homerange_scale_candidates)
sort(vif(model))

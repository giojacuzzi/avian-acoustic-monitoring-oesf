##############################################################################
# Finalize candidate set of variables for occurrence at the plot scale
#
##############################################################################
library(mapview)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggrepel)
theme_set(theme_minimal())

crs_m = 32610 # coordinate reference system EPSG:32610, UTM Zone 10N, (meters)

path_data_plot_scale = "data/cache/occurrence_covariates/data_plot_scale.rds"
message('Loading plot scale data from cache ', path_data_plot_scale)
data_plot_scale = readRDS(path_data_plot_scale)
data_plot_scale$stratum = factor(data_plot_scale$stratum)
levels(data_plot_scale$stratum) <- c(
  "stand initiation",
  "stem exclusion",
  "understory reinitiation",
  "old forest",
  "thinning"
)
data_plot_scale$stratum <- factor( # reorder levels so that thinning follows stem exclusion
  data_plot_scale$stratum,
  levels = c(
    "stand initiation",
    "stem exclusion",
    "thinning",
    "understory reinitiation",
    "old forest"
  )
)

path_rast_cover_clean = "data/cache/occurrence_covariates/rast_cover_clean.tif"
message('Loading raster cover data from cache ', path_rast_cover_clean)
rast_cover_clean = rast(path_rast_cover_clean)

watersheds = st_read('data/environment/GIS Data/watersheds/watersheds.shp') %>% janitor::clean_names() %>% st_transform(crs = crs_m)

mapview(rast_cover_clean,
        alpha.regions = 1.0) +
  # mapview(watersheds, col.regions = 'transparent', lwd = 3) +
  mapview(data_plot_scale, zcol = 'hs') +
  mapview(st_buffer(data_plot_scale, 100), col.regions = 'transparent', lwd = 2)

# Sites that were surveyed for habitat data in-person
table(data_plot_scale$hs)

# Age (mean and cv) [#]
hist(data_plot_scale$age_mean, breaks = seq(0, max(data_plot_scale$age_mean) + 10, by = 10))
summary(data_plot_scale$age_mean)
ggplot(data_plot_scale, aes(x = stratum, y = age_mean)) +
  geom_boxplot() + scale_y_continuous(breaks = seq(0, 440, by = 10)) + theme_minimal() +
  ggtitle('Predicted stand age [y]')
hist(data_plot_scale %>% filter(hs == TRUE) %>% pull(age_mean), breaks = seq(0, max(data_plot_scale$age_mean) + 10, by = 10))
summary(data_plot_scale %>% filter(hs == TRUE) %>% pull(age_mean))
ggplot(data_plot_scale %>% filter(hs == TRUE), aes(x = stratum, y = age_mean)) +
  geom_boxplot() + scale_y_continuous(breaks = seq(0, 440, by = 10)) + theme_minimal() +
  ggtitle("Predicted stand age [y]")

# Cover class [categorical]
table(data_plot_scale$stratum)
table(data_plot_scale %>% filter(hs == TRUE) %>% pull(stratum))
table(data_plot_scale$stage)
table(data_plot_scale %>% filter(hs == TRUE) %>% pull(stage))

# Thinning treatment [categorical]
table(data_plot_scale$thinning_treatment, useNA = 'ifany')
table(data_plot_scale %>% filter(hs == TRUE) %>% pull(thinning_treatment), useNA = 'ifany')

### PCA of all variables ##############################################################

# PCA of habitat survey variables
dpsc_hs = data_plot_scale %>% filter(hs == TRUE)  %>% select(where(is.numeric)) %>% st_drop_geometry()
dpsc_hs = dpsc_hs %>% select(-age_cv, -tree_all_evenness, -tree_gte10cm_evenness, -tree_lt10cm_evenness) # drop fields with NA values
dpsc_hs_scaled = scale(dpsc_hs)
pca = prcomp(dpsc_hs_scaled, center = TRUE, scale. = TRUE)
summary(pca)

# Loadings for first 3 PCs
round(pca$rotation[, 1:3], 2)

factoextra::fviz_eig(pca)
factoextra::fviz_pca_ind(pca, geom.ind = "point", col.ind = "cos2", repel = TRUE)
factoextra::fviz_pca_ind(pca, geom.ind = "point",
                         col.ind = st_drop_geometry(data_plot_scale)[data_plot_scale$hs == TRUE, 'stratum'],
                         addEllipses = TRUE, legend.title = "Group", repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", repel = TRUE)

### Comparison of habitat survey and remote sensing data ##############################################################

# Basal area [m2/ha] TODO: confirm if this is a mean value at local level
ggplot(data_plot_scale, aes(x = plot_ba_hs, y = plot_ba_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_ba_hs, data_plot_scale$plot_ba_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_ba_hs, data_plot_scale$plot_ba_rs, na.rm = TRUE)) +
  ggtitle('Basal area [m2/ha]')
cor(data_plot_scale$plot_ba_hs, data_plot_scale$plot_ba_rs, use = "pairwise.complete.obs")

# Tree density (large, dbh > 10 cm) [# trees/ha]
ggplot(data_plot_scale, aes(x = plot_treeden_gt10cmDbh_hs, y = plot_treeden_gt4inDbh_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_treeden_gt10cmDbh_hs, data_plot_scale$plot_treeden_gt4inDbh_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_treeden_gt10cmDbh_hs, data_plot_scale$plot_treeden_gt4inDbh_rs, na.rm = TRUE)) +
  ggtitle('Tree density (large, dbh > 10 cm) [# trees/ha]')
cor(data_plot_scale$plot_treeden_gt10cmDbh_hs, data_plot_scale$plot_treeden_gt4inDbh_rs, use = "pairwise.complete.obs")

# Tree density (small, dbh < 10 cm) [# trees/ha]
ggplot(data_plot_scale, aes(x = plot_treeden_lt10cmDbh_hs, y = plot_treeden_lt4inDbh_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_treeden_lt10cmDbh_hs, data_plot_scale$plot_treeden_lt4inDbh_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_treeden_lt10cmDbh_hs, data_plot_scale$plot_treeden_lt4inDbh_rs, na.rm = TRUE)) +
  ggtitle('Tree density (small, dbh < 10 cm) [# trees/ha]')
cor(data_plot_scale$plot_treeden_lt10cmDbh_hs, data_plot_scale$plot_treeden_lt4inDbh_rs, use = "pairwise.complete.obs")

# Total tree density (all sizes) [# trees/ha]
ggplot(data_plot_scale, aes(x = plot_treeden_all_hs, y = plot_treeden_all_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_treeden_all_hs, data_plot_scale$plot_treeden_all_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_treeden_all_hs, data_plot_scale$plot_treeden_all_rs, na.rm = TRUE)) +
  ggtitle('Tree density (total, all sizes) [# trees/ha]')
cor(data_plot_scale$plot_treeden_all_hs, data_plot_scale$plot_treeden_all_rs, use = "pairwise.complete.obs")

# Tree diameter (mean and SD) [cm]
ggplot(data_plot_scale, aes(x = plot_qmd_all_hs, y = plot_qmd_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_qmd_all_hs, data_plot_scale$plot_qmd_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_qmd_all_hs, data_plot_scale$plot_qmd_rs, na.rm = TRUE)) +
  ggtitle('Tree diameter (quadratic mean) [cm]')
cor(data_plot_scale$plot_qmd_all_hs, data_plot_scale$plot_qmd_rs, use = "pairwise.complete.obs")

# Tree height (mean and cv) [m]
ggplot(data_plot_scale, aes(x = plot_ht_hs, y = plot_htmax_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_ht_hs, data_plot_scale$plot_htmax_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_ht_hs, data_plot_scale$plot_htmax_rs, na.rm = TRUE)) +
  ggtitle('Tree height (mean) [m]')
cor(data_plot_scale$plot_ht_hs, data_plot_scale$plot_htmax_rs, use = "pairwise.complete.obs")
ggplot(data_plot_scale, aes(x = plot_ht_cv_hs, y = plot_htmax_cv_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_ht_cv_hs, data_plot_scale$plot_htmax_cv_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_ht_cv_hs, data_plot_scale$plot_htmax_cv_rs, na.rm = TRUE)) +
  ggtitle('Tree height (cv) [m]')
cor(data_plot_scale$plot_ht_cv_hs, data_plot_scale$plot_htmax_cv_rs, use = "pairwise.complete.obs")

# Density of all snags [# snags/ha]
ggplot(data_plot_scale, aes(x = plot_snagden_hs, y = plot_snagden_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_snagden_hs, data_plot_scale$plot_snagden_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_snagden_hs, data_plot_scale$plot_snagden_rs, na.rm = TRUE)) +
  ggtitle('Snag density [#/ha]')
cor(data_plot_scale$plot_snagden_hs, data_plot_scale$plot_snagden_rs, use = "pairwise.complete.obs")

# Downed wood volume [m3/ha]
ggplot(data_plot_scale, aes(x = plot_downvol_hs, y = plot_downvol_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_downvol_hs, data_plot_scale$plot_downvol_rs, na.rm = TRUE)) +
  ylim(0, max(data_plot_scale$plot_downvol_hs, data_plot_scale$plot_downvol_rs, na.rm = TRUE)) +
  ggtitle('Downed wood volume [m3/ha]')
cor(data_plot_scale$plot_downvol_hs, data_plot_scale$plot_downvol_rs, use = "pairwise.complete.obs")

# Tree density hs versus rs
# Keep HS, and consider keeping split between large and small (instead of "all")
ggplot(data_plot_scale %>% select(site, `stratum`, plot_treeden_lt10cmDbh_hs, plot_treeden_gt10cmDbh_hs, plot_treeden_all_hs) %>%
         tidyr::pivot_longer(cols = starts_with("plot_treeden"), names_to = "variable", values_to = "value"),
       aes(x = stratum, y = value, fill = variable)) +
  geom_boxplot() + coord_cartesian(ylim = c(0, 2500)) + theme_minimal() +
  ggtitle('Tree density (habitat surveys)')
ggplot(data_plot_scale %>% select(site, `stratum`, plot_treeden_lt4inDbh_rs, plot_treeden_gt4inDbh_rs, plot_treeden_all_rs) %>%
         tidyr::pivot_longer(cols = starts_with("plot_treeden"), names_to = "variable", values_to = "value"),
       aes(x = stratum, y = value, fill = variable)) +
  geom_boxplot() + theme_minimal() +
  ggtitle('Tree density (remote sensing)')

# QMD: Only keep "all"
ggplot(data_plot_scale %>% select(site, `stratum`, plot_qmd_lt10cmDbh_hs, plot_qmd_gt10cmDbh_hs, plot_qmd_all_hs) %>%
         tidyr::pivot_longer(cols = starts_with("plot_qmd"), names_to = "metric", values_to = "qmd"),
       aes(x = stratum, y = qmd, fill = metric)) +
  geom_boxplot() + theme_minimal() +
  ggtitle('Tree QMD (habitat surveys)')

### Collinearity analysis ##############################################################

# Start with all candidate variables
data_plot_scale_candidates = data_plot_scale %>% st_drop_geometry() %>% select(where(is.numeric))
length(names(data_plot_scale_candidates))

cor_matrix_plot_scale = cor(data_plot_scale_candidates, use = "pairwise.complete.obs", method = "pearson")
cor_matrix_plot_scale[lower.tri(cor_matrix_plot_scale, diag = TRUE)] = NA
(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix_plot_scale)), !is.na(Freq) & abs(Freq) >= 0.8))

# Drop redundant remote sensing variables in favor of more accurate habitat survey variables
vars_to_drop = c(
  'plot_ba_rs', 'plot_ba_sd_rs', 'plot_treeden_gt4inDbh_rs', 'plot_treeden_all_rs', 'plot_treeden_lt4inDbh_rs', 'plot_sdi_rs', 'plot_qmd_rs', 'plot_htmax_rs', 'plot_htmax_cv_rs', 'plot_snagden_rs', 'plot_downvol_rs'
)
data_plot_scale_candidates = data_plot_scale_candidates %>% select(-all_of(vars_to_drop))
cor_matrix_plot_scale = cor(data_plot_scale_candidates, use = "pairwise.complete.obs", method = "pearson")
cor_matrix_plot_scale[lower.tri(cor_matrix_plot_scale, diag = TRUE)] = NA
(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix_plot_scale)), !is.na(Freq) & abs(Freq) >= 0.8))

# Reduce list of candidate variables (highly correlated, less preferable, irrelevant, etc.)
vars_to_drop = c(
  # Age is not a variable of direct interest. Also, age_point is highly correlated with age_mean (0.99), which is a better representation of stand age. age_cv is unreliable.
  'age_mean', 'age_point', 'age_cv',
  # Basal area ignores tree count; stand structure is better represented by density and size as separate covariates.
  # Josh and Teddy recommend independent variables for 1) tree density and 2) tree size (QMD).
  # Dan, however, notes that integrative metrics that incorporate both tree density and (mean) tree size may do the best job in these analyses.
  'plot_ba_hs',
  # plot_treeden_lt10cmDbh_hs is highly correlated with that of all trees (0.98); favor breakout of tree sizes for interpretability.
  'plot_treeden_all_hs',
  # Quadratic mean diameter is highly correlated with that of all trees (0.99); favor all trees QMD for interpretability.
  'plot_qmd_gt10cmDbh_hs', 'plot_qmd_lt10cmDbh_hs',
  # Canopy cover and closure are highly correlated (0.99). Drop closure in favor of cover here because it is likely a more accurate measurement and is more widely used in landscape-scale analyses as being representative of landscape-scale structure.
  'plot_canopy_closure_rs',
  # Canopy cover and tree height are highly correlated (0.9). Drop height.
  'plot_ht_hs',
  # TODO: Canopy layers and canopy cover are highly correlated (0.85). Drop layers in favor of another vertical heterogeneity metric.
  'canopy_layers_rs',
  # Tree height, height-to-live-crown, and length-of-live-crown are highly correlated
  'plot_hlc_hs', 'plot_llc_hs', 'plot_lcr_hs',
  # TODO: Choose between understory cover and volume. We keep volume here to better represent vertical structure.
  'plot_understory_cover',
  # Individual tree species composition is better summarized with taxonomic diversity metrics
  'tree_gte10cm_density_psme', 'tree_gte10cm_density_thpl', 'tree_gte10cm_density_abam', 'tree_gte10cm_density_tshe', 'tree_gte10cm_density_alru', 'tree_gte10cm_density_pisi',
  'tree_all_density_thpl', 'tree_all_density_abam', 'tree_all_density_tshe', 'tree_all_density_alru', 'tree_all_density_pisi', 'tree_all_density_psme',
  # Most taxonomic diversity metrics are highly correlated between tree size classes; favor all sizes for interpretability.
  'tree_lt10cm_richness', 'tree_gte10cm_richness',
  'tree_lt10cm_evenness', 'tree_gte10cm_evenness',
  'tree_lt10cm_diversity', 'tree_gte10cm_diversity',
  # Furthermore, diversity and evenness are highly correlated; favor shannon diversity for interpretability.
  'tree_all_evenness', 'tree_all_richness',
  # Drop distance from intermittent streams; favor only consistent watercourses indicative of different riparian habitat.
  'dist_watercourses_all'
)
data_plot_scale_candidates = data_plot_scale_candidates %>% select(-all_of(vars_to_drop))
cor_matrix_plot_scale = cor(data_plot_scale_candidates, use = "pairwise.complete.obs", method = "pearson")
cor_matrix_plot_scale[lower.tri(cor_matrix_plot_scale, diag = TRUE)] = NA
(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix_plot_scale)), !is.na(Freq) & abs(Freq) >= 0.8))
names(data_plot_scale_candidates)
length(names(data_plot_scale_candidates))

corrplot::corrplot(cor_matrix_plot_scale, method = "color", tl.cex = 0.8)

# Inspect specific variables by stage / stratum
ggplot(data_plot_scale, aes(x = stage, y = age_mean)) +
  geom_boxplot() + ggtitle("Predicted stand age x structural stage")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_qmd_all_hs)) +
  geom_point() + ggtitle("Tree size (QMD) x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_qmd_all_hs)) +
  geom_boxplot() + ggtitle("Tree size (QMD) x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_treeden_gt10cmDbh_hs)) +
  geom_point() + ggtitle("Tree density (large) x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_treeden_gt10cmDbh_hs)) +
  geom_boxplot() + ggtitle("Tree density (large) x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_treeden_lt10cmDbh_hs)) +
  geom_point() + coord_cartesian(ylim = c(0,2750)) + ggtitle("Tree density (small) x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_treeden_lt10cmDbh_hs)) +
  geom_boxplot() + coord_cartesian(ylim = c(0,2750)) + ggtitle("Tree density (small) x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_ht_hs)) +
  geom_point() + ggtitle("Tree height x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_ht_hs)) +
  geom_boxplot() + ggtitle("Tree height x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_ht_cv_hs)) +
  geom_point() + ggtitle("Tree height coefficient of variation x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_ht_cv_hs)) +
  geom_boxplot() + ggtitle("Tree height coefficient of variation x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_canopy_cover_rs)) +
  geom_point() + ggtitle("Canopy cover x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_canopy_cover_rs)) +
  geom_boxplot() + ggtitle("Canopy cover x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_understory_vol)) +
  geom_point() + ggtitle("Understory volume x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_understory_vol)) +
  geom_boxplot() + ggtitle("Understory volume x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_snagden_hs)) +
  geom_point() + ggtitle("Snag density x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_snagden_hs)) +
  geom_boxplot() + ggtitle("Snag density x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_downvol_hs)) +
  geom_point() + ggtitle("Downed wood vol x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_downvol_hs)) +
  geom_boxplot() + ggtitle("Downed wood vol x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = plot_elev_rs)) +
  geom_point() + ggtitle("Elevation x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = plot_elev_rs)) +
  geom_boxplot() + ggtitle("Elevation x stratum")

ggplot(data_plot_scale, aes(x = age_mean, y = tree_all_diversity)) +
  geom_point() + ggtitle("Tree shannon diversity x stand age")
ggplot(data_plot_scale, aes(x = `stratum`, y = tree_all_diversity)) +
  geom_boxplot() + ggtitle("Tree shannon diversity x stratum")

# psme - Douglas fir
# thpl - Western redcedar
# abam - Pacific silver fir
# tshe - Western hemlock
# alru - Red alder
# pisi - Sitka spruce
#
# Stand initiation: small alru pioneers newly-disturbed habitat as an early-seral species, otherwise plantation psme and tshe are establishing
# Stem exclusion: alru is outcompeted by psme and tshe and mostly disappears, tree size is more uniform
# Understory reinit: tshe dominates, abam establishing, tree size less uniform
# Old forest: tshe is climax species, abam established
metric = "tree_all_density_" # "tree_all_density_", "tree_gte10cm_density_"
df_long = data_plot_scale %>% select(site, age_mean, `stratum`, starts_with(metric)) %>%
  pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
  mutate(species = gsub(metric, "", species))
ggplot(df_long, aes(x = species, y = density, fill = species)) +
  geom_boxplot() +
  facet_wrap(~ stratum, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1250)) +
  labs(title = "Tree species density by stratum")

# VIF analysis for multicollinearity
library(car)
model = lm(rep(1, nrow(data_plot_scale_candidates)) ~ ., data = data_plot_scale_candidates)
sort(vif(model))
# Consider dropping variable(s) with high VIF values
model = lm(rep(1, nrow(data_plot_scale_candidates)) ~ . -plot_canopy_cover_rs, data = data_plot_scale_candidates)
sort(vif(model))

# TODO: Outlier analysis and data sanity checking

### PCA of reduced variable set ##############################################################

data_plot_scale_candidates$site = data_plot_scale$site
data_plot_scale_candidates$stratum = data_plot_scale$stratum
data_plot_scale_candidates$hs = data_plot_scale$hs

# PCA of habitat survey variables
dpsc_hs = data_plot_scale_candidates %>% filter(hs == TRUE)  %>% select(where(is.numeric)) %>% st_drop_geometry()
dpsc_hs_scaled = scale(dpsc_hs)
pca = prcomp(dpsc_hs_scaled, center = TRUE, scale. = TRUE)
summary(pca)

factoextra::fviz_eig(pca)
factoextra::fviz_pca_ind(pca,
                         geom.ind = "point", 
                         col.ind = "cos2",
                         repel = TRUE)
factoextra::fviz_pca_ind(pca,
                         geom.ind = "point",
                         col.ind = data_plot_scale_candidates[data_plot_scale_candidates$hs == TRUE, 'stratum'],
                         addEllipses = TRUE,
                         legend.title = "Group",
                         repel = TRUE)
factoextra::fviz_pca_var(pca, 
                         col.var = "contrib",
                         repel = TRUE)


path_plot_scale_data       = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_homerange_scale_data  = "data/cache/occurrence_covariates/data_homerange_scale.rds"

data_plot_scale = readRDS(path_plot_scale_data)
data_homerange_scale = readRDS(path_homerange_scale_data)

library(tidyr)
library(sf)
library(ggrepel)
theme_set(theme_minimal())

#############################################################################################################
# Compare plot-scale data measured via habitat survey vs remote sensing

data_total = full_join(st_drop_geometry(data_plot_scale), data_homerange_scale[['plot']], by = 'site')

comparison_plot = function(x, y, s, title, x_lab = "Habitat survey", y_lab = "Remote sensing") {
  r = cor(x, y, use = "complete.obs", method = "pearson")
  ggplot() +
    geom_abline(slope = 1, color = "gray") +
    geom_point(aes(x = x, y = y)) +
    geom_text_repel(aes(x = x, y = y, label = s), size = 2) +
    xlim(0, max(x, y, na.rm = TRUE)) +
    ylim(0, max(x, y, na.rm = TRUE)) +
    labs(
      title = title, subtitle = paste0("(Pearson r = ", round(r, 2), ")"),
      x = x_lab, y = y_lab
    )
}

comparison_plot(data_total$plot_ba_hs, data_total$homerange_ba_mean, data_total$site, "Basal area [m2/ha]")

comparison_plot(data_total$plot_treeden_gt10cmDbh_hs, data_total$homerange_treeden_gt4in_dbh_mean, data_total$site, "Density (large trees, DBH > 10 cm) [# trees/ha]")

comparison_plot(data_total$plot_treeden_all_hs, data_total$homerange_treeden_all_mean, data_total$site, "Density (all tree sizes) [# trees/ha]")

comparison_plot(data_total$plot_qmd_all_hs, data_total$homerange_qmd_mean, data_total$site, "Quadratic mean diameter [cm]")

comparison_plot(data_total$plot_ht_hs, data_total$homerange_htmax_mean, data_total$site, "Height mean [cm]")

# TODO: plot_ht_cv_hs appears erroneous?
comparison_plot(data_total$plot_ht_cv_hs, data_total$homerange_htmax_cv, data_total$site, "Height CV [cm]")

comparison_plot(data_total$plot_snagden_hs, data_total$homerange_snagden_gt15dbh_mean, data_total$site, "Snag density [cm]")

comparison_plot(data_total$plot_downvol_hs, data_total$homerange_downvol_mean, data_total$site, "Downed wood volume [m3/ha]")

# Inspect specific variables by stage or stratum
ggplot(data_total, aes(x = stage, y = homerange_canopy_layers_mean)) +
  geom_boxplot()

#############################################################################################################
# Investigate plot-scale variables across strata

# factor `stratum` with labels
# 1 "stand initiation"
# 2 "stem exclusion"
# 3 "understory reinitiation"
# 4 "old forest"
# 5 "commercial thinning"
data_total$stratum = factor(
  data_total$stratum,
  levels = c(
    1, # stand initiation
    2, # stem exclusion
    5, # commercial thinning
    3, # understory reinitiation
    4  # old forest
  ),
  labels = c(
    "stand initiation",
    "stem exclusion",
    "commercial thinning",
    "understory reinitiation",
    "old forest"
  ),
  ordered = TRUE
)

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
metric = "plot_tree_all_density_" # "tree_all_density_", "tree_gte10cm_density_"
dps_long = data_total %>% select(site, `stratum`, starts_with(metric)) %>%
  pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
  mutate(species = gsub(metric, "", species))
ggplot(dps_long, aes(x = species, y = density, fill = species)) +
  geom_boxplot() +
  facet_wrap(~ stratum, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1250)) +  # TODO: confirm outliers
  labs(title = "Tree species density by stage") +
  theme_minimal()

ggplot(data_total, aes(x = stratum, y = homerange_age_mean)) +
  geom_boxplot() +
  labs(title = 'Stand age')

ggplot(data_total, aes(x = stratum, y = plot_ba_hs)) +
  geom_boxplot() +
  labs(title = 'Basal area (all live trees) [m2/ha]')

ggplot(data_total, aes(x = stratum, y = plot_treeden_gt10cmDbh_hs)) +
  geom_boxplot() +
  labs(title = 'Tree density (all live trees large, dbh > 10 cm) [# trees/ha]')
ggplot(data_total, aes(x = stratum, y = plot_treeden_lt10cmDbh_hs)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 3000)) + # TODO: check outliers
  labs(title = 'Tree density (all live trees small, dbh < 10 cm) [# trees/ha]')
ggplot(data_total, aes(x = stratum, y = plot_treeden_all_hs)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 3500)) + # TODO: check outliers
  labs(title = 'Total tree density (all sizes) [# trees/ha]')

ggplot(data_total, aes(x = stratum, y = plot_qmd_gt10cmDbh_hs)) +
  geom_boxplot() +
  labs(title = 'Tree quadratic mean diameter (large, dbh > 10cm) [cm]')
ggplot(data_total, aes(x = stratum, y = plot_qmd_lt10cmDbh_hs)) +
  geom_boxplot() +
  labs(title = 'Tree quadratic mean diameter (small, dbh < 10 cm) [cm]')
ggplot(data_total, aes(x = stratum, y = plot_qmd_all_hs)) +
  geom_boxplot() +
  labs(title = 'Tree quadratic mean diameter (all) [cm]')

ggplot(data_total, aes(x = stratum, y = plot_ht_hs)) +
  geom_boxplot() +
  labs(title = 'Tree height (mean) [m]')
ggplot(data_total, aes(x = stratum, y = plot_ht_cv_hs)) +
  geom_boxplot() +
  labs(title = 'Tree height (cv) [m]')

ggplot(data_total, aes(x = stratum, y = plot_hlc_hs)) +
  geom_boxplot() +
  labs(title = 'Tree height to live crown [m]')
ggplot(data_total, aes(x = stratum, y = plot_llc_hs)) +
  geom_boxplot() +
  labs(title = 'Tree length (depth) of live crown [m]')
ggplot(data_total, aes(x = stratum, y = plot_lcr_hs)) +
  geom_boxplot() +
  labs(title = 'Tree live crown ratio [#]')

ggplot(data_total, aes(x = stratum, y = plot_snagden_hs)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 200)) + # TODO: check outliers
  labs(title = 'Density of all snags [# snags/ha]')

ggplot(data_total, aes(x = stratum, y = plot_downvol_hs)) +
  geom_boxplot() +
  labs(title = 'Downed wood volume [m3/ha]')

ggplot(data_total, aes(x = stratum, y = plot_understory_cover)) +
  geom_boxplot() +
  labs(title = 'Understory vegetation cover [%]')
ggplot(data_total, aes(x = stratum, y = plot_understory_vol)) +
  geom_boxplot() +
  labs(title = 'Understory vegetation volume [m3/ha]')

ggplot(data_total, aes(x = stratum, y = homerange_canopy_cover_mean)) +
  geom_boxplot() +
  labs(title = 'Canopy cover (mean) [%]')
ggplot(data_total, aes(x = age_mean, y = homerange_canopy_cover_mean)) +
  geom_point() + geom_smooth() +
  labs(title = 'Canopy cover (mean) [%]')
ggplot(data_total, aes(x = stratum, y = homerange_canopy_cover_cv)) +
  geom_boxplot() +
  labs(title = 'Canopy cover (cv) [#]')
ggplot(data_total, aes(x = age_mean, y = homerange_canopy_cover_cv)) +
  geom_point() + geom_smooth() +
  labs(title = 'Canopy cover (cv) [#]')

ggplot(data_total, aes(x = stratum, y = homerange_canopy_layers_mean)) +
  geom_boxplot() +
  labs(title = 'Canopy layers (mean) [#]')

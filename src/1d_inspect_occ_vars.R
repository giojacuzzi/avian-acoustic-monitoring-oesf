path_plot_scale_data       = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_homerange_scale_data  = "data/cache/occurrence_covariates/data_homerange_scale.rds"

data_plot_scale = readRDS(path_plot_scale_data)
data_homerange_scale = readRDS(path_homerange_scale_data)

library(tidyr)
library(sf)
library(ggrepel)
library(patchwork)
theme_set(theme_minimal())

#############################################################################################################
# Compare plot-scale data measured via habitat survey vs remote sensing

data_total = full_join(st_drop_geometry(data_plot_scale), data_homerange_scale[['plot']], by = 'site')

comparison_plot = function(x, y, s, title, x_lab = "Habitat survey", y_lab = "Remote sensing") {
  r =   cor(x, y, use = "complete.obs", method = "pearson")
  rho = cor(x, y, use = "complete.obs", method = "spearman")
  ggplot() +
    geom_abline(slope = 1, color = "gray") +
    geom_point(aes(x = x, y = y), shape = 1) +
    geom_text_repel(aes(x = x, y = y, label = s), size = 2) +
    xlim(0, max(x, y, na.rm = TRUE)) +
    ylim(0, max(x, y, na.rm = TRUE)) +
    labs(
      title = title, subtitle = paste0("(r = ", round(r, 2), ", rho = ", round(rho, 2), ")"),
      x = x_lab, y = y_lab
    )
}

p_ba = comparison_plot(data_total$plot_ba_hs, data_total$homerange_ba_mean, data_total$site, "Basal area [m2/ha]")

p_tdl = comparison_plot((data_total$plot_treeden_gt10cmDbh_hs), (data_total$homerange_treeden_gt4in_dbh_mean), data_total$site, "Density (DBH > 10 cm) [#/ha]")

p_tda = comparison_plot(data_total$plot_treeden_all_hs, data_total$homerange_treeden_all_mean, data_total$site, "Density (all) [#/ha]")

p_qmd = comparison_plot(data_total$plot_qmd_all_hs, data_total$homerange_qmd_mean, data_total$site, "Quadratic mean diameter [cm]")

p_h = comparison_plot(data_total$plot_ht_hs, data_total$homerange_htmax_mean, data_total$site, "Height mean [cm]")

# TODO: plot_ht_cv_hs appears erroneous?
data_total$plot_ht_cv_hs = data_total$plot_ht_cv_hs/100
p_hcv = comparison_plot(data_total$plot_ht_cv_hs, data_total$homerange_htmax_cv, data_total$site, "Height CV [cm]")

p_ds = comparison_plot(data_total$plot_snagden_hs, data_total$homerange_snagden_gt15dbh_mean, data_total$site, "Snag density [#/ha]")

p_dwv = comparison_plot(data_total$plot_downvol_hs, data_total$homerange_downvol_mean, data_total$site, "Downed wood volume [m3/ha]")

p_1to1_HsRs = (p_ba  | p_tdl | p_tda) /
              (p_qmd | p_h   | p_hcv) /
              (p_ds  | p_dwv | plot_spacer()) + plot_annotation(title = "Fig 2")
p_1to1_HsRs

# Correlation heatmaps

survey_vars <- data_total %>%
  select(plot_ba_hs, plot_treeden_gt10cmDbh_hs, plot_treeden_all_hs, plot_qmd_all_hs, plot_ht_hs, plot_ht_cv_hs, plot_snagden_hs, plot_downvol_hs) #%>% log1p()

rs_vars <- data_total %>%
  select(homerange_ba_mean, homerange_treeden_gt4in_dbh_mean, homerange_treeden_all_mean, homerange_qmd_mean, homerange_htmax_mean, homerange_htmax_cv, homerange_snagden_gt15dbh_mean, homerange_downvol_mean) #%>% log1p()

method = "spearman"
cor_mat <- cor(survey_vars, rs_vars, method = method, use = "pairwise.complete.obs")

# Step 3: Reshape for ggplot heatmap
cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c("HS", "RS", "correlation")

# Step 4: Add a flag for matched pairs
cor_df <- cor_df %>%
  mutate(Matched = case_when(
    HS == "plot_ba_hs" & RS == "homerange_ba_mean" ~ TRUE,
    HS == "plot_treeden_gt10cmDbh_hs" & RS == "homerange_treeden_gt4in_dbh_mean" ~ TRUE,
    HS == "plot_treeden_all_hs" & RS == "homerange_treeden_all_mean" ~ TRUE,
    HS == "plot_qmd_all_hs" & RS == "homerange_qmd_mean" ~ TRUE,
    HS == "plot_ht_hs" & RS == "homerange_htmax_mean" ~ TRUE,
    HS == "plot_ht_cv_hs" & RS == "homerange_htmax_cv" ~ TRUE,
    HS == "plot_snagden_hs" & RS == "homerange_snagden_gt15dbh_mean" ~ TRUE,
    HS == "plot_downvol_hs" & RS == "homerange_downvol_mean" ~ TRUE,
    TRUE ~ FALSE
  ))

# Step 5: Plot heatmap
ggplot(cor_df, aes(x = HS, y = RS, fill = correlation)) +
  geom_tile(color = "white") +
  geom_tile(data = subset(cor_df, Matched), color = "black", size = 1.2, fill = NA) + # outline matched pairs
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab") +
  geom_text(aes(label = round(correlation, 2)), color = "black", size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste(method, "correlation"),
       y = "Remote sensing",
       x = "Habitat survey") +
  plot_annotation(title = "Fig 4")

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

p_age = ggplot(data_total, aes(x = stratum, y = homerange_age_mean, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Stand age')

p_ba = ggplot(data_total, aes(x = stratum, y = plot_ba_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Basal area [m2/ha]')

p_tdl = ggplot(data_total, aes(x = stratum, y = plot_treeden_gt10cmDbh_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Density (DBH > 10 cm) [# trees/ha]')
p_tds =ggplot(data_total, aes(x = stratum, y = plot_treeden_lt10cmDbh_hs, fill = stratum)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 3000)) + # TODO: check outliers
  labs(subtitle = 'Density (DBH < 10 cm) [# trees/ha]')
p_tda = ggplot(data_total, aes(x = stratum, y = plot_treeden_all_hs, fill = stratum)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 3500)) + # TODO: check outliers
  labs(subtitle = 'Density (all) [# trees/ha]')

p_qmdl = ggplot(data_total, aes(x = stratum, y = plot_qmd_gt10cmDbh_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'QMD (DBH > 10cm) [cm]')
p_qmds = ggplot(data_total, aes(x = stratum, y = plot_qmd_lt10cmDbh_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'QMD (DBH < 10 cm) [cm]')
p_qmda = ggplot(data_total, aes(x = stratum, y = plot_qmd_all_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'QMD (all) [cm]')

p_h = ggplot(data_total, aes(x = stratum, y = plot_ht_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Height (mean) [m]')
p_hcv = ggplot(data_total, aes(x = stratum, y = plot_ht_cv_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Height (cv) [m]')

p_hlc = ggplot(data_total, aes(x = stratum, y = plot_hlc_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Height to live crown [m]')
p_llc = ggplot(data_total, aes(x = stratum, y = plot_llc_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Length live crown [m]')
p_lcr = ggplot(data_total, aes(x = stratum, y = plot_lcr_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Live crown ratio [#]')

p_ds = ggplot(data_total, aes(x = stratum, y = plot_snagden_hs, fill = stratum)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 200)) + # TODO: check outliers
  labs(subtitle = 'Density all snags [# snags/ha]')

p_dwv = ggplot(data_total, aes(x = stratum, y = plot_downvol_hs, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Downed wood vol [m3/ha]')

p_uvc = ggplot(data_total, aes(x = stratum, y = plot_understory_cover, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Understory veg cover [%]')
p_uvv = ggplot(data_total, aes(x = stratum, y = plot_understory_vol, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Understory veg vol [m3/ha]')

p_cc = ggplot(data_total, aes(x = stratum, y = homerange_canopy_cover_mean, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Canopy cover (mean) [%]')
ggplot(data_total, aes(x = age_mean, y = homerange_canopy_cover_mean, fill = stratum)) +
  geom_point() + geom_smooth() +
  labs(subtitle = 'Canopy cover (mean) [%]')
p_cccv = ggplot(data_total, aes(x = stratum, y = homerange_canopy_cover_cv, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Canopy cover (cv) [#]')
ggplot(data_total, aes(x = age_mean, y = homerange_canopy_cover_cv, fill = stratum)) +
  geom_point() + geom_smooth() +
  labs(subtitle = 'Canopy cover (cv) [#]')

p_cl = ggplot(data_total, aes(x = stratum, y = homerange_canopy_layers_mean, fill = stratum)) +
  geom_boxplot() +
  labs(subtitle = 'Canopy layers (mean) [#]')


p_strata_structure = (p_age | p_ba | p_h | p_hcv) /
                     (p_hlc | p_llc | p_lcr | p_cl) /
                     (p_cc | p_cccv | p_tdl | p_tds) /
                     (p_tda | p_qmdl | p_qmds | p_qmda) /
                     (p_ds | p_dwv | p_uvc | p_uvv) + plot_annotation(title = "Fig 1")
p_strata_structure + plot_layout(guides = "collect") &  # collect legends
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

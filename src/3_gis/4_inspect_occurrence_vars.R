## 4_inspect_occurrence_vars.R ###############################################################################
# Inspect variables for occurrence
#
## INPUT:
path_plot_scale_data       = "data/cache/3_gis/3_calc_occurrence_vars/data_plot_scale_2020_clean_strata_4_sites.rds"
path_homerange_scale_data  = "data/cache/3_gis/3_calc_occurrence_vars/data_homerange_scale_2020_clean_strata_4_sites.rds"
###########################################################################################################

source("src/global.R")

data_plot_scale = readRDS(path_plot_scale_data)
data_homerange_scale = readRDS(path_homerange_scale_data)

# DEBUG: Distribution of median homerange scale cover proportions
debug_scale = "pileated woodpecker" # or "median"
hist(data_homerange_scale[[debug_scale]]$pcnt_standinit)
hist(data_homerange_scale[[debug_scale]]$pcnt_compex)
hist(data_homerange_scale[[debug_scale]]$pcnt_mature)
hist(data_homerange_scale[[debug_scale]]$pcnt_thin)

df_long = data_homerange_scale[[debug_scale]] %>%
  select(site, pcnt_standinit, pcnt_compex, pcnt_mature, pcnt_thin, pcnt_road_paved, pcnt_water) %>%
  pivot_longer(cols = starts_with("pcnt_"), names_to = "class", values_to = "percent")
site_order = data_homerange_scale[[debug_scale]] %>% arrange(desc(pcnt_compex)) %>% pull(site)
df_long = df_long %>% mutate(site = factor(site, levels = site_order))
ggplot(df_long, aes(x = percent, y = site, fill = class)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("darkgreen", "brown", "gray", "orange", "purple", "royalblue")) + labs(title = debug_scale)
# DEBUG

# Compare plot-scale data measured via habitat survey vs remote sensing -------------------------------

data_total = full_join(st_drop_geometry(data_plot_scale), data_homerange_scale[['plot']], by = 'site')

comparison_plot = function(d, x_lab, y_lab, y_conv = 1, s, title) {
  x = d[[x_lab]]
  y = d[[y_lab]] * y_conv
  r =   cor(x, y, use = "complete.obs", method = "pearson")
  rho = cor(x, y, use = "complete.obs", method = "spearman")
  ggplot() +
    geom_abline(slope = 1, color = "gray") +
    geom_point(aes(x = x, y = y), shape = 1) +
    # geom_text_repel(aes(x = x, y = y, label = s), size = 2) +
    xlim(0, max(x, y, na.rm = TRUE)) +
    ylim(0, max(x, y, na.rm = TRUE)) +
    labs(
      title = paste(x_lab, "x", y_lab), subtitle = paste0("(r = ", round(r, 2), ", rho = ", round(rho, 2), ")"), x = x_lab, y = y_lab
    )
}

# BA
p_ba_ba      = comparison_plot(data_total, "plot_ba_hs", "ba_mean",      conv_ft2peracre_to_m2perha, data_total$site); p_ba_ba
p_ba_ba_4    = comparison_plot(data_total, "plot_ba_hs", "ba_4_mean",    conv_ft2peracre_to_m2perha, data_total$site); p_ba_ba_4
p_ba_ba_6    = comparison_plot(data_total, "plot_ba_hs", "ba_6_mean",    conv_ft2peracre_to_m2perha, data_total$site); p_ba_ba_6
p_ba_ba_t100 = comparison_plot(data_total, "plot_ba_hs", "ba_t100_mean", conv_ft2peracre_to_m2perha, data_total$site); p_ba_ba_t100
# QMD all
p_qmd_qmd      = comparison_plot(data_total, "plot_qmd_all_hs", "qmd_mean",      conv_in_to_cm, data_total$site); p_qmd_qmd
p_qmd_qmd_6    = comparison_plot(data_total, "plot_qmd_all_hs", "qmd_6_mean",    conv_in_to_cm, data_total$site); p_qmd_qmd_6
p_qmd_qmd_t100 = comparison_plot(data_total, "plot_qmd_all_hs", "qmd_t100_mean", conv_in_to_cm, data_total$site); p_qmd_qmd_t100
# QMD > 10 cm DBH
p_qmd10_qmd     = comparison_plot(data_total, "plot_qmd_gt10cmDbh_hs", "qmd_mean",      conv_in_to_cm, data_total$site); p_qmd10_qmd
p_qmd10_qmd_6   = comparison_plot(data_total, "plot_qmd_gt10cmDbh_hs", "qmd_6_mean",    conv_in_to_cm, data_total$site); p_qmd10_qmd_6
p_qmd10_qmdt100 = comparison_plot(data_total, "plot_qmd_gt10cmDbh_hs", "qmd_t100_mean", conv_in_to_cm, data_total$site); p_qmd10_qmdt100
# Tree density all
p_td_ta = comparison_plot(data_total, "plot_treeden_all_hs", "tree_acre_mean",    conv_peracre_to_perha, data_total$site); p_td_ta
# Tree density > 10 cm DBH
p_td10_ta    = comparison_plot(data_total, "plot_treeden_gt10cmDbh_hs", "tree_acre_mean",    conv_peracre_to_perha, data_total$site); p_td10_ta
p_td10_ta_4  = comparison_plot(data_total, "plot_treeden_gt10cmDbh_hs", "tree_acre_4_mean",  conv_peracre_to_perha, data_total$site); p_td10_ta_4
p_td10_ta_6  = comparison_plot(data_total, "plot_treeden_gt10cmDbh_hs", "tree_acre_6_mean",  conv_peracre_to_perha, data_total$site); p_td10_ta_6
# Height
p_ht_htmax   = comparison_plot(data_total, "plot_ht_hs", "htmax_mean",   conv_ft_to_m, data_total$site); p_ht_htmax
p_ht_ht_t40  = comparison_plot(data_total, "plot_ht_hs", "ht_t40_mean",  conv_ft_to_m, data_total$site); p_ht_ht_t40
p_ht_ht_t100 = comparison_plot(data_total, "plot_ht_hs", "ht_t100_mean", conv_ft_to_m, data_total$site); p_ht_ht_t100
# Snag density
p_snag_snag_20 = comparison_plot(data_total, "plot_snagden_hs", "snag_acre_20_mean", conv_peracre_to_perha, data_total$site); p_snag_snag_20
p_snag_snag_30 = comparison_plot(data_total, "plot_snagden_hs", "snag_acre_30_mean", conv_peracre_to_perha, data_total$site); p_snag_snag_30
# Down wood
p_down_ddwm = comparison_plot(data_total, "plot_downvol_hs", "cfvol_ddwm_mean", conv_ft3peracre_to_m3perha, data_total$site); p_down_ddwm

# p_1to1_HsRs = (p_ba  | p_tdl | p_tda) /
#   (p_qmd | p_h   | p_hcv) /
#   (p_ds  | p_dwv | plot_spacer()) + plot_annotation(title = "Fig 2")
# p_1to1_HsRs

# Correlation heatmaps

survey_vars = data_total %>%
  select(plot_ba_hs, plot_qmd_all_hs, plot_qmd_gt10cmDbh_hs, plot_treeden_all_hs, plot_treeden_gt10cmDbh_hs, plot_ht_hs, plot_snagden_hs, plot_downvol_hs)

rsfris_vars = data_total %>%
  select(ba_mean, ba_4_mean, ba_6_mean, ba_t100_mean, qmd_mean, qmd_6_mean, qmd_t100_mean, qmd_mean, qmd_6_mean, qmd_t100_mean, tree_acre_mean, tree_acre_4_mean, tree_acre_6_mean, htmax_mean, ht_t40_mean, ht_t100_mean, snag_acre_20_mean, snag_acre_30_mean, cfvol_ddwm_mean)

method = "spearman"
cor_mat <- cor(survey_vars, rsfris_vars, method = method, use = "pairwise.complete.obs")

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
  geom_tile(data = subset(cor_df, Matched), color = "black", linewidth = 1.2, fill = NA) + # outline matched pairs
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab") +
  geom_text(aes(label = round(correlation, 2)), color = "black", size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste(method, "correlation"),
       y = "Remote sensing",
       x = "Habitat survey") +
  plot_annotation(title = "Fig 4")

# Best remote sensing proxies:
# BA                 -> ba_mean
# QMD all/>10"       -> qmd_6_mean or qmd_t100_mean
# Tree density > 10" -> tree_acre_4_mean or tree_acre_6_mean
# Height             -> ht_t100_mean

# Investigate plot-scale variables across strata -----------------------------------------------------

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
dps_long = data_total %>% select(site, stratum_4, starts_with(metric)) %>%
  pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
  mutate(species = gsub(metric, "", species))
ggplot(dps_long, aes(x = species, y = density, fill = species)) +
  geom_boxplot() +
  facet_wrap(~ stratum_4, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1250)) +  # TODO: confirm outliers
  labs(title = "Tree species density by stage") +
  theme_minimal()

# p_age = ggplot(data_total, aes(x = stratum_4, y = homerange_age_mean, fill = stratum_4)) +
#   geom_boxplot() +
#   labs(subtitle = 'Stand age')

p_ba = ggplot(data_total, aes(x = stratum_4, y = plot_ba_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Basal area [m2/ha]')

p_tdl = ggplot(data_total, aes(x = stratum_4, y = plot_treeden_gt10cmDbh_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Density (DBH > 10 cm) [# trees/ha]')
p_tds =ggplot(data_total, aes(x = stratum_4, y = plot_treeden_lt10cmDbh_hs, fill = stratum_4)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 3000)) + # TODO: check outliers
  labs(subtitle = 'Density (DBH < 10 cm) [# trees/ha]')
p_tda = ggplot(data_total, aes(x = stratum_4, y = plot_treeden_all_hs, fill = stratum_4)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 3500)) + # TODO: check outliers
  labs(subtitle = 'Density (all) [# trees/ha]')

p_qmdl = ggplot(data_total, aes(x = stratum_4, y = plot_qmd_gt10cmDbh_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'QMD (DBH > 10cm) [cm]')
p_qmds = ggplot(data_total, aes(x = stratum_4, y = plot_qmd_lt10cmDbh_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'QMD (DBH < 10 cm) [cm]')
p_qmda = ggplot(data_total, aes(x = stratum_4, y = plot_qmd_all_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'QMD (all) [cm]')

p_h = ggplot(data_total, aes(x = stratum_4, y = plot_ht_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Height (mean) [m]')
p_hcv = ggplot(data_total, aes(x = stratum_4, y = plot_ht_cv_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Height (cv) [m]')

p_hlc = ggplot(data_total, aes(x = stratum_4, y = plot_hlc_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Height to live crown [m]')
p_llc = ggplot(data_total, aes(x = stratum_4, y = plot_llc_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Length live crown [m]')
p_lcr = ggplot(data_total, aes(x = stratum_4, y = plot_lcr_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Live crown ratio [#]')

p_ds = ggplot(data_total, aes(x = stratum_4, y = plot_snagden_hs, fill = stratum_4)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 200)) + # TODO: check outliers
  labs(subtitle = 'Density all snags [# snags/ha]')

p_dwv = ggplot(data_total, aes(x = stratum_4, y = plot_downvol_hs, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Downed wood vol [m3/ha]')

p_uvc = ggplot(data_total, aes(x = stratum_4, y = plot_understory_cover, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Understory veg cover [%]')
p_uvv = ggplot(data_total, aes(x = stratum_4, y = plot_understory_vol, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Understory veg vol [m3/ha]')

p_cc = ggplot(data_total, aes(x = stratum_4, y = canopy_cover_mean, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Canopy cover (mean) [%]')
ggplot(data_total, aes(x = age_mean, y = canopy_cover_mean, fill = stratum_4)) +
  geom_point() + geom_smooth() +
  labs(subtitle = 'Canopy cover (mean) [%]')
p_cccv = ggplot(data_total, aes(x = stratum_4, y = canopy_cover_cv, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Canopy cover (cv) [#]')
ggplot(data_total, aes(x = age_mean, y = canopy_cover_cv, fill = stratum_4)) +
  geom_point() + geom_smooth() +
  labs(subtitle = 'Canopy cover (cv) [#]')

p_cl = ggplot(data_total, aes(x = stratum_4, y = canopy_layers_mean, fill = stratum_4)) +
  geom_boxplot() +
  labs(subtitle = 'Canopy layers (mean) [#]')


p_strata_structure = (p_ba | p_h | p_hcv) /
  (p_hlc | p_llc | p_lcr | p_cl) /
  (p_cc | p_cccv | p_tdl | p_tds) /
  (p_tda | p_qmdl | p_qmds | p_qmda) /
  (p_ds | p_dwv | p_uvc | p_uvv) + plot_annotation(title = "Fig 1")
p_strata_structure + plot_layout(guides = "collect") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")


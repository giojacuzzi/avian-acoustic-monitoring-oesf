## Calculate landscape matrix variables and test for their effects on richness, FDis, and old growth specialist presence
#
source("src/global.R")

t = 2020 # year: 2020, 2021, 2022, 2023
cover_classification = "clean_strata_4"

path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"

# Get response variables ====================================================

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)
message("Loading site key from ", path_site_key)
site_key = read_csv(path_site_key, show_col_types = FALSE) %>% mutate(site = str_to_lower(site))
message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
stages = model_data$stages

strata = as.factor(site_key$stratum[ match(sites, site_key$site) ])

species_traits = species_traits %>% filter(common_name %in% species)

# Get the expected occurrence for year t
z = msom$sims.list$z
z = z[ , , match(t , c(2020, 2021, 2022, 2023)), ]
z_mean <- apply(z, c(2, 3), mean)
rownames(z_mean) = sites
colnames(z_mean) = species

# Sum across species (dim 3) for each posterior draw and site
# Result: 6000 draws × 224 sites
richness_posterior <- apply(z, c(1, 2), sum)

# Summarize across draws (dim 1)
richness_mean  <- apply(richness_posterior, 2, mean)
richness_sd    <- apply(richness_posterior, 2, sd)
richness_lower <- apply(richness_posterior, 2, quantile, 0.025)
richness_upper <- apply(richness_posterior, 2, quantile, 0.975)

# Combine into a data frame
richness_df <- data.frame(
  site         = sites,
  mean         = richness_mean,
  sd           = richness_sd,
  lower_95     = richness_lower,
  upper_95     = richness_upper
)
stopifnot(richness_df$site == sites)
richness_df$strata = stages$stratum_4
SR = richness_df %>% select(mean, site, strata) %>% rename(SR = mean)

# Get FDis
# Create species-by-traits matrix
trait_matrix = species_traits %>% as.data.frame() %>%
  select(group_nest_ps,
         group_forage_substrate,
         group_diet,
         group_migrant,
         mass)
rownames(trait_matrix) = species_traits$common_name
# TODO: Subset species to those with at least one occurrence?

stopifnot(rownames(trait_matrix) == colnames(z_mean))  # must be TRUE

FDis = fdisp(gowdis(trait_matrix), z_mean)
FDis = stack(FDis$FDis)
colnames(FDis) = c("FDis", "site")
stopifnot(all(FDis$site == sites))
FDis$strata = stages$stratum_4

ggplot(richness_df, aes(x = strata, y = mean, fill = strata)) + geom_boxplot() +
  ggplot(FDis, aes(x = strata, y = FDis, fill = strata)) + geom_boxplot()

# Calculate landscape predictors ===========================================

path_rast_cover = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_", t, "_", cover_classification, ".tif")
path_site_cover_class = paste0("data/cache/3_gis/2_gen_cover_rasters/site_cover_class_", t, "_sf.rds")

sites_sf = read_rds(path_site_cover_class)

message('Loading raster cover data from cache ', path_rast_cover)
rast_cover = rast(path_rast_cover)

mapview(rast_cover) + mapview(sites_sf)

# For each site buffer:
# Fragment size — how much mature forest is nearby?
# Isolation — how far to nearest mature patch?
# Matrix proportion — what % is COMP EXCL?
# Edge density — interface between mature and matrix

library(landscapemetrics)
metrics <- sample_lsm(
  rast_cover,
  y = sites_sf,
  plot_id = sites_sf$site,
  shape = "circle",
  size = 5000,
  what = c(
    "lsm_c_ca",      # total class area (mature forest area in buffer)
    "lsm_c_np",      # number of patches
    "lsm_c_ed",      # edge density
    # "lsm_l_shdi",    # Shannon diversity index of landscape
    "lsm_c_enn_mn"   # mean Euclidean nearest neighbor distance (isolation)
    #     "lsm_p_area",    # area of focal patch (patch the point falls in)
    # "lsm_c_lpi",     # largest patch index (% buffer = largest patch)
    # "lsm_c_clumpy"   # clumpiness - how aggregated vs dispersed is mature forest
  )
)
levels(rast_cover)

# Create class lookup from your raster categories
class_lookup <- data.frame(
  class = levels(rast_cover)[[1]][,1],  # adjust based on levels() output
  stage = levels(rast_cover)[[1]][,2]  # adjust order
)

# Pivot to wide format: one row per site, one column per metric x stage
metrics_wide <- metrics %>%
  left_join(class_lookup, by = "class") %>%
  filter(stage != "other") %>%  # optionally drop non-forest
  pivot_wider(
    id_cols    = plot_id,
    names_from = c(metric, stage),
    values_from = value,
    names_sep  = "_"
  ) %>%
  rename(site = plot_id)

head(metrics_wide)

metrics_filtered <- metrics_wide %>%
  filter(site %in% sites) %>%
  slice(match(sites, site))
stopifnot(all(metrics_filtered$site == sites))

stop("DEBUGGY")

# Inspect
sites_sf_mature = sites_sf %>% filter(stratum_4 == "mature")
mapview(sites_sf_mature %>% left_join(mango %>% select(site, FDis, SR)), zcol = "SR")
mapview(rast_cover) + mapview(sites_sf_mature %>% left_join(mango %>% select(site, FDis, SR)), zcol = "SR")

# Combine data ======================

mango <- metrics_filtered %>%
  left_join(SR, by = c("site")) %>%
  left_join(FDis, by = c("site", "strata"))

# Quick sense check - do mature sites have more mature forest area nearby?
mango %>%
  group_by(strata) %>%
  summarise(
    mean_ca_mature   = mean(ca_mature,   na.rm = TRUE),
    mean_ca_compex   = mean(ca_compex,   na.rm = TRUE),
    mean_ca_standinit = mean(ca_standinit, na.rm = TRUE)
  )

# Replace NAs with 0 for class area and edge density
mango <- mango %>%
  mutate(across(starts_with("ca_"), ~replace_na(., 0)),
         across(starts_with("ed_"), ~replace_na(., 0)))

# Recheck summary - now 0s instead of NAs
mango %>%
  group_by(strata) %>%
  summarise(
    mean_ca_mature    = mean(ca_mature),
    mean_ca_compex    = mean(ca_compex),
    mean_ca_standinit = mean(ca_standinit),
    mean_ca_thin      = mean(ca_thin)
  )

stop("DEBUGGY")

# HYPOTHESIS TEST ============================================

library(lme4)

# Filter to mature sites only
mango_mature <- mango %>%
  filter(strata == "mature")

# Now the question is: within mature sites,
# does surrounding landscape predict FDis?
m_mature <- lm(FDis ~ ca_mature + ca_compex + ca_standinit + ca_thin, data = mango_mature)
summary(m_mature)

m_mature <- lm(SR ~ ca_mature + ca_compex + ca_standinit + ca_thin, data = mango_mature)
summary(m_mature)

# And for indicator species
m_creeper_mature <- lm(brown_creeper ~ ca_mature + ca_compex, 
                       data = mango_indicators %>% filter(strata == "mature"))
summary(m_creeper_mature)







# Does surrounding landscape context predict FDis 
# beyond local stage classification?

# Scale predictors for interpretability
mango_scaled <- mango %>%
  mutate(across(c(ca_mature, ca_compex, ca_standinit, ed_mature), scale))

# Base model - local stage only
m_base <- lmer(FDis ~ strata + (1|landscape_block), data = mango_scaled)

# Add landscape context
m_landscape <- lmer(FDis ~ strata + ca_mature + ca_compex + (1|landscape_block), 
                    data = mango_scaled)

# Compare
AIC(m_base, m_landscape)
summary(m_landscape)

# Key prediction: ca_mature should be positive
# (more mature forest nearby = higher FDis)
# ca_compex should be negative
# (more matrix = lower FDis, homogenization)

# # Does mature forest area in surrounding landscape predict FDis
# # beyond local stage classification?
# m1 <- lmer(FDis ~ strata + mature_area_500m + (1|landscape_block), 
#            data = fd_df)
# 
# # Does isolation predict which old-growth species are present?
# # Use brown creeper occupancy as indicator
# m2 <- lmer(brown_creeper_z ~ strata + dist_to_mature + matrix_prop + 
#              (1|landscape_block), data = site_df)

# MULTI-SCALE NOTES
# Extract metrics at multiple biologically relevant scales
# buffers <- c(100, 250, 500, 1000, 2000, 5000)
# 
# metrics_multiscale <- map(buffers, ~{
#   sample_lsm(lc, y = sites_sf, shape = "circle", size = .x,
#              what = c("lsm_c_ca", "lsm_c_ed", "lsm_l_shdi")) %>%
#     mutate(scale = .x)
# }) %>% bind_rows()
# 
# # Then for each metric, fit models at each scale and compare AIC
# # The scale with lowest AIC is the scale of effect
# scale_selection <- metrics_multiscale %>%
#   group_by(scale, metric) %>%
#   nest() %>%
#   mutate(
#     model = map(data, ~lm(FDis ~ value, data = left_join(.x, fd_df))),
#     AIC   = map_dbl(model, AIC)
#   )
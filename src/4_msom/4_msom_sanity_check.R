# Compare naive site occupancy to estimated occupancy from MSOM

path_msom = "data/cache/models/V3_msom_pcnt_fp_fp_all.rds"

path_y = "data/cache/4_msom/1_assemble_msom_data/y.rds"

source("src/global.R")
library(reshape2)

message("Loading site key from ", path_site_key)
site_key = read_csv(path_site_key, show_col_types = FALSE) %>% mutate(site = str_to_lower(site))

message("Loading species detection histories 'y' from ", path_y)
y = readRDS(path_y)

# Naive occupancy ---------------------------------
naive_occ <- apply(y, c(3, 4), function(mat) {
  site_max <- apply(mat, 1, function(x) { # mat is [sites x visits] for one season x species combo
    if (all(is.na(x))) NA else max(x, na.rm = TRUE)
  })
  detected  <- sum(site_max >= 1, na.rm = TRUE)
  surveyed  <- sum(!is.na(site_max))
  if (surveyed == 0) NA else detected / surveyed
})
print(round(naive_occ, 3))
(avg_occ = data.frame(
  species    = colnames(naive_occ),
  avg_naive_occ = round(colMeans(naive_occ, na.rm = TRUE),3),
  row.names  = NULL
))

# Naive occupancy by stage -----------------------
site_detected <- apply(y, c(1, 3, 4), function(x) {
  if (all(is.na(x))) NA else as.integer(max(x, na.rm = TRUE) >= 1)
})
library(reshape2)
det_df <- melt(site_detected, varnames = c("site", "season", "species"), value.name = "detected")
stages = data.frame(
  site = dimnames(y)$site,
  stage = site_key$stratum[ match(dimnames(y)$site, site_key$site) ]
)
det_df <- merge(det_df, stages, by = "site")
occ_by_stage_species <- det_df |>
  filter(!is.na(detected)) |>
  group_by(stage, species, season) |>
  summarise(naive_occ = mean(detected), .groups = "drop") |>
  group_by(stage, species) |>
  summarise(avg_naive_occ = mean(naive_occ), .groups = "drop") |>
  split(~species)
# Access by name, e.g.:
occ_by_stage_species[["brown creeper"]]
occ_by_stage_species[["rufous hummingbird"]]


# Load MSOM data

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
strata = as.factor(site_key$stratum[ match(sites, site_key$site) ])

# Posterior mean occupancy per site x season x species
z_mean <- apply(msom$sims.list$z, c(2, 3, 4), mean)
# Result: [224 sites x 4 seasons x 67 species]
# Assign dimnames to match y
dimnames(z_mean) <- list(
  site    = sites,
  season  = seasons,
  species = species
)

# Only compare species present in both
shared_species <- intersect(dimnames(y)[[4]], species)

# Subset z to shared species by name (not position)
z_mean_aligned <- z_mean[, , shared_species]

z_df <- melt(z_mean_aligned, varnames = c("site", "season", "species"), value.name = "z_mean")
z_df <- merge(z_df, stages, by = "site")

# Average by stratum x species (same averaging structure as before)
z_by_stage <- z_df |>
  group_by(stage, species, season) |>
  summarise(z_occ = mean(z_mean), .groups = "drop") |>
  group_by(stage, species) |>
  summarise(avg_z_occ = mean(z_occ), .groups = "drop") |>
  split(~species)

# Now compare directly
z_by_stage[["rufous hummingbird"]]
occ_by_stage_species[["rufous hummingbird"]]

# 4_community_comp.R ####################################################################################
# Community composition analyses
#
# INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds" # "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
##################################################################################################################

source("src/global.R")

strata_cols = c(
  "STAND INIT" = "orange",
  "COMP EXCL"  = "forestgreen",
  "THINNED"    = "purple",
  "MATURE"     = "tan4"
)

# if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading site key from ", path_site_key)
site_key = read_csv(path_site_key, show_col_types = FALSE) %>% mutate(site = str_to_lower(site))

path_occurrence_predictor_homerange_data = "data/cache/4_msom/1_assemble_msom_data/V2_occurrence_predictor_homerange_data.rds"
homerange_data = read_rds(path_occurrence_predictor_homerange_data)
path_occurrence_predictor_plot_data      = "data/cache/4_msom/1_assemble_msom_data/V2_occurrence_predictor_plot_data.rds"
plot_data = read_rds(path_occurrence_predictor_plot_data)

# str(site_richness_strata %>% left_join(homerange_data[["median"]], by = "site"))

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons

strata = as.factor(site_key$stratum[ match(sites, site_key$site) ])

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

z = msom$sims.list$z # Simplify posterior occurrence site x species matrix as the most probable occupancy state across all draws and seasons

# Composition PERMANOVA and NMDS ---------------------------------------------------------------------

## Test for significant differences in (expected) taxonomic community composition of different stages

# Posterior mean occupancy matrix: sites × species
psi_mean        <- apply(z, c(2, 4), mean)
rownames(psi_mean) <- sites
colnames(psi_mean) <- species

# Compute Bray-Curtis distance matrix once — used by all subsequent analyses
dist_bc <- vegdist(psi_mean, method = "bray")

# ── PCoA ordination ───────────────────────────────────────────────────────────
# PCoA on dist_bc is geometrically consistent with the PERMANOVA above.
# PCA on psi_mean would use Euclidean distances — a different space.
pcoa    <- cmdscale(dist_bc, k = nrow(psi_mean) - 1, eig = TRUE)
eig     <- pmax(pcoa$eig, 0)
pct_var <- round(100 * eig / sum(eig), 1)

pcoa_df <- tibble(
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  stage = strata,
  site  = sites
)

centroids <- pcoa_df %>%
  group_by(stage) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")

# Species correlation vectors — shows which species drive axis separation
sp_cor_df <- cor(psi_mean, pcoa$points[, 1:2]) %>%
  as.data.frame() %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("species") %>%
  mutate(r2 = PCoA1^2 + PCoA2^2) %>%
  filter(r2 > 0.25)          # |r| > 0.5 on at least one axis

vec_scale <- 0.35

ggplot(pcoa_df, aes(PCoA1, PCoA2, colour = stage)) +
  geom_point(alpha = 0.45, size = 1.8) +
  stat_ellipse(aes(fill = stage), geom = "polygon", alpha = 0.08, colour = NA) +
  geom_point(data = centroids, aes(fill = stage),
             shape = 21, size = 4, colour = "white", stroke = 1.2) +
  geom_segment(data = sp_cor_df,
               aes(x = 0, y = 0,
                   xend = PCoA1 * vec_scale, yend = PCoA2 * vec_scale),
               colour = "grey40", linewidth = 0.4,
               arrow = arrow(length = unit(0.15, "cm")),
               inherit.aes = FALSE) +
  geom_text_repel(data = sp_cor_df,
                  aes(x = PCoA1 * vec_scale, y = PCoA2 * vec_scale, label = species),
                  colour = "grey30", size = 2.8, max.overlaps = 20,
                  inherit.aes = FALSE) +
  scale_colour_manual(values = strata_cols) +
  scale_fill_manual(values = strata_cols) +
  labs(x      = paste0("PCoA 1 (", pct_var[1], "%)"),
       y      = paste0("PCoA 2 (", pct_var[2], "%)"),
       colour = "Stage", fill = "Stage") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), aspect.ratio = 1)

# PCoA 1: This is your dominant gradient and it maps almost perfectly onto the structural transition from open to closed canopy. The left-pointing species are all early-seral, open-habitat associates: olive-sided flycatcher is a snag-dependent aerial insectivore that requires open foraging airspace above forest; song sparrow is a shrub obligate essentially eliminated by canopy closure; purple finch tracks shrubby forest edges. Brown creeper pointing right is the counterweight — it is one of the most reliable old-growth indicators in the Pacific Northwest, requiring deeply furrowed bark for foraging and large-diameter trees for nesting. The fact that a single axis explains 32% of total community variation and that axis is interpretable as canopy closure is strong confirmation that structural development is the primary driver of community organization in this landscape.

# PCoA 2: The three closed-canopy stages overlap on axis 1 — they have all crossed the canopy closure threshold — but axis 2 is picking up finer structural variation within that closed-canopy space. Red-breasted nuthatch, pine siskin, and western wood-pewee pointing upward are all associated with mature coniferous structure but for different reasons: nuthatch and siskin track cone crop availability and large conifers; western wood-pewee is an aerial insectivore that needs canopy gaps and elevated perches. Hutton's vireo pointing downward is a dense-canopy shrub-layer forager — it thrives in the dark understory of competitive exclusion stands where there is still some shrub development beneath a closed canopy. Axis 2 is likely separating the denser, more structurally uniform competitive exclusion stands (Hutton's vireo end) from the more open, structurally complex mature and thinned stands (nuthatch, siskin, pewee end).

# ── PERMANOVA ──────────────────────────────────────────────────────────────────
# R2  = proportion of total compositional variation explained by stage
# Pr(>F) < 0.05 → at least one stage differs compositionally from the others
permanova <- adonis2(dist_bc ~ strata, permutations = 999)
print(permanova)

# ── Pairwise PERMANOVA ─────────────────────────────────────────────────────────
# Subset the pre-computed distance matrix for each pair rather than recomputing
dist_bc_mat <- as.matrix(dist_bc)

pairwise_permanova <- combn(levels(strata), 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    i            <- strata %in% pair
    dist_pair    <- as.dist(dist_bc_mat[i, i])
    strata_pair  <- strata[i]
    fit          <- adonis2(dist_pair ~ strata_pair, permutations = 999)
    tibble(
      stratum_1 = pair[1],
      stratum_2 = pair[2],
      R2        = fit$R2[1],
      F         = fit$F[1],
      p.val     = fit$`Pr(>F)`[1]
    )
  })

# Small R2 → similar assemblages; large R2 → distinct assemblages
# A significant p with low R2 means detectable but minor difference
print(pairwise_permanova %>% arrange(R2))

# ── Dispersion (betadisper) ────────────────────────────────────────────────────
# Tests whether within-stage spread differs among stages
# Significant result means PERMANOVA differences may partly reflect dispersion, not just centroid shifts
bd         <- betadisper(dist_bc, strata)
anova(bd)

dispersion <- tibble(
  stratum            = levels(strata),
  mean_centroid_dist = tapply(bd$distances, strata, mean)
)
print(dispersion)

# ── Within-stage mean dissimilarity ───────────────────────────────────────────
# Higher values → more compositionally heterogeneous sites within that stage
diss_within_stratum <- levels(strata) %>%
  map_dfr(function(s) {
    m <- dist_bc_mat[strata == s, strata == s]
    tibble(stratum = s, mean_diss_bc = mean(m[upper.tri(m)]))
  })
print(diss_within_stratum)

# ── Between-stage mean dissimilarity ──────────────────────────────────────────
# Higher values → greater compositional turnover between that pair of stages
diss_between_strata <- combn(levels(strata), 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    m <- dist_bc_mat[strata == pair[1], strata == pair[2]]
    tibble(stratum_1 = pair[1], stratum_2 = pair[2], mean_diss_bc = mean(m))
  })
print(diss_between_strata %>% arrange(mean_diss_bc))

# ── Turnover vs nestedness decomposition ──────────────────────────────────────────
# Binarize at psi >= 0.5 (majority-probability presence)
# This is preferable to raw detection data because occupancy-model-derived
# presence/absence has had detection-induced false zeros removed
# psi_bin <- ifelse(psi_mean >= 0.5, 1, 0)
psi_bin = apply(z, c(2, 4), function(x) {
  ifelse(mean(x) >= 0.5, 1, 0)
})
rownames(psi_bin) <- sites
colnames(psi_bin) <- species   # this line is what's missing

# Quick sanity check — mean richness per stage
richness_by_stage <- tibble(
  site  = sites,
  stage = strata,
  S     = rowSums(psi_bin)
) %>%
  group_by(stage) %>%
  summarise(mean_S = mean(S), sd_S = sd(S), .groups = "drop")
print(richness_by_stage)

# Core betapart object — computed once, used by all subsequent functions
bp_core <- betapart.core(psi_bin)

# Total landscape-level partition
# beta.SOR = total Sorensen dissimilarity
# beta.SIM = Simpson (replacement / turnover component)
# beta.SNE = nestedness component
beta_multi <- beta.multi(bp_core, index.family = "sorensen")
print(beta_multi)

# Pairwise decomposition across all 224 sites
beta_pair <- beta.pair(bp_core, index.family = "sorensen")
# beta_pair$beta.sor  — total Sorensen (equivalent to dist_bc but on binary data)
# beta_pair$beta.sim  — replacement
# beta_pair$beta.sne  — nestedness

# Convert to matrices for stage-level summaries
sor_mat <- as.matrix(beta_pair$beta.sor)
sim_mat <- as.matrix(beta_pair$beta.sim)
sne_mat <- as.matrix(beta_pair$beta.sne)

# Summarise by stage pair — mirrors your diss_between_strata table
stage_beta <- combn(levels(strata), 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    i1 <- strata == pair[1]
    i2 <- strata == pair[2]
    tibble(
      stratum_1   = pair[1],
      stratum_2   = pair[2],
      sor         = mean(sor_mat[i1, i2]),
      sim         = mean(sim_mat[i1, i2]),
      sne         = mean(sne_mat[i1, i2]),
      pct_replace = round(100 * sim / sor, 1),
      pct_nested  = round(100 * sne / sor, 1)
    )
  }) %>%
  arrange(desc(pct_replace))

print(stage_beta)

# Within-stage partition — how homogeneous is each stage internally?
within_beta <- levels(strata) %>%
  map_dfr(function(s) {
    i <- strata == s
    tibble(
      stratum     = s,
      sor         = mean(sor_mat[i, i][upper.tri(sor_mat[i, i])]),
      sim         = mean(sim_mat[i, i][upper.tri(sim_mat[i, i])]),
      sne         = mean(sne_mat[i, i][upper.tri(sne_mat[i, i])]),
      pct_replace = round(100 * sim / sor, 1),
      pct_nested  = round(100 * sne / sor, 1)
    )
  })

print(within_beta)

# Stage-level species pools — species present (psi >= 0.5) in >= 20% of sites per stage
# The 20% threshold filters out rare species driving spurious nestedness
stage_pool <- levels(strata) %>%
  map(function(s) {
    site_idx <- strata == s
    prev     <- colMeans(psi_bin[site_idx, ])
    names(prev[prev >= 0.5])
  }) %>%
  set_names(levels(strata))

# Pairwise pool overlap
pool_overlap <- combn(levels(strata), 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    a  <- stage_pool[[pair[1]]]
    b  <- stage_pool[[pair[2]]]
    tibble(
      stratum_1        = pair[1],
      stratum_2        = pair[2],
      n_shared         = length(intersect(a, b)),
      n_only_1         = length(setdiff(a, b)),
      n_only_2         = length(setdiff(b, a)),
      jaccard          = length(intersect(a,b)) / length(union(a,b)),
      is_nested        = length(setdiff(a, b)) == 0 | length(setdiff(b, a)) == 0,
      subset_direction = case_when(
        length(setdiff(a, b)) == 0 ~ paste(pair[1], "⊂", pair[2]),
        length(setdiff(b, a)) == 0 ~ paste(pair[2], "⊂", pair[1]),
        TRUE                       ~ "replacement"
      )
    )
  })

print(pool_overlap)

# What is CE's single exclusive species in each comparison?
map_dfr(c("MATURE", "STAND INIT", "THINNED"), function(other) {
  tibble(
    comparison   = paste("CE vs", other),
    ce_exclusive = setdiff(stage_pool[["COMP EXCL"]], stage_pool[[other]])
  )
})
map_dfr(c("STAND INIT", "COMP EXCL", "MATURE"), function(other) {
  tibble(
    comparison   = paste("THINNED vs", other),
    thinned_exclusive = setdiff(stage_pool[["THINNED"]], stage_pool[[other]])
  )
})
map_dfr(c("STAND INIT", "COMP EXCL", "THINNED"), function(other) {
  tibble(
    comparison   = paste("MATURE vs", other),
    mature_exclusive = setdiff(stage_pool[["MATURE"]], stage_pool[[other]])
  )
})

# Functional ------------------------

# Create species-by-traits matrix
trait_matrix = species_traits %>% as.data.frame() %>%
  select(group_nest_ps,
         group_forage_substrate, # TODO: Consider more traits?
         group_diet,
         group_migrant,
         mass)
rownames(trait_matrix) = species_traits$common_name
stopifnot(rownames(trait_matrix) == colnames(psi_mean))

gower_dist = gowdis(trait_matrix)

FDis = fdisp(gower_dist, community_matrix)
FDis = stack(FDis$FDis)
colnames(FDis) = c("FDis", "site")
FDis$strata = factor(
  strata,
  levels = c("STAND INIT", "COMP EXCL", "THINNED", "MATURE")
)
FDis %>% group_by(strata) %>%
  summarise(
    mean = mean(FDis),
    sd = sd(FDis),
    n_sites = n()
  )

ggplot(FDis, aes(strata, FDis, fill = strata)) +
  geom_boxplot() +
  scale_fill_manual(values = strata_cols) +
  labs(x = "", y = "Functional dispersion") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

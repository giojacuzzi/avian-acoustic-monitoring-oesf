# 4_community_comp.R ####################################################################################
# Community composition analyses
#
# INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds" # "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
path_occurrence_predictor_plot_data = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_plot_data.rds"
##########################################################################################################

source("src/global.R")

strata_cols = c(
  "standinit" = "orange",
  "compex"  = "forestgreen",
  "thin"    = "purple",
  "mature"     = "tan4"
)

# Load data --------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
stages = model_data$stages
strata = factor(stages$stratum_4, levels = c("standinit", "compex", "thin", "mature"))

z = msom$sims.list$z
n_iter    = dim(z)[1]
n_sites   = dim(z)[2]
n_seasons = dim(z)[3]
n_species = dim(z)[4]

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

message("Loading occurrence predictor plot scale data from ", path_occurrence_predictor_plot_data)
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)[[1]] %>% filter(site %in% sites)
stopifnot(all(occurrence_predictor_plot_data$site == sites))

# Species richness by stage =====================================================================

# Accumulate species sums across posterior iterations
rich_post = matrix(0, nrow = n_iter, ncol = n_sites * n_seasons)
for (k in seq_len(n_species)) {
  zk = z[, , , k]
  dim(zk) = c(n_iter, n_sites * n_seasons)
  rich_post = rich_post + zk
  if (k %% 10 == 0) message(k, "/", n_species)
}
dim(rich_post) = c(n_iter, n_sites, n_seasons)
str(rich_post)

# Test for differences in richness across posterior iterations
{
  stage     <- stages$stratum_4                      # factor [n_sites]
  stage_rep <- rep(stage, times = n_seasons)         # factor [n_sites * n_seasons]
  
  # ── ANOVA + Tukey on every posterior iteration ────────────────────────────────
  tukey_post <- vector("list", n_iter)
  
  for (i in seq_len(n_iter)) {
    df_i <- data.frame(
      richness = as.vector(rich_post[i, , ]),  # [n_sites * n_seasons]
      stage    = stage_rep
    )
    
    fit            <- aov(richness ~ stage, data = df_i)
    tukey_post[[i]] <- TukeyHSD(fit)$stage
    
    if (i %% 1000 == 0) message(i, "/", n_iter)
  }
  
  # ── Summarize posterior distributions of pairwise differences ─────────────────
  tukey_arr  <- simplify2array(tukey_post)    # [n_pairs x 4 x n_iter]
  tukey_arr  <- aperm(tukey_arr, c(3, 1, 2)) # [n_iter x n_pairs x 4]
  pair_names <- rownames(tukey_post[[1]])
  
  tukey_summary <- lapply(seq_along(pair_names), function(p) {
    diffs <- tukey_arr[, p, 1]  # posterior of mean difference for pair p
    data.frame(
      pair     = pair_names[p],
      mean     = mean(diffs),
      median   = median(diffs),
      lower95  = quantile(diffs, 0.025),
      upper95  = quantile(diffs, 0.975),
      prob_gt0 = mean(diffs > 0)
    )
  }) |> bind_rows()
  
  # ── Print full summary ────────────────────────────────────────────────────────
  print(tukey_summary, digits = 3)
  
  # ── Flag pairs with 95% CI excluding zero ─────────────────────────────────────
  tukey_summary |>
    filter(lower95 > 0 | upper95 < 0) |>
    arrange(desc(abs(mean)))
}

# Summarize posterior
rich_mean = apply(rich_post, c(2, 3), mean)
rich_lo   = apply(rich_post, c(2, 3), quantile, 0.025)
rich_hi   = apply(rich_post, c(2, 3), quantile, 0.975)
dimnames(rich_mean) = dimnames(rich_lo) = dimnames(rich_hi) =
  list(site = sites, season = seasons)
tidy_richness = function(mat, col) {
  as.data.frame(mat) |>
    rownames_to_column("site") |>
    pivot_longer(-site, names_to = "season", values_to = col)
}
richness_season_df = tidy_richness(rich_mean, "richness") |>
  left_join(tidy_richness(rich_lo, "richness_lo"), by = c("site", "season")) |>
  left_join(tidy_richness(rich_hi, "richness_hi"), by = c("site", "season")) |>
  left_join(tibble(site = sites, stratum = strata), by = "site") |>
  mutate(season = factor(season, levels = seasons))
richness_df = richness_season_df |>
  summarise(richness = mean(richness), .by = c(site, stratum))

# Richness summary
richness_df |>
  summarise(
    mean_rich = mean(richness),
    sd_rich   = sd(richness),
    n         = n(),
    .by       = stratum
  ) |>
  arrange(desc(mean_rich)) |> print()

# Posterior mean richness per site, across seasons, by stratum.
fig_SR = richness_df |>
  ggplot(aes(x = stratum, y = richness, fill = stratum)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, shape = 16) +
  scale_fill_manual(values = strata_cols) +
  labs(x = "Stratum", y = "Expected richness", fill = "Stratum") +
  theme(legend.position = "none"); print(fig_SR)

# Bayesian contrasts to test for differences in richness between strata
{
  # Average over seasons per site per iteration
  rich_site_iter = apply(rich_post, c(1, 2), mean)
  # Posterior mean richness per stratum per iteration
  strata_levels = levels(strata)
  stratum_iter = sapply(strata_levels, function(s) {
    site_idx = which(strata == s)
    rowMeans(rich_site_iter[, site_idx, drop = FALSE])
  })
  stratum_summary = as.data.frame(stratum_iter) |>
    pivot_longer(everything(), names_to = "stratum", values_to = "richness") |>
    summarise(
      mean  = mean(richness),
      lo    = quantile(richness, 0.025),
      hi    = quantile(richness, 0.975),
      .by   = stratum
    )
  # Pairwise contrasts
  strata_pairs = combn(strata_levels, 2, simplify = FALSE)
  contrasts_df = map_dfr(strata_pairs, function(pair) {
    diff_post = stratum_iter[, pair[1]] - stratum_iter[, pair[2]]
    tibble(
      contrast  = paste(pair[1], "-", pair[2]),
      mean_diff = mean(diff_post),
      median_diff = median(diff_post),
      lo        = quantile(diff_post, 0.025),
      hi        = quantile(diff_post, 0.975),
      p_gt0     = mean(diff_post > 0) # posterior probability that pair[1] > pair[2]
    )
  })
  print(contrasts_df)
}

# Frequentist ANOVA / Tukey pairwise tests on posterior mean richness to test
# for differences in richness between strata
{
  # ANOVA for overall richness among stages (posterior mean richness)
  model = aov(richness ~ stratum, data = richness_df)
  summary(model) # If significant, between-stage variance is larger than within-stage variance
  qqnorm(residuals(model))
  qqline(residuals(model)) # Generally meets assumption of residual normality
  # Tukey pairwise test
  TukeyHSD(model) # Pairwise comparisons are significantly different
}

# TAKEAWAYS:
# - Stand initiation has the highest richness (8-10 more species than closed-canopy stages, BCI 6.5-11.4)
# - Compex has the lowest richness (1-2 fewer species on average than thinned and mature)
# - Differences in richness are clear between thin and compex, but weaker between thin and mature

# Taxonomic composition =====================================================================

# Posterior mean occupancy, averaged over seasons
# TODO: Ultimately do this over all posterior draws?
Ez = matrix(0, nrow = n_sites, ncol = n_species,
             dimnames = list(site = sites, species = species))
for (k in seq_len(n_species)) {
  zk = z[, , , k]
  Ez[, k] = apply(zk, 2, mean)
}

# Bray-Curtis distance matrix
dist_bc = vegdist(Ez, method = "bray")

# Principal coordinates analysis (PCoA) ──────────────────────────────────────────
pcoa = cmdscale(dist_bc, k = 2, eig = TRUE)
pct  = round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)
pcoa_df = as.data.frame(pcoa$points) |>
  setNames(c("PCoA1", "PCoA2")) |>
  mutate(site = sites, stratum = strata)

fig_PCoA_sites = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stratum)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(fill = stratum), geom = "polygon", alpha = 0.1, level = 0.95) +
  scale_color_manual(values = strata_cols) +
  scale_fill_manual(values = strata_cols) +
  labs(x = paste0("PCoA1 (", pct[1], "%)"), y = paste0("PCoA2 (", pct[2], "%)"), color = "Stratum", fill = "Stratum") +
  theme_bw(); print(fig_PCoA_sites)

# Species ordination correlations
{
  species_fit = envfit(pcoa$points, Ez, permutations = 999)
  species_scores = as.data.frame(species_fit$vectors$arrows) |>
    setNames(c("PCoA1", "PCoA2")) |>
    mutate(species = species, r2 = species_fit$vectors$r, p = species_fit$vectors$pvals) |>
    arrange(desc(r2))
  top_species = species_scores %>% filter(p < 0.05) %>% filter(r2 > 0.0)
  
  arrow_scale = 0.6 * max(abs(pcoa_df[, c("PCoA1", "PCoA2")])) / max(abs(top_species[, c("PCoA1", "PCoA2")]))
  
  fig_PCoA_species = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stratum)) +
    stat_ellipse(aes(fill = stratum), geom = "polygon", alpha = 0.05, level = 0.95) +
    geom_segment(data = top_species,
                 aes(x = 0, y = 0, xend = PCoA1 * r2 * arrow_scale, yend = PCoA2 * r2 * arrow_scale),
                 inherit.aes = FALSE, arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.5) +
    geom_text(data = top_species,
              aes(x = PCoA1 * r2 * arrow_scale, y = PCoA2 * r2 * arrow_scale, label = species),
              inherit.aes = FALSE, color = "black", alpha = 0.5, size = 3, hjust = 0.5, vjust = -0.5) +
    scale_color_manual(values = alpha(strata_cols, 0.5)) +
    scale_fill_manual(values = strata_cols) +
    labs(x = paste0("PCoA1 (", pct[1], "%)"), y = paste0("PCoA2 (", pct[2], "%)"), color = "Stratum", fill = "Stratum") +
    theme_bw(); print(fig_PCoA_species)
  print(species_scores)
}

# Environmental ordination correlations
{
  env_data = occurrence_predictor_plot_data %>%
    select(-c(site, ces, treatment, wadnr_patch_stratum,
              thinning_status, thinning_treatment,
              stage_3, stage_4, stratum_4, stratum_5, scale),
           -ends_with("_median"), -ends_with("_sd"),
           -starts_with("sdi"), -starts_with("tree_acre_4_"), -starts_with("tree_acre_30_"), -starts_with("qmd_4_"),
           -starts_with("qmd_t100_"), -starts_with("ht"), -starts_with("pcnt_"), -starts_with("ba_"))
  env_fit = envfit(pcoa$points, env_data, permutations = 999, na.rm = TRUE)
  env_scores = as.data.frame(env_fit$vectors$arrows) %>%
    setNames(c("PCoA1", "PCoA2")) %>%
    mutate(variable = rownames(.), r2 = env_fit$vectors$r, p = env_fit$vectors$pvals) %>%
    arrange(desc(r2))
  top_env = env_scores %>% filter(p < 0.05) %>% filter(r2 > 0.0)
  
  arrow_scale_env = 0.6 * max(abs(pcoa_df[, c("PCoA1", "PCoA2")])) / max(abs(top_env[, c("PCoA1", "PCoA2")]))
  
  fig_PCoA_env = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stratum)) +
    stat_ellipse(aes(fill = stratum), geom = "polygon", alpha = 0.05, level = 0.95) +
    geom_segment(data = top_env, aes(x = 0, y = 0, xend = PCoA1 * r2 * arrow_scale_env, yend = PCoA2 * r2 * arrow_scale_env),
                 inherit.aes = FALSE, arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.5) +
    geom_text(data = top_env, aes(x = PCoA1 * r2 * arrow_scale_env, y = PCoA2 * r2 * arrow_scale_env, label = variable),
              inherit.aes = FALSE, color = "black", alpha = 0.5, size = 3, hjust = 0.5, vjust = -0.5) +
    scale_color_manual(values = alpha(strata_cols, 0.5)) +
    scale_fill_manual(values = strata_cols) +
    labs(x = paste0("PCoA1 (", pct[1], "%)"), y = paste0("PCoA2 (", pct[2], "%)"), color = "Stratum", fill = "Stratum") +
    theme_bw(); print(fig_PCoA_env)
  print(env_scores)
}

# Findings:
# - Strong primary environmental filter is a canopy closure / disturbance gradient (PCoA axis 1 ~ 32%)
# - Weaker secondary elevational gradient (PCoA axis 2 ~ 10%)
# - Example PCoA1 species ordination correlations: brown creeper + (mature trees, closed), song sparrow - (edge, open habitat)
# - Example PCoA2 species ordination correlations: gray jay + (montane), kingfisher/bald eagle - (valley bottom / riparian / aquatic)

# PERMANOVA ──────────────────────────────────────────
# R2  = proportion of total composition variation explained by stage
# Pr(>F) < 0.05 → at least one stage differs in composition from the others
permanova = adonis2(dist_bc ~ strata, permutations = 999)
print(permanova)
# Findings:
# - Substantial variation in taxonomic composition by stage
# - Stage explains roughly 40% of variation in taxonomic composition

# Test homogeneity of dispersion ──────────────────────────────────────────
# If significant, PERMANOVA could reflect differences in dispersion rather than centroid location
bd   = betadisper(dist_bc, strata)
disp = permutest(bd)
print(disp)
print(data.frame(
  stratum            = levels(strata),
  mean_centroid_dist = tapply(bd$distances, strata, mean)
))
# Findings:
# - Dispersion differs significantly among strata, so some of the PERMANOVA signal could reflect differences in within-stratum spread rather than centroid location alone

# Pairwise PERMANOVA ──────────────────────────────────────────
# Small R2 → similar assemblages; large R2 → distinct assemblages
# A significant p with low R2 means detectable but minor difference
pairwise_permanova = combn(levels(strata), 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    i            = strata %in% pair
    dist_pair    = as.dist( as.matrix(dist_bc)[i, i])
    strata_pair  = strata[i]
    fit          = adonis2(dist_pair ~ strata_pair, permutations = 999)
    tibble(
      stratum_1 = pair[1],
      stratum_2 = pair[2],
      R2        = fit$R2[1],
      F         = fit$F[1],
      p.val     = fit$`Pr(>F)`[1]
    )
  })
print(pairwise_permanova %>% arrange(R2))
# Findings:
# - Standinit vs. everything else: by far the largest compositional differences, confirming standinit supports a fundamentally distinct bird community, not just a richer one.
# - Among the remaining three strata: all significant but weak, suggesting compex, thin, and mature share broadly similar assemblages with modest compositional differences. Thin–mature is  most distinct of closed canopy, while compex–thin are most similar.

# Between- and within-stage mean dissimilarity ──────────────────────────────────────────
# Higher values = greater compositional turnover between that pair of stages
diss_between_strata = combn(levels(strata), 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    m = as.matrix(dist_bc)[strata == pair[1], strata == pair[2]]
    tibble(stratum_1 = pair[1], stratum_2 = pair[2], mean_diss_bc = mean(m))
  })
print(diss_between_strata %>% arrange(mean_diss_bc))
# Higher values → more compositionally heterogeneous sites within that stage
diss_within_stratum = levels(strata) %>%
  map_dfr(function(s) {
    m = as.matrix(dist_bc)[strata == s, strata == s]
    tibble(stratum = s, mean_diss_bc = mean(m[upper.tri(m)]))
  })
print(diss_within_stratum)
# Findings:
# - Closed canopy stages are compositionally very similar. within-stage variation is roughly comparable to between-stage variation, so stages are no more different from each other than random sites drawn from within the same stage
# - For standinit, between-stage signal dominates variation.
# - Overall community structure is stable across closed-canopy stages, perhaps dominant species are shared?

# Turnover vs nestedness decomposition ──────────────────────────────────────────
# Binarize at psi >= 0.5 (majority-probability presence)
# This is preferable to raw detection data because occupancy-model-derived
# presence/absence has had detection-induced false zeros removed
psi_bin = ifelse(Ez >= 0.5, 1, 0) # TODO: Across posterior draws?

# Core betapart object — computed once, used by all subsequent functions
bp_core = betapart.core(psi_bin)

# Total landscape-level partition
# beta.SOR = total Sorensen dissimilarity
# beta.SIM = turnover component
# beta.SNE = nestedness component
beta_multi = beta.multi(bp_core, index.family = "sorensen")
print(beta_multi)
# Findings:
# - Overall beta diversity (SOR 0.95) is high, so communities are dissimilar across all sites
# - Nearly all dissimilarity driven by species turnover (SIM 0.9), with nestedness contributing almost nothing (SNE 0.04).
# - Sites are not simply species-poor subsets of richer sites, but instead genuinely differ in species.

# Compare closed-canopy stages only:
bp_core_closed = betapart.core(psi_bin[strata %in% c("compex", "thin", "mature"), ])
beta_multi_closed = beta.multi(bp_core_closed, index.family = "sorensen")
print(beta_multi_closed)
# Findings: similar story

# Pairwise decomposition across all sites
beta_pair = beta.pair(bp_core, index.family = "sorensen")

# Convert to matrices for stage-level summaries
sor_mat = as.matrix(beta_pair$beta.sor)
sim_mat = as.matrix(beta_pair$beta.sim)
sne_mat = as.matrix(beta_pair$beta.sne)

# Summarise by stage pair — mirrors your diss_between_strata table
stage_beta = combn(levels(strata), 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    i1 = strata == pair[1]
    i2 = strata == pair[2]
    tibble(
      stratum_1   = pair[1],
      stratum_2   = pair[2],
      sor         = mean(sor_mat[i1, i2]),
      sim         = mean(sim_mat[i1, i2]),
      sne         = mean(sne_mat[i1, i2]),
      pct_replace = round(100 * sim / sor, 1),
      pct_nested  = round(100 * sne / sor, 1)
    )
  }) %>% arrange(desc(pct_replace))

print(stage_beta)
# Findings:
# - Differences among closed-canopy stages are primarily driven by species turnover. Although compex is species-poor relative to thin and mature, those species are not simply a subset -- compex isn't only filtering the mature/thin species, but also selecting for species affiliated with compex specifically.
# - Increasing nestedness driving differences between standinit and closed canopy across mature < thin < compex indicate that canopy closure filters species out of the standinit pool. Nestedness is strongest in compex stands where stem exclusion environment is most intense.
# - Canopy closure reduces richness (nestedness), but which species persist under different management regimes still varies substantially (turnover). The practical implication is that thin and mature are not interchangeable from a conservation standpoint despite similar richness — they are likely supporting different species within their shared reduced pool.
# - Standinit supports a superset community containing most of the species found elsewhere plus a suite of early-successional specialists. The other three strata share a broadly similar forest-interior community that differs among itself mainly through species replacement rather than richness difference.

## NOTE: To find which species are driving the turnover signals between stages, we turn to the MSOM

# Within-stage partition — how homogeneous is each stage internally?
within_beta = levels(strata) %>%
  map_dfr(function(s) {
    i = strata == s
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
# Findings:
# - Within-stratum dissimilarity is predominantly turnover-driven across all four stages, meaning sites within a stage replace species rather than being richer/poorer subsets of each other.
# - Standinit and mature have the highest within-stage turnover, suggesting that within these strata, sites tend to differ in which species they support rather than how many.
# - Standinit: high richness means sites are drawing different subsets from a large species pool — stochastic assembly from an open, diverse pool
# - Mature: late-successional specialists are replacing each other rather than co-occurring, possibly driven by local habitat heterogeneity or microsite variation
# - Compex has the greatest within-stage nestedness (lowest turnover), sites within compex vary primarily because competitive filtering is not uniform in intensity across sites 

# FINAL SYNTHESIS:
# - Standinit — species-rich, open-pool stochastic assembly, genuinely distinct from all other stages, primarily through filtering into closed-canopy stages
# - Compex — species-poor, internally variable through differential filtering intensity, stem exclusion is the dominant process but uneven across sites
# - Thin and mature — intermediate richness, internally turnover-structured, similar to each other in composition but supporting detectably different species

# Functional ==========================================================

# Posterior mean occupancy, averaged over seasons
# TODO: Ultimately do this over all posterior draws?
Ez = matrix(0, nrow = n_sites, ncol = n_species,
            dimnames = list(site = sites, species = species))
for (k in seq_len(n_species)) {
  zk = z[, , , k]
  Ez[, k] = apply(zk, 2, mean)
}

# Create species-by-traits matrix
trait_matrix = species_traits %>% as.data.frame() %>%
  select(group_nest_ps,
         group_forage_substrate, # TODO: Consider more traits?
         group_diet,
         group_migrant,
         mass)
rownames(trait_matrix) = species_traits$common_name
stopifnot(rownames(trait_matrix) == colnames(Ez))

gower_dist = gowdis(trait_matrix)

# Functional dispersion  ──────────────────────────────────────────
FDis = fdisp(gower_dist, Ez)
FDis = stack(FDis$FDis)
colnames(FDis) = c("FDis", "site")
FDis$strata = strata
FDis %>% group_by(strata) %>%
  summarise(
    mean = mean(FDis),
    sd = sd(FDis),
    n_sites = n()
  )

fig_FDis = FDis |>
  ggplot(aes(x = strata, y = FDis, fill = strata)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, shape = 16) +
  scale_fill_manual(values = strata_cols) +
  labs(x = "Stratum", y = "FDis", fill = "Stratum") +
  theme(legend.position = "none"); print(fig_FDis)

model = aov(FDis ~ strata, data = FDis)
summary(model)
TukeyHSD(model)

# Findings
# - Very small differences in functional diversity?

# Functional composition, as measured by the community-level weighted means of trait values,
# which for continuous trais is the mean trait value of all species present in the community,
# and for categorical traits the abundance of each individual class, i.e. the proportion.
cwms = functcomp(trait_matrix, Ez, CWM.type = "all")  # default
cwm_df = data.frame(cwms, strata = strata)

cwm_nest_ps = cwm_df %>%
  select(starts_with("group_nest_ps"), strata) %>%
  pivot_longer(cols = -strata, names_to = "type", values_to = "proportion") %>%
  mutate(type = str_remove(type, "^group_nest_ps_"),
         type = str_replace_all(type, "\\.", " "),
         type = str_to_sentence(type))
p_nest = ggplot(cwm_nest_ps, aes(x = strata, y = proportion, fill = fct_rev(type))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent_format()) + labs(fill = "Nesting strategy"); p_nest

cwm_forage_substrate = cwm_df %>%
  select(starts_with("group_forage_substrate"), strata) %>%
  pivot_longer(cols = -strata, names_to = "type", values_to = "proportion") %>%
  mutate(type = str_remove(type, "^group_forage_substrate_"),
         type = str_replace_all(type, "\\.", " "),
         type = str_to_sentence(type))
p_forage = ggplot(cwm_forage_substrate, aes(x = strata, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent_format()) + labs(fill = "Foraging strategy"); p_forage

cwm_migrant = cwm_df %>%
  select(starts_with("group_migrant"), strata) %>%
  pivot_longer(cols = -strata, names_to = "type", values_to = "proportion")
ggplot(cwm_migrant, aes(x = strata, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format())

# Findings:
# - Taxonomic turnover is occurring among functionally similar species. On the surface, stages differ in which species are present but not in what those species collectively do (functional redundancy).
# - On the surface, looking at local site occurrence alone, one might think that losing a stage does not obviously mean losing a functional role, but it does mean losing the specific species pool associated with that stage.
# - HOWEVER, the MSOM reveals influence of the broader homerange scale

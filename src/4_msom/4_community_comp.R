# 4_community_comp.R ####################################################################################
# Community composition analyses
#
# OUTPUT: 
out_cache_dir  = "data/cache/4_msom/4_community_comp"
#
# INPUT:
path_msom = "data/cache/models/V3_msom_pcnt_fp_fp_all.rds" # "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
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

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)

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

species_traits = species_traits %>% filter(common_name %in% species)

# Taxonomic and functional diversity -----------------------------------------------------------------

z = msom$sims.list$z # Simplify posterior occurrence site x species matrix as the most probable occupancy state across all draws and seasons
z_binary = apply(z, c(2, 4), function(x) {
  ifelse(mean(x) >= 0.5, 1, 0)
})
dim(z_binary)
rownames(z_binary) = sites
colnames(z_binary) = species

## Taxonomic diversity (species richness)

richness_per_site = rowSums(z_binary)
site_richness_strata = data.frame(
  site = names(richness_per_site),       # site IDs
  richness = as.numeric(richness_per_site),  # richness values
  strata = strata                        # stratum factor/vector
)
site_richness_strata$strata = factor(
  site_richness_strata$strata, levels = c("STAND INIT", "COMP EXCL", "THINNED", "MATURE")
)

ggplot(site_richness_strata %>% left_join(homerange_data[["median"]], by = "site") %>% left_join(plot_data, by = "site"), aes(x = age_mean, y = richness)) + geom_point() + geom_smooth()
ggplot(site_richness_strata %>% left_join(homerange_data[["median"]], by = "site") %>% left_join(plot_data, by = "site") %>% filter(age_mean < 25), aes(x = age_mean, y = richness)) + geom_point() + geom_smooth()

# ANOVA for overall richness among stages
model = aov(richness ~ strata, data = site_richness_strata)
summary(model) # If significant, between-stage variance is larger than within-stage variance
qqnorm(residuals(model))
qqline(residuals(model)) # Generally meets assumption of residual normality
# Tukey pairwise test
TukeyHSD(model) # Pairwise comparisons are significantly different

site_richness_strata %>% group_by(strata) %>%
  summarize(
    mean_richness = mean(richness),
    sd_richness = sd(richness),
    max_richness = max(richness),
    min_richness = min(richness),
    n_sites = n()
  )

p_SR = ggplot(site_richness_strata, aes(x = strata, y = richness, fill = strata)) +
  geom_boxplot() +
  scale_fill_manual(values = strata_cols) +
  labs(subtitle = "Taxonomic richness") +
  theme(legend.position = "bottom"); print(p_SR)

## Functional diversity (functional dispersion)

# Select traits
trait_matrix = species_traits %>% as.data.frame() %>% select(mass, group_migrant, group_nest_ps, group_diet)
rownames(trait_matrix) = species_traits$common_name

# Subset species to those with at least one occurrence
species_totals = colSums(z_binary)
present_species = names(species_totals[species_totals > 0])
community_matrix = z_binary[, present_species]
trait_matrix = trait_matrix[present_species, ] %>% mutate(across(where(is.character), as.factor))

# Compute functional diversity using FD
# "FDis is the mean distance in multidimensional trait space of individual species to the centroid of all species...
# [it] is the mulutivariate analogue of the weighted mean absolute deviation; this makes [it] unaffected by species
# richness by construction... 
# Functional dispersion (FDis; Laliberté and Legendre 2010) is computed from the uncorrected species-species
# distance matrix via fdisp. Axes with negatives eigenvalues are corrected following the approach of Anderson (2006).
# When all species have equal abundances (i.e. presence-absence data), FDis is simply the average distance to the
# centroid (i.e. multivariate dispersion) as originally described by Anderson (2006)...
# For unweighted presence-absence data, FDis can be used for a formal statistical test of differences in FD."
FDis = fdisp(gowdis(trait_matrix), community_matrix)
FDis = stack(FDis$FDis)
colnames(FDis) = c("FDis", "site")
FDis$strata = strata
FDis %>% group_by(strata) %>%
  summarise(
    mean = mean(FDis),
    sd = sd(FDis),
    n_sites = n()
  )

fd_results = dbFD(
  x = trait_matrix, a = community_matrix, corr = "sqrt" # "PCoA axes corresponding to negative eigenvalues are imaginary axes that cannot be represented in a Euclidean space, but simply ignoring these axes would lead to biased estimations of FD. Hence in dbFD one of four correction methods are used."
)
fd_df = data.frame(
  site = rownames(community_matrix),
  strata = strata,
  FRic = fd_results$FRic,
  FEve = fd_results$FEve, # 
  FDis = fd_results$FDis, # Mean distance of species to community centroid
  FDiv = fd_results$FDiv  # 
)
fd_df %>% group_by(strata) %>% summarise(mean = mean(FDis), sd = sd(FDis), n_sites = n())

fd_long = fd_df %>%
  pivot_longer(cols = c(FRic, FDis, FEve, FDiv), names_to = "metric", values_to = "value") %>%
  mutate(strata = factor(strata, levels = c("STAND INIT", "COMP EXCL", "THINNED", "MATURE")))

p_FD = ggplot(fd_long, aes(strata, value, fill = strata)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = strata_cols) +
  theme(legend.position = "none"); print(p_FD)

ggplot(site_richness_strata %>% left_join(FDis, by = "site") %>% left_join(plot_data, by = "site"), aes(x = age_mean, y = FDis)) + geom_point() + geom_smooth()

# ANOVA for overall functional dispersion among stages
model = aov(FDis ~ strata, data = FDis)
summary(model) # If significant, between-stage variance is larger than within-stage variance
qqnorm(residuals(model))
qqline(residuals(model)) # Generally meets assumption of residual normality
kruskal.test(FDis ~ strata, data = FDis) # Double-check with Kruskal for non-normal data
# Tukey pairwise test
TukeyHSD(model) # All pairwise comparisons except COMP EXCL-THINNED are significantly different

# Functional composition, as measured by the community-level weighted means of trait values,
# which for continuous trais is the mean trait value of all species present in the community,
# and for categorical traits the abundance of each individual class, i.e. the proportion.
cwms = functcomp(trait_matrix, community_matrix, CWM.type = "all")  # default
cwm_df = data.frame(cwms, strata = strata)
cwm_df$strata = factor(cwm_df$strata, levels = c("STAND INIT", "COMP EXCL", "THINNED", "MATURE"))

ggplot(cwm_df, aes(x = strata, y = mass, fill = strata)) +
  geom_boxplot() + scale_fill_manual(values = strata_cols)

cwm_nest_ps = cwm_df %>%
  select(starts_with("group_nest_ps"), strata) %>%
  pivot_longer(cols = -strata, names_to = "type", values_to = "proportion")
ggplot(cwm_nest_ps, aes(x = strata, y = proportion, fill = fct_rev(type))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent_format())

cwm_migrant = cwm_df %>%
  select(starts_with("group_migrant"), strata) %>%
  pivot_longer(cols = -strata, names_to = "type", values_to = "proportion")
ggplot(cwm_migrant, aes(x = strata, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format())

## Fourth-corner analysis
#
# The fourth-corner method measures and test relationships between
# species functional traits and environmental variables.
#
# Q: Do overal trait distributions differ among stages?
# Null: Distribution of nesting strategies is the same across all stages.
R = data.frame(strata = strata)
L = as.data.frame(community_matrix)
Q = as.data.frame(trait_matrix)
Q$group_migrant = factor(Q$group_migrant)
Q$group_nest_ps = factor(Q$group_nest_ps)
Q$group_diet = factor(Q$group_diet)
fourth = fourthcorner(R, L, Q, nrepet = 9999)
summary(fourth) # Nesting strategy overall differs among strata

# Which nesting strategies are associated with which stages?
# Q: Are particular nesting strategies overrepresented in a stage compared to all other stages?
dat_long = L %>% mutate(site = rownames(L), strata = R$strata) %>% pivot_longer(cols = -c(site, strata), names_to = "species", values_to = "occ") %>% filter(occ > 0)  # keep only present species
nest_df = data.frame(species = rownames(Q), group = Q$group_nest_ps)
dat_long = dat_long %>% left_join(nest_df, by = "species")
tab = table(dat_long$strata, dat_long$group)
results = list()
for (nest_type in colnames(tab)) {
  for (strata_level in rownames(tab)) {
    # counts in this strata
    a = tab[strata_level, nest_type]
    b = sum(tab[strata_level, ]) - a
    # counts in all other strata
    other_rows = setdiff(rownames(tab), strata_level)
    c = sum(tab[other_rows, nest_type])
    d = sum(tab[other_rows, ]) - c
    # Fisher's exact test for (overrepresented) count data
    mat = matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    p = fisher.test(mat, alternative = "greater")$p.value
    results = rbind(results, data.frame(
      strata = strata_level, nest_type = nest_type, p_value = p
    ))
  }
}
results$p_adj = p.adjust(results$p_value, method = "holm") # Adjust for multiple testing
(overrep_nest = results %>% filter(p_adj < 0.05))

# TODO: Explore relationship between species rarity (baseline occupancy) and stage

# Composition PERMANOVA and NMDS ---------------------------------------------------------------------

## Test for significant differences in (expected) community composition of different stages

# Get posterior mean occurrence probability site x species matrix across all posterior draws and seasons
str(z)
psi_mean = apply(z, c(2,4), mean)
rownames(psi_mean) = sites
colnames(psi_mean) = species

# PERMANOVA
# Model R2 is the % of variation in species composition across sites explained by strata as a predictor
# If Pr(>F) < 0.05, differences in composition between strata are statistically significant (at least one stratum differs from the others)
# If not, null hypothesis is true: there is no difference in species composition among strata
permanova = adonis2(psi_mean ~ strata, method = "bray")
permanova

# Pairwise PERMANOVA
strata_levels = unique(strata)
combs = combn(strata_levels, 2, simplify = FALSE)
pairwise_permanova = tibble()
for (pair in combs) {
  i = strata %in% pair
  psi_pair = psi_mean[i, ]
  strata_pair = strata[i]
  permanova_pair = adonis2(psi_pair ~ strata_pair, method = "bray", permutations = 999)
  pairwise_permanova = rbind(pairwise_permanova, tibble(
    stratum_1 = pair[1],
    stratum_2 = pair[2],
    R2        = permanova_pair$R2[1],
    F         = permanova_pair$F[1],
    p.val     = permanova_pair$`Pr(>F)`[1]
  ))
}
# Comparisons with small R2 host largely similar assemblages, while those
# with large R2 have very different assemblages.
# Even a small effect (low R2) can be highly significant (i.e. detectable) if groups are consistently different,
# but these differences are minor when compared to differences between strata with higher R2.
print(pairwise_permanova %>% arrange(R2))

# Dispersion (mean distance of sites to centroid)
dist_bc = vegdist(psi_mean, method = "bray")
bd = betadisper(dist_bc, strata)
anova(bd) # Does dispersion differ among strata? If so, PERMANOVA differences may reflect variation in within-stratum dispersion
dispersion = tibble(stratum = names(bd$group.distances), mean_centroid_dist = tapply(bd$distances, strata, mean))
print(dispersion) # Which strata have the most heterogeneous within-stratum assemblages / beta diversity (high values) and the most homogenous (low values)?

# Mean dissimilarity (beta diversity) among sites within the same stage
diss_within_stratum = tibble()
for (s in levels(strata)) {
  m = as.matrix(dist_bc)[which(strata == s), which(strata == s)]
  diss_within_stratum = rbind(diss_within_stratum, tibble(stratum = s, mean_diss_bc = mean(m[upper.tri(m)])))
}
print(diss_within_stratum) # Which strata has greater within-stage variation?

# Mean dissimilarity (beta diversity) among sites between stages
diss_between_strata = tibble()
combs = combn(levels(strata), 2, simplify = FALSE)
for (pair in combs) {
  m = as.matrix(dist_bc)[which(strata == pair[1]), which(strata == pair[2])]
  diss_between_strata = rbind(diss_between_strata, tibble(stratum_1=pair[1], stratum_2=pair[2], mean_diss_bc=mean(m)))
}
diss_between_strata # Which pairs of strata have the highest/lowest turnover?

# Turnover vs nestedness decomposition
within_turnover = tibble()
for (draw in 1:dim(z)[1]) {
  print(draw)
  for (t in seq(seasons)) {
    z_draw = z[draw, , t, ]
    beta_res = beta.pair(z_draw, index.family = "jaccard")
    for (s in levels(strata)) {
      i = which(strata == s)
      mat_turnover = as.matrix(beta_res$beta.jtu)[i, i]
      mat_nested   = as.matrix(beta_res$beta.jne)[i, i]
      mat_total    = as.matrix(beta_res$beta.jac)[i, i]
      within_turnover = rbind(
        within_turnover,
        tibble(draw = draw, season = t, stratum = s,
               turnover   = mean(mat_turnover[upper.tri(mat_turnover)]),
               nestedness = mean(mat_nested[upper.tri(mat_nested)]),
               mean_diss  = mean(mat_total[upper.tri(mat_total)]))
      )
    }
  }
}
within_turnover_summary = within_turnover %>% group_by(stratum) %>%
  summarise(
    total_mean      = mean(mean_diss, na.rm = TRUE),
    total_lower      = quantile(mean_diss, 0.025, na.rm = TRUE),
    total_upper      = quantile(mean_diss, 0.975, na.rm = TRUE),
    turnover_mean   = mean(turnover, na.rm = TRUE),
    turnover_lower   = quantile(turnover, 0.025, na.rm = TRUE),
    turnover_upper   = quantile(turnover, 0.975, na.rm = TRUE),
    nestedness_mean = mean(nestedness, na.rm = TRUE),
    nestedness_lower = quantile(nestedness, 0.025, na.rm = TRUE),
    nestedness_upper = quantile(nestedness, 0.975, na.rm = TRUE)
  ) %>% arrange(desc(total_mean))
print(within_turnover_summary)

between_turnover = tibble()
for (draw in 1:dim(z)[1]) {
  print(draw)
  for (t in seq(seasons)) {
    z_draw = z[draw, , t, ]
    beta_res = beta.pair(z_draw, index.family = "jaccard")
    combs = combn(levels(strata), 2, simplify = FALSE)
    for (pair in combs) {
      i = which(strata == pair[1])
      j = which(strata == pair[2])
      mat_turnover = as.matrix(beta_res$beta.jtu)[i, j]
      mat_nested   = as.matrix(beta_res$beta.jne)[i, j]
      mat_total    = as.matrix(beta_res$beta.jac)[i, j]
      between_turnover = rbind(
        between_turnover,
        tibble(stratum_1  = pair[1],
               stratum_2  = pair[2],
               turnover   = mean(mat_turnover),
               nestedness = mean(mat_nested),
               mean_diss  = mean(mat_total))
      )
    }
  }
}
between_turnover_summary = between_turnover %>% group_by(stratum_1, stratum_2) %>%
  summarise(
    total_mean      = mean(mean_diss, na.rm = TRUE),
    total_lower      = quantile(mean_diss, 0.025, na.rm = TRUE),
    total_upper      = quantile(mean_diss, 0.975, na.rm = TRUE),
    turnover_mean   = mean(turnover, na.rm = TRUE),
    turnover_lower   = quantile(turnover, 0.025, na.rm = TRUE),
    turnover_upper   = quantile(turnover, 0.975, na.rm = TRUE),
    nestedness_mean = mean(nestedness, na.rm = TRUE),
    nestedness_lower = quantile(nestedness, 0.025, na.rm = TRUE),
    nestedness_upper = quantile(nestedness, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(desc(total_mean))
print(between_turnover_summary)
# Between closed-canopy stages alone, nestedness is low (richness is similar) and turnover is
# more important, meaning that species replacement occurs along the developmental gradient.
# Between stand initiation and closed-canopy stages, nestedness is high (richness is dissimilar,
# as closed-canopy sites host reduced species subsets of open sites) and turnover is moderate,
# further reflecting species replacement across the developmental gradient.
# The higher overall richness of stand init sites can be attributed to a combination of
# nestedness (stand init includes several overlapping species) and turnover (novel early colonizers?)?
# https://www.sciencedirect.com/science/article/pii/S0378112711005779
between_turnover %>% arrange(desc(mean_diss))

# NMDS ----------------------------------------------------------------------------------

# Using posterior mean occupancy probabilities (site similarity in expected composition)
nmds_occprob_k2 = metaMDS(psi_mean, distance = "bray", trymax = 200, k=2)
nmds_occprob_k2

nmds_occprob_k3 = metaMDS(psi_mean, distance = "bray", trymax = 200, k=3)
nmds_occprob_k3

nmds_strata_cols = c("forestgreen", "tan4", "orange", "purple")
strata = site_key$stratum[ match(sites, site_key$site) ]
colors = strata_cols[strata]

ordiplot(nmds_occprob_k3,type="n")
ordiellipse(nmds_occprob_k3, groups = strata, draw = "polygon", kind = "sd", conf = 0.95, lwd = 2, col = NA, border = nmds_strata_cols, label = FALSE)
points(nmds_occprob_k3, display="sites",col=colors, pch = 16, cex=0.75)
# orditorp(nmds_occprob_k3,display="species",col="darkgray",air=0.1, cex=0.75)

# TODO: Visualize species vectors only for those that are consistent significant indicators (see below)

# Functional NMDS

# Select numeric trait columns only (exclude 'strata')
traits_numeric <- cwm_df %>% select(where(is.numeric))

# Calculate NMDS using Bray-Curtis distance
nmds_res <- metaMDS(traits_numeric, distance = "bray", k = 3, trymax = 100)
ordiplot(nmds_res, type = "n")
ordiellipse(nmds_res, groups = strata, draw = "polygon", kind = "sd", conf = 0.95, lwd = 2, col = NA, border = nmds_strata_cols, label = FALSE)
points(nmds_res, display="sites",col=colors, pch = 16, cex=0.75)
points(nmds_res, display="sites", pch = 16, cex=1)
ordihull(nmds_res, groups = strata, draw = "polygon", label = TRUE)

# Indicator species analysis -----------------------------------------------
# https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html

message("Conducting indicator species analysis")

## Indicator species analysis
nperm = 10000

# TODO: Instead, do this over posterior draws and summarize

occ = as.data.frame(z_binary)
rownames(occ) = sites
colnames(occ) = species
head(occ)

# Multi-stratum indicator analysis (De Cáceres, M., Legendre, P., Moretti, M. 2010)
message("Starting multi-stratum indicator analysis (", Sys.time(), ")"); t_start = Sys.time()
indval_multigroup = multipatt(occ, strata, duleg = FALSE, control = how(nperm=nperm))
message("Finished (", Sys.time() - t_start, ")")
# Indicator species associated with a single stratum are specialists of that stratum, while
# those associated with a combination of strata are generalists across those strata.
# Species with significant p-values are strongly associated with a stratum or combination of strata
# Specificity: "Component ‘A’ is sample estimate of the probability that the surveyed site belongs to the target site group given the fact that the species has been found. This conditional probability is called the specificity or positive predictive value of the species as indicator of the site group." Species with A = 1.0 occur exclusively in sites of a given stratum. In other words, species with high A values are largely restricted to a given stratum.
# Fidelity: "Component ‘B’ is sample estimate of the probability of finding the species in sites belonging to the site group. This second conditional probability is called the fidelity or sensitivity of the species as indicator of the target site group." Species with B = 1.0 appear in all sites belonging to a given stratum, i.e. it is widespread in the stratum.
summary(indval_multigroup, indvalcomp=TRUE)
# Show results for all species, regardless of test significance
# Note that of the X species in our dataset, only Y are selected for analysis. The remaining Z species
# "have their highest IndVal value for the set of all sites. In other words, those species occur in sites
# belonging to all groups. The association with the set of all sites cannot be statistically tested, because
# there is no external group for comparison."
summary(indval_singlegroup, alpha=1)
# The species that are not selected can be seen via $sign with NA p.value
indval_multigroup$sign

# NOTE: Species like brown creeper, pileated woodpecker, have similar IndVal values across multiple groups (e.g. THINNED and MATURE), causing them to jump between group assignments in single group analyses due to stocasticity between runs. Using a multi-group analysis is more appropriate.
# indval_singlegroup = multipatt(occ, strata, duleg = TRUE, control = how(nperm=nperm))

# Estimate richness by group -------------------------------------------------------------------------

samples = dim(z)[1]
J = dim(z)[2] # n sites
T = dim(z)[3] # n seasons
I = dim(z)[4] # n species

richness_results = list()
for (grouping in c("group_all", "group_migrant", "group_nest_ps", "group_size")) {
  message("Calculating richness for grouping: ", grouping)
  
  species_group = tibble(
    common_name = species_traits$common_name,
    group = species_traits[[grouping]],
    group_idx = as.integer(as.factor(group))
  ) %>%
  filter(common_name %in% species) %>%
  slice(match(species, common_name))
  
  G = length(unique(species_group$group))
  
  rich_group = array(NA, c(samples, J, T, G))
  for (g in unique(species_group$group)) {
    species_in_g = which(species %in% (species_group %>% filter(group == g) %>% pull(common_name)))
    gi = (species_group %>% filter(group == g) %>% pull(group_idx))[1]
    rich_group[ , , , gi] = apply(z[ , , , species_in_g, drop = FALSE], c(1,2,3), sum)
  }
  rich_group_mean  = apply(rich_group, c(2,3,4), mean)
  rich_group_lower = apply(rich_group, c(2,3,4), quantile, probs = 0.025)
  rich_group_upper = apply(rich_group, c(2,3,4), quantile, probs = 0.975)
  
  site_richness = apply(rich_group_mean, c(1,3), mean) # mean site richness per site averaged across years
  rownames(site_richness) = sites
  colnames(site_richness) = species_group %>% arrange(group_idx) %>% pull(group) %>% unique()
  
  year_richness = apply(rich_group_mean, c(2,3), mean) # mean site richness per year averaged across sites
  rownames(year_richness) = seasons
  colnames(year_richness) = species_group %>% arrange(group_idx) %>% pull(group) %>% unique()
  
  d = as_tibble(site_richness, rownames = "site") %>% left_join(site_key %>% select(site, stratum), by = "site")
  dl = d %>% pivot_longer(cols = where(is.numeric), names_to = "group", values_to = "richness") %>% 
    mutate(stratum = factor(stratum, levels = c("STAND INIT", "COMP EXCL", "THINNED", "MATURE")))
  
  p = ggplot(dl %>% filter(stratum %in% unique(strata)), aes(x = stratum, y = richness, fill = group)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # dodge so boxes for each group don't overlap
    scale_fill_brewer(palette = "Set2") +  # optional, nice color palette
    facet_wrap(~ group, scales = "free_y") +
    labs(title = grouping, x = "Site", y = "Species richness", fill = "Guild"); print(p)
  
  richness_results[[grouping]] = site_richness
}

source("src/global.R")
path_y                                   = "data/cache/4_msom/1_assemble_msom_data/y.rds"
path_trait_data                          = "data/cache/2_traits/1_agg_traits/trait_data.csv"

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)

message("Loading species detection histories 'y' from ", path_y)
y = readRDS(path_y)

species = dimnames(y)$species
seasons = dimnames(y)$season
surveys = dimnames(y)$survey
sites   = dimnames(y)$site

message("Total number of putative detections across all species: ", sum(y > 0, na.rm = TRUE))
message("Total number of putative detections per species:")
sort(apply(y, 4, function(x) sum(x > 0, na.rm = TRUE)), decreasing = TRUE) %>% print()

species_traits = species_traits %>% filter(common_name %in% species)

library(rotl)
library(ape)

sci_names <- species_traits$scientific_name  # already have these!

# Capitalize genus (Open Tree needs "Corvus brachyrhynchos" not "corvus brachyrhynchos")
sci_names_cap <- gsub("(^[a-z])", "\\U\\1", sci_names, perl = TRUE)
sci_names_cap[sci_names_cap == "Corthylio calendula"] <- "Regulus calendula"

# Match to Open Tree of Life
taxa <- tnrs_match_names(sci_names_cap, context_name = "Birds")

# Check for problems
taxa[taxa$approximate_match == TRUE | is.na(taxa$ott_id), ]

# Build tree
tree <- tol_induced_subtree(ott_ids = ott_id(taxa), label_format = "name")

# Clean tip labels
tree$tip.label <- gsub("_", " ", tree$tip.label)

plot(tree, cex = 0.6, no.margin = TRUE)

library(ggtree)

tree_tbl <- as_tibble(tree) |> as.treedata()

ggtree(tree_tbl, layout="fan", open.angle=180) +
  geom_tiplab(size = 2.5, fontface = "italic")

ggtree(tree_tbl, layout="roundrect") +
  geom_tiplab(size = 2.5, fontface = "italic")

name_lookup <- setNames(species_traits$common_name, species_traits$scientific_name)
name_lookup["regulus calendula"]       <- "ruby-crowned kinglet"
name_lookup["hesperiphona vespertina"] <- "evening grosbeak"
name_lookup["picoides villosus"]       <- "hairy woodpecker"
name_lookup["picoides pubescens"]      <- "downy woodpecker"

tree_common = tree
tree_common$tip.label <- str_to_sentence(name_lookup[tolower(tree$tip.label)])
ggtree(tree_common, layout="fan", open.angle=180) +
  geom_tiplab(size = 2.5, fontface = "italic")

####################

library(ggtree)
library(tidytree)
library(RColorBrewer)

# --- Build tip annotation data frame ---
# Map tip labels (common names) back to order
tip_data <- data.frame(
  label = tree_common$tip.label,
  stringsAsFactors = FALSE
) %>%
  left_join(
    species_traits %>%
      mutate(label = str_to_sentence(common_name)) %>%
      select(label, order),
    by = "label"
  )

# --- Color palette for orders ---
orders      <- sort(unique(tip_data$order))
n_orders    <- length(orders)
order_cols <- setNames(hue_pal()(length(orders)), orders)

# --- Plot with tip colors, then propagate to branches ---
p <- ggtree(tree_common, layout = "fan", open.angle = 180) %<+% tip_data +
  geom_tiplab(aes(color = order), size = 2.5, fontface = "italic") +
  geom_tippoint(aes(color = order), size = 1.5) +
  scale_color_manual(values = order_cols, name = "Order", na.value = "grey50") +
  theme(legend.position = "right")
p

# --- Fix: use a name that doesn't clash with ape internals ---
tree_grouped <- groupOTU(tree_common, order_groups, group_name = "tax_order")

# Check it worked - should see tax_order as a factor on the phylo
attr(tree_grouped, "tax_order") |> levels()

# --- Plot ---
library(scales)

p_branches <- ggtree(tree_grouped, aes(color = tax_order),
                     layout = "fan", open.angle = 180) %<+% tip_data +
  geom_tiplab(aes(color = tax_order), size = 2.5, fontface = "italic") +
  scale_color_manual(
    values = c("0" = "gray50", order_cols),
    name   = "Order",
    breaks = orders
  ) +
  theme(legend.position = "bottom")

p_branches


library(ggtext)

colored_labels <- setNames(
  paste0("<span style='color:", order_cols[orders], "'>", orders, "</span>"),
  orders
)

p_branches = ggtree(tree_grouped, aes(color = tax_order),
                     # layout = "fan", open.angle = 180
                    layout="roundrect", branch.length = "none"
                    ) %<+% tip_data +
  geom_tiplab(aes(color = tax_order), size = 2.5, fontface = "italic") +
  scale_color_manual(
    values = c("0" = "gray50", order_cols),
    name   = "Order",
    breaks = orders,
    labels = colored_labels
  ) +
  guides(color = guide_legend(
    override.aes  = list(alpha = 0),   # hide the line/point glyphs
    label.theme   = element_markdown(), # render HTML in labels
    nrow          = 2                   # wrap across 2 rows at bottom
  )) +
  xlim(0, 50) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text()
  ); p_branches

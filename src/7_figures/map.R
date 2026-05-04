source("src/global.R")

library(ggspatial)

path_msom = "data/cache/models/prefinal_msom_jags_nofp_all.rds"
message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)
(sites = model_data$sites)
rm(model_data)

pnt_sites = aru_sites %>% filter(site %in% sites)

rast_elev = rast("data/environment/elevation/elevation.tif") %>% project(crs_m_rast)
rast_elev = mask(rast_elev, vect(study_area))

rast_agg <- aggregate(rast_elev, fact = 5)  # 5x reduction
df <- as.data.frame(rast_agg, xy = TRUE)
colnames(df) <- c("x", "y", "z")

rast_cover = rast("data/cache/3_gis/2_gen_cover_rasters/rast_cover_2020_clean_strata_4.tif")

df_cover <- as.data.frame(rast_cover, xy = TRUE)
colnames(df_cover) <- c("x", "y", "stage")
df_cover$stage <- factor(df_cover$stage, levels = c("standinit", "compex", "thin", "mature", "water", "road_paved", "other"))

mask_outside <- st_difference(study_area, st_union(landscape_planning_units))
bbox <- st_bbox(st_buffer(pnt_sites, 1000))

fig_study_area = ggplot() +
  geom_raster(data = df_cover, aes(x, y, fill = stage), alpha = 0.95) +
  scale_fill_manual(
    values = c(
      colors_map,
      "lightblue1", "gray30", "lightblue1"),
    na.value = "transparent",
    name = "Forest Stage"
  ) +
  geom_relief(data = df, aes(x, y, z = z), sun.angle = 45, alpha = 0.2) +
  geom_sf(data = mask_outside, fill = "white", color = NA, alpha = 0.4, inherit.aes = FALSE) +
  geom_sf(data = st_union(landscape_planning_units), fill = NA, color = "gray30", linewidth = 0.25, inherit.aes = FALSE) +
  geom_sf(data = pnt_sites, shape = 21, color = "gray10", fill = "gray10", size = 1) +
  annotation_scale(location = "br", width_hint = 0.15) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_minimal,
                         height = unit(0.9, "cm"), width = unit(0.9, "cm")
  ) +
  coord_sf(crs = crs_m,
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  labs(x = "", y = "") +
  theme(
    # axis.text  = element_blank(),
    # axis.ticks = element_blank(),
    legend.position = "none"); print(fig_study_area)

ggsave("data/cache/figs/fig_study_area.pdf", fig_study_area,
       width = 9, height = 7)

# ---------------------------------------------------------------------

# wadnr_units_all = st_read("data/environment/GIS Data/WA_DNR_Units/WA_DNR_Units.shp", quiet = TRUE) %>%
#   janitor::clean_names() %>% st_transform(crs = crs_m)

wadnr_units_local = st_read("data/environment/GIS Data/WA_DNR_Units/WA_DNR_Units.shp", quiet = TRUE) %>%
  janitor::clean_names() %>% st_transform(crs = crs_m) %>% # select units that were sampled
  # mutate(jurisdic_2 = ifelse(jurisdic_2 == "Upper Clearwate", "Upper Clearwater", jurisdic_2)) %>%
  filter(jurisdic_2 %in% c("Willy - Huel", "Kalaloch", "Copper Mine", "Upper Clearwate", "Queets", "Goodman Creek", "Reade Hill", "Upper Sol Duc", "Dickodochtedar"))

wadnr_parcels_all = st_read("data/environment/GIS Data/WA_DNR_Managed_Land_Parcels/WA_DNR_Managed_Land_Parcels.shp", quiet = TRUE) %>%
  janitor::clean_names() %>% st_transform(crs = crs_m)

wadnr_ownership = st_intersection(wadnr_units_local, st_union(wadnr_parcels_all))

tribal_lands = st_read("/Users/giojacuzzi/repos/avian-acoustic-monitoring-oesf/data/environment/GIS Data/Tribal_Lands/Tribal_Lands.shp") %>%
  janitor::clean_names() %>% st_transform(crs = crs_m) %>% filter(
    name %in% c("Quinault Indian Reservation", "Hoh Indian Reservation", "Hoh Tribe Lands", "Quileute Indian Reservation")
  )

public_lands = st_read("/Users/giojacuzzi/repos/avian-acoustic-monitoring-oesf/data/environment/GIS Data/WA_Major_Public_Lands_(non-DNR)/WA_Major_Public_Lands_(non-DNR).shp", quiet = TRUE) %>%
  janitor::clean_names() %>% st_transform(crs = crs_m) %>% filter(
    name %in% c("Olympic National Forest", "Olympic National Park")
  )

mapview(wadnr_ownership, zcol = "jurisdic_2") + mapview(public_lands) + mapview(tribal_lands)

library(elevatr)
elev_rast <- get_elev_raster(
  locations = st_make_valid(public_lands) |> st_transform(4326) |> st_union() |> st_sf(),
  z = 9,
  clip = "bbox"
)
wa_counties_cb = counties(state = 'WA', cb = T) %>% st_transform(4326) %>% filter(
  NAME %in% c("Clallam", "Jefferson", "Grays Harbor")
)

counties_sp <- as(wa_counties_cb, "Spatial")   # raster needs sp, not sf

elev_cropped <- crop(elev_rast, counties_sp)    # clip to bounding box
elev_masked  <- mask(elev_cropped, counties_sp)

df_elev_all <- elev_masked |>
  rast() |>                          # RasterLayer → SpatRaster
  project(paste0("EPSG:", st_crs(landscape_planning_units)$epsg)) |>
  as.data.frame(xy = TRUE) |>
  setNames(c("x", "y", "z"))

fig_region = ggplot() +
  geom_relief(data = df_elev_all, aes(x, y, z = z), sun.angle = 45, alpha = 0.2) +
  geom_sf(data = public_lands %>% filter(name == "Olympic National Park"), fill = "tan4", alpha = 0.5, color = "gray30", linewidth = 0.25, inherit.aes = FALSE) +
  geom_sf(data = public_lands %>% filter(name == "Olympic National Forest"), fill = "forestgreen", alpha = 0.5, color = "gray30", linewidth = 0.25, inherit.aes = FALSE) +
  geom_sf(data = tribal_lands, fill = "skyblue", alpha = 0.5, color = "gray30", linewidth = 0.25, inherit.aes = FALSE) +
  geom_sf(data = wadnr_ownership, fill = "royalblue", alpha = 0.5, color = "gray30", linewidth = 0.25, inherit.aes = FALSE) +
  geom_sf(data = st_union(landscape_planning_units), fill = NA, color = "gray10", linewidth = 0.4, inherit.aes = FALSE) +
  coord_sf(crs = crs_m,
           xlim = c(bbox["xmin"] - 5000, bbox["xmax"] + 5000),
           ylim = c(bbox["ymin"] - 15000, bbox["ymax"] + 5000)) +
  labs(x = "", y = "") +
  theme(legend.position = "bottom"); print(fig_region)

ggsave("data/cache/figs/fig_region.pdf", fig_region,
       width = 4, height = 6)

stop()

# Home range site example
potential_sites = c("ap308i", "bp246i", "ca261i", "dz077i", "dz301i", "aa019i", "ap308i", "az024i", "cp226i",
                    "aa145i", "ba194i", "ba199i", "ca250i", "bp233i")
pnt = pnt_sites %>% filter(site == potential_sites[11]) # 8, 4, 11, 12, 13
pnt_plot_buffer = st_buffer(pnt, 100)
pnt_homerange_buffer = st_buffer(pnt, 1000)
bbox_homerange = st_bbox(pnt_homerange_buffer)
pnt_cover = crop(mask(rast_cover, pnt_homerange_buffer), pnt_homerange_buffer)

elev_df <- as.data.frame(crop(mask(rast_agg, pnt_homerange_buffer), pnt_homerange_buffer), xy = TRUE)
colnames(elev_df) <- c("x", "y", "z")

fig_homerange = ggplot() +
  geom_raster(data = pnt_cover, aes(x, y, fill = stage), alpha = 0.95) +
  scale_fill_manual(
    values = c(
      "thin" = colors_map[3],
      "standinit" = colors_map[1],
      "compex" = colors_map[2],
      "mature" = colors_map[4],
      "water" = "lightblue1",
      "road_paved" = "gray30",
      "other" = "lightblue1"),
    na.value = "transparent",
    name = "Forest Stage"
  ) +
  geom_relief(data = elev_df, aes(x, y, z = z), sun.angle = 45, alpha = 0.2) +
  geom_sf(data = pnt_homerange_buffer, color = "gray10", fill = NA, linewidth = 0.5) +
  geom_sf(data = pnt_plot_buffer, color = "gray10", fill = NA, linewidth = 0.5) +
  geom_sf(data = pnt, shape = 21, color = "gray10", fill = "gray10", size = 1) +
  coord_sf(crs = crs_m,
           xlim = c(bbox_homerange["xmin"], bbox_homerange["xmax"]),
           ylim = c(bbox_homerange["ymin"], bbox_homerange["ymax"])) +
  theme_void() +
  theme(legend.position = "right"); print(fig_homerange)

ggsave("data/cache/figs/fig_homerange.pdf", fig_homerange)
  
# rast_QMD_6 = load_raster("data/environment/rsfris_study_area/4.0/QMD_6.tif")
# pnt_qmd = crop(mask(rast_QMD_6, pnt_homerange_buffer), pnt_homerange_buffer)
# ggplot() +
#   geom_raster(data = pnt_qmd, aes(x, y, fill = Layer_1), alpha = 0.95) +
#   geom_sf(data = pnt_homerange_buffer, color = "gray10", fill = NA, linewidth = 0.5) +
#   geom_sf(data = pnt_plot_buffer, color = "gray10", fill = NA, linewidth = 0.5) +
#   geom_sf(data = pnt, shape = 21, color = "gray10", fill = "gray10", size = 1) +
#   coord_sf(crs = crs_m,
#            xlim = c(bbox_homerange["xmin"], bbox_homerange["xmax"]),
#            ylim = c(bbox_homerange["ymin"], bbox_homerange["ymax"])) +
#   theme_void() +
#   theme(legend.position = "bottom")

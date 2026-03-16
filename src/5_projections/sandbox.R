source("src/global.R")

# hcu <- st_read("data/environment/GIS Data/Habitat_Conservation_Plan_Units/Habitat_Conservation_Plan_Units.shp")
# units = st_read("data/environment/GIS Data/WA_DNR_Units/WA_DNR_Units.shp")
# parcels = st_read("data/environment/GIS Data/WA_DNR_Managed_Land_Parcels/WA_DNR_Managed_Land_Parcels.shp")
# hu = st_read("data/environment/GIS Data/harvest/T3uplandHarvestUnits/HarvestUnitsUpland.shp")
# u = st_read("data/environment/GIS Data/Watershed_Administrative_Units_-_Forest_Practices_Regulation/Watershed_Administrative_Units_-_Forest_Practices_Regulation.shp")

riu = st_read("data/environment/rs_fris/Polygon_RS-FRIS/Polygon_RS-FRIS.shp")

study_area_proj <- st_transform(study_area, st_crs(riu))
riu_crop <- riu[st_intersects(riu, study_area_proj, sparse = FALSE), ]

# UPDATE

path = "data/environment/oesf_rdeis_sof_act_20110419.gdb"
st_layers(path)

# State of the forest (developmental stages)

sof = st_read(dsn = path, layer = "sof")
str(sof)

study_area_proj = st_transform(study_area, st_crs(sof))
sof_decade_0 = sof %>% filter(DECADE == 0, ALTERNATIVE == "Landscape")
sof_decade_0_crop = sof_decade_0[st_intersects(sof_decade_0, study_area_proj, sparse = FALSE), ]
str(sof_decade_0_crop)

grp = factor(sof_decade_0_crop$OESF_SDS_GROUPED)
pal = setNames(rainbow(nlevels(grp)), levels(grp))
plot(st_geometry(sof_decade_0_crop), col = pal[grp], border = NA)
legend("topright", legend = names(pal), fill = pal, cex = 0.7)

# Activities (thinning and harvest)

act = st_read(dsn = path, layer = "act")
str(act)
study_area_proj = st_transform(study_area, st_crs(act))
act_decade_1 = act %>% filter(DECADE == 1, ALTERNATIVE == "Landscape")
act_decade_1_crop = act_decade_1[st_intersects(act_decade_1, study_area_proj, sparse = FALSE), ]
str(act_decade_1_crop)

mapview(act_decade_1_crop, zcol = "HARVEST_TYPE_NM")
grp = factor(act_decade_1_crop$HARVEST_TYPE_NM)
pal = setNames(rainbow(nlevels(grp)), levels(grp))
plot(st_geometry(act_decade_1_crop), col = pal[grp], border = NA)
legend("topright", legend = names(pal), fill = pal, cex = 0.7)

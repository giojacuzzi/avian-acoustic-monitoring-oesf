source("src/global.R")

# hcu <- st_read("data/environment/GIS Data/Habitat_Conservation_Plan_Units/Habitat_Conservation_Plan_Units.shp")
# units = st_read("data/environment/GIS Data/WA_DNR_Units/WA_DNR_Units.shp")
# parcels = st_read("data/environment/GIS Data/WA_DNR_Managed_Land_Parcels/WA_DNR_Managed_Land_Parcels.shp")
# hu = st_read("data/environment/GIS Data/harvest/T3uplandHarvestUnits/HarvestUnitsUpland.shp")
# u = st_read("data/environment/GIS Data/Watershed_Administrative_Units_-_Forest_Practices_Regulation/Watershed_Administrative_Units_-_Forest_Practices_Regulation.shp")

riu = st_read("data/environment/rs_fris/Polygon_RS-FRIS/Polygon_RS-FRIS.shp")

study_area_proj <- st_transform(study_area, st_crs(riu))
riu_crop <- riu[st_intersects(riu, study_area_proj, sparse = FALSE), ]

# ACT

st_layers("data/projections/FEIS_ACT_FGDB_20160201.gdb")

act_table <- st_read(
  "data/projections/FEIS_ACT_FGDB_20160201.gdb",
  layer = "FEIS_ACT_20150806"
)

act_poly <- st_read(
  "data/projections/FEIS_ACT_FGDB_20160201.gdb",
  layer = "R12302010_5TH"
)

names(act_poly)
names(act_table)

str(act_poly)
str(act_table)

study_area_proj <- st_transform(study_area, st_crs(act_poly))
act_poly_crop <- act_poly[st_intersects(act_poly, study_area_proj, sparse = FALSE), ]
mapview(act_poly_crop)

table(act_table$ALTERNATIVE)
length(unique(act_poly$REMSOFT_ID))

str(act_poly_crop)
str(act_table)

# Forestry activities during decades 1 and 2 under Pathway alternative
test = full_join(act_poly_crop, act_table %>%
                   filter(ALTERNATIVE == "Pathway") %>%
                   select(REMSOFT_ID, ACT_DECADE), by = "REMSOFT_ID")
mapview(test %>% filter(ACT_DECADE == 1), col.regions = "blue") +
  mapview(test %>% filter(ACT_DECADE == 2), col.regions = "red")

test = full_join(act_poly_crop, act_table %>%
                   filter(ALTERNATIVE == "Pathway") %>%
                   select(REMSOFT_ID, ACT_DECADE, HARVEST_TYPE_NM), by = "REMSOFT_ID")
mapview(test %>% filter(ACT_DECADE == 1), zcol = "HARVEST_TYPE_NM")

# Pathway

st_layers("data/projections/pathwayFGDB_20160106.gdb")

poly = st_read(
  "data/projections/pathwayFGDB_20160106.gdb",
  layer = "fcPATHWAY_picks_20160106"
)

names(poly)
str(poly)

study_area_proj <- st_transform(study_area, st_crs(poly))
poly_crop <- poly[st_intersects(poly, study_area_proj, sparse = FALSE), ]
mapview(poly_crop) + mapview(study_area_proj)

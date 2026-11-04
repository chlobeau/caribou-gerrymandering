require(TuktuTools)

load("./ms_varyingArea/results/summary_perm.rda")
load("./ms_varyingArea/results/summary_fidelity.rda")

load("./data/studyAreaBounds2025.rda")
load(file = "./data/nwt_lakes.rda")
load(file = "./data/varyingBounds.rda")

harvest$Harvest[harvest$Harvest == "Unsettled Land Claims Unit U"] <- "Unsettled Land Claims Zone U"

provinces <- st_read("./data/raw/studyAreaBounds/provinces/lpr_000a16a_e.shp") |>
  st_transform(st_crs(bounds))

# get coastline
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(st_crs(bounds))

renamePolys <- function(poly_sf, bound_name){
  poly_sf |>
    dplyr::select(bound_name) |>
    dplyr::rename(grouping_col = bound_name) |>
    mutate(bound = bound_name)
}

all <- rbind(renamePolys(bounds, "StudyArea"),
             renamePolys(admin, "Admin"),
             renamePolys(landscape, "Landscape"),
             renamePolys(rangeplan, "RangePlan"),
             renamePolys(harvest, "Harvest"),
             renamePolys(bounds2, "StudyArea2"))

used_polys = bind_rows(
  perm_poly_summary |>
    subset(included_perm) |>
    mutate(approach = "resample"),
  fidelity_poly_summary |>
    subset(included_fidelity) |>
    mutate(approach = "fidelity")) |> 
  dplyr::select(bound, grouping_col) |>
  mutate(grouping_col = ifelse(grouping_col == "Unsettled Land Claims Unit U",
                               "Unsettled Land Claims Zone U", grouping_col))

polys = all |>
  subset(grouping_col %in% used_polys$grouping_col)

# display fix: merge features that share a name (removes internal lines),
# then close the river gap in the two merged StudyArea2 areas
close_gap <- function(geom, d = 2000) {
  geom |>
    sf::st_buffer(d,  joinStyle = "MITRE", mitreLimit = 5) |>
    sf::st_buffer(-d, joinStyle = "MITRE", mitreLimit = 5)
}

polys <- polys |>
  dplyr::group_by(bound, grouping_col) |>
  dplyr::summarise(.groups = "drop")   # unions geometries within each group

fix <- polys$bound == "StudyArea2" & polys$grouping_col %in% c("Dehcho North & MVH", "Dehcho South & W. Liard")
sf::st_geometry(polys)[fix] <- close_gap(sf::st_geometry(polys)[fix])

bb <- st_bbox(polys |> subset(bound == "StudyArea")) |> st_as_sfc()

bb_inset <- st_bbox(c(xmin = -1900000, xmax = 1600000,
                      ymin = 340000,  ymax = 4100000),
                    crs = st_crs(polys)) |> st_as_sfc()
provinces_crop <- st_intersection(provinces, world)

# prep ----
## define colors ----

used_polys |>
  dplyr::select(bound, grouping_col) |>
  arrange(bound) |>
  unique()

zone_color_lookup <- c(
  # --- StudyArea / StudyArea2 shared identities ---
  "Cameron Hills"                      = "#CC79A7",
  "Mackenzie"                          = "#009E73",   # StudyArea + StudyArea2 (confirm same geom)
  "Pine Point - Buffalo Lake"          = "#0072B2",
  "Dehcho South"                       = "#56B4E9",
  "Dehcho South & W. Liard"            = "#56B4E9",   # merged identity, same color
  "Dehcho North"                       = "#E69F00",
  "Dehcho North & MVH"                 = "#E69F00",   # merged identity, same color
  "North Slave (TASR)"                 = "#D55E00",
  "MVH"                                = "#F7D999",
  "Hay River Lowlands"                 = "#F0E442",
  "West Liard River"                   = "#C3E4F6",
  
  # --- Admin ---
  "Dehcho"                             = "#7F3C8D",
  "South Slave"                        = "#11A579",
  "North Slave"                        = "#3969AC",
  
  # --- Landscape ---
  "East of Hay River"                               = "#0072B2",
  "Dehcho South - Hay River Lowlands"               = "#F0E442",
  "NE of Mackenzie River - S of Great Bear River"   = "#80BA5A",
  "W of Liard and  Mackenzie River - S of Keele River" = "#999999",
  
  # --- RangePlan ---
  "Southern NWT" = "#A5AA99",
  "Wekeezhii"    = "#D55E00",
  "Sahtu"        = "black",
  
  # --- Harvest ---
  "Unsettled Land Claims Zone U"   =  "#4D4D4D",  
  "Unsettled Land Claims Zone D"   = "#A6761D",  
  "Wek'èezhìi Boundary and Zone R" = "#E7298A",
  "Sahtu Dene AND Metis Settlement Area" = "black"
)

## define north/south order ----

# compute north-to-south order from polygon centroids
centroid_order <- polys |>
  mutate(centroid_y = st_coordinates(st_centroid(st_geometry(polys)))[, 2]) |>
  st_drop_geometry() |>
  dplyr::select(bound, grouping_col, centroid_y) |>
  unique()

## display labels for legends, where they differ from the raw polygon name ----
zone_label_lookup <- c(
  "Sahtu Dene AND Metis Settlement Area" = "Sahtu Dene & Metis Settlement Area",
  "Wek'èezhìi Boundary and Zone R"        = "Wek'èezhìi Boundary & Zone R"
)

# eccc polygons
eccc <- st_read("./data/raw/studyAreaBounds/ECCC_BorealCaribouRanges.shp") |>
  st_transform(st_crs(bounds))

# reference outline: full data-collection extent
sa_outline <- polys |>
  subset(bound == "StudyArea") |>
  st_buffer(5000) |>
  st_union() |>
  st_buffer(-5000) |>
  st_cast("POLYGON")

# keep outer ring only (drops the interior gap)
sa_outline <- st_sfc(lapply(st_geometry(sa_outline),
                            function(p) st_polygon(list(p[[1]]))),
                     crs = st_crs(polys)) |>
  st_union()

# BC and AB labels (lon/lat converted to the map CRS)
prov_labs <- data.frame(lab = c("BC", "AB"),
                        lon = c(-123, -116),
                        lat = c(59.3, 59.3)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(st_crs(bounds))
prov_xy <- st_coordinates(prov_labs)

# map all management boundaries ----

png("./ms_varyingArea/figs/boundary_schemes.png",
    width = 9, height = 9, units = "in", res = 600)
{
  layout(matrix(c(1, 2, 3,
                  4, 5, 6,
                  7, 7, 7), nrow = 3, byrow = TRUE),
         heights = c(3, 3, 1.5))
  
  par(mar = c(0.5, 0.5, 1.8, 0.5), oma = c(1, 1, 1, 1))
  
  for (g in c("StudyArea", "StudyArea2", "Landscape", 
              "RangePlan", "Harvest", "Admin")) {
    
    sub <- polys |> subset(bound == g)
    cols <- zone_color_lookup[sub$grouping_col]
    
    plot(bb, col = NA, border = NA)
    plot(st_geometry(provinces), col = "grey98", border = "grey50", lwd = 0.6, add = TRUE)
    plot(st_geometry(nwt_lakes), col = "white", border = alpha("grey", 0.5),
         add = TRUE)
    plot(st_geometry(sub), col = alpha(cols, 0.6), 
         border = "grey30", lwd = 0.8, add = TRUE)
    
    if (!g %in% c("StudyArea", "StudyArea2")) {
      plot(sa_outline, col = NA, border = "black", lty = 2, lwd = 1.2, add = TRUE)
    }
    text(prov_xy[, 1], prov_xy[, 2], prov_labs$lab, col = "grey40", font = 3, cex = 0.9)
    
    mtext(g, side = 3, line = 0.2, font = 2, cex = 1.2)
  }
  
  # inset panel 
  par(mar = c(0, 16, 0, 16))
  plot(bb_inset, col = NA, border = NA)
  plot(st_geometry(world), col = "grey98", border = "grey60", lwd = 0.5, add = TRUE)
  plot(st_geometry(provinces_crop), col = NA, border = "grey50", lwd = 0.6, add = TRUE)
  
  plot(st_geometry(eccc), col = alpha("#999933", 0.1), border = "#999933", add = TRUE)
  
  plot(bb, col = alpha("#882255", 0.2), border = "#882255", lwd = 2, add = TRUE)
  box(col = "grey70")
}
dev.off()

# map trajectories ----
load(file = "./data/nt_move_mort_1daily.rda")
load(file = "./ms_varyingArea/results/prop_assigned_all.rda")

move_sa <- merge(nt_resampled |> 
                   st_drop_geometry(),
                 prop_assigned_all |>
                   subset(bound == "StudyArea2"),
                 by = c("ID", "bioYear")) |>
  subset(!is.na(assigned) & assigned != "Calendar")

move_sa$cols <- zone_color_lookup[move_sa$assigned]

unique(move_sa$assigned)

move_h <- merge(nt_resampled |> 
                  st_drop_geometry(),
                prop_assigned_all |>

                  subset(bound == "Harvest"),
                 by = c("ID", "bioYear")) |>
  subset(!is.na(assigned))

move_h$assigned[move_h$assigned == "Unsettled Land Claims Unit U"] <- "Unsettled Land Claims Zone U"

move_h$cols <- zone_color_lookup[move_h$assigned]

unique(move_h$assigned)

# plot together
png("./ms_varyingArea/figs/trajectory_comparison.png",
    width = 10, height = 6, units = "in", res = 600)
{
  par(mfrow = c(1, 2), mar = rep(.2, 4), oma = c(.5, .5, 1.5, .5))
  
  # --- Panel a: StudyArea2 tracks ---
  plot(st_geometry(bb), col = NA, border = NA)
  plot(st_geometry(provinces_crop), col = "grey98", border = "grey50", lwd = 0.6, add = TRUE)
  plot(st_geometry(nwt_lakes), col = "white", border = alpha("grey", .5), add = TRUE)
  plot(st_geometry(polys |>
                     subset(bound == "StudyArea2")), 
       add = TRUE, col = NA, border = "black", lwd = 2)
  
  ddply(move_sa, c("ID", "bioYear"), function(df) {
    with(df |> arrange(Date), lines(x, y, col = alpha(cols, 0.5)))
    NULL
  })
  mtext("StudyArea2", side = 3, line = 0.2, font = 2, cex = 1.2)
  
  sa_zones <- centroid_order |>
    subset(bound == "StudyArea2") |>
    dplyr::arrange(dplyr::desc(centroid_y)) |>
    dplyr::pull(grouping_col)
  
  legend("topright", 
         legend = sa_zones, 
         pt.bg = zone_color_lookup[sa_zones], 
         pch = 22, pt.cex = 2, cex = .75,
         bg = "white")
  
  # --- Panel b: Harvest tracks ---
  plot(st_geometry(bb), col = NA, border = NA)
  plot(st_geometry(provinces_crop), col = "grey98", border = "grey50", lwd = 0.6, add = TRUE)
  plot(st_geometry(nwt_lakes), col = "white", border = alpha("grey", .5), add = TRUE)
  plot(st_geometry(polys |>
                     subset(bound == "Harvest")), 
       add = TRUE, col = NA, border = "black", lwd = 2)
  
  ddply(move_h |>
          subset(assigned != "Unsettled Land Claims Zone U"), 
        c("ID", "bioYear"), function(df) {
          with(df |> arrange(Date), lines(x, y, col = alpha(cols, 0.5)))
          NULL
          })
  # map U last so it shows up over D
  ddply(move_h |>
          subset(assigned == "Unsettled Land Claims Zone U"), 
        c("ID", "bioYear"), function(df) {
          with(df |> arrange(Date), lines(x, y, col = alpha(cols, 0.5)))
          NULL
        })
  mtext("Harvest", side = 3, line = 0.2, font = 2, cex = 1.2)
  
  harvest_zones <- centroid_order |>
    subset(bound == "Harvest") |>
    dplyr::arrange(dplyr::desc(centroid_y)) |>
    dplyr::pull(grouping_col)
  
  legend("topright", 
         legend = ifelse(harvest_zones %in% names(zone_label_lookup),
                         zone_label_lookup[harvest_zones],
                         harvest_zones),
         pt.bg = zone_color_lookup[harvest_zones], 
         pch = 22, pt.cex = 2, cex = .75,
         bg = "white")
}
dev.off()


# supporting ----

## map individual boundary scheme ----

for (g in c("StudyArea", "StudyArea2", "Landscape", 
            "RangePlan", "Harvest", "Admin")) {
  
  sub <- polys |> subset(bound == g)
  
  legend_names <- centroid_order |>
    subset(bound == g) |>
    dplyr::arrange(dplyr::desc(centroid_y)) |>
    dplyr::pull(grouping_col)
  legend_cols  <- zone_color_lookup[legend_names]
  
  cols <- zone_color_lookup[sub$grouping_col]
  
  png(paste0("./ms_varyingArea/figs/boundary_schemes/", g, ".png"),
      width = 7, height = 6, units = "in", res = 300)
  {
    par(mar = c(0.5, 0.5, 1.8, 0.5))
    
    plot(bb, col = NA, border = NA)
    plot(st_geometry(provinces), col = "grey98", border = "grey50", lwd = 0.6, add = TRUE)
    plot(st_geometry(nwt_lakes), col = "white", border = alpha("grey", 0.5), add = TRUE)
    plot(st_geometry(sub), col = alpha(cols, 0.6), 
         border = "grey30", lwd = 0.8, add = TRUE)
    mtext(g, side = 3, line = 0.2, font = 2, cex = 1.2)
    
    legend("topright",
           legend = legend_names,
           pt.bg  = alpha(legend_cols, 0.6),
           pch    = 22,
           pt.cex = 1.8,
           cex    = 0.85,
           bg     = "white")
    
    # --- inset map, bottom left ---
    par(fig = c(0, 0.3, 0, 0.21), mar = c(0, 0, 0, 0), new = TRUE)
    plot(bb_inset, col = NA, border = NA)
    plot(st_geometry(world), col = "grey98", border = "grey60", lwd = 0.5, add = TRUE)
    plot(st_geometry(provinces_crop), col = NA, border = "grey50", lwd = 0.6, add = TRUE)
    plot(bb, col = alpha("#7B2D8B", 0.2), border = "#7B2D8B", lwd = 2, add = TRUE)
    box(col = "grey70")
  }
  dev.off()
}

## tables ----

## combined bound-level summary table
# n_obs = polygon-years for Cross-validation/Permutations, animal-years for
# Boundary fidelity
bound_summary_all <- dplyr::bind_rows(
  perm_bound_summary,
  fidelity_bound_summary
) |>
  subset(bound != "Ecoregion") |>
  dplyr::select(method, bound, n_polygons, n_years, n_obs,
                min_per_polygon, max_per_polygon) |>
  dplyr::arrange(method, bound)

bound_summary_all

## detailed per-polygon table
# same polygon list/ordering as the existing `polygons_by_bound` table,
# with sample sizes sitting alongside each inclusion flag as their own
# columns rather than folded into a string
polygons_by_bound_detailed <- centroid_order |>
  dplyr::left_join(cv_poly_summary,       by = c("bound", "grouping_col")) |>
  dplyr::left_join(perm_poly_summary,     by = c("bound", "grouping_col")) |>
  dplyr::left_join(fidelity_poly_summary, by = c("bound", "grouping_col")) |>
  dplyr::mutate(
    included_cv       = tidyr::replace_na(included_cv, FALSE),
    included_perm     = tidyr::replace_na(included_perm, FALSE),
    included_fidelity = tidyr::replace_na(included_fidelity, FALSE)
  ) |>
  dplyr::arrange(bound, dplyr::desc(centroid_y)) |>
  dplyr::select(-centroid_y) |>
  dplyr::rename(Polygon = grouping_col)

polygons_by_bound_detailed

bound_summary_display <- bound_summary_all |>
  dplyr::arrange(method, bound) |>
  dplyr::select(
    "Method"          = method,
    "Boundary scheme" = bound,
    "n polygons"      = n_polygons,
    "n years"         = n_years,
    "n obs"           = n_obs,
    "min per polygon" = min_per_polygon,
    "max per polygon" = max_per_polygon
  )

readr::write_excel_csv(bound_summary_display,
                       "./ms_varyingArea/text/supporting_bound_summary.csv")

polygons_by_bound_display <- polygons_by_bound_detailed |>
  dplyr::select(
    bound,
    Polygon,
    "CV included"  = included_cv,
    "CV years"          = n_years_cv,
    "Permutations included"      = included_perm,
    "Permutation years" = n_years_perm,
    "Boundary fidelity included" = included_fidelity,
    "Individuals"       = n_individuals,
    "Animal-years"      = n_animal_years
  )

readr::write_excel_csv(polygons_by_bound_display,
                       "./ms_varyingArea/text/supporting_polygons_by_bound.csv")

# prop wetland ----
require(raster)
require(exactextractr)

eosd = raster("./data/raw/EOSD/EOSD_NWT1.tif")
crs(eosd)
crs(bounds)
eosd

bounds_t = st_transform(bounds, crs(eosd))

plot(eosd)
plot(bounds_t, add = TRUE)

ex = exact_extract(eosd, bounds_t)

counts = bind_rows(ex) |>
  count(value)

vals = counts |> subset(value != "NA")

sum(vals$n[vals$value >= 81 & vals$value <= 83]) / sum(vals$n)

# area size ----

tab = polys |>
  #subset(grouping_col %in% fidelity_poly_summary$grouping_col) |>
  group_by(bound) |>
  summarize(n = length(unique(grouping_col)),
            area_km2 = mean(st_area(geometry))/1e6,
            sd_km2 = sd(st_area(geometry))/1e6) |>
  st_drop_geometry()

sum(as.numeric(st_area(bounds)) / 1e6)     # total km2

tab |>taboundsb |>
  mutate(`Area, km2 (mean ± sd)` = paste0(
    format(round(as.numeric(area_km2)), big.mark = ","), " ± ",
    format(round(as.numeric(sd_km2)), big.mark = ","))) |>
  dplyr::select(bound, n, `Area, km2 (mean ± sd)`)

require(TuktuTools)
require(taigatools) #https://github.com/chlobeau/taigatools
require(glmmTMB)

load("./data/nt_move_mort_1daily.rda") # movement data are not publicly available
# the workflow needs a dataframe with one location per animal per day
# (a 1-daily resampled track) with the following columns:
# ID - individual animal id that matches ID in ids_annual
# bioYear - must match bioYear in ids_annual
# geometry - an sf POINT column with defined CRS

load("./ms_varyingArea/code/data_review/ids_annual.rda")
load("./ms_varyingArea/code/data_review/studyAreaBounds2025.rda")
load("./ms_varyingArea/code/data_review/varyingBounds.rda")
head(ids_annual)

ids_annual$StudyArea2 |> unique()

ids_annual = ids_annual |>
  subset(!(StudyArea2 %in% c("Calendar", "Maxhamish")))

# calculate prop assigned ----

# makes sure polygon names in ids_annual match polygon names in sf objects
make_names_tbl <- function(ids_col) {
  assigned_names <- unique(ids_col[!is.na(ids_col)])
  data.frame(
    assigned     = assigned_names,
    assigned_new = vctrs::vec_as_names(assigned_names, repair = "universal")
  )
}

# verify it matches existing names_tbl for StudyArea2
make_names_tbl(ids_annual$StudyArea2)

bounds_list <- list(
  StudyArea  = list(sf = bounds, col = "StudyArea"),
  StudyArea2 = list(sf = bounds2, col = "StudyArea2"),
  Landscape  = list(sf = landscape, col = "Landscape"),
  Admin      = list(sf = admin, col = "Admin"),
  RangePlan  = list(sf = rangeplan, col = "RangePlan"),
  Harvest    = list(sf = harvest, col = "Harvest"),
  Ecoregion  = list(sf = ecoregion, col = "Ecoregion")
)

eval = FALSE

if(eval){
  prop_assigned_all <- ldply(names(bounds_list), function(bound_name) {
    
    bounds_sf  <- bounds_list[[bound_name]]$sf
    bound_col  <- bounds_list[[bound_name]]$col
    names_tbl  <- make_names_tbl(ids_annual[[bound_col]])
    
    # merge locations with assigned polygon for this bound
    locs <- merge(nt_resampled,
                  ids_annual |> dplyr::select(ID, bioYear, all_of(bound_col)),
                  by = c("ID", "bioYear")) |>
      st_transform(st_crs(bounds_sf))
    
    pts_summary <- ddply(locs, "bioYear", function(df) {
      countPtsInPolygons(
        pts_sf  = st_as_sf(df),
        bounds = bounds_sf,
        bounds_col = bound_col,
        assigned_col = bound_col
      )
    })
    
    pts_summary <- merge(pts_summary, names_tbl,
                         by.x = "assigned", by.y = "assigned",
                         all.x = TRUE)
    
    # calculate prop_assigned
    pts_summary$prop_assigned <- NA
    for (i in 1:nrow(pts_summary)) {
      df <- pts_summary[i, ]
      if (!is.na(df$assigned_new)) {
        col_name <- paste0(df$assigned_new, "_days")
        if (col_name %in% names(df)) {
          pts_summary$prop_assigned[i] <- 
            as.numeric(df[[col_name]]) / df$total_days
        }
      }
    }
    
    pts_summary$bound <- bound_name
    pts_summary
  })
  
  save(prop_assigned_all, file = "./ms_varyingArea/results/prop_assigned_all.rda")
  
}else load(file = "./ms_varyingArea/results/prop_assigned_all.rda")

# clean data ----

prop_complete <- prop_assigned_all |>
  select(ID, bioYear, bound, assigned, prop_assigned) |>
  subset(!is.na(assigned))

# unique animal-years with complete prop_assigned across all six schemes
n_distinct(prop_complete$ID, prop_complete$bioYear)

# unique individuals
length(unique(prop_complete$ID))

table(prop_complete$bioYear)

prop_complete |>
  group_by(bound) |>
  dplyr::summarise(
    n_polygons = n_distinct(assigned),
    n_ids = n_distinct(ID)
  )

# save summary stats (for supporting information) ----

# per-polygon sample sizes 
fidelity_poly_summary <- prop_complete |>
  mutate(grouping_col = assigned) |>
  group_by(bound, grouping_col) |>
  dplyr::summarise(
    n_individuals  = n_distinct(ID),
    n_animal_years = n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(included_fidelity = TRUE)  # no sample-size filter applied

# bound-level summary 
fidelity_bound_summary <- prop_complete |>
  group_by(bound) |>
  dplyr::summarise(
    n_polygons      = n_distinct(assigned),
    n_years         = n_distinct(bioYear),
    n_obs           = n(),
    min_per_polygon = min(table(assigned)),
    max_per_polygon = max(table(assigned)),
    .groups = "drop"
  ) |>
  dplyr::mutate(method = "Boundary fidelity")

save(fidelity_poly_summary,
     fidelity_bound_summary, 
     file = "./ms_varyingArea/results/summary_fidelity.rda")

# model prop_assigned ~ boundary ----

fit_prop <- glmmTMB::glmmTMB(prop_assigned ~ bound-1 + (1|ID),
                             data = prop_complete,
                             family = "binomial")
summary(fit_prop)

## pairwise comparisons among boundary schemes ----
emmeans::emmeans(fit_prop, pairwise ~ bound)$contrasts |>
  as.data.frame() |>
  arrange(p.value)

## proportion of animal-years where animal crossed its boundary at all ----
prop_complete |>
  group_by(bound) |>
  dplyr::summarise(
    prop_crossing           = mean(prop_assigned < 1),
    mean_prop_when_crossing = mean(prop_assigned[prop_assigned < 1], na.rm = TRUE),
    n                       = n()
  ) |>
  arrange(prop_crossing)

# look at individuals that crossed
fit_crossed <- glmmTMB(prop_assigned ~ bound - 1 + (1 | ID),
                       data = prop_complete |>
                         subset(prop_assigned < 1),
                       family = binomial)

summary(fit_crossed)

# pairwise contrasts, Tukey-adjusted
emmeans::emmeans(fit_crossed, pairwise ~ bound)$contrasts |>
  as.data.frame() |>
  arrange(p.value)
# no significant differences between groups - leaving out of the ms

# summary table: mean/median proportion assigned per bound ----
prop_complete |>
  group_by(bound) |>
  dplyr::summarise(
    mean_prop   = round(mean(prop_assigned,   na.rm = TRUE), 3),
    median_prop = round(median(prop_assigned, na.rm = TRUE), 3),
    n           = n(),
    .groups     = "drop"
  ) |>
  arrange(desc(mean_prop))

# final formatted summary table
prop_complete |>
  group_by(bound) |>
  dplyr::summarise(
    mean_prop     = round(mean(prop_assigned, na.rm = TRUE), 3),
    pct_crossing  = round(mean(prop_assigned < 1, na.rm = TRUE) * 100, 1),
    #mean_crossing = round(mean(prop_assigned[prop_assigned < 1], na.rm = TRUE), 3),
    n             = n(),
    .groups = "drop"
  ) |>
  mutate(bound = factor(bound, 
                        levels = c("StudyArea", "StudyArea2", "Landscape", "Ecoregion",
                                   "RangePlan", "Admin", "Harvest"))) |>
  arrange(bound) |>
  dplyr::rename(
    `Management boundary scheme`    = bound,
    `Mean proportion assigned`      = mean_prop,
    `% animal-years crossing`       = pct_crossing,
    #`Mean proportion when crossing` = mean_crossing,
    `n`                             = n
  ) |> write.csv(file = "./ms_varyingArea/results/propAssignedCrossing.csv", row.names = FALSE)

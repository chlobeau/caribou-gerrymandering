require(TuktuTools)
require(taigatools) #https://github.com/chlobeau/taigatools
require(survival)
require(glmmTMB)

# remove locations from comp, survival, herd assignment data for review
# load("./data/cleanComp_sf.rda")
# cleanComp_sf = cleanComp_sf |>
#   st_drop_geometry() |>
#   dplyr::select(-c("x", "y"))
# save(cleanComp_sf, file = "./ms_varyingArea/code/data_review/cleanComp_sf.rda")

# load("./data/cleanSurv_20251217.rda")
# survival <- survival |>
#   dplyr::select(ID, startDate, endDate, Status, Sex, Died, days)
# save(survival, file = "./ms_varyingArea/code/data_review/cleanSurv_public.rda")

# load("./data/herd_assignments.rda")
# ids_annual <- ids_annual |>
#   dplyr::select(-c("centroid_x", "centroid_y"))
# save(ids_annual, file = "./ms_varyingArea/code/data_review/ids_annual.rda")

load("./ms_varyingArea/code/data_review/cleanComp_sf.rda")
load("./ms_varyingArea/code/data_review/ids_annual.rda")
load("./ms_varyingArea/code/data_review/cleanSurv_public.rda")

sample_size_threshold = 15

grouping_cols = c("StudyArea", "StudyArea2", "Landscape",
                  "RangePlan", "Admin", "Harvest")

# prep comp survey data
cleanComp_sf <- cleanComp_sf |>
  subset(!(StudyArea2 %in% c("Calendar", "Maxhamish")))
unique(cleanComp_sf$StudyArea2)

# prep annual ids data for survival
ids_annual <- ids_annual |>
  subset(!(StudyArea2 %in% c("Calendar", "Maxhamish")) &
           prop_assigned >= 0.8)

# first calculate monthly survival
monthly_morts <- ddply(survival |> 
                         subset(!is.na(endDate)) |>
                         subset(Sex == "Female"), 
                       "ID", setupMonthlySurvival)

monthly_morts <- addBioYear(monthly_morts |> mutate(Time = month_end), bioStart = 91) |>
  dplyr::select(-Time)

run_permutation <- function(i, grouping_cols) {
  
  message("Permutation ", i)
  set.seed(i)
  
  # shuffle survival
  ids_shuffled <- ids_annual
  years_ids <- unique(ids_annual$bioYear)
  ids_shuffled[grouping_cols] <- lapply(grouping_cols, function(g) {
    unlist(lapply(years_ids, function(yr) {
      idx <- ids_annual$bioYear == yr
      sample(ids_annual[[g]][idx])
    }))
  })
  
  surv_perm <- ldply(grouping_cols, function(g) {
    mm_g <- merge(monthly_morts, 
                  ids_shuffled[, c("ID", "bioYear", g)], 
                  by = c("ID", "bioYear"))
    monthly_km  <- ddply(mm_g |> subset(!is.na(get(g))),
                         c("Year", "month", g), getSurvTable)
    
    getAnnualSurvival(monthly_km, grouping_col = g)
  })
   
  # shuffle comp (no filter)
  comp_shuffled <- as.data.frame(cleanComp_sf)
  years_comp    <- unique(comp_shuffled$bioYear)
  comp_shuffled[grouping_cols] <- lapply(grouping_cols, function(g) {
    unlist(lapply(years_comp, function(yr) {
      idx <- comp_shuffled$bioYear == yr
      sample(comp_shuffled[[g]][idx])
    }))
  })
  comp_shuffled <- st_as_sf(comp_shuffled)
  
  rec_perm <- ldply(grouping_cols, function(g) {
    calculateRecruitment(comp_shuffled, grouping_col = g)
  })
  
  lambda_perm <- merge(
    rec_perm  |> dplyr::select(bound, grouping_col, Years,
                               CalfCow, RecruitedFemales, Rrm),
    surv_perm |> dplyr::select(bound, grouping_col, Years,
                               km_survival, km_var, n_survival),
    by = c("bound", "grouping_col", "Years")) |>
    subset(n_survival > sample_size_threshold) |> 
    mutate(
      lambda = km_survival / (1 - Rrm),
      w = n_survival / mean(n_survival, na.rm = TRUE)
    ) |>
    subset(!is.na(lambda))
  
  lambda_perm
}

fit_perm_models <- function(lambda_perm, i, grouping_cols){
  
  ldply(grouping_cols, function(g) {
    df <- lambda_perm |> subset(bound == g)
    if (nrow(df) < 5) return(NULL)
    
    fit_surv <- suppressWarnings(
      glmmTMB(km_survival ~ grouping_col + (1 | Years),
              data = df, weights = w, family = "binomial"))
    fit_rec  <- suppressWarnings(
      glmmTMB(CalfCow ~ grouping_col + (1 | Years),
              data = df, weights = w, family = gaussian))
    fit_lam  <- suppressWarnings(
      glmmTMB(lambda ~ grouping_col + (1 | Years),
              data = df |> subset(!is.na(lambda)),
              weights = w, family = gaussian))
    
    data.frame(
      bound = g,
      iter = i,
      r2_surv = performance::r2(fit_surv)$R2_marginal,
      r2_rec = performance::r2(fit_rec)$R2_marginal,
      r2_lambda = performance::r2(fit_lam)$R2_marginal
    )
  })
}

# resample in parallel ! ----
eval = FALSE
if(eval){
  r1 <- run_permutation(1, grouping_cols)
  r2 <- run_permutation(2, grouping_cols)
  identical(r1$km_survival, r2$km_survival)  # should be FALSE
  identical(r1, r2)
  
  require(parallel)
  
  cl <- makeCluster(20)
  
  clusterExport(cl, varlist = c(
    "run_permutation",
    "grouping_cols",
    "ids_annual",
    "cleanComp_sf",
    "monthly_morts",
    "sample_size_threshold",
    "calculateRecruitment",
    "getAnnualSurvival"
  ))
  
  clusterEvalQ(cl, {
    library(dplyr)
    library(plyr)
    library(sf)
    library(taigatools)
    library(survival)
    library(TuktuTools)
  })
  
  system.time(
    lambda_perm_list <- parLapply(cl, 1:499, run_permutation,
                                  grouping_cols = grouping_cols)
  )
  
  # test_parallel <- parLapply(cl, 1:6, run_permutation,
  #                            grouping_cols = grouping_cols) |>
  #   bind_rows()
  
  stopCluster(cl)
  
  save(lambda_perm_list, file = "./ms_varyingArea/results/lambda_perm_list.rda")
  
}else load(file = "./ms_varyingArea/results/lambda_perm_list.rda")

# fit models to resampled data ----

if(eval){
  
  cl <- makeCluster(20)
  
  clusterExport(cl, varlist = c(
    "fit_perm_models",
    "grouping_cols",
    "lambda_perm_list"
  ))
  
  clusterEvalQ(cl, {
    library(dplyr)
    library(plyr)
    library(glmmTMB)
    library(performance)
  })
  
  system.time(
    null_r2 <- parLapply(cl, seq_along(lambda_perm_list), function(i)
      fit_perm_models(lambda_perm_list[[i]], i, grouping_cols))
  )
  
  stopCluster(cl)
  
  null_r2 <- ldply(null_r2)
  save(null_r2, file = "./ms_varyingArea/results/null_r2_new.rda")
}else load(file = "./ms_varyingArea/results/null_r2_new.rda")

# calculate observed R2 ----

# compare structure of permuted vs observed data
# use first permutation as representative example

lambda_observed_all <- ldply(grouping_cols, function(g) {
  
  mm_g <- merge(monthly_morts,
                ids_annual[, c("ID", "bioYear", g)],
                by = c("ID", "bioYear"))
  
  monthly_km  <- ddply(mm_g |> subset(!is.na(get(g))),
                       c("Year", "month", g), getSurvTable)
  annual_surv <- getAnnualSurvival(monthly_km, grouping_col = g)
  
  rec <- calculateRecruitment(cleanComp_sf, grouping_col = g)
  
  merge(rec |>
          dplyr::select(grouping_col, Years,
                        CalfCow, RecruitedFemales, Rrm),
        annual_surv |>
          dplyr::select(grouping_col, Years,
                        km_survival, km_var, n_survival),
        by = c("grouping_col", "Years")) |>
    dplyr::mutate(bound = g,
                  lambda = km_survival / (1 - Rrm),
                  w = n_survival / mean(n_survival, na.rm = TRUE))
})

lambda_observed <- lambda_observed_all |>
  subset(n_survival > sample_size_threshold)

save(lambda_observed, file = "./ms_varyingArea/results/lambda_observed.rda")

observed_r2 <- fit_perm_models(lambda_observed, i = 0, grouping_cols = grouping_cols)

# save summary stats ----

# row counts
nrow(lambda_observed) # total polygon-year observations, all six bounds combined
nrow(lambda_perm_list[[1]])

# bounds represented
table(lambda_observed$bound)
table(lambda_perm_list[[1]]$bound)

# grouping_cols per bound
tapply(lambda_observed$grouping_col, lambda_observed$bound, function(x) length(unique(x)))
tapply(lambda_perm_list[[1]]$grouping_col, lambda_perm_list[[1]]$bound, function(x) length(unique(x)))

# years of data
table(lambda_observed$Years)

# mean rows per bound across all permutations
perm_df <- ldply(seq_along(lambda_perm_list), function(i) {
  df <- lambda_perm_list[[i]]
  df$iter <- i
  df
})

# mean rows and unique groups per bound across permutations
ddply(perm_df, "bound", function(df) {
  per_iter <- ddply(df, "iter", function(x) {
    data.frame(n_rows   = nrow(x),
               n_groups = length(unique(x$grouping_col)))
  })
  data.frame(mean_rows   = mean(per_iter$n_rows),
             mean_groups = mean(per_iter$n_groups))
})

ddply(lambda_observed, "bound", function(df) {
  data.frame(n_rows   = nrow(df),
             n_groups = length(unique(df$grouping_col)))
})

## per-polygon sample sizes
perm_poly_summary <- lambda_observed_all |>
  dplyr::select(bound, grouping_col) |>
  unique() |>
  dplyr::left_join(
    lambda_observed |>
      dplyr::group_by(bound, grouping_col) |>
      dplyr::summarise(n_years_perm = dplyr::n(), .groups = "drop"),
    by = c("bound", "grouping_col")
  ) |>
  dplyr::mutate(
    n_years_perm  = tidyr::replace_na(n_years_perm, 0),
    included_perm = n_years_perm > 0
  )

## bound-level summary
perm_bound_summary <- lambda_observed |>
  dplyr::group_by(bound) |>
  dplyr::summarise(
    n_polygons      = dplyr::n_distinct(grouping_col),
    n_years         = dplyr::n_distinct(Years),
    n_obs           = dplyr::n(),
    min_per_polygon = min(n_survival),
    max_per_polygon = max(n_survival),
    .groups = "drop"
  ) |>
  dplyr::mutate(method = "Permutations")

save(perm_poly_summary,
     perm_bound_summary, 
     file = "./ms_varyingArea/results/summary_perm.rda")

# compare R2 ----
# p-value = proportion of null R2 >= observed R2
# higher observed R2 than null = boundary explains more variation than chance
perm_summary <- null_r2 |>
  dplyr::group_by(bound) |>
  dplyr::summarise(
    null_mean_surv   = mean(r2_surv,   na.rm = TRUE),
    null_mean_rec    = mean(r2_rec,    na.rm = TRUE),
    null_mean_lambda = mean(r2_lambda, na.rm = TRUE),
    p_surv   = (sum(r2_surv   >= observed_r2$r2_surv[match(bound,   observed_r2$bound)], na.rm = TRUE) + 1) / (sum(!is.na(r2_surv))   + 1),
    p_rec    = (sum(r2_rec    >= observed_r2$r2_rec[match(bound,    observed_r2$bound)], na.rm = TRUE) + 1) / (sum(!is.na(r2_rec))    + 1),
    p_lambda = (sum(r2_lambda >= observed_r2$r2_lambda[match(bound, observed_r2$bound)], na.rm = TRUE) + 1) / (sum(!is.na(r2_lambda)) + 1)
  ) |>
  left_join(observed_r2, by = "bound")

perm_summary

# figure ----

bio_col   <- "#1b9e77"  # teal-green, colorblind safe
admin_col <- "#d95f02"  # orange, colorblind safe

bound_df <- data.frame(
  bound = c("StudyArea", "StudyArea2", "Landscape", "RangePlan", "Admin", "Harvest"),
  col   = c(bio_col, bio_col, bio_col, admin_col, admin_col, admin_col),
  grp   = c("Biological", "Biological", "Biological",
            "Administrative", "Administrative", "Administrative"),
  row.names = NULL
)

vital_list <- list(
  list(vr = "r2_surv",   p = "p_surv",   label = "Survival",
       ymax = 0.16),
  list(vr = "r2_rec",    p = "p_rec",    label = "Recruitment",
       ymax = 0.65),
  list(vr = "r2_lambda", p = "p_lambda", label = "Population growth rate (\u03bb)",
       ymax = 0.5)
)

png("./ms_varyingArea/figs/resampling_r2.png",
   width = 10, height = 4.5, units = "in", res = 600)

#win.graph(width = 10, height = 4.5)
{
  par(mfrow = c(1, 3),
      mar   = c(7, 2, 3, 1),
      oma   = c(2, 3, 0, 0),
      las   = 1,
      tcl   = 0.3,
      mgp   = c(3, 0.3, 0),
      bty   = "l")
  
  for (v in seq_along(vital_list)) {
    
    # merge observed and perm_summary with bound_df ordering
    obs  <- merge(observed_r2,  bound_df, by = "bound")
    ps   <- merge(perm_summary, bound_df, by = "bound")
    obs  <- obs[match(bound_df$bound, obs$bound), ]
    ps   <- ps[match(bound_df$bound, ps$bound), ]
    
    plot(0, type = "n",
         xlim = c(0.5, 6.5),
         ylim = c(0, vital_list[[v]]$ymax),
         xlab = "",
         ylab = if (v == 1) expression(R^2) else "",
         main = vital_list[[v]]$label,
         font.main = 2,
         yaxs = "i",
         xaxt = "n",
         xpd  = NA,
         cex.lab = 1.3)
    
    for (b in 1:6) {
      bound_name <- bound_df$bound[b]
      null_vals  <- null_r2[[vital_list[[v]]$vr]][null_r2$bound == bound_name]
      null_vals  <- null_vals[!is.na(null_vals)]
      
      # violin
      if (length(null_vals) > 10) {
        dens   <- density(null_vals, adjust = .6)
        dscale <- dens$y / max(dens$y) * 0.35
        polygon(c(b - dscale, rev(b + dscale)),
                c(dens$x, rev(dens$x)),
                col    = alpha(bound_df$col[b], 0.2),
                border = NA)
      }
      
      # boxplot overlay
      bxp <- boxplot.stats(null_vals)$stats
      segments(b, bxp[1], b, bxp[2], col = bound_df$col[b], lty = 1, lwd = 1)
      segments(b, bxp[4], b, bxp[5], col = bound_df$col[b], lty = 1, lwd = 1)
      rect(b - 0.08, bxp[2], b + 0.08, bxp[4],
           col = alpha(bound_df$col[b], 0.15), 
           border = bound_df$col[b], lwd = 1.2)
      segments(b - 0.08, bxp[3], b + 0.08, bxp[3],
               col = bound_df$col[b], lwd = 2)
      
      # observed point - filled if sig, open if not
      p_val  <- ps[[vital_list[[v]]$p]][b]
      pt_bg  <- if (!is.na(p_val) && p_val < 0.05) bound_df$col[b] else alpha("white", .6)
      points(b, obs[[vital_list[[v]]$vr]][b],
             pch = 21,
             bg  = pt_bg,
             col = bound_df$col[b],
             cex = 2,
             lwd = 2)
    }
    
    # x axis labels
    text(x = 1:6,
         y = par("usr")[3] - vital_list[[v]]$ymax * 0.015,
         labels = bound_df$bound,
         srt = 40,
         adj = 1,
         xpd = TRUE,
         cex = 1.1)
  }
  
  # legend
  par(fig = c(0, 1, 0, 1), 
      oma = c(0,0,0,0), 
      mar = c(0,22,0,22), # this is how I'm centering legends
      new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  legend("bottomleft",
         legend = c("Biological", "Administrative"),
         pch = 21,
         pt.bg = c(bio_col, admin_col),
         col = c(bio_col, admin_col),
         pt.cex = 2,
         bty = "n",
         horiz = TRUE,
         cex = 1.2,
         title = "Management boundary type",
         title.font = 2)
  
  legend("bottomright",
         legend = c("p < 0.05", "p \u2265 0.05"),
         pch = 21,
         pt.bg = c("#636363", "white"),
         pt.lwd = 2,
         col = "#636363",
         pt.cex = 2,
         bty = "n",
         horiz = TRUE,
         cex = 1.2,
         title = "Significance",
         title.font = 2)
}
dev.off()
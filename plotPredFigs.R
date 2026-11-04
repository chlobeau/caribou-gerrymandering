## Prediction figures: survival by year, recruitment by StudyArea2, lambda by Admin

require(dplyr)
require(plyr)
require(tidyr)
require(lubridate)
require(sf)
require(glmmTMB)
require(ggeffects)
require(scales)

# Prep recruitment ----

load("./ms_varyingArea/results/lambda_observed.rda")

d = lambda_observed |>
  subset(bound == "StudyArea2")

fit_r <- glmmTMB(CalfCow ~ grouping_col + (1 | Years),
                 weights = w, 
                 data = d, 
                 family = gaussian)

# Polygon ordering (north to south), StudyArea2 only

load("./data/varyingBounds.rda")

renamePolys <- function(poly_sf, bound_name) {
  poly_sf |>
    dplyr::select(bound_name) |>
    dplyr::rename(grouping_col = bound_name) |>
    mutate(bound = bound_name)
}

combined_polys <- renamePolys(bounds2, "StudyArea2") |>
    subset(grouping_col %in% unique(d$grouping_col))
           
centroid_order <- combined_polys |>
  mutate(centroid_y = st_coordinates(st_centroid(st_geometry(combined_polys)))[, 2]) |>
  st_drop_geometry() |>
  dplyr::select(bound, grouping_col, centroid_y)

poly_order_rec <- centroid_order |>
  subset(bound == "StudyArea2") |>
  dplyr::arrange(centroid_y) |>
  dplyr::pull(grouping_col)

# Fit models and get predictions ----

pred_rec <- as.data.frame(ggpredict(fit_r, terms = "grouping_col")) |>
  dplyr::rename(grouping_col = x)
pred_rec <- pred_rec[match(poly_order_rec, pred_rec$grouping_col), ]

# Figures ----

## recruitment by StudyArea2 ----
png("./ms_varyingArea/figs/pred_recruitment_studyarea2.png",
    width = 7, height = 5, units = "in", res = 600)
{
  par(mar = c(5.5, 11, 1, 1), bty = "l", tcl = 0.3, mgp = c(2, 0.4, 0))
  n <- length(poly_order_rec)
  
  plot(0, type = "n",
       xlim = c(0, max(c(pred_rec$conf.high, d$CalfCow), na.rm = TRUE) * 1.05),
       ylim = c(0.5, n + 0.5),
       xlab = "Recruitment (calf:cow ratio)", ylab = "",
       yaxt = "n")
  axis(2, at = 1:n, labels = poly_order_rec, las = 1)
  
  abline(v = 0.28, lty = 2, col = "grey50")
  
  for (i in seq_along(poly_order_rec)) {
    raw <- d$CalfCow[d$grouping_col == poly_order_rec[i]]
    points(raw, jitter(rep(i, length(raw)), amount = 0.12),
           pch = 19, col = alpha("#BBBBBB", .8))
  }
  
  segments(pred_rec$conf.low, 1:n, pred_rec$conf.high, 1:n, 
           col = "#AA3377", lwd = 1.5)
  points(pred_rec$predicted, 1:n, 
         pch = 19, col = "#AA3377", cex = 2.5)
  
  # legend drawn on a full-figure overlay, below the main plot and its axis title
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("bottom",
         legend = c("Predicted (95% CI)", "Raw polygon-year value"),
         pch    = 19,
         col    = c("#AA3377", "#BBBBBB"),
         pt.cex = c(2.2, 1),
         lwd = c(1.5, NA),
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.95)
}
dev.off()

require(TuktuTools)

load("./ms_varyingArea/results/lambda_observed.rda")
load(file = "./data/nwt_lakes.rda")
load(file = "./data/varyingBounds.rda")

provinces <- st_read("./data/raw/studyAreaBounds/provinces/lpr_000a16a_e.shp") |>
  st_transform(st_crs(bounds2))

bounds2 = bounds2 |>
  subset(!(StudyArea2 %in% c("Maxhamish", "Calendar")))

bounds2 <- bounds2 |>
  sf::st_buffer(2000,  joinStyle = "MITRE", mitreLimit = 5) |>
  sf::st_buffer(-2000, joinStyle = "MITRE", mitreLimit = 5)

require(RColorBrewer)
mypalette <- colorRampPalette(c("#E08214", "#FDB863", "#F2F2F2", "#B2ABD2", "#8073AC"))(101)

lambda_cols = data.frame(colors = mypalette,
                         values = seq(.92, 1.08, length.out = 101))

get_closest_color <- function(fit_value, lambda_df) {
  closest_index <- which.min(abs(lambda_df$values - fit_value))
  return(lambda_df$colors[closest_index])
}

fit_lambda_g <- function(g, sf_obj) {
  d <- lambda_observed |> subset(bound == g & !is.na(lambda))   # n > 15 data
  d$area <- factor(d$grouping_col)
  
  m <- glmmTMB(lambda ~ area + (1 | Years),
               weights = w, data = d, family = gaussian)
  
  pr <- predict(m, newdata = data.frame(area = levels(d$area), Years = NA, w = 1),
                se.fit = TRUE, allow.new.levels = TRUE) |> data.frame()
  pr$grouping_col <- levels(d$area)
  pr$predicted    <- pr$fit
  pr$conf.low     <- pr$fit - 1.96 * pr$se.fit
  pr$conf.high    <- pr$fit + 1.96 * pr$se.fit
  pr$col          <- sapply(pr$fit, get_closest_color, lambda_df = lambda_cols)
  
  # order polygons south to north so the north plots at the top
  cy <- data.frame(grouping_col = sf_obj[[g]],
                   cy = st_coordinates(st_centroid(st_geometry(sf_obj)))[, 2])
  pr <- merge(pr, cy) |> (\(x) x[order(x$cy), ])()
  
  list(map = merge(sf_obj, pr[, c("grouping_col", "col")], by.x = g, by.y = "grouping_col"),
       pr = pr, raw = d)
}

plot_map <- function(r) {
  par(mar = c(.5, .5, 3, .5), bty = "n")
  plot(st_geometry(bounds2), col = NA, border = NA)
  plot(st_geometry(nwt_lakes), col = "white", border = alpha("grey", .5), add = TRUE)
  plot(st_geometry(provinces), border = "grey", lty = 3, add = TRUE)
  plot(st_geometry(r$map), add = TRUE,
       col = r$map$col, border = "black")
  plot(st_geometry(bounds2), col = NA, border = "black", add = TRUE, lty = 3)
  title = names(r$map)[1]
  mtext(title, 3, cex = 1.5, line = .5)
}

plot_dots <- function(r, show_x = FALSE) {
  pr <- r$pr; n <- nrow(pr)
  par(mar = c(if (show_x) 6 else 2, 9, 1, 1),  
      bty = "l", tcl = 0.3, mgp = c(2, 0.4, 0))
  plot(0, type = "n", xlim = xr, ylim = c(0.5, n + 0.5), xaxt = "n", yaxt = "n",
       ylab = "", xlab = "")
  axis(1)
  if (show_x) mtext(expression(paste("Population growth rate (", lambda, ")")),
                    side = 1, line = 2.5, cex = 0.9)
  axis(2, at = 1:n, las = 1,
       labels = sapply(pr$grouping_col, function(x) paste(strwrap(x, 22), collapse = "\n")))
  abline(v = 1, lty = 2, col = "grey50")
  for (i in 1:n) {
    raw <- r$raw$lambda[r$raw$grouping_col == pr$grouping_col[i]]
    points(raw, jitter(rep(i, length(raw)), amount = 0.12),
           pch = 19, col = alpha("#BBBBBB", .8))
  }
  segments(pr$conf.low, 1:n, pr$conf.high, 1:n, lwd = 1.5, col = "grey30")
  points(pr$predicted, 1:n, pch = 21, bg = pr$col, lwd = 1.5, cex = 2.5)
}

lambda <- list("StudyArea2" = fit_lambda_g("StudyArea2", bounds2),
             #"Landscape" = fit_lambda_map("Landscape", landscape),
             "Admin" = fit_lambda_g("Admin", admin))

xr <- range(unlist(lapply(lambda, function(r) c(r$pr$conf.low, r$pr$conf.high, r$raw$lambda))))

for (nm in names(lambda)) plot_map(lambda[[nm]])

for (nm in names(lambda)) plot_dots(lambda[[nm]])

png("./ms_varyingArea/figs/lambda_maps_dots.png", 
    width = 7.5, height = 7.5, 
    units = "in", 
    res = 600)
layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1, 1.2))

plot_map(lambda[["StudyArea2"]])     
plot_dots(lambda[["StudyArea2"]], show_x = FALSE)          
plot_map(lambda[["Admin"]])
plot_dots(lambda[["Admin"]], show_x = TRUE) 

# legend overlay, restricted to the right column so it centers under the dot plots
par(fig = c(0.4545, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom",
       legend = c("Predicted (95% CI)", "Raw polygon-year value"),
       pch    = c(21, 19),
       pt.bg  = c("white", NA),
       col    = c("grey30", "#BBBBBB"),
       pt.cex = c(2.2, 1),
       lty    = c(1, NA),
       lwd    = c(1.5, NA),
       horiz  = TRUE,
       bty    = "n",
       cex    = 1)

dev.off()

lambda$StudyArea2$pr
lambda$Admin$pr

plot_dots_simple <- function(r) {
  pr <- r$pr; n <- nrow(pr)
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot(0, type = "n", xlim = xr, ylim = c(0.5, n + 0.5),
       axes = FALSE, xlab = "", ylab = "")
  abline(v = 1, lty = 2, col = "grey50", lwd = 2)
  segments(pr$conf.low, 1:n, pr$conf.high, 1:n, lwd = 3, col = "grey30")
  points(pr$predicted, 1:n, pch = 21, bg = pr$col, lwd = 2, cex = 3.5)
}

png("./ms_varyingArea/figs/lambda_schematic.png", 
    width = 7.5, height = 7.5, 
    units = "in", 
    res = 600)
layout(matrix(1:4, nrow = 2, byrow = TRUE))

plot_map(lambda[["StudyArea2"]])
plot_dots_simple(lambda[["StudyArea2"]])
plot_map(lambda[["Admin"]])
plot_dots_simple(lambda[["Admin"]])
dev.off()

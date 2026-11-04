require(TuktuTools)
require(glmmTMB)

# annual survival by individual ----
## for binomial models instead of kaplan meier estimates

load("./ms_varyingArea/code/data_review/ids_annual.rda")
load("./ms_varyingArea/code/data_review/cleanSurv_public.rda")

setupAnnualSurvival <- function(mort.df, bioStart = "04-01") {
  
  start_threshold <- as.Date(paste0(year(mort.df$startDate), "-", bioStart))
  end_threshold   <- as.Date(paste0(year(mort.df$endDate),   "-", bioStart))
  
  start_bioYear <- as.Date(paste0(ifelse(mort.df$startDate < start_threshold, 
                                         year(mort.df$startDate) - 1, 
                                         year(mort.df$startDate)), "-", bioStart))
  end_bioYear   <- as.Date(paste0(ifelse(mort.df$endDate < end_threshold,
                                         year(mort.df$endDate) - 1,
                                         year(mort.df$endDate)), "-", bioStart))
  
  year_sequence <- seq.Date(start_bioYear, end_bioYear %m+% years(1), by = "1 year")
  
  data.frame(
    year_start = year_sequence,
    year_end   = year_sequence %m+% years(1) - 1,
    ID         = mort.df$ID,
    Fate       = mort.df$Status,
    startDate  = mort.df$startDate,
    endDate    = mort.df$endDate
  ) |>
    mutate(
      end        = pmin(mort.df$endDate, year_end),
      DeadInYear = Fate == "Mort" & endDate <= year_end,
      CensInYear = Fate == "Cens" & endDate <= year_end,
      bioYear    = as.integer(substr(year_start, 1, 4)),
      duration   = as.numeric(end - startDate)
    )
}

yearly_morts <- survival |>
  subset(Sex == "Female" & !is.na(startDate) & !is.na(endDate)) |>
  ddply("ID", setupAnnualSurvival)

yearly_morts_bounds <- merge(yearly_morts, ids_annual,
                             by = c("ID", "bioYear")) |>
  na.omit()

surv_df <- yearly_morts_bounds |>
  subset(!CensInYear) |>
  mutate(
    survived = as.integer(!DeadInYear),
    bioYear  = as.integer(bioYear)
  )

length(unique(surv_df$ID))
nrow(surv_df)

# year random-intercept variance ----
## from the individual-based survival model

# Year-only model (no boundary scheme), individual-based survival

m_year_only <- glmmTMB(
  survived ~ 1 + (1 | bioYear),
  family = binomial,
  data = surv_df
)

summary(m_year_only)

performance::r2(m_year_only)$R2_conditional

# Pull variance components
vc <- VarCorr(m_year_only)
year_var <- as.numeric(vc$cond$bioYear)

year_var

# intraclass-correlation coefficient
?performance::icc
# ICC can be interpreted as 
# "the proportion of the variance explained by the grouping structure in the population"

year_icc <- performance::icc(m_year_only)
year_icc

# year and area random intercept ----

grouping_cols = c("StudyArea", "StudyArea2", "Landscape",
                  "RangePlan", "Admin", "Harvest")

year_area_vr <- function(response) {
  ldply(grouping_cols, function(g) {
    d <- lambda_observed |> subset(bound == g)
    d$y    <- d[[response]]
    d      <- d |> subset(!is.na(y))
    d$area <- factor(d$grouping_col)
    
    m_y <- suppressWarnings(glmmTMB(y ~ 1 + (1 | Years),
                                    weights = w, data = d, family = gaussian))
    m   <- suppressWarnings(glmmTMB(y ~ 1 + (1 | Years) + (1 | area),
                                    weights = w, data = d, family = gaussian))
    
    vc  <- VarCorr(m)$cond
    b   <- summary(m)$coefficients$cond[1, ]
    lrt <- anova(m_y, m)
    
    year_v <- as.numeric(vc$Years)
    area_v <- as.numeric(vc$area)
    total  <- year_v + area_v + sigma(m)^2
    
    data.frame(delineation = g,
               n_areas     = n_distinct(d$area),
               n_obs       = nrow(d),
               estimate    = as.numeric(b["Estimate"]),
               lower       = as.numeric(b["Estimate"] - 1.96 * b["Std. Error"]),
               upper       = as.numeric(b["Estimate"] + 1.96 * b["Std. Error"]),
               year_pct    = 100 * year_v / total,
               area_pct    = 100 * area_v / total,
               lrt_p = ifelse(is.na(lrt$`Pr(>Chisq)`[2]) & area_v < 1e-6,
                              1, lrt$`Pr(>Chisq)`[2]),
               r2_yr = performance::r2(m_y)$R2_conditional,
               r2_yrare = performance::r2(m)$R2_conditional)
  })
}

make_vr_table <- function(df, est_label) {
  data.frame(
    `Region delineation`      = df$delineation,
    `Polygons (n)`            = df$n_areas,
    est                       = sprintf("%.3f [%.3f-%.3f]",
                                        df$estimate, df$lower, df$upper),
    `Year (% of variance)`    = sprintf("%.1f", df$year_pct),
    `Polygon (% of variance)` = sprintf("%.1f", df$area_pct),
    `LRT p`                   = ifelse(is.na(df$lrt_p), "NA",
                                       ifelse(df$lrt_p < 0.001, "<0.001",
                                              sprintf("%.3f", df$lrt_p))),
    check.names = FALSE
  ) |>
    dplyr::rename(!!est_label := est)
}

## survival ----
year_var   # baseline survival m_year_only

surv_yeararea <- ldply(grouping_cols, function(g) {
  d <- surv_df
  d$area <- factor(d[[g]])
  
  m_y <- glmmTMB(survived ~ 1 + (1 | bioYear),
                 family = binomial, data = d)
  
  # random-intercept model: pooled survival and variance components
  m <- glmmTMB(survived ~ 1 + (1 | bioYear) + (1 | area),
               family = binomial, data = d)
  vc <- VarCorr(m)$cond
  b  <- summary(m)$coefficients$cond[1, ]   # intercept and SE (logit scale)
  ci <- plogis(b["Estimate"] + c(-1.96, 1.96) * b["Std. Error"])
  
  lrt <- anova(m_y, m)
  
  year_v <- as.numeric(vc$bioYear)
  area_v <- as.numeric(vc$area)
  total  <- year_v + area_v + pi^2 / 3   # logistic residual variance
  
  data.frame(delineation = g,
             n_areas     = n_distinct(d$area),
             estimate    = as.numeric(plogis(b["Estimate"])),
             lower       = as.numeric(ci[1]),
             upper       = as.numeric(ci[2]),
             year_pct    = 100 * year_v / total,
             area_pct    = 100 * area_v / total,
             lrt_p       = lrt$`Pr(>Chisq)`[2],
             r2_yr = performance::r2(m_y)$R2_conditional,
             r2_yrare = performance::r2(m)$R2_conditional)
})

## lambda / recruitment ----
# need lambda_observed from resampleVitalRates.R
load("./ms_varyingArea/results/lambda_observed.rda")
rec_yeararea <- year_area_vr("CalfCow")
lam_yeararea <- year_area_vr("lambda")

# ms tables ----

make_vr_table(surv_yeararea, "Survival [95% CI]") #|>
#write.table("clipboard", sep = "\t", row.names = FALSE, quote = FALSE)
make_vr_table(rec_yeararea, "Calves per cow [95% CI]")
make_vr_table(lam_yeararea, "Lambda [95% CI]")

fmt <- function(df) {
  data.frame(
    `Region delineation` = df$delineation,
    `Polygons (n)`       = as.character(df$n_areas),
    `Estimate [95% CI]`  = sprintf("%.3f [%.3f-%.3f]", df$estimate, df$lower, df$upper),
    `Year (%)`           = sprintf("%.1f", df$year_pct),
    `Polygon (%)`        = sprintf("%.1f", df$area_pct),
    `LRT p`              = ifelse(df$lrt_p < 0.001, "<0.001", sprintf("%.3f", df$lrt_p)),
    check.names = FALSE
  )
}

tab_stack <- rbind(
  c("Survival", rep("", 5)),
  as.matrix(fmt(surv_yeararea)),
  c("Recruitment (calves per cow)", rep("", 5)),
  as.matrix(fmt(rec_yeararea)),
  c("Lambda (\u03bb)", rep("", 5)),
  as.matrix(fmt(lam_yeararea))
)
colnames(tab_stack) <- names(fmt(surv_yeararea))

tab_stack
tab_stack |> write.table("clipboard", sep = "\t", row.names = FALSE, quote = FALSE)

# figure annual survival ----

png("./ms_varyingArea/figs/pred_survival_year_n.png",
    width = 8.5, height = 5, units = "in", res = 600)
{
  par(mar = c(4, 4, 1, 4), bty = "l", tcl = 0.3, mgp = c(2, 0.4, 0))   
  
  # strip for the bars sits below the data (0.62 to 0.70)
  ybase   <- 0.74
  ytop    <- 0.8
  scale_n <- (ytop - ybase) / 200     # bar height per individual-year (200 = top of right axis)
  
  plot(0, type = "n",
       xlim = range(pred_surv_year$bioYear) + c(-0.5, 0.5),
       ylim = c(ybase, 1.02),     
       yaxs = "i",                
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  mtext("Survival", side = 2, line = 2.5)
  
  axis(1, at = pred_surv_year$bioYear, labels = pred_surv_year$bioYear,
       las = 2, tick = FALSE)
  mtext("Biological year (Apr 01 - Mar 31)", side = 1, line = 2.8)
  
  # histogram bars
  rect(pred_surv_year$bioYear - 0.35, ybase,
       pred_surv_year$bioYear + 0.35, ybase + n_year * scale_n,
       col = adjustcolor("grey70", 0.6), border = NA)
  
  # left axis for survival only, right axis for counts
  axis(2, at = seq(0.8, 1.00, 0.05), las = 1)
  
  n_ticks <- c(0, 100, 200)
  axis(4, at = ybase + n_ticks * scale_n, labels = n_ticks, las = 1)
  mtext("Individuals", side = 4, line = 2, at = ybase + 100 * scale_n)
  
  abline(h = 0.85, lty = 2, col = "grey50")
  
  segments(pred_surv_year$bioYear, pred_surv_year$lo,
           pred_surv_year$bioYear, pred_surv_year$hi,
           col = '#228833', lwd = 1.5)
  points(pred_surv_year$bioYear, pred_surv_year$pred,
         pch = 19, col = '#228833', cex = 2)
  
  points(emp_surv_year$bioYear, emp_surv_year$mean_survived,
         pch = 22, 
         col = "black", 
         bg = alpha("white", .4),
         lwd = 1.5,
         cex = 2)
  
  legend("topright",
         inset  = c(-0.1, -0.05),
         xpd = NA,
         legend = c("Predicted (95% CI)", "Raw annual value"),
         pch    = c(21, 22),
         pt.bg  = c('#228833', "white"),
         col    = c('#228833', "black"),
         pt.cex = 2,
         pt.lwd = c(2.5, 1.5),
         lwd = c(1.5, NA),
         bty = "n")
}
dev.off()

# annual survival randomization ----
## removed from ms

n_perm <- 499
null_year_var <- numeric(n_perm)

eval = FALSE
if(eval){
  set.seed(123)  # reproducibility
  
  for (i in 1:n_perm) {
    
    perm_df <- surv_df
    
    # Shuffle bioYear labels across rows, ignoring true year structure,
    # while preserving the total number of observations per year
    # follows spatial randomization, keeps n individuals per year
    perm_df$bioYear_perm <- sample(perm_df$bioYear, replace = FALSE)
    
    m_null <- tryCatch(
      suppressWarnings(
        glmmTMB(survived ~ 1 + (1 | bioYear_perm),
                family = binomial, data = perm_df)),
      error = function(e) NULL
    )
    
    if (!is.null(m_null)) {
      vc_null <- VarCorr(m_null)
      null_year_var[i] <- as.numeric(vc_null$cond$bioYear_perm)
    } else {
      null_year_var[i] <- NA
    }
  }
  
  # no failures
  sum(is.na(null_year_var))
  
  save(null_year_var, year_var,
       file = "./ms_varyingArea/results/year_permutation.rda")
}else load(file = "./ms_varyingArea/results/year_permutation.rda")


# Permutation p-value: proportion of null year-variance values
# >= observed year variance
p_value_year <- (sum(null_year_var >= year_var) + 1) / (length(null_year_var) + 1)

year_var
p_value_year

table(surv_df$bioYear)

overall_intercept <- fixef(m_year_only)$cond[["(Intercept)"]]

pred_surv_year <- broom.mixed::tidy(m_year_only, effects = "ran_vals", conf.int = TRUE) |>
  subset(group == "bioYear") |>
  dplyr::transmute(
    bioYear = as.integer(level),
    pred    = plogis(overall_intercept + estimate),
    lo      = plogis(overall_intercept + conf.low),
    hi      = plogis(overall_intercept + conf.high)
  ) |>
  dplyr::arrange(bioYear)

n_year <- as.vector(table(surv_df$bioYear)[as.character(pred_surv_year$bioYear)])
sum(n_year)   # should equal the number of individual-years in the model (2144)


emp_surv_year <- surv_df |>
  dplyr::group_by(bioYear) |>
  dplyr::summarise(mean_survived = mean(survived), .groups = "drop")
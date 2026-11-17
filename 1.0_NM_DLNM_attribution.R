#===============================================================================
#  Unmanaged plastic waste and dengue burden in Sri Lanka
#  Analytical code: exposure–response estimation and backward attribution
#
#  Sections
#    0. Settings
#    1. Data preparation
#    2. Exposure scaling (rescaled within-district ECDF)
#    3. Stage 1: district-specific distributed lag non-linear models (DLNM)
#    4. Stage 2: multivariate random-effects meta-analysis
#    5. Pooled exposure–response curve (figure)
#    6. Backward attribution (main analysis and reference-percentile sensitivity)
#    7. Outputs and reproducibility checks
#
#  Downstream: the economic evaluation script (Part 1, E0) is run in the SAME R
#  session after this script and requires:
#    store, beta_pool, V_pool, draws, n_sim, make_ref_row, run_attribution,
#    cen_main, cen_sens, res_main, sens_all, df, outcome_from
#
#  ---------------------------------------------------------------------------
#  METHODS SUMMARY
#  ---------------------------------------------------------------------------
#  Data. Weekly notified dengue cases, district population (2024 census, held
#  constant over the study period) and the waste-specific Breteau Index
#  (BI_waste) for each district, 2023–2025. Outcomes are
#  modelled for 2024–2025 (104 weeks). 2023 is used only as a lag buffer, so
#  that each outcome week has a complete 12-week exposure history.
#
#  Exposure scaling. Entomological surveillance intensity and baseline
#  BI_waste differ across districts. Exposure was therefore expressed as a
#  within-district percentile, using the empirical cumulative distribution
#  function (ECDF) of each district's weekly series (2023–2025):
#
#      E_dt = 100 * (F_d(x_dt) - F_d(min x_d)) / (1 - F_d(min x_d))
#
#  where F_d is the district ECDF. This rescaling maps the district minimum to
#  0 and the maximum to 100. Without it, a series with many weeks at the minimum
#  (often zero) would start at F_d(min) > 0 instead of 0. Districts with a
#  constant series are excluded (rescaled value undefined).
#  Note: when many weeks are tied at the minimum, E_dt has a gap above 0. A
#  reference value of, for example, 10 on this scale is a fixed point on the
#  rescaled exposure axis. It is not the exposure value below which 10% of
#  weeks fall.
#
#  Stage 1 (district DLNMs). For each district d, a quasi-Poisson generalised
#  linear model was fitted to weekly cases in 2024–2025:
#
#      log E[Y_dt] = log(Pop_d) + cb(E_dt, lag 0–12) + ns(t, df = 8)
#
#    - Cross-basis cb: exposure dimension = natural cubic spline, knots at the
#      25th and 75th rescaled percentiles, boundary knots at 0 and 100 (3 df).
#      Lag dimension = natural cubic spline with one internal knot at lag 4
#      weeks, plus intercept (3 df). Total: 9 parameters.
#    - Long-term trend and seasonality: natural cubic spline of time with 8 df
#      over the 104 outcome weeks, fitted on 2024–2025 only.
#    - Offset: log district population (2024 census).
#    - Overdispersion: quasi-Poisson.
#  From each model we extracted (i) the full 9-parameter cross-basis
#  coefficients and covariance (used for attribution) and (ii) the overall
#  cumulative association over lags 0–12, reduced to 3 parameters with
#  reference 0 (used for the figure).
#
#  Stage 2 (meta-analysis). District coefficients were pooled with
#  multivariate random-effects meta-analysis (mixmeta, REML), separately for
#  the 3-parameter overall curve and the 9-parameter full cross-basis.
#
#  Backward attribution. For district d and week t, the attributable number is
#
#      AN_dt = n_dt * (1 - exp(-(w_dt - w_0)' beta))
#
#  where w_dt is the cross-basis row (exposure history over lags 0–12), w_0 is
#  the cross-basis row for exposure held at the reference percentile at all
#  lags, and beta is the pooled 9-parameter coefficient vector (REML mean).
#  The pooled coefficients, not district BLUPs, are applied to every district.
#  Weekly contributions are summed without truncation, so weeks below the
#  reference contribute negative AN. AF = 100 * AN / observed cases.
#  95% uncertainty intervals: 5,000 Monte Carlo draws of beta from
#  MVN(beta, V) (seed 2026). These draws are shared by all districts and
#  references. National intervals are percentiles of draw-wise national sums,
#  which preserves the correlation between districts.
#
#  Reference percentiles. Main analysis: 10th rescaled percentile (P10).
#  Sensitivity: P0, P25 and P50. P50 is reported to show why central
#  references are not used: they treat typical exposure as the counterfactual.
#===============================================================================


#===============================================================================
# 0. Settings
#===============================================================================
suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, tidyr, purrr, readr, ggplot2
  library(ISOweek)
  library(dlnm)
  library(splines)
  library(mixmeta)
})
# MASS is not attached (it masks dplyr::select); called as MASS::mvrnorm
select <- dplyr::select

# Paths: relative to the repository root
file_in <- file.path("data", "G1_weekly_series_v2_2_surinc_lag2023_v1_harmonised.csv")
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis settings
exposure_var <- "BI_waste"   # exposure column
outcome_from <- 2024             # first outcome year; 2023 = lag buffer only
lag_limit    <- c(0, 12)         # weeks
lags_to_meta <- c(0, 2, 4, 8, 10, 12)
time_df      <- 8                # df for ns(time), 2024–2025

# Cross-basis specification (used in Stage 1 AND attribution; must not differ)
argvar_cb <- list(fun = "ns", knots = c(25, 75), Boundary.knots = c(0, 100))
arglag_cb <- list(fun = "ns", knots = 4)

# Attribution settings
cen_main <- 10                   # main reference: P10
cen_sens <- c(0, 10, 25, 50)     # sensitivity references
n_sim    <- 5000                 # Monte Carlo draws
seed_mc  <- 2026


#===============================================================================
# 1. Data preparation
#===============================================================================
dat <- read_csv(file_in, show_col_types = FALSE)

# Population: one 2024 census value per district, applied to all years -------
pop_census <- dat %>%
  filter(!is.na(Population)) %>%
  group_by(District) %>%
  summarise(n_values   = n_distinct(Population),
            Population = first(Population), .groups = "drop")

stopifnot(all(pop_census$n_values == 1),
          setequal(pop_census$District, unique(dat$District)))

df <- dat %>%
  select(-Population) %>%
  left_join(select(pop_census, District, Population), by = "District") %>%
  transmute(District,
            Year       = as.integer(Year),
            Week_Num   = as.integer(Week_Num),
            Cases      = as.numeric(Cases),
            Population = as.numeric(Population),
            x          = as.numeric(.data[[exposure_var]])) %>%
  mutate(week_start = ISOweek2date(sprintf("%d-W%02d-1", Year, Week_Num))) %>%
  arrange(District, week_start) %>%
  group_by(District) %>%
  mutate(time = as.integer((week_start - min(week_start)) / 7) + 1) %>%
  ungroup()

stopifnot(!anyNA(df$Population), all(df$Population > 0))

# Data checks: continuous weekly series, no missing exposure, lag buffer present
data_check <- df %>%
  group_by(District) %>%
  summarise(n_weeks    = n(),
            gaps       = sum(diff(time) != 1),
            na_x       = sum(is.na(x)),
            n_buffer   = sum(Year < outcome_from),
            pct_at_min = round(100 * mean(x == min(x, na.rm = TRUE)), 1),
            .groups = "drop")
print(data_check, n = Inf)
stopifnot(all(data_check$gaps == 0), all(data_check$na_x == 0),
          all(data_check$n_buffer >= diff(lag_limit)))


#===============================================================================
# 2. Exposure scaling: rescaled within-district ECDF (min -> 0, max -> 100)
#===============================================================================
rescaled_ecdf <- function(x) {
  Fx    <- ecdf(x)
  F_min <- Fx(min(x))
  if (F_min >= 1) return(rep(NA_real_, length(x)))   # constant series
  100 * (Fx(x) - F_min) / (1 - F_min)
}


#===============================================================================
# 3. Stage 1: district-specific DLNMs
#===============================================================================
districts <- sort(unique(df$District))
n_dist    <- length(districts)
n_par     <- 9                                    # 3 (exposure) x 3 (lag)

coef_overall   <- matrix(NA, n_dist, 3, dimnames = list(districts, paste0("b", 1:3)))
vcov_overall   <- vector("list", n_dist)
coef_full      <- matrix(NA, n_dist, n_par, dimnames = list(districts, NULL))
vcov_full      <- vector("list", n_dist)
results_by_lag <- setNames(lapply(lags_to_meta, function(x) list()),
                           paste0("lag", lags_to_meta))
store          <- list()                          # inputs for attribution

for (i in seq_along(districts)) {
  sub <- df %>% filter(District == districts[i]) %>% arrange(time)
  sub$Exp_Pct <- rescaled_ecdf(sub$x)
  if (all(is.na(sub$Exp_Pct))) {
    message("Skipped (constant exposure): ", districts[i]); next
  }
  
  win <- sub$Year >= outcome_from
  y   <- ifelse(win, sub$Cases, NA)               # 2023 excluded from likelihood
  tb  <- matrix(NA, nrow(sub), time_df)
  tb[win, ] <- ns(sub$time[win], df = time_df)    # trend on 2024–2025 only
  
  cb <- crossbasis(sub$Exp_Pct, lag = lag_limit,
                   argvar = argvar_cb, arglag = arglag_cb)
  
  m <- glm(y ~ cb + tb + offset(log(sub$Population)),
           family = quasipoisson(), na.action = na.exclude)
  
  # Overall cumulative association (lag 0–12), reference 0: for the figure
  cr_all <- crossreduce(cb, m, type = "overall", cen = 0)
  coef_overall[i, ] <- coef(cr_all)
  vcov_overall[[i]] <- vcov(cr_all)
  
  # Lag-specific associations
  for (l in lags_to_meta) {
    cr_l <- crossreduce(cb, m, type = "lag", value = l, cen = 0)
    results_by_lag[[paste0("lag", l)]][[districts[i]]] <-
      list(coef = coef(cr_l), vcov = vcov(cr_l))
  }
  
  # Full cross-basis coefficients: for attribution
  id <- grep("^cb", names(coef(m)))
  stopifnot(length(id) == n_par)
  coef_full[i, ] <- coef(m)[id]
  vcov_full[[i]] <- vcov(m)[id, id]
  
  store[[districts[i]]] <- list(
    cb         = cb[win, , drop = FALSE],         # cross-basis, outcome weeks
    cb_full    = cb,                              # full cross-basis incl. 2023
    x_full     = sub$Exp_Pct,                     # rescaled exposure incl. 2023
    cases_full = ifelse(win, sub$Cases, NA),      # 2023 -> NA (not attributed)
    cases      = sub$Cases[win],
    year       = sub$Year[win],
    year_full  = sub$Year,
    pop        = sub$Population[win])
}


#===============================================================================
# 4. Stage 2: multivariate random-effects meta-analysis (REML)
#===============================================================================
keep <- !is.na(coef_overall[, 1])

meta_overall <- mixmeta(coef_overall[keep, ] ~ 1, S = vcov_overall[keep], method = "reml")
meta_full    <- mixmeta(coef_full[keep, ]    ~ 1, S = vcov_full[keep],    method = "reml")
summary(meta_overall)
summary(meta_full)

blups <- blup(meta_full, vcov = TRUE)             # district BLUPs (not used in AN)
names(blups) <- districts[keep]


#===============================================================================
# 5. Pooled overall exposure–response curve (lag 0–12), log-scale RR axis
#===============================================================================
b_pred  <- onebasis(0:100, fun = "ns", knots = argvar_cb$knots,
                    Boundary.knots = argvar_cb$Boundary.knots)
pred_ov <- crosspred(b_pred, coef = coef(meta_overall), vcov = vcov(meta_overall),
                     model.link = "log", cen = 0, by = 1)

# District-specific Stage-1 curves (same basis and reference)
regall <- apply(coef_overall[keep, ], 1, function(b) exp(b_pred %*% b))

show_districts <- FALSE                           # TRUE to overlay district curves
y_lim  <- c(0.5, 10)                              # log axis: lower limit must be > 0
y_tick <- c(0.5, 1, 2, 5, 10)

plot_overall <- function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  plot(pred_ov$predvar, pred_ov$allRRfit, type = "n", log = "y",
       ylim = y_lim, yaxt = "n", las = 1, bty = "l",
       xlab = expression(BI[waste] ~ "(district-specific rescaled percentile)"),
       ylab = "Relative risk (log scale)",
       main = "Pooled overall cumulative association, lag 0\u201312 weeks")
  axis(2, at = y_tick, labels = y_tick, las = 1)
  if (show_districts)
    matplot(0:100, regall, type = "l", col = grey(0.6), lty = 2, lwd = 0.8, add = TRUE)
  polygon(c(pred_ov$predvar, rev(pred_ov$predvar)),
          c(pred_ov$allRRlow, rev(pred_ov$allRRhigh)),
          col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)
  abline(h = 1, lty = 2, col = "grey40")
  lines(pred_ov$predvar, pred_ov$allRRfit, col = "darkblue", lwd = 3)
  leg <- c("Pooled estimate", "95% CI", if (show_districts) "District-specific")
  legend("topleft", legend = leg,
         col  = c("darkblue", NA, grey(0.6))[seq_along(leg)],
         lwd  = c(3, NA, 0.8)[seq_along(leg)],
         lty  = c(1, NA, 2)[seq_along(leg)],
         fill = c(NA, rgb(0.1, 0.1, 0.8, 0.2), NA)[seq_along(leg)],
         border = NA, bty = "n", cex = 0.9)
}

plot_overall()
pdf(file.path(out_dir, "Overall_RR_lag0_12.pdf"), 7, 6); plot_overall(); dev.off()

# Report values outside the plotting range (clipped)
if (any(pred_ov$allRRlow < y_lim[1] | pred_ov$allRRhigh > y_lim[2]))
  message("Pooled CI extends beyond y_lim; consider widening y_lim.")
if (show_districts) {
  clipped <- colnames(regall)[apply(regall, 2, function(r) any(r > y_lim[2] | r < y_lim[1]))]
  if (length(clipped)) message("District curves clipped: ", paste(clipped, collapse = ", "))
}


#===============================================================================
# 6. Backward attribution
#===============================================================================

# 6.1 Pooled full cross-basis coefficients (REML mean) ------------------------
beta_pool <- coef(meta_full)
V_pool    <- vcov(meta_full)

# 6.2 Monte Carlo draws: drawn ONCE, shared by all districts and references ---
set.seed(seed_mc)
draws <- MASS::mvrnorm(n_sim, mu = beta_pool, Sigma = V_pool)   # n_sim x 9

# 6.3 Reference cross-basis row: exposure held at `cen` over all lags (w_0) ---
make_ref_row <- function(cen) {
  cb_ref <- crossbasis(rep(cen, diff(lag_limit) + 1), lag = lag_limit,
                       argvar = argvar_cb, arglag = arglag_cb)
  as.numeric(cb_ref[nrow(cb_ref), ])
}

# 6.4 One district: AN_dt = n_dt * (1 - exp(-(w_dt - w_0)' beta)) -------------
attr_one_district <- function(s, ref_row) {
  cbc <- sweep(unclass(s$cb_full), 2, ref_row)               # w_dt - w_0
  ok  <- !is.na(s$cases_full) & complete.cases(cbc)          # 2024–25, full lag history
  X   <- cbc[ok, , drop = FALSE]
  n   <- s$cases_full[ok]
  
  eta_hat <- drop(X %*% beta_pool)
  eta_sim <- X %*% t(draws)                                  # n_weeks x n_sim
  
  list(n_weeks = sum(ok),
       cases   = sum(n),
       AN      = sum(n * (1 - exp(-eta_hat))),               # signed, not truncated
       AN_sim  = colSums(n * (1 - exp(-eta_sim))))
}

# 6.5 Attribution engine for one reference percentile -------------------------
run_attribution <- function(cen) {
  ref_row <- make_ref_row(cen)
  res     <- map(store, attr_one_district, ref_row = ref_row)
  
  tab_dist <- imap_dfr(res, function(a, d) tibble(
    District = d, n_weeks = a$n_weeks, cases = a$cases, AN = a$AN,
    AN_low   = unname(quantile(a$AN_sim, 0.025)),
    AN_high  = unname(quantile(a$AN_sim, 0.975))))
  
  nat_sim <- Reduce(`+`, map(res, "AN_sim"))                 # draw-wise national sum
  tab_nat <- tibble(
    District = "National", n_weeks = NA_integer_,
    cases = sum(tab_dist$cases), AN = sum(tab_dist$AN),
    AN_low  = unname(quantile(nat_sim, 0.025)),
    AN_high = unname(quantile(nat_sim, 0.975)))
  
  table <- bind_rows(tab_dist, tab_nat) %>%
    mutate(AF      = 100 * AN      / cases,
           AF_low  = 100 * AN_low  / cases,
           AF_high = 100 * AN_high / cases,
           cen     = cen)
  
  list(table = table, nat_sim = nat_sim, district_res = res)
}

# 6.6 Main analysis (P10) ------------------------------------------------------
res_main <- run_attribution(cen_main)

attr_table <- res_main$table %>%
  mutate(across(c(AN, AN_low, AN_high), round),
         across(c(AF, AF_low, AF_high), ~ round(.x, 1)))
print(attr_table, n = Inf)

# 6.7 Checks -------------------------------------------------------------------
# (a) every district uses all 104 outcome weeks
stopifnot(all(attr_table$n_weeks[attr_table$District != "National"] == 104))
# (b) districts with negative AF (exposure mostly below reference)
print(filter(attr_table, District != "National", AF < 0))
# (c) full pooled coefficients reproduce the reduced overall curve
cb_chk  <- crossbasis(rep(cen_main, diff(lag_limit) + 1), lag = lag_limit,
                      argvar = argvar_cb, arglag = arglag_cb)
cp_full <- crosspred(cb_chk, coef = beta_pool, vcov = V_pool,
                     model.link = "log", at = 0:100, cen = cen_main)
ref_ov  <- approx(pred_ov$predvar, pred_ov$allRRfit, xout = cen_main)$y
plot(cp_full, "overall", log = "y",
     xlab = "Exposure (rescaled percentile)", ylab = "RR (log scale)",
     main = paste0("Check: full vs reduced pooled curve, reference P", cen_main))
lines(pred_ov$predvar, pred_ov$allRRfit / ref_ov, col = "red", lty = 2)
legend("topleft", c("Full model", "Reduced overall (re-centred)"),
       col = c("black", "red"), lty = c(1, 2), bty = "n")

# 6.8 Sensitivity: alternative reference percentiles ---------------------------
sens_all <- map(cen_sens, run_attribution)
names(sens_all) <- paste0("p", cen_sens)

sens_nat <- map_dfr(sens_all, ~ filter(.x$table, District == "National")) %>%
  transmute(Reference = paste0("P", cen),
            Main      = cen == cen_main,
            AN = round(AN), AN_low = round(AN_low), AN_high = round(AN_high),
            AF = round(AF, 1), AF_low = round(AF_low, 1), AF_high = round(AF_high, 1),
            Reported  = sprintf("%.1f%% (95%% UI: %.1f\u2013%.1f); %s cases",
                                AF, AF_low, AF_high, format(AN, big.mark = ",")))
print(sens_nat)

sens_dist <- map_dfr(sens_all, "table") %>%
  mutate(est = sprintf("%.1f (%.1f to %.1f)", AF, AF_low, AF_high),
         Reference = paste0("P", cen)) %>%
  select(District, Reference, est) %>%
  pivot_wider(names_from = Reference, values_from = est)
print(sens_dist, n = Inf)


#===============================================================================
# 7. Outputs and reproducibility
#===============================================================================
write_csv(attr_table, file.path(out_dir, "attribution_P10_district_national.csv"))
write_csv(sens_nat,   file.path(out_dir, "attribution_sensitivity_national.csv"))
write_csv(sens_dist,  file.path(out_dir, "attribution_sensitivity_district.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# Objects required downstream (economic evaluation, Part 1 E0)
required_upstream <- c("store", "beta_pool", "V_pool", "draws", "n_sim",
                       "make_ref_row", "run_attribution", "cen_main", "cen_sens",
                       "res_main", "sens_all", "df", "outcome_from")
stopifnot(all(vapply(required_upstream, exists, logical(1))),
          nrow(draws) == n_sim, ncol(draws) == length(beta_pool),
          length(res_main$nat_sim) == n_sim)

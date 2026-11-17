

#===============================================================================
# NATIONAL ECONOMIC EVALUATION OF DENGUE ATTRIBUTABLE TO BI_waste
# Part 1: Link to attribution (Section 5) + economic and demographic inputs (B0)
#
# Run order: in the SAME R session, immediately after the attribution script
#            (Sections 0-5). No objects are re-estimated here.
# Price year: 2024 USD. Perspective: societal (programme cost reported separately).
# Uncertainty: propagated through the pooled Monte Carlo draws created in 5.2
#              (set.seed(2026); n_sim = 5000), shared by all components.
#===============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(readr)
})
select <- dplyr::select              # MASS (attached in Section 0) masks select

#-------------------------------------------------------------------------------
# E0. Link to attribution: required objects from Sections 3-5
#-------------------------------------------------------------------------------
#  store          list by district: cb_full, x_full, cases_full, cases, year, pop (Sec 3)
#  beta_pool      pooled 9-parameter cross-basis coefficients (REML)       (5.1)
#  V_pool         covariance of beta_pool                                  (5.1)
#  draws          n_sim x 9 Monte Carlo draws of beta_pool                 (5.2)
#  n_sim          number of draws (5000)                                   (5.0)
#  make_ref_row   reference cross-basis row at a given percentile          (5.3)
#  run_attribution attribution engine (district + draw-wise national)      (5.5)
#  cen_main       main reference percentile (10)                           (5.0)
#  cen_sens       reference scenarios reported (0, 10, 25, 50)             (5.0)
#  res_main       attribution at P10: table, nat_sim, district_res         (5.6)
#  sens_all       attribution at each reference in cen_sens                (5.8)
#  df             analysis data (Population = 2024 census)                 (Sec 1)
#  outcome_from   first outcome year (2024; 2023 = lag buffer)             (Sec 0)
required_upstream <- c("store", "beta_pool", "V_pool", "draws", "n_sim",
                       "make_ref_row", "run_attribution", "cen_main", "cen_sens",
                       "res_main", "sens_all", "df", "outcome_from")
missing_upstream <- required_upstream[!vapply(required_upstream, exists, logical(1))]
if (length(missing_upstream))
  stop("Run the attribution script (Sections 0-5) first. Missing: ",
       paste(missing_upstream, collapse = ", "))

# Draws must be the n_sim x 9 matrix from 5.2 (same draws used for the AN UIs)
draws <- as.matrix(draws)
if (ncol(draws) != length(beta_pool) && nrow(draws) == length(beta_pool)) draws <- t(draws)
if (nrow(draws) != n_sim || ncol(draws) != length(beta_pool))
  stop("`draws` must be ", n_sim, " x ", length(beta_pool),
       "; found ", nrow(draws), " x ", ncol(draws), ". Re-run Section 5.2.")
if (length(res_main$nat_sim) != n_sim)
  stop("res_main$nat_sim length differs from n_sim; re-run Section 5.6.")

#===============================================================================
# B0. Economic and demographic inputs
#   Productivity: human capital approach at observed earnings (LFS 2024)
#   Conversion (Turner et al., Value Health 2019;22:1026-32):
#     study-year USD -> LKR at study-year rate -> Sri Lanka CPI -> 2024 USD
#===============================================================================

#--- Age groups (mortality profiles; YLL computed at age-group mid-points) -----
age_groups <- c("< 1","1 - 4","5 - 9","10 - 14","15 - 19","20 - 24",
                "25 - 34","35 - 44","45 - 54","55 - 64","≥ 65")
age_mids   <- c(0.5, 2.5, 7, 12, 17, 22, 29.5, 39.5, 49.5, 59.5, 70)

#--- Observed dengue deaths by age. Source: Epidemiology Unit / NDCU, 2024-2025
study_demographics_2024 <- tibble(age_group = age_groups, age_mid = age_mids,
                                  observed_deaths = c(0, 1, 1, 0, 2, 3, 5, 6, 3, 2, 1))
study_demographics_2025 <- tibble(age_group = age_groups, age_mid = age_mids,
                                  observed_deaths = c(0, 0, 2, 0, 0, 3, 7, 7, 3, 3, 4))
mortality_profiles_by_year <- list("2024" = study_demographics_2024,
                                   "2025" = study_demographics_2025)
observed_deaths_lookup <- tibble(Year = c(2024, 2025), Observed_Deaths = c(24, 29))
stopifnot(all(map_dbl(mortality_profiles_by_year, ~ sum(.x$observed_deaths)) ==
                observed_deaths_lookup$Observed_Deaths))

#--- Hospitalisation and severity mix -------------------------------------------
# Source: NDCU NADSys surveillance database and NDCU hospital-based evaluations
# (unpublished).
p_hosp_primary <- 0.75               # proportion of reported cases hospitalised (main)
p_hosp_sa      <- 1.00               # reported scenario

prop_paed <- 0.203; prop_adult <- 0.797                          # age split, hospitalised
sev_paed_df  <- 0.650; sev_paed_dhf  <- 0.257; sev_paed_dss  <- 0.093
sev_adult_df <- 0.736; sev_adult_dhf <- 0.175; sev_adult_dss <- 0.089
stopifnot(abs(prop_paed + prop_adult - 1) < 1e-9,
          abs(sev_paed_df  + sev_paed_dhf  + sev_paed_dss  - 1) < 1e-9,
          abs(sev_adult_df + sev_adult_dhf + sev_adult_dss - 1) < 1e-9)

#--- Currency conversion to 2024 USD --------------------------------------------
# CPI: World Bank WDI FP.CPI.TOTL, Sri Lanka (2010 = 100)
# FX:  WDI PA.NUS.FCRF (LKR per USD, annual average)
cpi     <- c("2012" = 114.77, "2013" = 122.69, "2018" = 150.23,
             "2024" = 306.96, "2025" = 292.34)
fx_2013 <- 129.07
fx_2024 <- 302.12
fx_2025 <- 300.85  # Reference only, NOT used for conversion (CBSL reference rate; VERIFY primary source)
# factor = FX_study * (CPI_2024 / CPI_study) / FX_2024
conv_factor <- function(fx_study, cpi_from, cpi_to = cpi[["2024"]], fx_to = fx_2024)
  unname(fx_study * (cpi_to / cpi_from) / fx_to)
fx_18_24 <- conv_factor(1 / 0.0067, cpi[["2018"]])   # 2018 USD -> 2024 USD
fx_12_24 <- conv_factor(127.60,     cpi[["2012"]])   # 2012 USD -> 2024 USD
fx_13_24 <- conv_factor(fx_2013,    cpi[["2013"]])   # 2013 USD -> 2024 USD
lkr25_to_usd24 <- function(lkr) lkr * (cpi[["2024"]] / cpi[["2025"]]) / fx_2024
stopifnot(abs(fx_18_24 - 1.0094) < 1e-3, abs(fx_12_24 - 1.1297) < 1e-3,
          abs(fx_13_24 - 1.0688) < 1e-3)

#--- Provider cost per case, public sector (2024 USD) ---------------------------
# Source: Thalagala et al. 2016 (2012 USD). p = paediatric, a = adult
c_ward_p_df_mean  <- 215 * fx_12_24; c_ward_p_dhf_mean <- 293 * fx_12_24
c_icu_p_dhf_mean  <- 609 * fx_12_24; c_ward_a_df_mean  <- 196 * fx_12_24
c_ward_a_dhf_mean <- 255 * fx_12_24; c_icu_a_dhf_mean  <- 866 * fx_12_24
c_amb_min <- 39.20 * fx_12_24; c_amb_max <- 53.75 * fx_12_24
c_ambulatory_mean <- (c_amb_min + c_amb_max) / 2     # ambulatory: mid-point of range
c_amb_sa_shepard  <- 84 * fx_13_24                   # Shepard et al. 2016 (2013 USD)
provider_cv       <- 0.20                            # CV for provider-cost sampling

# Shepard et al. 2016 national unit costs (USD, as in source)
shepard_unit_costs_usd <- c(inpatient_provider = 307, inpatient_oop = 51,
                            inpatient_indirect = 6,   ambulatory_provider = 52,
                            ambulatory_oop = 29,      ambulatory_indirect = 29,
                            indirect_nonmed = 29,
                            death_adult = 95000,      death_child = 62000)

#--- Out-of-pocket cost per inpatient episode (2024 USD) ------------------------
# Source: Weerasinghe et al., BMC Health Serv Res 2022;22:657 (2018 USD)
oop_m  <- 16.29 * fx_18_24; oop_sd <- 14.02 * fx_18_24; n_src_oop <- 833
oop_amb_ratio <- 5.06 / 16.29        # ambulatory OOP = pre-admission share of inpatient OOP

#--- Productivity: days lost per episode. Source: Weerasinghe et al. 2022 --------
# Sampled as the mean: SE = SD / sqrt(n)
days_pt_m <- 21.51; days_pt_sd <- 41.71; n_days_pt <- 524   # patient work days
days_cg_m <- 2.91;  days_cg_sd <- 2.47;  n_days_cg <- 833   # caregiver days
p_emp     <- 524 / 994               # adults reporting missed work / adult respondents
ind_amb_ratio <- 48 / 85             # ambulatory : inpatient indirect cost (Shepard 2016)


#===============================================================================
# B1. Earnings: one source for both morbidity and mortality
#   Source: DCS Sri Lanka Labour Force Survey (LFS) Annual Report 2024, Table 4.6
#   Daily-paid workers stand in for informal workers, monthly-paid for formal.
#   Self-employed workers are valued at employee earnings (assumption).
#===============================================================================
earn_m_daily   <- 33523   # mean monthly earnings, daily-paid workers (LKR)
earn_m_monthly <- 58292   # mean monthly earnings, monthly-paid workers (LKR)
p_informal     <- 0.665   # informal employment share (LFS 2024, ILO definition)

# Weighted monthly earnings -> annual USD -> daily wage
earn_month_lkr  <- p_informal * earn_m_daily + (1 - p_informal) * earn_m_monthly  # 41,821 LKR
annual_earn_usd <- earn_month_lkr * 12 / fx_2024                                  # 1,661 USD/yr
work_days_yr    <- 260                                 # 5-day week (scenario: 312)
wage_usd        <- annual_earn_usd / work_days_yr      # 6.39 USD/day
wage_cv         <- 0.10                                # CV of shared wage multiplier

#--- Wage alternatives (USD/day) -------------------------------------------------
wi_2018 <- 100; wi_2024 <- 181.1   # CBSL Informal Private Sector Wage Rate Index (2018 = 100)
wage_usd_ipswri <- 2037 * wi_2024 / wi_2018 / fx_2024  # 12.21 upper: Weerasinghe wage x index
wage_usd_lfsdly <- (earn_m_daily / 26) / fx_2024       # 4.27 lower: daily-paid, 26 days/month
wage_usd_min    <- 700 / fx_2024                       # 2.32 floor: Minimum Wage Act No. 48/2024
wage_usd_cpi    <- 13.65 * fx_18_24                    # 13.78: 2018 wage, CPI-inflated only
real_wage_r     <- wage_usd_ipswri / wage_usd_cpi      # ~0.886: real wage change since 2018
stopifnot(abs(earn_month_lkr - 41821) < 1, abs(wage_usd - 6.39) < 0.01,
          abs(wage_usd_ipswri - 12.21) < 0.01, abs(real_wage_r - 0.886) < 0.005)

#===============================================================================
# B2. Valuing deaths: lost earnings over remaining life expectancy
#   Each year of life lost counts as earnings x employment-to-population ratio
#   at that age. No fixed retirement age; 0% real earnings growth.
#   Source: LFS 2024, Table 4.11 (both sexes)
#   Ages 55-64 and 65+: labour force participation (50.7%, 20.2%) x the
#   employed/labour-force ratio at ages 25-54 (64.2 / 66.4)
#===============================================================================
el_ratio_2554 <- 0.642 / 0.664
ep_ratio <- function(age) dplyr::case_when(
  age < 15 ~ 0,
  age < 25 ~ 0.185,
  age < 55 ~ 0.642,
  age < 65 ~ 0.507 * el_ratio_2554,     # 0.490
  TRUE     ~ 0.202 * el_ratio_2554)     # 0.195
stopifnot(abs(ep_ratio(70) - 0.195) < 0.001, abs(ep_ratio(60) - 0.490) < 0.001)

# Labour share of GDP, for the GDP-based alternative only
lab_share_lo <- 0.285   # ILO modelled estimate 2024
lab_share_hi <- 0.597   # DCS 2024

#===============================================================================
# B3. Disability weights and illness duration (for YLD)
#   YLD per case = disability weight x duration in years
#===============================================================================
daly_ambulatory_mean   <- 0.032; daly_ambulatory_sd   <- 0.006
daly_hospitalized_mean <- 0.036; daly_hospitalized_sd <- 0.007
dengue_duration_years  <- 14 / 365.25   # 14-day episode

#===============================================================================
# B4. GDP per capita (current USD), for the GDP-based alternative only
#   Source: World Bank WDI NY.GDP.PCAP.CD
#===============================================================================
gdp_lookup <- tibble(Year = c(2024, 2025), GDP_Capita_USD = c(4545.36, 5002.08))
gdp_cv     <- 0.10

#===============================================================================
# B5. Paediatric inpatient indirect cost (caregiver time)
#   Source: Fernando et al. 2021, Table 5 (2013 USD): DF $11, DHF/DSS $24,
#   overall $15. Converted to 2024 USD with the same factor as other 2013 inputs.
#   Review-reported values ($15/$33) are already price-adjusted and are not used,
#   to avoid converting twice.
#===============================================================================
ind_p_conversion   <- fx_13_24                       # 1.0688
ind_p_df_2024      <- 11.0 * ind_p_conversion        # 11.76
ind_p_dhf_2024     <- 24.0 * ind_p_conversion        # 25.65
ind_p_dss_2024     <- ind_p_dhf_2024                 # DHF and DSS reported together
ind_p_overall_2024 <- 15.0 * ind_p_conversion        # 16.03
ind_p_min <- ind_p_df_2024; ind_p_max <- ind_p_dhf_2024   # range if severity not split
stopifnot(abs(ind_p_conversion - 1.0688) < 1e-3,
          abs(ind_p_df_2024 - 11.7574) < 1e-3,
          abs(ind_p_dhf_2024 - 25.6525) < 1e-3,
          abs(ind_p_overall_2024 - 16.0328) < 1e-3)

#===============================================================================
# B6. Vector-control programme cost (reported separately from societal cost)
#   Source: NDCU / MOH expenditure records (LKR million)
#   2024: LKR / FX 2024. 2025: LKR 2025 -> 2024 prices via CPI, then / FX 2024.
#===============================================================================
programme_lookup <- tibble(
  Year              = c(2024, 2025),
  Operational_LKR_M = c(45.68 + 9.45 + 5.15 * 0.75, 67.01 + 18.90 + 1.22 * 0.75),
  Labor_LKR_M       = c((21176298 + 91466737 + 286853355) / 1e6,
                        (22290840 + 96280776 + 301950900) / 1e6)) %>%
  mutate(Total_Vector_Control_USD_M = if_else(
    Year == 2025,
    lkr25_to_usd24((Operational_LKR_M + Labor_LKR_M) * 1e6) / 1e6,
    (Operational_LKR_M + Labor_LKR_M) / fx_2024))

message(sprintf("Earnings: LKR %.0f/month | USD %.0f/yr | USD %.2f/day (real-wage r = %.3f)",
                earn_month_lkr, annual_earn_usd, wage_usd, real_wage_r))

#===============================================================================
# B7. Reference diagnostics: what the P10 reference means in raw BI per district
#   pct_tied_min           % of weeks at the district minimum (percentile 0)
#   BI_at_P10              raw BI_waste value at the 10th percentile
#   pct_wk_below_P10_2425  % of 2024-25 weeks below P10 (these push AN negative)
#===============================================================================
ref_diag <- imap_dfr(store, function(s, d) {
  x_raw <- df$x[df$District == d][order(df$time[df$District == d])]
  tibble(District  = d,
         pct_tied_min = round(100 * mean(s$x_full == 0), 1),
         BI_min       = min(x_raw),
         BI_at_P10    = min(x_raw[s$x_full >= 10]),
         BI_median    = median(x_raw),
         pct_wk_below_P10_2425 = round(100 * mean(s$x_full[s$year_full >= outcome_from] < 10), 1))
})
print(ref_diag, n = Inf)

#===============================================================================
# B8. Attributable cases by district, year and Monte Carlo draw
#   Same method as Section 5, split by year (2024, 2025). Each draw uses the
#   same pooled coefficients as Section 5 (draws, n_sim = 5000).
#   AN       = signed: weeks below the reference count as negative
#   AN_trunc = weekly AF floored at 0 (non-negative alternative)
#===============================================================================
stopifnot(is.matrix(draws), ncol(draws) == length(beta_pool))
n_sim     <- nrow(draws)
out_years <- c(2024, 2025)

# One district: weekly AF -> annual totals (point estimate and per draw)
attr_one_district_yr <- function(s, ref_row) {
  cbc <- sweep(unclass(s$cb_full), 2, ref_row)               # w_t - w_ref
  ok  <- !is.na(s$cases_full) & complete.cases(cbc)          # 2024-25 weeks only
  X  <- cbc[ok, , drop = FALSE]; n <- s$cases_full[ok]; yr <- s$year_full[ok]
  af_hat <- 1 - exp(-drop(X %*% beta_pool))                  # weeks
  af_sim <- 1 - exp(-(X %*% t(draws)))                       # weeks x n_sim
  list(cases        = rowsum(n, yr)[, 1],
       AN_hat       = rowsum(n * af_hat, yr)[, 1],
       AN_sim       = rowsum(n * af_sim, yr),                # years x n_sim
       AN_sim_trunc = rowsum(n * pmax(af_sim, 0), yr))
}

# All districts, long format: District x Year x sim
build_an_draws <- function(cen) {
  ref_row <- make_ref_row(cen)
  out <- imap_dfr(store, function(s, d) {
    r <- attr_one_district_yr(s, ref_row)
    map_dfr(rownames(r$AN_sim), function(y) tibble(
      District = d, Year = as.integer(y), sim = seq_len(n_sim),
      Cases  = unname(r$cases[y]), AN_hat = unname(r$AN_hat[y]),
      AN     = r$AN_sim[y, ],      AN_trunc = r$AN_sim_trunc[y, ]))
  }) %>% mutate(cen = cen)
  stopifnot(!anyNA(out$AN), all(out_years %in% out$Year),
            all(count(out, District, Year)$n == n_sim))
  out
}

an_draws <- build_an_draws(cen_main)                         # main: P10

# Attributable cases for any reference; af = "signed" (main) or "trunc"
national_an <- function(ref = cen_main, af = c("signed", "trunc")) {
  af <- match.arg(af)
  z  <- if (identical(as.numeric(ref), as.numeric(cen_main))) an_draws else build_an_draws(ref)
  z %>% mutate(AN_selected = if (af == "signed") AN else AN_trunc,
               af_type = af, reference_percentile = ref)
}

#--- Consistency check: must reproduce Section 5 exactly ------------------------
# (a) district point estimates, (b) national draw-by-draw totals
gate <- an_draws %>% filter(sim == 1) %>%
  group_by(District) %>% summarise(AN_bridge = sum(AN_hat), .groups = "drop") %>%
  left_join(res_main$table %>% select(District, AN_sec5 = AN), by = "District")
stopifnot(!anyNA(gate$AN_sec5), max(abs(gate$AN_bridge - gate$AN_sec5)) < 1e-6)

nat_chk <- an_draws %>% group_by(sim) %>% summarise(AN = sum(AN), .groups = "drop")
stopifnot(isTRUE(all.equal(nat_chk$AN, unname(res_main$nat_sim))))
message("Check passed: annual split reproduces Section 5 estimates and draws.")

#===============================================================================
# B9. Parameter draws: created once, shared by all analyses
#   Row i of P is paired with attribution draw i.
#   Costs: gamma (provider), log-normal (OOP), uniform (ambulatory range).
#   Days lost: gamma around the MEAN (SE = SD / sqrt(n)), not individual spread.
#   Wage and GDP: multipliers with mean 1; levels are applied later as scenarios.
#===============================================================================
set.seed(20261104)

# Log-normal for a mean cost, SE = SD / sqrt(n)
draw_mean_cost <- function(m, sd, n_src) {
  se <- sd / sqrt(n_src); s <- sqrt(log(1 + se^2 / m^2))
  rlnorm(n_sim, log(m) - s^2 / 2, s)
}
gam    <- function(m, cv) rgamma(n_sim, shape = 1 / cv^2, scale = m * cv^2)  # mean m, CV cv
gam_se <- function(m, sd, n) gam(m, (sd / sqrt(n)) / m)                      # gamma, SE of mean

# Severity mix (DF, DHF, DSS). 'mult' scales the combined DHF+DSS share
# (severity-mix scenario) and keeps the DHF:DSS ratio. The DF argument is not
# used; DF is calculated as 1 - severe.
severity_weights <- function(df, dhf, dss, mult = 1) {
  severe <- min(1, max(0, (dhf + dss) * mult))
  c(1 - severe, severe * dhf / (dhf + dss), severe * dss / (dhf + dss))
}

P <- tibble(
  c_amb     = runif(n_sim, c_amb_min, c_amb_max),          # ambulatory provider cost
  c_p_df    = gam(c_ward_p_df_mean,  provider_cv),         # paediatric: ward DF
  c_p_dhf   = gam(c_ward_p_dhf_mean, provider_cv),         #             ward DHF
  c_p_dss   = gam(c_icu_p_dhf_mean,  provider_cv),         #             ICU (DSS)
  c_a_df    = gam(c_ward_a_df_mean,  provider_cv),         # adult: ward DF
  c_a_dhf   = gam(c_ward_a_dhf_mean, provider_cv),         #        ward DHF
  c_a_dss   = gam(c_icu_a_dhf_mean,  provider_cv),         #        ICU (DSS)
  oop       = draw_mean_cost(oop_m, oop_sd, n_src_oop),    # inpatient OOP
  days_pt   = gam_se(days_pt_m, days_pt_sd, n_days_pt),    # patient work days lost
  days_cg   = gam_se(days_cg_m, days_cg_sd, n_days_cg),    # caregiver days lost
  ind_p     = rep(ind_p_overall_2024, n_sim),              # paediatric caregiver cost: fixed
  ind_p_df  = rep(ind_p_df_2024,  n_sim),                  #   (no SD from the study used)
  ind_p_dhf = rep(ind_p_dhf_2024, n_sim),
  ind_p_dss = rep(ind_p_dss_2024, n_sim),
  dw_amb    = pmax(rnorm(n_sim, daly_ambulatory_mean,   daly_ambulatory_sd),   0),
  dw_hosp   = pmax(rnorm(n_sim, daly_hospitalized_mean, daly_hospitalized_sd), 0),
  wage_mult = gam(1, wage_cv),                             # shared: morbidity + mortality
  gdp_mult  = gam(1, gdp_cv))

# Draws are complete, non-negative and centred on the input means
stopifnot(nrow(P) == n_sim, !anyNA(P), all(P >= 0),
          abs(mean(P$c_amb) / c_ambulatory_mean - 1) < 0.02,
          abs(mean(P$oop) / oop_m - 1) < 0.02,
          abs(mean(P$wage_mult) - 1) < 0.02,
          abs(mean(P$gdp_mult) - 1) < 0.02)


#===============================================================================
# B10. Life table and per-death constants
#   One death at age a is valued three ways:
#     yll_per_death  years of life lost (remaining life expectancy)
#     pv_work        present value of lost working years, weighted by the
#                    employment rate at each age (x annual earnings later)
#     pv_life        present value of all remaining years, unweighted
#                    (x GDP per capita later)
#   Each value is averaged over the observed age distribution of deaths.
#===============================================================================

#--- Step 1: reference life table ----------------------------------------------
# GBD reference life table: remaining life expectancy at exact age
# (to confirm against the GBD 2021 reference table)
gbd_le <- tibble(age = c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85),
                 le  = c(88.9, 88.0, 84.0, 79.0, 74.1, 69.1, 64.1, 59.2, 54.2, 49.3,
                         44.4, 39.6, 34.8, 30.1, 25.5, 21.1, 17.0, 13.1, 9.7))

# Linear interpolation at the age-group midpoints (e.g. 29.5, 70)
le_at <- function(a) approx(gbd_le$age, gbd_le$le, xout = a, rule = 2)$y

#--- Step 2: discounting --------------------------------------------------------
base_dr <- 0          # Main analysis: 0% discount rate for costs and YLL
sens_dr <- 0.03       # Sensitivity analysis: 3% applied jointly to costs and YLL

# Methodological references:
# Drummond MF, Sculpher MJ, Claxton K, Stoddart GL, Torrance GW.
# Methods for the Economic Evaluation of Health Care Programmes. 4th ed.
# Oxford University Press; 2015. Chapter 7: discounting.
# Publisher/book information:
# https://global.oup.com/academic/product/methods-for-the-economic-evaluation-of-health-care-programmes-9780199665884
# Book contents showing Chapter 7:
# https://media.repository.chds.harvard.edu/static/filer_public/b0/03/b003eca5-498b-49c2-8ebd-8d1587a78111/2015_drummond_toc_ch.pdf
#
# Tan-Torres Edejer T, Baltussen R, Adam T, et al., eds.
# Making Choices in Health: WHO Guide to Cost-Effectiveness Analysis.
# World Health Organization; 2003. See section 5, “Discounting.”
# WHO IRIS record and document:
# https://iris.who.int/handle/10665/42699
# Direct PDF:
# https://iris.who.int/bitstream/handle/10665/42699/9241546018.pdf


# Present value of t years at rate r (continuous annuity form; r = 0 -> t)
ann <- function(t, r) {
  stopifnot(length(r) == 1)
  if (r == 0) t else (1 - (1 + r)^(-t)) / r
}

#--- Step 3: present value of years lost for one death at age a ----------------
#   weight = "emp"  : each year weighted by the employment rate at that age
#   weight = "none" : every year counts fully
#   stop_age = Inf  : horizon = remaining life expectancy
#   stop_age = 60/65: horizon = years until retirement
#   The last part-year counts in proportion (frac).
pv_years <- function(a, le_rem, r, weight = c("emp", "none"), stop_age = Inf) {
  weight  <- match.arg(weight)
  horizon <- if (is.finite(stop_age)) max(0, stop_age - a) else le_rem
  if (horizon <= 0) return(0)
  t    <- seq.int(0, ceiling(horizon) - 1)          # years after death
  frac <- pmin(1, horizon - t)                      # last part-year
  w    <- if (weight == "emp") ep_ratio(a + t) else rep(1, length(t))
  sum(frac * w / (1 + r)^t)
}

# Checks: 10 years undiscounted = 10; discounting lowers the value
stopifnot(abs(pv_years(30, 10, 0, "none") - 10) < 1e-9,
          pv_years(30, 50, 0.03, "emp") < pv_years(30, 50, 0, "emp"))

#--- Step 4: average per death over the age distribution of deaths -------------
#   r_yll   discount rate for YLL (GBD main analysis: 0)
#   r_cost  discount rate for earnings
#   work_hz "profile" (to life expectancy), "stop60", "stop65"
age_consts <- function(prof, r_yll, r_cost, work_hz) {
  share    <- prof$observed_deaths / sum(prof$observed_deaths)
  a        <- prof$age_mid
  le       <- le_at(a)
  stop_age <- switch(work_hz, profile = Inf, stop60 = 60, stop65 = 65)
  c(yll_per_death = sum(share * ann(le, r_yll)),
    pv_work = sum(share * mapply(pv_years, a, le,
                                 MoreArgs = list(r = r_cost, weight = "emp",  stop_age = stop_age))),
    pv_life = sum(share * mapply(pv_years, a, le,
                                 MoreArgs = list(r = r_cost, weight = "none", stop_age = stop_age))))
}

#--- Step 5: per-death constants for one year and scenario ---------------------
get_kc <- function(yr, retire = "stop60", dr = base_dr, r_yll = 0) {
  stopifnot(retire %in% c("profile", "stop60", "stop65"), dr >= 0, r_yll >= 0)
  age_consts(mortality_profiles_by_year[[as.character(yr)]], r_yll, dr, retire)
}

# Quick look at the main constants
map_dfr(c(2024, 2025), ~ c(Year = .x, get_kc(.x))) %>% print()


#===============================================================================
# B11. Scenarios
#   Main grid = structural choices, all combinations:
#     work_hz   retirement horizon: stop60 (primary), profile (no retirement), stop65
#     p_hosp    share of cases hospitalised: primary, sensitivity
#     af_type   signed AN (primary) or AN with weekly AF floored at 0
#     discount  paired: 0% costs + 0% YLL (primary); 3% costs + 3% YLL
#   3 x 2 x 2 x 2 = 24 scenarios. Wage and mortality-valuation variants are
#   run on the primary scenario only (separate tables later).
#===============================================================================

#--- Step 1: discount rates (defined in part 4) --------------------------------
stopifnot(base_dr == 0, sens_dr == 0.03)

#--- Step 2: values held fixed across the main grid ----------------------------
scen_defaults <- tibble(
  wage_day      = wage_usd,          # USD per work day (LFS 2024)
  annual_earn   = annual_earn_usd,   # USD per year, for deaths
  mort_val      = "lfs",             # deaths valued at LFS earnings
  lab_share     = NA_real_,          # used only in the GDP-based variant
  variant       = "main",
  ind_amb_ratio = ind_amb_ratio,     # outpatient indirect cost ratio
  oop_amb_ratio = oop_amb_ratio,     # outpatient OOP cost ratio
  provider_type = "base")

#--- Step 3: build the grid ----------------------------------------------------
# Costs and YLL are always discounted at the same rate
discount_scenarios <- tibble(r_cost = c(base_dr, sens_dr),
                             r_yll  = c(base_dr, sens_dr))

scen <- expand_grid(work_hz = c("stop60", "profile", "stop65"),
                    p_hosp  = c(p_hosp_primary, p_hosp_sa),
                    af_type = c("signed", "trunc")) %>%
  cross_join(discount_scenarios) %>%
  cross_join(scen_defaults) %>%
  mutate(scen_id = row_number(),
         primary = work_hz == "stop60" & r_cost == base_dr & r_yll == base_dr &
           p_hosp == p_hosp_primary & af_type == "signed")

# Exactly one primary scenario
stopifnot(nrow(scen) == 24, sum(scen$primary) == 1)
scen_primary <- scen %>% filter(primary)

#--- Step 4: per-death values for each year x horizon x discount rate ----------
# Only 3 x 2 = 6 combinations per year matter here (p_hosp and af_type don't
# affect per-death values), so 12 rows in total.
const_tab <- expand_grid(Year = out_years,
                         distinct(scen, work_hz, r_cost, r_yll)) %>%
  mutate(v = pmap(list(Year, r_yll, r_cost, work_hz), function(y, ry, rc, wh)
    age_consts(mortality_profiles_by_year[[as.character(y)]], ry, rc, wh)),
    yll_per_death = map_dbl(v, "yll_per_death"),
    pv_work       = map_dbl(v, "pv_work"),
    pv_life       = map_dbl(v, "pv_life")) %>%
  select(-v)
print(const_tab, n = Inf)

# Checks: employment-weighted years < all years; average YLL per death is
# below life expectancy at birth (deaths are not all in infancy)
stopifnot(all(const_tab$pv_work < const_tab$pv_life),
          all(const_tab$yll_per_death < le_at(0.5) - 1))


#===============================================================================
# B12. Burden engine: one year x one scenario, vectors of length n_sim
#   Inputs
#     AN          attributable cases per draw (signed or trunc, chosen beforehand)
#     cases_tot   total reported cases that year (national)
#     s           one row of scen
#     deaths_obs  observed dengue deaths that year (national)
#     kc          one row of const_tab (yll_per_death, pv_work, pv_life)
#     prog_usd_m  national dengue programme cost that year, USD million
#     gdp_pc      GDP per capita, USD (used only for the GDP-based valuations)
#     omega       multiplier on AN (1 = reported cases; >1 = under-reporting)
#   Output: one row per draw. Costs are in USD million.
#===============================================================================

# Fixed inputs: age groups and severity shares each add up to 1
stopifnot(abs(prop_paed + prop_adult - 1) < 1e-9,
          abs(sum(severity_weights(NULL, sev_paed_dhf,  sev_paed_dss))  - 1) < 1e-9,
          abs(sum(severity_weights(NULL, sev_adult_dhf, sev_adult_dss)) - 1) < 1e-9)

burden_engine <- function(AN, cases_tot, s, deaths_obs, kc, prog_usd_m, gdp_pc, omega = 1) {
  stopifnot(length(AN) == n_sim, nrow(kc) == 1, nrow(s) == 1, cases_tot > 0)
  if (s$mort_val == "gdp_ls") stopifnot(!is.na(s$lab_share))
  
  #--- Step 1: split attributable cases by setting and age ---------------------
  AN   <- AN * omega
  af   <- AN / cases_tot              # attributable fraction of reported cases
  hosp <- AN * s$p_hosp               # hospitalised
  amb  <- AN - hosp                   # outpatient (ambulatory)
  pa   <- hosp * prop_paed            # hospitalised children
  ad   <- hosp * prop_adult           # hospitalised adults
  
  #--- Step 2: severity mix (DF / DHF / DSS) -----------------------------------
  # sev_mix_mult exists only in the severity-mix scenarios; 1 otherwise
  mix_mult <- if ("sev_mix_mult" %in% names(s)) s$sev_mix_mult else 1
  sev_p <- severity_weights(sev_paed_df,  sev_paed_dhf,  sev_paed_dss,  mix_mult)
  sev_a <- severity_weights(sev_adult_df, sev_adult_dhf, sev_adult_dss, mix_mult)
  
  #--- Step 3: provider (health system) cost -----------------------------------
  # Outpatient cost: drawn value, or the fixed Shepard value in that scenario
  c_amb_draw <- if (s$provider_type == "shepard") rep(c_amb_sa_shepard, n_sim) else P$c_amb
  provider <- amb * c_amb_draw +
    pa * (sev_p[1] * P$c_p_df + sev_p[2] * P$c_p_dhf + sev_p[3] * P$c_p_dss) +
    ad * (sev_a[1] * P$c_a_df + sev_a[2] * P$c_a_dhf + sev_a[3] * P$c_a_dss)
  
  #--- Step 4: out-of-pocket cost ----------------------------------------------
  # Outpatient OOP = inpatient OOP x ratio
  oop <- hosp * P$oop + amb * P$oop * s$oop_amb_ratio
  
  #--- Step 5: indirect cost of illness (human capital) ------------------------
  # Adults: patient days lost (x employment rate) + caregiver days, x daily wage.
  # Wage = scenario level x shared draw multiplier (mean 1).
  ind_a     <- (P$days_pt * p_emp + P$days_cg) * s$wage_day * P$wage_mult
  # Children: caregiver cost by severity (fixed values, no uncertainty)
  ind_p_sev <- sev_p[1] * P$ind_p_df + sev_p[2] * P$ind_p_dhf + sev_p[3] * P$ind_p_dss
  # Outpatients: age-weighted inpatient value x outpatient ratio
  ind_amb   <- s$ind_amb_ratio * (prop_adult * ind_a + prop_paed * P$ind_p)
  indirect  <- ad * ind_a + pa * ind_p_sev + amb * ind_amb
  
  #--- Step 6: health losses ---------------------------------------------------
  # Attributable deaths: observed deaths x attributable fraction (assumes the
  # same case fatality in attributable and other cases)
  deaths <- af * deaths_obs
  yld    <- (amb * P$dw_amb + hosp * P$dw_hosp) * dengue_duration_years
  yll    <- deaths * kc$yll_per_death
  
  #--- Step 7: productivity lost through deaths --------------------------------
  per_death <- switch(s$mort_val,
                      lfs    = s$annual_earn * P$wage_mult * kc$pv_work,              # main: LFS earnings x employment by age
                      gdp    = gdp_pc * P$gdp_mult * kc$pv_life,                      # Shepard: GDP per capita, every year
                      gdp_ls = gdp_pc * s$lab_share * P$gdp_mult * kc$pv_life,        # labour share of GDP per capita
                      stop("unknown mort_val: ", s$mort_val))
  mort_cost <- deaths * per_death
  
  #--- Step 8: output, one row per draw (costs in USD million) -----------------
  tibble(sim = seq_len(n_sim), AN, Cases = cases_tot, Deaths = deaths,
         YLD = yld, YLL = yll, DALY = yld + yll,
         Provider  = provider  / 1e6,
         OOP       = oop       / 1e6,
         Indirect  = indirect  / 1e6,
         Mortality = mort_cost / 1e6,
         Societal  = (provider + oop + indirect + mort_cost) / 1e6,  # programme cost excluded
         Programme_alloc = prog_usd_m * af,                          # reported separately
         YLL_per_death = kc$yll_per_death, PV_work = kc$pv_work)
}


#===============================================================================
# B13. Runner: year-specific results + combined total (Year = 0)
#   For each scenario x year:
#     1. add up attributable cases across districts within each draw (national)
#     2. pick signed or truncated AN, depending on the scenario
#     3. look up that year's deaths, programme cost, GDP per capita and
#        per-death values
#     4. run burden_engine()
#   Then add the years together within each draw to get Year = 0.
#
#   Uncertainty is carried through by keeping draws paired: draw k of 2024 is
#   added to draw k of 2025, and summaries (mean, 95% interval) are taken only
#   at the end. This is the standard Monte Carlo approach for probabilistic
#   analysis (Briggs, Claxton & Sculpher 2006, ch. 4). Attributable numbers
#   from the DLNM are additive across districts and years, so summing them
#   per draw is valid (Gasparrini & Leone 2014).
#
#   References
#     Briggs AH, Claxton K, Sculpher MJ. Decision Modelling for Health Economic
#       Evaluation. Oxford University Press; 2006.
#     Gasparrini A, Leone M. Attributable risk from distributed lag models.
#       BMC Med Res Methodol. 2014;14:55.
#===============================================================================

#--- Step 1: helpers and column sets -------------------------------------------
# One value for one year; stops if the year is missing, repeated or NA
lookup1 <- function(tab, col, yr) {
  v <- tab[[col]][tab$Year == yr]
  if (length(v) != 1 || is.na(v))
    stop("lookup1: need exactly one non-missing '", col, "' for Year ", yr)
  v
}

# Counts and costs that are added across years for the combined total
sum_cols <- c("AN", "Cases", "Deaths", "YLD", "YLL", "DALY", "Provider", "OOP",
              "Indirect", "Mortality", "Societal", "Programme_alloc")

#--- Step 2: runner -------------------------------------------------------------
run_burden <- function(draws_df, omega = 1, label = "National_waste", scen_df = scen) {
  stopifnot(all(c("Year", "sim", "AN", "AN_trunc", "Cases") %in% names(draws_df)),
            !anyDuplicated(scen_df$scen_id))
  
  # 2a. National totals per year and draw (sum over districts)
  agg <- draws_df %>%
    group_by(Year, sim) %>%
    summarise(AN = sum(AN), AN_trunc = sum(AN_trunc), Cases = sum(Cases),
              .groups = "drop")
  
  # 2b. Per-year inputs, prepared once (not inside the scenario loop)
  yr_data <- map(set_names(out_years), function(yr) {
    d <- agg %>% filter(Year == yr) %>% arrange(sim)
    # one row per draw, and reported cases are the same in every draw
    stopifnot(nrow(d) == n_sim, identical(d$sim, seq_len(n_sim)),
              n_distinct(d$Cases) == 1)
    list(d          = d,
         deaths_obs = lookup1(observed_deaths_lookup, "Observed_Deaths", yr),
         prog_usd_m = lookup1(programme_lookup, "Total_Vector_Control_USD_M", yr),
         gdp_pc     = lookup1(gdp_lookup, "GDP_Capita_USD", yr))
  })
  
  # 2c. Engine for every scenario x year
  by_year <- map_dfr(seq_len(nrow(scen_df)), function(k) {
    s <- scen_df[k, ]
    map_dfr(out_years, function(yr) {
      yd <- yr_data[[as.character(yr)]]
      kc <- const_tab %>% filter(Year == yr, work_hz == s$work_hz,
                                 r_cost == s$r_cost, r_yll == s$r_yll)
      if (nrow(kc) != 1)
        stop("No unique const_tab row for Year ", yr, ", scenario ", s$scen_id)
      
      burden_engine(
        AN         = if (s$af_type == "signed") yd$d$AN else yd$d$AN_trunc,
        cases_tot  = yd$d$Cases[1],
        s          = s,
        deaths_obs = yd$deaths_obs,
        kc         = kc,
        prog_usd_m = yd$prog_usd_m,
        gdp_pc     = yd$gdp_pc,
        omega      = omega) %>%
        mutate(Year = yr, scen_id = s$scen_id)
    })
  })
  
  # 2d. Combined total (Year = 0): add the years within each draw.
  #     YLL_per_death is the death-weighted average across years;
  #     PV_work has no single combined value, so it is left NA.
  combined <- by_year %>%
    group_by(scen_id, sim) %>%
    summarise(across(all_of(sum_cols), sum), .groups = "drop") %>%
    mutate(Year = 0, YLL_per_death = YLL / Deaths, PV_work = NA_real_)
  
  # 2e. Attach scenario settings and the analysis label
  out <- bind_rows(by_year, combined) %>%
    left_join(scen_df, by = "scen_id") %>%
    mutate(Analysis = label)
  
  #--- Step 3: checks ------------------------------------------------------------
  # Each scenario has n_sim rows for each year and for the combined total
  stopifnot(nrow(out) == nrow(scen_df) * n_sim * (length(out_years) + 1))
  
  # Combined = sum of years, for the main outputs
  chk <- out %>%
    group_by(scen_id, sim) %>%
    summarise(across(c(Deaths, DALY, Societal),
                     ~ sum(.x[Year != 0]) - sum(.x[Year == 0])),
              .groups = "drop")
  stopifnot(max(abs(as.matrix(chk[, c("Deaths", "DALY", "Societal")]))) < 1e-8)
  
  out
}

#===============================================================================
# B14. Main run + summaries
#   Each summary is a mean with a 95% percentile interval across the n_sim
#   draws (2.5th and 97.5th percentiles), the standard way of reporting
#   uncertainty from probabilistic analysis (Briggs, Claxton & Sculpher 2006,
#   ch. 4). Summaries are taken only at the end, after years and districts
#   have been added within each draw.
#
#   Year = 0 is the combined 2024 + 2025 total; it is labelled "Combined"
#   in the printed tables.
#
#   References
#     Briggs AH, Claxton K, Sculpher MJ. Decision Modelling for Health Economic
#       Evaluation. Oxford University Press; 2006.
#===============================================================================

#--- Step 1: summary helpers -----------------------------------------------------
qtl <- function(x, p) quantile(x, p, names = FALSE)
mlu <- list(M = mean, L = ~ qtl(.x, .025), U = ~ qtl(.x, .975))   # mean, lower, upper

year_label <- function(y) if_else(y == 0, "Combined", as.character(y))

#--- Step 2: main run (24 scenarios) ---------------------------------------------
sim_nat <- run_burden(an_draws)
stopifnot(nrow(sim_nat) == nrow(scen) * (length(out_years) + 1) * n_sim)

# Share of draws with negative AN (signed scenario only; see part 6)
neg_share <- sim_nat %>%
  filter(primary) %>%
  group_by(Year) %>%
  summarise(neg_AN_pct = 100 * mean(AN < 0), .groups = "drop") %>%
  mutate(Year = year_label(Year))
print(neg_share)

#--- Step 3: health-loss breakdown, primary scenario -------------------------------
#   AF_M        attributable fraction of reported cases (%)
#   YLL_per_death  ratio of means, so it is stable when signed deaths are near 0
#   YLD_pct     share of DALYs from illness rather than death
decomp <- sim_nat %>%
  filter(primary) %>%
  group_by(Year) %>%
  summarise(AN_M     = mean(AN),
            AF_M     = 100 * mean(AN / Cases),
            Deaths_M = mean(Deaths),
            YLD_M    = mean(YLD),
            YLL_M    = mean(YLL),
            DALY_M   = mean(DALY),
            DALY_L   = qtl(DALY, .025),
            DALY_U   = qtl(DALY, .975),
            YLL_per_death = mean(YLL) / mean(Deaths),
            YLD_pct  = 100 * mean(YLD) / mean(DALY),
            .groups = "drop") %>%
  mutate(Year = year_label(Year))
print(decomp, width = Inf)

#--- Step 4: cost components, primary scenario (USD million) -------------------------
#   Societal = Provider + OOP + Indirect + Mortality. The programme cost is
#   shown alongside but not included (see part 6).
comp <- sim_nat %>%
  filter(primary) %>%
  group_by(Year) %>%
  summarise(across(c(Provider, OOP, Indirect, Mortality, Societal, Programme_alloc), mlu),
            .groups = "drop") %>%
  mutate(Year = year_label(Year))
print(comp, width = Inf)

# Means add up exactly (percentile bounds do not, so only means are checked)
stopifnot(all(abs(comp$Provider_M + comp$OOP_M + comp$Indirect_M + comp$Mortality_M -
                    comp$Societal_M) < 1e-8))

#--- Step 5: all scenarios, by year and combined --------------------------------------
scen_summary <- sim_nat %>%
  group_by(Year, scen_id, work_hz, r_cost, r_yll, p_hosp, af_type, primary) %>%
  summarise(across(c(DALY, Indirect, Mortality, Societal), mlu),
            Societal_Md = median(Societal),
            .groups = "drop") %>%
  arrange(Year == 0, Year, scen_id) %>%          # 2024, 2025, then Combined
  mutate(Year = year_label(Year))
print(scen_summary, n = Inf, width = Inf)

#===============================================================================
# B15. Wage and mortality-valuation sensitivities (primary scenario only)
#   One input is changed at a time; everything else stays at the primary
#   values. All variants use the same AN draws and the same parameter draws
#   (P), so differences between variants come only from the changed input,
#   not from simulation noise ("common random numbers"; Stout & Goldie 2008).
#
#   Valuation of deaths:
#     lfs     annual LFS earnings x employment-weighted working years (main)
#     gdp     GDP per capita x all remaining life years, no employment
#             adjustment, as in Shepard et al. 2013
#     gdp_ls  GDP per capita x labour share of income (low / high bound)
#   Valuing lost output by earnings is the human-capital approach
#   (Drummond et al. 2015, ch. 4).
#
#   References
#     Stout NK, Goldie SJ. Keeping the noise down: common random numbers for
#       disease simulation modeling. Health Care Manag Sci. 2008;11:399-406.
#     Shepard DS, Undurraga EA, Halasa YA. Economic and disease burden of dengue
#       in Southeast Asia. PLoS Negl Trop Dis. 2013;7(2):e2055.
#     Drummond MF, Sculpher MJ, Claxton K, Stoddart GL, Torrance GW. Methods for
#       the Economic Evaluation of Health Care Programmes. 4th ed. Oxford
#       University Press; 2015.
#===============================================================================

#--- Step 1: build the variants -----------------------------------------------------
# Start from the primary scenario, reset the severity multiplier, then change
# the named inputs.
mk_var <- function(lab, ...) {
  scen_primary %>% mutate(variant = lab, sev_mix_mult = 1) %>% mutate(...)
}
main_label <- sprintf("Main: LFS %.2f/day, E/P profile", wage_usd)
usd <- function(x) sprintf("%.2f", x)

scen_val <- bind_rows(
  mk_var(main_label),
  # Hospital severity mix: DHF + DSS share scaled, DHF:DSS ratio kept (part 3)
  mk_var("Hospital severity mix: severe share -20% relative", sev_mix_mult = 0.8),
  mk_var("Hospital severity mix: severe share +20% relative", sev_mix_mult = 1.2),
  # Outpatient costs relative to inpatient
  mk_var("Indirect ambulatory ratio: 0.3", ind_amb_ratio = 0.3),
  mk_var("Indirect ambulatory ratio: 0.8", ind_amb_ratio = 0.8),
  mk_var("OOP ambulatory ratio: 0",        oop_amb_ratio = 0),
  mk_var("OOP ambulatory ratio: 1",        oop_amb_ratio = 1),
  mk_var("Ambulatory provider: Shepard 2013 unit cost", provider_type = "shepard"),
  # Daily wage: changes both illness days and deaths
  mk_var(paste("Wage: Weerasinghe x IPSWRI", usd(wage_usd_ipswri)),
         wage_day = wage_usd_ipswri, annual_earn = wage_usd_ipswri * work_days_yr),
  mk_var(paste("Wage: LFS daily-paid", usd(wage_usd_lfsdly)),
         wage_day = wage_usd_lfsdly, annual_earn = wage_usd_lfsdly * work_days_yr),
  mk_var(paste("Wage: legal minimum", usd(wage_usd_min)),
         wage_day = wage_usd_min,    annual_earn = wage_usd_min * work_days_yr),
  mk_var(paste("Wage: CPI-inflated", usd(wage_usd_cpi), "(old)"),
         wage_day = wage_usd_cpi,    annual_earn = wage_usd_cpi * work_days_yr),
  # Same annual earnings spread over 312 days: lower daily rate, so only the
  # illness-day cost changes; the cost of deaths stays the same
  mk_var("Working days 312", wage_day = annual_earn_usd / 312),
  # Valuation of deaths
  mk_var("Deaths: GDP pc, no E/P (Shepard)", mort_val = "gdp"),
  mk_var(sprintf("Deaths: GDP pc x %.1f%%", 100 * lab_share_lo),
         mort_val = "gdp_ls", lab_share = lab_share_lo),
  mk_var(sprintf("Deaths: GDP pc x %.1f%%", 100 * lab_share_hi),
         mort_val = "gdp_ls", lab_share = lab_share_hi)) %>%
  mutate(scen_id = 1000 + row_number())   # 1000+ so IDs never clash with scen

stopifnot(nrow(scen_val) == 16, !anyDuplicated(scen_val$variant))

#--- Step 2: run and summarise, by year and combined ------------------------------
sim_val <- run_burden(an_draws, label = "Valuation_SA", scen_df = scen_val)

val_sens <- sim_val %>%
  group_by(Year, scen_id, variant) %>%
  summarise(across(c(Indirect, Mortality, Societal), mlu), .groups = "drop") %>%
  arrange(Year == 0, Year, scen_id) %>%
  mutate(Year = year_label(Year))
print(val_sens, n = Inf, width = Inf)

#--- Step 3: check -----------------------------------------------------------------
# The main variant must reproduce the primary scenario exactly (same draws),
# for each year and combined.
main_chk <- val_sens %>%
  filter(variant == main_label) %>%
  select(Year, Societal_M) %>%
  inner_join(comp %>% select(Year, Societal_M_prim = Societal_M), by = "Year")
stopifnot(nrow(main_chk) == length(out_years) + 1,
          max(abs(main_chk$Societal_M - main_chk$Societal_M_prim)) < 1e-8)


#===============================================================================
# B16. National summary table (all scenarios) + Table 3 (primary scenario)
#   One row per scenario x year, with mean and 95% percentile interval for
#   each output (Briggs, Claxton & Sculpher 2006, ch. 4).
#
#   AF (attributable fraction, %) = 100 x AN / reported cases, calculated
#   within each draw and then summarised. Reported cases are fixed across
#   draws within a year, so this equals 100 x mean(AN) / Cases.
#
#   Year is kept as a number in tab_nat (0 = combined 2024 + 2025) so it can
#   be filtered directly; Table 3 shows it as "Combined".
#===============================================================================

#--- Step 1: summary of every scenario ------------------------------------------
scen_cols <- setdiff(names(scen), "scen_id")     # scenario settings (not outputs)
vars <- c("AF", setdiff(sum_cols, "Cases"))      # Cases is fixed; no interval needed

tab_nat <- sim_nat %>%
  mutate(AF = 100 * AN / Cases) %>%
  group_by(Analysis, scen_id, Year) %>%
  summarise(across(all_of(vars), mlu),
            Societal_Md = median(Societal),
            .groups = "drop") %>%
  left_join(scen, by = "scen_id")

# One row per scenario x (each year + combined)
stopifnot(nrow(tab_nat) == nrow(scen) * (length(out_years) + 1))

# Primary DALYs must match the decomposition from part 8, for every year
dal_chk <- tab_nat %>%
  filter(primary) %>%
  transmute(Year = year_label(Year), DALY_M) %>%
  inner_join(decomp %>% select(Year, DALY_M_dec = DALY_M), by = "Year")
stopifnot(nrow(dal_chk) == length(out_years) + 1,
          max(abs(dal_chk$DALY_M - dal_chk$DALY_M_dec)) < 1e-8)

#--- Step 2: Table 3 (primary scenario) -------------------------------------------
#   Columns: AN, AF, deaths, YLD, YLL, DALY and societal cost (USD million),
#   each as mean (_M), lower (_L) and upper (_U) 95% bound, plus the median
#   societal cost (Societal_Md).
tab3_vars <- c("AN", "AF", "Deaths", "YLD", "YLL", "DALY", "Societal")
tab3_prefixes <- paste0(tab3_vars, "_")

tab3 <- tab_nat %>%
  filter(primary) %>%
  select(Year, starts_with(tab3_prefixes, ignore.case = FALSE)) %>%
  arrange(Year == 0, Year) %>%                   # 2024, 2025, then Combined
  mutate(Year = year_label(Year))

stopifnot(
  nrow(tab3) == length(out_years) + 1,
  !any(scen_cols %in% names(tab3)),              # no scenario settings leaked in
  all(c(outer(tab3_vars, c("_M", "_L", "_U"), paste0)) %in% names(tab3)),
  all(tab3$DALY_L <= tab3$DALY_M & tab3$DALY_M <= tab3$DALY_U)
)
print(tab3, width = Inf)

#===============================================================================
# B17. Reference-period sensitivity (P0 / P10 / P25 / P50), primary scenario
#   Attributable cases depend on the reference (expected) level of cases that
#   observed cases are compared against. The main analysis uses P<cen_main>;
#   here the same analysis is repeated with the other reference levels.
#   A higher reference level means more cases are counted as expected, so
#   fewer are attributable: AF should fall from P0 to P50.
#
#   Only the reference level changes. Reported cases, deaths, costs and
#   parameter draws are the same in every run, so differences come only from
#   AN (one-way sensitivity analysis; Briggs, Claxton & Sculpher 2006, ch. 4).
#===============================================================================

#--- Step 1: AN draws for each reference level ----------------------------------
ref_levels <- c(0, 10, 25, 50)
ref_names  <- paste0("P", ref_levels)
stopifnot(cen_main %in% ref_levels)

# Reuse the main draws for the main level, so that run reproduces sim_nat
ref_draws <- set_names(
  map(ref_levels, ~ if (.x == cen_main) an_draws else build_an_draws(.x)),
  ref_names)

# Reported cases must be identical across reference levels (only AN changes)
cases_by_ref <- map(ref_draws, ~ .x %>%
                      filter(sim == 1, Year %in% out_years) %>%
                      arrange(District, Year) %>%
                      pull(Cases))
stopifnot(all(map_lgl(cases_by_ref, ~ identical(.x, cases_by_ref[[1]]))))

#--- Step 2: run and summarise, by year and combined ---------------------------
sim_ref <- imap_dfr(ref_draws, ~ run_burden(.x, label = .y, scen_df = scen_primary))

ref_sens_econ <- sim_ref %>%
  mutate(AF = 100 * AN / Cases) %>%
  group_by(Analysis, Year) %>%
  summarise(across(c(AF, DALY, Societal), mlu),
            Societal_Md = median(Societal),
            Programme_alloc_M = mean(Programme_alloc),
            .groups = "drop") %>%
  arrange(Year == 0, Year, factor(Analysis, levels = ref_names))

print(ref_sens_econ %>% mutate(Year = year_label(Year)), n = Inf, width = Inf)

#--- Step 3: checks -----------------------------------------------------------------
# The main reference level must reproduce the primary scenario, every year
ref_chk <- ref_sens_econ %>%
  filter(Analysis == paste0("P", cen_main)) %>%
  select(Year, DALY_M, Societal_M) %>%
  inner_join(tab_nat %>% filter(primary) %>%
               select(Year, DALY_M_nat = DALY_M, Societal_M_nat = Societal_M),
             by = "Year")
stopifnot(nrow(ref_chk) == length(out_years) + 1,
          max(abs(ref_chk$DALY_M - ref_chk$DALY_M_nat)) < 1e-8,
          max(abs(ref_chk$Societal_M - ref_chk$Societal_M_nat)) < 1e-8)

# AF should fall as the reference level rises (P0 -> P50), in every year
af_order <- ref_sens_econ %>%
  group_by(Year) %>%
  arrange(factor(Analysis, levels = ref_names), .by_group = TRUE) %>%
  summarise(monotone = !is.unsorted(rev(AF_M)), .groups = "drop")
if (!all(af_order$monotone))
  warning("AF does not fall from P0 to P50 in: ",
          paste(year_label(af_order$Year[!af_order$monotone]), collapse = ", "),
          ". Inspect ref_sens_econ.")



#===============================================================================
# B18. Tornado analysis (pointer)
#   The one-way (tornado) sensitivity analysis is kept in a separate script:
#   tornado_oneway_sensitivity.R. Run it after this script has built `comp`
#   (part 8). It uses Year = "Combined" for the 2024 + 2025 total.
#===============================================================================

#===============================================================================
# B19. Expected-range check on the primary results (soft check)
#   Purpose: catch wiring errors after a revision (wrong argument, unit error,
#   dropped component). Gives a warning, not an error, because the ranges are
#   judgement-based, not identities. The identities (components sum to the
#   societal total; 2024 + 2025 = combined) are hard checks in parts 8-9.
#
#   Ranges are for the primary scenario, combined 2024 + 2025, million USD.
#   Revision note: the indirect upper bound was widened from 3.2 to 3.5 after
#   p_emp changed from 524/1064 to 524/994 (adult denominator). Patient work
#   loss rises by about 7% (1064/994), so indirect is now ~3.22 (was ~3.0).
#   Update these ranges whenever an input changes on purpose.
#===============================================================================

#--- Step 1: combined results vs expected ranges --------------------------------
p0 <- comp %>% filter(Year == "Combined")
stopifnot(nrow(p0) == 1)

chk_tab <- tibble(
  Component = c("Provider+OOP", "Indirect", "Mortality", "Societal"),
  Actual    = c(p0$Provider_M + p0$OOP_M, p0$Indirect_M,
                p0$Mortality_M, p0$Societal_M),
  Exp_lo    = c(11.5, 2.8, 0.0, 14.5),
  Exp_hi    = c(12.5, 3.5, 0.6, 16.0)) %>%
  mutate(OK = Actual >= Exp_lo & Actual <= Exp_hi)
print(chk_tab)

if (!all(chk_tab$OK))
  warning("Primary results (combined 2024+2025) outside expected ranges: ",
          paste(chk_tab$Component[!chk_tab$OK], collapse = ", "),
          ". Check wiring before reporting.")

#--- Step 2: each year's share of the combined total -----------------------------
# The two years have similar case numbers, so each should be roughly half.
# A share near 0% or 100% means a year was dropped or counted twice.
yr_shares <- comp %>%
  filter(Year != "Combined") %>%
  select(Year, Provider_M, OOP_M, Indirect_M, Mortality_M, Societal_M) %>%
  pivot_longer(-Year, names_to = "Component", values_to = "Value") %>%
  group_by(Component) %>%
  mutate(Share_pct = 100 * Value / sum(Value)) %>%
  ungroup()

chk_year <- yr_shares %>%
  mutate(OK = Share_pct > 30 & Share_pct < 70)
print(chk_year, n = Inf)

if (!all(chk_year$OK))
  warning("Unusual year split (outside 30-70% of combined) for: ",
          paste(unique(paste(chk_year$Component[!chk_year$OK],
                             chk_year$Year[!chk_year$OK])), collapse = ", "),
          ". Check the year inputs (cases, GDP, wages) before reporting.")


#===============================================================================
# PABI eight-sentinel-district pilot extension (append-only; 5 Nov 2026)
#   Run AFTER Section 5 attribution and the final Version 2 national block
#   (B0-B4) in the locked PABI_econmic_Nov-5.R.
#   Pilot blocks parallel national BRIDGES B0-B4. Locked national economics,
#   draws, functions, scenario grids and table conventions are REUSED, not
#   redefined. No model is refitted and no draws are resampled.
#===============================================================================

#-------------------------------------------------------------------------------
# PL1. Preconditions: every upstream object must exist and match the lock
#-------------------------------------------------------------------------------
pilot_required <- c(
  # Section 5 attribution
  "store", "draws", "beta_pool", "make_ref_row", "cen_main", "ref_levels",
  "res_main", "df",
  # National B0 economic and demographic inputs
  "programme_lookup", "gdp_lookup", "observed_deaths_lookup",
  "mortality_profiles_by_year", "age_groups", "age_mids", "age_consts",
  "const_tab", "fx_2024", "fx_2025", "lkr25_to_usd24",
  # National B1-B4 engine, scenarios and summary helpers
  "P", "n_sim", "out_years", "scen", "scen_primary", "scen_val",
  "burden_engine", "lookup1", "sum_cols", "mlu"
)
pilot_missing <- pilot_required[!vapply(pilot_required, exists, logical(1),
                                        envir = .GlobalEnv, inherits = TRUE)]
if (length(pilot_missing))
  stop("Run attribution and the locked national pipeline first; missing: ",
       paste(pilot_missing, collapse = ", "))

if (!identical(as.numeric(cen_main), 10))
  stop("Required main reference is cen_main = 10.")
if (!identical(as.numeric(ref_levels), c(0, 10, 25, 50)))
  stop("Locked national `ref_levels` must be c(0, 10, 25, 50).")
if (!is.matrix(draws) || nrow(draws) != n_sim || nrow(draws) != 5000L ||
    ncol(draws) != length(beta_pool))
  stop("Shared draws must be 5000 x length(beta_pool); do not resample.")
if (!is.list(res_main) || is.null(res_main$district_res))
  stop("`res_main` must contain named district_res from Section 5.")
if (!is.list(store) || is.null(names(store)))
  stop("`store` must be a named district list.")

message("PL1 preconditions passed: ", length(pilot_required), " objects found; ",
        length(store), " districts in store; cen_main = P", cen_main, "; ",
        n_sim, " shared draws.")

#===============================================================================
# BRIDGE B0 (pilot): population and allocation weight only
#   National economics (unit costs, FX/CPI, GDP, DALY constants) are reused.
#===============================================================================
# Pilot district populations (fixed inputs for the eight sentinel districts)
pilot_population <- c(
  Jaffna = 594751, Monaragala = 527585, Kandy = 1461895, Galle = 1097372,
  Ratnapura = 1145423, Gampaha = 2436142, Kalutara = 1305784,
  Kurunegala = 1768156
)
pilot_districts        <- names(pilot_population)
pilot_population_total <- sum(pilot_population)

stopifnot(pilot_population_total == 10337108L,
          all(pilot_districts %in% names(store)),
          all(pilot_districts %in% names(res_main$district_res)))

# National population from df: exactly one positive value per district
if (!all(c("District", "Population") %in% names(df)))
  stop("National `df` needs District and Population to allocate the operational budget.")

pilot_pop_rows <- unique(df[!is.na(df$District) & !is.na(df$Population) &
                              df$Population > 0, c("District", "Population")])
pilot_pop_n <- dplyr::count(pilot_pop_rows, District, name = "n_pop")
if (any(pilot_pop_n$n_pop != 1L) || nrow(pilot_pop_rows) != nrow(pilot_pop_n))
  stop("`df` must provide exactly one unique positive Population per national District.")

pilot_national_population_total <- sum(pilot_pop_rows$Population)
if (!is.finite(pilot_national_population_total) ||
    pilot_national_population_total < pilot_population_total)
  stop("Invalid national population total in df$Population.")

# Cross-check: hard-coded pilot populations vs df (warning only, not a stop)
pilot_pop_check <- tibble::tibble(District = pilot_districts,
                                  Pop_pilot = unname(pilot_population)) %>%
  dplyr::left_join(pilot_pop_rows, by = "District") %>%
  dplyr::mutate(Diff_pct = round(100 * (Pop_pilot - Population) / Population, 2))
if (any(is.na(pilot_pop_check$Population)) || any(abs(pilot_pop_check$Diff_pct) > 0.5, na.rm = TRUE)) {
  warning("Pilot populations differ from df$Population by more than 0.5% (or are missing):")
  print(pilot_pop_check)
}

# Share of the national population covered by the pilot (used to allocate budget)
pilot_population_weight <- pilot_population_total / pilot_national_population_total

message("Pilot B0: ", length(pilot_districts), " districts; population ",
        format(pilot_population_total, big.mark = ","), " (",
        round(100 * pilot_population_weight, 1), "% of national ",
        format(pilot_national_population_total, big.mark = ","), ").")

#===============================================================================
# BRIDGE B1-B2 (pilot): backward AN draws by year (shared draws; no refit)
#   Section 5 `store` fields cb_full / cases_full / year_full feed a year split
#   parallel to locked attr_one_district_yr / build_an_draws, with the SAME
#   `draws`. At cen_main, 2024 + 2025 must reconcile with res_main$district_res.
#
#   eta_{d,t,k} = (w_{d,t} - w_0)' beta_k          for draw k
#   AF_{d,t,k}  = 1 - exp(-eta_{d,t,k})
#   AN_{d,y,k}  = sum over weeks t in year y of n_{d,t} * AF_{d,t,k}
#===============================================================================

#--- B1. Year-split attribution for one district --------------------------------
pilot_attr_by_year <- function(s, ref_row) {
  stopifnot(all(c("cb_full", "cases_full", "year_full") %in% names(s)))
  
  cbc <- sweep(unclass(s$cb_full), 2, ref_row)        # w_{d,t} - w_0
  ok  <- !is.na(s$cases_full) & complete.cases(cbc)   # same row filter as Section 5
  X   <- cbc[ok, , drop = FALSE]
  n   <- s$cases_full[ok]
  yr  <- s$year_full[ok]
  
  if (!length(n) || !all(c(2024, 2025) %in% yr))
    stop("Pilot district lacks valid 2024/2025 rows.")
  
  af_sim <- 1 - exp(-(X %*% t(draws)))                # n_weeks x n_sim
  
  list(cases        = rowsum(n, yr)[, 1],                    # cases by year
       AN_sim       = rowsum(n * af_sim, yr),                # years x n_sim (signed)
       AN_sim_trunc = rowsum(n * pmax(af_sim, 0), yr))       # protective weeks set to 0 (SA)
}

#--- Helper: pull the Section 5 per-draw AN vector for one district -------------
#    Accepts the same shapes the original block handled.
pilot_s5_target <- function(target) {
  if (is.list(target) && !is.data.frame(target)) {
    if (!is.null(target$AN_sim))      target <- target$AN_sim
    else if (!is.null(target$AN))     target <- target$AN
  }
  if (is.data.frame(target)) {
    if ("AN_sim" %in% names(target))                 target <- target$AN_sim
    else if (all(c("sim", "AN") %in% names(target))) target <- target$AN[order(target$sim)]
  }
  if (is.matrix(target) && nrow(target) != n_sim && ncol(target) == n_sim)
    target <- colSums(target)                         # years x draws -> draws
  as.numeric(target)
}

#--- B2. Long draw table (District x Year x sim) and reconciliation ------------
make_pilot_an_draws <- function(cen) {
  ref <- make_ref_row(cen)
  
  res <- purrr::set_names(pilot_districts) %>%
    purrr::map(~ pilot_attr_by_year(store[[.x]], ref))
  
  out <- purrr::imap_dfr(res, function(r, d) {
    purrr::map_dfr(rownames(r$AN_sim), function(y) tibble::tibble(
      District = d,
      Year     = as.integer(y),
      sim      = seq_len(n_sim),
      Cases    = unname(r$cases[y]),
      AN       = as.numeric(r$AN_sim[y, ]),
      AN_trunc = as.numeric(r$AN_sim_trunc[y, ])))
  }) %>%
    dplyr::mutate(cen = cen)
  
  # Structural checks
  stopifnot(!anyNA(out$AN), !anyNA(out$AN_trunc),
            all(c(2024, 2025) %in% out$Year),
            all(dplyr::count(out, District, Year)$n == n_sim),
            setequal(unique(out$District), pilot_districts))
  
  # Reconciliation with Section 5 (main reference only)
  if (isTRUE(all.equal(as.numeric(cen), as.numeric(cen_main)))) {
    for (d in pilot_districts) {
      target <- pilot_s5_target(res_main$district_res[[d]])
      got    <- colSums(res[[d]]$AN_sim)              # 2024 + 2025, per draw
      if (length(target) != n_sim)
        stop("Section 5 target for ", d, " is not a length-", n_sim, " draw vector.")
      tol <- 1e-8 * max(1, abs(target))
      if (max(abs(got - target)) > tol)
        stop("Pilot annual draws do not reconcile to res_main$district_res[[", d, "]].")
    }
    message("B1-B2: ", length(pilot_districts), " pilot districts reconcile with ",
            "Section 5 (", n_sim, " draws, P", cen, ").")
  }
  out
}

pilot_an_draws <- make_pilot_an_draws(cen_main)

#--- Quick summary: mean AN and 95% UI by district and year ---------------------
pilot_an_summary <- pilot_an_draws %>%
  dplyr::group_by(District, Year, Cases) %>%
  dplyr::summarise(AN_mean = mean(AN),
                   AN_low  = quantile(AN, 0.025),
                   AN_high = quantile(AN, 0.975),
                   .groups = "drop") %>%
  dplyr::mutate(AF_pct = round(100 * AN_mean / Cases, 1))
print(pilot_an_summary, n = Inf)


#-------------------------------------------------------------------------------
# B2b (pilot). Plastic material-share draws
#   For district d:  plastic_pos = sum(pl_A + pl_B + pl_AB)
#                    other_pos   = sum(A + B + AB)
#   Posterior share  s_d ~ Beta(plastic_pos + 1, other_pos + 1)   (uniform prior)
#   One independent share per simulation k = 1..n_sim, so it lines up row-wise
#   with the shared AN draws. Isolated seed; caller's RNG state is restored.
#-------------------------------------------------------------------------------

#--- Input path (default, or user override via `pilot_plastic_csv`) ------------
pilot_plastic_csv_default <- "~/Documents/Heidelberg/Plastics_Indicator/sri_lanka_indicator/data/pl_focussed_rev_4.csv"
pilot_plastic_csv_path <- path.expand(
  if (exists("pilot_plastic_csv", envir = .GlobalEnv, inherits = TRUE))
    pilot_plastic_csv else pilot_plastic_csv_default
)
if (!file.exists(pilot_plastic_csv_path))
  stop("Plastic CSV not found; set `pilot_plastic_csv`: ", pilot_plastic_csv_path)

#--- Read (blank headers become "...9" etc. via name_repair) -------------------
pilot_plastic_raw <- readr::read_csv(pilot_plastic_csv_path,
                                     show_col_types = FALSE,
                                     name_repair = "unique")

# Safety net: name any remaining empty or NA column names
pilot_col_names <- names(pilot_plastic_raw)
pilot_bad_names <- is.na(pilot_col_names) | pilot_col_names == ""
if (any(pilot_bad_names))
  names(pilot_plastic_raw)[pilot_bad_names] <- paste0("unnamed_", which(pilot_bad_names))

#--- Required columns ----------------------------------------------------------
pilot_plastic_outcomes <- c("pl_Wet", "pl_Dry", "pl_A", "pl_B", "pl_AB",
                            "Wet", "Dry", "A", "B", "AB")
pilot_missing_cols <- setdiff(c("District", pilot_plastic_outcomes), names(pilot_plastic_raw))
if (length(pilot_missing_cols))
  stop("Plastic input is missing: ", paste(pilot_missing_cols, collapse = ", "),
       ". It must contain District and all 10 outcomes.")

#--- Clean: trim names, drop all-empty rows, coerce outcomes to numeric --------
pilot_plastic_clean <- pilot_plastic_raw %>%
  dplyr::mutate(District = trimws(District)) %>%
  dplyr::filter(!dplyr::if_all(dplyr::all_of(pilot_plastic_outcomes), is.na))

# Warn if any non-blank entry is not numeric (it would otherwise become 0 silently)
pilot_nonnum <- vapply(pilot_plastic_outcomes, function(v) {
  x <- pilot_plastic_clean[[v]]
  sum(!is.na(x) & is.na(suppressWarnings(as.numeric(x))))
}, integer(1))
if (any(pilot_nonnum > 0))
  warning("Non-numeric plastic entries set to 0: ",
          paste(names(pilot_nonnum)[pilot_nonnum > 0], pilot_nonnum[pilot_nonnum > 0],
                sep = "=", collapse = ", "))

pilot_plastic_clean <- pilot_plastic_clean %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(pilot_plastic_outcomes),
                              ~ tidyr::replace_na(suppressWarnings(as.numeric(.x)), 0)))

# Warn about CSV districts that don't match a pilot name (likely spelling)
pilot_csv_unmatched <- setdiff(unique(pilot_plastic_clean$District),
                               c(pilot_districts, NA, ""))
if (length(pilot_csv_unmatched))
  message("Note: CSV districts not in the pilot set (ignored): ",
          paste(pilot_csv_unmatched, collapse = ", "),
          ". Check spelling if any should be a pilot district.")

#--- District counts ------------------------------------------------------------
pilot_plastic_counts <- pilot_plastic_clean %>%
  dplyr::filter(District %in% pilot_districts) %>%
  dplyr::group_by(District) %>%
  dplyr::summarise(plastic_pos = sum(pl_A + pl_B + pl_AB),
                   other_pos   = sum(A + B + AB),
                   .groups = "drop")

pilot_missing_districts <- setdiff(pilot_districts, pilot_plastic_counts$District)
if (length(pilot_missing_districts))
  warning("No plastic survey rows for: ", paste(pilot_missing_districts, collapse = ", "),
          " (they will use a flat Beta(1,1) share).")

# Complete to all 8 pilot districts, in pilot order
pilot_plastic_counts <- tibble::tibble(District = pilot_districts) %>%
  dplyr::left_join(pilot_plastic_counts, by = "District") %>%
  dplyr::mutate(dplyr::across(c(plastic_pos, other_pos), ~ tidyr::replace_na(.x, 0)))

if (any(pilot_plastic_counts$plastic_pos < 0 | pilot_plastic_counts$other_pos < 0))
  stop("Plastic counts must be non-negative.")

#--- Beta draws in an isolated RNG stream --------------------------------------
pilot_draw_plastic_share <- function(counts, n, seed = 20261105) {
  # Save the caller's RNG state and restore it on exit, even if an error occurs
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  
  set.seed(seed)
  out <- matrix(NA_real_, nrow = n, ncol = nrow(counts),
                dimnames = list(NULL, counts$District))
  for (j in seq_len(nrow(counts)))          # same order as before -> identical draws
    out[, j] <- stats::rbeta(n, counts$plastic_pos[j] + 1, counts$other_pos[j] + 1)
  out
}

pilot_rng_before <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
pilot_plastic_share <- pilot_draw_plastic_share(pilot_plastic_counts, n_sim)
stopifnot(identical(pilot_rng_before,
                    if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL),
          identical(dim(pilot_plastic_share), c(as.integer(n_sim), length(pilot_districts))),
          !anyNA(pilot_plastic_share),
          all(pilot_plastic_share >= 0 & pilot_plastic_share <= 1))

#--- Summary: posterior mean and 95% interval by district -----------------------
pilot_plastic_summary <- pilot_plastic_counts %>%
  dplyr::mutate(
    n_pos      = plastic_pos + other_pos,
    share_mean = (plastic_pos + 1) / (n_pos + 2),          # analytic Beta mean
    share_low  = apply(pilot_plastic_share, 2, quantile, 0.025)[District],
    share_high = apply(pilot_plastic_share, 2, quantile, 0.975)[District]
  )
print(pilot_plastic_summary, n = Inf)

message("B2b: plastic-share draws ready (", n_sim, " x ", length(pilot_districts),
        "); RNG state restored.")

#-------------------------------------------------------------------------------
# B2c (pilot). Plastic-attributable AN draws
#   AN_plastic[d, y, k]       = AN[d, y, k]       * s_d[k]
#   AN_plastic_trunc[d, y, k] = AN_trunc[d, y, k] * s_d[k]
#   s_d[k] = pilot_plastic_share[k, d]; same share for both years within draw k
#   (district property), independent of beta draws across k.
#-------------------------------------------------------------------------------

#--- Pre-checks on the lookup keys ---------------------------------------------
stopifnot(exists("pilot_an_draws"), exists("pilot_plastic_share"),
          identical(colnames(pilot_plastic_share), pilot_districts),
          nrow(pilot_plastic_share) == n_sim)

pilot_col_idx <- match(pilot_an_draws$District, pilot_districts)
if (anyNA(pilot_col_idx))
  stop("Districts in pilot_an_draws not found in pilot_plastic_share: ",
       paste(unique(pilot_an_draws$District[is.na(pilot_col_idx)]), collapse = ", "))
if (any(pilot_an_draws$sim < 1L | pilot_an_draws$sim > n_sim |
        pilot_an_draws$sim != round(pilot_an_draws$sim)))
  stop("`sim` in pilot_an_draws must be integers in 1..", n_sim, ".")

#--- Combine -------------------------------------------------------------------
pilot_plastic_draws <- pilot_an_draws %>%
  dplyr::mutate(
    plastic_share    = pilot_plastic_share[cbind(as.integer(sim), pilot_col_idx)],
    AN_plastic       = AN * plastic_share,
    AN_plastic_trunc = AN_trunc * plastic_share
  )
rm(pilot_col_idx)

#--- Post-checks ---------------------------------------------------------------
stopifnot(
  nrow(pilot_plastic_draws) == nrow(pilot_an_draws),
  !anyNA(pilot_plastic_draws$plastic_share),
  !anyNA(pilot_plastic_draws$AN_plastic),
  all(pilot_plastic_draws$plastic_share >= 0 & pilot_plastic_draws$plastic_share <= 1),
  all(abs(pilot_plastic_draws$AN_plastic) <= abs(pilot_plastic_draws$AN) + 1e-9),
  all(pilot_plastic_draws$AN_plastic_trunc >= 0)
)

# Same share for 2024 and 2025 within each district-draw
pilot_share_check <- pilot_plastic_draws %>%
  dplyr::group_by(District, sim) %>%
  dplyr::summarise(n_share = dplyr::n_distinct(plastic_share), .groups = "drop")
if (any(pilot_share_check$n_share != 1L))
  stop("Plastic share differs between years within a district-draw.")
rm(pilot_share_check)

#--- Summary 1: district x year ------------------------------------------------
pilot_plastic_an_summary <- pilot_plastic_draws %>%
  dplyr::group_by(District, Year, Cases) %>%
  dplyr::summarise(
    AN_mean          = mean(AN),
    share_mean       = mean(plastic_share),
    AN_plastic_mean  = mean(AN_plastic),
    AN_plastic_low   = quantile(AN_plastic, 0.025),
    AN_plastic_high  = quantile(AN_plastic, 0.975),
    AN_plastic_trunc_mean = mean(AN_plastic_trunc),
    .groups = "drop"
  ) %>%
  dplyr::mutate(AF_plastic_pct = round(100 * AN_plastic_mean / Cases, 2))
print(pilot_plastic_an_summary, n = Inf)

#--- Summary 2: pilot total by year (sum within draw, THEN quantiles) ----------
pilot_plastic_total_draws <- pilot_plastic_draws %>%
  dplyr::group_by(Year, sim) %>%
  dplyr::summarise(Cases = sum(Cases), AN = sum(AN),
                   AN_plastic = sum(AN_plastic),
                   AN_plastic_trunc = sum(AN_plastic_trunc), .groups = "drop")

pilot_plastic_total_summary <- pilot_plastic_total_draws %>%
  dplyr::group_by(Year, Cases) %>%
  dplyr::summarise(
    AN_mean         = mean(AN),
    AN_plastic_mean = mean(AN_plastic),
    AN_plastic_low  = quantile(AN_plastic, 0.025),
    AN_plastic_high = quantile(AN_plastic, 0.975),
    AN_plastic_trunc_mean = mean(AN_plastic_trunc),
    .groups = "drop"
  ) %>%
  dplyr::mutate(AF_plastic_pct = round(100 * AN_plastic_mean / Cases, 2))
print(pilot_plastic_total_summary)

message("B2c: plastic-attributable AN draws ready (", nrow(pilot_plastic_draws),
        " rows = ", length(pilot_districts), " districts x ",
        length(unique(pilot_plastic_draws$Year)), " years x ", n_sim, " draws).")


#-------------------------------------------------------------------------------
# Pilot mortality profile (national-share fallback until pilot age counts supplied)
#   Default: national age shares for year y, rescaled to pilot observed deaths.
#   Override: set `pilot_age_profile_input` BEFORE running, as a NAMED LIST
#     list(`2024` = <df>, `2025` = <df>), each with columns
#     age_group, age_mid, observed_deaths summing to the pilot deaths (14 / 16).
#-------------------------------------------------------------------------------
pilot_deaths <- tibble::tibble(Year = c(2024L, 2025L), Observed_Deaths = c(14, 16))

pilot_age_profile_override <- NULL
if (exists("pilot_age_profile_input", envir = .GlobalEnv, inherits = TRUE))
  pilot_age_profile_override <- pilot_age_profile_input

pilot_age_profiles <- stats::setNames(vector("list", 2), c("2024", "2025"))
for (yr in c(2024L, 2025L)) {
  target_deaths <- pilot_deaths$Observed_Deaths[pilot_deaths$Year == yr]
  if (is.null(pilot_age_profile_override)) {
    # Fallback: national age shape, pilot death count
    base <- mortality_profiles_by_year[[as.character(yr)]]
    if (is.null(base) || !"observed_deaths" %in% names(base) ||
        !is.finite(sum(base$observed_deaths)) || sum(base$observed_deaths) <= 0)
      stop("National mortality profile for ", yr, " is missing or has zero deaths.")
    prof <- base %>%
      dplyr::mutate(observed_deaths = observed_deaths / sum(observed_deaths) * target_deaths)
  } else {
    # Override: pilot-specific profile, used as supplied
    prof <- pilot_age_profile_override[[as.character(yr)]]
    if (is.null(prof) || !all(c("age_group", "age_mid", "observed_deaths") %in% names(prof)))
      stop("pilot_age_profile_input must be a named list with elements `2024` and `2025`, ",
           "each containing age_group, age_mid, observed_deaths.")
    if (abs(sum(prof$observed_deaths) - target_deaths) > 1e-8)
      stop("Age-profile deaths must sum to pilot observed deaths in ", yr,
           " (", target_deaths, ").")
  }
  pilot_age_profiles[[as.character(yr)]] <- prof
}
rm(yr, prof, target_deaths)
if (exists("base", inherits = FALSE)) rm(base)

pilot_age_profile_assumption <- if (is.null(pilot_age_profile_override)) {
  "Assumption: national age shares rescaled to pilot deaths (14/16); replace via pilot_age_profile_input"
} else {
  "Pilot-specific age profile supplied via pilot_age_profile_input"
}

#-------------------------------------------------------------------------------
# Pilot mortality constants
#   Reuse locked `age_consts`; only the death profile is pilot-specific.
#   All discount/horizon settings come from locked `scen`.
#-------------------------------------------------------------------------------
pilot_const_tab <- tidyr::expand_grid(
  Year = c(2024L, 2025L),
  dplyr::distinct(scen, work_hz, r_cost, r_yll)
) %>%
  dplyr::mutate(
    v = purrr::pmap(list(Year, r_yll, r_cost, work_hz), function(y, ry, rc, wh)
      age_consts(pilot_age_profiles[[as.character(y)]], ry, rc, wh)),
    yll_per_death = purrr::map_dbl(v, "yll_per_death"),
    pv_work       = purrr::map_dbl(v, "pv_work"),
    pv_life       = purrr::map_dbl(v, "pv_life")
  ) %>%
  dplyr::select(-v)

stopifnot(
  nrow(pilot_const_tab) == 2L * nrow(dplyr::distinct(scen, work_hz, r_cost, r_yll)),
  all(is.finite(pilot_const_tab$yll_per_death)),
  all(is.finite(pilot_const_tab$pv_work)),
  all(is.finite(pilot_const_tab$pv_life))
)

message(pilot_age_profile_assumption)
message("Pilot mortality constants ready: ", nrow(pilot_const_tab),
        " rows (2 years x ", nrow(pilot_const_tab) / 2, " scenario settings).")


#-------------------------------------------------------------------------------
# Pilot programme costs (LKR -> USD using locked annual FX)
#
# METHOD
#   Programme_total[y] (LKR M) = Operational_national[y] x w_pop   (top-down)
#                              + Pilot_labour                       (bottom-up)
#
#   (1) Operational: locked national operational budget for year y, allocated to
#       the pilot by population share w_pop = pilot_population_weight.
#       Assumes operational spending is proportional to population.
#   (2) Labour: annual pilot staff cost from supplied salaries and staffing.
#       Each term = monthly salary (LKR) x staff-equivalents x 12 months.
#       Part-time cadres are converted to full-time equivalents as
#       headcount x (8 days / 28 days), i.e. 8 working days per 4-week month.
#       National labour is NOT reallocated, so labour is not double-counted.
#       The same nominal LKR labour cost is applied in 2024 and 2025.
#   (3) Currency/price year: all outputs are in 2024 USD (millions).
#       2024: LKR / fx_2024 (locked 2024 annual average rate).
#       2025: lkr25_to_usd24() (locked helper: 2025 LKR -> 2024 USD).
#
#   Programme cost is reported separately and is EXCLUDED from the societal total
#   (consistent with the national engine, where it is attribution-allocated as
#   Programme_alloc = Programme_USD_M x AF).
#-------------------------------------------------------------------------------

#--- Input checks --------------------------------------------------------------
if (!exists("pilot_population_weight"))
  stop("`pilot_population_weight` (pilot share of national population) is not defined.")
stopifnot(is.numeric(pilot_population_weight), length(pilot_population_weight) == 1L,
          is.finite(pilot_population_weight),
          pilot_population_weight > 0, pilot_population_weight <= 1)
if (!all(c("Year", "Operational_LKR_M") %in% names(programme_lookup)))
  stop("programme_lookup must contain Year and Operational_LKR_M.")

#--- Pilot labour (annual LKR) -------------------------------------------------
#   Cadre 1:   8 full-time staff        x 84,435 LKR/month x 12
#   Cadre 2: 158 staff x (8/28) FTE     x 50,781 LKR/month x 12
#   Cadre 3: 667 staff x (8/28) FTE     x 37,725 LKR/month x 12
pilot_labor_components_lkr <- c(
  cadre_1 = 84435 * 8 * 12,
  cadre_2 = 50781 * ((158 * 8) / 28) * 12,
  cadre_3 = 37725 * ((667 * 8) / 28) * 12
)
pilot_labor_lkr <- sum(pilot_labor_components_lkr)   # ~121.89 LKR M per year

#--- Programme cost by year ----------------------------------------------------
pilot_programme <- programme_lookup %>%
  dplyr::filter(Year %in% c(2024, 2025)) %>%
  dplyr::transmute(
    Year = as.integer(Year),
    # (1) top-down operational allocation by population share
    Operational_Pilot_LKR_M = Operational_LKR_M * pilot_population_weight,
    # (2) bottom-up pilot labour
    Pilot_Labor_LKR_M = pilot_labor_lkr / 1e6,
    Programme_total   = Operational_Pilot_LKR_M + Pilot_Labor_LKR_M,   # LKR M
    # (3) convert to 2024 USD M
    Programme_USD_M = dplyr::if_else(
      Year == 2025,
      lkr25_to_usd24(Programme_total * 1e6) / 1e6,
      Programme_total / fx_2024)
  )

#--- Output checks -------------------------------------------------------------
if (nrow(pilot_programme) != 2L || !setequal(pilot_programme$Year, c(2024L, 2025L)))
  stop("Locked programme_lookup must have exactly one 2024 row and one 2025 row.")
if (!all(is.finite(pilot_programme$Programme_USD_M)) ||
    any(pilot_programme$Programme_USD_M <= 0))
  stop("Pilot programme USD values are non-finite or non-positive; check fx_2024 / lkr25_to_usd24().")

print(pilot_labor_components_lkr / 1e6)   # labour components, LKR M
print(pilot_programme)
message("Pilot programme costs ready (2024 USD M): ",
        paste0(pilot_programme$Year, " = ",
               round(pilot_programme$Programme_USD_M, 3), collapse = "; "))


##===============================================================================
# BRIDGE B3: pilot engine (parallel; call locked burden_engine and reuse P/scen)
#
# METHOD
#   For each scenario s and year y (2024, 2025), per draw k:
#     AN_pilot[y,k] = sum over districts of AN[d,y,k]
#                     (all-waste or plastic; signed or truncated per s$af_type)
#     burden_engine(AN_pilot, observed pilot cases, observed pilot deaths,
#                   pilot age constants, pilot programme USD, GDP per capita)
#   Combined 2024-25 (Year = 0): additive outcomes summed within draw.
#   Districts are summed within draw before any quantiles, so correlation
#   across districts and years is preserved.
#===============================================================================

# Additive outcomes summed to build the combined 2024-25 period
pilot_sum_cols <- c("AN", "Cases", "Deaths", "YLD", "YLL", "DALY", "Provider", "OOP",
                    "Indirect", "Mortality", "Societal", "Programme_total", "Programme_alloc")

run_pilot_economics <- function(x, analysis, scenario_rows, plastic = FALSE) {
  
  #--- (1) Choose the AN series: all-waste or plastic-specific -----------------
  if (plastic) {
    use_an    <- "AN_plastic"
    use_trunc <- "AN_plastic_trunc"
  } else {
    use_an    <- "AN"
    use_trunc <- "AN_trunc"
  }
  
  #--- (2) Aggregate districts to the pilot, within each Year x draw -----------
  agg <- x %>% dplyr::group_by(Year, sim) %>%
    dplyr::summarise(AN       = sum(.data[[use_an]]),
                     AN_trunc = sum(.data[[use_trunc]]),
                     Cases    = sum(Cases), .groups = "drop")
  
  #--- (3) Run the locked engine for each scenario x year ----------------------
  by_year <- purrr::map_dfr(seq_len(nrow(scenario_rows)), function(k) {
    ss <- scenario_rows[k, , drop = FALSE]          # one scenario row
    purrr::map_dfr(c(2024L, 2025L), function(yr) {
      
      # Pilot draws for this year, ordered by sim
      d <- agg %>% dplyr::filter(Year == yr) %>% dplyr::arrange(sim)
      # Checks: all draws present, in order; observed cases fixed across draws
      stopifnot(nrow(d) == n_sim, identical(d$sim, seq_len(n_sim)),
                dplyr::n_distinct(d$Cases) == 1)
      
      # Pilot age constants matching this year and the scenario's discount/horizon
      kc <- pilot_const_tab %>% dplyr::filter(Year == yr, work_hz == ss$work_hz,
                                              r_cost == ss$r_cost, r_yll == ss$r_yll)
      if (nrow(kc) != 1L) stop("Pilot age-constant lookup is not unique.")
      
      # Pilot programme cost (2024 USD M) and observed pilot deaths
      p_cost <- pilot_programme %>% dplyr::filter(Year == yr)
      obs    <- pilot_deaths$Observed_Deaths[pilot_deaths$Year == yr]
      
      # "signed" uses the main AN; any other label (e.g. "trunc") uses truncated AN
      if (ss$af_type == "signed") {
        pilot_AN <- d$AN
      } else {
        pilot_AN <- d$AN_trunc
      }
      
      # Locked national engine with pilot-specific inputs
      result <- burden_engine(AN = pilot_AN,
                              cases_tot  = d$Cases[1],
                              s          = ss,
                              deaths_obs = obs,
                              kc         = kc,
                              prog_usd_m = p_cost$Programme_USD_M,
                              gdp_pc     = lookup1(gdp_lookup, "GDP_Capita_USD", yr))
      
      # Tag year/scenario and record the full pilot programme cost
      result %>% dplyr::mutate(Year = yr, scen_id = ss$scen_id,
                               Programme_total = p_cost$Programme_USD_M)
    })
  })
  
  #--- (4) Combined 2024-25 (Year = 0) -----------------------------------------
  # Additive outcomes summed within draw; YLL per death recomputed as a ratio;
  # PV_work is not meaningful over two years, so it is NA.
  combined <- by_year %>% dplyr::group_by(scen_id, sim) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(pilot_sum_cols), sum), .groups = "drop") %>%
    dplyr::mutate(Year = 0L, YLL_per_death = YLL / Deaths, PV_work = NA_real_)
  
  #--- (5) Attach scenario labels, analysis flags and per-100k rates -----------
  out <- dplyr::bind_rows(by_year, combined) %>%
    dplyr::left_join(scenario_rows, by = "scen_id") %>%
    dplyr::mutate(Analysis = analysis, Plastic_specific = plastic,
                  AN_per100k    = AN    / pilot_population_total * 1e5,
                  Cases_per100k = Cases / pilot_population_total * 1e5)
  
  #--- (6) Additivity check: 2024 + 2025 DALYs must equal combined, per draw ---
  chk <- out %>% dplyr::group_by(scen_id, sim) %>%
    dplyr::summarise(delta = sum(DALY[Year != 0]) - sum(DALY[Year == 0]), .groups = "drop")
  stopifnot(max(abs(chk$delta)) < 1e-8)
  
  out
}

#--- Main grid (locked scen) ---------------------------------------------------
pilot_sim_all_waste <- run_pilot_economics(pilot_an_draws, "Pilot all-waste", scen, FALSE)
pilot_sim_plastic   <- run_pilot_economics(pilot_plastic_draws, "Pilot plastic-specific", scen, TRUE)

#--- Locked wage/mortality valuation sensitivity grid (exact scenario labels) --
pilot_val_all_waste <- run_pilot_economics(pilot_an_draws,
                                           "Pilot all-waste valuation sensitivity", scen_val, FALSE)
pilot_val_plastic   <- run_pilot_economics(pilot_plastic_draws,
                                           "Pilot plastic-specific valuation sensitivity", scen_val, TRUE)


#===============================================================================
# Pilot reference-level sensitivity
#
# Reference sensitivity follows locked `ref_levels` = c(0,10,25,50); each
# reference uses the same shared coefficient draws, economic P, and beta shares.
# Only the reference level changes between analyses, so differences reflect
# the counterfactual choice alone.
#===============================================================================

#--- (1) Pilot all-waste AN draws for each reference level ---------------------
# The main reference reuses the existing draws (keeps the primary result
# identical); other levels are regenerated with make_pilot_an_draws(ref).
# Names P0, P10, P25, P50 become the Analysis labels downstream.
pilot_ref_an <- stats::setNames(purrr::map(ref_levels, function(ref) {
  if (ref == cen_main) {
    pilot_an_draws
  } else {
    make_pilot_an_draws(ref)
  }
}), paste0("P", ref_levels))

#--- (2) Plastic-specific AN for each reference level --------------------------
# Draw-specific district plastic share via matrix indexing
# [row = sim, column = district position]; AN and truncated AN scaled by it.
pilot_ref_plastic <- purrr::map(pilot_ref_an, ~ .x %>%
                                  dplyr::mutate(plastic_share    = pilot_plastic_share[cbind(sim, match(District, pilot_districts))],
                                                AN_plastic       = AN       * plastic_share,
                                                AN_plastic_trunc = AN_trunc * plastic_share))

#--- (3)-(4) All-waste: run engine (primary scenario), summarise 2024-25 -------
# AF (%) per draw; mean/UI via mlu; median societal; mean programme allocation.
pilot_ref_sens_all_waste <- purrr::imap_dfr(pilot_ref_an,
                                            ~ run_pilot_economics(.x, .y, scen_primary, FALSE)) %>%
  dplyr::filter(Year == 0) %>%
  dplyr::mutate(AF = 100 * AN / Cases) %>%
  dplyr::group_by(Analysis) %>%
  dplyr::summarise(dplyr::across(c(AF, DALY, Societal), mlu),
                   Societal_Md       = stats::median(Societal),
                   Programme_alloc_M = mean(Programme_alloc),
                   .groups = "drop")

#--- (3)-(4) Plastic-specific: same pipeline using plastic AN ------------------
pilot_ref_sens_plastic <- purrr::imap_dfr(pilot_ref_plastic,
                                          ~ run_pilot_economics(.x, .y, scen_primary, TRUE)) %>%
  dplyr::filter(Year == 0) %>%
  dplyr::mutate(AF = 100 * AN / Cases) %>%
  dplyr::group_by(Analysis) %>%
  dplyr::summarise(dplyr::across(c(AF, DALY, Societal), mlu),
                   Societal_Md       = stats::median(Societal),
                   Programme_alloc_M = mean(Programme_alloc),
                   .groups = "drop")

#===============================================================================
# BRIDGE B4: pilot summaries + Table 3 (locked primary table layout)
#
# METHOD
#   Per draw k:  AF_k = 100 * AN_k / Cases
#   Per Analysis x scenario x Year: mean and UI via mlu() for
#     AN, AF, Deaths, YLD, YLL, DALY, Societal; plus median Societal.
#   Keep the locked primary scenario only.
#   Rows: Year 0 (2024-25 combined), 2024, 2025.
#   Column families and order match the locked national `tab3`.
#===============================================================================
# Three rows: combined period (Year = 0), 2024, and 2025; same seven outcome
# families, column order, and M/L/U summary columns as locked national `tab3`.

make_pilot_tab3 <- function(x) {
  tab <- x %>%
    #--- (1) Per-draw attributable fraction (%) ---------------------------------
  dplyr::mutate(AF = 100 * AN / Cases) %>%
    #--- (2) Summarise draws per analysis x scenario x period -------------------
  #    mlu() -> <var>_M, <var>_L, <var>_U ; median societal cost added
  dplyr::group_by(Analysis, scen_id, Year) %>%
    dplyr::summarise(dplyr::across(c(AN, AF, Deaths, YLD, YLL, DALY, Societal), mlu),
                     Societal_Md = stats::median(Societal), .groups = "drop") %>%
    #--- (3) Keep the locked primary scenario only ------------------------------
  dplyr::left_join(scen_primary, by = "scen_id") %>%
    dplyr::filter(primary) %>%
    #--- (4) Locked column layout: Year, then the seven outcome families --------
  #    (case-sensitive prefix match, so the order follows the locked tab3)
  dplyr::select(Year, dplyr::starts_with(
    c("AN_", "AF_", "Deaths_", "YLD_", "YLL_", "DALY_", "Societal_"),
    ignore.case = FALSE)) %>%
    #--- (5) Row order: 0 (combined), 2024, 2025 --------------------------------
  dplyr::arrange(Year)
  tab
}

pilot_tab3_all_waste <- make_pilot_tab3(pilot_sim_all_waste)
pilot_tab3_plastic   <- make_pilot_tab3(pilot_sim_plastic)

# These are the two requested primary Table-3-style outputs; the valuation and
# reference sensitivity objects are retained separately in pilot_outputs.

#--- (6) Structural check: exactly one row per period --------------------------
stopifnot(nrow(pilot_tab3_all_waste) == 3L, nrow(pilot_tab3_plastic) == 3L,
          setequal(pilot_tab3_all_waste$Year, c(0L, 2024L, 2025L)),
          setequal(pilot_tab3_plastic$Year,   c(0L, 2024L, 2025L)))

#===============================================================================
# Presentation tables (gt)
#===============================================================================
# Required gt tables. Scenario label is explicitly tied to the locked primary row.
if (!requireNamespace("gt", quietly = TRUE))
  stop("Install package `gt` to create the required presentation tables.")

# Subtitle states the primary assumptions, built directly from scen_primary
pilot_scenario_label <- with(scen_primary,
                             sprintf("Primary: signed AN; %.0f%% hospitalized; cost discount %.0f%%; YLL discount %.0f%%; work horizon %s",
                                     100 * p_hosp, 100 * r_cost, 100 * r_yll, work_hz))

make_pilot_gt <- function(x, title) {
  gt::gt(x %>% dplyr::select(Year, dplyr::everything())) %>%
    gt::tab_header(title = title, subtitle = pilot_scenario_label) %>%
    # Column labels: M = mean, L/U = uncertainty bounds; units where relevant
    gt::cols_label(Year = "Year (0 = combined)",
                   AN_M = "AN_M", AN_L = "AN_L", AN_U = "AN_U",
                   AF_M = "AF_M (%)", AF_L = "AF_L (%)", AF_U = "AF_U (%)",
                   Deaths_M = "Deaths_M", Deaths_L = "Deaths_L", Deaths_U = "Deaths_U",
                   YLD_M = "YLD_M", YLD_L = "YLD_L", YLD_U = "YLD_U",
                   YLL_M = "YLL_M", YLL_L = "YLL_L", YLL_U = "YLL_U",
                   DALY_M = "DALY_M", DALY_L = "DALY_L", DALY_U = "DALY_U",
                   Societal_M = "Societal_M (USD M)", Societal_L = "Societal_L (USD M)",
                   Societal_U = "Societal_U (USD M)", Societal_Md = "Societal_Md (USD M)") %>%
    # All numeric columns shown to 2 decimals
    gt::fmt_number(columns = dplyr::where(is.numeric), decimals = 2)
}

pilot_gt_all_waste <- make_pilot_gt(pilot_tab3_all_waste, "PABI sentinel pilot — all-waste attribution")
pilot_gt_plastic   <- make_pilot_gt(pilot_tab3_plastic,   "PABI sentinel pilot — plastic-specific attribution")


#-------------------------------------------------------------------------------
# Output objects (national table names are not overwritten)
#
# All pilot results go into one named list for saving and reloading.
# Names inside the list drop the "pilot_" prefix where the list context
# already makes the meaning clear (e.g. sim_all_waste, tab3_plastic).
#-------------------------------------------------------------------------------
pilot_outputs <- list(
  #--- Population and weighting inputs -----------------------------------------
  pilot_population                = pilot_population,
  pilot_population_total          = pilot_population_total,
  pilot_national_population_total = pilot_national_population_total,
  pilot_population_weight         = pilot_population_weight,
  
  #--- Economic inputs ---------------------------------------------------------
  pilot_labor_lkr                 = pilot_labor_lkr,
  pilot_programme                 = pilot_programme,
  
  #--- Age profiles and the assumption used to build them ----------------------
  pilot_age_profiles              = pilot_age_profiles,
  pilot_age_profile_assumption    = pilot_age_profile_assumption,
  
  #--- Plastic attribution inputs (counts, per-draw shares) --------------------
  pilot_plastic_counts            = pilot_plastic_counts,
  pilot_plastic_share_draws       = pilot_plastic_share,
  
  #--- Attributable-number draws (district x year x sim) -----------------------
  pilot_an_draws                  = pilot_an_draws,
  pilot_plastic_draws             = pilot_plastic_draws,
  
  #--- Per-draw economic simulations, main scenario grid (B3) ------------------
  sim_all_waste                   = pilot_sim_all_waste,
  sim_plastic                     = pilot_sim_plastic,
  
  #--- Primary Table-3 summaries and gt tables (B4) ----------------------------
  tab3_all_waste                  = pilot_tab3_all_waste,
  tab3_plastic                    = pilot_tab3_plastic,
  gt_all_waste                    = pilot_gt_all_waste,
  gt_plastic                      = pilot_gt_plastic,
  
  #--- Wage/mortality valuation sensitivity grid (B3) --------------------------
  valuation_all_waste             = pilot_val_all_waste,
  valuation_plastic               = pilot_val_plastic,
  
  #--- Primary scenario description (used as the gt subtitle) ------------------
  scenario_label                  = pilot_scenario_label)

#--- Save to disk (disabled). File names are tagged with the main reference ---
#    level so that runs at different P levels do not overwrite each other.
#saveRDS(pilot_outputs, file.path(output_dir, sprintf("PABI_pilot_economic_outputs_P%d.rds", cen_main)))
#gt::gtsave(pilot_gt_all_waste, file.path(output_dir, sprintf("PABI_pilot_tab3_all_waste_P%d.html", cen_main)))
#gt::gtsave(pilot_gt_plastic,   file.path(output_dir, sprintf("PABI_pilot_tab3_plastic_P%d.html",   cen_main)))

#--- Console report: full tables, with no row or column truncation -------------
print(pilot_tab3_all_waste, n = Inf, width = Inf)
print(pilot_tab3_plastic,   n = Inf, width = Inf)

# Record the age-profile assumption alongside the results in the log
message("Pilot outputs saved. Age-profile assumption: ", pilot_age_profile_assumption)


#===============================================================================
# Unified requested results table (national BI_waste vs pilot plastic attribution)
#
# METHOD
#   Primary scenario only, for both national (sim_nat) and pilot (plastic) draws.
#   For each metric and year: compute per draw, then report
#     mean (2.5th–97.5th percentile) across draws.
#   Column order: National 2024 | Pilot 2024 | National 2025 | Pilot 2025.
#   Per-capita rates use each geography's own population total.
#   Vector-control programme allocation is an alternative scenario and is
#   NOT included in societal cost.
#===============================================================================
if (!requireNamespace("gt", quietly = TRUE))
  stop("Install package `gt` to render the merged results table.")

#--- Formatting helpers --------------------------------------------------------
# Mean with 95% UI from the draw distribution: "mean (lower–upper)" + optional suffix
fmt_ui <- function(x, digits = 2, suffix = "") {
  q <- stats::quantile(x, c(.025, .975), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)%s"),
          mean(x, na.rm = TRUE), q[1], q[2], suffix)
}
# Fixed integers with thousands separators (used for population)
fmt_const <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

#--- Primary-scenario draws ----------------------------------------------------
primary_nat   <- sim_nat %>% dplyr::filter(primary)            # national, BI_waste
primary_pilot <- pilot_sim_plastic %>% dplyr::filter(primary)  # pilot, plastic-specific

#--- Context denominators ------------------------------------------------------
# National population: one value per district (stop if a district has conflicting values)
national_population <- unique(df[!is.na(df$District) & !is.na(df$Population) & df$Population > 0,
                                 c("District", "Population")])
if (anyDuplicated(national_population$District))
  stop("National population is not unique by district; cannot render table.")
national_pop_total  <- sum(national_population$Population)
national_district_n <- length(unique(names(store)))   # districts with fitted models
pilot_pop_total     <- pilot_population_total
pilot_district_n    <- length(pilot_districts)

#--- Row builders --------------------------------------------------------------
# One formatted row for a metric. `metric(z, which)` returns a per-draw vector;
# `which` is "national" or "pilot" so the metric can pick the right denominator.
make_measure <- function(label, section, metric, digits = 2, suffix = "") {
  cells <- list()
  for (yr in c(2024L, 2025L)) {
    n <- primary_nat   %>% dplyr::filter(Year == yr) %>% dplyr::arrange(sim)
    p <- primary_pilot %>% dplyr::filter(Year == yr) %>% dplyr::arrange(sim)
    if (!nrow(n) || !nrow(p)) stop("Missing primary results for year ", yr)
    nvec <- metric(n, "national")
    pvec <- metric(p, "pilot")
    # Accumulate in the requested order: 2024 national/pilot, then 2025 national/pilot.
    cells[[paste0("National_", yr)]] <- fmt_ui(nvec, digits, suffix)
    cells[[paste0("Pilot_", yr)]]    <- fmt_ui(pvec, digits, suffix)
  }
  tibble::as_tibble(c(list(Section = section, Measure = label), cells), .name_repair = "minimal")
}

# Fixed descriptor row: the same national/pilot text in both years
context_row <- function(label, section, national, pilot) {
  tibble::tibble(Section = section, Measure = label,
                 National_2024 = national, Pilot_2024 = pilot,
                 National_2025 = national, Pilot_2025 = pilot)
}

# Thin wrapper around make_measure (kept for readability of the row list)
metric_row <- function(label, section, fun, digits = 2, suffix = "") {
  make_measure(label, section, function(z, which) fun(z, which), digits, suffix)
}

#--- Assemble rows -------------------------------------------------------------
rows <- dplyr::bind_rows(
  #--- Context ---
  tibble::tibble(Section = "Context", Measure = "Year",
                 National_2024 = "2024", Pilot_2024 = "2024",
                 National_2025 = "2025", Pilot_2025 = "2025"),
  context_row("Analytical framework", "Context", "BI_waste attribution", "Plastic attribution"),
  context_row("Health districts (n)", "Context", as.character(national_district_n), as.character(pilot_district_n)),
  context_row("Population*", "Context", fmt_const(national_pop_total), fmt_const(pilot_pop_total)),
  
  #--- Epidemiological (AF computed per draw) ---
  metric_row("Total reported cases", "Epidemiological", function(z, w) z$Cases, 0),
  metric_row("Attributable fraction (AF)", "Epidemiological", function(z, w) 100 * z$AN / z$Cases, 2, "%"),
  metric_row("Attributable cases", "Epidemiological", function(z, w) z$AN, 0),
  metric_row("DALYs", "Epidemiological", function(z, w) z$DALY, 1),
  metric_row("Deaths (n)", "Epidemiological", function(z, w) z$Deaths, 2),
  
  #--- Per-capita: each geography uses its own population denominator ---
  metric_row("Cases", "Per-capita (per 100k)", function(z, w) z$Cases / (if (w == "national") national_pop_total else pilot_pop_total) * 1e5, 1),
  metric_row("DALYs", "Per-capita (per 100k)", function(z, w) z$DALY / (if (w == "national") national_pop_total else pilot_pop_total) * 1e5, 2),
  metric_row("Deaths", "Per-capita (per 100k)", function(z, w) z$Deaths / (if (w == "national") national_pop_total else pilot_pop_total) * 1e5, 3),
  
  #--- Societal cost components and total (USD M) ---
  metric_row("Provider treatment cost", "Societal costs attributable to waste-associated dengue (USD M)", function(z, w) z$Provider, 2),
  metric_row("Household OOP", "Societal costs attributable to waste-associated dengue (USD M)", function(z, w) z$OOP, 2),
  metric_row("Indirect morbidity productivity loss", "Societal costs attributable to waste-associated dengue (USD M)", function(z, w) z$Indirect, 2),
  metric_row("Premature mortality productivity loss (HCA)", "Societal costs attributable to waste-associated dengue (USD M)", function(z, w) z$Mortality, 2),
  metric_row("Total societal costs", "Societal costs attributable to waste-associated dengue (USD M)", function(z, w) z$Societal, 2),
  
  #--- Alternative scenario: AF-allocated vector-control resources ---
  metric_row("Vector-control resources corresponding to the attributable fraction†", "Alternative vector-control resource-allocation scenario", function(z, w) z$Programme_alloc, 2),
  metric_row("Total including vector-control programme resources†", "Alternative vector-control resource-allocation scenario", function(z, w) z$Societal + z$Programme_alloc, 2)
) %>%
  # Fix the section display order
  dplyr::mutate(Section = factor(Section, levels = c(
    "Context", "Epidemiological", "Per-capita (per 100k)",
    "Societal costs attributable to waste-associated dengue (USD M)",
    "Alternative vector-control resource-allocation scenario"))) %>%
  dplyr::arrange(Section)

#--- gt presentation table -----------------------------------------------------
merged_gt <- gt::gt(rows, rowname_col = "Measure", groupname_col = "Section") %>%
  gt::tab_header(title = "Dengue burden and societal cost estimates",
                 subtitle = "National BI_waste attribution and pilot-district plastic attribution; 2024 and 2025") %>%
  gt::cols_label(National_2024 = "National 2024 — Estimate (95% UI)",
                 Pilot_2024 = "Pilot districts 2024 — Estimate (95% UI)",
                 National_2025 = "National 2025 — Estimate (95% UI)",
                 Pilot_2025 = "Pilot 2025 — Estimate (95% UI)") %>%
  # Bold the two total rows
  gt::tab_style(style = gt::cell_text(weight = "bold"),
                locations = gt::cells_stub(rows = Measure %in% c("Total societal costs",
                                                                 "Total including vector-control programme resources†"))) %>%
  gt::tab_source_note(gt::md("* Population: source-script population inputs. † Vector-control is an alternative resource-allocation scenario; it is not included in societal costs."))

#--- Save: HTML always; PNG attempted (needs webshot2), failures silenced -----
#gt::gtsave(merged_gt, file.path(output_dir, "ful_ecominic_evalaution_results.html"))
#try(gt::gtsave(merged_gt, file.path(output_dir, "ful_ecominic_evalaution_results.png")), silent = TRUE)
message("Merged results table written to ", output_dir)


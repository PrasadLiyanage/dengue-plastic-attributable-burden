#===============================================================================
#===============================================================================
# STEP 1. PILOT DISTRICT DEMOGRAPHIC AND SPATIAL WEIGHTING FRAMEWORK
#===============================================================================
#
# This section defines the demographic and spatial framework for the pilot
# plastic-specific economic burden analysis.
#
# The pilot analysis evaluates the contribution of plastic-associated breeding
# habitats using detailed entomological partitioning data collected from
# selected districts.
#
# Unlike the national analysis, which uses aggregated waste-related exposure
# indicators, this module incorporates material-specific information to
# estimate the proportion of dengue burden potentially attributable to
# plastic-associated Aedes breeding habitats.
#
#------------------------------------------------------------------------------
# 1. Pilot district selection
#------------------------------------------------------------------------------
#
# The analysis includes eight sentinel districts representing diverse
# environmental, climatic, and urbanisation contexts:
#
#   • Jaffna
#   • Monaragala
#   • Kandy
#   • Galle
#   • Ratnapura
#   • Gampaha
#   • Kalutara
#   • Kurunegala

#------------------------------------------------------------------------------
# 2. Population denominator
#------------------------------------------------------------------------------
#
# District population estimates are used to:
#
#   • describe the population represented by the pilot dataset
#   • calculate population-weighted burden estimates
#   • scale pilot findings relative to the national population denominator
#
# No individual-level demographic information is included.
#
#
#------------------------------------------------------------------------------
# 3. Spatial population weighting
#------------------------------------------------------------------------------
#
# The pilot districts collectively represent a proportion of the national
# population.
#
# The population coverage fraction is calculated as:
#
#       Spatial population weight =
#
#       Population represented by pilot districts /
#       Total national population
#
#
# This weighting factor allows interpretation of pilot district estimates
# within the broader national context while preserving the distinction between
# pilot observations and national surveillance estimates.
#
#
#------------------------------------------------------------------------------
# 4. Purpose within the plastic partitioning model
#------------------------------------------------------------------------------
#
# The spatial weighting framework provides the foundation for integrating:
#
#   • plastic-specific entomological indicators
#   • district-level exposure estimates
#   • attributable dengue burden
#   • economic consequences
#
# The pilot model is intended as a validation and scaling framework for the
# plastic-specific indicator rather than a replacement for the national
# waste-based exposure model.
#
#
#------------------------------------------------------------------------------
# 5. Reproducibility and data protection
#------------------------------------------------------------------------------
#
# This section contains only aggregated district-level parameters.
#
# Confidential surveillance records and individual household/container-level
# observations are not included in the public repository.
#
# The framework can be reproduced using the provided aggregated inputs and
# adapted to other settings where material-specific entomological surveillance
# data are available.
#
#===============================================================================
#===============================================================================

library(tidyverse)
library(MASS)
library(dlnm)
library(gt)
library(scales)

set.seed(123)
n_sim <- 5000

#------------------------------------------------------------------------------
# STEP 1: DEMOGRAPHIC & SPATIAL LOOKUPS
#------------------------------------------------------------------------------
national_population_denominator <- 21781800

pilot_districts <- c(
  "Jaffna", "Monaragala", "Kandy", "Galle",
  "Ratnapura", "Gampaha", "Kalutara", "Kurunegala"
)

master_lookup <- tibble(
  District = c("Gampaha", "Kalutara", "Kandy", "Galle", "Jaffna", "Monaragala", "Ratnapura", "Kurunegala"),
  Population = c(2436142, 1305784, 1461895, 1097372, 594751, 527585, 1145423, 1768156)
)
pilot_population_total <- sum(master_lookup$Population) # 10,337,108
spatial_pop_weight <- pilot_population_total / national_population_denominator # ~0.474575

#------------------------------------------------------------------------------
#===============================================================================
# STEP 2. AGE-STRUCTURED MORTALITY AND CLINICAL SEVERITY FRAMEWORK
#===============================================================================
#
# This section defines the temporal demographic and clinical structure used for
# mortality valuation and disability-adjusted life year (DALY) estimation in
# the pilot plastic-specific economic burden model.
#
# The framework integrates:
#
#   1. Age-specific dengue mortality profiles
#   2. Annual mortality variation
#   3. Age-dependent life expectancy calculations
#   4. Clinical severity distributions
#
# These aggregated parameters allow attributable dengue deaths to be translated
# into years of life lost (YLL) and economic productivity losses.
#
#------------------------------------------------------------------------------
# 1. Age-specific mortality profiles
#------------------------------------------------------------------------------
#
# Mortality distributions are specified separately for each study year
# (2024 and 2025) to account for temporal variation in dengue mortality
# patterns.
#
# Each demographic matrix contains:
#
#   age_group:
#       Standard epidemiological age categories.
#
#   age_mid:
#       Midpoint age assigned to each age category, used for:
#          • remaining life expectancy calculations
#          • productivity loss estimation
#          • discounted mortality valuation
#
#   observed_deaths:
#       Number of observed dengue deaths distributed across age groups.
#
# The age-specific mortality profile is used to allocate simulated attributable
# deaths proportionally across age groups:
#
#       Age-specific attributable deaths =
#
#       Total attributable deaths ×
#       (observed deaths in age group / total observed deaths)
#
#
#------------------------------------------------------------------------------
# 2. Annual mortality calibration
#------------------------------------------------------------------------------
#
# Total observed dengue deaths are provided separately for each study year.
#
# These annual totals are used to:
#
#   • calibrate mortality simulations
#   • generate stochastic death counts using Poisson sampling
#   • preserve observed temporal mortality patterns
#
# Mortality uncertainty is incorporated during Monte Carlo simulation rather
# than treating observed deaths as fixed values.
#
#
#------------------------------------------------------------------------------
# 3. Clinical severity stratification
#------------------------------------------------------------------------------
#
# Attributable dengue cases are separated into:
#
#   • hospitalized cases
#   • ambulatory cases
#
# Hospitalized cases are further stratified by:
#
#   • paediatric population
#   • adult population
#   • clinical severity category
#
#
# The severity categories include:
#
#   DF:
#       Dengue fever
#
#   DHF:
#       Dengue haemorrhagic fever
#
#   DSS:
#       Dengue shock syndrome / severe dengue requiring intensive management
#
#
# Paediatric and adult severity distributions are applied as follows:
#
#   Paediatric hospitalized cases:
#       DF, DHF, DSS proportions
#
#   Adult hospitalized cases:
#       DF, DHF, DSS proportions
#
#
# These clinical proportions determine:
#
#   • healthcare utilisation patterns
#   • provider cost estimates
#   • disability duration assumptions
#   • disease severity weighting
#
#
#------------------------------------------------------------------------------
# 4. Role in economic burden estimation
#------------------------------------------------------------------------------
#
# The demographic and clinical framework provides the link between
# epidemiological attribution and economic outcomes.
#
# It enables estimation of:
#
#   Epidemiological outcomes:
#
#      • attributable deaths
#      • years of life lost (YLL)
#      • disability-adjusted life years (DALYs)
#
#
#   Economic outcomes:
#
#      • mortality-related productivity losses
#      • healthcare costs by severity category
#      • household and societal economic burden
#
#
#------------------------------------------------------------------------------
# 5. Data confidentiality and reproducibility
#------------------------------------------------------------------------------
#
# Only aggregated age-specific mortality profiles and clinical distributions
# are included.
#
# No individual patient records, hospital identifiers, or confidential
# surveillance observations are required for reproduction of this component.
#
#===============================================================================
#------------------------------------------------------------------------------
study_demographics_2024 <- tibble(
  age_group = c("< 1", "1 - 4", "5 - 9", "10 - 14", "15 - 19", "20 - 24", "25 - 34", "35 - 44", "45 - 54", "55 - 64", "≥ 65"),
  age_mid   = c(0.5, 2.5, 7.0, 12.0, 17.0, 22.0, 29.5, 39.5, 49.5, 59.5, 70.0),
  observed_deaths = c(0, 1, 1, 0, 2, 3, 5, 6, 3, 2, 1)
)

study_demographics_2025 <- tibble(
  age_group = c("< 1", "1 - 4", "5 - 9", "10 - 14", "15 - 19", "20 - 24", "25 - 34", "35 - 44", "45 - 54", "55 - 64", "≥ 65"),
  age_mid   = c(0.5, 2.5, 7.0, 12.0, 17.0, 22.0, 29.5, 39.5, 49.5, 59.5, 70.0),
  observed_deaths = c(0, 0, 2, 0, 0, 3, 7, 7, 3, 3, 4)
)

mortality_profiles_by_year <- list("2024" = study_demographics_2024, "2025" = study_demographics_2025)
observed_deaths_lookup     = tibble(Year = c(2024, 2025), Observed_Deaths = c(14, 16))

# Clinical Stratification
prop_paed  <- 0.203; prop_adult <- 0.797
sev_paed_df   <- 0.650; sev_paed_dhf  <- 0.257; sev_paed_dss  <- 0.093
sev_adult_df  <- 0.736; sev_adult_dhf <- 0.175; sev_adult_dss <- 0.089

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# STEP 3: FINANCIAL & UNIFIED PROGRAMMATIC BASELINES
#------------------------------------------------------------------------------
#
# Defines economic parameters used in the Monte Carlo cost-of-illness model.
#
# Components include:
#
#   1. Healthcare costs:
#      Inflation- and exchange-rate-adjusted costs for dengue management,
#      including ambulatory care, ward admission, and intensive care for
#      paediatric and adult dengue severity categories.
#
#   2. Household and productivity costs:
#      Out-of-pocket expenditure and indirect productivity losses are modeled
#      using log-normal distributions to account for right-skewed cost data.
#
#   3. Disability valuation:
#      Disability weights and illness duration parameters are used for DALY
#      estimation (YLD component).
#
#   4. Mortality valuation:
#      GDP per capita is used for human-capital valuation of productivity
#      losses due to premature mortality.
#
#   5. Vector-control programme costs:
#      Operational costs are population-weighted to the eight pilot districts,
#      while human-resource costs are calculated specifically for the pilot
#      surveillance area.
#
#   6. Cost uncertainty:
#      Healthcare costs are sampled during Monte Carlo simulation using gamma
#      distributions parameterized by mean cost and coefficient of variation.
#
# All monetary values are standardized to 2024 USD.
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
usd_rate <- 300.00
ex_2012  <- 131.0; ex_2018 <- 162.5; ex_2024 <- 300.9
# Inflation multipliers for converting historical costs to 2024 values
inf_12_24 <- 2.67; inf_18_24 <- 2.10

#------------------------------------------------------------------------------
# Inflation-adjusted healthcare provider costs
#------------------------------------------------------------------------------

c_ward_p_df_mean  <- (215 * ex_2012 * inf_12_24) / ex_2024
c_ward_p_dhf_mean <- (293 * ex_2012 * inf_12_24) / ex_2024
c_icu_p_dhf_mean  <- (609 * ex_2012 * inf_12_24) / ex_2024
c_ward_a_df_mean  <- (196 * ex_2012 * inf_12_24) / ex_2024
c_ward_a_dhf_mean <- (255 * ex_2012 * inf_12_24) / ex_2024
c_icu_a_dhf_mean  <- (866 * ex_2012 * inf_12_24) / ex_2024
c_ambulatory_mean <- (84  * ex_2012 * inf_12_24) / ex_2024

# Coefficient of variation for healthcare cost uncertainty
provider_cv        <- 0.20

#------------------------------------------------------------------------------
# Household out-of-pocket and productivity loss parameters
#------------------------------------------------------------------------------

oop_m <- (16.29 * ex_2018 * inf_18_24) / ex_2024
# Log-normal parameters for household costs

oop_mu <- log(oop_m) - 0.5 * log(1 + (14.02^2 / 16.29^2))
oop_sigma <- sqrt(log(1 + (14.02^2 / 16.29^2)))

# Adult productivity loss

ind_a <- (303.99 * ex_2018 * inf_18_24) / ex_2024
ind_adult_mu <- log(ind_a) - 0.5 * log(1 + (569.77^2 / 303.99^2))
ind_adult_sigma <- sqrt(log(1 + (569.77^2 / 303.99^2)))

# Paediatric caregiver productivity loss range
ind_p_min <- 15.0; ind_p_max <- 33.0

#------------------------------------------------------------------------------
# DALY parameters
#------------------------------------------------------------------------------

daly_ambulatory_mean   <- 0.032; daly_ambulatory_sd   <- 0.006
daly_hospitalized_mean <- 0.036; daly_hospitalized_sd <- 0.007

# Average dengue illness duration (years)
dengue_duration_years  <- 14 / 365.25

#------------------------------------------------------------------------------
# Macroeconomic parameters for mortality valuation
#------------------------------------------------------------------------------

gdp_per_capita_usd_2024 <- 4516.00

# GDP uncertainty parameter
gdp_cv<- 0.10


#------------------------------------------------------------------------------
# Pilot-specific vector-control programme costing
#------------------------------------------------------------------------------

# Human resource costs are estimated for the eight pilot districts.
# Salaries are annualized and converted to million USD.

labor_pilot_annual_lkr <- (84435 * 8 * 12) + (50781 * (158 * 8 / 28) * 12) + (37725 * (667 * 8 / 28) * 12)
labor_pilot_usd_m      <- (labor_pilot_annual_lkr / 1e6) / usd_rate


# Operational programme expenditure is scaled according to
# the population represented by the pilot districts.

programme_lookup <- tibble(
  Year = c(2024, 2025),
  Operational_LKR_M = c(45.68 + 9.45 + (5.15 * 0.75), 67.01 + 18.90 + (1.22 * 0.75))
) %>%
  dplyr::mutate(
    Total_Vector_Control_USD_M = ((Operational_LKR_M / usd_rate) * spatial_pop_weight) + labor_pilot_usd_m
  )


#------------------------------------------------------------------------------
# Gamma distribution sampler for uncertain healthcare costs
#------------------------------------------------------------------------------

draw_gamma_cost <- function(mean_value, cv) {
  shape <- 1 / (cv^2); scale <- mean_value / shape
  rgamma(1, shape = shape, scale = scale)
}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# STEP 4: CLEAN FIELD DATA & BUILD ENTOMOLOGICAL BETA PRIOR FALLBACKS
#------------------------------------------------------------------------------
#
# This section prepares the pilot entomological field dataset collected
# specifically to quantify the relative contribution of plastic and non-plastic
# water-holding containers to Aedes breeding in the eight pilot districts.
#
# The field dataset contains container-level observations collected during a
# dedicated survey and is cleaned by:
#
#   • removing records with no observed entomological information,
#   • replacing missing outcome values with zero, and
#   • generating a complete analytical dataset for subsequent modelling.
#
# To improve estimation stability, empirical Beta prior parameters are derived
# from the complete pilot survey. These global priors are subsequently used as
# fallback distributions for districts or container categories with sparse or
# zero observations, reducing instability in Monte Carlo sampling while
# preserving observed district-level information.
#
# NOTE:
# The original pilot field dataset is not included in the public repository.
# Researchers wishing to reproduce the complete analysis should replace the
# placeholder input below with their own pilot field dataset collected using
# an equivalent sampling protocol.
#------------------------------------------------------------------------------

# Placeholder for pilot entomological field dataset
new_pl <- read_csv("path/to/pilot_plastic_container_survey.csv")
outcomes <- c("pl_Wet", "pl_Dry", "pl_A", "pl_B", "pl_AB", "Wet", "Dry", "A", "B", "AB")

cleaned_pl <- new_pl %>%
  filter(!if_all(all_of(outcomes), is.na)) %>%
  mutate(across(all_of(outcomes), ~replace_na(., 0)))

#===============================================================================
# Empirical Beta prior parameters derived from dedicated entomological
# partitioning surveillance conducted across the eight sentinel districts.
#==============================================================================
# The pilot survey quantified the relative contribution of plastic and
# non-plastic discarded habitats to Aedes larval productivity by recording
# Aedes-positive larval habitats according to container material category.
#
# For each district i, the plastic-associated habitat contribution was defined
# probabilistically as:
#
# omega_plastic,i =
# Plastic habitats infested with Aedes larvae /
# All waste-related habitats infested with Aedes larvae
#
# District-specific uncertainty in omega_plastic,i was represented using
# Beta distributions parameterized by the observed number of plastic and
# non-plastic Aedes-positive habitats:
#
# omega_plastic,i ~ Beta(alpha_i, beta_i)
#
# where:
#   alpha_i = number of plastic-associated Aedes-positive habitats
#   beta_i  = number of non-plastic-associated Aedes-positive habitats
#
# The resulting omega_plastic,i estimates were used to partition the
# waste-attributable burden into plastic-specific burden:
#
# AF_plastic,i = AF_waste,i × omega_plastic,i
#
# District-specific observations were retained where available, while the
# pooled empirical Beta distribution below was used as an informative
# fallback for sparse strata or districts with insufficient observations.
#
# The pooled prior represents uncertainty around the overall contribution
# of plastic-associated breeding habitats in the pilot surveillance dataset.

global_a1 <- 674   # Plastic-associated Aedes-positive habitats
global_a2 <- 351   # Non-plastic-associated Aedes-positive habitats

#===============================================================================
# GENERATE ENTOMOLOGICAL PARTITIONING PROFILES (LAPLACE SMOOTHED)
#===============================================================================
ento_partitioning_profiles <- cleaned_pl %>%
  dplyr::filter(District %in% pilot_districts) %>%
  dplyr::group_by(District) %>%
  dplyr::summarise(
    plastic_pos = sum(pl_A + pl_B + pl_AB, na.rm = TRUE),
    other_pos   = sum(A + B + AB, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Clean Bayesian Laplace Smoothing (adds 1 pseudo-count to both channels)
    beta_alpha = as.numeric(plastic_pos) + 1,
    beta_beta  = as.numeric(other_pos) + 1,
    
    # Track raw totals exclusively for sample size diagnostics
    total_observed_positive = plastic_pos + other_pos
  )

# Combine DLNM relative risk thresholds with field profiles
pred_ov <- crosspred(b_pred, coef = coef(meta_overall), vcov = vcov(meta_overall), model.link = "log", cen = 0)
rr_lookup <- data.frame(Exp_Pct = round(pred_ov$predvar), RR = pred_ov$allRRfit)

ar_data <- final_analysis_data %>%
  dplyr::filter(District %in% pilot_districts) %>%
  dplyr::group_by(District) %>%
  dplyr::mutate(
    Exp_Pct = ecdf(.data[[exposure_var]])(.data[[exposure_var]]) * 100,
    Exp_Pct_round = round(Exp_Pct)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(rr_lookup, by = c("Exp_Pct_round" = "Exp_Pct")) %>%
  dplyr::left_join(ento_partitioning_profiles %>% dplyr::select(District, beta_alpha, beta_beta), by = "District")

years <- sort(unique(ar_data$Year))
simulated_gdp_vector <- rgamma(n_sim, shape = 1/(gdp_cv^2), scale = gdp_per_capita_usd_2024 * (gdp_cv^2))
coef_sim_matrix      <- MASS::mvrnorm(n_sim, coef(meta_overall), vcov(meta_overall))
vcov_zero_identity   <- vcov(meta_overall) * 0

policy_scenarios <- expand.grid(
  Retirement = c(60, 65),
  Discount   = c("Main_3.5" = 0.035, "Undiscounted_0" = 0.000, "High_5.0" = 0.050),
  stringsAsFactors = FALSE
)

final_compiled_results <- list()

#------------------------------------------------------------------------------
#===============================================================================
# GENERATE DISTRICT-SPECIFIC ENTOMOLOGICAL PARTITIONING PROFILES
#===============================================================================
#
# This section integrates pilot entomological survey data with the DLNM-derived
# exposure–response function to estimate the fraction of attributable dengue
# burden associated with plastic breeding habitats.
#
# The workflow consists of:
#
#   1. District-specific entomological partitioning
#      Positive Aedes breeding sites are aggregated into:
#
#         • plastic-associated containers
#         • other container types
#
#   2. Bayesian smoothing
#      Laplace smoothing (adding one pseudo-count to each category) is applied
#      to generate Beta distribution parameters for every district. These
#      distributions represent uncertainty in the proportion of Aedes-positive
#      breeding occurring in plastic containers and avoid unstable estimates in
#      districts with small sample sizes.
#
#   3. DLNM integration
#      The district-specific Beta parameters are merged with the DLNM-derived
#      relative risk estimates and attributable fraction dataset, linking
#      epidemiological attribution with observed entomological productivity.
#
#   4. Monte Carlo preparation
#      Simulation objects are initialized, including:
#
#         • annual study periods
#         • GDP uncertainty distributions
#         • multivariate sampling of DLNM coefficients
#         • policy scenarios for retirement age and discount rate
#
# These inputs form the basis of the subsequent probabilistic economic burden
# simulation for plastic-associated dengue transmission.
#===============================================================================
#------------------------------------------------------------------------------
for (s in seq_len(nrow(policy_scenarios))) {
  retire_age <- policy_scenarios$Retirement[s]
  disc_rate  <- policy_scenarios$Discount[s]
  
  yearly_summary_accumulator <- list()
  
  for (yr in years) {
    data_year          <- ar_data %>% dplyr::filter(Year == yr)
    observed_deaths_yr <- observed_deaths_lookup %>% filter(Year == yr) %>% pull(Observed_Deaths)
    program_base_usd_m <- programme_lookup %>% filter(Year == yr) %>% pull(Total_Vector_Control_USD_M)
    
    active_demographics <- mortality_profiles_by_year[[as.character(yr)]]
    cohort_years_lost   <- pmax(retire_age - active_demographics$age_mid, 0)
    
    Provider_sim  <- numeric(n_sim); Program_sim   <- numeric(n_sim)
    Household_sim <- numeric(n_sim); Indirect_sim  <- numeric(n_sim); Mortality_sim <- numeric(n_sim)
    YLD_sim       <- numeric(n_sim); YLL_sim       <- numeric(n_sim)
    Cases_sim     <- numeric(n_sim); AF_sim        <- numeric(n_sim)
    Deaths_sim    <- numeric(n_sim); Societal_sim  <- numeric(n_sim)
    
    for (i in seq_len(n_sim)) {
      pred_sim_loop <- crosspred(b_pred, coef = coef_sim_matrix[i, ], vcov = vcov_zero_identity, model.link = "log", cen = 0)
      rr_lookup_sim <- data.frame(Exp_Pct = round(pred_sim_loop$predvar), RR_sim = pred_sim_loop$allRRfit)
      
      loop_ar_data <- data_year %>%
        dplyr::left_join(rr_lookup_sim, by = c("Exp_Pct_round" = "Exp_Pct")) %>%
        dplyr::mutate(
          RR_sim = dplyr::if_else(is.na(RR_sim), 1, RR_sim),
          AF_waste = dplyr::if_else(RR_sim > 1, (RR_sim - 1) / RR_sim, 0),
          
          # Row-by-row beta draw using the local entomological vector productivity weights
          pl_draw = rbeta(dplyr::n(), shape1 = beta_alpha, shape2 = beta_beta),
          AF_plastic = AF_waste * pl_draw,
          AN_plastic = AF_plastic * Cases
        )
      
      an_cases <- sum(loop_ar_data$AN_plastic, na.rm = TRUE)
      mean_af  <- ifelse(sum(loop_ar_data$Cases) > 0, an_cases / sum(loop_ar_data$Cases), 0)
      
      hosp_cases <- an_cases * 0.25
      amb_cases  <- an_cases * 0.75
      
      Cases_sim[i] <- an_cases
      AF_sim[i]    <- mean_af * 100
      
      paed_cases <- hosp_cases * prop_paed; adult_cases <- hosp_cases * prop_adult
      p_df  <- paed_cases * sev_paed_df;  p_dhf  <- paed_cases * sev_paed_dhf;  p_dss  <- paed_cases * sev_paed_dss
      a_df  <- adult_cases * sev_adult_df; a_dhf  <- adult_cases * sev_adult_dhf; a_dss  <- adult_cases * sev_adult_dss
      
      c_w_p_df  <- draw_gamma_cost(c_ward_p_df_mean, provider_cv)
      c_w_p_dhf <- draw_gamma_cost(c_ward_p_dhf_mean, provider_cv)
      c_i_p_dss <- draw_gamma_cost(c_icu_p_dhf_mean, provider_cv)
      c_w_a_df  <- draw_gamma_cost(c_ward_a_df_mean, provider_cv)
      c_w_a_dhf <- draw_gamma_cost(c_ward_a_dhf_mean, provider_cv)
      c_i_a_dss <- draw_gamma_cost(c_icu_a_dhf_mean, provider_cv)
      c_amb     <- draw_gamma_cost(c_ambulatory_mean, provider_cv)
      
      Provider_sim[i] <- (
        (amb_cases * c_amb) + (p_df * c_w_p_df) + (p_dhf * c_w_p_dhf) + (p_dss * c_i_p_dss) +
          (a_df * c_w_a_df) + (a_dhf * c_w_a_dhf) + (a_dss * c_i_a_dss)
      ) / 1e6
      
      prog_draw      <- draw_gamma_cost(program_base_usd_m, 0.20)
      Program_sim[i] <- prog_draw * mean_af
      
      Household_sim[i] <- (an_cases * rlnorm(1, oop_mu, oop_sigma)) / 1e6
      Indirect_sim[i]  <- ((adult_cases * rlnorm(1, ind_adult_mu, ind_adult_sigma)) + (paed_cases * runif(1, ind_p_min, ind_p_max))) / 1e6
      
      observed_deaths_draw <- rpois(1, observed_deaths_yr)
      attributable_deaths_drawn <- mean_af * observed_deaths_draw
      Deaths_sim[i] <- attributable_deaths_drawn
      
      cohort_death_vector <- attributable_deaths_drawn * (active_demographics$observed_deaths / observed_deaths_yr)
      
      cohort_annuity_factors <- ifelse(disc_rate == 0 | cohort_years_lost <= 0,
                                       cohort_years_lost,
                                       (1 - (1 + disc_rate)^(-cohort_years_lost)) / disc_rate)
      Mortality_sim[i] <- sum(cohort_death_vector * (simulated_gdp_vector[i] * cohort_annuity_factors)) / 1e6
      
      d_amb_weight  <- max(rnorm(1, daly_ambulatory_mean, daly_ambulatory_sd), 0)
      d_hosp_weight <- max(rnorm(1, daly_hospitalized_mean, daly_hospitalized_sd), 0)
      YLD_sim[i]    <- ((amb_cases * d_amb_weight) + (hosp_cases * d_hosp_weight)) * dengue_duration_years
      
      cohort_le_remaining    <- pmax(76.0 - active_demographics$age_mid, 0)
      cohort_yll_multipliers <- ifelse(disc_rate == 0,
                                       cohort_le_remaining,
                                       (1 - exp(-disc_rate * cohort_le_remaining)) / disc_rate)
      YLL_sim[i]             <- sum(cohort_death_vector * cohort_yll_multipliers)
      
      Societal_sim[i] <- Provider_sim[i] + Program_sim[i] + Household_sim[i] + Indirect_sim[i] + Mortality_sim[i]
    }
    
    total_daly_vector <- YLD_sim + YLL_sim
    
    yearly_summary_accumulator[[as.character(yr)]] <- tibble(
      Year = yr, Retirement_Age = retire_age, Discount_Rate = disc_rate,
      
      Cases_M   = round(mean(Cases_sim)),  
      Cases_L   = round(quantile(Cases_sim, 0.025)),  
      Cases_U   = round(quantile(Cases_sim, 0.975)),
      
      AF_M      = mean(AF_sim), AF_L = quantile(AF_sim, 0.025), AF_U = quantile(AF_sim, 0.975),
      DALY_M    = round(mean(total_daly_vector)), DALY_L = round(quantile(total_daly_vector, 0.025)), DALY_U = round(quantile(total_daly_vector, 0.975)),
      Deaths_M  = round(mean(Deaths_sim)), Deaths_L = round(quantile(Deaths_sim, 0.025)), Deaths_U = round(quantile(Deaths_sim, 0.975)),
      
      Cases_Per_100k_M  = (mean(Cases_sim) / pilot_population_total) * 100000,
      Cases_Per_100k_L  = (quantile(Cases_sim, 0.025) / pilot_population_total) * 100000,
      Cases_Per_100k_U  = (quantile(Cases_sim, 0.975) / pilot_population_total) * 100000,
      
      DALY_Per_100k_M   = (mean(total_daly_vector) / pilot_population_total) * 100000,
      DALY_Per_100k_L   = (quantile(total_daly_vector, 0.025) / pilot_population_total) * 100000,
      DALY_Per_100k_U   = (quantile(total_daly_vector, 0.975) / pilot_population_total) * 100000,
      
      Deaths_Per_100k_M = (mean(Deaths_sim) / pilot_population_total) * 100000,
      Deaths_Per_100k_L = (quantile(Deaths_sim, 0.025) / pilot_population_total) * 100000,
      Deaths_Per_100k_U = (quantile(Deaths_sim, 0.975) / pilot_population_total) * 100000,
      
      Prov_M  = mean(Provider_sim),  Prov_L = quantile(Provider_sim, 0.025),  Prov_U = quantile(Provider_sim, 0.975),
      Prog_M  = mean(Program_sim),   Prog_L = quantile(Program_sim, 0.025),   Prog_U = quantile(Program_sim, 0.975),
      OOP_M   = mean(Household_sim), OOP_L  = quantile(Household_sim, 0.025), OOP_U  = quantile(Household_sim, 0.975),
      Ind_M   = mean(Indirect_sim),  Ind_L  = quantile(Indirect_sim, 0.025),  Ind_U  = quantile(Indirect_sim, 0.975),
      Mort_M  = mean(Mortality_sim), Mort_L = quantile(Mortality_sim, 0.025), Mort_U = quantile(Mortality_sim, 0.975),
      
      Societal_M = mean(Societal_sim), Societal_L = quantile(Societal_sim, 0.025), Societal_U = quantile(Societal_sim, 0.975)
    )
  }
  final_compiled_results[[paste0("R_", retire_age, "_D_", disc_rate)]] <- bind_rows(yearly_summary_accumulator)
}


#------------------------------------------------------------------------------
# STEP 6: COMPILE FINAL SIMULATION OUTPUT
#------------------------------------------------------------------------------
#
# Combine the summarized results from all Monte Carlo simulations into a single
# output table. The resulting dataset contains annual estimates of attributable
# dengue burden, health outcomes, and economic costs for each policy scenario
# (retirement age and discount rate), together with their corresponding
# uncertainty intervals. This table forms the primary analytical output used for
# subsequent result summaries, tables, and figures.
#------------------------------------------------------------------------------
master_plastic_economic_output <- bind_rows(final_compiled_results)


# Unmanaged plastic waste and dengue burden: analytical code

This repository contains the R code for the analyses in the manuscript **"Unmanaged plastic waste and dengue burden"** (currently under peer review).

**Authors:** Prasad Liyanage, Sakunthala Janaki, D S Anoja F Dheerasinghe, Preshila Samaraweera, Jagath Amarasekera, Sampath Weerakoon, Costas Velis, Marina Treskova, Joacim Rocklöv

**Corresponding author:** Professor Joacim Rocklöv, Umeå University (joacim.rocklov@umu.se)

> The code is shared to support peer review and reproducibility. It may be updated during the review process. The final version will be archived with a permanent DOI on publication.

---

## Overview

The study extends the WHO-recommended Breteau Index to identify Aedes-positive breeding habitats formed by discarded non-biodegradable waste (BI_waste). It then quantifies the association between BI_waste and dengue incidence using routine surveillance data from 28 districts in Sri Lanka (2024–2025).

The analysis has two linked components:

1. **Exposure–response modelling and attribution.** This component estimates the dengue burden attributable to waste-associated habitats at the district and national levels. It also estimates the burden attributable to plastic-associated habitats in sentinel districts.
2. **Cost-of-illness evaluation.** This component converts the attributable burden into economic costs at the national level for BI_waste and for the sentinel districts for plastics.

---

## Repository structure

```
.
├── 1.0_NM_DLNM_attribution.R   # Two-stage DLNM and attributable burden
├── 2.0_NM_attribution_COI.R    # National and sentinel-district cost-of-illness
├── .gitignore
└── README.md
```

---

## Run order

The scripts must be run **in the same R session, in this order**. The second script uses objects created by the first and does not re-estimate them.

```r
source("1.0_NM_DLNM_attribution.R")   # Two-stage DLNM and attributable burden
source("2.0_NM_attribution_COI.R")    # Cost-of-illness evaluation
```

---

## Pipeline description

### 1. `1.0_NM_DLNM_attribution.R`: exposure–response and attribution

**First stage (district level)**
- Builds district-specific distributed lag non-linear models (DLNMs) of weekly dengue incidence on BI_waste over a 12-week lag period.
- Adjusts for rainfall, temperature and surveillance intensity.
- Before modelling, rescales exposure within each district with an empirical cumulative distribution function (ECDF) to percentiles. This makes the exposure scale comparable across districts with different absolute BI_waste levels.

**Second stage (pooling)**
- Pools district-specific cumulative exposure–response coefficients with a multivariate meta-analysis.
- Derives best linear unbiased predictions (BLUPs) for each district, which borrow strength from the pooled estimate.

**Attribution**
- Calculates attributable fractions and numbers using **backward attribution**.
- Defines the counterfactual as low within-district exposure (the reference percentile is set in the script). The script maps reference percentiles back to the original BI_waste scale and reports point estimates for each percentile.
- Obtains uncertainty intervals by Monte Carlo simulation from the coefficient distributions.
- For the eight sentinel districts, apportions the waste-attributable burden to plastic-associated habitats using the plastic share of Aedes-positive waste containers. These shares come from a material-partitioning survey conducted from March to June 2026.

**Outputs:** cumulative exposure–response curves (relative risk on a log scale), and attributable fractions and cases by district, for sentinel districts and nationally.

### 2. `2.0_NM_attribution_COI.R`: cost-of-illness evaluation

- Takes the attributable case estimates from Script 1 as inputs.
- Estimates the cost of illness at the **national level** and for the **sentinel districts**.
- Distributes cases across severity categories using the national dengue severity mix.
- Uses the following cost components:
  - direct medical costs
  - direct non-medical costs
  - indirect (productivity) costs, including caregiver costs for paediatric cases.
- Reports all costs in **2024 USD** from a **societal perspective**. Vector-control programme costs are reported separately.
- Propagates uncertainty from the attribution stage into the economic estimates.

---

## Data availability

The dengue surveillance and entomological data come from the National Dengue Control Unit, Ministry of Health, Sri Lanka. **These data are not publicly available** because of data-sharing restrictions. Researchers can request access from the National Dengue Control Unit, subject to its approval.

No surveillance data are included in this repository. Users with authorised access should place input files in a local `data/` folder, which Git ignores.

Published economic parameters (such as unit costs and GDP) are documented with their sources within the code.

---

## Software requirements

- R (version 4.3 or later recommended)
- Required packages are loaded at the top of each script. Install any missing packages with `install.packages()` before running the scripts.

---

## Citation

If you use this code, please cite the associated article once it is published. Citation details and the archive DOI will be added here.

---

## Licence

This code is released under the MIT Licence. See `LICENSE` for details.

---

## Contact

For questions about the code, please open an issue in this repository or contact the corresponding author.

# 🪸 AR Attractor–Producer 🐠

Using repeated fish assemblage surveys to evaluate how shipwreck artificial reefs alter fish abundance, diversity, and life-stage structure relative to nearby natural reefs around Koh Tao, Thailand.

This repository combines the former **AR-attractor** and **AR-producer** workflows into one reproducible analysis for the Global Reef shipwreck project.

---

## 🐠 Overview

Artificial reefs may increase local fish abundance by attracting fish from surrounding habitats, supporting recruitment and growth, or both. This project uses paired artificial-reef and natural-reef surveys to test how fish assemblages change through time following shipwreck deployment.

The analysis includes three site pairs:

- **Aow Mao:** newly deployed wreck with pre- and post-deployment surveys
- **No Name:** newly deployed wreck with pre- and post-deployment surveys
- **Sattakut:** established wreck used as a post-deployment reference pair

Abundance surveys span July 2023 to late 2025. Size-structured surveys span March 2024 to October 2025.

---

## 🎯 Objectives

- Quantify temporal changes in total fish abundance on artificial and natural reefs.
- test whether abundance trajectories differ among reef types, site pairs, and feeding guilds.
- quantify changes in Shannon diversity through time and across deployment periods.
- compare juvenile and adult abundance between artificial and natural reefs.
- test whether juvenile-to-adult composition varies among site pairs.
- evaluate species-specific juvenile and adult patterns for focal damselfishes and snappers.
- provide evidence relevant to the attraction-versus-production role of shipwreck artificial reefs.

---

## 📊 Final analyses

### Feeding-guild abundance

A negative-binomial generalized linear mixed model is used to estimate temporal abundance trajectories:

```r
Count ~ feeding_guild * type * pair * date_s + (1 | survey_id)
```

The final model uses `nbinom2` with feeding-guild-specific zero inflation. A GAMM sensitivity analysis is evaluated using a temporal holdout, but the GLMM is retained because nonlinear smoothing does not materially improve predictive performance.

### Total abundance

Survey-level total abundance is modelled as:

```r
Count ~ type * pair * date_s
```

using a negative-binomial model.

### Shannon diversity

Shannon diversity is modelled as:

```r
shannon ~ type * (date_s + period) + pair
```

using a Gaussian model with reef-type-specific residual dispersion.

### Community life-stage structure

Juvenile and adult counts are reconstructed across all survey × species × life-stage combinations, including genuine zeroes, and modelled as:

```r
stage_count ~ life_stage * type * pair + (1 | Species) + (1 | survey_id)
```

using `nbinom2` with reef-type-specific dispersion.

### Species-specific life-stage models

Five focal species are modelled separately using:

```r
stage_count ~ life_stage * type * pair + (1 | survey_id)
```

The focal taxa are:

- *Neopomacentrus cyanomos*
- *Pomacentrus alexanderae*
- *Lutjanus russellii*
- *Lutjanus vitta*
- *Lutjanus argentimaculatus*

---

## 📁 Repository structure

```text
AR-attractor-producer/
├── 00_SETUP.R              # metadata, packages, paths, shared lookups and helpers
├── 01A_CLEAN_abund.R       # abundance-data cleaning and survey aggregation
├── 01B_CLEAN_size.R        # size-bin cleaning and probabilistic life stages
├── 02_EFFORT.R             # survey effort and dataset summaries
├── 03_EXPLORE.R            # formal pre-model exploration and data checks
├── 04_MODELS.R             # final models, diagnostics and sensitivity analyses
├── 05_RESULTS.R            # emmeans, contrasts and manuscript result tables
├── 06_PLOTS.R              # final manuscript figures
├── 99_scratch.R            # temporary checks and development code
├── size-fixes/             # size-bin corrections and maturity lookup
├── data/
│   ├── raw/                # source datasets
│   ├── interim/            # optional intermediate data
│   └── processed/          # cleaned analysis datasets
├── docs/                   # manuscript and project documentation
└── outputs/
    └── Analysis_YYYY.MM.DD/
        ├── diagnostics/    # model summaries and session information
        ├── exploration/    # exploration report and survey timeline
        ├── figures/        # final manuscript figures
        ├── model_objects/  # optional retained model objects
        └── tables/         # effort and final results tables
```

Only final cleaned datasets, manuscript figures, model summaries, and publication-relevant tables are written to disk. Temporary candidate models and intermediate calculations remain in the R session unless explicitly retained.

---

## ⚙️ Running the analysis

1. Clone the repository.

   ```bash
   git clone https://github.com/global-reef/AR-attractor-producer.git
   ```

2. Open the R project in RStudio.

3. Run scripts in order:

   ```text
   00_SETUP.R
   01A_CLEAN_abund.R
   01B_CLEAN_size.R
   02_EFFORT.R
   03_EXPLORE.R
   04_MODELS.R
   05_RESULTS.R
   06_PLOTS.R
   ```

4. Outputs are written to a date-stamped folder defined by `analysis_date` in `00_SETUP.R`.

The scripts are designed to be run within the same R session because downstream scripts use objects created earlier in the workflow.

---

## 📦 Main dependencies

- `dplyr`, `tidyr`, `readr`, `stringr`, `forcats`, `lubridate`, `purrr`, `tibble`
- `glmmTMB`, `lme4`, `lmerTest`, `mgcv`
- `emmeans`, `broom.mixed`, `performance`, `DHARMa`
- `vegan`
- `ggplot2`, `patchwork`, `ggeffects`
- `here`

`00_SETUP.R` checks that all required packages are installed before the analysis begins.

---

## 📈 Main outputs

- survey effort and temporal coverage summaries
- feeding-guild and total-abundance trajectories
- GLMM-versus-GAMM temporal holdout comparison
- Shannon diversity trajectories and reef-type contrasts
- community juvenile-to-adult abundance ratios
- species-specific juvenile proportions and life-stage contrasts
- model diagnostics and temporal autocorrelation tests
- publication-ready figures and CSV result tables

---

## 📝 Notes

- Surveys were conducted by Global Reef researchers around Koh Tao, Thailand.
- Artificial and natural reefs are analysed as paired systems rather than as interchangeable replicate sites.
- Aow Mao and No Name include pre- and post-deployment observations; Sattakut is post-deployment only.
- The analyses evaluate patterns consistent with attraction, redistribution, recruitment, and juvenile enrichment, but do not directly measure fish production as new biomass.

---

## License

This repository and its data are private and are not licensed for redistribution. For collaboration inquiries, contact [scarlett@global-reef.com](mailto:scarlett@global-reef.com).

**Affiliation:** [Global Reef](https://global-reef.com), Koh Tao, Thailand

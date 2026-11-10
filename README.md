# FONR *Ericameria* Survival

Analysis and supporting materials for a field experiment examining survival of *Ericameria ericoides* seedlings at the University of California, Santa Cruz Fort Ord Natural Reserve (UCSC-FONR), California.

## Project overview

This study examines how seedling survival varies with:

- **Patch type:** Dieback Patch vs. Legacy Patch
- **Herbivore access:** Caged vs. Uncaged
- **Transplant cohort:** 2021 vs. 2022

Seedlings were transplanted into replicated plots across five study sites in long-unburned maritime chaparral.

The experiment tests whether local patch structure and herbivore access interact to influence shrub seedling survival and recruitment.

Because the two cohorts were established in different years and seasons, differences between cohorts are interpreted as combined interannual and seasonal differences rather than as isolated rainfall effects.

## Experimental design

The experiment included:

- 5 study sites
- 2 patch types
- 2 herbivore treatments
- 3 replicate plots per patch × treatment combination
- 5 seedlings per plot
- 300 seedlings per cohort

### Patch types

**Dieback Patch**  
Relatively enclosed canopy gaps associated with shrub dieback and surrounded by existing shrub cover.

**Legacy Patch**  
Broader, more open patches within the shrubland.

### Herbivore treatments

**Caged**  
Seedlings protected from herbivore access.

**Uncaged**  
Seedlings exposed to natural herbivore access.

## Repository structure

```text
fonr-ericameria-survival/
├── README.md
├── LICENSE
├── .gitignore
├── data/
├── scripts/
├── figures/
│   ├── manuscript/
│   └── supplement/
├── maps/
├── docs/
└── output/
```

Folder contents may change as the manuscript is finalized.

## Data

The main datasets contain repeated survival observations for two transplant cohorts.

### 2021 cohort

The first cohort was planted in April 2021.

Primary files include:

```text
data/erer_cohort_1.csv
data/erer_cohort_1_all_surveys.csv
```

The main analysis is restricted to the period during which the herbivore-exclusion treatment remained active.

The all-surveys dataset includes later monitoring after cages were removed.

### 2022 cohort

The second cohort was planted in January 2022.

Primary file:

```text
data/erer_cohort_2.csv
```

### Plot and spatial data

Additional files contain plot coding and site-location information, including:

```text
data/Coh1 plot codes.csv
data/Coh2 plot codes.csv
data/ERER_plot_locations.csv
```

These files are used to reconcile experimental treatments, plot identities, and study-site locations.

### Herbivore damage dataset

```text
data/survival_damage_uncaged_site_patch_cohort.csv
```

This cleaned dataset summarizes uncaged seedlings by site × patch type × cohort for the supplementary analysis of recorded herbivore damage versus survival.

Each row represents one site × patch type × cohort combination containing 15 seedlings.

## Analysis

The primary R analysis is contained in the `scripts/` directory.

The analysis includes:

### GEE survival analysis

Generalized estimating equations (GEE) with a binomial error distribution and logit link are used to evaluate seedling survival through time.

The model evaluates effects of:

- patch type
- herbivore treatment
- cohort
- interactions among these factors
- time since planting
- seasonal terms

Repeated observations are clustered by individual seedling.

Planned contrasts are reported as **odds ratios (ORs)**.

### Patch contrasts

Patch-effect contrasts compare:

```text
Legacy Patch / Dieback Patch
```

Interpretation:

- OR < 1: lower survival odds in Legacy Patches
- OR > 1: higher survival odds in Legacy Patches

### Herbivore contrasts

Herbivore-effect contrasts compare:

```text
Uncaged / Caged
```

Interpretation:

- OR < 1: lower survival odds under herbivore access

### Time-to-event sensitivity analysis

Cox proportional hazards models are used as a supplementary sensitivity analysis using a **last confirmed alive** survival definition.

These analyses are reported as hazard ratios (HRs) and are used to evaluate whether the main conclusions are robust to a time-to-event formulation.

### Herbivore damage and survival

For uncaged seedlings, percent survival is compared with percent showing recorded herbivore damage using separate simple linear regressions for the 2021 and 2022 cohorts.

Each point represents one site × patch-type combination.

Missing and uprooted seedlings showing evidence consistent with herbivore removal are included in the herbivore-damage classification.

## Figures

Current manuscript-supporting figures include:

### Main manuscript

- Study-site overview map
- Observed seedling survival trajectories

### Supplementary material

- GAM-smoothed survival trajectories
- Herbivore damage versus survival
- Detailed site maps / UAV imagery, if retained

Current R-generated figure files include:

```text
Fig_1_observed_survival_curves.png
Fig_S1_GAM_smoothed_survival_curves.png
Fig_S3_damage_vs_survival.png
```

PDF versions are also generated where appropriate.

## Reproducibility

The project is organized to use project-relative paths so that the repository can be cloned and run on another computer without requiring machine-specific drive locations.

For example:

```r
readr::read_csv("data/erer_cohort_1.csv")
```

rather than:

```r
readr::read_csv("D:/GitHub/fonr-ericameria-survival/data/erer_cohort_1.csv")
```

Required R packages are listed or loaded near the beginning of the analysis script.

The primary analysis script is maintained in:

```text
scripts/
```

## Manuscript

This repository supports an in-progress manuscript examining the effects of patch structure and herbivore access on *Ericameria ericoides* seedling survival in long-unburned maritime chaparral.

Current working title:

> **Patch structure and herbivory shape shrub seedling survival across cohorts in long-unburned maritime chaparral**

The title and target journal may change prior to submission.

## Contributors

- **Jon Detka** — analysis, data processing, visualization, UAV/GIS support, and manuscript development
- **Laurel Fox** — project development, field study design, ecological interpretation, and manuscript development

Additional contributors and affiliations can be added as the manuscript is finalized.

## Licensing

Unless otherwise noted:

- **Code** in this repository is released under the **MIT License**. See `LICENSE`.
- **Data, figures, and documentation** are made available under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**.

If this repository is made public before manuscript publication, users are encouraged to contact the authors before reusing unpublished data or analysis products.

## Citation

A formal manuscript citation and repository citation will be added following publication.

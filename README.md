# Thesis Code Overview

This repository contains the code and data pipeline for estimating heterogeneous labor-market responses to oil supply news shocks. The main empirical objects are local-projection impulse response functions (IRFs) by occupation, industry, and occupation-industry cells, plus related significance tests, plots, and variance-decomposition exercises.

## Main Entry Scripts

- `main_ind_4.jl`: Runs the four-group industry analysis. It prepares the industry panel, estimates local projections, runs significance tests, produces IRF plots, and saves the industry-level OilIntensity cross-sectional regression results.
- `main_occ.jl`: Runs the occupation-level analysis across the main outcome variants.
- `main_occ_ind.jl`: Runs the occupation-by-industry panel analysis and plots IRFs by industry group.

## Folder Structure

### `data/`

Stores the raw and intermediate input data used by the project. This includes CPS fixed-width files, macroeconomic series, oil supply news shocks, crosswalk files, and O*NET task/skill data.

Key contents:

- CPS data files used to construct labor-market panels.
- `VARdata.xlsx`, `FEDFUNDS.csv`, `INDPRO.csv`, `T10Y3M.csv`, `UNRATE.csv`, and `GDPC1.csv` for macro controls.
- `oilSupplyNewsShocks_2025M06.xlsx` for the oil supply news shock series.
- Occupation crosswalk and task-intensity files.
- `data/ONET/` contains O*NET measures such as abilities, knowledge, skills, and work activities.

### `src/`

Contains all source code. The subfolders separate the project by empirical module.

### `src/ind/`

Industry-level analysis. This folder builds industry panels, estimates local projections, runs significance tests, plots industry IRFs, merges IRF coefficients, and runs cross-sectional regressions using industry OilIntensity.

The current four-group industry workflow is mainly organized as:

- `01_data_prep_ind_4.jl`: Builds the industry panel and constructs OilShare/OilIntensity exposure measures.
- `02_lp_ind_4_estimation.jl`: Estimates local projections and runs cross-sectional OilIntensity regressions.
- `03_significance_tests_ind_4.jl`: Performs significance testing for the industry IRFs.
- `04_plots_ind_4.jl`: Produces IRF and persistence-summary plots.
- `05_merge_irf_betas_ind_4.jl`: Merges IRF coefficient outputs across outcomes.
- `06_oilshare_no_mining.jl`: Additional OilShare/OilIntensity robustness code excluding mining.

### `src/occ/`

Occupation-level analysis. This folder prepares occupation panels, estimates occupation-specific IRFs, runs significance tests, creates plots, and links occupation responses to task and O*NET measures.

The revised occupation workflow is mainly organized as:

- `01_data_prep_revised.jl`: Builds occupation-level labor-market panels.
- `02_lp_occ_estimation_revised.jl`: Estimates local projections by occupation group.
- `03_significance_tests_revised.jl`: Runs occupation-level significance tests.
- `04_plots_revised.jl`: Produces occupation-level IRF plots.
- Additional scripts analyze IRF correlations, clustering, and task/O*NET links.

### `src/occ_ind/`

Occupation-by-industry analysis. This folder builds panels at the occupation-industry level, estimates local projections, plots IRFs, and performs variance-decomposition exercises to compare the relative roles of occupation and industry heterogeneity.

Key files include:

- `01_data_prep_occ_ind.jl`: Constructs the occupation-by-industry panel.
- `02_lp_estimation_occ_ind.jl`: Estimates local projections for occupation-industry cells.
- `03_lp_plots.jl`: Plots occupation-industry IRFs.
- `05_shapley_variance_decomposition.jl` and related files: Decompose variation in responses across occupation and industry dimensions.
- `07_plot_occ_ind_4_by_occupation.jl`: Plots four-group industry results by occupation.

### `src/macro/`

Macro-level scripts used to estimate or plot responses of aggregate variables such as CPI, GDP, and industrial production.

### `src/data_processing/`

Auxiliary data-processing utilities. These scripts help prepare occupation crosswalks and wide-format output sheets used elsewhere in the project.

### `result/`

Stores generated outputs from the empirical pipeline, including IRF estimates, merged coefficient files, plots, exposure measures, and wide-format summary tables.

Important subfolders:

- `result/ind_4/`: Outputs from the four-group industry analysis, including merged IRF coefficients, OilIntensity by industry group, and OilShare by occupation group.
- `result/occ/`: Occupation-level outputs, including wide-format outcome files and figures.
- `result/occ_ind_4/`: Outputs from occupation-by-industry analyses.
- `result/macro/`: Macro-response outputs.
- `result/mapping/`: Mapping files and crosswalk outputs used to connect occupations, industries, and task measures.

## Outcome Variants

Most workflows are run for several outcome variants:

- `:hourly_rate`: log real hourly wage.
- `:hours`: log hours worked.
- `:income`: weekly income or log weekly income, depending on the module.
- `:income_share_var`: income share.
- `:inequality`: within-group 75-25 log income ratio.
- `:median`: median log income or wage.
- `:unemployment`: unemployment rate.
- `:employment`: log employment.

## Typical Workflow

1. Run the relevant main script from the repository root.
2. The script loads the data-preparation, estimation, testing, and plotting modules for that empirical level.
3. Outputs are written to the corresponding folder under `result/`.



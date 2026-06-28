# Code Map

Last updated: 2026-06-28.

## Verified Code Structure

Verified:
- `code/00_setup`: path setup and package installation helpers.
- `code/01_data_prep`: Stata and Python data-preparation pipeline.
- `code/02_empirical`: reduced-form bridge-effect regressions.
- `code/03_model`: Julia model, solver, calibration, and plotting code.
- `code/sandbox`: exploratory exposure construction.

## Setup

Verified:
- `code/00_setup/data_paths.do` defines global paths for `data`, `output`, and subfolders.
- `code/00_setup/install_stata_pkgs.do` installs Stata packages including `reghdfe`, `estout`, `coefplot`, `winsor2`, `_gwtmean`, and `geoinpoly`.
- `code/00_setup/install_python_pkgs.py` checks/installs `numpy`, `pandas`, `geopandas`, `osmnx`, and optional `pandana`.
- `code/00_setup/install_julia_pkgs.jl` activates and instantiates the Julia project.
- `Project.toml` includes Julia packages for Stata file IO, optimization, regressions, Excel/CSV/JLD2, and plotting.

Inferred:
- Stata handles data cleaning and reduced-form regressions.
- Python handles geospatial shortest-path calculation.
- Julia handles structural model solution and calibration.

Uncertain:
- Whether package installation scripts should be run automatically in production, given network and system dependency risks.

## Data Preparation

### `01_prep_cied.do`

Verified:
- Loads CIED raw files for 1998-2014.
- Keeps Qingdao observations and manufacturing firms.
- Constructs firm age, export variables, firm type, wage per worker, initial 2011 wage, industry code, county/town assignment by coordinates.
- Saves `data/cied/temp/cied_qingdao_07_14.dta`.
- Saves 2011 CIED firm sample to `data/ctsd/processed/firm_cied_qingdao.dta`.

Inferred:
- CIED supplements CTSD, especially around 2011 model/sample construction.

Uncertain:
- Whether CIED source URL in comments is sufficient citation for the paper.

### `02_prep_ctsd.do`

Verified:
- Loads CTSD raw files for 2007-2020.
- Harmonizes variable names across changing annual survey schemas.
- Cleans firm IDs/names and industry codes.
- Keeps manufacturing sectors.
- Constructs revenue, intermediate inputs, real capital by perpetual inventory, labor, wage, export, ownership type, and price-index-adjusted variables.
- Creates translog production-function inputs `q`, `k`, `l`, and `m`.
- Saves markdown estimation data to `data/regression/temp/markdown_est_07_20.dta`.
- Creates Qingdao firm panel datasets with and without social-security payments included in wages.

Inferred:
- This is the heaviest and most central cleaning script.
- The no-security branch is for robustness or an older wage definition.

Uncertain:
- A sector-observation counting block appears inconsistent after `keep ind ind_name n_firm`; confirm whether it is dead code.

### `03_prep_market_access.do`

Verified:
- Builds a model sample using CTSD + CIED around 2010/2011.
- Builds a regression sample for 2007-2020.
- Exports firm coordinates to Excel for Python shortest-path calculations.
- Prepares Qingdao town government coordinates.
- Prepares 2010 Qingdao town working-age population.
- Calls `03_market_access_func.py`.
- Merges with-bridge and no-bridge travel-time matrices.
- Computes `dma = dzj - dzj_prime` and `dln_ma = ln(dma)`.
- Saves regression data to `data/regression/processed/regression_qingdao_07_20.dta`.

Inferred:
- `dzj` is no-bridge travel time and `dzj_prime` is with-bridge travel time.
- Positive `dma` means the bridge/tunnel reduces travel time.

Uncertain:
- The comment says `using QGIS` in Step 6, but current implementation imports Python-generated CSVs.

### `03_market_access_func.py`

Verified:
- Loads or downloads the Qingdao OSM drive network.
- Applies Qingdao road-speed standards by road class.
- Reads `bridge_selected.csv` with selected OSM IDs.
- Marks bridge/tunnel segments by exact OSM ID fingerprint.
- Calculates shortest paths from firms to town governments with `pandana`.
- Produces with-bridge and no-bridge commute-time CSVs for no-security, model, and 2007-2020 samples.
- Exports `output/logs/check_blocked_segments_exact.csv` for blocked-segment verification.

Inferred:
- Removing bridge/tunnel edges is preferred over assigning extreme travel times because of `pandana` contraction-hierarchy issues.

Uncertain:
- Whether selected bridge/tunnel OSM IDs should include only water-crossing segments or also approach viaducts.

### `04_estimate_markdown.do` and `04_markdown_func.do`

Verified:
- Estimates production functions by 12 broad manufacturing sectors.
- Uses a translog production function and DLW/ACF-style GMM-IV procedure.
- Computes material-input output elasticity `theta_m_tl`, labor elasticity `theta_l_tl`, markup `mu_DLW_TL`, and labor markdown variable `md_TL`.
- Winsorizes markdown by industry-year.
- Saves firm-level markdowns and sector-level weighted markdowns.

Inferred:
- The labor markdown object is central to the labor-market-power interpretation.

Uncertain:
- The exact theoretical definition of `md_TL = (1 / mu_DLW_TL) * theta_l_tl / alpha_l` needs to be documented in the paper to avoid confusion.

## Empirical Code

### `bridge_effect.do`

Verified:
- Merges regression sample with Qingdao markdown estimates.
- Constructs `ln_employ`, `ln_md`, and `ln_wage`.
- Defines `post`, `big_ma`, `ln_w0`, `ln_age`, firm type, and town-year cluster.
- Runs employment, markdown, and wage regressions.
- Runs event-study plots and heterogeneity by initial wage.
- Exports estimates to `output/tables/bridge_effect_estimates.csv` when run.

Inferred:
- The empirical design is a DID/event-study around the bridge opening with market-access intensity and initial-wage heterogeneity.

Uncertain:
- Current `output/tables` is empty, so this script has not produced visible current tables in the present workspace state.

## Model Code

### `main.jl`

Verified:
- Orchestrates package loading, Stata scripts, and model execution.
- Loads `data/model/processed/firm_qingdao_model_10.dta`.
- Builds arrays for population `l`, observed employment `l_j_data`, travel times `d` and `d′`, and observed wage `w_j_data`.
- Estimates `η`, `θ`, and `σ` by matching model moments to `β_target`, with `α` fixed.
- Uses firm-level `ind_agg` as the nested-logit sector label `s(j)`.
- Solves baseline and counterfactual equilibria.
- Exports wages and model figures.

Inferred:
- This is the intended main structural script.

Uncertain:
- It currently expects variables not present in the current model dataset.
- Running the whole script would rerun large Stata data preparation and may be expensive.

### `functions.jl`

Verified:
- Defines `RunStata`, summary/check helpers, nested-logit worker choice, a generic convergence routine, `SolveModel`, `SolveZfromW`, `ComputeModelMoments`, and `ObjectiveFunction`.
- Current `SolveModel` uses sector-nested worker choice when `firm_sector` is supplied and keeps the one-level logit as the `σ = 1` compatibility case.

Inferred:
- This file is the authoritative current structural solver.

Uncertain:
- Whether convergence normalization `x ./ mean(x)` is the final normalization, given wage bill normalization inside the update rule.

### `calibration.jl`

Verified:
- Contains a template SMM structure with undefined `β1_data` and `β2_data`.

Inferred:
- It is a placeholder or early calibration sketch.

Uncertain:
- Whether it should remain part of the active workflow.

### `model_output.jl`

Verified:
- References tariff/trade-cost `.jld2` files that are absent from `data/model/raw`.
- Calls a different `SolveModel` signature than the one in `functions.jl`.

Inferred:
- It is legacy or copied from a trade-model project.

Uncertain:
- Whether it has any role in the current bridge paper.

## Sandbox

### `code/sandbox/draft.do`

Verified:
- Constructs `firm_exposure_qingdao.dta` using 2012-2013 Qingdao business registrations and with-bridge travel times to Huangdao/Jimo new-entry activity.
- `bridge_effect.do` uses `firm_exposure_qingdao.dta` in an auxiliary regression block.

Inferred:
- This is an exploratory control or mechanism variable, not yet integrated into the main pipeline.

Uncertain:
- Whether exposure should be promoted into `code/01_data_prep` or kept as exploratory.

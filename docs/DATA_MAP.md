# Data Map

Last updated: 2026-07-11.

The repository's `data/` directory is ignored by git but exists locally. This map records what is currently on disk and what the scripts appear to expect.

## Directory Overview

Verified:
- `data/cied`: China Industrial Enterprise Database raw/temp/processed folders.
- `data/ctsd`: China Tax Survey Database raw/temp/processed folders.
- `data/geo`: road networks, town boundaries, government coordinates, population, firm coordinates, maps, and shapefiles.
- `data/model`: model temp and processed firm data.
- `data/regression`: regression temp and processed firm panel data.
- `data/busin`: Qingdao business registration data for 2012 and 2013.

Inferred:
- `raw` contains source or near-source files.
- `temp` contains intermediate generated files.
- `processed` contains analysis-ready outputs.

Uncertain:
- Some generated outputs required by scripts are absent from the current workspace.

## Key Current Datasets

### Model sample

Verified:
- `data/model/processed/firm_qingdao_model_10.dta`
  - Current shape: 701,696 rows, 14 variables.
  - Variables currently present: `id`, `town`, `dzj`, `dzj_prime`, `pop`, `year`, `employ`, `ind_code2`, `longitude`, `latitude`, `wage_inital`, `ind_agg`, `county`, `firm_town`.
  - `town` is the worker-origin / residence town in the commute matrix; `firm_town` is the firm's coordinate-derived location town.
  - This is the active model input loaded by `code/03_model/main.jl`.
  - It contains `ind_agg`, used as the sector nest `s(j)` in the nested-logit model code.

Uncertain:
- Whether `data/model/processed/firm_qingdao_model_11.dta` should also be used in a later model robustness or validation exercise.

### Regression sample

Verified:
- `data/regression/processed/regression_qingdao_07_20.dta`
  - Current shape: 48,020 rows, 35 variables.
  - Key variables include `id`, `year`, `wage_total`, `employ`, `revenue`, `export`, `ind_code2`, `firm_type`, `longitude`, `latitude`, `wage_inital`, `age`, `county`, `town`, `dzj`, `dzj_prime`, `pop`, `dma`, `dln_ma`, `md_TL`, `md_TL_w`, `ln_employ`, `ln_md`, `ln_wage`, `post`, `big_ma`.
  - Years span 2007-2020.

Inferred:
- This is the main reduced-form regression dataset.

Uncertain:
- Current file includes generated regression variables, so it may already have been modified by `bridge_effect.do`.

### CTSD prepared Qingdao sample

Verified:
- `data/ctsd/temp/ctsd_qingdao_07_20.dta`
  - Current shape: 32,600 rows, 16 variables.
  - Key variables include `id`, `year`, `sdid`, `wage_total`, `employ`, `revenue`, `export`, `ind_code2`, `firm_type`, `longitude`, `latitude`, `wage_inital_2010`, `wage_inital_2011`, `age`, `export_intensity`, `export_bool`.

Inferred:
- This is the CTSD Qingdao manufacturing firm panel used for regression and model sample assembly.

### CIED prepared Qingdao sample

Verified:
- `data/cied/temp/cied_qingdao_07_14.dta`
  - Current shape: 27,684 rows, 16 variables.
  - Key variables include `group`, `year`, `revenue`, `export`, `employ`, `wage_total`, `county`, `longitude`, `latitude`, `age`, `export_bool`, `firm_type`, `export_intensity`, `wage_inital`, `ind_code2`, `town`.
- `data/ctsd/processed/firm_cied_qingdao.dta`
  - Current shape: 4,377 rows, 16 variables.
  - Contains 2011 CIED firms used in model/sample construction.

Inferred:
- CIED helps supplement firm observations where CTSD coverage or matching is incomplete.

### Markdown outputs

Verified:
- `data/ctsd/processed/result_markdown_est_07_20.dta`
  - Large firm-year markdown output for all relevant CTSD data.
- `data/ctsd/processed/result_markdown_est_07_20_sector.dta`
  - Current shape: 12 rows, 3 variables: `ind`, `ind_name`, `markdown_wt_j`.
- `data/ctsd/processed/result_markdown_est_qingdao_07_20.dta`
  - Current shape: 12,499 rows, 4 variables: `id`, `year`, `md_TL`, `md_TL_w`.

Inferred:
- `md_TL_w` is the winsorized markdown used in reduced-form analysis.

### Geography and market access

Verified:
- `data/geo/raw/bridge_selected.csv`
  - 18 rows of selected OSM IDs for bridge/tunnel segment blocking.
- `output/logs/check_blocked_segments_exact.csv`
  - 18 blocked OSM segments exported by Python verification.
  - Includes rows named Jiaozhou Bay Tunnel and Jiaozhou Bay Bridge.
- `data/geo/raw/osm_qingdao.graphml`
  - OSMnx graph used by Python shortest-path code.
- `data/geo/raw/osm_qingdao.gpkg`
  - 103,884 road features in current GeoPackage table.
- `data/geo/raw/govern_qingdao.xlsx`
  - 141 town government coordinate rows.
- `data/geo/processed/pop_census_qingdao_2010.xlsx`
  - 172 rows of Qingdao town population data, including `15-64岁`.
- `data/geo/processed/pop_census_qingdao_2010.dta`
  - Generated Stata version used by market-access script.
- `data/geo/temp/firm_list_qingdao_model.xlsx`
  - 7,231 firm coordinate rows for model-sample routing.
- `data/geo/raw/firm_qingdao_07_20.xlsx`
  - Firm coordinate export for routing regression sample.

Inferred:
- Population weights are 2010 working-age population by town.
- Government points serve as town-level commuting destinations.

Uncertain:
- The final commute-time CSVs referenced by scripts are not present in the current data listing.

### CTSD-CIED matching

Verified:
- `data/model/temp/match_ctsd_cied.dta`
  - Current shape: 1,093,445 rows, variables `sdid`, `group`.
- `data/regression/temp/match_cied_ctsd_07_14.dta`
  - Current shape: 1,093,445 rows, variables including `sdid`, `gqid`, `group`, `年份`, `企业匹配唯一标识码`, `组织机构代码`, `企业名称`.

Inferred:
- These files map CTSD observations to CIED panel identifiers.

Uncertain:
- The matching algorithm/source is not documented in the scripts read here.

### Business registration data

Verified:
- `data/busin/BusinRegis_QingDao_12.dta`
  - 78,450 rows; variables `newgcid`, longitude, latitude, county, year.
- `data/busin/BusinRegis_QingDao_13.dta`
  - 120,173 rows; same structure.
- `code/sandbox/draft.do` uses these to construct Huangdao/Jimo exposure variables.

Inferred:
- These are auxiliary data for entry/exposure mechanisms.

Uncertain:
- They are not part of the main orchestrated pipeline in `main.jl`.

## Missing Or Stale Generated Outputs

Verified:
- `output/tables` is currently empty.
- `output/figures` is currently empty.
- `data/model/raw` is currently empty.
- `code/03_model/model_output.jl` references absent `.jld2` files under `data/model/raw`.

Uncertain:
- Whether output files are intentionally omitted, not yet generated, or cleaned.

## Data Source Claims

Verified from local files only:
- `01_prep_cied.do` contains a CIED source URL in a comment.
- `03_market_access_func.py` downloads/loads OSM data for `Qingdao, China`.
- Road speed assignment is described in code as based on Qingdao Central Urban Area Road Network Planning.
- Historical notes mention OSM, Ma Lin transportation data, Steven Davis roads, and a GDP grid source, but these are not all active in current code.

Uncertain:
- Final paper-ready data-source citations need user confirmation and external verification.

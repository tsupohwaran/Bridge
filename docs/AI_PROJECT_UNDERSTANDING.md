# AI Project Understanding

Last updated: 2026-06-05.

This document is an AI-readable map of the project as understood from the current repository. It is not a substitute for the paper. Statements are explicitly classified as `Verified`, `Inferred`, or `Uncertain`.

## Sources Read

Verified:
- Code: `code/00_setup`, `code/01_data_prep`, `code/02_empirical`, `code/03_model`, `code/sandbox`.
- Manuscript scaffold: `manuscript/paper/main.tex`.
- Draft model notes: `manuscript/draft/ModelBridge.lyx`, `manuscript/draft/Model.lyx`, `manuscript/draft/Labor reallocation (general setup).lyx`.
- Historical notes/results: selected files in `manuscript/notes` and `manuscript/archive`.
- Current data inventory under `data/` and current output inventory under `output/`.

Uncertain:
- I did not verify external data sources online.
- I did not rerun the data pipeline or model.
- I treated historical Word/PDF notes as evidence of project thinking, not as final paper results.

## Project Topic

Verified:
- The LaTeX title is `Commuting Cost, Labor Market Power and Misallocation: Evidence from Chinese Bridge`.
- `code/03_model/main.jl` is headed `Cross-Sea Bridge, Labor Market Power and Welfare`.
- `code/02_empirical/bridge_effect.do` is headed `Effect of Cross-sea Bridge on Firms`.
- Geography and data preparation are centered on Qingdao: firm coordinates, town governments, Qingdao census population, OSM Qingdao road network, and selected bridge/tunnel segments.
- Historical notes mention Jiaozhou Bay Bridge, Jiaozhou Bay Tunnel, Red Island interchange, and a 2011-06-30 opening date.

Inferred:
- The project studies how a Qingdao cross-sea transport improvement changed commuting costs / market access and thereby affected firm employment, wages, and labor-market power.
- The economic mechanism is likely labor reallocation: lower commuting costs alter which firms workers can reach, changing firm-specific labor supply elasticities and markdowns.
- The bridge shock is probably used as both reduced-form evidence and a calibration/counterfactual input for the quantitative model.

Uncertain:
- The final paper's exact research question is not written in the LaTeX body.
- It is not yet clear whether the treatment should be described as bridge only, bridge plus tunnel, or a broader cross-bay road-access package.
- It is not yet clear which historical result tables are final.

## Empirical Design

Verified:
- `03_prep_market_access.do` computes firm-to-town travel times with and without bridge/tunnel segments, then constructs `dma = dzj - dzj_prime` and `dln_ma = ln(dma)` for the regression sample.
- In that script, `dzj_prime` comes from `commute_time_qingdao*_python.csv` and `dzj` comes from `commute_time_no_bridge_qingdao*_python.csv`.
- `bridge_effect.do` defines `post = (year > 2011)`.
- `bridge_effect.do` defines `big_ma = 0` when `dln_ma <= -1` or missing, and `big_ma = 1` otherwise.
- Outcomes used in `bridge_effect.do` include `ln_employ`, `ln_md`, and `ln_wage`.
- Main controls include age plus export status by year, firm type by year, industry by year, firm fixed effects, and year fixed effects.
- The empirical script runs event-study style interactions for 2010-2013 and DID-style `big_ma × post` regressions.
- The user confirmed that the preferred event timing treats 2011 as the treatment year, so the empirical `post` definition should be `year >= 2011`.
- The user confirmed that the treatment should be interpreted as firm-level accessibility improvement induced by bridge construction, not as a simple Huangdao / Jiaozhou Bay location indicator.

Inferred:
- The reduced-form design compares firms with larger versus smaller market-access improvements after the 2011 bridge/tunnel opening.
- Initial wage heterogeneity is central: regressions interact treatment with `ln_w0_diff`, the firm's initial wage relative to the treated-group mean.
- The model calibration moments are meant to mimic a regression of employment changes on `BigMA`, wage heterogeneity, and their interaction.
- Because SUTVA / general-equilibrium spillovers can affect less exposed firms too, the DID coefficients should be interpreted as relative changes for more affected firms compared with less affected firms, not as the absolute total effect of the bridge.
- This relative-effect limitation motivates the structural general-equilibrium model.

Uncertain:
- Local empirical scripts still need to be synchronized with the confirmed treatment timing (`year >= 2011`) and the model-side `BigMA = 1[dMA >= 0.5]` definition.
- The preferred balanced/unbalanced panel definition is unclear.
- `bridge_effect.do` uses `repeat != 5` while the overall data span is 2007-2020; this likely targets a 5-year window after filters, but the exact intended sample window needs confirmation.

## Structural Model

Verified:
- Draft notes and Julia code define workers by residence/town `z` and firms by `j`.
- Workers choose firms based on wages and commuting time.
- The current Julia model uses:
  - choice probabilities `π_zj`;
  - residence population `l_z`;
  - firm labor `l_j = sum_z π_zj l_z`;
  - individual/location-firm labor supply elasticity `ε_zj = θ(1 - π_zj)`;
  - firm elasticity `ε_j` as a worker-share-weighted average;
  - firm wage condition `w_j = α z_j l_j^(α-1) ε_j / (1 + ε_j)`;
  - markdown/labor market power object `ν_j = 1 + 1 / ε_j`.
- `main.jl` fixes `α = 0.4` and estimates or sets commuting/preference parameters `η` and `θ`.
- `main.jl` uses `SolveFirmPrimitivesFromData` to back out firm amenity `a_j` from observed employment and firm productivity `z_j` from observed wages before solving counterfactuals.
- The final paper should call `a_j` firm amenity, interpret `z` as a region/street/town, and call `ν_j = 1 + 1/ε_j` markdown with `ν_j > 1`.
- The final model should continue using the current one-level logit structure.
- Current model-side `bigMA` is fixed as `1[dMA >= 0.5]`, where `dMA` is measured in minutes.

Inferred:
- `η` is the commuting-cost elasticity in worker utility.
- `θ` controls preference dispersion / responsiveness to wage and commute differences.
- The model counterfactual is intended to move from no-bridge commuting times `d` to with-bridge times `d′`.
- The model is designed to rationalize why market-access improvements affect high- and low-wage firms differently.
- The user wants the main mechanism to be that traffic integration weakens labor-market power on average, while high-wage / large firms with large accessibility gains may experience increased labor-market power / markdown.

Uncertain:
- Draft notes include a nested structure with parameter `σ`, but this is historical relative to the current one-level final model.
- `main.jl` expects `pop`, `dzj`, `dzj_prime`, and `wage_inital`, but the current `data/model/processed/firm_qingdao_model.dta` on disk has only 9 variables and lacks these fields.
- `β_target = [-0.092242, 0.0939565]` comes from empirical results, but local empirical code may not yet be updated to reproduce those values.

## Current Repository State

Verified:
- `.gitignore` excludes `data/`, `*.dta`, `output/logs/`, and common temporary files.
- `data/` exists locally and contains large raw and intermediate datasets.
- `output/tables` and `output/figures` currently contain no files; `output/logs` contains `check_blocked_segments_exact.csv`.
- `data/model/raw` currently has no files, but `code/03_model/model_output.jl` references missing `.jld2` trade-cost files there.
- `model_output.jl` calls a `SolveModel` signature that does not match the current `functions.jl`, so it appears legacy or from another project branch.

Inferred:
- `main.jl` is a high-level orchestration script, but the full end-to-end pipeline is expensive and not in a currently completed generated-output state.
- The current repository is closer to a working research workspace than a reproducible public replication package.

Uncertain:
- Which generated files should be treated as canonical current outputs.
- Whether missing commute-time CSVs and model variables are expected because the user has not run the latest pipeline, or because the current scripts are inconsistent.

## Historical Result Notes

Verified:
- `manuscript/notes/实证结果汇总20260410.pdf` contains tables for employment, markdown, and wage DID estimates.
- `manuscript/notes/结果汇总20260111.docx` says the project switched from QGIS travel-time calculation to Python, reports high QGIS/Python correlation, and says model regression coefficients still had sign/scale mismatch at that time.
- `manuscript/archive/Bridge&LandPrice_250810.docx` records a separate land-price/background exploration and policy-timing concerns.
- `manuscript/archive/Size&Md_QD_251025.docx` records descriptive work on firm size and markdown.

Inferred:
- The project has evolved from land-price and market-access explorations toward a labor-market-power mechanism.
- Some old documents are concept notes, not current paper sections.

Uncertain:
- Whether any historical result note should be migrated into the main paper.
- Whether the current empirical findings are meant to be employment down on average but less negative / more positive for high-wage treated firms, markdown down, and wages up, as some notes suggest.

## One-Sentence Working Summary

Inferred:
- This appears to be a Qingdao cross-sea bridge/tunnel project asking whether commuting-cost reductions changed firm labor-market power and misallocation by reallocating workers across firms with different initial wages.

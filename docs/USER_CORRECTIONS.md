# User Corrections

Last updated: 2026-08-02.

Use this file to record corrections from the project owner. When this file changes, update the other AI understanding documents accordingly.

## Confirmed Corrections

- `code/01_data_prep/02_prep_ctsd.do` Step 4 should retain the deflated investment variable `invest_nomi` alongside `capital_real` in `ctsd_07_20_step4.dta`.
- `ComputeModelMoments` should estimate calibration regression moments from an appended two-period firm panel: baseline equilibrium as period 0 and counterfactual equilibrium as period 1. The target specification is the two-period DID analogue of `reghdfe lnl i1.bigMA#i1.post c.w_diff#i1.bigMA#i1.post, a(id year c.w_diff#year)`, with `post = 1` for the counterfactual period. `β₁` is the coefficient on `1.bigMA#1.post`; `β₂` is the coefficient on `c.w_diff#1.bigMA#1.post`.
- The worker utility function should include a firm-level non-pecuniary amenity. Current implementation uses additive firm shifter `a_j` in utility, normalizes `mean(a_j) = 0`, solves `a_j` from observed firm employment, and solves `z_j` from observed wages.
- Model solver calls should explicitly pass `a_j` in `vars`; do not silently default missing amenities to zero.
- In the final paper, call `a_j` firm amenity.
- The final paper model should use the nested logit distribution from `manuscript/draft/Model.lyx`. This supersedes the earlier one-level-logit preference for the final paper text. In the final paper, the nest index `s` should denote a sector, not a street/town.
- In the model code, use `ind_agg` as the firm sector nest `s(j)`.
- For now, use the same three calibration moments to back out `η`, `θ`, and the nested-logit parameter `σ`, with `α` fixed.
- In paper text, `z` should be interpreted strictly as a region, specifically street/town.
- In the paper, call `ν_j = 1 + 1 / ε_j` markdown, with `ν_j > 1`.
- The final `bigMA` threshold should be fixed at `0.5` minutes.
- `β_target` comes from empirical results. The local empirical code may not yet be updated to reproduce the corresponding target values; the user will update it later.
- For now, do not add external labor-supply-elasticity moments or other moments to discipline `θ`.
- The empirical treatment should be interpreted as firm-level accessibility improvement induced by the bridge construction fact. Being located in Huangdao / Jiaozhou Bay is not itself equivalent to `BigMA`.
- The preferred event timing treats 2011 as the treatment year. Empirical `post` should be defined as `year >= 2011`.
- SUTVA / general-equilibrium spillovers are acknowledged. The DID coefficients identify relative changes of more affected firms compared with less affected firms, not the absolute total effect of the bridge. This limitation is one reason the paper needs a general-equilibrium model.
- The preferred model mechanism is that traffic integration weakens labor-market power on average, but high-wage / large firms with large accessibility improvements may experience increased labor-market power / markdown.
- After routing, when multiple government-place names map to the same census `town`, aggregate `dzj` and `dzj_prime` to their fixed arithmetic means by `id town`. Do not use forced duplicate dropping, which can retain different source-town travel times for firms at identical coordinates.

## Pending Correction Slots

- Research question: traffic integration, labor-market power, and reallocation/misallocation; exact final wording still pending.
- Treatment definition: firm-level accessibility improvement from bridge construction, with `BigMA` based on the confirmed 0.5-minute market-access threshold; not a simple Huangdao / Jiaozhou Bay location dummy.
- Main sample:
- Preferred data sources:
- Preferred empirical specification: two-period DID moments in `ComputeModelMoments` using appended baseline/counterfactual firm data and absorbed `id`, `year`, and `c.w_diff#year` effects; empirical `post` should use `year >= 2011`.
- Calibration targets: `β_target` comes from empirical results; local empirical scripts may be updated later.
- Model version: nested logit over sector nests, following `manuscript/draft/Model.lyx`, with `s` denoting `ind_agg`. The current Julia calibration should estimate `η`, `θ`, and `σ` from the same three reduced-form moments while keeping `α` fixed.
- Results to treat as current:
- Files/scripts to treat as deprecated:
- Terminology preferences: `a_j` = firm amenity; `z` = region/street/town; `ν_j` = markdown (`>1`).

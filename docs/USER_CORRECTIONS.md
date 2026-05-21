# User Corrections

Last generated: 2026-05-07.

Use this file to record corrections from the project owner. When this file changes, update the other AI understanding documents accordingly.

## Confirmed Corrections

- `ComputeModelMoments` should estimate calibration regression moments from an appended two-period firm panel: baseline equilibrium as period 0 and counterfactual equilibrium as period 1. The target specification is the two-period DID analogue of `reghdfe lnl i1.bigMA#i1.post c.w_diff#i1.bigMA#i1.post, a(id year c.w_diff#year)`, with `post = 1` for the counterfactual period. `β₁` is the coefficient on `1.bigMA#1.post`; `β₂` is the coefficient on `c.w_diff#1.bigMA#1.post`.
- The worker utility function should include a firm-level non-pecuniary amenity. Current implementation uses additive firm shifter `a_j` in utility, normalizes `mean(a_j) = 0`, solves `a_j` from observed firm employment, and solves `z_j` from observed wages.
- Model solver calls should explicitly pass `a_j` in `vars`; do not silently default missing amenities to zero.
- In the final paper, call `a_j` firm amenity.
- The final model should continue using the current one-level logit structure.
- In paper text, `z` should be interpreted strictly as a region, specifically street/town.
- In the paper, call `ν_j = 1 + 1 / ε_j` markdown, with `ν_j > 1`.
- The final `bigMA` threshold should be fixed at `0.5` minutes.
- `β_target` comes from empirical results. The local empirical code may not yet be updated to reproduce the corresponding target values; the user will update it later.
- For now, do not add external labor-supply-elasticity moments or other moments to discipline `θ`.

## Pending Correction Slots

- Research question:
- Treatment definition:
- Main sample:
- Preferred data sources:
- Preferred empirical specification: two-period DID moments in `ComputeModelMoments` using appended baseline/counterfactual firm data and absorbed `id`, `year`, and `c.w_diff#year` effects.
- Calibration targets: `β_target` comes from empirical results; local empirical scripts may be updated later.
- Model version: current one-level logit.
- Results to treat as current:
- Files/scripts to treat as deprecated:
- Terminology preferences: `a_j` = firm amenity; `z` = region/street/town; `ν_j` = markdown (`>1`).

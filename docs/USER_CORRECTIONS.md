# User Corrections

Last updated: 2026-05-20.

Use this file to record corrections from the project owner. When this file changes, update the other AI understanding documents accordingly.

## Confirmed Corrections

- `ComputeModelMoments` should estimate calibration regression moments from an appended two-period firm panel: baseline equilibrium as period 0 and counterfactual equilibrium as period 1. The target specification is the two-period DID analogue of `reghdfe lnl i1.bigMA#i1.post c.w_diff#i1.bigMA#i1.post, a(id year c.w_diff#year)`, with `post = 1` for the counterfactual period. `β₁` is the coefficient on `1.bigMA#1.post`; `β₂` is the coefficient on `c.w_diff#1.bigMA#1.post`.
- The worker utility function should include a firm-level non-pecuniary amenity. Current implementation uses additive firm shifter `a_j` in utility, normalizes `mean(a_j) = 0`, solves `a_j` from observed firm employment, and solves `z_j` from observed wages.
- Model solver calls should explicitly pass `a_j` in `vars`; do not silently default missing amenities to zero.
- Worker utility should use commuting time in minutes directly, `-η d_zj`, rather than log commuting time, `-η ln(d_zj)`.

## Pending Correction Slots

- Research question:
- Treatment definition:
- Main sample:
- Preferred data sources:
- Preferred empirical specification: two-period DID moments in `ComputeModelMoments` using appended baseline/counterfactual firm data and absorbed `id`, `year`, and `c.w_diff#year` effects.
- Calibration targets:
- Model version:
- Results to treat as current:
- Files/scripts to treat as deprecated:
- Terminology preferences:

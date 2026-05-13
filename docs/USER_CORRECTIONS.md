# User Corrections

Last generated: 2026-05-07.

Use this file to record corrections from the project owner. When this file changes, update the other AI understanding documents accordingly.

## Confirmed Corrections

- `ComputeModelMoments` should estimate calibration regression moments from an appended two-period firm panel: baseline equilibrium as period 0 and counterfactual equilibrium as period 1. The target specification is the two-period DID analogue of `reghdfe lnl i1.bigMA#i1.post c.w_diff#i1.bigMA#i1.post, a(id year c.w_diff#year)`, with `post = 1` for the counterfactual period. `β₁` is the coefficient on `1.bigMA#1.post`; `β₂` is the coefficient on `c.w_diff#1.bigMA#1.post`.

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

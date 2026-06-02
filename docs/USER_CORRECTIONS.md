# User Corrections

Last generated: 2026-05-07.

Use this file to record corrections from the project owner. When this file changes, update the other AI understanding documents accordingly.

## Confirmed Corrections

- `ComputeModelMoments` should estimate calibration regression moments from an appended two-period firm panel: baseline equilibrium as period 0 and counterfactual equilibrium as period 1. The current target specification should match the two-regression design in `code/02_empirical/calculate_calibration_4_moments.do`: regression 1 uses the full two-year sample and reports `1.BIG#1.post` and `1.BIG#1.post#c.demean_lnw0`; regression 2 restricts to `BIG==1` and reports `1.post#c.lndma` and `1.post#c.demean_lnw0#c.lndma`. The model moment order is `[BIG#post, BIG#post#demean_lnw0, lndma#post among BIG==1, demean_lnw0#lndma#post among BIG==1]`.
- The worker utility function should include a firm-level non-pecuniary amenity. Current implementation uses additive firm shifter `a_j` in utility, normalizes `mean(a_j) = 0`, solves `a_j` from observed firm employment, and solves `z_j` from observed wages.
- Model solver calls should explicitly pass `a_j` in `vars`; do not silently default missing amenities to zero.
- In the final paper, call `a_j` firm amenity.
- The final model should continue using the current one-level logit structure.
- In paper text, `z` should be interpreted strictly as a region, specifically street/town.
- In the paper, call `ν_j = 1 + 1 / ε_j` markdown, with `ν_j > 1`.
- The final model-side `BIG` threshold should match `calculate_calibration_4_moments.do`: `BIG = 1[dMA > 0.5]`.
- `β_target` should come from the four reported empirical coefficients in `calculate_calibration_4_moments.do`.
- As of the current two-regression `calculate_calibration_4_moments.do` output, the calibration target order and values are `[BIG#post, BIG#post#demean_lnw0, lndma#post among BIG==1, demean_lnw0#lndma#post among BIG==1] = [-0.0912855, 0.1068467, -0.0302357, 0.0512688]`.
- For now, do not add external labor-supply-elasticity moments or other moments to discipline `θ`.

## Pending Correction Slots

- Research question:
- Treatment definition:
- Main sample:
- Preferred data sources:
- Preferred empirical specification: four two-period DID moments in `ComputeModelMoments` using appended baseline/counterfactual firm data. Match the current two-regression `calculate_calibration_4_moments.do` design: full-sample `BIG` and `BIG × demean_lnw0` moments, plus treated-only `lndma` and `lndma × demean_lnw0` moments.
- Calibration targets: `β_target` should contain the four reported coefficients from `calculate_calibration_4_moments.do`.
- Model version: current one-level logit.
- Results to treat as current:
- Files/scripts to treat as deprecated:
- Terminology preferences: `a_j` = firm amenity; `z` = region/street/town; `ν_j` = markdown (`>1`).

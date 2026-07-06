# Wage-Ladder Poaching GE Solver

Last updated: 2026-07-05.

This document records the solution logic for the new diagnostic model in `code/03_model/diagnose_poaching_groups.jl`. It separates what is implemented from the economic interpretation and remaining open questions.

## Status

Verified:
- The active diagnostic script is `code/03_model/diagnose_poaching_groups.jl`.
- The diagnostic extends the sector-nested logit model with baseline employer groups defined as `ind_agg × baseline wage decile`.
- The script now runs fast amenity inversion by default, with `RUN_INVERSION=1` and `FAST_AMENITY=1`.
- The script now runs both static post-bridge labor transitions and full general equilibrium wage fixed points by default, with `RUN_GE=1`.
- The default diagnostic grid is `rho ∈ {0, 0.25, 0.5, 1}` and `kappa = 0`.
- The current output table is `output/tables/poaching_group_beta_diagnostics.csv`.

Inferred:
- This extension is designed to represent a wage-ladder poaching channel: workers initially attached to lower-wage employer groups are more responsive to firms above their baseline wage group.
- The purpose is to see whether the model can generate the positive control-firm wage slope in employment changes that appears in the reduced-form moments.

Uncertain:
- The extension is not yet the final calibrated model. The parameter `rho` is diagnostic and has not yet been estimated from a formal moment condition.
- The extension uses model-implied baseline origin and employer-group masses, not observed worker-level flows.

## Core Objects

Verified:
- Residence locations are indexed by `z`.
- Firms are indexed by `j`.
- Sector nests are indexed by `s(j)` and implemented with firm-level `ind_agg`.
- Baseline employer groups are indexed by `g` and are constructed as

```text
g(j) = (ind_agg(j), wage_decile(j)).
```

- The wage decile is computed from normalized baseline log wages.
- Group-level mean log wage is employment-weighted:

```text
ell_g = sum_{j:g(j)=g} L_j_data log(w_j_data) / sum_{j:g(j)=g} L_j_data.
```

Inferred:
- `g` is an approximation to workers' baseline labor-market segment. It is not a true worker type observed in the data.

## Fast Amenity Inversion

Verified:
- The old amenity inversion used the full primitive recovery routine at each iteration and was slow.
- The fast inversion solves only the employment-matching problem during iterations and computes elasticities only once at the end.
- Let `r_j = theta * a_j`. For fixed observed wages, baseline travel times, and parameters, each iteration computes nested-logit shares from

```text
log Q_zj = log(w_j_data) + a_j - eta log(d_zj).
```

- It then computes model employment

```text
L_j(r) = sum_z l_z pi_zj(r).
```

- The update is a damped log employment correction:

```text
r_tilde_j = r_j + log(L_j_data) - log(L_j(r))
r_tilde   = r_tilde - mean(r_tilde)
r_new     = damp * r + (1 - damp) * r_tilde
r_new     = r_new - mean(r_new)
```

- Current default settings are:

```text
AMENITY_TOL     = 1e-5
AMENITY_MAXITER = 2000
AMENITY_DAMP    = 0.6
```

- In the verified default run, fast inversion converged in 779 iterations with final employment gap `9.8855e-6`.

Inferred:
- The speed gain comes from avoiding repeated computation of origin-firm elasticities inside the inversion loop. Since the inversion only needs to match employment, probabilities are enough until the final pass.

Uncertain:
- The best damping parameter may differ if the model is later estimated over a wider parameter space.

## Baseline Group Masses

Verified:
- After baseline amenities are inverted, the script computes baseline choice probabilities `pi_zj`.
- The share of firm `j`'s baseline employment coming from origin `z` is

```text
gamma_zj = pi_zj * l_z / L_j.
```

- The mass of workers in origin `z` attached to baseline employer group `g` is

```text
M_zg = sum_{j:g(j)=g} L_j_data * gamma_zj.
```

- The implementation checks that group masses and total mass are preserved. In the default run:

```text
group_mass_gap = 1.04e-17
total_mass_gap = 1.11e-16
```

Inferred:
- `M_zg` is the key state variable that lets the post-bridge choice problem distinguish workers by both residence and baseline employer segment.

Uncertain:
- Because `M_zg` is model-implied, any misspecification in baseline amenities or baseline choice probabilities affects the inferred worker-group composition.

## Post-Bridge Choice With Poaching

Verified:
- For a worker living in `z` and attached to baseline group `g`, the post-bridge systematic utility from firm `j` is

```text
log Q_zgj =
    log(w_j)
    + a_j
    - eta log(d'_zj)
    + rho * max(log(w_j) - ell_g, 0)
    - kappa * 1{ind_agg(j) != ind_g}.
```

- `rho` governs wage-ladder poaching.
- `kappa` governs an optional cross-sector penalty. The current default diagnostic sets `kappa = 0`.
- Given `Q_zgj`, the script applies the same sector-nested logit formula used in the baseline model.

Inferred:
- When `rho > 0`, firms above a worker group's baseline wage become more attractive over and above the direct log wage term.
- This makes low-wage baseline groups especially responsive to high-wage firms whose commute access improves after the bridge.

Uncertain:
- Whether the final paper should keep `kappa = 0`, estimate it, or omit it entirely remains open.

## Static Transition

Verified:
- The static transition holds wages fixed at baseline observed wages.
- For each `rho, kappa`, the script computes group-specific post-bridge choice probabilities and post labor:

```text
L'_j = sum_z sum_g M_zg * pi_zgj(d', w_data, a).
```

- It then estimates the first-difference diagnostic regression:

```text
Delta log L_j =
    gamma_L * w_diff_j
    + beta_1 * BigMA_j
    + beta_2 * BigMA_j * w_diff_j
    + error_j.
```

- The output column `beta_labor_wdiff_post` is `gamma_L`, the control-firm wage slope.

Inferred:
- Static transition isolates the worker-reallocation force before equilibrium wage feedback dampens or amplifies it.

## Full GE Solver

Verified:
- Full GE updates wages using the firm first-order condition:

```text
w_j = alpha * A_j * L_j^(alpha - 1) * epsilon_j / (1 + epsilon_j).
```

- For `rho = 0` and `kappa = 0`, the script calls the standard `SolveModel` because the model collapses to the sector-nested logit baseline.
- For `rho > 0` or `kappa > 0`, the script uses a custom log-wage fixed point.
- In each GE iteration:

```text
1. Normalize current wages.
2. Compute pi_zgj, L_j, and epsilon_j from the group-specific nested logit.
3. Compute updated wages from the firm FOC.
4. Normalize updated wages.
5. Compute the max absolute log wage update.
6. Apply a damped and step-clipped log-wage update.
```

- The own-wage derivative of the wage-ladder utility is

```text
lambda_gj = 1 + rho * 1{log(w_j) > ell_g}.
```

- Therefore the origin-group-firm elasticity is

```text
epsilon_zgj =
    lambda_gj * theta *
    [1/sigma + (1 - 1/sigma) * pi_zgj|s - pi_zgj].
```

- Firm-level elasticity is the employment-share-weighted average across origins and baseline employer groups:

```text
epsilon_j =
    sum_z sum_g (M_zg * pi_zgj / L_j) * epsilon_zgj.
```

- Current default GE settings are:

```text
GE_TOL      = 1e-4
GE_MAXITER  = 300
GE_DAMP     = 0.95
GE_LOG_STEP = 0.30
```

- The script uses rho-continuation: for a fixed `kappa`, the GE solution at one `rho` is used as the initial wage vector for the next `rho`.

Inferred:
- Rho-continuation is important for stability because the poaching term makes the wage fixed point more nonlinear.
- The full-GE slope is smaller than the static slope because wages adjust and partly offset pure reallocation pressure.

Uncertain:
- The current solver is adequate for diagnostics but may be too slow for a large calibration loop. Further optimization or a more formal nonlinear solver may be needed before estimating `rho`.

## Verified Diagnostic Results

Verified:
- The default run uses:

```text
eta = 0.9
theta = 10.0
sigma = 0.25
alpha = 0.8
wage groups = 10
kappa = 0
baseline mode = fast_inverted_amenities
```

- Fast amenity inversion converges:

```text
amenity_converged = true
amenity_iterations = 779
amenity_gap = 9.8855e-6
```

- Full GE converges for all default `rho` values:

| rho | static beta_labor_wdiff_post | full-GE beta_labor_wdiff_post | GE iterations | GE gap |
|---:|---:|---:|---:|---:|
| 0.00 | -0.405378 | -0.040654 | 239 | 9.644e-5 |
| 0.25 | 1.198043 | 0.324993 | 112 | 9.740e-5 |
| 0.50 | 1.303201 | 0.575035 | 243 | 9.924e-5 |
| 1.00 | 1.190585 | 0.870615 | 116 | 9.749e-5 |

Inferred:
- The sign problem is not solved by full GE alone: with `rho = 0`, the full-GE control-firm wage slope is still negative.
- A positive wage-ladder channel can reverse the sign, even after GE wage adjustment.

Uncertain:
- The empirical value of `rho` should be disciplined by a formal moment, most naturally the control-firm wage slope `beta_labor_wdiff_post` or an equivalent moment.

## Relation To Paper Text

Verified:
- The paper draft now documents the extension in `manuscript/paper/main.tex`.
- The working notes now document the mechanism and diagnostic results in `manuscript/notes/results_260701.tex`.

Inferred:
- The extension should be presented as a diagnostic or candidate model extension until `rho` is formally estimated.

Uncertain:
- The final paper may choose between three options:
  1. keep the current sector-nested model as the main model and present wage-ladder poaching as a diagnostic;
  2. add `rho` as a fourth structural parameter and estimate it with an additional moment;
  3. replace the sector-only nest structure with the baseline employer-group extension as the main worker-choice model.

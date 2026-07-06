# Model Notes

Last updated: 2026-07-05.

This document summarizes the structural model as read from current Julia code, draft notes, and user corrections. It separates verified code status from remaining open interpretation issues. The latest update is that the baseline sector-nested logit model is now accompanied by a diagnostic wage-ladder poaching extension, documented in `docs/WAGE_LADDER_POACHING_GE.md`.

## Current Implemented Julia Model

Verified:
- Implemented in `code/03_model/functions.jl` and called by `code/03_model/main.jl`.
- Current Julia code implements the sector-nested logit worker-choice model.
- `main.jl` reads firm sector nests from `ind_agg` in the model sample and passes them as `firm_sector`.
- The code keeps the one-level logit as the `σ = 1` compatibility case when no sector nest is supplied.
- Current code indices:
  - `z`: residence / town / worker origin location.
  - `s`: sector nest, implemented with firm-level `ind_agg`.
  - `j`: firm.
- Current code inputs:
  - `l_z`: residence population, normalized to sum to 1.
  - `d_zj`: commuting time from residence `z` to firm `j`.
  - `w_j`: firm wage.
  - `a_j`: firm-level non-pecuniary amenity, called firm amenity in the final paper, implemented as an additive utility shifter and inverted from observed firm employment.
  - `z_j`: firm productivity, backed out from observed wages.
  - `α`: decreasing returns parameter, set to `0.8` in the current `main.jl` calibration block.
  - `η`: commuting-cost elasticity.
  - `θ`: preference dispersion / responsiveness parameter.
  - `σ`: nested-logit sector correlation parameter.

Verified:
- Current Julia code corresponds to the same worker indirect utility as the paper:

```text
U_izj = ln(w_j) + a_j - η ln(d_zj) + (1/θ) ξ_ij
```

- With sector nests, `WorkerChoice` computes `π_zs`, `π_zj|s`, and `π_zj = π_zs × π_zj|s` using log-sum-exp normalization.
- Current Julia code computes:

```text
ε_zj = θ [1/σ + (1 - 1/σ) π_zj|s - π_zj]
γ_zj = π_zj l_z / l_j
ε_j = sum_z γ_zj ε_zj
ν_j = 1 + 1 / ε_j
```

Uncertain:
- The sector-only nested-logit simulations still need final paper-facing validation, especially because the wage-ladder poaching extension is now under diagnostic testing.

## Target Manuscript Model: Sector-Nested Logit

Verified:
- The user instructed on 2026-06-28: use the nested logit distribution in `manuscript/draft/Model.lyx`, but define `s` as sector.
- In the final paper, `z` should remain a residence region, specifically street/town.
- In the final paper, `a_j` should be called firm amenity.
- In the final paper, `ν_j = 1 + 1 / ε_j` should be called markdown, with `ν_j > 1`.

Verified:
- Target manuscript indices:
  - `z ∈ Z`: residence location / street / town.
  - `s ∈ S`: sector nest.
  - `J(s)`: firms in sector `s`.
  - `j ∈ J(s)`: firm `j` belonging to sector `s`.
- Worker indirect utility should keep the firm amenity term:

```text
U_izj = ln(w_j) + a_j - η ln(d_zj) + (1/θ) ξ_ij
```

- The draft nested-logit distribution is:

```text
F(ξ_i1, ..., ξ_iJ)
= exp{ - sum_{s∈S} [ sum_{j∈J(s)} exp(-ξ_ij / σ) ]^σ }.
```

Inferred:
- Sector nests capture correlation in workers' idiosyncratic preferences across firms in the same sector, consistent with sector-specific skills, search, and outside options.
- `σ = 1` collapses the model to the one-level logit compatibility case.
- Smaller `σ` implies stronger within-sector correlation in preference shocks and makes within-sector alternatives more closely related.

Uncertain:
- The empirical and economic interpretation of the `ind_agg` sector aggregation should be explained carefully in the paper.

## Worker Choice In The Target Model

Verified:
- Define the systematic attractiveness index:

```text
q_zj = w_j exp(a_j) d_zj^(-η).
```

- For a worker living in `z`, the probability of choosing sector `s` is:

```text
π_zs =
{ sum_{k∈J(s)} q_zk^(θ/σ) }^σ
/
sum_{r∈S} { sum_{k∈J(r)} q_zk^(θ/σ) }^σ.
```

- Conditional on choosing sector `s`, the probability of choosing firm `j ∈ J(s)` is:

```text
π_zj|s =
q_zj^(θ/σ)
/
sum_{k∈J(s)} q_zk^(θ/σ).
```

- The unconditional firm choice probability is:

```text
π_zj = π_zs × π_zj|s.
```

- Firm labor is:

```text
L_j = sum_z l_z π_zj.
```

Inferred:
- Lower commuting time raises `q_zj`, increasing the attractiveness of firm `j` to workers in affected origins.
- Higher firm amenity raises a firm's attractiveness independently of wages and commuting time.
- The nested structure lets a commuting shock change both sector choice probabilities `π_zs` and within-sector firm choice probabilities `π_zj|s`.

## Labor Supply Elasticity And Markdown

Verified:
- For `j ∈ J(s)`, the origin-firm labor supply elasticity implied by the sector-nested logit is:

```text
ε_zj
= ∂ ln π_zj / ∂ ln w_j
= θ [ 1/σ + (1 - 1/σ) π_zj|s - π_zj ].
```

- When `σ → 1`, this collapses to the one-level-logit elasticity:

```text
ε_zj → θ (1 - π_zj).
```

- Firm `j`'s workforce share from origin `z` is:

```text
γ_zj = π_zj l_z / sum_{z'} π_z'j l_z'.
```

- Firm `j`'s aggregate labor supply elasticity is:

```text
ε_j = sum_z γ_zj ε_zj.
```

- The model markdown is:

```text
ν_j = 1 + 1 / ε_j > 1.
```

Inferred:
- A bridge changes `π_zs`, `π_zj|s`, `π_zj`, `γ_zj`, and `ε_j`, so it can change firm markdowns even without changing productivity or amenities.
- If a firm draws many workers whose outside options are weak or highly correlated within its sector, its effective labor supply elasticity can be lower and its markdown higher.
- Holding `a_j` fixed in counterfactuals, amenities enter markdowns through worker choice probabilities and worker-origin composition.

## Firm Problem

Verified:
- Draft notes define production as:

```text
Y_j(L_j) = A_j L_j^α
```

- The wage first-order condition is:

```text
w_j = α A_j L_j^(α - 1) ε_j / (1 + ε_j).
```

- Product prices are normalized in the paper text unless the user restores price heterogeneity.

Inferred:
- This is a monopsonistic labor-market condition: the wage is marginal product times an elasticity adjustment.
- Lower `ε_j` means a larger wedge between marginal product and wage.
- Firm amenities do not enter the wage first-order condition directly; they affect equilibrium wages through `π_zj`, `L_j`, and `ε_j`.

Uncertain:
- It is unclear whether `α` is fixed globally or should vary by firm/sector in the final paper.

## Bridge Counterfactual

Verified:
- Market-access Stata code labels no-bridge travel time as `dzj` and with-bridge travel time as `dzj_prime`.
- It computes `dma = dzj - dzj_prime`.
- Positive `dma` means travel time is lower with the bridge/tunnel.
- `main.jl` currently solves a baseline using `d` and a counterfactual using `d′`.

Inferred:
- Baseline is no-bridge commuting times; counterfactual is with-bridge commuting times.
- The target nested model should compare the same two commuting networks, but worker reallocation would operate through sector choice and within-sector firm choice.
- The model predicts changes in employment, wages, and markdowns from the commuting-time change alone, holding inferred productivity and firm amenities fixed.

Uncertain:
- The label "baseline" may need clarification because observed wages are from 2010/2011 and the bridge opened in 2011.

## Calibration / Moment Matching

Verified:
- Current `main.jl` sets:

```text
α = 0.8
moment_order = ["labor_bigMA", "labor_bigMA_wdiff", "wage_bigMA"]
σ_bounds = [0.1, 1.0]
```

- Current code reads the three calibration targets from `output/tables/calibration_moments.csv`.
- Current code optimizes over `[η, θ, σ]` using `Optim.NelderMead()` with fixed `α`.
- Current code inverts firm amenities and productivity from baseline observed employment and wages:

```text
Find a_j such that L_j_data = sum_z π_zj(w_data, a, d) l_z
Normalize mean_j(a_j) = 0
A_j = w_j_data (1 + ε_j) / [α L_j_data^(α - 1) ε_j]
```

- Current `ComputeModelMoments` solves the model before and after the travel-time change, appends baseline and counterfactual firm outcomes as two periods, constructs `post`, `bigMA`, and baseline wage heterogeneity `w_diff`, and estimates the two-period DID moment analogue of:

```text
lnL_jt ~ bigMA_j × post_t + w_diff_j × bigMA_j × post_t
          + id_j fixed effects + year_t fixed effects + w_diff_j × year_t slopes
```

- Current model moments are:

```text
β1 = coefficient on bigMA × post
β2 = coefficient on w_diff × bigMA × post
```

- The calibration targets come from empirical results written to `output/tables/calibration_moments.csv`. Local empirical code may still be updated later if the preferred target estimates change.
- The final `bigMA` threshold should be fixed at `0.5` minutes.

Inferred:
- The calibration tries to make model-generated labor reallocation match reduced-form employment heterogeneity by market-access treatment and initial wage.
- The third wage moment now helps discipline `σ` jointly with `η` and `θ`; the model remains exactly identified in moment count but may still need sensitivity checks.

Uncertain:
- The preferred economic interpretation and reporting normalization for inverted firm amenities need confirmation.
- Whether three moments are sufficiently informative for stable separate identification of `η`, `θ`, and `σ` remains an empirical diagnostics question.

## Wage-Ladder Poaching Extension

Verified:
- The candidate extension is implemented in `code/03_model/diagnose_poaching_groups.jl`.
- A detailed solution note is recorded in `docs/WAGE_LADDER_POACHING_GE.md`.
- The extension constructs baseline employer groups as:

```text
g(j) = ind_agg(j) × baseline_wage_decile(j)
```

- It infers origin-by-group baseline masses from model-implied baseline origin shares:

```text
M_zg = sum_{j:g(j)=g} L_j_data * gamma_zj,
gamma_zj = pi_zj * l_z / L_j.
```

- Post-bridge utility for a worker in baseline group `g` choosing firm `j` includes:

```text
rho * max(log(w_j) - group_mean_log_wage_g, 0)
```

- The full-GE solver recomputes group-specific nested-logit probabilities, firm labor, firm labor supply elasticities, and wages until the log-wage fixed point converges.
- The fast amenity inversion in this diagnostic matched baseline firm employment with final gap about `9.89e-6` in 779 iterations.
- In the verified diagnostic run, full GE converged for `rho ∈ {0, 0.25, 0.5, 1}` with `kappa = 0`.

Inferred:
- A positive `rho` gives workers initially attached to low-wage employer groups an additional pull toward higher-wage firms.
- This mechanism can make low-wage control firms lose more workers after the bridge, which turns the control-firm wage slope in employment changes positive.
- The full-GE diagnostic shows that `rho = 0` leaves `beta_labor_wdiff_post` slightly negative, while `rho > 0` makes it positive.

Uncertain:
- `rho` has not yet been formally calibrated.
- It remains open whether this extension should be the final model, a robustness mechanism, or a diagnostic explanation for the sector-only model's sign problem.

## Mechanism As Currently Understood

Verified:
- The user wants the core mechanism to be that traffic integration weakens labor-market power on average.
- The user also wants the model to allow high-wage / large firms with large accessibility improvements to experience increased labor-market power / markdown.

Inferred:
- The bridge reduces commuting costs for some residence-firm pairs.
- Workers reallocate across sectors and across firms within sectors as commuting times change.
- High-wage firms may attract more distant workers after the bridge, changing their worker-origin composition.
- Worker-origin composition and within-sector substitution patterns affect firm labor supply elasticity.
- Changes in firm labor supply elasticity affect labor-market power / markdowns and wages.
- The empirical heterogeneity by initial wage is intended to discipline this mechanism.

Uncertain:
- Current simulations still need to verify whether the nested model output matches the intended average and high-wage-firm markdown patterns.

# Model Notes

Last updated: 2026-05-21.

This document summarizes the structural model as read from current Julia code and draft notes. It distinguishes the current implemented model from older or richer draft formulations.

## Current Implemented Julia Model

Verified:
- Implemented in `code/03_model/functions.jl` and called by `code/03_model/main.jl`.
- Indices:
  - `z`: residence / town / worker origin location.
  - `j`: firm.
- Inputs:
  - `l_z`: residence population, normalized to sum to 1.
  - `d_zj`: commuting time from residence `z` to firm `j`.
  - `w_j`: firm wage.
  - `a_j`: firm-level non-pecuniary amenity, called firm amenity in the final paper, implemented as an additive utility shifter and inverted from observed firm employment.
  - `z_j`: firm productivity, backed out from observed wages.
  - `α`: decreasing returns parameter, set to `0.4` in `main.jl`.
  - `η`: commuting-cost elasticity.
  - `θ`: preference dispersion / responsiveness parameter.

## Worker Choice

Verified:
- Current Julia code corresponds to worker utility:

```text
U_izj = ln(w_j) + a_j - η ln(d_zj) + (1/θ) ε_ij
```

- Current Julia code implies:

```text
π_zj = exp(θ [ln(w_j) + a_j - η ln(d_zj)])
       / sum_k exp(θ [ln(w_k) + a_k - η ln(d_zk)])
```

- Equivalently:

```text
π_zj = (w_j exp(a_j) d_zj^(-η))^θ
       / sum_k (w_k exp(a_k) d_zk^(-η))^θ
```

- Firm labor is:

```text
l_j = sum_z π_zj l_z
```

Inferred:
- Lower commuting time raises the attractiveness of a firm for workers in affected origins.
- Higher firm amenity raises a firm's attractiveness independently of wages and commuting time.
- Higher `θ` makes worker allocation more sensitive to wage/commuting differences.
- Higher `η` makes commuting time more important.

Verified:
- In paper text, `z` should be interpreted strictly as a region, specifically street/town.
- In the final paper, `a_j` should be called firm amenity.

## Labor Supply Elasticity And Markdown

Verified:
- Current code computes:

```text
ε_zj = θ (1 - π_zj)
γ_zj = π_zj l_z / l_j
ε_j = sum_z γ_zj ε_zj
ν_j = 1 + 1 / ε_j
```

- Draft notes interpret `γ_zj` as the share of firm `j`'s workers from origin `z`.
- Draft notes say firm markdown/labor-market power depends on the composition of workers by origin and outside options.

Inferred:
- If a firm draws many workers whose outside options are weak or hard to reach, its effective labor supply elasticity is lower and its labor-market power is higher.
- A bridge changes `π_zj`, `γ_zj`, and `ε_j`, so it can change firm markdowns even without changing productivity or amenities.
- Holding `a_j` fixed in counterfactuals, the wage elasticity formula is unchanged; amenities enter through worker choice probabilities and worker-origin composition.

Verified:
- In the paper, `ν_j = 1 + 1/ε_j` should be called markdown, with `ν_j > 1`.

## Firm Problem

Verified:
- Draft notes define production as:

```text
Y_j(l_j) = z_j l_j^α
```

- Current Julia wage update is:

```text
w_j = α z_j l_j^(α - 1) ε_j / (1 + ε_j)
```

- Wages are normalized in the solver.

Inferred:
- This is a monopsonistic labor-market condition: the wage is marginal product times an elasticity adjustment.
- Lower `ε_j` means a larger wedge between marginal product and wage.
- Firm amenities do not enter the wage first-order condition directly in the current code; they affect equilibrium wages through `π_zj`, `l_j`, and `ε_j`.

Uncertain:
- Product prices `P_j` appear in draft notes but are normalized away in current code.
- It is unclear whether `α` is fixed globally or should vary by firm/sector in the final paper.

## Bridge Counterfactual

Verified:
- Market-access Stata code labels no-bridge travel time as `dzj` and with-bridge travel time as `dzj_prime`.
- It computes `dma = dzj - dzj_prime`.
- Positive `dma` means travel time is lower with the bridge/tunnel.
- `main.jl` solves a baseline using `d` and a counterfactual using `d′`.

Inferred:
- Baseline is no-bridge commuting times; counterfactual is with-bridge commuting times.
- The model predicts changes in employment, wages, and markdowns from the commuting-time change alone, holding inferred productivity fixed.

Uncertain:
- The label "baseline" may need clarification because observed wages are from 2010/2011 and the bridge opened in 2011.

## Calibration / Moment Matching

Verified:
- `main.jl` sets:

```text
α = 0.4
β_target = [-0.06990708, -0.02787439, 0.02603689, 0.04537406]
x0 = [8.0, 0.75]
```

- It optimizes over `[η, θ]` using `Optim.NelderMead()`.
- For each candidate `[η, θ]`, `functions.jl` inverts firm amenities and productivity from baseline observed employment and wages:

```text
Find a_j such that l_j_data = sum_z π_zj(w_data, a, d) l_z
Normalize mean_j(a_j) = 0
z_j = w_j_data (1 + ε_j) / [α l_j_data^(α - 1) ε_j]
```

- `main.jl` uses `a_j = 0` only as the starting guess for the amenity inversion.
- `ComputeModelMoments` solves the model before and after the travel-time change, appends baseline and counterfactual firm outcomes as two periods, constructs `post`, `BIG`, demeaned baseline log wage (`demean_lnw0`), and demeaned log market-access intensity (`lndma`), and estimates the two-period DID moment analogue of `code/02_empirical/calculate_calibration_4_moments.do`:

```text
lnl_jt ~ BIG_j × post_t
         + lndma_j × BIG_j × post_t
         + demean_lnw0_j × BIG_j × post_t
         + demean_lnw0_j × lndma_j × post_t
         + lndma_j × post_t
         + id_j fixed effects + year_t fixed effects
         + demean_lnw0_j × year_t slopes
```

- As in Stata, `lndma × post` is collinear with the treatment and treated-intensity terms and is omitted; the model moments are the four reported coefficients:

```text
β1 = coefficient on BIG × post
β2 = coefficient on lndma × BIG × post
β3 = coefficient on demean_lnw0 × BIG × post
β4 = coefficient on demean_lnw0 × lndma × post
```

- `β_target` should contain the four reported coefficients from `calculate_calibration_4_moments.do`.
- The final `bigMA` threshold should be fixed at `0.5` minutes.
- For now, no external labor-supply-elasticity moment or other moment should be added to discipline `θ`.

Inferred:
- The calibration tries to make model-generated labor reallocation match reduced-form employment heterogeneity by market-access treatment and initial wage.

Uncertain:
- The preferred economic interpretation and reporting normalization for inverted firm amenities need confirmation.

## Draft Model Variants

Verified:
- `ModelBridge.lyx` and `Model.lyx` include a richer nested preference structure with a parameter `σ` for within-town/within-street firm preference correlation.
- These drafts include:

```text
π_zj = π_zs × π_zj|s
ε_zj = θ [1/σ + (1 - 1/σ) π_zj|s - π_zj]
```

- Current `functions.jl` does not implement `σ`.

Inferred:
- The current Julia code corresponds to the special or simplified non-nested case.

Verified:
- The final model should continue using the current one-level logit structure.

## Mechanism As Currently Understood

Inferred:
- The bridge reduces commuting costs for some firm-origin pairs.
- Workers reallocate toward firms made more attractive by shorter commutes and/or higher wages.
- Workers also reallocate toward higher-amenity firms if nonzero amenities are supplied.
- High-wage firms may attract more distant workers after the bridge, changing their worker-origin composition.
- Worker-origin composition affects firm labor supply elasticity.
- Changes in firm labor supply elasticity affect labor-market power / markdowns and wages.
- The empirical heterogeneity by initial wage is intended to discipline this mechanism.

Uncertain:
- The sign of the markdown effect for high-wage versus low-wage firms is not fully settled across notes and current simulations.

## Implementation Warnings

Verified:
- Current `data/model/processed/firm_qingdao_model.dta` lacks variables needed by `main.jl`.
- `model_output.jl` appears incompatible with the current model functions.
- `calibration.jl` is a placeholder with undefined empirical moments.
- `functions.jl` now solves firm amenities from `l_j_data` and firm productivity from `w_j_data` for each calibration parameter guess.
- Direct calls to `SolveModel` and `SolveZfromW` now require `a_j` in `vars`; there is no zero-amenity fallback inside the solver.

Inferred:
- Run `main.jl` only after confirming generated market-access files and model input variables are current.

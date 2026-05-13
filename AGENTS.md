# AI Project Guide

Last generated: 2026-05-07.

This repository is a quantitative spatial / labor-market-power project about a Qingdao cross-sea bridge transport shock. The current LaTeX paper body is only a scaffold, so AI agents should rely first on code, data structure, and draft notes rather than inventing narrative claims.

## Mandatory Reading Order

1. `docs/AI_PROJECT_UNDERSTANDING.md`
2. `docs/CODE_MAP.md`
3. `docs/DATA_MAP.md`
4. `docs/MODEL_NOTES.md`
5. `docs/OPEN_QUESTIONS.md`
6. `docs/USER_CORRECTIONS.md`

## Grounding Rules

- Separate all understanding into `Verified`, `Inferred`, and `Uncertain`.
- Do not invent the research question, identifying assumptions, data sources, calibration targets, or results.
- Treat `manuscript/notes` and `manuscript/archive` as historical working notes unless the user says they are authoritative.
- Treat current generated data and output files as possibly stale. Several scripts expect files or variables not present in the current working tree.
- Before editing code, inspect the relevant script and check whether a similarly named no-security, sandbox, or legacy file exists.

## Current High-Level Understanding

Verified:
- The paper title in `manuscript/paper/main.tex` is `Commuting Cost, Labor Market Power and Misallocation: Evidence from Chinese Bridge`.
- Main code comments describe the project as `Cross-Sea Bridge, Labor Market Power and Welfare` and `Effect of Cross-sea Bridge on Firms`.
- The data pipeline uses Qingdao firm data, town-level geography, OSM road networks, bridge/tunnel road segments, CTSD, and CIED.
- The current structural model uses workers in residence locations `z`, firms `j`, commuting time `d_zj`, wages `w_j`, worker choice probabilities, firm labor supply elasticities, and a monopsony markdown term.

Inferred:
- The likely empirical setting is the 2011 opening of the Jiaozhou Bay Bridge / tunnel system as a shock to commuting costs and market access for Qingdao manufacturing firms.
- The project seems to connect reduced-form employment, wage, and markdown effects with a quantitative spatial model calibrated to empirical moments.

Uncertain:
- The exact final research question, preferred sample, target empirical moments, and final identification strategy still need user confirmation.

## Edit Discipline

- Preserve user work. Do not modify data, raw manuscripts, or code unless the user explicitly asks.
- This documentation pass created only the approved documentation files.
- If new user corrections arrive, update `docs/USER_CORRECTIONS.md` first, then revise the relevant understanding docs.

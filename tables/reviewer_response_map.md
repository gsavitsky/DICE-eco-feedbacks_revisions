# Reviewer comments: status against the corrected model (2026-11-20)

Status key:
- **Done**: the analysis exists.
- **Text**: needs writing or rewording only.
- **Open**: needs a decision or new work.

Unless noted, numbers are from the corrected model with GAMS controls and standard DICE settings.

## Reviewer 1 (carbon-cycle / permafrost expert; minor revision)

| # | Comment | Status | Answer / evidence |
|---|---|---|---|
| 1 | Add a gradual-thaw-only (Hector) version | Done | Abrupt off vs on: abrupt thaw adds ≤0.01 °C by 2100 and 0.09–0.12 °C by 2300. Report the gradual-only configuration alongside. |
| 2 | Permafrost feedback too high (440 PgC in atm); compare with Schuur 2022 (149 PgC), AR6 (≤240 PgC by 2100) | Done + Text | Corrected: about 105 PgC released by 2100 and 38 PgC/°C (AR6 18, range 3.1–41; Schädel 63). This is inside the assessed ranges. Explain: the old 1035 PgC pool and step abrupt term (see R3.2). |
| 3 | Graph of DICE CO2 emissions | Done | `fig_emissions.png` |
| s | Line numbers; S-curve wording; "negative feedbacks" terminology; explain that negative emissions cause the cooling; CO2 notation; give ppm; second axis on Fig 2b | Text / Done | ppm axis is already in `Visualizations_corrected.R`. Use "positive/negative (reinforcing/balancing)" on first use. |

## Reviewer 2 (design / presentation)

| # | Comment | Status | Answer / evidence |
|---|---|---|---|
| B | Shorten; engage normative and justice aspects (Rivadeneira & Carton 2022; Keen 2021) | Text | The new framing does this: late, partly locked-in warming falls on later generations, and the global welfare metric hides local Arctic impacts. |
| 1 | Reframe around damage-function insensitivity | Open (respond, do not adopt) | Corrected results and R3.5 undermine the "structural" claim: with a quadratic damage function, a larger coefficient only rescales damages. We reframe around timing (late, locked-in warming) instead. |
| 2 | Factorial: show the distribution of the off→on change; is each factor clearly balancing or reinforcing? Too many significant figures | Done | New factorial, `factorial_off_to_on_distribution.csv`. Signs are consistent in 32/32 contexts at the peak for every factor. Fertilisation flips sign by 2200 in 8/32 contexts (it releases carbon when CO2 falls during the removal phase). |
| 3 | No rationale for 400 years; stop at 2200 | Done | Simulations end in 2300 and the main reported year is 2200. The factorial ends in 2200 because after 2200 the DICE control path's sustained net removal pushes every Lenton run below 280 ppm. |
| 4 | Fig 7 is indistinguishable; plot a shorter window | Text / figure | — |
| 5 | Jumps in abatement cost at 2200 and 2300 | Text | Caused by the DICE 2023 μ upper-bound schedule (1.1 → 1.05 around 2205, → 1.0 around 2305), not by the feedbacks. |
| 6–7 | Howard–Sterner is only a parameter change, not structural | Done (removed) | The damage and discounting cells are dropped from the main paper. |
| s | §1.2 fertilisation should be "balancing" | Text | — |

## Reviewer 3 (recommends reject and resubmit)

| # | Comment | Status | Answer / evidence |
|---|---|---|---|
| 1(i) | Why Lenton rather than FRIDA-Clim or Lade et al.? | Text + Done | Lenton gives process-level attribution that DFAIR cannot. DFAIR has only two aggregate terms (irC, irT). Lenton separates fertilisation, photosynthesis, plant and soil respiration, and ocean buffer. See `factorial_effects.csv` and `lenton_grouped_effects.csv`. Still need to discuss FRIDA-Clim and Lade. |
| 1(ii) | Add permafrost on top of DFAIR | Done | "DICE 2023 + permafrost" in all runs. Permafrost effect: DFAIR 38 vs Lenton 39 PgC/°C. ZEC: +0.10 vs +0.11 °C. |
| 1(iii) | Contribution beyond Kessler 2017 and Wirths 2018 | Text | Timing and locked-in thaw (about half of warming to 2100 comes from carbon thawed before 2020); a literature-bounded SCC range; robustness to overshoot; process attribution. |
| 1 | FaIR has no major gap except permafrost | Concede | DICE-C ≈ DICE 2023 (peaks 2.56 vs 2.63 °C), which agrees with the reviewer. |
| 2a | Abrupt thaw misattributed; step at 2 °C; compare a smoother version | Done | Now smooth Armstrong McKay: 1.5 °C, 200 yr, +50%. Threshold test at 2.3 °C also run. |
| 2b | 1035 PgC pool includes non-frozen soil; Woodard uses 865 | Done | 865 PgC pre-industrial, spun up to 2020 (691.5 frozen + 150.4 thawed). |
| 2c | Magnitude about 80 GtC/°C vs AR6 18 and Schädel 63; report gradual-only | Done | Now 38 PgC/°C; gradual-only reported (see R1.1). |
| 3a | "Policy without full information" wording | Text | Use "policy that disregards the feedbacks". |
| 3b | Do not report beyond 2300 | Done | Simulations end in 2300. |
| 3c | Feedbacks do not "undermine" mitigation; they make it more valuable | Done + Text | Agreed. The SCC rises 2–28% with permafrost, about 12% at Hector strength. Remove the "undermine" claim. |
| 4a | DICE-C cooler than DICE 2023 for 50 years is implausible | Done | That was the old Lenton bug. Corrected DICE-C tracks DICE 2023 to 2080 (`ba1`, `ba2`). |
| 4b | Lenton worse at emulating ESMs | Done + Text | C4MIP at 4×CO2 vs CMIP6 (Arora 2020 Table A1, mean ± s.d.). Land β: Lenton 0.69 vs 0.97 ± 0.40. Ocean β: 0.86 vs 0.79 ± 0.07. Land γ: −32 vs −45 ± 51 PgC/K. All within the CMIP6 spread. Ocean γ: 0 (variant iii, no ocean warming) vs −17 ± 5; ocean warming (Lenton's variant iv) is now an optional process. DFAIR total β 1.84 / γ −97 vs CMIP6 total 1.76 / −62 ± 51, also within. |
| 4c | ScenarioMIP comparison invalid (concentration-driven) | Done | Replace with the emissions-driven RCMIP SSP check (`ssp_check.py`) and the ZEC test against MacDougall. |
| 5 | Welfare percentages depend on utility normalisation; use consumption-equivalent or SCC; Howard–Sterner description | Done / removed | Report the SCC (`tab_scc_full`). Damage-function cells dropped. |
| 6 | Tipping-point definitions; intro promises (policy adherence, social dynamics) the paper doesn't keep; the paper is about permafrost, not ecological feedbacks generally | Text | Rewrite §1.2 and §1.4; the title should name permafrost. |
| m | Minor points: fertilisation wording; missing citations; Dietz +8.4% for permafrost alone; Wirths 2018; explain STELLA → Python, fixed controls; "pre-indistrial"; STP terminology; "continuous-time"; "three reasons"; SCM abbreviation | Text | Rename the Lenton configuration away from "SCM" (for example "Lenton carbon cycle"). |

## DFAIR feedback strengths (new, 2026-11-20)

Switching each term off, relative to default, permafrost off:

| Term switched off | Peak | 2100 | 2200 |
|---|---|---|---|
| irC (uptake saturation) | −0.33 °C | −0.31 °C | −0.29 °C |
| irT (temperature) | −0.21 °C | −0.19 °C | −0.17 °C |
| Both | −0.60 °C | −0.57 °C | −0.48 °C |

Effect of background feedback strength on permafrost warming, with irC and irT each scaled from 0× to 2×:

| | Peak | 2200 |
|---|---|---|
| Range | 0.09–0.24 °C | 0.19–0.57 °C |
| Default | 0.18 °C | 0.40 °C |

Lenton grouped (full model):

| Group | Peak | 2200 |
|---|---|---|
| Concentration (ocean + fertilisation) | −0.10 °C | +0.36 °C |
| Temperature (photosynthesis + plant + soil respiration) | +0.10 °C | +0.08 °C |

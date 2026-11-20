# Chapter 2 (DICE-C / DICE-CP): updated code, figures and notes (17 Nov 2026; updated 20 Nov 2026)

## 0. Update, 20 Nov 2026 (read first; supersedes sections 1-4 where they differ)

### Model and set-up changes
| Change | Detail |
|---|---|
| Lenton variant | **(iv)**, Lenton's ocean-warming version (surface boxes warm with global T, his Eq. 16) with his K_M = 145. Was (iii), his control (fixed ocean T, K_M = 120). `LENTON_VARIANT=iii` restores it. 2020 state: `lenton_km145_iv_state2020.json` (418.6 ppm). |
| Horizon | **2020-2200** everywhere (simulations, welfare, SCC). After 2200 the DICE 2023 control path's sustained net removal (μ > 1 from ~2120) pushes every Lenton run below 280 ppm, and some feedback-switch combinations reach zero CO2 (`check_2300_physical.py`). |
| Feedback switches | `cfg['fbx']` in `ch2_models.py`: ocean buffer, ocean warming, CO2 fertilisation, photosynthesis T response, plant and soil respiration T responses. "Off" = the process keeps its pre-industrial response from 1800 (own spin-up per combination). `cfg['KM']` overrides fertilisation strength. |
| Other model options | prescribed emissions `E_fn` (GtC/yr or (fossil, land-use)), non-CO2 forcing `F_nonco2` (constant or function of year), reversible abrupt thaw, running-mean Q10 (sensitivity tests only). Defaults reproduce earlier runs exactly. |
| Chapter 3 optimiser | `code/dice_pyomo_corrected.py` uses variant (iv): 2-D chemistry fit (DIC x ocean warming, error 1.6e-6) and temperature-dependent solubility. Forward model matches `ch2_models.py` (peak T within 7e-4 C). VM runs pending: `chapter3/VM_VARIANT_IV_RUNS.md`. |

### Evaluation (variant iv)
- Lenton's own experiment (iv): 1018 ppm / 5.72 K (his 1020 / 5.7). C4MIP (1%/yr to 4xCO2): land β 0.81, ocean β 0.86, land γ −33, ocean γ −21 PgC/K, all within CMIP6 ±1 s.d. (Arora et al. 2020, Table A1: 0.97±0.40, 0.79±0.07, −45±51, −17±5). DFAIR total β 1.84, γ −97 (CMIP6 total 1.75, −62).
- Emissions-driven SSPs (`ssp_check_v2.py`, RCMIP v5.1.0): CO2 rise 2020-2100 in SSP2-4.5 is 184 (Lenton + permafrost), 193 (DFAIR + permafrost) vs 188 (MAGICC); permafrost releases 75-128 PgC by 2100 across SSPs.
- Hector validation: frozen pool matches Hector's own update rule within ~1% (RCP2.6-8.5, 2000-2100). Zero-emissions commitment from permafrost: +0.10 C after 50 yr on both carbon cycles (MacDougall: 0.06, 0.02-0.14).
- Figure: `R/vis_corrected/figures/fig_model_validation.png`.

### Main results (GAMS controls, standard settings, 2020-2200)
| | SCC 2020 | Peak T | T 2100 / 2200 |
|---|---|---|---|
| DICE 2023 | $45.0 | 2.63 | 2.59 / 1.81 |
| DICE 2023 + permafrost | $48.9 | 2.81 | 2.73 / 2.21 |
| DICE-C | $45.6 | 2.59 | 2.56 / 1.53 |
| DICE-CP | $49.9 | 2.77 | 2.71 / 1.92 |

- Permafrost (Hector strength): +0.15 C in 2100, +0.39 C in 2200; about half of the 2100 warming is from carbon thawed before 2020. SCC +9.4% (Lenton), +8.8% (DFAIR); gradual-only +8.3% (`ch2_gradual_only.csv`).
- **SCC by carbon-cycle process (Shapley values, 128-run factorial, Lenton), 2020 pulse:** all feedbacks at pre-industrial $42.8; ocean buffer +13.8, ocean warming +2.2, CO2 fertilisation −14.2, photosynthesis T −12.0, plant respiration +6.2, soil respiration +7.8, permafrost +3.3; all on $49.9. 2080 pulse: $189 → $231, same pattern (`shapley_values_p2080.csv`).
- Interactions: only the ocean-buffer pairs are sign-consistent across all 32 backgrounds (× fertilisation −7.9, × photosynthesis −6.9, × soil +4.5, × plant +3.6, × permafrost +1.5).
- DFAIR (2 aggregate terms + permafrost): $24.3 → uptake saturation +13.3, temperature term +8.6, permafrost +2.7 → $48.9.
- Across CMIP6 land β (K_M 100 / 145 / 275, each spun up; 2020 CO2 then differs from observed): fertilisation −10.5 / −14.2 / −21.7, permafrost 3.4 / 3.3 / 3.0, **net feedback +12.6 / +7.1 / −3.5** (sign depends on fertilisation strength).
- DFAIR scaled to the CMIP6 β/γ ±1 s.d. corners: SCC $34-59, permafrost premium 7.7-9.2% (`dfair_cmip6_range.csv`).
- Zero time preference at 2200: DICE 2023 SCC $472 (was $1,630 to 2300); permafrost premium +19-21%.

### New scripts (code/permafrost_revision_test/) and outputs (R/vis_corrected/)
| Script | Output |
|---|---|
| `ch2_final_runs_parallel.py`, `ch2_tables.py` | main grid, sweep, tables (now 2200) |
| `ch2_warming_time.py` | `ch2_permafrost_warming_time(series).csv` |
| `ch2_gradual_only.py` | `ch2_gradual_only.csv` |
| `ch2_feedback_factorial.py` | `feedback_factorial/factorial_*`, `fert_sweep_*` (6-process factorial; fertilisation-strength sweep) |
| `ch2_scc_shapley.py` (env `SHAPLEY_KM`, `SHAPLEY_TAG`, `SHAPLEY_PULSES`, `SHAPLEY_SKIP_DFAIR`) | `feedback_factorial/shapley_runs*.csv`, `shapley_values*.csv`, `dfair_cmip6_range.csv` |
| `ch2_shapley_interactions.py` | `shapley_interactions.csv` (incl. equal-weight, all-on, range, sign consistency), `shapley_off_to_on.csv` |
| `ch2_dfair_shapley.py` | `dfair_shapley_values.csv`, `dfair_shapley_interactions.csv` |
| `ch2_dfair_feedbacks.py` | DFAIR irC/irT scaling, `c4mip_compare.csv`, `lenton_grouped_effects.csv` |
| `ssp_check_v2.py` | `ssp_check_v2.csv` (RCMIP inputs in `hist_data/` are World-level SSP extracts) |
| `overshoot_checks.py` (`ZEC_ONLY=1` for the zero-emissions test) | zero-emissions and overshoot checks (overshoot results not used in Paper 2) |
| `export_permafrost_before_after.py` | `permafrost_before_after.csv` (original pools from a Python replication of the STELLA module, within 0.008 C of STELLA) |
| `km_beta_map.py`, `lenton_iv_evaluation.py`, `check_2300_physical.py` | K_M → land β map; C4MIP β/γ for (iii)/(iv); 2300 feasibility check |
| `run_all_2200.sh`, `run_paper2_extras.sh` | reruns everything in order |

R figure scripts (`R/`): `Visualizations_before_after.R` (`figures/before_after/ba1-ba4`), `Visualizations_scc_waterfall.R` (args: values file, pulse years, output), `Visualizations_scc_interactions.R`, `Visualizations_scc_decomposition.R` (Lenton vs DFAIR), `Visualizations_permafrost_timing.R`, `Visualizations_validation.R`. Python scripts write to `~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected`; edit the `OUT` paths to write elsewhere.

### Still to do
- `tables/methods_revised.tex` still describes the 2300 horizon and does not yet describe the factorial / Shapley method.
- Chapter 3 re-solve with variant (iv) on the VM; `chapter3/` results below are still variant (iii) to 2420.
- `tables/reviewer_response_map.md` predates the variant (iv) switch and the Shapley analysis.

---


## 1. What is in this folder

| Folder | Contents |
|---|---|
| `code/` | Python model and data. `permafrost_revision_test/ch2_models.py` is the model; `ch2_final_runs_parallel.py` runs every Chapter 2 experiment; `ch2_tables.py` writes the LaTeX tables. `dice_pyomo_corrected.py` is the Chapter 3 re-optimisation with the same physics. |
| `R/` | `Visualizations_corrected.R` (the original figure script pointed at the corrected outputs) and `Visualizations_permafrost_hooks.R` (SCC-bounds, welfare-range and emissions figures). `R/vis_corrected/` holds the current outputs (CSV) and `figures/`. |
| `tables/` | LaTeX tables generated from the outputs, the revised Methods section (`methods_revised.tex`), the parameter table (`tab_lenton_params.tex`) and new bibliography entries (`new_references.bib`). |
| `chapter3/` | Chapter 3 re-optimisation results (corrected physics): optimal controls and trajectories (`results/`), outcomes table (`ch3_corrected_outcomes.md`/`.csv`), figure data (`figdata_corrected/`), figures (`figures_ch3_corrected/`, with `caption_notes.md`), the scripts that make them (`ch3_outcomes_corrected.py`, `ch3_figures_corrected.R`), and the VM instructions (`VM_CORRECTED_RUNS.md`). The optimisation model itself is `code/dice_pyomo_corrected.py`. |

### Model configurations
All configurations share the DICE 2023 economy, two-box temperature model, forcing, damages and welfare. They are driven by the DICE 2023 cost–benefit optimal controls (GAMS). Simulations run 2020–2300.

| Configuration | Carbon cycle | Permafrost |
|---|---|---|
| DICE 2023 | DFAIR | none |
| DFAIR + permafrost | DFAIR | Hector v2.3pf + smooth abrupt thaw |
| DICE-C | Lenton (2000) | none |
| DICE-CP | Lenton (2000) | Hector v2.3pf + smooth abrupt thaw |

Each is run under standard / Howard–Sterner damage × standard / zero-time-preference discounting (16 runs). There is also a permafrost-strength sweep (respiration-rate multiplier m = 0.1–50 on both carbon cycles, all four damage × discounting cells, with the SCC in every cell) and a permafrost-pool "ceiling" test.

### How to run
Python 3.11 with `numpy`, `scipy`, `pandas`, `openpyxl` and `matplotlib`. The Chapter 3 script additionally needs `pyomo`, `ipopt` and HSL MA27.

```bash
cd code/permafrost_revision_test
python3 ch2_final_runs_parallel.py     # ~233 runs, ~5 min on 8 cores -> R/vis_corrected/
python3 ch2_tables.py                  # -> tables/*.tex
cd ../../R
Rscript Visualizations_corrected.R     # -> R/vis_corrected/figures/
Rscript Visualizations_permafrost_hooks.R
```

- Importing `dice_stella_validation.py` runs its own checks against the old STELLA export (`FINALresults.xlsx`) and takes about 1 minute. The module is used only for parameters, exogenous drivers and the GAMS controls.
- The feedback factorial is switched off by default (`RUN_FACTORIAL=1` enables it; see the to-do list).

A single run from Python:
```python
import ch2_models as M, permafrost_revision_sim as prs
d = M.run(dict(carbon='lenton', perma=dict(m=1.0, abrupt='smooth'),
               a2=prs.DAMAGE['std'], rartp=prs.DISC['std']), nyears=280)
# DICE-CP: peak warming 2.736 C, permafrost release 2020-2100 = 105.0 PgC, U_acc(2300) = 12045.464
```

---

## 2. Mistakes found in the submitted model (STELLA `DICE 2023 5_12.stmx` and its Python/Pyomo copies)

### Lenton carbon cycle
1. **Ocean interior disconnected (units).** The five ocean-mixing flows were multiplied by `3.15576e7*12e-15`, a mol→PgC conversion, although the pools are already in PgC. Lenton writes his equations in moles. Interior mixing was about 1e-13 PgC/yr, so intermediate and deep carbon never changed and ocean uptake after 2025 was about 0.4 PgC/yr (observed ~2.9).
2. **Incomplete circulation.** The model lacked warm→cold surface and cold→deep advection, and had cold→intermediate where Lenton's Fig. 1 has intermediate→cold. There was no biological pump and no nitrate or alkalinity; a Revelle power law stood in for the carbonate chemistry.
3. **Deep-ocean pool** set to 23,860 PgC; Lenton's Fig. 1 has 26,830.
4. **Initial-state shock.** The atmosphere started at its 2020 value (887 PgC) while land and ocean started at pre-industrial values. About 50 PgC left the atmosphere in 2020–21 and CO2 fell 28 PgC in two years. This caused the early DICE-C cooling that Reviewer 3 flagged.
5. **Land use.** Deforestation divided GtCO2 by 2.12 instead of 3.667, and used a separate land-emission series (2.6) inconsistent with the economy's (5.9). Land-use carbon was also not removed from vegetation (Lenton Eq. 17).
6. **Terrestrial CO2 switch (feedback factorial).** With the switch off, CO2 was set to 0.0004165 in a ppmv equation, so photosynthesis went to zero and land released ~2,000 PgC. This produced the "dominant" fertilisation feedback in the old Fig. 4.

### Permafrost module (after Woodard et al. 2021, Hector v2.3pf)
7. **Missing 0.02/yr respiration rate.** Hector respires thawed carbon at (thawed − static) × 0.02 × Q10^(T200/10). The STELLA model had no 0.02, so all thawed carbon reached the atmosphere within a few years: 304 PgC by 2100, about 95 PgC/°C against AR6's 18. This was checked in the v2.3pf source (Zenodo 10.5281/zenodo.4876800, `simpleNbox.cpp`). Hector's static carbon is bookkeeping equal to f_static × thawed, so the static-fraction treatment itself was not the error; Eq. 12 as printed in the paper is misleading on this point.
8. **Initial stock.** 1035 PgC includes non-frozen soil. Hector uses 865 PgC pre-industrial, which gives about 730 PgC frozen in 2010. The model also started with no thawed carbon in 2020.
9. **Frozen fraction.** Used ln(1 + T_HL) and σ = 0.99. Hector uses ln(T_HL) and σ = 0.986.
10. **T200** was a 200-year *delay* of global temperature, with a synthetic ramp before 2220. Hector uses a 200-year *running mean* of T_HL.
11. **Abrupt thaw.** A step to +40% of gradual thaw at 2 °C. Armstrong McKay et al. (2022) give a threshold of ~1.5 °C (1.0–2.3), a timescale of ~200 years and about +50% of gradual-thaw emissions.

---

## 3. Changes and rationale

| Change | Rationale |
|---|---|
| Rebuilt Lenton's four-box ocean: Fig. 1 circulation, Knox & McElroy (1984) biology (C:N 7, CaCO3:Corg 0.25, 30% remineralised in intermediate water, α = 0.01), and full carbonate chemistry (Weiss; Goyet & Poisson; Dickson). Mean alkalinity is tuned to 2326 µmol/kg to give Lenton's 280 ppm steady state. | Reproduce Lenton (2000) as published. |
| Lenton's own calibration: K_M = 120, k_D = 0.27, variant (iii) (his control, no ocean-surface warming); land use removed from vegetation. | Faithful to Lenton. A re-tune to K_M = 209 under Global Carbon Budget emissions was tried and rejected as a departure from Lenton. |
| Carbon cycle spun up from the 1800 steady state with Lenton's Table 3 emissions to 1990, then GCB emissions to 2019 (2020 CO2 = 419 ppm). | Removes the initial-state shock; land and ocean start consistent with history. |
| Permafrost made Hector-faithful: 0.02/yr rate, ln(w_f T) with σ = 0.986, 200-year running-mean T200 (sticky Q10), refreeze limited to the thawed pool. | Match the source model. |
| Permafrost pools spun up 1800–2020 from 865 PgC, giving 692 PgC frozen and 150 PgC thawed in 2020. | The 2020 thawed pool is consistent with Hector's history and with the modern stock (~730 PgC in 2010). |
| Smooth Armstrong McKay abrupt thaw: +50% of cumulative gradual thaw, onset w(T) = ½[1 + tanh((T − 1.5)/0.5)], 200-yr timescale, irreversible. | Removes the global discontinuity at 2 °C. Hector deliberately omits abrupt thaw and Woodard et al. say it may underestimate the feedback, so this is not double counting. |
| New configuration: DFAIR + permafrost. | Reviewer 3: separates the permafrost effect from the carbon-cycle swap. |
| Python forward model (RK4, Δt = 0.05 yr) replaces STELLA for Chapter 2 results. | Before the fixes it reproduced the STELLA runs to within 0.01 °C. STELLA is no longer needed; the `.stmx` is superseded. |
| Simulations end in 2300; welfare and SCC are measured at 2300. Figures report to 2300. | DICE 2023's control schedule changes after 2300 (end of net-negative emissions). Addresses Reviewers 2 and 3 on the horizon. |
| SCC by pulse method in every run. | Reviewer 3 asked for the SCC; percentage welfare changes hide policy-relevant effects. |

### Evaluation of the corrected model (`tables/tab_model_evaluation.tex`)
- **Lenton's preindustrial state (Fig. 1):** reproduced to within about 2%.
- **Lenton's extended business-as-usual experiment,** run with his inputs (Lenton's value in brackets):
  - 1990 CO2 352.7 ppm (353);
  - land becomes a carbon source in 2113 (~2120);
  - peak CO2 984 ppm (985);
  - peak warming 5.53 K (5.5).
- **C4MIP feedback metrics:**
  - land β 0.69 (CMIP6 1.22 ± 0.40);
  - ocean β 0.86 (0.79 ± 0.07);
  - land γ −32 PgC/K (−58 ± 38).
- **Permafrost at Hector strength:** 105 PgC released over 2020–2100, about 39 PgC/°C. This matches Hector's own runs and sits near the top of AR6's range (18, 3.1–41).

---

## 4. Main results (GAMS controls, 2020–2300)

| Configuration | Peak T (°C, year) | 2020 SCC, std discounting ($/tCO2) | SCC, HS damage | Welfare, std damage |
|---|---|---|---|---|
| DICE 2023 | 2.63 (2113) | 48.2 | 140.7 | — |
| DFAIR + permafrost | 2.81 (2117) | 53.9 (+11.7%) | 158.0 (+12.3%) | −0.10% |
| DICE-C | 2.56 (2108) | 47.6 (−1.3%) | 138.5 (−1.5%) | +0.05% |
| DICE-CP | 2.74 (2112) | 53.3 (+10.5%) | 155.7 (+10.7%) | −0.04% |

**Carbon-cycle structure** (DICE-C vs DICE 2023) changes little: peak −0.07 °C, SCC −1.3%. Under zero time preference the difference is larger, with an SCC 43% lower.

**Permafrost strength sweep:**
- Across 5–89 PgC/°C, permafrost raises the 2020 SCC by +1.6% to +28%. It is about +12% at Hector strength, bracketing Dietz et al. (2021)'s +8.4%.
- The result is nearly identical for DFAIR and Lenton, and for the two damage functions. Howard–Sterner triples SCC levels but leaves the permafrost percentage unchanged.

**Thresholds:**
- SCC +10% requires about 31–33 PgC/°C.
- SCC +25% requires about 81–83 PgC/°C.

**Discounting sensitivity (zero time preference):**
- The permafrost SCC effect is several times larger: about +62% at Hector strength.
- +10% is reached at about 5 PgC/°C.

**Welfare** changes stay below 1% in all cases.

**Upper end of the sweep:** above m = 8 (~89 PgC/°C), the spin-up has already respired most of the thawed carbon before 2020. That would imply more than 100 PgC of pre-2020 permafrost emissions, which history does not support, so the sweep is cut at m = 8. Within this module the maximum is set by how much carbon thaws (pool size × frozen-fraction parameters), not by the respiration rate (`ch2_permafrost_ceiling.csv`).

## 4b. Chapter 3: re-optimised policy (corrected physics, 17 Nov 2026)

`code/dice_pyomo_corrected.py` uses the same physics as the Chapter 2 model.
- **Ocean:** 4 DIC states. Alkalinity and nitrate stay exactly steady, so biology enters as constant source and sink terms. The carbonate chemistry is fitted to within 4e-6.
- **Permafrost:** Hector-faithful, with the 2020 pools from the spin-up.
- **Controls:** annual μ and s, with the GAMS bounds (including TATM ≥ 0.5 °C).
- **Time step:** one RK4 step per year is numerically unstable for the corrected surface ocean, so solves use 2 and 5 substeps per year.

**Convergence:** every solve ended "optimal". For each model, two starts (GAMS controls and the previous optimum) and two resolutions (2 and 5 substeps) give the same objective:
- DICE-C: 6449.4698 (padded objective);
- DICE-CP: 6437.3011.

**Consistency:** the Chapter 2 forward model driven by the optimal controls reproduces the optimised temperatures to within 0.001 °C. Reporting window is 2020–2300; the optimisation solves to 2470 to avoid end-of-horizon effects.

| Each model at its own optimum | DICE 2023 | DICE-C | DICE-CP |
|---|---|---|---|
| Peak T (year) | 2.76 (2117) | 2.73 (2116) | 2.83 (2115) |
| μ in 2050 / 2100 | 0.359 / 0.795 | 0.356 / 0.763 | 0.375 / 0.811 |
| Year μ first reaches 1 | 2119 | 2125 | 2119 |
| Cumulative industrial CO2, 2020–2100 | 760 GtC | 777 GtC | 738 GtC |
| 2020 SCC (2019$/tCO2) | 48.3 | 48.3 | 53.0 (+9.8%) |
| Welfare vs DICE 2023 | — | +0.05% | −0.05% |
| Welfare gain from re-optimising vs GAMS controls | +0.010% | +0.014% | +0.005% |

What the results show:
- **The policy response to permafrost is modest.** DICE-CP abates a little more mid-century, offsetting about 20% of the ~105 PgC of permafrost release by 2100. It stays ~0.07 °C warmer than DICE 2023 at peak.
- **Re-optimising raises peak warming relative to the GAMS controls in every model** (+0.10 to +0.18 °C). This is the annual- vs five-year-grid effect already seen for DICE 2023, so configurations should be compared at their own optima.
- **Permafrost refreezes as the climate cools.** After ~2130 carbon returns to the frozen pool, and the thawed pool is empty by ~2255.
- **DICE-C reaches the inherited 0.5 °C lower bound on temperature near 2290.**

To rebuild the outputs:
```bash
cd chapter3
python3 ch3_outcomes_corrected.py
Rscript ch3_figures_corrected.R
```

---

## 5. To do

### Before submitting the draft dissertation (~3 weeks)
1. **Results, Discussion, Abstract and Conclusions:** rewrite with the new numbers and framing (centre: bounding the SCC of permafrost carbon; DICE-C vs DICE 2023 as a short subsection; discounting as a sensitivity).
2. **Methods:** insert `methods_revised.tex`, `tab_lenton_params.tex` and the new bibliography entries, then compile and fix cross-references. Remove claims about a "200-year lagged response" and feedback-driven "tipping points".
3. **Tables:** replace `summary.tex` and the supplement tables with the generated ones; add the SCC, thresholds and model-evaluation tables; update the initial-conditions table.
4. **Figures:** swap in `R/vis_corrected/figures/` and update every caption with the new numbers. The CMIP6 comparison is kept (Brian); the zero-TP comparison is kept (Dale).
5. **Feedback factorial (old Fig. 4 and Table):** rerun with the corrected model, or drop it.
   - The new switch definitions freeze each dependence at its 2020 value, which keeps every run mass-consistent.
   - Some switch combinations currently drive CO2 to zero and crash. This needs debugging before the 64 runs can be produced.
6. **Chapter 3:** re-solved on 17 Nov (see section 4b and `chapter3/`). Next, rewrite the chapter text with the new numbers. All Chapter 3 numbers in the current draft use the faulty Lenton ocean and the old permafrost module.
7. **Consistency check:** every number in the text, tables and captions traced to the CSVs in `R/vis_corrected/`.
8. **Code and data:** commit this code. Mark the STELLA `.stmx` as superseded in the public repository.

### Before submitting the revised paper (Communications Sustainability, ~12 weeks)
1. **Point-by-point response letter,** plus marked-up and clean manuscripts.
2. **Framing and title** around the SCC bounds. Engage the literature Reviewer 3 cited: Kessler 2017, Wirths et al. 2018, Dietz et al. 2021a/b, FRIDA-Clim (Wells et al. 2026), Lade et al. 2018. Add Reviewer 2's normative point on discounting and welfare.
3. **Reviewer 1 items:**
   - permafrost magnitude now matches Hector and AR6 (done);
   - gradual-only run (available);
   - emissions figure (done);
   - ppm axis (done);
   - "balancing" vs "reinforcing", CO2 capitalisation, line numbers, describing the cooling mechanism in the introduction.
4. **Reviewer 3 items:**
   - DFAIR + permafrost (done); SCC (done); horizon cut at 2300 (done); early-cooling explanation (done);
   - damage coefficient only rescales results: concede and reword;
   - Howard–Sterner description: "upper end of non-catastrophic", plus a note on collapsing their functional form to a quadratic;
   - consumption-equivalent welfare (optional; the SCC partly answers this);
   - tipping-point definition;
   - rename "SCM";
   - STP terminology;
   - "continuous-time" wording;
   - "three reasons" vs five subsections;
   - missing citations;
   - Dietz 25% vs 8.4%;
   - published version of reference [42].
5. **CMIP6/ScenarioMIP figure:** Reviewer 3 objects because ScenarioMIP runs are concentration-driven. Consider replacing or supplementing it with an emissions-driven comparison. An SSP emissions-driven check against MAGICC concentrations exists in `ssp_check.py`; it needs re-running with the final K_M = 120 configuration and the RCMIP data (not included; download from the RCMIP v5.1.0 site).
6. **Methane:** permafrost CH4 is treated as carbon at CO2 forcing. State this as a limitation, or add CH4 forcing (Woodard et al. report CH4 adds 24–29% to the permafrost temperature effect).
7. **Data availability:** update the public `DICE_eco_feedbacks` repository with the Python model and the outputs.

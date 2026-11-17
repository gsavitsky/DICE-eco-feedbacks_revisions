"""
Chapter 3 outcomes with the corrected physics (2026-11-17). Supersedes chapter4_outcomes_table.py for the
corrected model. Reporting window 2020-2300 (as Chapter 2); the optimisations themselves solve to 2470.

Optimised solutions (Pyomo, 5 substeps per year; 2-substep runs and a second start agree):
  DICE 2023 : dice2023_annual_optimal_controls.json + trajectory_dice2023.json (unchanged model)
  DICE-C    : dicec_corr_ss5.json  + traj_dicec_corr_ss5.json   (dice_pyomo_corrected.py)
  DICE-CP   : dicecp_corr_ss5.json + traj_dicecp_corr_ss5.json
GAMS-control runs and all SCCs: permafrost_revision_test/ch2_models.py forward model (same physics), RK4 dt=0.05,
U_acc to 2300.

Outputs: ch3_corrected_outcomes.csv/.md, figdata_corrected/*.csv (same layout as figdata/, read by
ch3_figures_corrected.R)
"""
import json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "code", "permafrost_revision_test"))   # bundle layout
import ch2_models as M
import permafrost_revision_sim as prs
os.chdir(HERE)

NY = 280                                   # report 2020-2300
YRS = np.arange(2020, 2021 + NY)
GTCO2_PER_GTC = 3.666
FILES = {"DICE 2023": ("results/dice2023_annual_optimal_controls.json", "results/trajectory_dice2023.json"),
         "DICE-C": ("results/dicec_corr_ss5.json", "results/traj_dicec_corr_ss5.json"),
         "DICE-CP": ("results/dicecp_corr_ss5.json", "results/traj_dicecp_corr_ss5.json")}
CARBON = {"DICE 2023": dict(carbon='dfair'), "DICE-C": dict(carbon='lenton'),
          "DICE-CP": dict(carbon='lenton', perma=dict(m=1.0, abrupt='smooth'))}
CFGS = list(FILES)
STD = dict(a2=prs.DAMAGE['std'], rartp=prs.DISC['std'])

miuup_p = np.zeros(81); pv = np.arange(1, 82)
miuup_p[0], miuup_p[1] = 0.05, 0.10
miuup_p[pv > 2] = 0.12 * (pv[pv > 2] - 1)
miuup_p[pv > 8] = 0.85 + 0.05 * (pv[pv > 8] - 8)
miuup_p[pv > 11], miuup_p[pv > 20], miuup_p[pv > 37], miuup_p[pv > 57] = 1.0, 1.1, 1.05, 1.0
MIUUP = np.array([miuup_p[min(1 + n // 5, 80) - 1] for n in range(400)])


def load(cfg):
    c, t = (json.load(open(f)) for f in FILES[cfg])
    d = {k: np.array(t[k][:NY + 1], float) for k in t if k != "years" and isinstance(t[k], list) and len(t[k]) > NY}
    d["TATM"] = d["Tbox1"] + d["Tbox2"]
    if "Ca" not in d:
        d["Ca"] = 588.0 + d["R0"] + d["R1"] + d["R2"] + d["R3"]
    d["miu"], d["s"] = np.array(c["miu"], float)[:400], np.array(c["s"], float)[:400]
    return d


def sim(cfg, ctrl=None, pulse=None):
    c = dict(CARBON[cfg], **STD)
    if ctrl is not None:
        c['ctrl'] = ctrl
    return M.run(c, pulse=pulse, nyears=NY)


def scc(cfg, ctrl, year):
    b = sim(cfg, ctrl); p = sim(cfg, ctrl, (year, 1.0)); i = year - 2020
    mu = [1000 * prs.scale1_ * b['cons'][j][1] * b['cons'][j][0] ** (-prs.eta_) for j in (i, i + 1)]
    return -1000 * (p['Uacc'][-1] - b['Uacc'][-1]) / (0.5 * (mu[0] + mu[1])), b


D = {c: load(c) for c in CFGS}
ix = lambda y: int(y - 2020)

# forward model: optimised controls (consistency + SCC) and GAMS controls (comparison)
F_opt, F_gams, SCC = {}, {}, {}
for c in CFGS:
    ctrl = (D[c]["miu"], D[c]["s"])
    s20, F_opt[c] = scc(c, ctrl, 2020); s50, _ = scc(c, ctrl, 2050)
    g20, F_gams[c] = scc(c, None, 2020); g50, _ = scc(c, None, 2050)
    SCC[c] = dict(opt2020=s20, opt2050=s50, gams2020=g20, gams2050=g50)
    print(f"{c}: forward model with optimised controls peak T {F_opt[c]['TATM'].max():.3f} vs Pyomo {D[c]['TATM'].max():.3f}; "
          f"SCC2020 opt {s20:.1f} / GAMS controls {g20:.1f}", flush=True)


def first_year(mask):
    w = np.where(mask)[0]
    return int(YRS[w[0]]) if len(w) else None


rows = []
def add(metric, vals, unit="", prec=2):
    r = dict(metric=metric, unit=unit)
    for c in CFGS:
        v = vals.get(c)
        r[c] = "" if v is None else (f"{v:.{prec}f}" if isinstance(v, float) else str(v))
    rows.append(r)

g = lambda c, k, y: float(D[c][k][ix(y)])
sig = np.array([prs.exo(float(y))[2] for y in YRS]); aL = np.array([prs.exo(float(y))[1] for y in YRS])
L = np.array([prs.exo(float(y))[0] for y in YRS])
def eind(c):        # industrial CO2 GtCO2/yr at year start, that year's control
    yg = aL * (L / 1000.0) ** 0.7 * D[c]["K"] ** 0.3
    return sig[:NY] * yg[:NY] * (1 - D[c]["miu"][:NY])

# policy
add("Year mu first reaches 1.0", {c: first_year(D[c]["miu"][:NY] >= 1.0 - 1e-6) for c in CFGS}, "year", 0)
for y in (2030, 2050, 2075, 2100):
    add(f"mu {y}", {c: g(c, "miu", y) for c in CFGS}, "", 3)
add("Cumulative industrial CO2 2020-2100", {c: float(eind(c)[:80].sum() / GTCO2_PER_GTC) for c in CFGS}, "GtC", 0)
# warming
pk = {c: int(np.argmax(D[c]["TATM"])) for c in CFGS}
add("Peak TATM", {c: float(D[c]["TATM"][pk[c]]) for c in CFGS}, "C", 3)
add("Peak TATM year", {c: int(YRS[pk[c]]) for c in CFGS}, "year", 0)
for y in (2050, 2100, 2200, 2300):
    add(f"TATM {y}", {c: g(c, "TATM", y) for c in CFGS}, "C", 2)
for thr in (1.5, 2.0, 2.5):
    add(f"Years above {thr} C (2020-2300)", {c: int(np.sum(D[c]["TATM"][:NY] > thr)) for c in CFGS}, "years", 0)
add("Peak atmospheric CO2", {c: float(D[c]["Ca"].max() / 2.1295) for c in CFGS}, "ppm", 0)
# permafrost
cp = D["DICE-CP"]
for y in (2050, 2100, 2200, 2300):
    add(f"Frozen permafrost carbon {y}", {"DICE-CP": g("DICE-CP", "Cperm", y)}, "PgC", 0)
add("Thawed permafrost carbon peak", {"DICE-CP": float(cp["thawedC"].max())}, "PgC", 1)
add("Permafrost release to atmosphere 2020-2100 (forward model)", {"DICE-CP": float(F_opt["DICE-CP"]['release'][80])}, "PgC", 0)
# welfare and SCC
U = {c: g(c, "Uacc", 2300) for c in CFGS}
add("U_acc(2300), optimised", U, "", 2)
add("U_acc(2300) % vs DICE 2023 optimised", {c: 100 * (U[c] / U["DICE 2023"] - 1) for c in CFGS if c != "DICE 2023"}, "%", 3)
add("U_acc(2300) gain from re-optimising vs GAMS controls (forward model)",
    {c: 100 * (F_opt[c]['Uacc'][-1] / F_gams[c]['Uacc'][-1] - 1) for c in CFGS}, "%", 3)
add("SCC 2020 at own optimum", {c: SCC[c]['opt2020'] for c in CFGS}, "$/tCO2", 1)
add("SCC 2050 at own optimum", {c: SCC[c]['opt2050'] for c in CFGS}, "$/tCO2", 1)
add("SCC 2020 at GAMS controls", {c: SCC[c]['gams2020'] for c in CFGS}, "$/tCO2", 1)
add("SCC 2020 at own optimum, % vs DICE 2023", {c: 100 * (SCC[c]['opt2020'] / SCC['DICE 2023']['opt2020'] - 1) for c in CFGS if c != "DICE 2023"}, "%", 1)
# optimised vs GAMS controls (policy effect)
for c in CFGS:
    add(f"{c}: peak TATM, optimised minus GAMS controls", {c: float(D[c]["TATM"].max() - F_gams[c]['TATM'].max())}, "C", 3)

tab = pd.DataFrame(rows)
tab.to_csv("ch3_corrected_outcomes.csv", index=False)
md = ["| " + " | ".join(tab.columns) + " |", "|" + "---|" * len(tab.columns)] + \
     ["| " + " | ".join(str(x) for x in r) + " |" for r in tab.astype(str).values]
open("ch3_corrected_outcomes.md", "w").write("# Chapter 3 outcomes, corrected physics (2020-2300)\n\n" + "\n".join(md) + "\n")
pd.set_option("display.width", 220, "display.max_colwidth", 70)
print(tab.to_string(index=False))

# ---------------- figure data (same layout as figdata/) ----------------
FD = "figdata_corrected"; os.makedirs(FD, exist_ok=True)
ts = []
for c in CFGS:
    d = D[c]; e = eind(c)
    for i, y in enumerate(YRS):
        ts.append(dict(config=c, year=int(y), miu=float(d["miu"][i]) if i < NY else np.nan,
                       s=float(d["s"][i]) if i < NY else np.nan, miu_ceiling=float(MIUUP[i]) if i < NY else np.nan,
                       eind_gtco2=float(e[i]) if i < NY else np.nan, TATM=float(d["TATM"][i]), Ca=float(d["Ca"][i]),
                       K=float(d["K"][i]), Uacc=float(d["Uacc"][i]),
                       Cperm=float(d["Cperm"][i]) if "Cperm" in d else np.nan,
                       thawedC=float(d["thawedC"][i]) if "thawedC" in d else np.nan))
pd.DataFrame(ts).to_csv(f"{FD}/ch4_timeseries.csv", index=False)
gp, op = F_gams["DICE-CP"], D["DICE-CP"]
pd.DataFrame(dict(year=YRS, TATM_ch3=gp['TATM'], TATM_ch4_model_gams_controls=gp['TATM'], TATM_ch4_optimized=op["TATM"],
                  Cperm_ch4_model_gams_controls=gp['Cperm'], Cperm_ch4_optimized=op["Cperm"],
                  thawedC_ch4_model_gams_controls=gp['thawed'], thawedC_ch4_optimized=op["thawedC"],
                  Uacc_ch3=gp['Uacc'], Uacc_ch4_model_gams_controls=gp['Uacc'], Uacc_ch4_optimized=op["Uacc"])
             ).to_csv(f"{FD}/ch4_ch3_vs_ch4_dicecp.csv", index=False)
pd.DataFrame(dict(year=YRS, dice2023_opt=D["DICE 2023"]["Uacc"], dicec_opt=D["DICE-C"]["Uacc"], dicec_gams=F_gams["DICE-C"]['Uacc'],
                  dicecp_opt=op["Uacc"], dicecp_gams=gp['Uacc'])).to_csv(f"{FD}/ch4_uacc_opt_vs_gams.csv", index=False)
pd.DataFrame(dict(year=YRS, dice2023_opt=D["DICE 2023"]["TATM"], dicec_opt=D["DICE-C"]["TATM"], dicec_gams=F_gams["DICE-C"]['TATM'],
                  dicecp_opt=op["TATM"], dicecp_gams=gp['TATM'])).to_csv(f"{FD}/ch4_tatm_opt_vs_gams.csv", index=False)
RUNS = [("DICE-C", 2, "GAMS-controls start", "results/dicec_corr_ss2_gams.json", "results/traj_dicec_corr_ss2_gams.json"),
        ("DICE-C", 2, "previous-optimum start", "results/dicec_corr_ss2_alt.json", "results/traj_dicec_corr_ss2_alt.json"),
        ("DICE-C", 5, "warm start from the 2-substep solution", "results/dicec_corr_ss5.json", "results/traj_dicec_corr_ss5.json"),
        ("DICE-CP", 2, "GAMS-controls start", "results/dicecp_corr_ss2_gams.json", "results/traj_dicecp_corr_ss2_gams.json"),
        ("DICE-CP", 2, "previous-optimum start", "results/dicecp_corr_ss2_alt.json", "results/traj_dicecp_corr_ss2_alt.json"),
        ("DICE-CP", 5, "warm start from the 2-substep solution", "results/dicecp_corr_ss5.json", "results/traj_dicecp_corr_ss5.json")]
rr = []
for cfg, ss, start, cf, tf in RUNS:
    ctl, tr = json.load(open(cf)), json.load(open(tf))
    for i in range(NY):
        rr.append(dict(config=cfg, substeps=ss, start=start, year=2020 + i, miu=float(ctl["miu"][i]), s=float(ctl["s"][i]),
                       TATM=float(tr["Tbox1"][i] + tr["Tbox2"][i]), Uacc2420=float(tr["Uacc"][NY]), final=(cf == FILES[cfg][0])))
pd.DataFrame(rr).to_csv(f"{FD}/ch4_runs_resolution_multistart.csv", index=False)
json.dump(SCC, open(f"{FD}/ch3_scc.json", "w"), indent=1)
print(f"wrote {FD}/")

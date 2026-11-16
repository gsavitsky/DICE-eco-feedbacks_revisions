"""
Hector-faithful permafrost test (2026-11-16). FORWARD SIMULATION ONLY, GAMS controls, no optimisation.

Reuses permafrost_revision_sim.py (imported, not modified) and swaps in a permafrost module that follows
the Hector v3 source code (JGCRI/hector main, src/simpleNbox-runtime.cpp, rh_ftpa_co2 /
compute_pf_thaw_refreeze / slowparameval), which is what Woodard et al. (2021) results are consistent with:

  respiration   J = thawedC * (1 - fstatic) * 0.02 * Q10^(T200/10)       (CO2 + CH4 together)
                -> the 0.02 /yr base rate is MISSING in the Chapter 2 STELLA model (effective rate 1 /yr)
  frozen frac   f = 1 - Phi( (ln(wf*T) - mu) / sigma ),  mu = 1.67, sigma = 0.986, wf = 2
                (Chapter 2 uses ln(1 + wf*T), sigma = 0.99)
  T200          200-yr running mean of wf*T (Hector uses land temperature * wf; here global T * wf),
                and the Q10 effect is 'sticky' (never declines), as in Hector
  pool          730 PgC in 2020 (Hector's modern pool, Woodard Table 2; 865 PgC is Hector's 1745 value)
  refreeze      unchanged limiter (thawed pool can only give back what it holds)

Pre-2020 temperature history for the running mean: T(y) = 1.25*(exp(k(y-1850))-1)/(exp(170k)-1),
k = 0.02, zero before 1850 (a simple reconstruction; documented assumption).

Configurations:
  DICE-CP old + 0.02 rate      : Chapter 2 module, only the respiration base rate added (isolates the main fix)
  DICE-CP Hector gradual-only  : all Hector-faithful changes, no abrupt thaw
  DICE-CP Hector + smooth AMcK : all Hector-faithful changes + smooth Armstrong McKay abrupt term
Outputs: *_hector.csv / .json / .png in this folder (no existing file is overwritten).
"""
import csv, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import permafrost_revision_sim as prs  # import only; its __main__ block does not run

OUT = HERE
SQ2PI = math.sqrt(2 * math.pi)
WF, MU, SIG_H, RATE = 2.0, 1.67, 0.986, 0.02
CUR = {}


def dff_dT_hector(T):
    if T <= 1e-9:
        return 0.0
    z = (math.log(WF * T) - MU) / SIG_H
    return -math.exp(-0.5 * z * z) / SQ2PI / (SIG_H * T)


def perma_flows_cfg(T, dT, Cperm, thawedC, A, G, T200, mode, Sst=0.0, inert=False):
    cfg = CUR
    dff = dff_dT_hector(T) if cfg.get('hector_ff') else prs.dff_dT(T)
    gradual_raw = -dff * dT * Cperm
    lim = -thawedC / prs.REFREEZE_TAU
    gradual = 0.5 * (gradual_raw + lim + math.sqrt((gradual_raw - lim) ** 2 + 1e-9))
    if mode == 'step':
        abrupt = 0.4 * max(0.0, gradual_raw) if T > 2.0 else 0.0
    elif mode == 'smooth':
        abrupt = prs.smax0(prs.AB_FRAC * prs.w_onset(T) * prs.smax0(G) - A) / prs.AB_TAU
    else:
        abrupt = 0.0
    tot = gradual + abrupt
    J = (1 - prs.fstatic_) * thawedC * cfg.get('rate', 1.0) * prs.Q10_ ** (T200 / 10.0)
    return (-tot, tot - J, abrupt, gradual, J, J, 0.0)


prs.perma_flows = perma_flows_cfg  # rhs() looks this up as a module global


def T_hist(y, k=0.02):
    if y < 1850:
        return 0.0
    return 1.25 * (math.exp(k * (y - 1850)) - 1) / (math.exp(170 * k) - 1)


def run_h(cfg, dt=prs.DT):
    CUR.clear(); CUR.update(cfg)
    n_sub = int(round(1 / dt)); lag = int(round(200 / dt))
    S = np.array([0.1477, 1.099454, 295.0, 633.5, 0.0, 0.518, prs.P['Ca0'], prs.P['Cv0'], prs.P['Cs0'],
                  prs.P['Col0'], prs.P['Coh0'], prs.P['Coi0'], prs.P['Cod0'], 0.0,
                  cfg['C0'] if cfg['perma'] else 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    running = cfg.get('t200') == 'runmean'
    if running:  # history buffer of wf*T at substep resolution, 1820-2020
        hist = [WF * T_hist(2020 - 200 + j * dt) for j in range(lag)]
        buf = list(hist); csum = float(np.sum(hist)); t200_max = 0.0
    else:
        buf = []
    out, minS = [], np.full(len(S), np.inf)
    for kk in range(prs.N * n_sub + 1):
        y = 2020 + kk * dt
        if kk % n_sub == 0:
            out.append(S.copy())
        if kk == prs.N * n_sub:
            break
        if running:
            T200 = csum / lag
            t200_max = max(t200_max, T200); T200 = t200_max   # sticky Q10, as in Hector
        else:
            buf.append(S[0] + S[1])
            T200 = buf[kk - lag] if kk >= lag else prs.t200_init(y)
        k1 = prs.rhs(S, y, T200, cfg); k2 = prs.rhs(S + dt / 2 * k1, y + dt / 2, T200, cfg)
        k3 = prs.rhs(S + dt / 2 * k2, y + dt / 2, T200, cfg); k4 = prs.rhs(S + dt * k3, y + dt, T200, cfg)
        S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        minS = np.minimum(minS, S)
        if running:
            new = WF * (S[0] + S[1]); csum += new - buf[kk]; buf.append(new)
    A = np.array(out)
    d = {nm: A[:, i] for i, nm in enumerate(prs.NAMES)}
    d['TATM'] = d['Tbox1'] + d['Tbox2']
    d['min_thawedC'] = float(minS[15]); d['min_Cperm'] = float(minS[14])
    return d


MODELS = {
    'DICE-C': dict(perma=False, C0=0.0, abrupt='none'),
    'DICE-CP old': dict(perma=True, C0=1035.0, abrupt='step'),
    'DICE-CP old + 0.02 rate': dict(perma=True, C0=1035.0, abrupt='step', rate=RATE),
    'DICE-CP Hector gradual-only': dict(perma=True, C0=730.0, abrupt='none', rate=RATE, hector_ff=True, t200='runmean'),
    'DICE-CP Hector + smooth AMcK': dict(perma=True, C0=730.0, abrupt='smooth', rate=RATE, hector_ff=True, t200='runmean'),
}

if __name__ == "__main__":
    t0 = time.time()
    ST = prs.stella_xlsx()
    yrs = np.arange(2020, 2421)
    cells = [(dam, disc) for disc in ('std', 'zeroTP') for dam in ('std', 'HS')]
    res = {}
    for name, m in MODELS.items():
        for dam, disc in cells:
            cfg = dict(m, a2=prs.DAMAGE[dam], rartp=prs.DISC[disc])
            res[(name, dam, disc)] = run_h(cfg)
            print(f"  {name:30s} {dam:3s} {disc:6s} U_acc {res[(name, dam, disc)]['Uacc'][-1]:14.2f}", flush=True)

    base = {'std': ST[1]['Uacc'][-1], 'zeroTP': ST[2]['Uacc'][-1]}          # DICE 2023 std damage, STELLA
    stella_old = {(dam, disc): ST[prs.stella_run(13, dam, disc)]['Uacc'][-1] for dam, disc in cells}
    rows = []
    for (name, dam, disc), d in res.items():
        T = d['TATM']; i = int(np.argmax(T))
        first = lambda th: (int(yrs[np.where(T > th)[0][0]]) if np.any(T > th) else None)
        last = lambda th: (int(yrs[np.where(T > th)[0][-1]]) if np.any(T > th) else None)
        u = d['Uacc'][-1]
        # STELLA-anchored: for permafrost models, STELLA old value + (python new - python old)
        u_anch = u
        if MODELS[name]['perma']:
            u_anch = stella_old[(dam, disc)] + (u - res[('DICE-CP old', dam, disc)]['Uacc'][-1])
        elif name == 'DICE-C':
            u_anch = ST[prs.stella_run(9, dam, disc)]['Uacc'][-1]
        row = dict(model=name, damage=dam, disc=disc, peak_T=T[i], peak_year=int(yrs[i]),
                   T2100=T[80], T2200=T[180], T2300=T[280], T2420=T[400],
                   first_above_2=first(2.0), first_above_3=first(3.0), last_above_2=last(2.0), last_above_1p5=last(1.5),
                   years_above_1p5=int(np.sum(T > 1.5)), years_above_2=int(np.sum(T > 2.0)),
                   uacc_2420_py=u, uacc_2420_anchored=u_anch,
                   loss_pct_vs_DICE2023_std=100 * (u_anch / base[disc] - 1))
        if MODELS[name]['perma']:
            dc = res[('DICE-C', dam, disc)]
            rel = d['R'][80]
            row.update(release_2100=rel, release_2200=d['R'][180], release_2300=d['R'][280], release_2420=d['R'][400],
                       thaw_2100=d['G'][80] + d['A'][80], abrupt_2100=d['A'][80],
                       PgC_per_C_2100=rel / T[80],
                       Ca_excess_over_DICEC_2100=d['Ca'][80] - dc['Ca'][80],
                       Ca_excess_per_C_2100=(d['Ca'][80] - dc['Ca'][80]) / T[80],
                       dT_vs_DICEC_2100=T[80] - dc['TATM'][80], dTpeak_vs_DICEC=T[i] - float(np.max(dc['TATM'])),
                       min_thawedC=d['min_thawedC'],
                       mass_err_max=float(np.max(np.abs(d['Cperm'] + d['thawedC'] + d['R'] - MODELS[name]['C0']))))
        rows.append(row)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(OUT, "summary_hector.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()})
    with open(os.path.join(OUT, "timeseries_std_hector.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["model", "year", "TATM", "Ca", "Cperm", "thawedC", "A", "G", "release_cum"])
        for name in MODELS:
            d = res[(name, 'std', 'std')]
            for j, y in enumerate(yrs):
                w.writerow([name, y] + [f"{d[k][j]:.6f}" for k in ("TATM", "Ca", "Cperm", "thawedC", "A", "G", "R")])

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for name in MODELS:
        d = res[(name, 'std', 'std')]
        ax[0].plot(yrs, d['TATM'], label=name)
        if MODELS[name]['perma']:
            ax[1].plot(yrs, d['R'], label=name)
    ax[0].set_ylabel("TATM (°C)"); ax[1].set_ylabel("Cumulative permafrost release (PgC)")
    for a in ax: a.set_xlabel("Year"); a.grid(alpha=.3)
    ax[0].legend(fontsize=8)
    fig.suptitle("Hector-faithful permafrost test, GAMS controls (standard damage, standard discounting)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "permafrost_hector_comparison.png"), dpi=150)
    print(f"done in {time.time() - t0:.0f}s")

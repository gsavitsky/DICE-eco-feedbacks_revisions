"""
SCC decomposition by carbon-cycle process, Shapley values (forward model, GAMS controls, standard DICE settings).

Lenton factorial over 7 processes (ch2_models cfg['fbx'], 'off' = pre-industrial response from 1800, own spin-up):
  ocean  : carbonate-chemistry buffer saturation       oceanT : ocean warming (solubility/chemistry; Lenton variant iv)
  fert   : CO2 fertilisation                             photoT : temperature response of photosynthesis
  vegT   : temperature response of plant respiration    soilT  : temperature response of soil respiration
  perma  : permafrost (Hector strength, smooth abrupt thaw)
2^7 = 128 combinations. SCC (2019$/tCO2) of a 1 GtCO2 pulse in 2020 and in 2100, welfare to 2200 (the factorial
window; after 2200 the control path's sustained net removal drives Lenton CO2 below pre-industrial).
Shapley value of process i: phi_i = sum over subsets S of the others of |S|!(n-|S|-1)!/n! [v(S+i) - v(S)];
the phi_i add exactly to v(all on) - v(all off).

DFAIR robustness: scale irC / irT so DFAIR's emergent C4MIP total beta and gamma span the CMIP6 mean +- 1 s.d.
(Arora et al. 2020 land + ocean, s.d. combined assuming independence), and report the permafrost SCC premium.
Outputs (vis_corrected/feedback_factorial/): shapley_runs.csv, shapley_values.csv, dfair_cmip6_range.csv.
"""
import os, sys, itertools, math, time
import numpy as np, pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected/feedback_factorial")
PROC = ['ocean', 'oceanT', 'fert', 'photoT', 'vegT', 'soilT', 'perma']
NY = int(os.environ.get('SHAPLEY_NY', '180'))      # welfare horizon: 180 -> 2200, 280 -> 2300
PULSES = tuple(int(x) for x in os.environ.get('SHAPLEY_PULSES', '2020,2100').split(','))
BETA = (0.97 + 0.78, math.hypot(0.40, 0.07)); GAMMA = (-45.07 - 17.21, math.hypot(50.59, 4.95))  # CMIP6 4xCO2, Arora 2020 Table A1


def scc_of(M, b, p, year):
    prs = M.prs; i = year - 2020
    mu = [1000 * prs.scale1_ * b['cons'][j][1] * b['cons'][j][0] ** (-prs.eta_) for j in (i, i + 1)]
    return -1000 * (p['Uacc'][-1] - b['Uacc'][-1]) / (0.5 * (mu[0] + mu[1]))


def job(args):
    sys.path.insert(0, HERE)
    import ch2_models as M
    tag, cfg, irs = args
    C0, T0 = M.SC.P['irC'], M.SC.P['irT']
    if irs:
        M.SC.P['irC'], M.SC.P['irT'] = C0 * irs[0], T0 * irs[1]
    try:
        b = M.run(cfg, nyears=NY)
        out = dict(peak_T=b['TATM'].max(), T2100=b['TATM'][80], T2200=b['TATM'][180])
        for y in PULSES:
            out[f'SCC{y}'] = scc_of(M, b, M.run(cfg, pulse=(y, 1.0), nyears=NY), y)
    finally:
        M.SC.P['irC'], M.SC.P['irT'] = C0, T0
    return tag, out


def dfair_c4mip(sC, sT):
    sys.path.insert(0, HERE)
    import ch2_models as M, ch2_dfair_feedbacks as D
    C0, T0 = M.SC.P['irC'], M.SC.P['irT']
    M.SC.P['irC'], M.SC.P['irT'] = C0 * sC, T0 * sT
    try:
        ub, dppm, dT = D.c4mip_dfair(False); uc, _, _ = D.c4mip_dfair(True)
    finally:
        M.SC.P['irC'], M.SC.P['irT'] = C0, T0
    return ub / dppm, (uc - ub) / dT


def solve_scales(target):
    """irC scale for beta, then irT scale for gamma at that irC scale (beta does not depend on irT)."""
    from scipy.optimize import brentq
    b, g = target
    sC = brentq(lambda s: dfair_c4mip(s, 1.0)[0] - b, 0.0, 6.0, xtol=1e-3)
    sT = brentq(lambda s: dfair_c4mip(sC, s)[1] - g, 0.0, 6.0, xtol=1e-3)
    return dict(beta_target=b, gamma_target=g, irC_scale=sC, irT_scale=sT)


def shapley(v, n):
    phi = np.zeros(n)
    for i in range(n):
        others = [j for j in range(n) if j != i]
        for k in range(n):
            w = math.factorial(k) * math.factorial(n - k - 1) / math.factorial(n)
            for S in itertools.combinations(others, k):
                on = set(S)
                key_without = tuple(1 if j in on else 0 for j in range(n))
                key_with = tuple(1 if (j in on or j == i) else 0 for j in range(n))
                phi[i] += w * (v[key_with] - v[key_without])
    return phi


if __name__ == "__main__":
    import permafrost_revision_sim as prs
    import ch2_models as M
    t0 = time.time()
    std = dict(a2=0.003467, rartp=prs.DISC['std']); PF = dict(m=1.0, abrupt='smooth')
    jobs = []
    for c in itertools.product((0, 1), repeat=7):
        fbx = {k: v for k, v in zip(PROC[:6], c[:6])}
        cfg = dict(std, carbon='lenton', fbx=fbx, KM=float(os.environ.get('SHAPLEY_KM', M.KM)))
        if c[6]: cfg['perma'] = PF
        jobs.append((('fac',) + c, cfg, None))
    with Pool(8) as pool:
        SKIP = os.environ.get('SHAPLEY_SKIP_DFAIR') == '1'
        corners = [] if SKIP else [(b, g) for b in (BETA[0] - BETA[1], BETA[0] + BETA[1]) for g in (GAMMA[0] + GAMMA[1], GAMMA[0] - GAMMA[1])]
        if not SKIP: corners.append(BETA[:1] + GAMMA[:1])
        sc = pool.map(solve_scales, corners) if corners else []
        if not SKIP: sc.append(dict(beta_target=np.nan, gamma_target=np.nan, irC_scale=1.0, irT_scale=1.0))   # DICE 2023 default
        for k, s in enumerate(sc):
            for m in (0, 1):
                cfg = dict(std, carbon='dfair', **({'perma': PF} if m else {}))
                jobs.append((('dfair', k, m), cfg, (s['irC_scale'], s['irT_scale'])))
        res = dict(pool.map(job, jobs, chunksize=1))
        emer = pool.starmap(dfair_c4mip, [(s['irC_scale'], s['irT_scale']) for s in sc])
    print(f"runs done [{time.time() - t0:.0f}s]", flush=True)

    # ---- Lenton factorial + Shapley ----
    R = pd.DataFrame([dict(zip(PROC, k[1:]), **v) for k, v in res.items() if k[0] == 'fac'])
    TAG = os.environ.get("SHAPLEY_TAG", "")
    R.to_csv(os.path.join(OUT, f"shapley_runs{TAG}.csv"), index=False)
    metrics = [f'SCC{y}' for y in PULSES] + ['peak_T', 'T2100', 'T2200']
    rows = []
    for mtr in metrics:
        v = {tuple(int(x) for x in r[PROC]): r[mtr] for _, r in R.iterrows()}
        phi = shapley(v, len(PROC))
        allon, alloff = v[(1,) * 7], v[(0,) * 7]
        for p, ph in zip(PROC, phi):
            rows.append(dict(metric=mtr, process=p, shapley=ph, all_off=alloff, all_on=allon, check_sum=phi.sum(),
                             check_total=allon - alloff))
    S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, f"shapley_values{TAG}.csv"), index=False)

    # ---- DFAIR across the CMIP6 beta/gamma range ----
    drows = []
    for k, s in enumerate(sc):
        b0, b1 = res[('dfair', k, 0)], res[('dfair', k, 1)]
        drows.append(dict(s, beta_emergent=emer[k][0], gamma_emergent=emer[k][1],
                          SCC2020_noP=b0['SCC2020'], SCC2020_P=b1['SCC2020'], premium_pct=100 * (b1['SCC2020'] / b0['SCC2020'] - 1),
                          perma_dpeak=b1['peak_T'] - b0['peak_T'], perma_dT2200=b1['T2200'] - b0['T2200']))
    DR = pd.DataFrame(drows)
    if not SKIP: DR.to_csv(os.path.join(OUT, "dfair_cmip6_range.csv"), index=False)

    pd.set_option('display.width', 220)
    print(S.pivot(index='process', columns='metric', values='shapley').loc[PROC, metrics].round(3).to_string())
    print(S.groupby('metric')[['all_off', 'all_on', 'check_sum', 'check_total']].first().round(3).to_string())
    print(DR.round(3).to_string(index=False))
    print(f"[{time.time() - t0:.0f}s]")

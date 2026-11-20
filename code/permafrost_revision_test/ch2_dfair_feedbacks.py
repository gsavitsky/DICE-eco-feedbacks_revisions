"""
DFAIR feedback strengths (forward model, GAMS controls, standard DICE settings, 2020-2200), and a like-for-like
comparison with the Lenton factorial.

DFAIR's carbon-cycle feedbacks enter only through iIRF100 = r0 + irC * (cumulative uptake) + irT * T.
 1. Scale irC (uptake saturation) and irT (climate-carbon feedback) by S in SCALES, permafrost off / Hector strength.
 2. Emergent C4MIP beta (PgC/ppm) and gamma (PgC/K), 1%/yr CO2 to 4x from pre-industrial, total land + ocean:
    DFAIR (compatible emissions from the prescribed concentration path) and Lenton (variant and K_M from ch2_models), vs CMIP6
    (Arora et al. 2020, 4xCO2: land + ocean beta 1.76, gamma -62.3 PgC/K).
 3. Grouped Lenton effects from the factorial: concentration feedbacks (ocean buffer + fertilisation) and temperature
    feedbacks (photosynthesis, plant and soil respiration), switched together.
Outputs: vis_corrected/feedback_factorial/dfair_scaling_runs.csv, dfair_scaling_timeseries.csv, c4mip_compare.csv,
lenton_grouped_effects.csv.
"""
import os, sys, math, time
import numpy as np, pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected/feedback_factorial")
SCALES = [0.0, 0.5, 1.0, 1.5, 2.0]
NY = 180
YRS = np.arange(2020, 2020 + NY + 1)


def dfair_job(args):
    sys.path.insert(0, HERE)
    import ch2_models as M
    sC, sT, m = args
    C0, T0 = M.SC.P['irC'], M.SC.P['irT']
    M.SC.P['irC'], M.SC.P['irT'] = C0 * sC, T0 * sT
    try:
        cfg = dict(carbon='dfair', a2=0.003467, rartp=M.prs.DISC['std'])
        if m:
            cfg['perma'] = dict(m=1.0, abrupt='smooth')
        d = M.run(cfg, nyears=NY)
    finally:
        M.SC.P['irC'], M.SC.P['irT'] = C0, T0
    return args, d['TATM'], d['CO2ppm']


def c4mip_dfair(coupled, years=140, dt=0.05, TCR=1.8):
    import ch2_models as M
    SC = M.SC; Meq = SC.Meq
    R = np.zeros(4); cumE = 0.0
    for k in range(int(years / dt)):
        t = (k + 1) * dt
        ppm_new = 280 * 1.01 ** t; MATnew = Meq * ppm_new / 280.0
        T = TCR * math.log2(ppm_new / 280) if coupled else 0.0
        a = SC.alpha_of(cumE - R.sum(), T)
        decay = R / (np.array(SC.TAU) * a)
        # implicit-enough: choose E so that sum(R) reaches MATnew - Meq after this step
        E = ((MATnew - Meq) - (R - dt * decay).sum()) / (dt * sum(SC.ESH))
        R = R + dt * (np.array(SC.ESH) * E - decay); cumE += E * dt
    dppm = 280 * 1.01 ** years - 280
    return cumE - R.sum(), dppm, TCR * math.log2(1.01 ** years)


def lenton_c4mip(coupled):
    sys.path.insert(0, HERE)
    import lenton_c4mip_check as C
    import ch2_models as M
    dl, do, dppm, dT = C.run(M.KM, coupled, ocean_warming=M.OW)
    return coupled, dl, do, dppm, dT


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(sC, sT, m) for sC in SCALES for sT in SCALES for m in (0, 1)]
    with Pool(8) as pool:
        lres = pool.map_async(lenton_c4mip, [False, True])
        res = {a: (T, C) for a, T, C in pool.map(dfair_job, jobs, chunksize=1)}
        lc = {c: r for c, *r in lres.get()}
    rows, tsl = [], []
    for (sC, sT, m), (T, C) in res.items():
        T0 = res[(sC, sT, 0)][0]
        rows.append(dict(irC_scale=sC, irT_scale=sT, perma=m, peak_T=T.max(), T2100=T[80], T2150=T[130], T2200=T[180],
                         CO2_2100=C[80], CO2_2200=C[180], perma_dpeak_T=T.max() - T0.max(), perma_dT2100=T[80] - T0[80],
                         perma_dT2200=T[180] - T0[180]))
        tsl.append(pd.DataFrame(dict(year=YRS, TATM=T, CO2ppm=C, irC_scale=sC, irT_scale=sT, perma=m)))
    D = pd.DataFrame(rows).sort_values(['perma', 'irC_scale', 'irT_scale'])
    D.to_csv(os.path.join(OUT, "dfair_scaling_runs.csv"), index=False)
    pd.concat(tsl).to_csv(os.path.join(OUT, "dfair_scaling_timeseries.csv"), index=False)

    # C4MIP beta / gamma
    import ch2_models as M
    ub, dppm, dT = c4mip_dfair(False); uc, _, _ = c4mip_dfair(True)
    lb = lc[False]; lcp = lc[True]
    C4 = pd.DataFrame([
        dict(model='DFAIR (DICE 2023)', beta_total=ub / dppm, gamma_total=(uc - ub) / dT),
        dict(model=f'Lenton ({M.VARIANT}, K_M={M.KM:g})', beta_total=(lb[0] + lb[1]) / dppm, gamma_total=((lcp[0] + lcp[1]) - (lb[0] + lb[1])) / dT,
             beta_land=lb[0] / dppm, beta_ocean=lb[1] / dppm, gamma_land=(lcp[0] - lb[0]) / dT, gamma_ocean=(lcp[1] - lb[1]) / dT),
        dict(model='CMIP6 (Arora 2020) mean, 4xCO2', beta_total=0.97 + 0.78, gamma_total=-45.07 - 17.21, beta_land=0.97,
             beta_ocean=0.78, gamma_land=-45.07, gamma_ocean=-17.21)])
    C4.to_csv(os.path.join(OUT, "c4mip_compare.csv"), index=False)

    # Lenton grouped effects from the factorial (contribution in the full model)
    F = pd.read_csv(os.path.join(OUT, "factorial_runs.csv"))
    grp = {'concentration (ocean buffer + fertilisation)': ['ocean', 'fert'],
           'temperature (photosynthesis + plant + soil respiration)': ['photoT', 'vegT', 'soilT']}
    FB = ['ocean', 'fert', 'photoT', 'vegT', 'soilT']
    g = []
    for p in (0, 1):
        on = F[(F[FB] == 1).all(axis=1) & (F.perma == p)].iloc[0]
        for name, cols in grp.items():
            rest = [c for c in FB if c not in cols]
            off = F[(F[cols] == 0).all(axis=1) & (F[rest] == 1).all(axis=1) & (F.perma == p)].iloc[0]
            g.append(dict(perma=p, group=name, **{k: on[k] - off[k] for k in ('peak_T', 'T2100', 'T2200')}))
    G = pd.DataFrame(g); G.to_csv(os.path.join(OUT, "lenton_grouped_effects.csv"), index=False)

    pd.set_option('display.width', 200)
    base = D[(D.irC_scale == 1) & (D.irT_scale == 1) & (D.perma == 0)].iloc[0]
    print("DFAIR: effect of switching each term off (scale 0) vs default, permafrost off")
    for name, q in (('irC off', (D.irC_scale == 0) & (D.irT_scale == 1)), ('irT off', (D.irC_scale == 1) & (D.irT_scale == 0)),
                    ('both off', (D.irC_scale == 0) & (D.irT_scale == 0))):
        r = D[q & (D.perma == 0)].iloc[0]
        print(f"  {name:9s}: peak {base.peak_T - r.peak_T:+.2f}  2100 {base.T2100 - r.T2100:+.2f}  2200 {base.T2200 - r.T2200:+.2f}")
    print(D[D.perma == 1].pivot(index='irC_scale', columns='irT_scale', values='perma_dpeak_T').round(3).to_string())
    print(D[D.perma == 1].pivot(index='irC_scale', columns='irT_scale', values='perma_dT2200').round(3).to_string())
    print(C4.round(2).to_string(index=False))
    print(G.round(3).to_string(index=False))
    print(f"[{time.time() - t0:.0f}s]")

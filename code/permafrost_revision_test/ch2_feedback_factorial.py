"""
Feedback factorial and CO2-fertilisation sweep for DICE-C / DICE-CP (forward model, GAMS controls, standard DICE
damages and discounting, 2020-2200).
Ends in 2200: after 2200 the GAMS control path's sustained net removal (MIU > 1 from ~2120) drives every
configuration, including the default DICE-C, below pre-industrial CO2, so the late period reflects the control path
rather than the feedbacks (and one no-sink combination reaches zero CO2).

Switches (ch2_models cfg['fbx']) are referenced to pre-industrial: an 'off' process keeps its pre-industrial response
through the 1800-2020 spin-up and after (fert: pre-industrial CO2 in photosynthesis; photoT/vegT/soilT: zero
temperature anomaly in that process; ocean: surface pCO2 linear about the pre-industrial state, i.e. no buffer
saturation). Each combination has its own 2020 carbon state; the 2020 climate state is the same for all.

 1. 2^6 factorial: ocean, fert, photoT, vegT, soilT, permafrost (Hector strength, smooth abrupt thaw). K_M from ch2_models.
 2. Fertilisation sweep: K_M in KMS (land beta from km_beta_map.csv) x permafrost off / m = 1 / m = 4; all else on.
Outputs (vis_corrected/feedback_factorial/): factorial_runs.csv, factorial_effects.csv, factorial_timeseries.csv,
fert_sweep_runs.csv, fert_sweep_timeseries.csv.
"""
import os, sys, itertools, time
import numpy as np, pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected/feedback_factorial")
FB = ['ocean', 'fert', 'photoT', 'vegT', 'soilT', 'perma']
KMS = [40, 80, 120, 145, 160, 220, 300, 400]
NY = 180
YRS = np.arange(2020, 2020 + NY + 1)


def job(args):
    sys.path.insert(0, HERE)
    import ch2_models as M
    tag, fbx, km, m = args
    cfg = dict(carbon='lenton', a2=0.003467, rartp=M.prs.DISC['std'], fbx=fbx, KM=float(km))
    if m:
        cfg['perma'] = dict(m=float(m), abrupt='smooth')
    try:
        d = M.run(cfg, nyears=NY)
    except ValueError:
        return tag, None
    keep = {k: d[k] for k in ('TATM', 'CO2ppm', 'release', 'land', 'ocean', 'atm')}
    return tag, keep


def summarise(d):
    T = d['TATM']; i = int(np.argmax(T))
    r = dict(peak_T=T[i], peak_year=int(YRS[i]), CO2_2020=d['CO2ppm'][0], CO2_peak=d['CO2ppm'].max(),
             CO2_min=d['CO2ppm'].min())
    for y in (2100, 2150, 2200):
        j = y - 2020
        r[f'T{y}'] = T[j]; r[f'CO2_{y}'] = d['CO2ppm'][j]; r[f'release_{y}'] = d['release'][j]
        r[f'dland_{y}'] = d['land'][j] - d['land'][0]; r[f'docean_{y}'] = d['ocean'][j] - d['ocean'][0]
    return r


def ts(tag, d, extra):
    return pd.DataFrame(dict(year=YRS, TATM=d['TATM'], CO2ppm=d['CO2ppm'], release=d['release'],
                             land_C=d['land'], ocean_C=d['ocean'], **{k: v for k, v in extra.items()}))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True); t0 = time.time()
    combos = list(itertools.product((0, 1), repeat=6))
    import ch2_models as M0
    jobs = [(('fac',) + c, {k: v for k, v in zip(FB[:5], c[:5])}, M0.KM, 1.0 if c[5] else 0) for c in combos]
    jobs += [(('fert', km, m), {}, km, m) for km in KMS for m in (0, 1.0, 4.0)]
    raw = os.path.join(OUT, "raw.pkl")
    if os.environ.get('REUSE') and os.path.exists(raw):
        import pickle; res = pickle.load(open(raw, 'rb'))
    else:
        with Pool(8) as pool:
            res = dict(pool.map(job, jobs, chunksize=1))
    print(f"runs done [{time.time() - t0:.0f}s]", flush=True)
    bad = [k for k, v in res.items() if v is None]
    print("FAILED (CO2 driven to zero):", bad, flush=True)
    import pickle; pickle.dump(res, open(os.path.join(OUT, "raw.pkl"), "wb"))
    if bad:
        sys.exit(0)

    # ---- factorial ----
    rows, tsl = [], []
    for c in combos:
        d = res[('fac',) + c]; r = dict(zip(FB, c)); r.update(summarise(d)); rows.append(r)
        tsl.append(ts(None, d, dict(zip(FB, c))))
    F = pd.DataFrame(rows); F.to_csv(os.path.join(OUT, "factorial_runs.csv"), index=False)
    TS = pd.concat(tsl); TS.to_csv(os.path.join(OUT, "factorial_timeseries.csv"), index=False)
    outs = ['peak_T', 'T2100', 'T2150', 'T2200', 'CO2_2100', 'CO2_2200']
    eff = []
    allon = F[(F[FB] == 1).all(axis=1)].iloc[0]
    for f in FB:                                   # main effects (2^k contrast) and contribution within the full model
        e = dict(term=f, kind='main effect (mean on - mean off)')
        e.update({o: F.loc[F[f] == 1, o].mean() - F.loc[F[f] == 0, o].mean() for o in outs}); eff.append(e)
        oth = [g for g in FB if g != f]
        one_off = F[(F[f] == 0) & (F[oth] == 1).all(axis=1)].iloc[0]
        e = dict(term=f, kind='contribution in full model (all on - this off)')
        e.update({o: allon[o] - one_off[o] for o in outs}); eff.append(e)
    for f, g in itertools.combinations(FB, 2):     # two-way interactions (2^k contrast)
        sgn = (2 * F[f] - 1) * (2 * F[g] - 1)
        e = dict(term=f'{f} x {g}', kind='interaction')
        e.update({o: (F[o] * sgn).sum() / (len(F) / 2) for o in outs}); eff.append(e)
    E = pd.DataFrame(eff); E.to_csv(os.path.join(OUT, "factorial_effects.csv"), index=False)
    # main effects over time
    me = []
    for f in FB:
        a = TS[TS[f] == 1].groupby('year')[['TATM', 'CO2ppm']].mean(); b = TS[TS[f] == 0].groupby('year')[['TATM', 'CO2ppm']].mean()
        me.append((a - b).assign(term=f).reset_index())
    pd.concat(me).to_csv(os.path.join(OUT, "factorial_main_effects_time.csv"), index=False)

    # ---- fertilisation sweep ----
    beta = pd.read_csv(os.path.join(HERE, "km_beta_map.csv")).set_index('KM')['land_beta']
    rows, tsl = [], []
    for km in KMS:
        b = float(np.interp(km, beta.index, beta.values))
        for m in (0, 1.0, 4.0):
            d = res[('fert', km, m)]; r = dict(KM=km, land_beta=b, perma_m=m); r.update(summarise(d)); rows.append(r)
            tsl.append(ts(None, d, dict(KM=km, land_beta=b, perma_m=m)))
    S = pd.DataFrame(rows)
    base = S[S.perma_m == 0].set_index('KM')
    for col in ('peak_T', 'T2100', 'T2150', 'T2200'):
        S['perma_d' + col] = S.apply(lambda r: r[col] - base.loc[r.KM, col], axis=1)
    S.to_csv(os.path.join(OUT, "fert_sweep_runs.csv"), index=False)
    pd.concat(tsl).to_csv(os.path.join(OUT, "fert_sweep_timeseries.csv"), index=False)
    pd.set_option('display.width', 220)
    print(E[E.kind != 'interaction'].round(3).to_string(index=False))
    print(E[E.kind == 'interaction'].assign(a=lambda x: x.T2200.abs()).sort_values('a', ascending=False).head(8).drop(columns='a').round(3).to_string(index=False))
    print(S[['KM', 'land_beta', 'perma_m', 'CO2_2020', 'peak_T', 'T2100', 'T2200', 'CO2_min', 'perma_dpeak_T', 'perma_dT2100', 'perma_dT2200']].round(2).to_string(index=False))
    print(f"[{time.time() - t0:.0f}s]")

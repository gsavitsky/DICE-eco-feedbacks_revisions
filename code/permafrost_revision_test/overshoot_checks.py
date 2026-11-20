"""
Overshoot / irreversibility checks (forward model only, 2020-2300).

 0. Hector validation: gradual permafrost module forced with Hector's own RCP temperatures (Woodard et al. 2021
    archive, default parameters, 2000-2100), frozen fraction compared with Hector's f_frozen.
 1. Zero-emissions commitment (ZEC): CO2 emissions prescribed (non-CO2 forcing held at 2020), emitted until
    cumulative emissions reach 1000 / 2000 PgC (DICE accounting, 633.5 PgC in 2020), then zero.
    Four models: DFAIR, DFAIR + permafrost, Lenton, Lenton + permafrost.
 2. Overshoot grid: reference path (net-zero industrial CO2 by 2050) and overshoot paths with peak +0.25/+0.5/+1.0 C
    above the reference peak, lasting ~50/100 years above it, built from delayed abatement (2025-2095) followed by net
    removal (MIU > 1). Each path is matched on cumulative net anthropogenic CO2 by 2300 *within each model*.
    Penalty = difference-in-differences: [X_os - X_ref]_with permafrost - [X_os - X_ref]_without.
 3. SCC of a 2080 tonne on every model x path, standard discounting, same pulse method throughout.
Permafrost variants: Q10 on sticky max (Hector) or running mean; abrupt thaw irreversible or reversible.
Outputs in overshoot_checks/.
"""
import os, sys, json, math, time
import numpy as np
import pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "overshoot_checks")
NY = 280
YRS = np.arange(2020, 2020 + NY + 1)
STD = dict(a2=0.003467)
VARIANTS = {'hector': dict(), 'runQ10': dict(q10='running'), 'revAbrupt': dict(abrupt_rev=True),
            'runQ10_revAbrupt': dict(q10='running', abrupt_rev=True)}
MODELS = {'DFAIR': dict(carbon='dfair'), 'Lenton': dict(carbon='lenton')}
OS_CASES = [(dp, D) for dp in (0.25, 0.5, 1.0) for D in (50, 100)]
_M = None


def M():
    global _M
    if _M is None:
        sys.path.insert(0, HERE)
        import ch2_models as m
        _M = m
    return _M


def cfg_of(model, variant):
    m = M()
    c = dict(MODELS[model], a2=STD['a2'], rartp=m.prs.DISC['std'])
    if variant is not None:
        c['perma'] = dict(m=1.0, abrupt='smooth', **VARIANTS[variant])
    return c


def base_controls():
    m = M()
    s = np.array([m.prs.exo(float(y))[11] for y in YRS])
    return s


def hump(t, a, b):
    x = np.clip((t - a) / (b - a), 0, 1)
    return np.sin(np.pi * x) ** 2


def miu_path(a, b, L):
    t = YRS.astype(float)
    ref = np.minimum(1.0, 0.05 + 0.95 * (t - 2020) / 30.0)
    return np.clip(ref - a * hump(t, 2025, 2095), 0, None) + b * hump(t, 2095, 2095 + L)


def run_path(cfg, miu, S, pulse=None):
    c = dict(cfg, ctrl=(miu, S))
    return M().run(c, pulse=pulse, nyears=NY)


def duration_above(T, thr):
    """years with T > thr, with linear interpolation at the crossings (continuous in the parameters)."""
    d = 0.0
    for i in range(len(T) - 1):
        a, b = T[i] - thr, T[i + 1] - thr
        if a > 0 and b > 0: d += 1
        elif a > 0 >= b: d += a / (a - b)
        elif b > 0 >= a: d += b / (b - a)
    return d


def match_b(cfg, a, L, S, ctot_ref, b0=0.3):
    """removal amplitude b so cumulative net CO2 (Ctot) in 2300 equals the reference (secant)."""
    f = lambda b: run_path(cfg, miu_path(a, b, L), S)['Ctot'][-1] - ctot_ref
    x0, x1 = 0.0, b0; f0, f1 = f(x0), f(x1)
    for _ in range(12):
        if abs(f1) < 0.05: break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0); x0, f0, x1 = x1, f1, x2; f1 = f(x1)
    return x1, f1


def tune_case(job):
    """tune (a, L) on DFAIR without permafrost so the overshoot path has peak Tref_pk + dp and duration D."""
    dp, D = job
    from scipy.optimize import least_squares
    cfg = cfg_of('DFAIR', None); S = base_controls()
    ref = run_path(cfg, miu_path(0, 0, 100), S); Tpk = ref['TATM'].max(); ctot = ref['Ctot'][-1]
    cache = {}

    def resid(p):
        a, L = float(p[0]), float(p[1])
        b, err = match_b(cfg, a, L, S, ctot)
        d = run_path(cfg, miu_path(a, b, L), S)
        cache[(a, L)] = b
        return [(d['TATM'].max() - Tpk - dp) / 0.05, (duration_above(d['TATM'], Tpk) - D) / 5.0]
    a0 = min(1.0, 0.45 * dp + 0.25); L0 = 40.0 + D
    r = least_squares(resid, [a0, L0], bounds=([0.0, 20.0], [1.0, 250.0]), diff_step=[0.02, 0.02], max_nfev=40)
    a, L = r.x; b, _ = match_b(cfg, a, L, S, ctot)
    d = run_path(cfg, miu_path(a, b, L), S)
    return dict(dp_target=dp, D_target=D, a=a, L=L, b=b, dp=d['TATM'].max() - Tpk, D=duration_above(d['TATM'], Tpk),
                Tref_pk=Tpk, miu_min=float(miu_path(a, b, L).min()), miu_max=float(miu_path(a, b, L).max()), ok=bool(r.success))


def run_combo(job):
    """one model x variant x path: retune b to match cumulative CO2 within this model, run base and 2080 pulse."""
    model, variant, case = job
    cfg = cfg_of(model, variant); S = base_controls()
    ref = run_path(cfg, miu_path(0, 0, 100), S)
    if case is None:
        d, b = ref, 0.0
    else:
        b, _ = match_b(cfg, case['a'], case['L'], S, ref['Ctot'][-1], b0=case['b'])
        d = run_path(cfg, miu_path(case['a'], b, case['L']), S)
    miu = miu_path(0, 0, 100) if case is None else miu_path(case['a'], b, case['L'])
    p = run_path(cfg, miu, S, pulse=(2080, 1.0))
    import permafrost_revision_sim as prs
    i = 2080 - 2020
    mu = [1000 * prs.scale1_ * d['cons'][j][1] * d['cons'][j][0] ** (-prs.eta_) for j in (i, i + 1)]
    scc = -1000 * (p['Uacc'][-1] - d['Uacc'][-1]) / (0.5 * (mu[0] + mu[1]))
    tag = 'ref' if case is None else f"os{case['dp_target']}_{case['D_target']}"
    ts = pd.DataFrame(dict(year=YRS, TATM=d['TATM'], CO2ppm=d['CO2ppm'], Ctot=d['Ctot'], release=d['release'],
                           Cperm=d['Cperm'], thawed=d['thawed'], miu=miu))
    ts['model'], ts['variant'], ts['path'] = model, variant or 'none', tag
    return dict(model=model, variant=variant or 'none', path=tag, b=b, Tpk=d['TATM'].max(), T2100=d['TATM'][80],
                T2200=d['TATM'][180], T2300=d['TATM'][-1], Ctot2300=d['Ctot'][-1], release2300=d['release'][-1],
                Cperm2300=d['Cperm'][-1], Uacc=d['Uacc'][-1], SCC2080=scc), ts


def zec_job(job):
    model, variant, total = job
    m = M()
    cfg = cfg_of(model, variant)
    x0 = m.init(cfg)
    F20 = sum(m.prs.exo(2020.0)[k] for k in (4, 5))                 # Fmisc + FGHGexo in 2020
    E20, g, cum0 = 11.0, 0.01, 633.5
    # emissions grow 1 %/yr from 11 GtC/yr until cumulative = total, then zero
    n = math.log(1 + g * (total - cum0) / E20) / math.log(1 + g); tz = 2020 + n
    Efn = lambda y: E20 * (1 + g) ** (y - 2020) if y < tz else 0.0
    c = dict(cfg, E_fn=Efn, F_nonco2=F20)
    d = m.run(c, nyears=NY)
    i0 = int(round(tz - 2020))
    T = d['TATM']
    return dict(model=model, variant=variant or 'none', total=total, t_zero=tz, T_zero=T[i0],
                ZEC25=T[i0 + 25] - T[i0], ZEC50=T[i0 + 50] - T[i0], ZEC100=T[i0 + 100] - T[i0],
                release_50=d['release'][i0 + 50] - d['release'][i0]), \
        pd.DataFrame(dict(year=YRS, TATM=T, CO2ppm=d['CO2ppm'], release=d['release'], model=model,
                          variant=variant or 'none', total=total))


def hector_validation():
    m = M(); phf = m.phf
    H = pd.read_csv(os.path.join(HERE, "hist_data", "hector_woodard_default_rcp.csv"))
    rows = []
    pm = dict(m=1.0, abrupt='none')
    for rcp, g in H.groupby('rcp'):
        g = g.sort_values('year'); yy = g.year.values.astype(float); TT = g.Tgav.values
        Tf = lambda y: float(np.interp(y, yy, TT)); dT = lambda y: float(np.interp(y, yy[1:] - 0.5, np.diff(TT)))
        f0 = g.f_frozen.values[0]
        S = np.array([f0 * 865.0, (1 - f0) * 865.0 * 0.5, 0.0, 0.0]); dt = 0.05
        out = [S[0] / 865.0]
        for k in range(int(100 / dt)):
            y = 2000 + k * dt
            fn = lambda S, y: m.perma_flows(Tf(y), dT(y), S[0], S[1], S[2], S[3], 0.0, pm)[0][:4]
            k1 = fn(S, y); k2 = fn(S + dt / 2 * k1, y + dt / 2); k3 = fn(S + dt / 2 * k2, y + dt / 2); k4 = fn(S + dt * k3, y + dt)
            S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            if (k + 1) % int(1 / dt) == 0: out.append(S[0] / 865.0)
        # Hector's frozen pool from its own annual rule (permafrost_c *= 1 - (f_old - f_new), simpleNbox-runtime.cpp)
        hp = [f0 * 865.0]
        for fo_, fn_ in zip(g.f_frozen.values[:-1], g.f_frozen.values[1:]):
            hp.append(hp[-1] * (1 - (fo_ - fn_)))
        from statistics import NormalDist
        for y, fh, fo, T, h in zip(g.year.values, g.f_frozen.values, out, TT, hp):
            rows.append(dict(rcp=rcp, year=y, Tgav=T, f_frozen_hector=fh,
                             f_diag_module=1 - NormalDist().cdf((math.log(phf.WF * T) - phf.MU) / phf.SIG_H),
                             pool_hector_rule=h, pool_module=fo * 865.0))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True); t0 = time.time()
    if os.environ.get('ZEC_ONLY') != '1':
        hv = hector_validation(); hv.to_csv(os.path.join(OUT, "hector_validation.csv"), index=False)
    else:
        hv = pd.read_csv(os.path.join(OUT, "hector_validation.csv"))
    e = hv.groupby('rcp').apply(lambda g: pd.Series(dict(ffrac_max_err=(g.f_diag_module - g.f_frozen_hector).abs().max(),
                                                         pool2100_hector=g.pool_hector_rule.iloc[-1], pool2100_module=g.pool_module.iloc[-1])))
    print("Hector validation\n", e.to_string(), flush=True)
    with Pool(8) as pool:
        zjobs = [(mo, v, tot) for tot in (1000, 2000) for mo in MODELS for v in (None, 'hector', 'revAbrupt')]
        zr = pool.map(zec_job, zjobs)
        pd.DataFrame([r[0] for r in zr]).to_csv(os.path.join(OUT, "zec_summary.csv"), index=False)
        pd.concat([r[1] for r in zr]).to_csv(os.path.join(OUT, "zec_timeseries.csv"), index=False)
        print("ZEC\n", pd.DataFrame([r[0] for r in zr]).round(3).to_string(), f"\n[{time.time()-t0:.0f}s]", flush=True)
        if os.environ.get('ZEC_ONLY') == '1':
            sys.exit(0)
        cases = pool.map(tune_case, OS_CASES)
        pd.DataFrame(cases).to_csv(os.path.join(OUT, "overshoot_paths.csv"), index=False)
        print("paths\n", pd.DataFrame(cases).round(3).to_string(), f"\n[{time.time()-t0:.0f}s]", flush=True)
        jobs = [(mo, v, c) for mo in MODELS for v in [None] + list(VARIANTS) for c in [None] + cases]
        res = pool.map(run_combo, jobs)
    R = pd.DataFrame([r[0] for r in res]); R.to_csv(os.path.join(OUT, "overshoot_runs.csv"), index=False)
    pd.concat([r[1] for r in res]).to_csv(os.path.join(OUT, "overshoot_timeseries.csv"), index=False)
    # difference-in-differences
    rows = []
    for mo in MODELS:
        nop = R[(R.model == mo) & (R.variant == 'none')].set_index('path')
        for v in VARIANTS:
            wp = R[(R.model == mo) & (R.variant == v)].set_index('path')
            for pth in wp.index:
                if pth == 'ref': continue
                row = dict(model=mo, variant=v, path=pth)
                for k in ('Tpk', 'T2100', 'T2200', 'T2300', 'release2300', 'SCC2080'):
                    row['DiD_' + k] = (wp.loc[pth, k] - wp.loc['ref', k]) - (nop.loc[pth, k] - nop.loc['ref', k])
                row['SCC2080_ref_withP'] = wp.loc['ref', 'SCC2080']; row['SCC2080_os_withP'] = wp.loc[pth, 'SCC2080']
                row['SCC2080_ref_noP'] = nop.loc['ref', 'SCC2080']; row['SCC2080_os_noP'] = nop.loc[pth, 'SCC2080']
                rows.append(row)
    DD = pd.DataFrame(rows); DD.to_csv(os.path.join(OUT, "overshoot_did.csv"), index=False)
    print("DiD\n", DD.round(3).to_string(), f"\n[{time.time()-t0:.0f}s]", flush=True)

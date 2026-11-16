"""
SSP emissions-driven check, 2020-2100 (2026-11-16).
Same CO2 emissions (RCMIP v5.1.0, MAGICC Fossil & Industrial + AFOLU) and same non-CO2 forcing
(RCMIP total ERF - CO2 ERF) into each carbon cycle, with the DICE 2023 two-box climate:
  * DFAIR (DICE 2023 carbon cycle, from the DICE 2020 state)
  * corrected Lenton (iii) and (iv), from their 2020 historical states, with and without the
    Hector-faithful permafrost module (730 PgC, gradual)
Reference: RCMIP SSP CO2 concentrations (MAGICC, emissions-consistent; Meinshausen et al. 2020).
Compared as the CO2 increase 2020->2050/2100 (ppm), so the 2020 baseline offset does not matter.
Output: ssp_check.csv
"""
import json, math, os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lenton_rebuild as LR
import lenton_full as LF
import permafrost_hector_faithful_sim as phf

P = LR.P
HD = os.path.join(HERE, "hist_data")
SSPS = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']


def series(fn, var, scen):
    d = pd.read_csv(os.path.join(HD, fn))
    r = d[(d.Region == 'World') & (d.Variable == var) & (d.Scenario == scen)]
    yrcols = [c for c in d.columns if c.isdigit()]
    s = r[yrcols].iloc[0].astype(float); s.index = s.index.astype(int)
    return s.interpolate()


def get(scen):
    ff = series('rcmip-emissions.csv', 'Emissions|CO2|MAGICC Fossil and Industrial', scen) / 3664.0
    lu = series('rcmip-emissions.csv', 'Emissions|CO2|MAGICC AFOLU', scen) / 3664.0
    erf = series('rcmip-radiative-forcing.csv', 'Effective Radiative Forcing', scen)
    eco2 = series('rcmip-radiative-forcing.csv', 'Effective Radiative Forcing|Anthropogenic|CO2', scen)
    conc = series('rcmip-concentrations.csv', 'Atmospheric Concentrations|CO2', scen)
    return ff, lu, (erf - eco2), conc


Teq1, Teq2, d1, d2, fco22x, Meq = P['Teq1'], P['Teq2'], P['d1'], P['d2'], P['fco22x'], P['Meq']
ESH = [0.2173, 0.224, 0.2824, 0.2763]
TAU = [P['tau0'], P['tau1'], P['tau2'], P['tau3']]


def alpha_of(cacc, T):
    target = 32.4 + P['irC'] * cacc + P['irT'] * T
    f = lambda a: sum(a * e * t * (1 - math.exp(-100 / (a * t))) for e, t in zip(ESH, TAU)) - target
    lo, hi = 0.1, 100.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)


def run_dfair(ff, lu, nonco2, dt=0.1):
    Tb = [0.1477, 1.099454]; R = [150.093, 102.698, 39.534, 6.1865]; Ctot = 633.5
    out = {}
    for y in range(2020, 2101):
        MAT = Meq + sum(R); out[y] = (MAT, Tb[0] + Tb[1])
        for j in range(int(1 / dt)):
            E = ff[y] + lu[y]
            MAT = Meq + sum(R); T = Tb[0] + Tb[1]
            a = alpha_of(Ctot - (MAT - Meq), T)
            F = fco22x * math.log(MAT / Meq) / math.log(2) + nonco2[y]
            Tb = [Tb[0] + dt * (Teq1 * F - Tb[0]) / d1, Tb[1] + dt * (Teq2 * F - Tb[1]) / d2]
            R = [r + dt * (e * E - r / (t * a)) for r, e, t in zip(R, ESH, TAU)]
            Ctot += dt * E
    return out


def run_lenton(var, ff, lu, nonco2, perma=False, dt=0.1):
    KM = json.load(open(os.path.join(HERE, "lenton_full.json")))[var]['KM']
    H0 = pd.read_csv(os.path.join(HERE, f"lenton_full_hist_{var}.csv"), index_col=0).loc[2020]
    keys = [f"{l}_{b}" for l in ('DIC', 'TA', 'NO3') for b in 'lhid']
    S = np.array([H0['Ca'], H0['Cv'], H0['Cs'], H0['Dcum']] + [H0[k] for k in keys] + [730.0 if perma else 0.0, 0.0, 0.0, 0.0])
    Tb = np.array([0.1477, 1.099454])
    cfg = dict(rate=phf.RATE, hector_ff=True, abrupt='none'); phf.CUR.clear(); phf.CUR.update(cfg)
    lag = int(200 / dt)
    buf = [phf.WF * phf.T_hist(1820 + j * dt) for j in range(lag)]; csum = sum(buf); t200m = 0.0; kk = 0
    out = {}

    def f(S, Tb, y, T200):
        Ca, Cv, Cs, D = S[:4]; X = S[4:16]; Cp, Th = S[16], S[17]
        T = Tb.sum()
        F = fco22x * math.log(Ca / Meq) / math.log(2) + nonco2[y]
        dTb = np.array([(Teq1 * F - Tb[0]) / d1, (Teq2 * F - Tb[1]) / d2])
        dX, asg, _, _ = LF.ocean_rhs(Ca, X, T if var == 'iv' else 0.0)
        Ph, Rp, Tu, Rs = LR.land(Ca, Cv, Cs, T, D, KM)
        if perma:
            dCp, dTh, _, _, _, Fm, _ = phf.perma_flows_cfg(T, dTb.sum(), Cp, Th, 0, 0, T200, 'none')
        else:
            dCp = dTh = Fm = 0.0
        dS = np.concatenate([[-Ph + Rp + Rs - asg + ff[y] + lu[y] + Fm, Ph - Rp - Tu - lu[y], Tu - Rs, lu[y]], dX, [dCp, dTh, 0, 0]])
        return dS, dTb

    for y in range(2020, 2101):
        out[y] = (S[0], Tb.sum())
        for j in range(int(1 / dt)):
            T200 = max(t200m, csum / lag); t200m = T200
            k1 = f(S, Tb, y, T200); k2 = f(S + dt / 2 * k1[0], Tb + dt / 2 * k1[1], y, T200)
            k3 = f(S + dt / 2 * k2[0], Tb + dt / 2 * k2[1], y, T200); k4 = f(S + dt * k3[0], Tb + dt * k3[1], y, T200)
            S = S + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); Tb = Tb + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            new = phf.WF * Tb.sum(); csum += new - buf[kk]; buf.append(new); kk += 1
    return out


if __name__ == "__main__":
    rows = []
    for scen in SSPS:
        ff, lu, nonco2, conc = get(scen)
        res = {'DFAIR (DICE 2023)': run_dfair(ff, lu, nonco2)}
        for var in ('iii', 'iv'):
            res[f'Lenton ({var})'] = run_lenton(var, ff, lu, nonco2)
            res[f'Lenton ({var}) + permafrost'] = run_lenton(var, ff, lu, nonco2, perma=True)
        ref = dict(model='RCMIP/MAGICC concentrations', scenario=scen, ppm2020=conc[2020],
                   d2050=conc[2050] - conc[2020], d2100=conc[2100] - conc[2020])
        rows.append(ref)
        for name, o in res.items():
            p = {y: o[y][0] * LR.PPM for y in o}
            rows.append(dict(model=name, scenario=scen, ppm2020=p[2020], d2050=p[2050] - p[2020], d2100=p[2100] - p[2020],
                             T2100=o[2100][1]))
        print(scen, " | ".join(f"{r['model']}: +{r['d2100']:.0f}" for r in rows if r['scenario'] == scen), flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "ssp_check.csv"), index=False)

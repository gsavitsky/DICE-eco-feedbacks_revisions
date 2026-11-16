"""
Full Lenton (2000) ocean with Knox & McElroy (1984, JGR 89 D3, Table 1) biology and chemistry (2026-11-16).

Tracers in each of 4 boxes (warm surface l, cold surface h, intermediate i, deep d): DIC, total
alkalinity (TA), nitrate (N). All are advected by Lenton's Fig. 1 circulation (Table 2 values):
  i->l (kT+kwi), l->i (kwi), l->h (kT), i->h (kU), h->d (kT+kU), d->i (kT+kU+kid), i->d (kid)
Biology (K&M Table 1):
  warm surface: all nitrate supplied is consumed (NO3 = 0 in WSW)
  cold surface: fraction alpha = 0.01 of incoming nitrate consumed
  C:N = 7 (organic), CaCO3:Corg = 0.25
  remineralisation: 30% of Corg (and N) in intermediate, 70% deep; CaCO3 dissolves 100% deep
  TA: CaCO3 formation removes 2 eq/mol, nitrate uptake adds 1 eq/mol (reversed on remineralisation)
Inventories: ocean-average TA 2365, DIC 2250, NO3 33 umol/kg (K&M Table 1); total ocean carbon set
  to Lenton Fig. 1 (37,740 GtC).
Air-sea (Lenton eq. 7) and carbonate chemistry as in lenton_rebuild.py (K0 Weiss 1974, K1/K2 Goyet &
  Poisson 1989, KB Dickson 1990, S = 34.7). T: warm 22 C, cold 2 C (+ dT if ocean warming, Lenton eq. 16).
Land: Lenton eqs. 1-6, 17 (as lenton_rebuild.py).
Steps: (1) spin up the ocean-atmosphere to steady state (no land, no forcing) and compare with Fig. 1;
       (2) land pools at equilibrium; (3) re-tune K_M to 1990 = 353 ppm under GCB emissions (Lenton 4.3),
       variants (iii) no ocean warming and (iv) ocean warming; report 2020 and sinks.
"""
import json, math, os, sys
import numpy as np
import pandas as pd
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lenton_rebuild as LR

P = LR.P; SY = LR.SY; RHO = LR.RHO; V = LR.V; A = LR.A
kT, kU, kwi, kid = LR.kT, LR.kU, LR.kwi, LR.kid
B = 'lhid'
MASS = {b: V[b] * RHO for b in B}                       # kg per box
ALPHA, CN, RAIN, F_ORG_I = 0.01, 7.0, 0.25, 0.3
Q = {('i', 'l'): kT + kwi, ('l', 'i'): kwi, ('l', 'h'): kT, ('i', 'h'): kU,
     ('h', 'd'): kT + kU, ('d', 'i'): kT + kU + kid, ('i', 'd'): kid}
IDX = {b: k for k, b in enumerate(B)}


def advect(c):
    """c: concentrations (mol/kg) per box -> d(total mol)/dt (mol/yr) per box."""
    d = np.zeros(4)
    for (src, dst), q in Q.items():
        f = q * SY * RHO * c[IDX[src]]
        d[IDX[src]] -= f; d[IDX[dst]] += f
    return d


def biology(N):
    """N: nitrate concentrations -> (dDIC, dTA, dN) in mol/yr, and export diagnostics."""
    Nin_l = Q[('i', 'l')] * SY * RHO * N[2]                                   # all consumed
    Nin_h = (Q[('l', 'h')] * SY * RHO * N[0] + Q[('i', 'h')] * SY * RHO * N[2])
    Ex = {'l': Nin_l, 'h': ALPHA * Nin_h}
    dD, dA, dN = np.zeros(4), np.zeros(4), np.zeros(4)
    for b, n in Ex.items():
        corg = CN * n; caco3 = RAIN * corg
        k = IDX[b]
        dD[k] -= corg + caco3; dA[k] += -2 * caco3 + n; dN[k] -= n
        dD[2] += F_ORG_I * corg; dN[2] += F_ORG_I * n; dA[2] -= F_ORG_I * n
        dD[3] += (1 - F_ORG_I) * corg + caco3; dN[3] += (1 - F_ORG_I) * n; dA[3] += 2 * caco3 - (1 - F_ORG_I) * n
    return dD, dA, dN, Ex


def ocean_rhs(Ca, X, dT=0.0):
    """X = [DIC(4), TA(4), N(4)] total mol. Returns dX (mol/yr), net air-sea (GtC/yr) and diagnostics."""
    Dm, Am, Nm = X[0:4], X[4:8], X[8:12]
    m = np.array([MASS[b] for b in B])
    cD, cA, cN = Dm / m, Am / m, Nm / m
    pa = Ca / P['ka'] * P['fgtm']
    fl = {}
    for b, T in (('l', LR.Tl0 + dT), ('h', LR.Th0 + dT)):
        k = IDX[b]
        ps = LR.pco2_surface(cD[k], cA[k], T)
        g = P['ks'] * A[b] * RHO * LR.K0(T)             # mol/yr/atm
        fl[b] = g * (pa - ps); fl[b + '_in'] = g * pa; fl[b + '_out'] = g * ps; fl[b + '_pco2'] = ps
    dD, dA, dN, Ex = biology(cN)
    dD = dD + advect(cD); dA = dA + advect(cA); dN = dN + advect(cN)
    dD[0] += fl['l']; dD[1] += fl['h']
    as_gtc = (fl['l'] + fl['h']) * 12e-15
    return np.concatenate([dD, dA, dN]), as_gtc, fl, Ex


def gtc(Dm):
    return Dm * 12e-15


# ---------------- coupled model ----------------
def rhs(S, T, eff, eluc, cfg):
    Ca, Cv, Cs, Dcum = S[0], S[1], S[2], S[3]
    X = S[4:]
    dX, as_gtc, fl, _ = ocean_rhs(Ca, X, T if cfg['ocean_warming'] else 0.0)
    Photo, PlantR, Turn, SoilR = LR.land(Ca, Cv, Cs, T, Dcum, cfg['KM'])
    dCa = -Photo + PlantR + SoilR - as_gtc + eff + eluc
    d = np.concatenate([[dCa, Photo - PlantR - Turn - eluc, Turn - SoilR, eluc], dX])
    return d, dict(ocean=as_gtc, nep=Photo - PlantR - SoilR)


def rk4(f, S, dt, *a):
    k1, _ = f(S, *a); k2, _ = f(S + dt / 2 * k1, *a); k3, _ = f(S + dt / 2 * k2, *a); k4, _ = f(S + dt * k3, *a)
    return S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def oa_rhs(Y, dT=0.0):  # ocean + atmosphere only (spin-up)
    dX, as_gtc, _, _ = ocean_rhs(Y[0], Y[1:], dT)
    return np.concatenate([[-as_gtc], dX]), None


if __name__ == "__main__":
    out = {}
    # ---- (1) spin-up ----
    m = np.array([MASS[b] for b in B])
    Mtot = m.sum()
    Ctot_mol = 37740.0 / 12e-15
    cD0 = Ctot_mol / Mtot
    # uniform initial concentrations at the K&M ocean averages; biology redistributes them (inventories conserved)
    def spin(TA_avg, years=3000, dt=0.25):
        Y = np.concatenate([[596.0], cD0 * m, TA_avg * m, 33e-6 * m])
        hs = []
        for yr in range(years):
            for _ in range(int(1 / dt)):
                Y = rk4(oa_rhs, Y, dt)
            if yr % 500 == 0:
                hs.append((yr, Y[0] * LR.PPM))
        return Y, hs

    Y, hist_spin = spin(2365e-6)
    print("With K&M average alkalinity 2365 umol/kg: steady state %.1f ppm" % (Y[0] * LR.PPM))
    # Lenton tuned the steady state to exactly 280 ppm (Sect. 3.1.3). His inventories are not published;
    # the least-constrained one is mean alkalinity, so it is tuned here (total carbon fixed at Fig. 1).
    TA_fit = brentq(lambda ta: spin(ta, years=1500)[0][0] * LR.PPM - 280.0, 2300e-6, 2365e-6, xtol=1e-9)
    Y, hist_spin = spin(TA_fit)
    print("Tuned mean alkalinity: %.1f umol/kg (K&M 2365)" % (TA_fit * 1e6))
    out['TA_avg_tuned_umolkg'] = TA_fit * 1e6
    dY, _ = oa_rhs(Y)
    _, as_gtc, fl, Ex = ocean_rhs(Y[0], Y[1:])
    pools = gtc(Y[1:5])
    conc = {lab: (Y[1 + 4 * j:5 + 4 * j] / m * 1e6).round(1).tolist() for j, lab in enumerate(['DIC', 'TA', 'NO3'])}
    print("Spin-up CO2 (ppm) every 500 yr:", [round(p, 1) for _, p in hist_spin])
    print("Steady state: Ca %.1f GtC = %.1f ppm (Lenton 596 / 280); max |dC/dt| %.2e GtC/yr"
          % (Y[0], Y[0] * LR.PPM, np.max(np.abs(gtc(dY[1:5])))))
    print("  ocean C pools (GtC) warm/cold/int/deep:", pools.round(0), " Lenton Fig.1: 730 140 10040 26830")
    print("  concentrations (umol/kg):", conc)
    print("  gross air-sea (GtC/yr): warm in %.1f out %.1f (57.0/57.9); cold in %.1f out %.1f (19.1/18.2)"
          % (fl['l_in'] * 12e-15, fl['l_out'] * 12e-15, fl['h_in'] * 12e-15, fl['h_out'] * 12e-15))
    print("  surface pCO2 (ppm): warm %.1f cold %.1f" % (fl['l_pco2'] * 1e6, fl['h_pco2'] * 1e6))
    print("  export (GtC/yr incl. CaCO3): warm %.2f cold %.3f" % (Ex['l'] * CN * (1 + RAIN) * 12e-15, Ex['h'] * CN * (1 + RAIN) * 12e-15))
    out['steady_state'] = dict(Ca=Y[0], ppm=Y[0] * LR.PPM, pools=pools.tolist(), conc=conc)
    Yss = Y.copy()

    # ---- (2)+(3) historical with K_M tuning ----
    E = pd.read_csv(os.path.join(HERE, "hist_data", "world_co2_1850_2023.csv"), index_col=0)

    def land_eq(KM, Ca):
        Photo, _, _, _ = LR.land(Ca, 1.0, 1.0, 0.0, 0.0, KM)
        rP = P['kr'] * P['kA'] * math.exp(-P['Ea'] / (P['R'] * 288))
        Cv = Photo / (rP + P['kt'])
        rS = P['ksr'] * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 - 227.13))))
        return Cv, P['kt'] * Cv / rS

    def hist(KM, ow, end=2020, keep=False, dt=0.25):
        cfg = dict(ocean_warming=ow, KM=KM)
        Cv, Cs = land_eq(KM, Yss[0])
        S = np.concatenate([[Yss[0], Cv, Cs, 0.0], Yss[1:]])
        rows = []
        for y in range(1850, end + 1):
            eff = float(E.loc[y, 'co2']) / 3664.0; eluc = float(E.loc[y, 'land_use_change_co2']) / 3664.0
            if keep:
                _, f = rhs(S, LR.T_hist(y), eff, eluc, cfg)
                rows.append(dict(year=y, ppm=S[0] * LR.PPM, ocean_sink=f['ocean'], nep=f['nep'], Ca=S[0], Cv=S[1],
                                 Cs=S[2], Dcum=S[3], **{f"{lab}_{b}": S[4 + 4 * j + IDX[b]]
                                                        for j, lab in enumerate(['DIC', 'TA', 'NO3']) for b in B}))
            if y == end:
                break
            for j in range(int(round(1 / dt))):
                S = rk4(rhs, S, dt, LR.T_hist(y + j * dt), eff, eluc, cfg)
        return pd.DataFrame(rows).set_index('year') if keep else S[0] * LR.PPM

    for tag, ow in (('iii', False), ('iv', True)):
        KM = brentq(lambda km: hist(km, ow, end=1990) - 353.0, 20, 800, xtol=0.1)
        H = hist(KM, ow, keep=True)
        H.to_csv(os.path.join(HERE, f"lenton_full_hist_{tag}.csv"))
        r = dict(KM=KM, ppm1990=H.loc[1990, 'ppm'], ppm2020=H.loc[2020, 'ppm'],
                 growth_1960_1990=(H.loc[1990, 'ppm'] - H.loc[1960, 'ppm']) / 30,
                 ocean_1980s=H.loc[1980:1989, 'ocean_sink'].mean(), nep_1980s=H.loc[1980:1989, 'nep'].mean(),
                 ocean_2010s=H.loc[2010:2019, 'ocean_sink'].mean(), nep_2010s=H.loc[2010:2019, 'nep'].mean())
        out[tag] = r
        print(f"\nvariant ({tag}) ocean warming={ow}: K_M = {KM:.1f}")
        print(f"   1990 {r['ppm1990']:.1f} | 2020 {r['ppm2020']:.1f} ppm (obs ~412) | growth 1960-90 {r['growth_1960_1990']:.2f} ppm/yr")
        print(f"   ocean sink 1980s {r['ocean_1980s']:.2f} (IPCC 2.0+-0.8) 2010s {r['ocean_2010s']:.2f} (GCB 2.6-2.9)"
              f" | NEP 1980s {r['nep_1980s']:.2f} 2010s {r['nep_2010s']:.2f} (GCB ~3.2)")
    json.dump(out, open(os.path.join(HERE, "lenton_full.json"), "w"), indent=1, default=float)
    np.save(os.path.join(HERE, "lenton_full_steadystate.npy"), Yss)

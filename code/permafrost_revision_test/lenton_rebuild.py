"""
Rebuild of the Lenton (2000, Tellus 52B) carbon cycle from the paper alone (2026-11-16).
Knox & McElroy (1984) not yet available: biological export and alkalinity are BACKED OUT from
Lenton's own published pre-industrial steady state (Fig. 1), so the model sits exactly on it.

Ocean (Fig. 1, Table 2), pools in GtC, concentrations C/V, water fluxes in m3/s:
  advection (thermohaline):  intermediate -> warm (kT), warm -> cold (kT), intermediate -> cold (kU),
                             cold -> deep (kT+kU), deep -> intermediate (kT+kU)
  two-way mixing:            warm <-> intermediate (kwi), intermediate <-> deep (kid)
  particles (fixed at their pre-industrial values, backed out of the steady state):
                             warm export Pw -> fraction fi to intermediate, (1-fi) to deep;
                             cold export Pc -> intermediate
  air-sea:  F = ks*A*rho*K0*(pCO2a - pCO2s)   (Lenton eq. 7)
  surface pCO2 from full carbonate chemistry (DIC + total alkalinity, K1/K2 Goyet & Poisson 1989,
  KB Dickson 1990, K0 Weiss 1974, S = 34.7). Surface alkalinities chosen so pCO2s gives the Fig. 1
  gross fluxes (warm 57.9 out / 57.0 in, cold 18.2 out / 19.1 in); held constant (fixed biology).
  Optional ocean warming (Lenton eq. 16): Tl = 295.15 + dT, Th = 275.15 + dT.
Land: Lenton eqs. 1-6, 17 (land use removed from vegetation; photosynthesis on Cv(s) = Cv0 - kD*D).
Validation: Lenton's historical run (Sect. 4): 1990 CO2 ~353 ppmv, 1980s ocean sink ~2.0 GtC/yr,
surface-deep exchange ~100 GtC/yr.
"""
import math, os, sys, json
import numpy as np
import pandas as pd
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    sys.path.insert(0, os.path.dirname(HERE))
    os.chdir(os.path.dirname(HERE))
    import dice_stella_validation as dsv
P = dict(dsv.P_)
SY = 3.15576e7
S_ = 34.7
RHO = P['rho']
GT2MOL = 1e15 / 12.0
PPM = 1 / P['ka'] * P['fgtm'] * 1e6            # ppm per GtC
V = dict(l=P['Vl'], h=P['Vh'], i=P['Vi'], d=P['Vd'])
kT, kU, kwi, kid = P['kT'], P['kU'], P['kwi'], P['kid']
A = dict(l=P['aO'] * P['Al'], h=P['aO'] * P['Ah'])
FIG1 = dict(Ca=596.0, l=730.0, h=140.0, i=10040.0, d=26830.0, Cv=550.0, Cs=1500.0)
GROSS = dict(l=(57.9, 57.0), h=(18.2, 19.1))     # (out, in) GtC/yr


# ---------------- carbonate chemistry ----------------
def K0(T):
    return math.exp(-60.2409 + 93.4517 * (100 / T) + 23.3585 * math.log(T / 100)
                    + S_ * (0.023517 - 0.023656 * (T / 100) + 0.0047036 * (T / 100) ** 2))

def K1K2(T):
    pK1 = 812.27 / T + 3.356 - 0.00171 * S_ * math.log(T) + 0.000091 * S_ ** 2
    pK2 = 1450.87 / T + 4.604 - 0.00385 * S_ * math.log(T) + 0.000182 * S_ ** 2
    return 10 ** -pK1, 10 ** -pK2

def KB(T):
    s = S_
    return math.exp((-8966.90 - 2890.53 * s ** .5 - 77.942 * s + 1.728 * s ** 1.5 - 0.0996 * s ** 2) / T
                    + 148.0248 + 137.1942 * s ** .5 + 1.62142 * s
                    - (24.4344 + 25.085 * s ** .5 + 0.2474 * s) * math.log(T) + 0.053105 * s ** .5 * T)

TB = 0.000416 * S_ / 35.0

def pco2_surface(DIC, ALK, T):
    """DIC, ALK in mol/kg -> pCO2 in atm (Park 1969 style: ALK = carbonate + borate alkalinity)."""
    k1, k2 = K1K2(T); kb = KB(T)
    f = lambda H: DIC * (k1 * H + 2 * k1 * k2) / (H * H + k1 * H + k1 * k2) + TB * kb / (kb + H) - ALK
    H = brentq(f, 1e-11, 1e-5, xtol=1e-16)
    return DIC * H * H / (H * H + k1 * H + k1 * k2) / K0(T)

def dic(Cgt, box):
    return Cgt * GT2MOL / (V[box] * RHO)


# ---------------- steady-state construction from Fig. 1 ----------------
Tl0, Th0 = P['Tl'], P['Th']
pa0 = FIG1['Ca'] / P['ka'] * P['fgtm']              # atm
pCO2s_target = {b: pa0 * GROSS[b][0] / GROSS[b][1] for b in ('l', 'h')}
ALK = {}
for b, T in (('l', Tl0), ('h', Th0)):
    D = dic(FIG1[b], b)
    ALK[b] = brentq(lambda a: pco2_surface(D, a, T) - pCO2s_target[b], 1.5e-3, 3.5e-3, xtol=1e-14)


def airsea(Ca, Cl, Ch, Tl, Th):
    pa = Ca / P['ka'] * P['fgtm']
    out = {}
    for b, C, T in (('l', Cl, Tl), ('h', Ch, Th)):
        ps = pco2_surface(dic(C, b), ALK[b], T)
        g = P['ks'] * A[b] * RHO * K0(T) * 12e-15        # GtC/yr per atm
        out[b] = g * (pa - ps)
        out[b + '_gross_in'], out[b + '_gross_out'] = g * pa, g * ps
    return out


def circulation(Cl, Ch, Ci, Cd):
    c = {b: x / V[b] for b, x in zip('lhid', (Cl, Ch, Ci, Cd))}
    F = {}
    F['i>l'] = (kT + kwi) * c['i'] * SY; F['l>i'] = kwi * c['l'] * SY
    F['l>h'] = kT * c['l'] * SY;          F['i>h'] = kU * c['i'] * SY
    F['h>d'] = (kT + kU) * c['h'] * SY;   F['d>i'] = (kT + kU + kid) * c['d'] * SY
    F['i>d'] = kid * c['i'] * SY
    dl = F['i>l'] - F['l>i'] - F['l>h']
    dh = F['l>h'] + F['i>h'] - F['h>d']
    di = F['l>i'] + F['d>i'] - F['i>l'] - F['i>h'] - F['i>d']
    dd = F['h>d'] + F['i>d'] - F['d>i']
    return np.array([dl, dh, di, dd]), F


# back out particle export so Fig. 1 is an exact steady state
as0 = airsea(FIG1['Ca'], FIG1['l'], FIG1['h'], Tl0, Th0)
circ0, F0 = circulation(FIG1['l'], FIG1['h'], FIG1['i'], FIG1['d'])
Pw = circ0[0] + as0['l']           # export needed from warm box
Pc = circ0[1] + as0['h']           # export needed from cold box
# intermediate:  circ_i + fi*Pw + Pc = 0  ->  fi
fi = -(circ0[2] + Pc) / Pw
BIO = dict(Pw=Pw, Pc=Pc, fi=fi)

# land: Lenton eqs 1-6; rate normalisations from the STELLA/Lenton parameter set
def land(Ca, Cv, Cs, T, Dcum, KM):
    ppm = Ca * PPM
    Cvs = P['Cv0'] - P['kD'] * Dcum if 'kD' in P else P['Cv0'] - 0.27 * Dcum
    # kMM normalises pre-industrial photosynthesis to Lenton's value (1.478 at KM=120, 1.578 at KM=145)
    kMM = 1.478 * (251.0 / (120 + 251.0)) / (251.0 / (KM + 251.0))
    Photo = P['kp'] * Cvs * kMM * ((ppm - P['kc']) / (KM + (ppm - P['kc']))) * ((15 + T) ** 2 * (25 - T) / 5625)
    PlantR = P['kr'] * Cv * P['kA'] * math.exp(-P['Ea'] / (P['R'] * (288 + T)))
    Turn = P['kt'] * Cv
    SoilR = P['ksr'] * Cs * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 + T - 227.13))))
    return Photo, PlantR, Turn, SoilR


def rhs(S, T, e_ff, e_luc, cfg):
    Ca, Cv, Cs, Cl, Ch, Ci, Cd, Dcum = S
    dT = T if cfg['ocean_warming'] else 0.0
    a = airsea(Ca, Cl, Ch, Tl0 + dT, Th0 + dT)
    circ, _ = circulation(Cl, Ch, Ci, Cd)
    Photo, PlantR, Turn, SoilR = land(Ca, Cv, Cs, T, Dcum, cfg['KM'])
    dCa = -Photo + PlantR + SoilR - a['l'] - a['h'] + e_ff + e_luc
    dCv = Photo - PlantR - Turn - e_luc
    dCs = Turn - SoilR
    dCl = circ[0] + a['l'] - BIO['Pw']
    dCh = circ[1] + a['h'] - BIO['Pc']
    dCi = circ[2] + BIO['fi'] * BIO['Pw'] + BIO['Pc']
    dCd = circ[3] + (1 - BIO['fi']) * BIO['Pw']
    return np.array([dCa, dCv, dCs, dCl, dCh, dCi, dCd, e_luc]), dict(
        ocean=a['l'] + a['h'], land=Photo - PlantR - SoilR - e_luc, nep=Photo - PlantR - SoilR)


def rk4(S, T, eff, eluc, cfg, dt):
    k1, _ = rhs(S, T, eff, eluc, cfg); k2, _ = rhs(S + dt / 2 * k1, T, eff, eluc, cfg)
    k3, _ = rhs(S + dt / 2 * k2, T, eff, eluc, cfg); k4, _ = rhs(S + dt * k3, T, eff, eluc, cfg)
    return S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def T_hist(y, kk=0.02):
    return 0.0 if y < 1850 else 1.25 * (math.exp(kk * (y - 1850)) - 1) / (math.exp(170 * kk) - 1)


if __name__ == "__main__":
    print("Pre-industrial (Fig. 1) checks")
    print("  surface pCO2 targets (ppm): warm %.1f, cold %.1f" % (pCO2s_target['l'] * 1e6, pCO2s_target['h'] * 1e6))
    print("  backed-out surface alkalinity (umol/kg): warm %.0f, cold %.0f; surface DIC: warm %.0f, cold %.0f"
          % (ALK['l'] * 1e6, ALK['h'] * 1e6, dic(730, 'l') * 1e6, dic(140, 'h') * 1e6))
    print("  model gross air-sea (GtC/yr): warm in %.1f out %.1f (Fig.1 57.0/57.9); cold in %.1f out %.1f (19.1/18.2)"
          % (as0['l_gross_in'], as0['l_gross_out'], as0['h_gross_in'], as0['h_gross_out']))
    print("  net air-sea: warm %.2f, cold %.2f GtC/yr" % (as0['l'], as0['h']))
    print("  water-borne fluxes (GtC/yr):", {k: round(v, 1) for k, v in F0.items()})
    surf_to_deep = F0['l>i'] + F0['h>d'] + BIO['Pw'] + BIO['Pc']
    print("  surface -> interior carbon transfer: %.1f GtC/yr (Lenton ~100); deep->surface upwelling i>l + i>h: %.1f"
          % (surf_to_deep, F0['i>l'] + F0['i>h']))
    print("  backed-out biological export: warm %.2f, cold %.2f GtC/yr; fraction of warm export remineralised in intermediate %.3f"
          % (BIO['Pw'], BIO['Pc'], BIO['fi']))
    S0 = np.array([FIG1['Ca'], FIG1['Cv'], FIG1['Cs'], FIG1['l'], FIG1['h'], FIG1['i'], FIG1['d'], 0.0])
    for cfg in (dict(name='Lenton (iii): no ocean warming, KM=120', ocean_warming=False, KM=120),
                dict(name='Lenton (iv): ocean warming, KM=145', ocean_warming=True, KM=145)):
        d0, _ = rhs(S0, 0.0, 0.0, 0.0, cfg)
        print("\n" + cfg['name'])
        print("  pre-industrial tendencies (GtC/yr):", np.round(d0[:7], 4))
        # 200-yr control run (no forcing)
        S = S0.copy()
        for _ in range(int(200 / 0.1)):
            S = rk4(S, 0.0, 0.0, 0.0, cfg, 0.1)
        print("  drift after 200 yr control: Ca %+.2f GtC, ocean pools %s" % (S[0] - S0[0], np.round(S[3:7] - S0[3:7], 2)))
        # historical run from steady state with GCB emissions (fossil + LUC) 1850-2020
        E = pd.read_csv(os.path.join(HERE, "hist_data", "world_co2_1850_2023.csv"), index_col=0)
        S = S0.copy(); rows = []; dt = 0.1
        for y in range(1850, 2021):
            eff = float(E.loc[y, 'co2']) / 3664.0
            eluc = float(E.loc[y, 'land_use_change_co2']) / 3664.0
            _, f = rhs(S, T_hist(y), eff, eluc, cfg)
            rows.append(dict(year=y, ppm=S[0] * PPM, ocean_sink=f['ocean'], land_sink_net=f['land'],
                             nep=f['nep'], eff=eff, eluc=eluc))
            for j in range(int(1 / dt)):
                S = rk4(S, T_hist(y + j * dt), eff, eluc, cfg, dt)
        H = pd.DataFrame(rows).set_index('year')
        dec80 = H.loc[1980:1989].mean()
        print("  1990 CO2: %.1f ppm (Lenton/obs 353); 2020: %.1f ppm (obs ~414)" % (H.loc[1990, 'ppm'], H.loc[2020, 'ppm']))
        print("  1980s mean: ocean sink %.2f GtC/yr (IPCC 2.0+-0.8), NEP %.2f, net land incl. LUC %.2f"
              % (dec80.ocean_sink, dec80.nep, dec80.land_sink_net))
        print("  2010s mean: ocean sink %.2f (GCB ~2.6-2.9), NEP %.2f (GCB land sink ~3.2)"
              % (H.loc[2010:2019].ocean_sink.mean(), H.loc[2010:2019].nep.mean()))
        tag = 'iv' if cfg['ocean_warming'] else 'iii'
        H.to_csv(os.path.join(HERE, f"lenton_rebuild_hist_{tag}.csv"))
        json.dump(dict(zip(['Ca', 'Cv', 'Cs', 'Cl', 'Ch', 'Ci', 'Cd', 'Dcum'], map(float, S))),
                  open(os.path.join(HERE, f"lenton_rebuild_state2021_{tag}.json"), "w"), indent=1)
    json.dump(dict(ALK_umolkg={k: v * 1e6 for k, v in ALK.items()}, BIO=BIO), open(os.path.join(HERE, "lenton_rebuild_backedout.json"), "w"), indent=1)

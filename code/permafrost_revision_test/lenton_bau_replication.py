"""
Replication of Lenton (2000) Sect. 5.1 'extended business as usual' with his own inputs (2026-11-16):
  * emissions: Table 3 (1800-1990, linear between 1800/1850/1947/1980/1990) and Table 4 IS92a
    (1990-2100, total and land-use), then linear decline to zero at 2200 (Table 5 scenario 1);
    land-use part follows Table 3/4 and is set to 0 after 2100
  * model version (iii): K_M = 120 ppmv, k_D = 0.27, land use removed from vegetation, no ocean warming;
    and (iv): K_M = 145, ocean surface warms with Ts (eq. 16)
  * Lenton's own grey-atmosphere climate (eqs. 10-15): Fd = (1-A) S/4 (1 + 3/4 tau),
    tau = 1.73 pCO2^0.263 + 0.0126 PH2O^0.503 + 0.0231, PH2O = H P0 exp(-L/(R Ts)),
    c dTs/dt = (Fd - sigma Ts^4) aE  (H tuned to give 288.15 K at 280 ppmv)
  * ocean/land = rebuilt Lenton with Knox & McElroy biology (lenton_full.py), started 1800 at steady state
Targets from Lenton's text: 1990 353 ppmv; 1990 NEP 1.98 & ocean 1.98 GtC/yr (iii); NEP peak ~2.5 GtC/yr
in 2030, NEP < 0 by 2120, max source ~1.1 GtC/yr ~2150, land loses 69 GtC while a source; ocean sink peak
~6.8 GtC/yr in 2130, 4.9 in 2060; CO2 peak ~985 ppmv in 2170; warming peak +5.5 K in 2180;
(iv): CO2 1020 ppmv, +5.7 K, ocean 4.1 GtC/yr in 2060.
"""
import math, os, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lenton_rebuild as LR
import lenton_full as LF

P = LR.P
S0, A0, SIG = 1368.0, 0.225, 5.67e-8
P0, Lh, Rg, c_heat, aE = 1.4e11, 43655.0, 8.314, 4.69e23, 5.101e14

def tau(ppm, Ts, H):
    return 1.73 * (ppm * 1e-6) ** 0.263 + 0.0126 * (H * P0 * math.exp(-Lh / (Rg * Ts))) ** 0.503 + 0.0231

def net_flux(ppm, Ts, H):
    return (1 - A0) * S0 / 4 * (1 + 0.75 * tau(ppm, Ts, H)) - SIG * Ts ** 4

H_fit = brentq(lambda h: net_flux(280.0, 288.15, h), 0.3, 0.9)

def interp(pts, y):
    xs, vs = zip(*pts)
    return float(np.interp(y, xs, vs))

TOT = [(1800, .35), (1850, .4), (1947, 2.0), (1980, 6.7), (1990, 7.4), (1995, 7.9), (2000, 8.4), (2005, 9.2), (2010, 9.9),
       (2015, 10.6), (2020, 11.4), (2025, 12.2), (2050, 14.5), (2075, 16.3), (2100, 20.3), (2200, 0.0)]
LUC = [(1800, .35), (1850, .35), (1947, .8), (1980, 1.5), (1990, 1.4), (2005, 1.36), (2010, 1.32), (2015, 1.28), (2020, 1.24),
       (2025, 1.2), (2050, .9), (2075, .45), (2100, 0.0), (2200, 0.0)]

def run(KM, ow, end=2400, dt=0.1):
    Yss = np.load(os.path.join(HERE, "lenton_full_steadystate.npy"))
    # land steady state at 280 ppm
    Photo, _, _, _ = LR.land(Yss[0], 1.0, 1.0, 0.0, 0.0, KM)
    rP = P['kr'] * P['kA'] * math.exp(-P['Ea'] / (P['R'] * 288))
    Cv = Photo / (rP + P['kt']); Cs = P['kt'] * Cv / (P['ksr'] * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 - 227.13)))))
    S = np.concatenate([[Yss[0], Cv, Cs, 0.0, 288.15], Yss[1:]])
    def f(S, y):
        Ca, Cv, Cs, D, Ts = S[:5]; X = S[5:]
        dT = Ts - 288.15; ppm = Ca * LR.PPM
        etot, eluc = interp(TOT, y), interp(LUC, y)
        dX, asg, _, _ = LF.ocean_rhs(Ca, X, dT if ow else 0.0)
        Ph, Rp, Tu, Rs = LR.land(Ca, Cv, Cs, dT, D, KM)
        dTs = net_flux(ppm, Ts, H_fit) * aE / c_heat * LR.SY
        return np.concatenate([[-Ph + Rp + Rs - asg + etot, Ph - Rp - Tu - eluc, Tu - Rs, eluc, dTs], dX]), (Ph - Rp - Rs, asg)
    rows = []
    for y in range(1800, end + 1):
        _, (nep, oc) = f(S, y)
        rows.append((y, S[0] * LR.PPM, S[4] - 288.15, nep, oc, S[1] + S[2]))
        for j in range(int(1 / dt)):
            t = y + j * dt
            k1 = f(S, t)[0]; k2 = f(S + dt / 2 * k1, t + dt / 2)[0]; k3 = f(S + dt / 2 * k2, t + dt / 2)[0]; k4 = f(S + dt * k3, t + dt)[0]
            S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return np.array(rows)

if __name__ == "__main__":
    print(f"grey-atmosphere relative humidity tuned to 288.15 K: H = {H_fit:.4f} (Lenton 0.6177)")
    for lab, KM, ow in (("(iii) K_M=120, no ocean warming", 120.0, False), ("(iv) K_M=145, ocean warming", 145.0, True)):
        R = run(KM, ow)
        yr, ppm, T, nep, oc, land = R.T
        i = lambda y: int(y - 1800)
        ip, it = int(np.argmax(ppm)), int(np.argmax(T))
        src = np.where((yr > 2050) & (nep < 0))[0]
        lost = (land[src[0]] - land[src[-1]]) if len(src) else 0.0
        print(f"\n{lab}")
        print(f"  1990: {ppm[i(1990)]:.1f} ppm (353) | NEP {nep[i(1990)]:.2f} (1.98) | ocean {oc[i(1990)]:.2f} (1.98) | dT {T[i(1990)]:.2f} K (1860-1990 ~0.6: {T[i(1990)] - T[i(1860)]:.2f})")
        print(f"  NEP peak {nep.max():.2f} in {int(yr[np.argmax(nep)])} (~2.5 in 2030) | NEP<0 from {int(yr[src[0]]) if len(src) else None} (~2120)"
              f" | max source {-nep.min():.2f} in {int(yr[np.argmin(nep)])} (~1.1, ~2150) | land C lost while source {lost:.0f} GtC (69)")
        print(f"  ocean sink 2060 {oc[i(2060)]:.2f} (iii 4.9 / iv 4.1) | peak {oc.max():.2f} in {int(yr[np.argmax(oc)])} (~6.8, 2130)")
        print(f"  CO2 peak {ppm[ip]:.0f} ppm in {int(yr[ip])} (iii 985/2170, iv 1020) | warming peak {T[it]:.2f} K in {int(yr[it])} (iii 5.5/2180, iv 5.7)")
        np.savetxt(os.path.join(HERE, f"lenton_bau_{'iii' if not ow else 'iv'}.csv"), R, delimiter=",",
                   header="year,ppm,dT,NEP,ocean_sink,land_C", comments="")

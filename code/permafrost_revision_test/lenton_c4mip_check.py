"""
C4MIP-style feedback check of the rebuilt Lenton carbon cycle (2026-11-16).
Concentration-driven 1%/yr CO2 from pre-industrial steady state (280 ppm) to 4x (140 yr), as in
Arora et al. (2020, Biogeosciences 17, 4173): beta (PgC/ppm) from the biogeochemically coupled run
(T = 0); gamma (PgC/K) from the difference with a fully coupled run, with warming prescribed as
T = TCR * log2(CO2/280), TCR = 1.8 K (CMIP6 mean ~1.8). Values at year 140 (4xCO2).
CMIP6 (Arora et al. 2020, Table A1 and abstract, 4xCO2, BGC-COU approach): land beta 0.97 +- 0.40, ocean beta
0.79 +- 0.07 PgC/ppm; land gamma -45.1 +- 50.6, ocean gamma -17.2 +- 5.0 PgC/K (11-model mean +- sample s.d.).
(Corrected 2026-11-20: earlier text quoted land beta 1.22, which is the 2xCO2 value, and land gamma -58.4.)
Runs for K_M = 208.6 (re-tuned) and 120 (Lenton's original).
"""
import math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lenton_rebuild as LR
import lenton_full as LF

Yss = np.load(os.path.join(HERE, "lenton_full_steadystate.npy"))
P = LR.P


def land_eq(KM, Ca):
    Photo, _, _, _ = LR.land(Ca, 1.0, 1.0, 0.0, 0.0, KM)
    rP = P['kr'] * P['kA'] * math.exp(-P['Ea'] / (P['R'] * 288))
    Cv = Photo / (rP + P['kt'])
    rS = P['ksr'] * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 - 227.13))))
    return Cv, P['kt'] * Cv / rS


def run(KM, coupled, ocean_warming=False, years=140, dt=0.1, TCR=1.8):
    Cv, Cs = land_eq(KM, Yss[0]); X = Yss[1:].copy()
    land0, ocean0 = Cv + Cs, LF.gtc(X[0:4]).sum()
    for k in range(int(years / dt)):
        t = k * dt
        def f(st, tt):
            ppm = 280 * 1.01 ** tt; Ca = ppm / LR.PPM
            T = TCR * math.log2(ppm / 280) if coupled else 0.0
            Cv_, Cs_ = st[0], st[1]
            Photo, PlantR, Turn, SoilR = LR.land(Ca, Cv_, Cs_, T, 0.0, KM)
            dX, _, _, _ = LF.ocean_rhs(Ca, st[2:], T if ocean_warming else 0.0)
            return np.concatenate([[Photo - PlantR - Turn, Turn - SoilR], dX])
        st = np.concatenate([[Cv, Cs], X])
        k1 = f(st, t); k2 = f(st + dt / 2 * k1, t + dt / 2); k3 = f(st + dt / 2 * k2, t + dt / 2); k4 = f(st + dt * k3, t + dt)
        st = st + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Cv, Cs, X = st[0], st[1], st[2:]
    dppm = 280 * 1.01 ** years - 280
    return (Cv + Cs - land0), (LF.gtc(X[0:4]).sum() - ocean0), dppm, 1.8 * math.log2(1.01 ** years)


if __name__ == "__main__":
    print("CMIP6 (Arora 2020, 4xCO2): land beta 0.97+-0.40, ocean beta 0.79+-0.07 PgC/ppm; "
          "land gamma -45+-51, ocean gamma -17+-5 PgC/K")
    for KM in (208.6, 120.0):
        dl_b, do_b, dppm, dT = run(KM, False)
        dl_c, do_c, _, _ = run(KM, True)
        print(f"K_M={KM:6.1f}: land beta {dl_b / dppm:.2f}  ocean beta {do_b / dppm:.2f} PgC/ppm | "
              f"land gamma {(dl_c - dl_b) / dT:.1f}  ocean gamma {(do_c - do_b) / dT:.1f} PgC/K "
              f"(no ocean warming in variant iii) | uptake by 4x: land {dl_c:.0f} ocean {do_c:.0f} PgC")

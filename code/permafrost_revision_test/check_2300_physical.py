"""Do all 128 Lenton factorial combinations (variant from ch2_models) stay physical to 2300 under the GAMS controls?"""
import os, sys, itertools
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PROC = ['ocean', 'oceanT', 'fert', 'photoT', 'vegT', 'soilT', 'perma']


def job(c):
    import ch2_models as M
    cfg = dict(carbon='lenton', a2=0.003467, rartp=M.prs.DISC['std'], fbx=dict(zip(PROC[:6], c[:6])), KM=float(M.KM))
    if c[6]: cfg['perma'] = dict(m=1.0, abrupt='smooth')
    try:
        d = M.run(cfg, nyears=280)
    except (ValueError, OverflowError):
        return c, None, None, None
    return c, float(d['CO2ppm'].min()), int(2020 + np.argmax(d['CO2ppm'] < 280)) if (d['CO2ppm'] < 280).any() else None, float(d['TATM'].min())


if __name__ == "__main__":
    with Pool(8) as p:
        R = p.map(job, list(itertools.product((0, 1), repeat=7)), chunksize=1)
    bad = [r for r in R if r[1] is None]
    ok = [r for r in R if r[1] is not None]
    print(f"crashed: {len(bad)} of 128", [dict(zip(PROC, r[0])) for r in bad][:6])
    below = [r for r in ok if r[2] is not None]
    print(f"below 280 ppm before 2300: {len(below)} of {len(ok)}; earliest year {min((r[2] for r in below), default=None)}")
    print(f"lowest CO2 {min(r[1] for r in ok):.0f} ppm; lowest T {min(r[3] for r in ok):.2f} C")
    allon = [r for r in ok if all(r[0])][0]; print(f"all on: min CO2 {allon[1]:.0f} ppm, below 280 from {allon[2]}, min T {allon[3]:.2f}")

"""
SSP emissions-driven check, 2020-2100, on the current Chapter 2 model (ch2_models: Lenton variant from LENTON_VARIANT,
default iv with K_M = 145; spun-up permafrost; 2026-11-20). Replaces ssp_check.py (which used the dropped K_M re-tune).
Same RCMIP v5.1.0 CO2 emissions (MAGICC Fossil & Industrial + AFOLU) and non-CO2 forcing (total ERF - CO2 ERF) into
DFAIR and Lenton, each with and without permafrost (Hector strength, smooth abrupt thaw), DICE 2023 two-box climate.
Reference: RCMIP SSP CO2 concentrations (MAGICC, emissions-consistent). Compared as the CO2 increase 2020->2050/2100.
Output: ssp_check_v2.csv
"""
import os, sys
import pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SSPS = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
CONFIGS = {'DFAIR (DICE 2023)': ('dfair', False), 'DFAIR + permafrost': ('dfair', True),
           'Lenton': ('lenton', False), 'Lenton + permafrost': ('lenton', True)}


def job(a):
    import ch2_models as M, ssp_check as O
    scen, name = a
    carbon, pf = CONFIGS[name]
    ff, lu, non, _ = O.get(scen)
    cfg = dict(carbon=carbon, a2=0.003467, rartp=M.prs.DISC['std'], E_fn=lambda y: (float(ff[int(y)]), float(lu[int(y)])),
               F_nonco2=lambda y: float(non[int(y)]))
    if pf: cfg['perma'] = dict(m=1.0, abrupt='smooth')
    d = M.run(cfg, nyears=80)
    c = d['CO2ppm']
    return dict(model=f"{name}" + (f" [Lenton {M.VARIANT}]" if carbon == 'lenton' else ''), scenario=scen, ppm2020=c[0],
                d2050=c[30] - c[0], d2100=c[80] - c[0], T2100=d['TATM'][80], release2100=d['release'][80])


if __name__ == "__main__":
    import ssp_check as O
    with Pool(8) as p:
        rows = p.map(job, [(s, n) for s in SSPS for n in CONFIGS], chunksize=1)
    for s in SSPS:
        conc = O.get(s)[3]
        rows.append(dict(model='RCMIP/MAGICC concentrations', scenario=s, ppm2020=conc[2020], d2050=conc[2050] - conc[2020],
                         d2100=conc[2100] - conc[2020]))
    D = pd.DataFrame(rows).sort_values(['scenario', 'model']); D.to_csv(os.path.join(HERE, "ssp_check_v2.csv"), index=False)
    print(D.pivot(index='model', columns='scenario', values='d2100').round(0).to_string())
    print(D[D.model.str.contains('permafrost')].pivot(index='model', columns='scenario', values='release2100').round(0).to_string())

"""
Gradual-thaw-only permafrost configuration (Reviewer 1, comment 1): Hector v2.3pf gradual thaw without the abrupt-thaw
term, alongside the full module (gradual + smooth Armstrong McKay abrupt thaw), on both carbon cycles.
GAMS controls, standard DICE settings, 2020-2200, SCC of a 1 GtCO2 pulse in 2020 (welfare to 2200).
Output: vis_corrected/ch2_gradual_only.csv
"""
import os, sys
import pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected/ch2_gradual_only.csv")
CASES = [(c, a) for c in ('dfair', 'lenton') for a in (None, 'none', 'smooth')]


def job(case):
    import ch2_models as M, ch2_scc_shapley as S
    carbon, ab = case
    cfg = dict(carbon=carbon, a2=0.003467, rartp=M.prs.DISC['std'])
    if ab: cfg['perma'] = dict(m=1.0, abrupt=ab)
    b = M.run(cfg, nyears=180); p = M.run(cfg, pulse=(2020, 1.0), nyears=180)
    T = b['TATM']
    return case, dict(peak_T=T.max(), T2100=T[80], T2150=T[130], T2200=T[180], release_2100=b['release'][80],
                      release_2200=b['release'][180], PgC_per_C_2100=b['release'][80] / T[80] if ab else float('nan'),
                      SCC2020=S.scc_of(M, b, p, 2020))


if __name__ == "__main__":
    with Pool(6) as p:
        res = dict(p.map(job, CASES, chunksize=1))
    rows = []
    for c, a in CASES:
        base = res[(c, None)]; r = res[(c, a)]
        lab = {'dfair': 'DFAIR', 'lenton': 'Lenton'}[c] + {None: '', 'none': ' + permafrost, gradual only', 'smooth': ' + permafrost, gradual + abrupt'}[a]
        rows.append(dict(configuration=lab, **r, dT_peak=r['peak_T'] - base['peak_T'], dT_2100=r['T2100'] - base['T2100'],
                         dT_2200=r['T2200'] - base['T2200'], SCC_pct=100 * (r['SCC2020'] / base['SCC2020'] - 1)))
    D = pd.DataFrame(rows); D.to_csv(OUT, index=False)
    print(D[['configuration', 'peak_T', 'T2200', 'release_2100', 'release_2200', 'PgC_per_C_2100', 'dT_peak', 'dT_2200', 'SCC2020', 'SCC_pct']].round(2).to_string(index=False))

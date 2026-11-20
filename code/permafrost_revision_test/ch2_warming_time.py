"""
Permafrost warming over time across the strength sweep (forward model, GAMS controls, 2020-2200).
For each carbon cycle x damage x strength m: warming added by permafrost in 2100/2150/2200 and at the peak, carbon
released, and the part due to carbon already thawed in 2020 ("committed": the same run with no new thaw after 2020,
i.e. gradual and abrupt thaw switched off, respiration of the 2020 thawed pool unchanged).
Writes vis_corrected/ch2_permafrost_warming_time.csv. Physical results do not depend on discounting (fixed controls).
"""
import os, sys
import numpy as np
import pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected")
MS = [0.1, 0.17, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
DAMAGE = {'standard': 0.003467, 'Howard-Sterner': 0.01}
YEARS = (2100, 2150, 2200)


def job(args):
    sys.path.insert(0, HERE)
    import ch2_models as M
    carbon, dmg, m, committed = args
    cfg = dict(carbon=carbon, a2=DAMAGE[dmg], rartp=M.prs.DISC['std'])
    if m is not None:
        cfg['perma'] = dict(m=m, abrupt='smooth')
        M.perma_2020_state(cfg['perma'])                       # spin-up with thaw on, cached before any switch-off
    f0, w0 = M.phf.dff_dT_hector, M.prs.w_onset
    try:
        if m is not None and committed:
            M.phf.dff_dT_hector = lambda T: 0.0; M.prs.w_onset = lambda T: 0.0
        d = M.run(cfg, nyears=180)
    finally:
        M.phf.dff_dT_hector, M.prs.w_onset = f0, w0
    return args, d['TATM'], d['release']


if __name__ == "__main__":
    jobs = [(c, dm, None, False) for c in ('dfair', 'lenton') for dm in DAMAGE]
    jobs += [(c, dm, m, k) for c in ('dfair', 'lenton') for dm in DAMAGE for m in MS for k in (False, True)]
    with Pool(8) as pool:                              # the committed switch is restored after every job
        res = {a: (T, R) for a, T, R in pool.map(job, jobs, chunksize=1)}
    rows = []
    for c in ('dfair', 'lenton'):
        for dm in DAMAGE:
            T0 = res[(c, dm, None, False)][0]
            for m in MS:
                T, R = res[(c, dm, m, False)]; Tc, Rc = res[(c, dm, m, True)]
                r = dict(carbon='DFAIR' if c == 'dfair' else 'Lenton', damage=dm, m=m,
                         PgC_per_C_2100=R[80] / T[80], dT_peak=T.max() - T0.max())
                for y in YEARS:
                    i = y - 2020
                    r[f'dT_{y}'] = T[i] - T0[i]; r[f'release_{y}'] = R[i]
                    r[f'dT_committed_{y}'] = Tc[i] - T0[i]; r[f'release_committed_{y}'] = Rc[i]
                    r[f'committed_share_{y}'] = (Tc[i] - T0[i]) / (T[i] - T0[i])
                rows.append(r)
    W = pd.DataFrame(rows); W.to_csv(os.path.join(OUT, "ch2_permafrost_warming_time.csv"), index=False)
    ts = []                                                    # time series for the timing figure (standard damage)
    for c in ('dfair', 'lenton'):
        T0 = res[(c, 'standard', None, False)][0]
        for m in MS:
            T, _ = res[(c, 'standard', m, False)]; Tc, _ = res[(c, 'standard', m, True)]
            ts.append(pd.DataFrame(dict(year=np.arange(2020, 2020 + len(T)), carbon='DFAIR' if c == 'dfair' else 'Lenton', m=m,
                                        dT_total=T - T0, dT_committed=Tc - T0)))
    pd.concat(ts).to_csv(os.path.join(OUT, "ch2_permafrost_warming_timeseries.csv"), index=False)
    print(W[W.damage == 'standard'][['carbon', 'm', 'PgC_per_C_2100', 'dT_peak', 'dT_2100', 'dT_2150', 'dT_2200',
                                     'committed_share_2100', 'committed_share_2200']].round(2).to_string(index=False))

"""
Shapley decomposition of the SCC for DFAIR (DICE 2023 carbon cycle): factors irC (uptake saturation), irT (temperature
dependence of uptake), perma (Hector permafrost). 'Off' = coefficient set to 0 (baseline IRF r0 and reservoirs kept).
2^3 = 8 runs, SCC of 1 GtCO2 pulses in 2020 and 2100, welfare to 2200, GAMS controls, standard DICE settings.
Output: vis_corrected/feedback_factorial/dfair_shapley_values.csv, dfair_shapley_interactions.csv
"""
import os, sys, itertools, math
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected/feedback_factorial")
F = ['irC', 'irT', 'perma']


def job(c):
    import ch2_models as M, ch2_scc_shapley as S
    std = dict(a2=0.003467, rartp=M.prs.DISC['std'])
    cfg = dict(std, carbon='dfair', **({'perma': dict(m=1.0, abrupt='smooth')} if c[2] else {}))
    return c, S.job((c, cfg, (float(c[0]), float(c[1]))))[1]


if __name__ == "__main__":
    with Pool(8) as p:
        res = dict(p.map(job, list(itertools.product((0, 1), repeat=3))))
    n = 3; rows, irows = [], []
    for mtr in ('SCC2020', 'SCC2100', 'peak_T', 'T2100', 'T2200'):
        v = {k: r[mtr] for k, r in res.items()}
        key = lambda on: tuple(1 if k in on else 0 for k in range(n))
        for i in range(n):
            oth = [k for k in range(n) if k != i]; phi = 0.0
            for s in range(n):
                for S in itertools.combinations(oth, s):
                    w = math.factorial(len(S)) * math.factorial(n - len(S) - 1) / math.factorial(n)
                    phi += w * (v[key(set(S) | {i})] - v[key(set(S))])
            rows.append(dict(metric=mtr, process=F[i], shapley=phi, all_off=v[(0, 0, 0)], all_on=v[(1, 1, 1)]))
        for i, j in itertools.combinations(range(n), 2):
            oth = [k for k in range(n) if k not in (i, j)]; tot = 0.0
            for s in range(len(oth) + 1):
                for S in itertools.combinations(oth, s):
                    w = math.factorial(len(S)) * math.factorial(n - len(S) - 2) / math.factorial(n - 1); S = set(S)
                    tot += w * (v[key(S | {i, j})] - v[key(S | {i})] - v[key(S | {j})] + v[key(S)])
            irows.append(dict(metric=mtr, pair=f"{F[i]} x {F[j]}", i=F[i], j=F[j], interaction=tot))
    V = pd.DataFrame(rows); V.to_csv(os.path.join(OUT, "dfair_shapley_values.csv"), index=False)
    I = pd.DataFrame(irows); I.to_csv(os.path.join(OUT, "dfair_shapley_interactions.csv"), index=False)
    print(V[V.metric.isin(['SCC2020', 'SCC2100'])].round(2).to_string(index=False))
    print(I[I.metric == 'SCC2020'].round(2).to_string(index=False))

"""
Pairwise Shapley interaction indices and off->on distributions for the SCC (and temperature) from the 128-run Lenton
factorial (shapley_runs.csv, ch2_scc_shapley.py; welfare to 2200). No model runs.
Shapley interaction index (Grabisch 1997) for processes i, j:
  I_ij = sum over S in N\\{i,j} of |S|! (n-|S|-2)! / (n-1)! * [v(S+i+j) - v(S+i) - v(S+j) + v(S)]
Positive: the two together raise the metric more than the sum of each alone (they reinforce each other).
Outputs (vis_corrected/feedback_factorial/): shapley_interactions.csv, shapley_off_to_on.csv
"""
import os, itertools, math
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected/feedback_factorial")
PROC = ['ocean', 'oceanT', 'fert', 'photoT', 'vegT', 'soilT', 'perma']
R = pd.read_csv(os.path.join(OUT, "shapley_runs.csv"))
n = len(PROC)
METRICS = ['SCC2020', 'SCC2100', 'peak_T', 'T2100', 'T2200']


def val(mtr):
    return {tuple(int(x) for x in r[PROC]): r[mtr] for _, r in R.iterrows()}


def key(on):
    return tuple(1 if k in on else 0 for k in range(n))


rows = []
for mtr in METRICS:
    v = val(mtr)
    for i, j in itertools.combinations(range(n), 2):
        oth = [k for k in range(n) if k not in (i, j)]
        tot = 0.0
        for s in range(len(oth) + 1):
            w = math.factorial(s) * math.factorial(n - s - 2) / math.factorial(n - 1)
            for S in itertools.combinations(oth, s):
                S = set(S)
                tot += w * (v[key(S | {i, j})] - v[key(S | {i})] - v[key(S | {j})] + v[key(S)])
        rows.append(dict(metric=mtr, pair=f"{PROC[i]} x {PROC[j]}", i=PROC[i], j=PROC[j], interaction=tot))
I = pd.DataFrame(rows)
# robustness: the 2x2 contrast in each of the 2^(n-2) settings of the other processes (equal-weight mean, range, value
# with all others on = the full model, all others off, and whether the sign is the same in every setting)
rob = []
for mtr in METRICS:
    v = val(mtr)
    for i, j in itertools.combinations(range(n), 2):
        oth = [k for k in range(n) if k not in (i, j)]; ds = []
        for bg in itertools.product((0, 1), repeat=len(oth)):
            S = {o for o, x in zip(oth, bg) if x}
            ds.append(v[key(S | {i, j})] - v[key(S | {i})] - v[key(S | {j})] + v[key(S)])
        ds = np.array(ds)
        rob.append(dict(metric=mtr, pair=f"{PROC[i]} x {PROC[j]}", equal_weight=ds.mean(), others_all_on=ds[-1],
                        others_all_off=ds[0], min=ds.min(), max=ds.max(), sign_consistent=bool((ds > 0).all() or (ds < 0).all())))
I = I.merge(pd.DataFrame(rob), on=['metric', 'pair'])
I.to_csv(os.path.join(OUT, "shapley_interactions.csv"), index=False)

rows = []
for mtr in METRICS:
    v = val(mtr)
    for i, p in enumerate(PROC):
        d = np.array([v[k[:i] + (1,) + k[i + 1:]] - v[k[:i] + (0,) + k[i + 1:]]
                      for k in v if k[i] == 0])
        sgn = np.sign(np.median(d))
        rows.append(dict(metric=mtr, process=p, mean=d.mean(), median=np.median(d), p5=np.percentile(d, 5),
                         p95=np.percentile(d, 95), min=d.min(), max=d.max(), same_sign=f"{int((np.sign(d) == sgn).sum())}/{len(d)}"))
D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, "shapley_off_to_on.csv"), index=False)

pd.set_option('display.width', 200)
for mtr in ('SCC2020', 'T2200'):
    print(f"--- {mtr}: largest pairwise interactions")
    x = I[I.metric == mtr].assign(a=lambda t: t.interaction.abs()).sort_values('a', ascending=False)
    print(x[['pair', 'interaction']].head(8).round(3).to_string(index=False))
    print(f"--- {mtr}: off->on change across the 64 settings of the other processes")
    print(D[D.metric == mtr][['process', 'mean', 'min', 'max', 'same_sign']].round(2).to_string(index=False))

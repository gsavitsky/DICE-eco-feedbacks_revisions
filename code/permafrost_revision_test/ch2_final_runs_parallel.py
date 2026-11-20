"""
Chapter 2 final runs, parallel (2026-11-17). Supersedes ch2_final_runs.py + ch2_sectionB_runs.py.
Simulations bounded to 2200 (CH2_NYEARS=180; welfare and SCC measured at 2200). Columns *_end refer to the last year.
Changes vs the 2026-11-16 runs: permafrost pools in 2020 come from a 1800-2020 spin-up from Hector's 865 PgC
(ch2_models.perma_2020_state); the feedback factorial freezes each dependence at its 2020 value (mass-consistent).

Outputs in DICE2023 Rproj/vis_corrected/:
  Welfare sensitivity/ (STELLA-layout CSVs, 16 runs), feedback sensitivity/ (64-run factorial),
  ch2_runs_summary.csv, ch2_permafrost_sweep.csv, ch2_permafrost_ceiling.csv, ch2_feedback_main_effects.csv,
  ch2_emissions_std.csv
"""
import itertools, os, shutil, sys, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RPROJ = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj")
OUT = os.path.join(RPROJ, "vis_corrected")
PERMA = dict(m=1.0, abrupt='smooth')
DAMAGE = {0: 0.003467, 1: 0.01}
SW = ['ocean', 'fert', 'photoT', 'soilT', 'vegT', 'perma']
COLS = {'ocean': 'Lenton Carbon.Ocean CO2 switch', 'fert': 'Lenton Carbon.Terr CO2 switch',
        'photoT': 'Lenton Carbon.photo temp switch', 'soilT': 'Lenton Carbon.temp-resp switch soil',
        'vegT': 'Lenton Carbon.temp-resp switch veg', 'perma': 'Permafrost dynamics.permafrost switch'}
KEEP = ('TATM', 'Uacc', 'release', 'thaw', 'cons', 'Y', 'Ygross', 'Cabate', 'Dfrac', 'MAT', 'CO2ppm', 'Ctot', 'thawed', 'Cperm')


def work(job):
    import ch2_models as M
    tag, cfg, pulse = job
    d = M.run(cfg, pulse=pulse, nyears=int(os.environ.get('CH2_NYEARS', '180')))   # simulations bounded to 2200
    return tag, {k: d[k] for k in KEEP}


def scc_from(b, p, year, dE=1.0):
    import permafrost_revision_sim as prs
    i = year - 2020
    mu = [1000 * prs.scale1_ * b['cons'][j][1] * b['cons'][j][0] ** (-prs.eta_) for j in (i, i + 1)]
    return -1000 * (p['Uacc'][-1] - b['Uacc'][-1]) / dE / (0.5 * (mu[0] + mu[1]))


if __name__ == "__main__":
    import permafrost_revision_sim as prs
    DISC = {0: prs.DISC['std'], 1: prs.DISC['zeroTP']}
    t0 = time.time()
    BASES = {1: ('DICE 2023', dict(carbon='dfair'), 0, 0), 5: ('DFAIR + permafrost', dict(carbon='dfair', perma=PERMA), 0, 1),
             9: ('DICE-C', dict(carbon='lenton'), 1, 0), 13: ('DICE-CP', dict(carbon='lenton', perma=PERMA), 1, 1)}
    jobs = []
    for b, (lab, c, _, _) in BASES.items():
        for hs in (0, 1):
            for zt in (0, 1):
                r = b + 2 * hs + zt; cfg = dict(c, a2=DAMAGE[hs], rartp=DISC[zt])
                jobs += [(('grid', r, None), cfg, None), (('grid', r, 2020), cfg, (2020, 1.0)), (('grid', r, 2050), cfg, (2050, 1.0))]
    MS = [0.1, 0.17, 0.25, 0.5, 1, 2, 4, 8, 16, 50]
    for carbon in ('dfair', 'lenton'):
        for hs in (0, 1):
            for zt in (0, 1):
                base = dict(carbon=carbon, a2=DAMAGE[hs], rartp=DISC[zt])
                jobs += [(('sw0', carbon, hs, zt, None), base, None), (('sw0', carbon, hs, zt, 2020), base, (2020, 1.0))]
                for m in MS:
                    cfg = dict(base, perma=dict(m=m, abrupt='smooth'))
                    jobs += [(('sw', carbon, hs, zt, m, None), cfg, None), (('sw', carbon, hs, zt, m, 2020), cfg, (2020, 1.0))]
    std = dict(a2=DAMAGE[0], rartp=DISC[0])
    for C0 in (740.0, 865.0, 991.0):          # Woodard et al. (2021) range for the pre-industrial pool
        for m in (1, 50, 500):
            jobs.append((('cap', C0, m), dict(std, carbon='lenton', perma=dict(m=m, abrupt='smooth', C0pre=C0)), None))
    RUN_FACTORIAL = os.environ.get('RUN_FACTORIAL') == '1'   # factorial postponed (some combinations unstable)
    for combo in (itertools.product([0, 1], repeat=6) if RUN_FACTORIAL else []):
        cfg = dict(std, carbon='lenton', fb=dict(zip(SW, combo)))
        if combo[5]:
            cfg['perma'] = PERMA
        jobs.append((('fact', combo), cfg, None))
    print(f"{len(jobs)} runs on {os.cpu_count()} cores", flush=True)
    with Pool(os.cpu_count()) as pool:
        R = dict(pool.map(work, jobs, chunksize=1))
    print(f"runs done in {time.time() - t0:.0f}s", flush=True)
    NY = int(os.environ.get('CH2_NYEARS', '180')); yrs = np.arange(2020, 2021 + NY)

    # ---- grid ----
    rows = []
    for b, (lab, c, lent, perm) in BASES.items():
        for hs in (0, 1):
            for zt in (0, 1):
                r = b + 2 * hs + zt; d = R[('grid', r, None)]; T = d['TATM']
                rows.append(dict(Run=r, config=lab, damage='Howard-Sterner' if hs else 'standard',
                                 discounting='zero time preference' if zt else 'standard', lenton=lent, permafrost=perm,
                                 peak_T=T.max(), peak_year=int(yrs[T.argmax()]), T2100=T[80], T2150=T[130], T2200=T[180],
                                 CO2ppm_2100=d['CO2ppm'][80], release_2100=d['release'][80], release_end=d['release'][-1],
                                 thawed_2020=d['thawed'][0], Cperm_2020=d['Cperm'][0], Uacc_end=d['Uacc'][-1],
                                 SCC2020=scc_from(d, R[('grid', r, 2020)], 2020), SCC2050=scc_from(d, R[('grid', r, 2050)], 2050)))
    Sm = pd.DataFrame(rows).sort_values('Run')
    Sm.to_csv(os.path.join(OUT, "ch2_runs_summary.csv"), index=False)
    print(Sm[['Run', 'config', 'damage', 'discounting', 'peak_T', 'peak_year', 'release_2100', 'Uacc_end', 'SCC2020']].round(3).to_string(index=False))
    WS = os.path.join(OUT, "Welfare sensitivity")
    os.makedirs(os.path.join(WS, "by_variable"), exist_ok=True)
    VARS = {'Climate_Atmospheric_temperature_anomaly_TATM_': 'TATM', 'Climate_atmospheric_carbon': 'MAT',
            'Economy_Accumulated_Utility_Uacc_': 'Uacc', 'Economy_Output_Y_': 'Y',
            'Economy_cost_of_emissions_abatement_Cabate_': 'Cabate',
            'Economy_damages_as_fraction_of_gross_output_Dfrac_': 'Dfrac', 'Economy_gross_output_Ygross_': 'Ygross'}
    for fname, key in VARS.items():
        df = pd.DataFrame({'Years': yrs})
        for r in range(1, 17):
            df[f'Run_{r}'] = R[('grid', r, None)][key]
        df.to_csv(os.path.join(WS, "by_variable", fname + ".csv"), index=False)
    pd.DataFrame([dict(Run=r, **{'Climate.Lenton?': BASES[r - (r - 1) % 4][2], 'Economy.altdam': ((r - 1) % 4) // 2,
                                 'Economy.altdisc': (r - 1) % 2, 'Permafrost dynamics.permafrost switch': BASES[r - (r - 1) % 4][3]})
                  for r in range(1, 17)]).to_csv(os.path.join(WS, "run_configurations.csv"), index=False)
    shutil.copy(os.path.join(RPROJ, "vis", "Welfare sensitivity", "GAMS_reference.csv"), os.path.join(WS, "GAMS_reference.csv"))
    E = pd.DataFrame({'Years': yrs[:-1]})
    for b, (lab, *_r) in BASES.items():
        E[lab] = np.diff(R[('grid', b, None)]['Ctot']) * 3.666
    E.to_csv(os.path.join(OUT, "ch2_emissions_std.csv"), index=False)

    # ---- sweep ----
    sw = []
    for carbon in ('dfair', 'lenton'):
        for hs in (0, 1):
            for zt in (0, 1):
                b0 = R[('sw0', carbon, hs, zt, None)]; s0 = scc_from(b0, R[('sw0', carbon, hs, zt, 2020)], 2020)
                for m in MS:
                    bb = R[('sw', carbon, hs, zt, m, None)]; s = scc_from(bb, R[('sw', carbon, hs, zt, m, 2020)], 2020); T = bb['TATM']
                    sw.append(dict(carbon='DFAIR' if carbon == 'dfair' else 'Lenton', damage='Howard-Sterner' if hs else 'standard',
                                   discounting='zero time preference' if zt else 'standard', m=m,
                                   PgC_per_C_2100=bb['release'][80] / T[80], release_2100=bb['release'][80], thaw_2100=bb['thaw'][80],
                                   release_end=bb['release'][-1], thawed_2020=bb['thawed'][0],
                                   dT_peak=T.max() - b0['TATM'].max(), dT_2100=T[80] - b0['TATM'][80],
                                   dT_2150=T[130] - b0['TATM'][130], dT_2200=T[180] - b0['TATM'][180], welfare_pct=100 * (bb['Uacc'][-1] / b0['Uacc'][-1] - 1),
                                   SCC=s, SCC_noperma=s0, SCC_pct=100 * (s / s0 - 1)))
    W = pd.DataFrame(sw); W.to_csv(os.path.join(OUT, "ch2_permafrost_sweep.csv"), index=False)
    print(W[(W.damage == 'standard')].pivot_table(index=['carbon', 'm'], columns='discounting',
                                                    values=['PgC_per_C_2100', 'SCC_pct', 'welfare_pct']).round(2).to_string())

    # ---- ceiling ----
    cap = []
    for C0 in (740.0, 865.0, 991.0):
        for m in (1, 50, 500):
            d = R[('cap', C0, m)]; T = d['TATM']
            cap.append(dict(C0_preindustrial=C0, m=m, Cperm_2020=d['Cperm'][0], thawed_2020=d['thawed'][0], thaw_2100=d['thaw'][80],
                            release_2100=d['release'][80], PgC_per_C_2100=d['release'][80] / T[80], T2100=T[80]))
    pd.DataFrame(cap).to_csv(os.path.join(OUT, "ch2_permafrost_ceiling.csv"), index=False)
    print(pd.DataFrame(cap).round(1).to_string(index=False))

    if not RUN_FACTORIAL:
        print(f'done {time.time() - t0:.0f}s (factorial skipped)'); sys.exit(0)
    # ---- factorial ----
    fs = os.path.join(OUT, "feedback sensitivity", "by_variable"); os.makedirs(fs, exist_ok=True)
    Tt = pd.DataFrame({'Years': yrs}); Cc = pd.DataFrame({'Years': yrs}); cfg_rows = []
    for r, combo in enumerate(itertools.product([0, 1], repeat=6), start=1):
        d = R[('fact', combo)]
        Tt[f'Run_{r}'] = d['TATM']; Cc[f'Run_{r}'] = d['MAT']
        cfg_rows.append({'Run': r, 'Climate.Lenton?': 1, **{COLS[k]: v for k, v in zip(SW, combo)}})
    Tt.to_csv(os.path.join(fs, "Climate_Atmospheric_temperature_anomaly_TATM_.csv"), index=False)
    Cc.to_csv(os.path.join(fs, "Lenton_Carbon_Atmospheric_Carbon.csv"), index=False)
    cdf = pd.DataFrame(cfg_rows); cdf.to_csv(os.path.join(OUT, "feedback sensitivity", "run_configurations.csv"), index=False)
    cdf['peak_T'] = [Tt[f'Run_{i}'][:281].max() for i in cdf.Run]; cdf['peak_C'] = [Cc[f'Run_{i}'][:281].max() for i in cdf.Run]
    eff = []
    for k in SW:
        on = cdf[cdf[COLS[k]] == 1].sort_values([COLS[j] for j in SW if j != k]); off = cdf[cdf[COLS[k]] == 0].sort_values([COLS[j] for j in SW if j != k])
        dT = on.peak_T.values - off.peak_T.values
        eff.append(dict(feedback=k, dT_peak_mean=dT.mean(), dT_peak_min=dT.min(), dT_peak_max=dT.max(),
                        dC_peak_mean=(on.peak_C.values - off.peak_C.values).mean()))
    pd.DataFrame(eff).to_csv(os.path.join(OUT, "ch2_feedback_main_effects.csv"), index=False)
    print(pd.DataFrame(eff).round(3).to_string(index=False))
    print(f"done {time.time() - t0:.0f}s")

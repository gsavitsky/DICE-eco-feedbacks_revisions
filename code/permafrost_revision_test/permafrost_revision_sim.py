"""
Permafrost-revision test for Chapter 2 (2026-11-16).  FORWARD SIMULATION ONLY, no optimisation.

Based on ../fixed_gams_controls_check.py (not modified): Chapter 2 models driven by the GAMS
controls (controls.json, linear interpolation as in STELLA), RK4 dt = 0.01, T200 from a fine
history buffer, 2020-2420.

The RHS below is a numpy copy of dice_pyomo_dicecp_reform_nb.dicecp_rhs (analytic gradual-thaw
derivative + refreeze limiter, tau = 2 yr), with three differences:
  * the frozen fraction and its derivative use the exact erf/normal-pdf form (STELLA NORMALCDF),
    not the tanh approximation the Pyomo model needs;
  * the old abrupt term is the STELLA hard step (IF TATM > 2 THEN 0.4*gradual), not the logistic
    surrogate used in Pyomo;
  * switches for damage (standard a2 = 0.003467 / Howard-Sterner a2 = 0.01) and discounting
    (standard rho = 0.001, beta_clim = 0.5 / zero-TP rho = beta_clim = 0, precautionary kept),
    exactly as the STELLA altdam / altdisc switches.
Three bookkeeping states are added: A (cumulative abrupt transfer Cperm -> thawedC), G (cumulative
net gradual transfer), R (cumulative release thawedC -> atmosphere).

Permafrost configurations:
  DICE-C                 : permafrost off
  DICE-CP old            : C0 = 1035, step abrupt (+40% of gradual rate when TATM > 2 C)  [Chapter 2]
  DICE-CP 865+step       : C0 = 865,  step abrupt                                         [effect of (a)]
  DICE-CP new-smooth     : C0 = 865,  smooth Armstrong McKay abrupt (1.5 C, 200 yr, +50%) [(a)+(b)]
  DICE-CP gradual-only   : C0 = 865,  no abrupt term                                      [(c)]
"""
import contextlib, io, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT)
os.chdir(PARENT)  # dsv reads controls.json etc. relative to its own folder; keep cwd safe
OUT = HERE
STELLA_XLSX = os.path.expanduser(
    "~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/STELLA results/FINALresults.xlsx")
DT = float(os.environ.get("PF_DT", "0.01"))
N = 400

print("Importing dice_stella_validation (runs its own checks on import, ~1 min)...", flush=True)
with contextlib.redirect_stdout(io.StringIO()):
    import dice_stella_validation as dsv
print("  done.", flush=True)
P = dsv.P_
miu_of, s_of = dsv.miu_of, dsv.s_of

# ---------------- constants (as in dice_pyomo_dicecp_reform_nb.py) ----------------
gamma_, a3_, theta2_ = 0.3, P['a3'], P['theta2']
fco22x_, Meq_ = P['fco22x'], P['Meq']
Teq1_, Teq2_, d1_, d2_ = P['Teq1'], P['Teq2'], P['d1'], P['d2']
Fc1_, Fc2_ = P['Fc1'], P['Fc2']
dK_, sd_cons_, eta_, pi_ = P['dK'], P['sd_cons'], P['eta'], P['pi_']
scale1_ = P['scale1']
Q10_, fstatic_ = P['Q10'], P['fstatic']
REFREEZE_TAU = 2.0


def _bcoef(T):
    return math.exp(-60.2409 + 93.4517 * (100 / T) + 23.3585 * math.log(T / 100)
                    + 34.7 * (0.023517 - 0.023656 * (T / 100) + 0.0047036 * (T / 100) ** 2))


bw_, bc_ = _bcoef(P['Tl']), _bcoef(P['Th'])


def exo_at(y):  # identical to dice_pyomo_dicecp_reform_nb.exo_at
    L = float(dsv.interp_annual('L', y)); aL = float(dsv.interp_annual('aL', y))
    sig = float(dsv.interp_annual('sigma', y)); eland = float(dsv.interp_annual('eland', y))
    Fmisc = float(dsv.interp_annual('F_Misc', y))
    CO2EGHGabateB = float(dsv.interp_annual('CO2E_GHGabateB', y))
    FGHGexo = float(dsv.fghgexo(y))
    pback = P['pback2050'] * math.exp(-0.01 * (y - 2050)) if y - 2020 <= 30 else P['pback2050'] * math.exp(-0.001 * (y - 2050))
    remiss = P['remiss0'] + ((P['remiss80'] - P['remiss0']) / 80) * (y - 2020) if y <= 2100 else P['remiss80']
    elandL = P['eland0_L'] * (1 - P['deland_L']) ** (y - 2015)
    return (L, aL, sig, eland, Fmisc, FGHGexo, CO2EGHGabateB, pback, remiss, elandL,
            float(miu_of(y)), float(s_of(y)))


_EXO = {}
def exo(y):
    k = round(y, 6)
    v = _EXO.get(k)
    if v is None:
        v = _EXO[k] = exo_at(y)
    return v


# ---------------- permafrost pieces ----------------
SQ2PI = math.sqrt(2 * math.pi)

def frozen_fraction(T):
    z = (math.log(1 + 2 * T) - 1.67) / 0.99
    return 1 - 0.5 * (1 + math.erf(z / math.sqrt(2)))

def dff_dT(T):
    z = (math.log(1 + 2 * T) - 1.67) / 0.99
    return -math.exp(-0.5 * z * z) / SQ2PI * (2.0 / 0.99) / (1 + 2 * T)

def smax0(x, eps=1e-6):
    return 0.5 * (x + math.sqrt(x * x + eps))

# Smooth Armstrong McKay (2022) abrupt thaw: parameters
AB_T0, AB_W, AB_TAU, AB_FRAC = 1.5, 0.5, 200.0, 0.5

def w_onset(T):
    return 0.5 * (1 + math.tanh((T - AB_T0) / AB_W))


def perma_flows(T, dT, Cperm, thawedC, A, G, T200, mode, Sst=0.0, inert=False):
    """Returns (dCperm, dthawed, dA, dG, dR, Fmelt, dSst). mode in {'step','smooth','none'}.
    inert=True (DIAGNOSTIC only): Hector v2.3pf Eq. 12/14 bookkeeping -- static carbon Sst accumulates as
    f_static x thaw, is never respired, refreezes in proportion; respiration acts on (thawedC - Sst)."""
    gradual_raw = -dff_dT(T) * dT * Cperm                       # Cperm -> thawedC, + when warming
    lim = -thawedC / REFREEZE_TAU                                # refreeze limiter (nb model)
    gradual = 0.5 * (gradual_raw + lim + math.sqrt((gradual_raw - lim) ** 2 + 1e-9))
    if mode == 'step':                                           # STELLA Eq. 17
        abrupt = 0.4 * max(0.0, gradual_raw) if T > 2.0 else 0.0
    elif mode == 'smooth':                                       # new, see REPORT.md
        target = AB_FRAC * w_onset(T) * smax0(G)
        abrupt = smax0(target - A) / AB_TAU
    else:
        abrupt = 0.0
    tot = gradual + abrupt
    if inert:
        J = (thawedC - Sst) * Q10_ ** (T200 / 10.0)
        dS = fstatic_ * tot if tot > 0 else (Sst / thawedC * tot if thawedC > 1e-12 else 0.0)
    else:
        J = (1 - fstatic_) * thawedC * Q10_ ** (T200 / 10.0)     # Eq. 18, thawedC -> atmosphere
        dS = 0.0
    return (-tot, tot - J, abrupt, gradual, J, J, dS)


# ---------------- full RHS ----------------
def rhs(S, y, T200, cfg):
    (Tb1, Tb2, K, Ctot, Uacc, FGHG, Ca, Cv, Cs, Col, Coh, Coi, Cod, Ldef,
     Cperm, thawedC, A, G, R, Sst) = S
    L, aL, sig, eland, Fmisc, FGHGexo, CO2EGHGabateB, pback, remiss, elandL, miu, s = exo(y)
    T = Tb1 + Tb2
    Ygross = aL * (L / 1000.0) ** (1 - gamma_) * K ** gamma_
    Dfrac = cfg['a2'] * T ** a3_
    theta1 = sig * remiss * pback / theta2_ / 1000.0
    Cabate = Ygross * theta1 * miu ** theta2_
    Y = Ygross * (1 - Dfrac) - Cabate
    I = s * Y; Cons = Y - I
    ECo2 = (sig * Ygross + eland) * (1 - miu)
    EIND = sig * Ygross * (1 - miu)
    F = fco22x_ * math.log(Ca / Meq_) / math.log(2) + Fmisc + FGHGexo
    dTb1 = (Teq1_ * F - Tb1) / d1_; dTb2 = (Teq2_ * F - Tb2) / d2_
    dFGHG = Fc1_ * CO2EGHGabateB * (1 - miu) + math.log(Fc2_) * FGHG
    dK = I - dK_ * K
    dCtot = ECo2 / 3.666
    Vc = min(sd_cons_ ** 2 * (y - 2020), sd_cons_ ** 2 * 235)
    rpre = -0.5 * Vc * eta_ ** 2
    STP = 1.0 / (1 + cfg['rartp']) ** (y - 2020) / (1 + rpre) ** (y - 2020)
    c = Cons * 1000.0 / L
    U = (c ** (1 - eta_) - 1) / (1 - eta_) - 1
    dUacc = U * L * STP * scale1_
    # Lenton carbon cycle
    Cvs = P['Cv0'] - Ldef * P['kD']
    pCO2 = Ca / P['ka'] * P['fgtm']; ppm = pCO2 * 1e6
    Photo = (P['kp'] * Cvs * P['kMM']) * ((ppm - P['kc']) / (P['KM'] + (ppm - P['kc']))) * ((15 + T) ** 2 * (25 - T) / 5625)
    PlantR = P['kr'] * Cv * P['kA'] * math.exp(-P['Ea'] / (P['R'] * (288 + T)))
    Turn = P['kt'] * Cv
    SoilR = P['ksr'] * Cs * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 + T - 227.13))))
    pCO2w = P['pCO2w0'] * (Col / P['Col0']) ** P['beta']
    pCO2c = P['pCO2c0'] * (Coh / P['Coh0']) ** P['beta']
    Faw = P['ks'] * P['aO'] * P['Al'] * bw_ * (pCO2 - pCO2w) * 12e-15 * P['rho']
    Fac = P['ks'] * P['aO'] * P['Ah'] * bc_ * (pCO2 - pCO2c) * 12e-15 * P['rho']
    k = 3.15576e7 * 12e-15
    WarmInt = P['kwi'] * (Col / P['Vl']) * k; IntDeep = P['kid'] * (Coi / P['Vi']) * k
    ColdInt = P['kU'] * (Coh / P['Vh']) * k; DeepInt = P['kid'] * (Cod / P['Vd']) * k
    IntWarm = (P['kT'] + P['kwi']) * (Coi / P['Vi']) * k
    GHGemis = EIND * 0.2727 + eland * (1 - miu) * 0.2727
    Defor = elandL / 2.12
    if cfg['perma']:
        dCp, dTh, dA, dG, dR, Fmelt, dSs = perma_flows(T, dTb1 + dTb2, Cperm, thawedC, A, G, T200, cfg['abrupt'],
                                                        Sst, cfg.get('inert', False))
    else:
        dCp = dTh = dA = dG = dR = Fmelt = dSs = 0.0
    dCa = PlantR + SoilR - Faw - Fac + GHGemis + Fmelt - Photo
    return np.array([dTb1, dTb2, dK, dCtot, dUacc, dFGHG, dCa, Photo - PlantR - Turn, Turn - SoilR,
                     IntWarm + Faw - WarmInt, Fac - ColdInt, WarmInt + ColdInt + DeepInt - IntDeep - IntWarm,
                     IntDeep - DeepInt, Defor, dCp, dTh, dA, dG, dR, dSs])


NAMES = ['Tbox1', 'Tbox2', 'K', 'Ctot', 'Uacc', 'FGHG', 'Ca', 'Cv', 'Cs', 'Col', 'Coh', 'Coi', 'Cod',
         'Ldef', 'Cperm', 'thawedC', 'A', 'G', 'R', 'Sstatic']


def t200_init(y):
    return (1.25 / (math.exp(0.01 * 200) - 1.0)) * (math.exp(0.01 * (y - 2020)) - 1.0)


def run(cfg, dt=DT):
    n_sub = int(round(1 / dt)); lag = int(round(200 / dt))
    S = np.array([0.1477, 1.099454, 295.0, 633.5, 0.0, 0.518, P['Ca0'], P['Cv0'], P['Cs0'], P['Col0'],
                  P['Coh0'], P['Coi0'], P['Cod0'], 0.0, cfg['C0'] if cfg['perma'] else 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    buf, out, minS = [], [], np.full(len(S), np.inf)
    for kk in range(N * n_sub + 1):
        y = 2020 + kk * dt
        if kk % n_sub == 0:
            out.append(S.copy())
        buf.append(S[0] + S[1])
        if kk == N * n_sub:
            break
        T200 = buf[kk - lag] if kk >= lag else t200_init(y)
        k1 = rhs(S, y, T200, cfg); k2 = rhs(S + dt / 2 * k1, y + dt / 2, T200, cfg)
        k3 = rhs(S + dt / 2 * k2, y + dt / 2, T200, cfg); k4 = rhs(S + dt * k3, y + dt, T200, cfg)
        S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        minS = np.minimum(minS, S)
    A = np.array(out)
    d = {nm: A[:, i] for i, nm in enumerate(NAMES)}
    d['TATM'] = d['Tbox1'] + d['Tbox2']
    d['min_thawedC'] = float(minS[15]); d['min_Cperm'] = float(minS[14])
    return d


# ---------------- configurations ----------------
RARTP_STD = math.exp(0.001 + 0.5 * pi_) - 1
DAMAGE = {'std': 0.003467, 'HS': 0.01}
DISC = {'std': RARTP_STD, 'zeroTP': 0.0}
MODELS = {
    'DICE-C':               dict(perma=False, C0=0.0, abrupt='none'),
    'DICE-CP old':          dict(perma=True, C0=1035.0, abrupt='step'),
    'DICE-CP 865+step':     dict(perma=True, C0=865.0, abrupt='step'),
    'DICE-CP new-smooth':   dict(perma=True, C0=865.0, abrupt='smooth'),
    'DICE-CP gradual-only': dict(perma=True, C0=865.0, abrupt='none'),
    # DIAGNOSTIC (not requested): gradual-only with Hector's inert static fraction
    'DIAG gradual-only+inert-static': dict(perma=True, C0=865.0, abrupt='none', inert=True),
}


def stella_xlsx():
    import openpyxl
    R = list(openpyxl.load_workbook(STELLA_XLSX, read_only=True, data_only=True)["Results"].iter_rows(values_only=True))
    h = R[0]
    def col(run, key):
        k = [i for i, x in enumerate(h) if x and str(x).startswith(f"Run {run}:") and key in str(x)][0]
        return np.array([float(r[k]) for r in R[1:] if r[0] is not None])
    out = {}
    for run in [1, 2, 3, 4, 9, 10, 11, 12, 13, 14, 15, 16]:
        out[run] = dict(TATM=col(run, "Atmospheric temperature"), Uacc=col(run, "Uacc"), Ca=col(run, "atmospheric carbon"))
    return out


# STELLA run numbers for each (model, damage, disc)
def stella_run(base, dam, disc):
    return base + (0 if dam == 'std' else 2) + (0 if disc == 'std' else 1)


if __name__ == "__main__":
    t0 = time.time()
    ST = stella_xlsx()
    results = {}
    cells = [(dam, disc) for disc in ('std', 'zeroTP') for dam in ('std', 'HS')]
    for name, m in MODELS.items():
        for dam, disc in cells:
            cfg = dict(m, a2=DAMAGE[dam], rartp=DISC[disc])
            ts = time.time()
            results[(name, dam, disc)] = run(cfg)
            print(f"  {name:22s} damage={dam:3s} disc={disc:6s}  U_acc(2420)={results[(name, dam, disc)]['Uacc'][-1]:14.2f}"
                  f"  ({time.time() - ts:.0f}s)", flush=True)
    yrs = np.arange(2020, 2421)

    # ---- time series CSV (standard damage, standard discounting; plus Uacc for all cells) ----
    import csv
    with open(os.path.join(OUT, "timeseries_std.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["model", "year", "TATM", "Ca", "Cperm", "thawedC", "abrupt_cum_A", "gradual_cum_G",
                     "release_cum_R", "Uacc"])
        for name in MODELS:
            d = results[(name, 'std', 'std')]
            for i, y in enumerate(yrs):
                wr.writerow([name, y] + [f"{d[k][i]:.6f}" for k in ("TATM", "Ca", "Cperm", "thawedC", "A", "G", "R", "Uacc")])
    with open(os.path.join(OUT, "uacc_timeseries_all_cells.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["model", "damage", "discounting", "year", "TATM", "Uacc"])
        for (name, dam, disc), d in results.items():
            for i, y in enumerate(yrs):
                wr.writerow([name, dam, disc, y, f"{d['TATM'][i]:.6f}", f"{d['Uacc'][i]:.6f}"])

    # ---- reproduction check vs STELLA ----
    print("\nReproduction vs STELLA (FINALresults.xlsx):")
    repro = []
    for name, base in (('DICE-C', 9), ('DICE-CP old', 13)):
        for dam, disc in cells:
            r = stella_run(base, dam, disc); d = results[(name, dam, disc)]
            dT = float(np.max(np.abs(d['TATM'] - ST[r]['TATM'])))
            dCa = float(np.max(np.abs(d['Ca'] - ST[r]['Ca'])))
            row = dict(model=name, damage=dam, disc=disc, stella_run=r, uacc_py=d['Uacc'][-1],
                       uacc_stella=ST[r]['Uacc'][-1], uacc_diff=d['Uacc'][-1] - ST[r]['Uacc'][-1],
                       max_abs_dTATM=dT, max_abs_dCa=dCa)
            repro.append(row)
            print(f"  {name:12s} {dam:3s} {disc:6s} run {r:2d}: U_acc py {row['uacc_py']:14.2f} STELLA {row['uacc_stella']:14.2f}"
                  f" diff {row['uacc_diff']:+8.2f}  max|dTATM| {dT:.4f} C  max|dCa| {dCa:.2f} GtC")
    with open(os.path.join(OUT, "reproduction_check.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(repro[0])); wr.writeheader(); wr.writerows(repro)

    # ---- summary CSV ----
    def tstats(T):
        i = int(np.argmax(T))
        above = lambda th: int(np.sum(T > th))
        def first(th):
            j = np.where(T > th)[0]; return int(yrs[j[0]]) if len(j) else None
        def last(th):
            j = np.where(T > th)[0]; return int(yrs[j[-1]]) if len(j) else None
        return dict(peak_T=T[i], peak_year=int(yrs[i]), T2100=T[80], T2200=T[180], T2300=T[280], T2420=T[400],
                    first_above_1p5=first(1.5), first_above_2=first(2.0), first_above_2p5=first(2.5), first_above_3=first(3.0),
                    years_above_1p5=above(1.5), years_above_2=above(2.0), years_above_2p5=above(2.5),
                    last_above_2=last(2.0), last_above_1p5=last(1.5))
    rows = []
    base_std = {'std': ST[1]['Uacc'][-1], 'zeroTP': ST[2]['Uacc'][-1]}   # DICE 2023 standard damage (STELLA runs 1, 2)
    # DICE 2023 reference rows (STELLA; unaffected by the permafrost module)
    for dam, disc in cells:
        r = stella_run(1, dam, disc)
        row = dict(model='DICE 2023 (STELLA run %d)' % r, damage=dam, disc=disc, source='STELLA FINALresults.xlsx')
        row.update(tstats(ST[r]['TATM']))
        row.update(uacc_2420=ST[r]['Uacc'][-1], loss_pct_vs_DICE2023_std=100 * (ST[r]['Uacc'][-1] / base_std[disc] - 1))
        rows.append(row)
    for (name, dam, disc), d in results.items():
        row = dict(model=name, damage=dam, disc=disc, source='python forward sim, GAMS controls')
        row.update(tstats(d['TATM']))
        row.update(uacc_2420=d['Uacc'][-1], loss_pct_vs_DICE2023_std=100 * (d['Uacc'][-1] / base_std[disc] - 1))
        if MODELS[name]['perma']:
            dc = results[('DICE-C', dam, disc)]
            Rel2100 = d['R'][80]; T2100 = d['TATM'][80]
            row.update(C0=MODELS[name]['C0'], thaw_2100=d['G'][80] + d['A'][80], abrupt_2100=d['A'][80],
                       release_2100=Rel2100, release_2200=d['R'][180], release_2300=d['R'][280], release_2420=d['R'][400],
                       Cperm_loss_2420=MODELS[name]['C0'] - d['Cperm'][400], abrupt_2420=d['A'][400],
                       thawedC_peak=float(np.max(d['thawedC'])), thawedC_peak_year=int(yrs[np.argmax(d['thawedC'])]),
                       PgC_per_C_2100_vs_preind=Rel2100 / T2100,
                       PgC_per_C_2100_vs_2020=Rel2100 / (T2100 - d['TATM'][0]),
                       Ca_excess_over_DICEC_2100=d['Ca'][80] - dc['Ca'][80],
                       Ca_excess_per_C_2100=(d['Ca'][80] - dc['Ca'][80]) / T2100,
                       dT_vs_DICEC_2100=T2100 - dc['TATM'][80],
                       min_thawedC=d['min_thawedC'], min_Cperm=d['min_Cperm'],
                       mass_err_max=float(np.max(np.abs(d['Cperm'] + d['thawedC'] + d['R'] - MODELS[name]['C0']))))
        rows.append(row)
    keys = []
    for r_ in rows:
        for k_ in r_:
            if k_ not in keys: keys.append(k_)
    with open(os.path.join(OUT, "summary_all_cells.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        for r_ in rows:
            wr.writerow({k_: (f"{v:.10g}" if isinstance(v, float) else v) for k_, v in r_.items()})
    json.dump({f"{a}|{b}|{c}": {k_: (v.tolist() if isinstance(v, np.ndarray) else v) for k_, v in d.items()}
               for (a, b, c), d in results.items()}, open(os.path.join(OUT, "results_all.json"), "w"))
    print(f"\nWrote CSVs to {OUT}  (total {time.time() - t0:.0f}s)")

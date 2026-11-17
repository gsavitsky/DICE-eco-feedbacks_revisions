"""
Chapter 3 re-optimisation with the corrected Chapter 2 physics (2026-11-16).

DICE_CONFIG = 'dicec'  : DICE 2023 economy + Lenton (2000) carbon cycle
DICE_CONFIG = 'dicecp' : + Hector v2.3pf-faithful permafrost with smooth Armstrong McKay abrupt thaw
(DICE 2023 itself is unchanged; its existing solution stands.)

Physics = permafrost_revision_test/ch2_models.py ('lenton' carbon cycle, K_M = 120, variant iii):
  * land: Lenton eqs. 1-6 & 17; land-use CO2 removed from vegetation, C_v(s) = C_v0 - k_D * Dcum
  * ocean: Lenton Fig. 1 four-box circulation in GtC (fluxes = Q * C/V, no mol conversion), Knox & McElroy
    biology. With fixed circulation and biology, alkalinity and nitrate are exactly steady (verified over
    1850-2020), so only the 4 DIC pools are states and biological export is a constant source/sink per box.
  * surface pCO2: exact carbonate chemistry tabulated as a function of each surface box's DIC (TA, T, S fixed)
    and fitted with a polynomial in log space (fit error reported at start-up; DIC states bounded to the fit range)
  * permafrost (dicecp): 0.02/yr respiration, 730 PgC in 2020, ln(wf*T) frozen fraction (tanh normal-CDF
    approximation, as in the previous Pyomo models), refreeze limiter, T200 = 200-yr running mean of wf*T
    (pre-2020 values from the same historical reconstruction as ch2_models.py). Hector's 'sticky' Q10 is
    omitted: tested in ch2_models.py, it changes nothing on these paths (U_acc identical to 1e-3).
  * initial state: 2020 Lenton state from the 1800-2020 historical run (lenton_km120_state2020.json); permafrost
    pools from the 1800-2020 permafrost spin-up (perma_state2020.json)
Numerics: annual decision variables MIU/S, annual RK4 defect constraints with SUBSTEPS per year (env), as before.

Env: DICE_CONFIG, DICE_SUBSTEPS (default 1), DICE_WARM_START (controls.json), DICE_OUTPUT_FILE,
     DICE_TRAJ_FILE, DICE_HORIZON_PAD (50), DICE_LINEAR_SOLVER (ma27), DICE_HSLLIB, DICE_CHECK_ONLY (1 = no solve)
"""
import contextlib, io, json, math, os, sys, time
import numpy as np
import pyomo.environ as pyo

HERE = os.path.dirname(os.path.abspath(__file__))
PRT = os.path.join(HERE, "permafrost_revision_test")
CONFIG = os.environ.get("DICE_CONFIG", "dicecp")
assert CONFIG in ("dicec", "dicecp")
PERMA = CONFIG == "dicecp"
SUB = int(os.environ.get("DICE_SUBSTEPS", "1"))
WARM_START_FILE = os.environ.get("DICE_WARM_START", "controls.json")
OUTPUT_FILE = os.environ.get("DICE_OUTPUT_FILE", f"{CONFIG}_corrected_ss{SUB}_controls.json")
HORIZON_PAD = int(os.environ.get("DICE_HORIZON_PAD", "50"))
N_TOTAL = 400 + HORIZON_PAD
LINEAR_SOLVER = os.environ.get("DICE_LINEAR_SOLVER", "ma27")
HSLLIB = os.environ.get("DICE_HSLLIB", os.path.join(HERE, "coinhsl_extracted",
                        "CoinHSL-archive.v2023.11.17.x86_64-apple-darwin-libgfortran5", "lib", "libcoinhsl.dylib"))
if LINEAR_SOLVER != "mumps":
    _d = os.path.dirname(HSLLIB)
    os.environ["DYLD_LIBRARY_PATH"] = _d + os.pathsep + os.environ.get("DYLD_LIBRARY_PATH", "")
    os.environ["LD_LIBRARY_PATH"] = _d + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")

with contextlib.redirect_stdout(io.StringIO()):
    sys.path.insert(0, HERE)
    import dice_stella_validation as dsv
P = dsv.P_

# ---------------- DICE 2023 constants and GAMS bounds (identical to the previous Chapter 3 scripts) ----------------
gamma_, a3_, theta2_ = 0.300, P['a3'], P['theta2']
fco22x_, Meq_ = P['fco22x'], P['Meq']
Teq1_, Teq2_, d1_, d2_ = P['Teq1'], P['Teq2'], P['d1'], P['d2']
Fc1_, Fc2_, dK_param = P['Fc1'], P['Fc2'], P['dK']
sd_cons_, eta_, pi_param = P['sd_cons'], P['eta'], P['pi_']
scale1_, scale2_ = P['scale1'], P['scale2']
pback2050_, remiss0_, remiss80_ = P['pback2050'], P['remiss0'], P['remiss80']
rartp_cont = math.exp(0.001 + 0.5 * pi_param) - 1
A2 = float(os.environ.get("DICE_A2", "0.003467"))

NPER = 81
pv = np.arange(1, NPER + 1)
miuup = np.zeros(NPER); miuup[0], miuup[1] = 0.05, 0.10
for i in range(NPER):
    if pv[i] > 2: miuup[i] = 0.12 * (pv[i] - 1)
    if pv[i] > 8: miuup[i] = 0.85 + 0.05 * (pv[i] - 8)
    if pv[i] > 11: miuup[i] = 1.0
    if pv[i] > 20: miuup[i] = 1.1
    if pv[i] > 37: miuup[i] = 1.05
    if pv[i] > 57: miuup[i] = 1.0


def exo_at(y):
    L = float(dsv.interp_annual('L', y)); aL = float(dsv.interp_annual('aL', y))
    sig = float(dsv.interp_annual('sigma', y)); eland = float(dsv.interp_annual('eland', y))
    Fmisc = float(dsv.interp_annual('F_Misc', y)); CO2EGHGabateB = float(dsv.interp_annual('CO2E_GHGabateB', y))
    FGHGexo = float(dsv.fghgexo(y))
    pback = pback2050_ * math.exp(-0.01 * (y - 2050)) if y - 2020 <= 30 else pback2050_ * math.exp(-0.001 * (y - 2050))
    remiss = remiss0_ + ((remiss80_ - remiss0_) / 80) * (y - 2020) if y <= 2100 else remiss80_
    return dict(L=L, aL=aL, sig=sig, eland=eland, Fmisc=Fmisc, FGHGexo=FGHGexo, CO2EGHGabateB=CO2EGHGabateB,
                pback=pback, remiss=remiss)


_EXO = {}
def exo(y):
    k = round(y, 6)
    if k not in _EXO:
        _EXO[k] = exo_at(y)
    return _EXO[k]


# ---------------- Lenton carbon cycle constants ----------------
SY = 3.15576e7
RHO = P['rho']
GT2MOL = 1e15 / 12.0
V = dict(l=P['Vl'], h=P['Vh'], i=P['Vi'], d=P['Vd'])
kT, kU, kwi, kid = P['kT'], P['kU'], P['kwi'], P['kid']
AREA = dict(l=P['aO'] * P['Al'], h=P['aO'] * P['Ah'])
KM = 120.0
kMM = 1.478 * (251.0 / (120 + 251.0)) / (251.0 / (KM + 251.0))
ka_fgtm = P['fgtm'] / P['ka']                     # atm per GtC
L20 = json.load(open(os.path.join(PRT, "lenton_km120_state2020.json")))
MASS = {b: V[b] * RHO for b in 'lhid'}


def _chem_consts(T, S=34.7):
    K0 = math.exp(-60.2409 + 93.4517 * (100 / T) + 23.3585 * math.log(T / 100)
                  + S * (0.023517 - 0.023656 * (T / 100) + 0.0047036 * (T / 100) ** 2))
    K1 = 10 ** -(812.27 / T + 3.356 - 0.00171 * S * math.log(T) + 0.000091 * S ** 2)
    K2 = 10 ** -(1450.87 / T + 4.604 - 0.00385 * S * math.log(T) + 0.000182 * S ** 2)
    KB = math.exp((-8966.90 - 2890.53 * S ** .5 - 77.942 * S + 1.728 * S ** 1.5 - 0.0996 * S ** 2) / T
                  + 148.0248 + 137.1942 * S ** .5 + 1.62142 * S - (24.4344 + 25.085 * S ** .5 + 0.2474 * S) * math.log(T)
                  + 0.053105 * S ** .5 * T)
    return K0, K1, K2, KB, 0.000416 * S / 35.0


def _pco2_exact(dic, ta, T):
    K0, K1, K2, KB, TB = _chem_consts(T)
    lo, hi = 1e-11, 1e-5
    f = lambda H: dic * (K1 * H + 2 * K1 * K2) / (H * H + K1 * H + K1 * K2) + TB * KB / (KB + H) - ta
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if f(lo) * f(mid) <= 0: hi = mid
        else: lo = mid
    H = math.sqrt(lo * hi)
    return dic * H * H / (H * H + K1 * H + K1 * K2) / K0


# polynomial fit of ln(pCO2_s) in x = DIC/DIC20 - 1 for each surface box (TA and T fixed)
XLO, XHI, DEG = -0.10, 0.06, 10
CHEM = {}
for b, T in (('l', P['Tl']), ('h', P['Th'])):
    D20 = L20[f'DIC_{b}']                                  # mol
    ta = L20[f'TA_{b}'] / MASS[b]
    xs = np.linspace(XLO, XHI, 400)
    ys = np.array([math.log(_pco2_exact(D20 * (1 + x) / MASS[b], ta, T)) for x in xs])
    us = (2 * xs - (XLO + XHI)) / (XHI - XLO)               # scaled to [-1, 1]; Chebyshev fit, power-basis evaluation
    coef = np.polynomial.chebyshev.cheb2poly(np.polynomial.chebyshev.chebfit(us, ys, DEG))
    err = np.max(np.abs(np.exp(np.polynomial.polynomial.polyval(us, coef) - ys) - 1))
    CHEM[b] = dict(D20_gtc=D20 * 12e-15, coef=[float(c) for c in coef], relerr=err,
                   g=P['ks'] * AREA[b] * RHO * _chem_consts(T)[0] * 12e-15)   # GtC/yr per atm


def pco2_surface(b, dic_gtc):
    c = CHEM[b]; x = dic_gtc / c['D20_gtc'] - 1.0
    u = (2 * x - (XLO + XHI)) / (XHI - XLO)
    poly = 0.0
    for k in reversed(range(len(c['coef']))):
        poly = poly * u + c['coef'][k]
    return pyo.exp(poly)


# constant biological source/sink (GtC/yr) per box, from the steady nitrate field (Knox & McElroy Table 1)
def _bio_const():
    N = {b: L20[f'NO3_{b}'] / MASS[b] for b in 'lhid'}
    Nin_l = (kT + kwi) * SY * RHO * N['i']
    Nin_h = kT * SY * RHO * N['l'] + kU * SY * RHO * N['i']
    d = dict(l=0.0, h=0.0, i=0.0, d=0.0)
    for b, n in (('l', Nin_l), ('h', 0.01 * Nin_h)):
        corg = 7.0 * n; caco3 = 0.25 * corg
        d[b] -= (corg + caco3) * 12e-15
        d['i'] += 0.3 * corg * 12e-15
        d['d'] += (0.7 * corg + caco3) * 12e-15
    return d
BIO = _bio_const()
QF = {('i', 'l'): kT + kwi, ('l', 'i'): kwi, ('l', 'h'): kT, ('i', 'h'): kU, ('h', 'd'): kT + kU, ('d', 'i'): kT + kU + kid, ('i', 'd'): kid}

# ---------------- permafrost constants ----------------
WF, MU_PF, SIG_PF = 2.0, 1.67, 0.986
FSTATIC, Q10, RATE = P['fstatic'], P['Q10'], 0.02
REFREEZE_TAU = 2.0
AB_T0, AB_W, AB_TAU, AB_FRAC = 1.5, 0.5, 200.0, 0.5


def T_hist(y, kk=0.02):
    return 0.0 if y < 1850 else 1.25 * (math.exp(kk * (y - 1850)) - 1) / (math.exp(170 * kk) - 1)


def _phi(z):
    return 0.5 * (1.0 + pyo.tanh(math.sqrt(2.0 / math.pi) * (z + 0.044715 * z ** 3)))


def dff_dT(T):
    c = math.sqrt(2.0 / math.pi)
    T = 0.5 * (T + pyo.sqrt(T ** 2 + 1e-4))      # smooth positive floor (~0.005 C): keeps log() defined if a
                                                 # warm-start path cools to ~0 C; no effect for T > ~0.05 C
    z = (pyo.log(WF * T) - MU_PF) / SIG_PF
    u = c * (z + 0.044715 * z ** 3)
    return -0.5 * (1 - pyo.tanh(u) ** 2) * c * (1 + 3 * 0.044715 * z ** 2) / (SIG_PF * T)


def smax0(x, eps=1e-6):
    return 0.5 * (x + pyo.sqrt(x ** 2 + eps))


STATE_NAMES = ['Tbox1', 'Tbox2', 'K', 'Ctot', 'Uacc', 'FGHG', 'Ca', 'Cv', 'Cs', 'Dcum', 'Dl', 'Dh', 'Di', 'Dd']
if PERMA:
    STATE_NAMES += ['Cperm', 'thawedC', 'Aab', 'Ggr']
IC = dict(Tbox1=0.1477, Tbox2=1.099454, K=295.0, Ctot=633.5, Uacc=0.0, FGHG=0.518,
          Ca=L20['Ca'], Cv=L20['Cv'], Cs=L20['Cs'], Dcum=L20['Dcum'],
          Dl=L20['DIC_l'] * 12e-15, Dh=L20['DIC_h'] * 12e-15, Di=L20['DIC_i'] * 12e-15, Dd=L20['DIC_d'] * 12e-15,
          **{k: v for k, v in json.load(open(os.path.join(PRT, "perma_state2020.json"))).items() if k != 'note'})
# permafrost pools in 2020 from the 1800-2020 spin-up (865 PgC pre-industrial, Hector), as in ch2_models.perma_2020_state


def rhs(st, miu, s, y, ex, T200, a2=A2):
    st = dict(zip(STATE_NAMES, st))
    Tb1, Tb2, K = st['Tbox1'], st['Tbox2'], st['K']
    TATM = Tb1 + Tb2
    Yg = ex['aL'] * (ex['L'] / 1000.0) ** (1 - gamma_) * K ** gamma_
    miu_safe = 0.5 * (miu + pyo.sqrt(miu ** 2 + 1e-12))
    Cab = Yg * ex['sig'] * ex['remiss'] * ex['pback'] / theta2_ / 1000.0 * miu_safe ** theta2_
    Y = Yg * (1 - a2 * TATM ** a3_) - Cab
    I = s * Y; Cons = Y - I
    ECo2 = (ex['sig'] * Yg + ex['eland']) * (1 - miu)
    Ca = st['Ca']
    F = fco22x_ * (pyo.log(Ca / Meq_) / math.log(2)) + ex['Fmisc'] + ex['FGHGexo']
    dTb1 = (Teq1_ * F - Tb1) / d1_; dTb2 = (Teq2_ * F - Tb2) / d2_
    Vc = min(sd_cons_ ** 2 * (y - 2020), sd_cons_ ** 2 * 235)
    STP = 1.0 / (1 + rartp_cont) ** (y - 2020) / (1 - 0.5 * Vc * eta_ ** 2) ** (y - 2020)
    c = Cons * 1000.0 / ex['L']
    c_safe = 0.5 * (c + pyo.sqrt(c ** 2 + 1e-12))
    U = (c_safe ** (1 - eta_) - 1) / (1 - eta_) - 1
    # land
    e_ff = ex['sig'] * Yg * (1 - miu) * 0.2727
    e_luc = ex['eland'] * (1 - miu) * 0.2727
    ppm = Ca * ka_fgtm * 1e6
    Cvs = P['Cv0'] - P['kD'] * st['Dcum']
    Photo = P['kp'] * Cvs * kMM * ((ppm - P['kc']) / (KM + ppm - P['kc'])) * ((15 + TATM) ** 2 * (25 - TATM) / 5625)
    PlantR = P['kr'] * st['Cv'] * P['kA'] * pyo.exp(-P['Ea'] / (P['R'] * (288 + TATM)))
    Turn = P['kt'] * st['Cv']
    SoilR = P['ksr'] * st['Cs'] * P['kB'] * pyo.exp(math.log(2.7) * (-(308.56 / (288 + TATM - 227.13))))
    # ocean (GtC)
    C = dict(l=st['Dl'], h=st['Dh'], i=st['Di'], d=st['Dd'])
    dO = dict(l=BIO['l'], h=BIO['h'], i=BIO['i'], d=BIO['d'])
    for (src, dst), q in QF.items():
        f = q * SY * C[src] / V[src]
        dO[src] = dO[src] - f; dO[dst] = dO[dst] + f
    pa = Ca * ka_fgtm
    Fl = CHEM['l']['g'] * (pa - pco2_surface('l', C['l']))
    Fh = CHEM['h']['g'] * (pa - pco2_surface('h', C['h']))
    dO['l'] = dO['l'] + Fl; dO['h'] = dO['h'] + Fh
    out = [dTb1, dTb2, I - dK_param * K, ECo2 / 3.666, U * ex['L'] * STP * scale1_,
           Fc1_ * ex['CO2EGHGabateB'] * (1 - miu) - (-math.log(Fc2_)) * st['FGHG']]
    Fmelt = 0.0
    if PERMA:
        raw = -dff_dT(TATM) * (dTb1 + dTb2) * st['Cperm']
        lim = -st['thawedC'] / REFREEZE_TAU
        grad = 0.5 * (raw + lim + pyo.sqrt((raw - lim) ** 2 + 1e-9))
        w = 0.5 * (1 + pyo.tanh((TATM - AB_T0) / AB_W))
        ab = smax0(AB_FRAC * w * smax0(st['Ggr']) - st['Aab']) / AB_TAU
        Fmelt = RATE * (1 - FSTATIC) * st['thawedC'] * pyo.exp(math.log(Q10) * T200 / 10.0)
        perm = [-(grad + ab), grad + ab - Fmelt, ab, grad]
    out += [-Photo + PlantR + SoilR - Fl - Fh + e_ff + e_luc + Fmelt, Photo - PlantR - Turn - e_luc, Turn - SoilR, e_luc,
            dO['l'], dO['h'], dO['i'], dO['d']]
    if PERMA:
        out += perm
    return out


def t200_value(Tpath_func, n):
    """200-yr running mean of WF*T over years n-200 .. n-1 (annual), pre-2020 from T_hist."""
    tot = 0.0
    for k in range(n - 200, n):
        tot = tot + WF * (T_hist(2020 + k) if k < 0 else Tpath_func(k))
    return tot / 200.0


# ---------------- numpy forward simulation (warm start and checks) ----------------
def load_controls(path):
    d = json.load(open(path))
    miu = d["miu"]["ypts"] if isinstance(d["miu"], dict) else d["miu"]
    s = d["s"]["ypts"] if isinstance(d["s"], dict) else d["s"]
    return np.array(miu), np.array(s)


def controls_fn(path, step=True):
    mu, s = load_controls(os.path.join(HERE, path))
    if len(mu) == 400:
        f = lambda t, a: float(a[min(max(int(t) - 2020, 0), 399)])
        return (lambda t: f(t, mu)), (lambda t: f(t, s))
    yrs = np.linspace(2020, 2420, 81)
    if step:   # annual step, value at the start of each year (as the Pyomo decision variables)
        return (lambda t: float(np.interp(math.floor(t), yrs, mu))), (lambda t: float(np.interp(math.floor(t), yrs, s)))
    return (lambda t: float(np.interp(t, yrs, mu))), (lambda t: float(np.interp(t, yrs, s)))


def forward(miu_f, s_f, nyears=N_TOTAL, sub=1, record=None):
    """RK4 with `sub` substeps per year, T200 updated once per year (as in the Pyomo model).
    record: optional dict filled with every substep's start state, RK4 slopes and stage states (warm start)."""
    S = np.array([IC[n] for n in STATE_NAMES], dtype=float)
    traj = [S.copy()]; Ts = []
    h = 1.0 / sub
    for n in range(nyears):
        Ts.append(S[0] + S[1])
        T200 = t200_value(lambda k: Ts[k], n) if PERMA else 0.0
        for j in range(sub):
            y = 2020 + n + j * h
            mu, sv = miu_f(y), s_f(y)
            k1 = np.array(rhs(S, mu, sv, y, exo(y), T200))
            k2 = np.array(rhs(S + h / 2 * k1, mu, sv, y + h / 2, exo(y + h / 2), T200))
            k3 = np.array(rhs(S + h / 2 * k2, mu, sv, y + h / 2, exo(y + h / 2), T200))
            k4 = np.array(rhs(S + h * k3, mu, sv, y + h, exo(y + h), T200))
            if record is not None:
                record[(n, j)] = dict(S=S.copy(), K=(k1, k2, k3, k4), H=(S + h / 2 * k1, S + h / 2 * k2, S + h * k3))
            S = S + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        traj.append(S.copy())
    return np.array(traj)


print(f"[{CONFIG}] chemistry fit max rel. error in pCO2: warm {CHEM['l']['relerr']:.2e}, cold {CHEM['h']['relerr']:.2e} "
      f"(DIC fit range {XLO:+.0%}..{XHI:+.0%} of 2020)", flush=True)

if os.environ.get("DICE_CHECK_ONLY") == "1" and __name__ == "__main__":
    # forward checks: annual RK4 (sub = 1, 2, 5) vs fine reference (sub = 20), GAMS controls
    mf, sf = controls_fn(WARM_START_FILE, step=True)
    ref = forward(mf, sf, nyears=400, sub=20)
    mfl, sfl = controls_fn(WARM_START_FILE, step=False)
    refl = forward(mfl, sfl, nyears=400, sub=20)
    _T = refl[:, 0] + refl[:, 1]
    print(f"fine reference, linearly interpolated controls (as Chapter 2): peak T {_T.max():.4f} ({2020 + _T.argmax()}), "
          f"U_acc(2420) {refl[400, STATE_NAMES.index('Uacc')]:.3f}")
    iT = [STATE_NAMES.index('Tbox1'), STATE_NAMES.index('Tbox2')]; iU = STATE_NAMES.index('Uacc'); iC = STATE_NAMES.index('Ca')
    Tref = ref[:, iT[0]] + ref[:, iT[1]]
    print(f"fine reference (sub=20): peak T {Tref.max():.4f} ({2020 + Tref.argmax()}), U_acc(2420) {ref[400, iU]:.3f}, "
          f"Ca2100 {ref[80, iC]:.1f}")
    for sub in (1, 2, 5):
        t0 = time.time()
        try:
            with np.errstate(all='ignore'):
                tr = forward(mf, sf, nyears=400, sub=sub)
        except (OverflowError, ValueError) as e:
            print(f"sub={sub}: UNSTABLE ({e.__class__.__name__}) -- annual RK4 step too large for the surface-ocean timescale")
            continue
        T = tr[:, iT[0]] + tr[:, iT[1]]
        print(f"sub={sub}: max|dT| {np.max(np.abs(T - Tref)):.2e}  dU_acc {tr[400, iU] - ref[400, iU]:+.4f}  "
              f"max|dCa| {np.max(np.abs(tr[:, iC] - ref[:, iC])):.3f}  ({time.time() - t0:.1f}s)")
    lo = {b: (1 + XLO) * CHEM[b]['D20_gtc'] for b in 'lh'}; hi = {b: (1 + XHI) * CHEM[b]['D20_gtc'] for b in 'lh'}
    for b, nm in (('l', 'Dl'), ('h', 'Dh')):
        x = ref[:, STATE_NAMES.index(nm)]
        print(f"DIC_{b} range along path: {x.min() / CHEM[b]['D20_gtc'] - 1:+.3f} .. {x.max() / CHEM[b]['D20_gtc'] - 1:+.3f} (fit {XLO:+.2f}..{XHI:+.2f})")
    np.save(os.path.join(PRT, f"pyomo_corrected_{CONFIG}_gamsref.npy"), ref)
    sys.exit(0)

# ---------------- Pyomo model ----------------
print(f"[{CONFIG}] building warm start from {WARM_START_FILE} ...", flush=True)
mf0, sf = controls_fn(WARM_START_FILE, step=True)
period_of = lambda n: min(1 + n // 5, 80)
miuup_annual = np.array([miuup[period_of(n) - 1] for n in range(N_TOTAL)])
mf = lambda t: min(max(mf0(t), 0.0), float(miuup_annual[min(int(t) - 2020, N_TOTAL - 1)]))   # warm start obeys the MIU bounds
REC = {}
warm = forward(mf, sf, nyears=N_TOTAL, sub=SUB, record=REC)

m = pyo.ConcreteModel()
m.N = pyo.RangeSet(0, N_TOTAL)
m.MIU = pyo.Var(range(N_TOTAL), bounds=lambda m, n: (0, float(miuup_annual[n])),
                initialize=lambda m, n: min(max(mf(2020 + n), 0.0), float(miuup_annual[n])))
m.S = pyo.Var(range(N_TOTAL), bounds=(0, 1), initialize=lambda m, n: sf(2020 + n))
for n in range(N_TOTAL):
    if period_of(n) > 37:
        m.S[n].fix(0.28)

bounds = dict(Tbox1=(-2, 25), Tbox2=(-2, 25), K=(1, None), Ctot=(0, None), Uacc=(None, None), FGHG=(None, None),
              Ca=(300, None), Cv=(1, None), Cs=(1, None), Dcum=(None, None),
              Dl=((1 + XLO) * CHEM['l']['D20_gtc'], (1 + XHI) * CHEM['l']['D20_gtc']),
              Dh=((1 + XLO) * CHEM['h']['D20_gtc'], (1 + XHI) * CHEM['h']['D20_gtc']),
              Di=(1, None), Dd=(1, None), Cperm=(0, 1300.0), thawedC=(-0.5, 1300.0), Aab=(-1, None), Ggr=(None, None))
SV = {}
for k, nm in enumerate(STATE_NAMES):
    v = pyo.Var(m.N, bounds=bounds[nm], initialize=lambda m, n, _k=k: float(warm[n][_k]))
    setattr(m, nm, v); SV[nm] = v
    v[0].fix(IC[nm])
m.TATM = pyo.Var(m.N, bounds=(float(os.environ.get('DICE_TATM_LO', '0.5')), 20), initialize=lambda m, n: float(warm[n][0] + warm[n][1]))
m.tatm_def = pyo.Constraint(m.N, rule=lambda m, n: m.TATM[n] == SV['Tbox1'][n] + SV['Tbox2'][n])

print(f"[{CONFIG}] building RK4 defects ({N_TOTAL} yr x {SUB} substeps x {len(STATE_NAMES)} states) ...", flush=True)
defs = []
h = 1.0 / SUB
m.SUBST = pyo.Var(range(N_TOTAL), range(1, SUB), range(len(STATE_NAMES)),
                  initialize=lambda m, n, j, k: float(REC[(n, j)]['S'][k]))
m.KS = pyo.Var(range(N_TOTAL), range(SUB), range(1, 5), range(len(STATE_NAMES)),
               initialize=lambda m, n, j, st, k: float(REC[(n, j)]['K'][st - 1][k]))
m.HS = pyo.Var(range(N_TOTAL), range(SUB), range(1, 4), range(len(STATE_NAMES)),
               initialize=lambda m, n, j, st, k: float(REC[(n, j)]['H'][st - 1][k]))
NS = len(STATE_NAMES)
for n in range(N_TOTAL):
    T200 = t200_value(lambda k: m.TATM[k], n) if PERMA else 0.0
    for j in range(SUB):
        y = 2020 + n + j * h
        S0 = [SV[nm][n] for nm in STATE_NAMES] if j == 0 else [m.SUBST[n, j, k] for k in range(NS)]
        S1 = [SV[nm][n + 1] for nm in STATE_NAMES] if j == SUB - 1 else [m.SUBST[n, j + 1, k] for k in range(NS)]
        Ks = [[m.KS[n, j, st, k] for k in range(NS)] for st in range(1, 5)]
        Hs = [[m.HS[n, j, st, k] for k in range(NS)] for st in range(1, 4)]
        r1 = rhs(S0, m.MIU[n], m.S[n], y, exo(y), T200)
        defs += [(Ks[0][k], r1[k]) for k in range(NS)]
        defs += [(Hs[0][k], S0[k] + 0.5 * h * Ks[0][k]) for k in range(NS)]
        r2 = rhs(Hs[0], m.MIU[n], m.S[n], y + h / 2, exo(y + h / 2), T200)
        defs += [(Ks[1][k], r2[k]) for k in range(NS)]
        defs += [(Hs[1][k], S0[k] + 0.5 * h * Ks[1][k]) for k in range(NS)]
        r3 = rhs(Hs[1], m.MIU[n], m.S[n], y + h / 2, exo(y + h / 2), T200)
        defs += [(Ks[2][k], r3[k]) for k in range(NS)]
        defs += [(Hs[2][k], S0[k] + h * Ks[2][k]) for k in range(NS)]
        r4 = rhs(Hs[2], m.MIU[n], m.S[n], y + h, exo(y + h), T200)
        defs += [(Ks[3][k], r4[k]) for k in range(NS)]
        defs += [(S1[k], S0[k] + h / 6 * (Ks[0][k] + 2 * Ks[1][k] + 2 * Ks[2][k] + Ks[3][k])) for k in range(NS)]
m.idx = pyo.RangeSet(0, len(defs) - 1)
m.defs = pyo.Constraint(m.idx, rule=lambda m, i: defs[i][0] == defs[i][1])
m.UTILITY = pyo.Objective(expr=SV['Uacc'][N_TOTAL] + scale2_, sense=pyo.maximize)
print(f"  {len(defs)} defect constraints", flush=True)

if __name__ == "__main__":
    solver = pyo.SolverFactory("ipopt")
    solver.options["max_iter"] = int(os.environ.get("DICE_MAX_ITER", "25000"))
    solver.options["mu_strategy"] = "adaptive"
    solver.options["max_cpu_time"] = 200000
    if LINEAR_SOLVER != "mumps":
        solver.options["linear_solver"] = LINEAR_SOLVER
        solver.options["hsllib"] = HSLLIB
    t0 = time.time()
    res = solver.solve(m, tee=True)
    print("Solver status:", res.solver.status, res.solver.termination_condition, f"({time.time() - t0:.0f}s)")
    print(f"Objective (padded): {pyo.value(m.UTILITY):.6f}; U_acc(2420)+scale2: {pyo.value(SV['Uacc'][400]) + scale2_:.6f}")
    T = np.array([pyo.value(m.TATM[n]) for n in range(401)])
    print(f"peak T {T.max():.4f} ({2020 + T.argmax()}); T2100 {T[80]:.4f}; T2300 {T[280]:.4f}")
    out = dict(config=CONFIG, substeps=SUB, miu=[pyo.value(m.MIU[n]) for n in range(400)],
               s=[pyo.value(m.S[n]) for n in range(400)], objective=pyo.value(m.UTILITY))
    json.dump(out, open(os.path.join(HERE, OUTPUT_FILE), "w"))
    tf = os.environ.get("DICE_TRAJ_FILE")
    if tf:
        json.dump({"years": list(range(2020, 2021 + N_TOTAL)),
                   **{nm: [pyo.value(SV[nm][n]) for n in range(N_TOTAL + 1)] for nm in STATE_NAMES}},
                  open(os.path.join(HERE, tf), "w"))
    print("saved", OUTPUT_FILE)

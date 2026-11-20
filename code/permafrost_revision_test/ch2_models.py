"""
Unified Chapter 2 forward model (2026-11-16): DICE 2023 economy/climate/welfare (as in
permafrost_revision_sim.rhs) with a choice of carbon cycle and permafrost module.

  carbon = 'dfair'  : DICE 2023 DFAIR (four reservoirs, state-dependent alpha), from the DICE 2020 state
  carbon = 'lenton' : full Lenton (2000) model, variant (iv) (ocean surface warms with global temperature, Lenton eq. 16)
                      with Lenton's own K_M = 145 and k_D = 0.27 (LENTON_VARIANT=iii: his control, no ocean warming, K_M = 120),
                      Knox & McElroy ocean biology/chemistry; 2020 state from a run starting at the
                      1800 steady state with Lenton's Table 3 emissions to 1990 and GCB emissions 1990-2019
  perma  = None or dict(m=multiplier, abrupt='none'|'smooth'):
                      Hector v2.3pf-faithful module (730 PgC in 2020, ln(wf T) frozen fraction,
                      200-yr running-mean T200 with sticky Q10); release rate = m * 0.02/yr (m = 1 is Hector;
                      m is the strength knob for the literature-range sweep). Permafrost carbon enters the
                      atmosphere (Lenton) or is added to emissions entering the DFAIR reservoirs.
Controls: GAMS optimal MIU and S (controls.json), as in Chapter 2. Optional CO2 pulse for the SCC.
"""
import json, math, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import permafrost_revision_sim as prs
import permafrost_hector_faithful_sim as phf
import lenton_rebuild as LR
import lenton_full as LF
import lenton_bau_replication as BAU
import ssp_check as SC

P = prs.P
VARIANT = os.environ.get('LENTON_VARIANT', 'iv')
assert VARIANT in ('iii', 'iv')
OW = VARIANT == 'iv'                 # ocean warming
KM = 145.0 if OW else 120.0
ESH, TAU = SC.ESH, SC.TAU
XKEYS = [f"{l}_{b}" for l in ('DIC', 'TA', 'NO3') for b in 'lhid']


def lenton_2020_state():
    fn = os.path.join(HERE, "lenton_km120_state2020.json" if not OW else "lenton_km145_iv_state2020.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    E = pd.read_csv(os.path.join(HERE, "hist_data", "world_co2_1850_2023.csv"), index_col=0)
    Yss = np.load(os.path.join(HERE, "lenton_full_steadystate.npy"))
    Ph, _, _, _ = LR.land(Yss[0], 1.0, 1.0, 0.0, 0.0, KM)
    rP = P['kr'] * P['kA'] * math.exp(-P['Ea'] / (P['R'] * 288))
    Cv = Ph / (rP + P['kt']); Cs = P['kt'] * Cv / (P['ksr'] * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 - 227.13)))))
    S = np.concatenate([[Yss[0], Cv, Cs, 0.0], Yss[1:]]); cfg = dict(ocean_warming=OW, KM=KM); dt = 0.1
    for y in range(1800, 2020):
        if y < 1990:
            et, el = BAU.interp(BAU.TOT, y), BAU.interp(BAU.LUC, y)
        else:
            el = E.loc[y, 'land_use_change_co2'] / 3664.0; et = E.loc[y, 'co2'] / 3664.0 + el
        for j in range(int(1 / dt)):
            S = LF.rk4(LF.rhs, S, dt, LR.T_hist(y + j * dt), et - el, el, cfg)
    st = dict(Ca=S[0], Cv=S[1], Cs=S[2], Dcum=S[3], **{k: S[4 + n] for n, k in enumerate(XKEYS)})
    json.dump({k: float(v) for k, v in st.items()}, open(fn, "w"), indent=1)
    return st


L20 = lenton_2020_state()
_fl20 = LF.ocean_rhs(L20['Ca'], np.array([L20[k] for k in XKEYS]), 0.0)[2]
PS20 = {'l': _fl20['l_pco2'], 'h': _fl20['h_pco2']}            # 2020 surface-ocean pCO2 (atm), for the factorial
# ---------------- feedback switches referenced to pre-industrial (factorial / fertilisation sweep) ----------------
# cfg['fbx'] = {ocean, fert, photoT, vegT, soilT: 0/1, default 1; oceanT: 0/1, default 1 in variant iv, 0 in iii},
# cfg['KM'] = Michaelis-Menten constant (fertilisation strength).
# 'Off' means the process keeps its pre-industrial response throughout, including the 1800-2020 spin-up:
#   fert   : photosynthesis sees pre-industrial CO2;   photoT/vegT/soilT : temperature anomaly 0 in that process;
#   ocean  : surface pCO2 linear in DIC and alkalinity about the pre-industrial state (constant buffer factor,
#            i.e. no buffer saturation).  Each switch combination gets its own 2020 state from its own spin-up.
_YSS = np.load(os.path.join(HERE, "lenton_full_steadystate.npy"))
CA_PI = float(_YSS[0])
_MBOX = np.array([LF.MASS[b] for b in LF.B])
_C0 = {b: (_YSS[1 + LF.IDX[b]] / _MBOX[LF.IDX[b]], _YSS[5 + LF.IDX[b]] / _MBOX[LF.IDX[b]]) for b in 'lh'}
_LIN = {}
for _b, _T in (('l', LR.Tl0), ('h', LR.Th0)):
    _d, _a = _C0[_b]; _h = 1e-6
    _p0 = LR.pco2_surface(_d, _a, _T)
    _LIN[_b] = (_p0, (LR.pco2_surface(_d * (1 + _h), _a, _T) - LR.pco2_surface(_d * (1 - _h), _a, _T)) / (2 * _d * _h),
                (LR.pco2_surface(_d, _a * (1 + _h), _T) - LR.pco2_surface(_d, _a * (1 - _h), _T)) / (2 * _a * _h))


def ocean_rhs_lin(Ca, X, dT=0.0):
    """LF.ocean_rhs with surface pCO2 linearised in DIC and alkalinity about the pre-industrial state (no buffer
    saturation); the temperature dependence of solubility and chemistry (dT > 0) is kept exact."""
    Dm, Am, Nm = X[0:4], X[4:8], X[8:12]
    cD, cA, cN = Dm / _MBOX, Am / _MBOX, Nm / _MBOX
    pa = Ca / P['ka'] * P['fgtm']
    fl = {}
    for b, T0 in (('l', LR.Tl0), ('h', LR.Th0)):
        k = LF.IDX[b]; p0, gD, gA = _LIN[b]; T = T0 + dT
        ps = p0 + gD * (cD[k] - _C0[b][0]) + gA * (cA[k] - _C0[b][1])
        if dT:
            ps += LR.pco2_surface(_C0[b][0], _C0[b][1], T) - p0          # warming shifts the chemistry
        g = P['ks'] * LF.A[b] * LF.RHO * LR.K0(T)
        fl[b] = g * (pa - ps)
    dD, dA, dN, _ = LF.biology(cN)
    dD = dD + LF.advect(cD); dA = dA + LF.advect(cA); dN = dN + LF.advect(cN)
    dD[0] += fl['l']; dD[1] += fl['h']
    return np.concatenate([dD, dA, dN]), (fl['l'] + fl['h']) * 12e-15


def lenton_fluxes(Ca, Cv, Cs, Dc, X, T, fbx, km):
    on = lambda k: fbx.get(k, 1)
    Ph, _, _, _ = LR.land(Ca if on('fert') else CA_PI, Cv, Cs, T if on('photoT') else 0.0, Dc, km)
    _, Rp, Tu, _ = LR.land(Ca, Cv, Cs, T if on('vegT') else 0.0, Dc, km)
    _, _, _, Rs = LR.land(Ca, Cv, Cs, T if on('soilT') else 0.0, Dc, km)
    dTo = T if fbx.get('oceanT', int(OW)) else 0.0    # ocean warming: default follows the Lenton variant
    if on('ocean'):
        dX, asg, _, _ = LF.ocean_rhs(Ca, X, dTo)
    else:
        dX, asg = ocean_rhs_lin(Ca, X, dTo)
    return Ph, Rp, Tu, Rs, dX, asg


_L20_CACHE = {}


def lenton_2020_state_fb(fbx, km):
    """2020 Lenton state for a feedback-switch combination and K_M: same spin-up as lenton_2020_state (1800 steady
    state, Lenton Table 3 emissions to 1990, GCB emissions 1990-2019, historical warming), with the switches applied."""
    key = (tuple(sorted(fbx.items())), float(km))
    if key in _L20_CACHE:
        return _L20_CACHE[key]
    E = pd.read_csv(os.path.join(HERE, "hist_data", "world_co2_1850_2023.csv"), index_col=0)
    Ph, _, _, _ = LR.land(_YSS[0], 1.0, 1.0, 0.0, 0.0, km)
    rP = P['kr'] * P['kA'] * math.exp(-P['Ea'] / (P['R'] * 288))
    Cv = Ph / (rP + P['kt']); Cs = P['kt'] * Cv / (P['ksr'] * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 - 227.13)))))
    S = np.concatenate([[_YSS[0], Cv, Cs, 0.0], _YSS[1:]]); dt = 0.1

    def f(S, T, eff, eluc):
        Ca, Cv, Cs, Dc = S[:4]
        Ph, Rp, Tu, Rs, dX, asg = lenton_fluxes(Ca, Cv, Cs, Dc, S[4:], T, fbx, km)
        return np.concatenate([[-Ph + Rp + Rs - asg + eff + eluc, Ph - Rp - Tu - eluc, Tu - Rs, eluc], dX])
    for y in range(1800, 2020):
        if y < 1990:
            et, el = BAU.interp(BAU.TOT, y), BAU.interp(BAU.LUC, y)
        else:
            el = E.loc[y, 'land_use_change_co2'] / 3664.0; et = E.loc[y, 'co2'] / 3664.0 + el
        for j in range(int(1 / dt)):
            T = LR.T_hist(y + j * dt); a = (T, et - el, el)
            k1 = f(S, *a); k2 = f(S + dt / 2 * k1, *a); k3 = f(S + dt / 2 * k2, *a); k4 = f(S + dt * k3, *a)
            S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    st = dict(Ca=S[0], Cv=S[1], Cs=S[2], Dcum=S[3], **{k: S[4 + n] for n, k in enumerate(XKEYS)})
    _L20_CACHE[key] = {k: float(v) for k, v in st.items()}
    return _L20_CACHE[key]


# state layout: 0-5 economy/climate [Tb1,Tb2,K,Ctot_econ,Uacc,FGHG]; carbon block; permafrost block [Cp,Th,A,G,R]
NE = 6


_PF_CACHE = {}


def perma_2020_state(pm):
    """Permafrost pools in 2020 from a spin-up 1800-2020 (Hector-faithful module, same equations and strength m as the
    run), starting from Hector's pre-industrial pool (865 PgC, Woodard et al. 2021) with no thawed carbon, forced by the
    historical warming reconstruction T_hist. Pre-2020 permafrost emissions are not added to the atmosphere: they are
    part of the observed historical carbon budget to which the carbon cycles are initialised.
    Returns [C_perm, C_thawed, A, G] in 2020."""
    key = (round(pm.get('m', 1.0), 6), pm.get('abrupt'), pm.get('C0pre', 865.0), bool(pm.get('abrupt_rev')))
    if key in _PF_CACHE:
        return _PF_CACHE[key]
    kk, dt = 0.02, 0.05
    Th_ = lambda y: phf.T_hist(y)
    dTdt = lambda y: 0.0 if y < 1850 else 1.25 * kk * math.exp(kk * (y - 1850)) / (math.exp(170 * kk) - 1)
    # 200-yr running mean of wf*T (T = 0 before 1850), evaluated annually
    yrs = np.arange(1600, 2021)
    Ta = np.array([Th_(y) for y in yrs])
    t200 = {int(y): phf.WF * Ta[(yrs >= y - 200) & (yrs < y)].mean() for y in range(1800, 2021)}
    S = np.array([pm.get('C0pre', 865.0), 0.0, 0.0, 0.0])
    f = lambda S, y, T200: perma_flows(Th_(y), dTdt(y), S[0], S[1], S[2], S[3], T200, pm)[0][:4]
    for k in range(int(220 / dt)):
        y = 1800 + k * dt; T200 = t200[int(y)]
        k1 = f(S, y, T200); k2 = f(S + dt / 2 * k1, y + dt / 2, T200); k3 = f(S + dt / 2 * k2, y + dt / 2, T200); k4 = f(S + dt * k3, y + dt, T200)
        S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    _PF_CACHE[key] = [float(v) for v in S]
    return _PF_CACHE[key]


def init(cfg):
    eco = [0.1477, 1.099454, 295.0, 633.5, 0.0, 0.518]
    if cfg['carbon'] == 'dfair':
        carb = [150.093, 102.698, 39.534, 6.1865, 633.5]                 # R0..R3, cumulative emissions
    else:
        L = lenton_2020_state_fb(cfg.get('fbx', {}), cfg.get('KM', KM)) if ('fbx' in cfg or 'KM' in cfg) else L20
        carb = [L['Ca'], L['Cv'], L['Cs'], L['Dcum']] + [L[k] for k in XKEYS]
    perm = (perma_2020_state(cfg['perma']) + [0.0]) if cfg.get('perma') else [0.0] * 5
    return np.array(eco + carb + perm, dtype=float)


def alpha_of(cacc, T):
    return SC.alpha_of(cacc, T)


def perma_flows(T, dT, Cp, Th, A, G, T200, pm):
    dff = phf.dff_dT_hector(T)
    raw = -dff * dT * Cp
    lim = -Th / prs.REFREEZE_TAU
    grad = 0.5 * (raw + lim + math.sqrt((raw - lim) ** 2 + 1e-9))
    if pm.get('abrupt') == 'smooth':
        tgt = prs.AB_FRAC * prs.w_onset(T) * prs.smax0(G)
        ab = ((tgt - A) if pm.get('abrupt_rev') else prs.smax0(tgt - A)) / prs.AB_TAU
    else:
        ab = 0.0
    tot = grad + ab
    if pm.get('abrupt_rev'):            # reversible abrupt thaw: refreeze limited by the thawed pool, as for gradual
        tot = 0.5 * (tot + lim + math.sqrt((tot - lim) ** 2 + 1e-9)); ab = tot - grad
    J = (1 - prs.fstatic_) * Th * 0.02 * pm.get('m', 1.0) * prs.Q10_ ** (T200 / 10.0)
    return np.array([-tot, tot - J, ab, grad, J]), J


def rhs(S, y, T200, cfg, pulse_gtc=0.0):
    Tb1, Tb2, K, Ctot, Uacc, FGHG = S[:NE]
    L, aL, sig, eland, Fmisc, FGHGexo, CO2EGHGabateB, pback, remiss, elandL, miu, s = prs.exo(y)
    if cfg.get('ctrl') is not None:            # optional annual control path (e.g. a re-optimised solution), step per year
        _m, _s = cfg['ctrl']; _i = min(max(int(y) - 2020, 0), len(_m) - 1); miu, s = float(_m[_i]), float(_s[_i])
    T = Tb1 + Tb2
    Yg = aL * (L / 1000.0) ** (1 - prs.gamma_) * K ** prs.gamma_
    Dfrac = cfg['a2'] * T ** prs.a3_
    Cab = Yg * sig * remiss * pback / prs.theta2_ / 1000.0 * miu ** prs.theta2_
    Y = Yg * (1 - Dfrac) - Cab
    I = s * Y; Cons = Y - I
    ECo2 = (sig * Yg + eland) * (1 - miu)
    Vc = min(prs.sd_cons_ ** 2 * (y - 2020), prs.sd_cons_ ** 2 * 235)
    STP = 1.0 / (1 + cfg['rartp']) ** (y - 2020) / (1 + -0.5 * Vc * prs.eta_ ** 2) ** (y - 2020)
    c = Cons * 1000.0 / L
    U = (c ** (1 - prs.eta_) - 1) / (1 - prs.eta_) - 1
    pi = len(S) - 5
    if cfg.get('perma'):
        # temperature tendency needed for gradual thaw: compute forcing first
        pass
    # carbon cycle and forcing
    if cfg['carbon'] == 'dfair':
        R = S[NE:NE + 4]; Cum = S[NE + 4]
        MAT = SC.Meq + R.sum()
    else:
        Ca = S[NE]; MAT = Ca
    Efn = cfg.get('E_fn')                       # optional prescribed CO2 emissions (GtC/yr), e.g. zero-emissions test
    Fnon = (cfg['F_nonco2'](y) if callable(cfg['F_nonco2']) else cfg['F_nonco2']) if 'F_nonco2' in cfg else Fmisc + FGHGexo
    F = prs.fco22x_ * math.log(MAT / prs.Meq_) / math.log(2) + Fnon
    dTb1 = (prs.Teq1_ * F - Tb1) / prs.d1_; dTb2 = (prs.Teq2_ * F - Tb2) / prs.d2_
    if cfg.get('perma'):
        dperm, Fmelt = perma_flows(T, dTb1 + dTb2, S[pi], S[pi + 1], S[pi + 2], S[pi + 3], T200, cfg['perma'])
        if not cfg.get('fb', {}).get('perma', 1):                 # factorial: permafrost release to the atmosphere off
            dperm = dperm.copy(); dperm[1] += Fmelt; dperm[4] = 0.0; Fmelt = 0.0
    else:
        dperm, Fmelt = np.zeros(5), 0.0
    if cfg['carbon'] == 'dfair':
        E = (sum(np.atleast_1d(Efn(y))) if Efn else ECo2 / 3.667) + Fmelt + pulse_gtc
        a = alpha_of(Cum - (MAT - SC.Meq), T)
        dcarb = np.concatenate([np.array(ESH) * E - R / (np.array(TAU) * a), [E]])
    else:
        Ca, Cv, Cs, Dc = S[NE:NE + 4]; X = S[NE + 4:NE + 16]
        e_ff = sig * Yg * (1 - miu) * 0.2727; e_luc = eland * (1 - miu) * 0.2727
        if Efn:                              # prescribed emissions: GtC/yr, or (fossil, land-use) GtC/yr
            _e = np.atleast_1d(Efn(y)); e_ff, e_luc = float(_e[0]), (float(_e[1]) if len(_e) > 1 else 0.0)
        if 'fbx' in cfg or 'KM' in cfg:      # feedback switches referenced to pre-industrial, and/or K_M override
            Ph, Rp, Tu, Rs, dX, asg = lenton_fluxes(Ca, Cv, Cs, Dc, X, T, cfg.get('fbx', {}), cfg.get('KM', KM))
            dcarb = np.concatenate([[-Ph + Rp + Rs - asg + e_ff + e_luc + Fmelt + pulse_gtc, Ph - Rp - Tu - e_luc, Tu - Rs, e_luc], dX])
            deco = np.array([dTb1, dTb2, I - prs.dK_ * K, ECo2 / 3.666, U * L * STP * prs.scale1_,
                             prs.Fc1_ * CO2EGHGabateB * (1 - miu) + math.log(prs.Fc2_) * FGHG])
            return np.concatenate([deco, dcarb, dperm]), (c, STP, Y, Yg, Cab, Dfrac, MAT)
        fb = cfg.get('fb', {})          # feedback switches (factorial); default all on. 'Off' freezes the dependence at 2020:
        T20 = 0.1477 + 1.099454         # temperature responses use the 2020 anomaly; fertilisation uses 2020 atmospheric CO2;
        dX, asg, fl, _ = LF.ocean_rhs(Ca, X, T if OW else 0.0)   # (factorial below: surface-ocean pCO2 held at 2020 values)
        if not fb.get('ocean', 1):
            pa = Ca * P['fgtm'] / P['ka']
            for j, b in ((0, 'l'), (1, 'h')):
                Fj = fl[b + '_in'] / pa * (pa - PS20[b])
                dX = dX.copy(); dX[j] += Fj - fl[b]; fl = dict(fl); fl[b] = Fj
            asg = (fl['l'] + fl['h']) * 12e-15
        Ph, _, _, _ = LR.land(Ca if fb.get('fert', 1) else L20['Ca'], Cv, Cs, T if fb.get('photoT', 1) else T20, Dc, KM)
        _, Rp, Tu, _ = LR.land(Ca, Cv, Cs, T if fb.get('vegT', 1) else T20, Dc, KM)
        _, _, _, Rs = LR.land(Ca, Cv, Cs, T if fb.get('soilT', 1) else T20, Dc, KM)
        dcarb = np.concatenate([[-Ph + Rp + Rs - asg + e_ff + e_luc + Fmelt + pulse_gtc, Ph - Rp - Tu - e_luc, Tu - Rs, e_luc], dX])
    deco = np.array([dTb1, dTb2, I - prs.dK_ * K, ECo2 / 3.666, U * L * STP * prs.scale1_,
                     prs.Fc1_ * CO2EGHGabateB * (1 - miu) + math.log(prs.Fc2_) * FGHG])
    return np.concatenate([deco, dcarb, dperm]), (c, STP, Y, Yg, Cab, Dfrac, MAT)


def run(cfg, dt=0.05, pulse=None, nyears=400):
    """pulse = (year, GtCO2) released evenly through that year."""
    S = init(cfg); n = int(round(1 / dt)); lag = int(round(200 / dt))
    buf = [phf.WF * phf.T_hist(1820 + j * dt) for j in range(lag)]; csum = float(sum(buf)); t2m = 0.0
    out, cons = [], []
    pi = len(S) - 5
    for kk in range(nyears * n + 1):
        y = 2020 + kk * dt
        if kk % n == 0:
            _, diag = rhs(S, y, 0.0, cfg)
            out.append(S.copy()); cons.append(diag)
        if kk == nyears * n:
            break
        if (cfg.get('perma') or {}).get('q10') == 'running':   # sensitivity: Q10 on the running mean (not sticky)
            T200 = csum / lag
        else:
            T200 = max(t2m, csum / lag); t2m = T200
        pg = (pulse[1] / 3.667) if (pulse and pulse[0] <= y < pulse[0] + 1) else 0.0
        k1 = rhs(S, y, T200, cfg, pg)[0]; k2 = rhs(S + dt / 2 * k1, y + dt / 2, T200, cfg, pg)[0]
        k3 = rhs(S + dt / 2 * k2, y + dt / 2, T200, cfg, pg)[0]; k4 = rhs(S + dt * k3, y + dt, T200, cfg, pg)[0]
        S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        new = phf.WF * (S[0] + S[1]); csum += new - buf[kk]; buf.append(new)
    A = np.array(out)
    C = np.array(cons)
    d = dict(TATM=A[:, 0] + A[:, 1], Uacc=A[:, 4], release=A[:, pi + 4], cons=C[:, :2],
             thaw=A[:, pi + 2] + A[:, pi + 3] - A[0, pi + 2] - A[0, pi + 3], thawed=A[:, pi + 1], Cperm=A[:, pi], Ctot=A[:, 3], Y=C[:, 2], Ygross=C[:, 3], Cabate=C[:, 4], Dfrac=C[:, 5], MAT=C[:, 6])
    d['CO2ppm'] = (SC.Meq + A[:, NE:NE + 4].sum(axis=1)) * LR.PPM if cfg['carbon'] == 'dfair' else A[:, NE] * LR.PPM
    if cfg['carbon'] == 'lenton':                     # carbon stocks (PgC) for flux attribution
        d['land'] = A[:, NE + 1] + A[:, NE + 2]; d['ocean'] = A[:, NE + 4:NE + 8].sum(axis=1) * 12e-15; d['atm'] = A[:, NE]
    return d


def scc(cfg, year=2020, dE=1.0, base=None):
    """$ per tCO2 (2019$): -1000 * dW/dE / dW/dC, W = U_acc(2420); dW/dC = 1000*scale1*STP*c^-eta over the year."""
    b = base or run(cfg)
    p = run(cfg, pulse=(year, dE))
    dW = (p['Uacc'][-1] - b['Uacc'][-1]) / dE
    i = year - 2020
    mu = [1000 * prs.scale1_ * b['cons'][j][1] * b['cons'][j][0] ** (-prs.eta_) for j in (i, i + 1)]
    return -1000 * dW / (0.5 * (mu[0] + mu[1])), b

"""
Faithful Python reimplementation of the STELLA DICE-C / DICE-CP system.
Purpose: VALIDATION ONLY -- confirm the equation extraction reproduces
FINALresults.xlsx before translating to Mathematica.

Time: continuous, RK4, DT=0.01, 2020..2420 (matches STELLA sim_specs).
Exogenous drivers (pop, TFP, sigma, eland, non-CO2) reconstructed from the
DICE2023 / Nonco2 GAMS formulas on the 5-year period grid, then linearly
interpolated to annual -- this is the defensible, hard-coded approach.

Three configurations via switches:
  DICE 2023 : Lenton?=0 (DFAIR carbon), permafrost off
  DICE-C    : Lenton?=1, permafrost off
  DICE-CP   : Lenton?=1, permafrost on
"""
import numpy as np
import math

# ============================================================
# EXOGENOUS DRIVERS from GAMS formulas (5-year period grid)
# period index p = 1..81  <->  year = 2020 + 5*(p-1)
# ============================================================

def build_exogenous():
    # --- DICE2023-b-4-3-10 parameters ---
    gama=0.300; pop1=7752.9; popadj=0.145; popasym=10825.0
    q1=135.7; AL1=5.84; gA1=0.066; delA=0.0015
    gsigma1=-0.015; delgsig=0.96; asymgsig=-0.005
    e1=37.56; miu1=0.05
    tstep=5
    # --- Nonco2 parameters ---
    eland0_nc=5.9; deland_nc=0.1
    F_Misc2020=-0.054; F_Misc2100=0.265
    ECO2eGHGB2020=9.96; ECO2eGHGB2100=15.5
    emissrat2020=1.40; emissrat2100=1.21

    P=81
    pv=np.arange(1,P+1)   # t.val

    # population: L("1")=pop1; loop L(t+1)=L(t)*(popasym/L(t))**popadj
    L=np.zeros(P); L[0]=pop1
    for i in range(P-1):
        L[i+1]=L[i]*(popasym/L[i])**popadj

    # TFP: gA(t)=gA1*exp(-delA*5*(t.val-1)); aL("1")=AL1; aL(t+1)=aL(t)/(1-gA(t))
    gA=gA1*np.exp(-delA*5*(pv-1))
    aL=np.zeros(P); aL[0]=AL1
    for i in range(P-1):
        aL[i+1]=aL[i]/(1-gA[i])

    # sigma: sig1=e1/(q1*(1-miu1)); gsig(t)=min(gsigma1*delgsig**(t.val-1),asymgsig); sigma(t+1)=sigma(t)*exp(5*gsig(t))
    sig1=e1/(q1*(1-miu1))
    gsig=np.minimum(gsigma1*delgsig**(pv-1), asymgsig)
    sigma=np.zeros(P); sigma[0]=sig1
    for i in range(P-1):
        sigma[i+1]=sigma[i]*np.exp(5*gsig[i])

    # eland (Nonco2): eland0*(1-deland)**(t.val-1)   [GtCO2/yr]
    eland=eland0_nc*(1-deland_nc)**(pv-1)

    # non-CO2 piecewise-linear-to-2100-then-flat (period 16 = year 2095; 17=2100)
    def pw(v2020,v2100):
        out=np.where(pv<=16, v2020+((v2100-v2020)/16)*(pv-1), v2100)
        return out
    F_Misc=pw(F_Misc2020,F_Misc2100)
    CO2E_GHGabateB=pw(ECO2eGHGB2020,ECO2eGHGB2100)
    emissrat=pw(emissrat2020,emissrat2100)
    sigmatot=sigma*emissrat

    years5=2020+5*(pv-1)
    return dict(years5=years5, L=L, aL=aL, sigma=sigma, eland=eland,
                F_Misc=F_Misc, CO2E_GHGabateB=CO2E_GHGabateB,
                emissrat=emissrat, sigmatot=sigmatot)

EXO=build_exogenous()

def interp_annual(name, t):
    return np.interp(t, EXO['years5'], EXO[name])

# ============================================================
# STELLA parameters (from DICE_2023_5_12.stmx)
# ============================================================
P_ = dict(
    # Climate
    Teq1=0.324, Teq2=0.44, d1=236.0, d2=4.07, fco22x=3.93, Meq=588.0,
    # Economy
    dK=0.0819, gamma=0.3, a1=0.0, a3=2.0, theta2=2.6,
    eta=0.95, pi_=0.05, lam_land=0.02, Eland0=5.9,
    remiss0=1.4, remiss80=1.21, pback2050=515.0,
    scale1=0.00891061, scale2=-6275.91, sd_cons=0.01,
    Fc1=0.00955, Fc2=0.861, sigma_scale=1.0,
    ECO2eGHGB2020=9.96, ECO2eGHGB2100=15.5,
    # DFAIR
    eps0=0.2173,eps1=0.224,eps2=0.2824,eps3=0.2763,
    tau0=1e6,tau1=394.4,tau2=36.53,tau3=4.304,
    IRF0=32.4,irC=0.019,irT=4.165,
    # Lenton
    kp=0.184,kr=0.092,kt=0.092,ksr=0.0337,KM=120.0,Ea=54830.0,
    kMM=1.478,kA=8.7039e9,kT=6.5e7,kU=4.87e7,kwi=1.25e7,kid=2e8,
    ks=1752.0,kc=29.0,kB=157.072,kD=0.27,R=8.314,Th=275.15,Tl=295.15,
    Al=0.85,Ah=0.15,aO=3.6e14,beta=11.0,pCO2w0=286e-6,pCO2c0=270e-6,
    rho=1027.0,Vl=3.06e16,Vh=5.4e15,Vi=3.6e17,Vd=9.7e17,
    eland0_L=2.6,deland_L=0.022009596,ka=1.773e20,fgtm=8.3259e13,
    Ca0=887.0,Cv0=550.0,Cs0=1500.0,Col0=730.0,Coh0=140.0,Coi0=10040.0,Cod0=23860.0,
    # Permafrost
    Q10=2.0,F_CH4=0.023,fstatic=0.74,
)


from scipy.optimize import brentq
def solve_alpha(IRFtarget):
    def f(a):
        s=0.0
        for eps,tau in [(P_["eps0"],P_["tau0"]),(P_["eps1"],P_["tau1"]),(P_["eps2"],P_["tau2"]),(P_["eps3"],P_["tau3"])]:
            s+=a*eps*tau*(1-np.exp(-100/(a*tau)))
        return s-IRFtarget
    try:
        return brentq(f,0.01,100.0,xtol=1e-10)
    except Exception:
        return 0.45


import os as _os
_HERE=_os.path.dirname(_os.path.abspath(__file__))
import json as _json
_lut=_json.load(open(_os.path.join(_HERE,"alpha_lut.json")))
_lut_x=np.linspace(_lut["xmin"],_lut["xmax"],len(_lut["ypts"]))
_lut_y=np.array(_lut["ypts"])
def alpha_lut(irf):
    return np.interp(irf,_lut_x,_lut_y)

_em=_json.load(open(_os.path.join(_HERE,"exomat.json")))
_em_x=np.linspace(_em["xmin"],_em["xmax"],len(_em["ypts"]))
_em_y=np.array(_em["ypts"])
def exomat(y): return np.interp(y,_em_x,_em_y)
_fg=_json.load(open(_os.path.join(_HERE,"fghgexo.json")))
_fg_x=np.linspace(_fg["xmin"],_fg["xmax"],len(_fg["ypts"]))
_fg_y=np.array(_fg["ypts"])
def fghgexo(y): return np.interp(y,_fg_x,_fg_y)


# DFAIR alpha lookup (IRF100 -> alpha), xscale 10..95
_alpha_x=np.linspace(10,95,13)
_alpha_y=np.array([0.0000549,0.0000563,0.0000578,0.0000594,0.0000609,0.0000626,
                   0.0000642,0.000066,0.0000677,0.0000696,0.0000715,0.0000734,0.0])
# NOTE: last y=0 is a STELLA artifact; clamp usage below.

# Control graphical functions (optimal MIU and S paths), xscale 2020..2420
_miu_x=np.linspace(2020,2420,13)
_miu_y=np.array([0.05,0.1,0.24,0.27468931,0.30992069,0.34890274,0.39194427,
                 0.42714712,0.46449982,0.50401703,0.54570128,0.58954338,0.63552])
_s_x=np.linspace(2020,2420,13)
_s_y=np.array([0.25636878,0.25373313,0.25114651,0.24956814,0.24842159,0.24760613,
               0.247041,0.24678246,0.24676827,0.24698039,0.24740432,0.24740432,0.24740432])
# NOTE: last two _s_y padded (STELLA table truncated at 0); using plateau value.
import json as _cj
_ctrl=_cj.load(open(_os.path.join(_HERE,"controls.json")))
_miu_xx=np.linspace(_ctrl["miu"]["xmin"],_ctrl["miu"]["xmax"],len(_ctrl["miu"]["ypts"]))
_miu_yy=np.array(_ctrl["miu"]["ypts"])
_s_xx=np.linspace(_ctrl["s"]["xmin"],_ctrl["s"]["xmax"],len(_ctrl["s"]["ypts"]))
_s_yy=np.array(_ctrl["s"]["ypts"])
def miu_of(t): return np.interp(t,_miu_xx,_miu_yy)
def s_of(t): return np.interp(t,_s_xx,_s_yy)



print("Exogenous drivers built. Sample at 2020/2100/2300:")
for nm in ['L','aL','sigma','eland','F_Misc','CO2E_GHGabateB']:
    vals=[interp_annual(nm,y) for y in (2020,2100,2300)]
    print(f"  {nm:16s}: {vals[0]:.4f}  {vals[1]:.4f}  {vals[2]:.4f}")

# ============================================================
# DYNAMICAL CORE
# State vector depends on config. We carry a superset and only
# evolve the active carbon module.
#
# Common states: Tbox1, Tbox2, K, Ctot, Uacc, FGHG
# DFAIR states:  R0,R1,R2,R3
# Lenton states: Ca,Cv,Cs,Col,Coh,Coi,Cod,Ldef
# Permafrost:    Cperm, thawedC   (+ history buffers for DELAY)
# ============================================================

def frozen_fraction_of(TATM):
    # Eq. 15: frozen_fraction = NORMALCDF(LN(1+T_HL), 1e6, 1.67, 0.99), T_HL=2*TATM
    T_HL=2.0*TATM
    if 1.0+T_HL<=0: return 0.0
    z=(np.log(1.0+T_HL)-1.67)/0.99
    return 1.0-0.5*(1.0+math.erf(z/np.sqrt(2.0)))   # upper tail

def make_deriv(config):
    lenton = (config in ('DICE-C','DICE-CP'))
    perma  = (config == 'DICE-CP')
    p=P_

    # index layout
    # 0 Tbox1,1 Tbox2,2 K,3 Ctot,4 Uacc,5 FGHG
    # 6..9  R0..R3           (DFAIR)
    # 10 Ca,11 Cv,12 Cs,13 Col,14 Coh,15 Coi,16 Cod,17 Ldef  (Lenton)
    # 18 Cperm,19 thawedC    (permafrost)

    def econ_and_emissions(y, K, TATM, FGHG):
        L=interp_annual('L',y); aL=interp_annual('aL',y); sig=interp_annual('sigma',y)*p['sigma_scale']
        eland=interp_annual('eland',y)  # GtCO2/yr (Nonco2 eland0=5.9,deland=0.1)
        miu=miu_of(y); s=s_of(y)
        Ygross=aL*(L/1000.0)**(1-p['gamma'])*(K**p['gamma'])
        Dfrac=p['a1']*TATM + (0.003467)*TATM**p['a3']    # altdam=0
        # backstop
        if y-2020<=30: pback=p['pback2050']*np.exp(-0.01*(y-2050))
        else:          pback=p['pback2050']*np.exp(-0.001*(y-2050))
        remiss = p['remiss0']+((p['remiss80']-p['remiss0'])/80)*(y-2020) if y<=2100 else p['remiss80']
        sigmatot=sig*remiss
        theta1=(sigmatot*pback)/p['theta2']/1000.0
        Cabate=Ygross*theta1*miu**p['theta2']
        Ynet=Ygross*(1-Dfrac)
        Y=Ynet-Cabate
        I=s*Y; Cons=Y-I
        # emissions
        ECo2=(sig*Ygross+eland)*(1-miu)   # GtCO2/yr
        EIND=(sig*Ygross)*(1-miu)
        return dict(L=L,Ygross=Ygross,Dfrac=Dfrac,Y=Y,I=I,Cons=Cons,ECo2=ECo2,EIND=EIND,miu=miu)

    def deriv(y, S, hist):
        Tbox1,Tbox2,K,Ctot,Uacc,FGHG = S[0],S[1],S[2],S[3],S[4],S[5]
        TATM=Tbox1+Tbox2

        e=econ_and_emissions(y,K,TATM,FGHG)
        ECo2=e['ECo2']; EIND=e['EIND']; miu=e['miu']

        # ---- carbon module -> atmospheric_carbon (GtC) ----
        if lenton:
            Ca,Cv,Cs,Col,Coh,Coi,Cod,Ldef = S[10:18]
            atmC=Ca
        else:
            R0,R1,R2,R3=S[6:10]
            atmC=exomat(y)   # DICE2023 (Lenton?=0): atmospheric carbon = exo MAT playback

        # ---- climate forcing ----
        Fmisc = p['Fc1']*0+ ( (-0.054)+((0.265-(-0.054))/(2100-2020))*(y-2020) if y<2100 else 0.265 )
        FGHGexo=fghgexo(y)   # exogenous non-CO2 forcing lookup
        F=p['fco22x']*(np.log(atmC/p['Meq'])/np.log(2))+Fmisc+FGHGexo
        dTbox1=(p['Teq1']*F-Tbox1)/p['d1']
        dTbox2=(p['Teq2']*F-Tbox2)/p['d2']

        # ---- FGHG stock (Flow 7) ----
        CO2E_GHGabateB=interp_annual('CO2E_GHGabateB',y)
        dFGHG=p['Fc1']*CO2E_GHGabateB*(1-miu)-(-np.log(p['Fc2']))*FGHG

        # ---- capital, Ctot, utility ----
        dK=e['I']-p['dK']*K
        dCtot=ECo2/3.666
        # utility
        Vc=min(p['sd_cons']**2*(y-2020), p['sd_cons']**2*235)
        rpre=-0.5*Vc*p['eta']**2
        rho=0.001; betaclim=0.5   # altdisc=0
        rartp=np.exp(rho+betaclim*p['pi_'])-1
        # GAMS: RR(t) = RR1(t)*(1+rprecaut(t))**(-tstep*(t.val-1)) -- negative
        # exponent, precautionary term offsets base decay. Confirmed against
        # raw .gms source and GAMSresults.csv (actual Nordhaus output).
        STP=1/(1+rartp)**(y-2020)/(1+rpre)**(y-2020)
        c=e['Cons']*1000/e['L']
        U=(c**(1-p['eta'])-1)/(1-p['eta'])-1
        dUacc=U*e['L']*STP*p['scale1']

        dS=np.zeros_like(S)
        dS[0]=dTbox1; dS[1]=dTbox2; dS[2]=dK; dS[3]=dCtot; dS[4]=dUacc; dS[5]=dFGHG

        # ---- carbon dynamics ----
        if not lenton:
            # DFAIR
            Cacc=Ctot-(atmC-p['Meq'])
            IRF100=p['IRF0']+p['irC']*Cacc+p['irT']*TATM
            alpha=alpha_lut(IRF100)
            for j,(eps,tau) in enumerate([(p['eps0'],p['tau0']),(p['eps1'],p['tau1']),
                                          (p['eps2'],p['tau2']),(p['eps3'],p['tau3'])]):
                Rj=S[6+j]
                dS[6+j]=(eps*ECo2/3.667)-((1-np.exp(-5/(tau*alpha)))/5)*Rj
        else:
            Ca,Cv,Cs,Col,Coh,Coi,Cod,Ldef=S[10:18]
            anomaly1=TATM
            surfT=288+anomaly1
            anomaly=anomaly1  # photo_temp_switch=1
            Cvs=p['Cv0']-Ldef*p['kD']
            pCO2=Ca/p['ka']*p['fgtm']
            pCO2ppmv=pCO2*1e6  # Terr switch=1
            # land fluxes
            Photo=(p['kp']*Cvs*p['kMM'])*((pCO2ppmv-p['kc'])/(p['KM']+(pCO2ppmv-p['kc'])))*((15+anomaly)**2*(25-anomaly)/5625)
            PlantR=p['kr']*Cv*p['kA']*np.exp(-p['Ea']/(p['R']*surfT))
            Turn=p['kt']*Cv
            SoilR=p['ksr']*Cs*p['kB']*2.7**(-1*(308.56/(surfT-227.13)))
            # ocean (WITH 12e-15 factor -- ORIGINAL/CORRECT STELLA form)
            bw=np.exp(-60.2409+93.4517*(100/p['Tl'])+23.3585*np.log(p['Tl']/100)+34.7*(0.023517-0.023656*(p['Tl']/100)+0.0047036*(p['Tl']/100)**2))
            bc=np.exp(-60.2409+93.4517*(100/p['Th'])+23.3585*np.log(p['Th']/100)+34.7*(0.023517-0.023656*(p['Th']/100)+0.0047036*(p['Th']/100)**2))
            pCO2w=p['pCO2w0']*(Col/p['Col0'])**p['beta']
            pCO2c=p['pCO2c0']*(Coh/p['Coh0'])**p['beta']
            Faw=p['ks']*p['aO']*p['Al']*bw*(pCO2-pCO2w)*12e-15*p['rho']
            Fac=p['ks']*p['aO']*p['Ah']*bc*(pCO2-pCO2c)*12e-15*p['rho']
            AirSeaLL=-Faw; AirSeaHL=-Fac
            WarmInt=p['kwi']*(Col/p['Vl'])*3.15576e7*12e-15
            IntDeep=p['kid']*(Coi/p['Vi'])*3.15576e7*12e-15
            ColdInt=p['kU']*(Coh/p['Vh'])*3.15576e7*12e-15
            DeepInt=p['kid']*(Cod/p['Vd'])*3.15576e7*12e-15
            IntWarm=(p['kT']+p['kwi'])*(Coi/p['Vi'])*3.15576e7*12e-15
            elandL=p['eland0_L']*(1-p['deland_L'])**(y-2015)
            Defor=elandL/2.12
            GHGemis=EIND*0.2727
            # Land-use change emissions also flow into the SCM atmosphere (paper 2.3:
            # "Anthropogenic CO2 emissions, comprising both industrial emissions and
            # land-use change emissions ... flow into the SCM atmospheric carbon pool").
            # STELLA's GHG-emissions ghost (EIND) carries total anthropogenic CO2; the
            # land-use piece did not serialize in the XML. eland is GtCO2/yr; 0.2727=12/44.
            elandCO2=interp_annual('eland',y)*(1-miu)
            GHGemis=GHGemis+elandCO2*0.2727

            # ---- permafrost module (Eqs. 15-18), DICE-CP only ----
            Fmelt=0.0
            dCperm=0.0; dthawed=0.0
            if perma:
                Cperm=S[18]; thawedC=S[19]
                # Eq. 15 (STELLA): frozen_fraction = NORMALCDF(LN(1+T_HL), 1e6, 1.67, 0.99)
                #   = upper-tail prob = 1 - Phi((ln(1+T_HL)-1.67)/0.99),  T_HL = 2*T_atm
                #   -> frozen fraction DECREASES as warming rises (correct physics)
                f_frozen=frozen_fraction_of(TATM)
                # Eq. 16 (STELLA source, DICE 2023 5_12.stmx line 2247):
                #   delta_Cperm = (DELAY(frozen_fraction,1) - frozen_fraction) * Cperm
                #   DELAY(x,1) is a literal 1-year pipeline delay (x(t-1)), not an
                #   instantaneous derivative -- STELLA's DELAY() ~= transport delay.
                #   NOT clamped in STELLA -- only the Cperm/thawed_C *stocks* carry <non_negative/>,
                #   not this flow. So gradual thaw goes negative (re-freeze) when climate cools.
                ff_delay1=hist.get('ff_delay1', f_frozen)
                gradual_thaw=(ff_delay1-f_frozen)*Cperm
                # Eq. 17 (STELLA line 2296): abrupt_thaw = IF Anomaly>2 THEN delta_Cperm*0.4 ELSE 0
                #   this flow DOES carry <non_negative/> in STELLA -- only the abrupt piece is
                #   barred from going negative, gradual thaw is not.
                abrupt=max(0.0, 0.4*gradual_thaw) if TATM>2.0 else 0.0
                total_thaw=gradual_thaw+abrupt        # Cperm -> thawedC (>=0)
                # Eq. 18: Fthawedatmo = labile*Q10^(T200/10), labile=(1-fstatic)*thawedC
                Clabile=(1.0-p['fstatic'])*thawedC
                T200=hist.get('T200',1.25)
                Jthaw_atm=Clabile*p['Q10']**(T200/10.0)
                Fmelt=Jthaw_atm                       # RH_pf+RH_CH4 both to atmosphere
                dCperm=-total_thaw                    # Cperm outflows: delta_Cperm+abrupt
                dthawed=total_thaw-Jthaw_atm          # thawedC in: thaw; out: Fthawedatmo

            dCa=PlantR+SoilR+AirSeaLL+AirSeaHL+GHGemis+Fmelt-Photo
            dCv=Photo-PlantR-Turn
            dCs=Turn-SoilR
            dCol=IntWarm+Faw-WarmInt          # inflow intwarm; out: AirSea(LL)=-Faw so +Faw; -WarmInt
            dCoh=Fac-ColdInt
            dCoi=WarmInt+ColdInt+DeepInt-IntDeep-IntWarm
            dCod=IntDeep-DeepInt
            dLdef=Defor
            dS[10:18]=[dCa,dCv,dCs,dCol,dCoh,dCoi,dCod,dLdef]
            if perma:
                dS[18]=dCperm; dS[19]=dthawed
        return dS, TATM, atmC

    return deriv, lenton, perma

print("\nDeriv builder ready.")

# ============================================================
# RK4 STEPPER (no-permafrost configs first)
# ============================================================
def simulate(config, t0=2020, t1=2420, dt=0.01):
    deriv, lenton, perma = make_deriv(config)
    N=int(round((t1-t0)/dt))
    S=np.zeros(20)
    # init common
    S[0]=0.1477; S[1]=1.099454; S[2]=295.0; S[3]=633.5; S[4]=0.0; S[5]=0.518
    # init DFAIR
    S[6]=150.093; S[7]=102.698; S[8]=39.534; S[9]=6.1865
    # init Lenton
    S[10]=P_['Ca0']; S[11]=P_['Cv0']; S[12]=P_['Cs0']; S[13]=P_['Col0']
    S[14]=P_['Coh0']; S[15]=P_['Coi0']; S[16]=P_['Cod0']; S[17]=0.0
    # init permafrost
    S[18]=1035.0; S[19]=0.0

    ff0=frozen_fraction_of(S[0]+S[1])  # frozen_fraction at IC; DELAY(...,1) with no
                                        # explicit initial value defaults to this
    hist=dict(T200=1.099454+0.1477, ff_delay1=ff0)
    # T200: 200-yr delayed TATM. Ring buffer of past TATM at each step.
    lag_steps=int(round(200.0/dt))
    tatm_buf=[]  # append current TATM each step; T200 = value lag_steps ago
    # frozen_fraction: 1-yr delayed. Ring buffer of past frozen_fraction at each step.
    lag_steps_ff=int(round(1.0/dt))
    ff_buf=[]
    out_years=[]; out_TATM=[]; out_atmC=[]
    rec_every=int(round(1.0/dt))  # record annually
    for k in range(N+1):
        y=t0+k*dt
        # current TATM for buffer + T200 lookup
        TATM_now=S[0]+S[1]
        tatm_buf.append(TATM_now)
        if y>=t0+200:
            # DELAY(Anomaly, 200, 1.25): value 200 yr ago
            hist['T200']=tatm_buf[k-lag_steps]
        else:
            # STELLA analytic ramp: (1.25/(exp(0.01*200)-1))*(exp(0.01*(t-2020))-1)
            hist['T200']=(1.25/(np.exp(0.01*200)-1.0))*(np.exp(0.01*(y-2020))-1.0)
        if perma:
            ff_buf.append(frozen_fraction_of(TATM_now))
            # DELAY(frozen_fraction, 1): value 1 yr ago, else defaults to initial value
            hist['ff_delay1']=ff_buf[k-lag_steps_ff] if y>=t0+1 else ff0
        if k%rec_every==0:
            _,TATM,atmC=deriv(y,S,hist)
            out_years.append(y); out_TATM.append(TATM); out_atmC.append(atmC)
        if k==N: break
        # RK4 (permafrost delay quantities held fixed across the 4 substeps)
        k1,_,_=deriv(y,S,hist)
        k2,_,_=deriv(y+dt/2,S+dt/2*k1,hist)
        k3,_,_=deriv(y+dt/2,S+dt/2*k2,hist)
        k4,_,_=deriv(y+dt,S+dt*k3,hist)
        S=S+dt/6*(k1+2*k2+2*k3+k4)
        if perma:
            # STELLA <non_negative/> on Cperm and thawed_C stocks (not on the delta_Cperm
            # flow itself -- gradual thaw can go negative, but the stocks it feeds can't).
            # MASS-CONSERVING non-negative stock (fixed 2026-09-25): real STELLA clips the
            # (negative) delta_Cperm flow to what thawed_C can supply, so re-freezing never
            # exceeds the pile and Cperm cannot grow from nothing (verified against a STELLA
            # output: Cperm stays flat once thawed_C ~ 0 while the raw flow would be ~ -0.15/yr).
            # A bare "clamp thawed_C to 0" left the excess refreeze in Cperm (created carbon).
            # Post-step equivalent: if thawed_C went negative by d, only (refreeze - d) was
            # available, so remove d from Cperm as well.
            if S[19]<0.0:
                S[18]+=S[19]   # S[19] < 0: take the un-supplied refreeze back out of Cperm
                S[19]=0.0
            if S[18]<0.0: S[18]=0.0
    return np.array(out_years),np.array(out_TATM),np.array(out_atmC)

# ---- validate DICE2023 and DICE-C ----
import openpyxl
wb=openpyxl.load_workbook(_os.path.join(_HERE,'FINALresults.xlsx'),data_only=True,read_only=True)
ws=wb['Results']; rows=list(ws.iter_rows(values_only=True)); H=rows[0]; Dr=rows[1:]
def ci(run,var): return H.index(f'Run {run}: {var}')
xlyears=[r[0] for r in Dr]
def ref(run,var,y):
    i=xlyears.index(y); return Dr[i][ci(run,var)]

for config,run in [('DICE2023',1),('DICE-C',9),('DICE-CP',13)]:
    yr,tatm,atmc=simulate(config)
    print(f"\n=== {config} (Run {run}) ===")
    print(f"{'year':>6} {'TATM_sim':>10} {'TATM_ref':>10} {'atmC_sim':>10} {'atmC_ref':>10}")
    for y in [2020,2021,2025,2050,2100,2200,2300]:
        i=list(yr).index(y)
        rt=ref(run,'Climate.Atmospheric temperature anomaly (TATM)',y)
        ra=ref(run,'Climate.atmospheric carbon',y)
        print(f"{y:>6} {tatm[i]:>10.4f} {rt:>10.4f} {atmc[i]:>10.2f} {ra:>10.2f}")

# ============================================================
# DIAGNOSTIC: first-step budget for DICE2023
# ============================================================
print("\n\n===== DIAGNOSTIC: DICE2023 first-step =====")
deriv,lenton,perma=make_deriv('DICE2023')
S=np.zeros(20)
S[0]=0.1477;S[1]=1.099454;S[2]=295.0;S[3]=633.5;S[4]=0.0;S[5]=0.518
S[6]=150.093;S[7]=102.698;S[8]=39.534;S[9]=6.1865
hist=dict(Fthawedatmo=0.0,dCperm=0.0,dthawed=0.0)
# recompute pieces at y=2020
y=2020.0
L=interp_annual('L',y);aL=interp_annual('aL',y);sig=interp_annual('sigma',y)
eland=interp_annual('eland',y);miu=miu_of(y);s=s_of(y)
K=S[2];TATM=S[0]+S[1]
Ygross=aL*(L/1000)**(0.7)*(K**0.3)
print(f"Ygross={Ygross:.4f} (ref YGROSS 2020 ~ {ref if False else '134.16 from GAMS'})")
ECo2=(sig*Ygross+eland)*(1-miu)
print(f"sigma={sig:.5f} eland={eland:.4f} miu={miu:.4f} s={s:.5f}")
print(f"ECo2={ECo2:.4f} GtCO2/yr  (=> {ECo2/3.667:.4f} GtC/yr into reservoirs)")
Ctot=S[3];atmC=P_['Meq']+S[6]+S[7]+S[8]+S[9]
Cacc=Ctot-(atmC-P_['Meq'])
IRF100=P_['IRF0']+P_['irC']*Cacc+P_['irT']*TATM
print(f"atmC={atmC:.4f} Cacc={Cacc:.4f} IRF100={IRF100:.4f}")
alpha=np.interp(IRF100,_alpha_x[:-1],_alpha_y[:-1])
print(f"alpha={alpha:.8e}  (alpha_x range {_alpha_x[0]}..{_alpha_x[-2]})")
# reservoir 0 deriv
eps,tau=P_['eps0'],P_['tau0']
term_in=eps*ECo2/3.667
term_out=((1-np.exp(-5/(tau*alpha)))/5)*S[6]
print(f"R0: in={term_in:.5f} out={term_out:.5f} dR0={term_in-term_out:.5f}")
eps,tau=P_['eps3'],P_['tau3']
term_in=eps*ECo2/3.667
term_out=((1-np.exp(-5/(tau*alpha)))/5)*S[9]
print(f"R3: in={term_in:.5f} out={term_out:.5f} dR3={term_in-term_out:.5f}")

print("\n===== Solve DFAIR alpha implicitly (GAMS irfeqlhs=irfeqrhs) =====")
from scipy.optimize import brentq
def IRF_of_alpha(alpha):
    s=0
    for eps,tau in [(0.2173,1e6),(0.224,394.4),(0.2824,36.53),(0.2763,4.304)]:
        s+=alpha*eps*tau*(1-np.exp(-100/(alpha*tau)))
    return s
# target IRF100 at 2020 = irf0 + irC*Cacc + irT*TATM
Cacc=334.9885; TATM=1.247154
IRFtarget=32.4+0.019*Cacc+4.165*TATM
print(f"IRFtarget={IRFtarget:.4f}")
alpha_solved=brentq(lambda a: IRF_of_alpha(a)-IRFtarget, 0.01, 100)
print(f"alpha_solved={alpha_solved:.6f}")
# check reservoir decay with this alpha
tau0=1e6
print(f"exp(-5/(tau0*alpha))={np.exp(-5/(tau0*alpha_solved)):.6f}")

# ============================================================
# DENSE DIVERGENCE TRACE: Run 1 (DICE2023), every year 2020-2032
# ============================================================
print("\n\n===== DENSE TRACE: DICE2023 Run1 annual =====")
yr,tatm,atmc=simulate('DICE2023')
print(f"{'year':>6} {'atmC_sim':>11} {'atmC_ref':>11} {'d%':>7} {'TATM_sim':>10} {'TATM_ref':>10}")
for y in range(2020,2033):
    i=list(yr).index(y)
    ra=ref(1,'Climate.atmospheric carbon',y)
    rt=ref(1,'Climate.Atmospheric temperature anomaly (TATM)',y)
    dpct=100*(atmc[i]-ra)/ra
    print(f"{y:>6} {atmc[i]:>11.4f} {ra:>11.4f} {dpct:>7.3f} {tatm[i]:>10.5f} {rt:>10.5f}")

# Also check: what is atmospheric carbon actually in ref = MAT = Meq + sum(R)?
# ref 'atmospheric carbon' for Run1 (Lenton?=0) should equal Meq+R0+R1+R2+R3

# ============================================================
# DICE-C first-step Lenton flux budget
# ============================================================
print("\n\n===== DICE-C Lenton flux budget at 2020 =====")
p=P_
y=2020.0
Ca=887.0;Cv=550.0;Cs=1500.0;Col=730.0;Coh=140.0;Coi=10040.0;Cod=23860.0;Ldef=0.0
TATM=1.247154
# economy for EIND
L=interp_annual('L',y);aL=interp_annual('aL',y);sig=interp_annual('sigma',y)
K=295.0;miu=miu_of(y)
Ygross=aL*(L/1000)**0.7*K**0.3
EIND=(sig*Ygross)*(1-miu)
anomaly1=TATM;surfT=288+anomaly1;anomaly=anomaly1
Cvs=p['Cv0']-Ldef*p['kD']
pCO2=Ca/p['ka']*p['fgtm'];pCO2ppmv=pCO2*1e6
Photo=(p['kp']*Cvs*p['kMM'])*((pCO2ppmv-p['kc'])/(p['KM']+(pCO2ppmv-p['kc'])))*((15+anomaly)**2*(25-anomaly)/5625)
PlantR=p['kr']*Cv*p['kA']*np.exp(-p['Ea']/(p['R']*surfT))
Turn=p['kt']*Cv
SoilR=p['ksr']*Cs*p['kB']*2.7**(-1*(308.56/(surfT-227.13)))
bw=np.exp(-60.2409+93.4517*(100/p['Tl'])+23.3585*np.log(p['Tl']/100)+34.7*(0.023517-0.023656*(p['Tl']/100)+0.0047036*(p['Tl']/100)**2))
bc=np.exp(-60.2409+93.4517*(100/p['Th'])+23.3585*np.log(p['Th']/100)+34.7*(0.023517-0.023656*(p['Th']/100)+0.0047036*(p['Th']/100)**2))
pCO2w=p['pCO2w0']*(Col/p['Col0'])**p['beta']
pCO2c=p['pCO2c0']*(Coh/p['Coh0'])**p['beta']
Faw=p['ks']*p['aO']*p['Al']*bw*(pCO2-pCO2w)*12e-15*p['rho']
Fac=p['ks']*p['aO']*p['Ah']*bc*(pCO2-pCO2c)*12e-15*p['rho']
GHGemis=EIND*0.2727
print(f"pCO2ppmv={pCO2ppmv:.4f}")
print(f"Photo={Photo:.4f}  PlantR={PlantR:.4f}  SoilR={SoilR:.4f}  Turn={Turn:.4f}")
print(f"Faw={Faw:.4f}  Fac={Fac:.4f}  (AirSeaLL=-Faw={-Faw:.4f}, AirSeaHL=-Fac={-Fac:.4f})")
print(f"GHGemis={GHGemis:.4f}  EIND={EIND:.4f}")
dCa=PlantR+SoilR+(-Faw)+(-Fac)+GHGemis+0-Photo
print(f"dCa/dt = {dCa:.4f}  => Ca(2021)~{Ca+dCa:.3f} (ref 874.51, so need dCa~-12.49)")
print(f"  breakdown: +PlantR{PlantR:.2f} +SoilR{SoilR:.2f} -Faw{Faw:.2f} -Fac{Fac:.2f} +GHG{GHGemis:.2f} -Photo{Photo:.2f}")

# ---- DICE-CP permafrost internal diagnostic ----
print("\n===== DICE-CP permafrost trajectory =====")
deriv,lenton,perma=make_deriv('DICE-CP')
N=int(round((2420-2020)/0.01)); S=np.zeros(20)
S[0]=0.1477;S[1]=1.099454;S[2]=295.0;S[3]=633.5;S[4]=0.0;S[5]=0.518
S[6]=150.093;S[7]=102.698;S[8]=39.534;S[9]=6.1865
S[10]=887;S[11]=550;S[12]=1500;S[13]=730;S[14]=140;S[15]=10040;S[16]=23860;S[17]=0
S[18]=1035.0;S[19]=0.0
lag=int(round(200/0.01)); buf=[]; hist={}
snap={2020,2050,2100,2150,2200,2300,2420}
for k in range(N+1):
    y=2020+k*0.01; T=S[0]+S[1]; buf.append(T)
    hist['T200']=buf[k-lag] if y>=2220 else (1.25/(np.exp(2.0)-1))*(np.exp(0.01*(y-2020))-1)
    if abs(y-round(y))<1e-6 and int(round(y)) in snap:
        print(f"{int(round(y))}: Cperm={S[18]:.1f} thawedC={S[19]:.1f} T200={hist['T200']:.3f} TATM={T:.3f}")
    if k==N: break
    k1,_,_=deriv(y,S,hist);k2,_,_=deriv(y+0.005,S+0.005*k1,hist)
    k3,_,_=deriv(y+0.005,S+0.005*k2,hist);k4,_,_=deriv(y+0.01,S+0.01*k3,hist)
    S=S+0.01/6*(k1+2*k2+2*k3+k4)

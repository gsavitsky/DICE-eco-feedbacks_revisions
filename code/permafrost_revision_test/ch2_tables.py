"""
Chapter 2 LaTeX tables from the final corrected-model outputs (2026-11-17).
Reporting convention (user decision 2026-11-17): simulations bounded to 2020-2300; welfare and SCC use U_acc to 2300; the
DICE objective horizon); time-series quantities are reported to 2300 (peak values found within 2020-2300).
Inputs: DICE2023 Rproj/vis_corrected/{ch2_runs_summary.csv, ch2_permafrost_sweep.csv, Welfare sensitivity/by_variable}
Outputs: ../../tables/*.tex (bundle-relative)
"""
import json, math, os
import numpy as np
import pandas as pd

VC = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "R", "vis_corrected"))
OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables"))
os.makedirs(OUT, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
S = pd.read_csv(os.path.join(VC, "ch2_runs_summary.csv"))
W = pd.read_csv(os.path.join(VC, "ch2_permafrost_sweep.csv"))
BV = os.path.join(VC, "Welfare sensitivity", "by_variable")
TATM = pd.read_csv(os.path.join(BV, "Climate_Atmospheric_temperature_anomaly_TATM_.csv"))
DFRAC = pd.read_csv(os.path.join(BV, "Economy_damages_as_fraction_of_gross_output_Dfrac_.csv"))
MAT = pd.read_csv(os.path.join(BV, "Climate_atmospheric_carbon.csv"))
CONF = ['DICE 2023', 'DFAIR + permafrost', 'DICE-C', 'DICE-CP']
BASE = {'DICE 2023': 1, 'DFAIR + permafrost': 5, 'DICE-C': 9, 'DICE-CP': 13}
PPM = 1 / 2.1295            # ppm per GtC (Lenton k_a, f_gtm)


def stp_std(y):
    pi_, eta, sd = 0.5 * 0.0 + 0.0, 0.95, 0.0
    return None


# standard-discounting STP exactly as in the model (rartp = exp(0.001 + 0.5*pi) - 1, precautionary term)
import sys
sys.path.insert(0, HERE)
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    sys.path.insert(0, os.path.dirname(HERE))
    import dice_stella_validation as dsv
P = dsv.P_
RARTP = math.exp(0.001 + 0.5 * P['pi_']) - 1


def STP(y):
    Vc = min(P['sd_cons'] ** 2 * (y - 2020), P['sd_cons'] ** 2 * 235)
    return 1.0 / (1 + RARTP) ** (y - 2020) / (1 - 0.5 * Vc * P['eta'] ** 2) ** (y - 2020)


def fmt(x, d=2, sign=False):
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s.replace('-', '$-$')


def write(name, body):
    open(os.path.join(OUT, name), "w").write(body)
    print("wrote", name)


yrs = TATM.Years.values
H = yrs <= 2300

# ---------- 1. summary table (standard damage & discounting) ----------
rows = []
ref = None
for c in CONF:
    r = BASE[c]
    T = TATM[f'Run_{r}'].values[H]; i = int(np.argmax(T)); yp = int(yrs[i])
    D = DFRAC[f'Run_{r}'].values[H][i] * 100
    rec = dict(c=c, T=T[i], y=yp, D=D, stp=STP(yp), T2100=T[80], T2200=T[180], T2300=T[280],
               co2=MAT[f'Run_{r}'].values[80] * PPM)
    if ref is None:
        ref = rec
    rows.append(rec)
lines = []
for k, r in enumerate(rows):
    dT = '---' if k == 0 else fmt(r['T'] - ref['T'], 2, True)
    dD = '---' if k == 0 else fmt(r['D'] - ref['D'], 2, True)
    lines.append(f"{k + 1} -- {r['c']} & {r['T']:.2f} & {r['y']} & {dT} & {r['D']:.2f}\\% & {dD} & {r['stp']:.4f} & "
                 f"{r['T2100']:.2f} / {r['T2200']:.2f} / {r['T2300']:.2f} & {r['co2']:.0f} \\\\")
write("tab_summary.tex", r"""\begin{table}[H]
\centering
\caption{Peak temperature (T\textsubscript{peak}, 2020--2300) and its year, change relative to DICE~2023, damage fraction of gross output at the peak ($D_\text{frac}$), its change relative to DICE~2023, social time preference factor at the peak, temperature in 2100, 2200 and 2300, and atmospheric CO\textsubscript{2} in 2100, for the four configurations under the DICE~2023 cost--benefit optimal controls, standard damage and standard discounting.}
\small
\setlength{\tabcolsep}{3pt}
\begin{tabular}{lcccccccc}
\textbf{Configuration} & \textbf{T\textsubscript{peak}} & \textbf{Peak} & \textbf{$\Delta$T\textsubscript{peak}} & \textbf{$D_\text{frac}$} & \textbf{$\Delta D_\text{frac}$} & \textbf{STP} & \textbf{T\textsubscript{2100/2200/2300}} & \textbf{CO\textsubscript{2}, 2100} \\
 & \textbf{(\textdegree C)} & \textbf{year} & \textbf{(\textdegree C)} & \textbf{at peak} & \textbf{(pp)} & \textbf{at peak} & \textbf{(\textdegree C)} & \textbf{(ppm)} \\
\hline
""" + "\n".join(lines) + r"""
\end{tabular}
\label{tab:summary}
\end{table}
""")

# ---------- 2. welfare and SCC (main text) ----------
def cell(c, dam, disc):
    r = BASE[c] + (2 if dam == 'HS' else 0) + (1 if disc == 'zero' else 0)
    return S[S.Run == r].iloc[0]


def ref_cell(dam, disc):
    return cell('DICE 2023', dam, disc)


lines = []
for c in CONF:
    parts = []
    for disc in ('std', 'zero'):
        for dam in ('std', 'HS'):
            x = cell(c, dam, disc); base = cell('DICE 2023', 'std', disc)
            parts.append('---' if (c == 'DICE 2023' and dam == 'std') else fmt(100 * (x.Uacc_2300 / base.Uacc_2300 - 1), 2, True) + r'\%')
    scc_std = cell(c, 'std', 'std').SCC2020; scc_hs = cell(c, 'HS', 'std').SCC2020
    pz = []
    for dam in ('std', 'HS'):
        x, b = cell(c, dam, 'zero'), ref_cell(dam, 'zero')
        pz.append('---' if c == 'DICE 2023' else fmt(100 * (x.SCC2020 / b.SCC2020 - 1), 0, True) + r'\%')
    pst = []
    for dam, v in (('std', scc_std), ('HS', scc_hs)):
        b = ref_cell(dam, 'std').SCC2020
        pst.append(f"{v:.1f}" + ('' if c == 'DICE 2023' else f" ({fmt(100 * (v / b - 1), 1, True)}\\%)"))
    lines.append(f"{c} & {parts[0]} & {parts[1]} & {parts[2]} & {parts[3]} & {pst[0]} & {pst[1]} & {pz[0]} & {pz[1]} \\\\")
write("tab_welfare_scc.tex", r"""\begin{table}[H]
\centering
\caption{Welfare and social cost of carbon (SCC) by configuration. Welfare is the change in $U_\text{acc}$ (to 2300) relative to DICE~2023 with standard damage under the same discounting convention. The 2020 SCC is in 2019 US\$ per tCO\textsubscript{2} under standard discounting (percentage change relative to DICE~2023 with the same damage function in parentheses); under zero time preference only the percentage change is shown (dollar values in Table~\ref{tab:scc_full}). Std: standard damage; HS: Howard--Sterner damage.}
\label{tab:welfare_scc}
\small
\setlength{\tabcolsep}{3pt}
\begin{tabular}{lcccccccc}
\toprule
 & \multicolumn{4}{c}{Welfare change (\%)} & \multicolumn{2}{c}{SCC 2020, standard disc.} & \multicolumn{2}{c}{SCC 2020, zero TP} \\
\cmidrule(lr){2-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}
 & \multicolumn{2}{c}{Standard disc.} & \multicolumn{2}{c}{Zero TP} & & & & \\
Configuration & Std & HS & Std & HS & Std & HS & Std & HS \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ---------- 3. supplement: U_acc and SCC dollars, all 16 runs ----------
lines = []
for c in CONF:
    for dam in ('std', 'HS'):
        a, z = cell(c, dam, 'std'), cell(c, dam, 'zero')
        ba, bz = cell('DICE 2023', 'std', 'std'), cell('DICE 2023', 'std', 'zero')
        la = '---' if (c == 'DICE 2023' and dam == 'std') else fmt(100 * (a.Uacc_2300 / ba.Uacc_2300 - 1), 2, True) + r'\%'
        lz = '---' if (c == 'DICE 2023' and dam == 'std') else fmt(100 * (z.Uacc_2300 / bz.Uacc_2300 - 1), 2, True) + r'\%'
        lines.append(f"{c if dam == 'std' else ''} & {'Standard' if dam == 'std' else 'Howard--Sterner'} & "
                     f"{a.Uacc_2300:,.1f} & {la} & {a.SCC2020:,.1f} & {a.SCC2050:,.1f} & "
                     f"{z.Uacc_2300:,.0f} & {lz} & {z.SCC2020:,.0f} & {z.SCC2050:,.0f} \\\\".replace(',', '{,}'))
    lines.append(r"\addlinespace")
write("tab_scc_full.tex", r"""\begin{table}[H]
\centering
\caption{Accumulated welfare ($U_\text{acc}$ to 2300), welfare change relative to DICE~2023 with standard damage under the same discounting convention, and social cost of carbon in 2020 and 2050 (2019 US\$ per tCO\textsubscript{2}) for all sixteen runs. Under zero time preference the discount factor rises over time and the SCC is dominated by damages in later centuries.}
\label{tab:scc_full}
\footnotesize
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llrrrrrrrr}
\toprule
 & & \multicolumn{4}{c}{Standard discounting} & \multicolumn{4}{c}{Zero time preference} \\
\cmidrule(lr){3-6}\cmidrule(lr){7-10}
Configuration & Damage & $U_\text{acc}$ & Change & SCC\textsubscript{2020} & SCC\textsubscript{2050} & $U_\text{acc}$ & Change & SCC\textsubscript{2020} & SCC\textsubscript{2050} \\
\midrule
""" + "\n".join(lines[:-1]) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ---------- 4. thresholds from the sweep ----------
def need(g, metric, target):
    # above m = 8 the 2020 thawed pool is depleted by the spin-up and release no longer rises
    g = g[g.m <= 8].sort_values('PgC_per_C_2100'); x = g.PgC_per_C_2100.values
    y = g[metric].values * (-1 if metric == 'welfare_pct' else 1); t = target * (-1 if metric == 'welfare_pct' else 1)
    return np.interp(t, y, x) if y.max() >= t else np.nan

CAP = W[(W.carbon == 'Lenton') & (W.m <= 8)].PgC_per_C_2100.max()
targets = [('SCC_pct', 10, r'SCC\textsubscript{2020} $+10\%$'), ('SCC_pct', 25, r'SCC\textsubscript{2020} $+25\%$'),
           ('dT_peak', 0.1, r'Peak warming $+0.1$\textdegree C'), ('dT_peak', 0.25, r'Peak warming $+0.25$\textdegree C'),
           ('welfare_pct', -0.1, r'Welfare $-0.1\%$'), ('welfare_pct', -0.25, r'Welfare $-0.25\%$')]
cells = [('standard', 'standard'), ('Howard-Sterner', 'standard'), ('standard', 'zero time preference'), ('Howard-Sterner', 'zero time preference')]
lines = []
for metric, tg, lab in targets:
    vals = []
    for dam, disc in cells:
        v = need(W[(W.carbon == 'Lenton') & (W.damage == dam) & (W.discounting == disc)], metric, tg)
        vals.append(f'$>{CAP:.0f}$' if np.isnan(v) else f"{v:.0f}")
    lines.append(f"{lab} & " + " & ".join(vals) + r" \\")
write("tab_thresholds.tex", r"""\begin{table}[H]
\centering
\caption{Permafrost feedback strength (PgC released per \textdegree C of warming by 2100) required to reach a given effect in DICE-CP, by damage function and discounting convention, interpolated from the permafrost-strength sweep. For comparison, AR6 assesses 18 (3.1--41) PgC\,\textdegree C$^{-1}$ \citep{canadellGlobalCarbonOther2021} and \citet{schadelPermafrostWildfireCarbon2026} estimate ${\sim}63$ PgC\,\textdegree C$^{-1}$ including abrupt thaw and fire. With Hector's pool, the module's maximum is ${\sim}""" + f"{CAP:.0f}" + r"""$ PgC\,\textdegree C$^{-1}$; ``$>""" + f"{CAP:.0f}" + r"""$'' means the effect is not reached within the attainable range.}
\label{tab:thresholds}
\small
\begin{tabular}{lcccc}
\toprule
 & \multicolumn{2}{c}{Standard discounting} & \multicolumn{2}{c}{Zero time preference} \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Effect & Std damage & HS damage & Std damage & HS damage \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ---------- 5. initial conditions (supplement replacement) ----------
L20 = json.load(open(os.path.join(HERE, "lenton_km120_state2020.json")))
PF = json.load(open(os.path.join(HERE, "perma_state2020.json")))
g = lambda k: L20[k] * 12e-15
init = [
    (r"\multicolumn{4}{@{}l}{\textbf{Lenton carbon cycle pools}} \\", None),
    (r"$C_{a,0}$ & Atmospheric carbon & " + f"{L20['Ca']:,.1f}".replace(',', '{,}') + r" & PgC ($" + f"{L20['Ca'] * PPM:.1f}" + r"$ ppmv) \\", None),
    (r"$C_{v,0}$ & Vegetation carbon & " + f"{L20['Cv']:.1f}" + r" & PgC \\", None),
    (r"$C_{s,0}$ & Soil carbon & " + f"{L20['Cs']:,.1f}".replace(',', '{,}') + r" & PgC \\", None),
    (r"$D_{0}$ & Cumulative land-use emissions since 1800 & " + f"{L20['Dcum']:.1f}" + r" & PgC \\", None),
    (r"DIC$_{l,0}$ & Warm surface ocean carbon & " + f"{g('DIC_l'):.1f}" + r" & PgC \\", None),
    (r"DIC$_{h,0}$ & Cold surface ocean carbon & " + f"{g('DIC_h'):.1f}" + r" & PgC \\", None),
    (r"DIC$_{i,0}$ & Intermediate ocean carbon & " + f"{g('DIC_i'):,.0f}".replace(',', '{,}') + r" & PgC \\", None),
    (r"DIC$_{d,0}$ & Deep ocean carbon & " + f"{g('DIC_d'):,.0f}".replace(',', '{,}') + r" & PgC \\", None),
    (r"\multicolumn{4}{@{}l}{\textbf{Permafrost pools}} \\", None),
    (r"$C_{\text{perm},0}$ & Frozen permafrost carbon (spin-up from 865 PgC in 1800) & " + f"{PF['Cperm']:.1f}" + r" & PgC \\", None),
    (r"$C_{\text{thawed},0}$ & Thawed permafrost carbon & " + f"{PF['thawedC']:.1f}" + r" & PgC \\", None),
]
body = "\n".join(t for t, _ in init)
write("tab_init_conditions_carbon.tex", r"""% Replaces the "Simple carbon model (Lenton) pools" and "Permafrost pools" blocks of tab:init_conditions in Supplement.tex.
% DFAIR, climate and economy initial conditions are unchanged.
""" + body + "\n")

# ---------- 6. model evaluation ----------
write("tab_model_evaluation.tex", r"""\begin{table}[H]
\centering
\caption{Evaluation of the Lenton carbon cycle implementation against \citet{lentonLandOceanCarbon2000} (pre-industrial steady state, Fig.~1; extended business-as-usual experiment, Sect.~5.1, forced with his emissions and coupled to his energy balance) and against CMIP6 carbon-cycle feedback metrics for a 1\%\,yr$^{-1}$ CO\textsubscript{2} increase to 4$\times$CO\textsubscript{2} \citep{aroraCarbonConcentrationCarbon2020}.}
\label{tab:model_evaluation}
\small
\begin{tabular}{lcc}
\toprule
Quantity & Reference & This study \\
\midrule
\multicolumn{3}{@{}l}{\textit{Pre-industrial steady state (Lenton 2000, Fig.~1)}} \\
Atmospheric CO\textsubscript{2} (ppmv) & 280 & 280.0 \\
Ocean carbon, warm/cold/intermediate/deep (PgC) & 730 / 140 / 10{,}040 / 26{,}830 & 716 / 139 / 9{,}997 / 26{,}888 \\
Air--sea flux, warm surface in/out (PgC\,yr$^{-1}$) & 57.0 / 57.9 & 56.9 / 58.4 \\
Air--sea flux, cold surface in/out (PgC\,yr$^{-1}$) & 19.1 / 18.2 & 19.0 / 17.5 \\
\addlinespace
\multicolumn{3}{@{}l}{\textit{Extended business as usual (Lenton 2000, Sect.~5.1)}} \\
CO\textsubscript{2} in 1990 (ppmv) & 353 & 352.7 \\
Land / ocean sink in 1990 (PgC\,yr$^{-1}$) & 1.98 / 1.98 & 1.98 / 1.99 \\
Peak land sink (PgC\,yr$^{-1}$, year) & ${\sim}2.5$ (2030) & 2.43 (2028) \\
Land becomes a carbon source (year) & ${\sim}2120$ & 2113 \\
Carbon lost while land is a source (PgC) & 69 & 68 \\
Ocean sink in 2060; peak (PgC\,yr$^{-1}$, year) & 4.9; ${\sim}6.8$ (2130) & 4.89; 6.75 (2132) \\
Peak CO\textsubscript{2} (ppmv, year) & 985 (2170) & 984 (2174) \\
Peak warming (K, year) & 5.5 (2180) & 5.53 (2184) \\
\addlinespace
\multicolumn{3}{@{}l}{\textit{Carbon-cycle feedbacks, 1\%\,yr$^{-1}$ to 4$\times$CO\textsubscript{2} (CMIP6 mean $\pm$ s.d.)}} \\
Land $\beta$ (PgC\,ppm$^{-1}$) & $1.22 \pm 0.40$ & 0.69 \\
Ocean $\beta$ (PgC\,ppm$^{-1}$) & $0.79 \pm 0.07$ & 0.86 \\
Land $\gamma$ (PgC\,K$^{-1}$) & $-58 \pm 38$ & $-32$ \\
Ocean $\gamma$ (PgC\,K$^{-1}$) & $-17 \pm 5$ & 0 (fixed ocean surface temperature) \\
\bottomrule
\end{tabular}
\end{table}
""")

# ---------- 7. feedback attribution ----------
if os.path.exists(os.path.join(VC, "ch2_feedback_main_effects.csv")):   # factorial postponed
    FE = pd.read_csv(os.path.join(VC, "ch2_feedback_main_effects.csv")).set_index('feedback')
    NAMES = [('fert', r'Terrestrial CO\textsubscript{2} fertilisation', r'Photosynthetic uptake increases with atmospheric CO\textsubscript{2}', r'Eq.~\ref{eq:photosynthesis}', r'CO\textsubscript{2} in photosynthesis held at its 2020 value'),
             ('photoT', 'Photosynthesis temperature response', 'Temperature scaling of gross photosynthesis', r'Eq.~\ref{eq:photosynthesis}', r'$\Delta T$ in photosynthesis held at its 2020 value'),
             ('vegT', 'Plant respiration temperature response', 'Arrhenius temperature dependence of plant respiration', r'Sect.~\ref{eq:land}', r'$\Delta T$ in plant respiration held at its 2020 value'),
             ('soilT', 'Soil respiration temperature response', 'Lloyd--Taylor temperature dependence of soil respiration', r'Sect.~\ref{eq:land}', r'$\Delta T$ in soil respiration held at its 2020 value'),
             ('ocean', 'Ocean carbonate buffering', r'Surface-ocean $p$CO\textsubscript{2} rises nonlinearly with dissolved carbon, weakening uptake', r'Eqs.~\ref{eq:oceanflux}--\ref{eq:carbonate}', r'Surface-ocean $p$CO\textsubscript{2} held at its 2020 values'),
             ('perma', 'Permafrost carbon release', 'Gradual and abrupt thaw and respiration of thawed carbon', r'Eqs.~\ref{eq:frozenfrac}--\ref{eq:abruptthaw}', 'Permafrost module removed')]
    lines = []
    for k, n1, n2, ref, off in NAMES:
        e = FE.loc[k]
        lines.append(f"{n1} & {n2} & {ref} & {off} & {fmt(e.dT_peak_mean, 2, True)} ({fmt(e.dT_peak_min, 2, True)} to {fmt(e.dT_peak_max, 2, True)}) \\\\\n\\addlinespace")
    write("tab_feedbacks.tex", r"""\begin{table}[H]
    \centering
    \caption{Feedbacks in the Lenton carbon cycle and permafrost module switched in the full factorial analysis, how each is disengaged, and its effect on peak warming (2020--2300): mean difference between runs with the feedback on and off across the 32 combinations of the other five (range in parentheses). Standard damage and discounting, DICE~2023 optimal controls.}
    \label{tab:feedbacks}
    \footnotesize
    \begin{tabular}{>{\raggedright\arraybackslash}p{2.8cm} >{\raggedright\arraybackslash}p{3.6cm} p{1.6cm} >{\raggedright\arraybackslash}p{3.2cm} >{\raggedright\arraybackslash}p{2.6cm}}
    \toprule
    \textbf{Feedback} & \textbf{Description} & \textbf{Reference} & \textbf{Disengagement} & \textbf{$\Delta$T\textsubscript{peak} (\textdegree C)} \\
    \midrule
    """ + "\n".join(lines) + r"""
    \bottomrule
    \end{tabular}
    \end{table}
    """)

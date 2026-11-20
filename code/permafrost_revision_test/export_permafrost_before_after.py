"""Permafrost pools over time, original submission vs corrected (standard damages and discounting, GAMS controls).
Original: Python replication of the STELLA DICE-CP ('DICE-CP old': 1035 PgC frozen in 2020, step abrupt thaw), checked
against the STELLA temperature output. Corrected: ch2_models (Hector-faithful module, spin-up from 865 PgC in 1800),
on the Lenton carbon cycle (DICE-CP) and on DFAIR (DICE 2023 + permafrost).
Writes vis_corrected/permafrost_before_after.csv: year, version, config, frozen (PgC), thawed (PgC),
released (cumulative to the atmosphere since 2020, PgC), frozen_frac (frozen / pre-industrial frozen pool), TATM."""
import os
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis_corrected/permafrost_before_after.csv")
import ch2_models as M
prs = M.prs
yrs = np.arange(2020, 2201)
rows = []
d = prs.run(dict(prs.MODELS['DICE-CP old'], a2=0.003467, rartp=prs.DISC['std']))
old = pd.DataFrame(dict(year=np.arange(2020, 2421), version='Original submission', config='DICE-CP', frozen=d['Cperm'],
                        thawed=d['thawedC'], released=d['R'], frozen_frac=d['Cperm'] / 1035.0, TATM=d['TATM']))
rows.append(old[old.year <= 2200])
for label, carbon in (('DICE-CP', 'lenton'), ('DICE 2023 + permafrost', 'dfair')):
    r = M.run(dict(carbon=carbon, a2=0.003467, rartp=prs.DISC['std'], perma=dict(m=1.0, abrupt='smooth')), nyears=180)
    rows.append(pd.DataFrame(dict(year=yrs, version='Corrected', config=label, frozen=r['Cperm'], thawed=r['thawed'],
                                  released=r['release'], frozen_frac=r['Cperm'] / 865.0, TATM=r['TATM'])))
D = pd.concat(rows); D.to_csv(OUT, index=False)
st = pd.read_csv(os.path.expanduser("~/Desktop/DICE_eco_feedbacks/DICE2023 Rproj/vis/welfare sensitivity/by_variable/"
                                    "Climate_Atmospheric_temperature_anomaly_TATM_.csv")).set_index('Years')['Run_13']
o = old.set_index('year')
print('replication vs STELLA DICE-CP TATM, max |diff| 2020-2200:', float((o.TATM - st.reindex(o.index)).abs().loc[:2200].max()))
print(D[D.year.isin([2020, 2100, 2150, 2200])].round(3).to_string(index=False))

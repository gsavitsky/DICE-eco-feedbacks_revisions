# VM runs: Chapter 3 re-optimisation with Lenton variant (iv) (2026-11-20)

The model now uses Lenton's ocean-warming version (iv), with his own K_M = 145. This is the default in `dice_pyomo_corrected.py`; set `LENTON_VARIANT=iii` for the previous runs.
- The surface-ocean chemistry fit is now two-dimensional: dissolved carbon × ocean warming from −0.5 to +5 °C. Its maximum error is 1.6e-6.
- CO2 solubility in the air–sea flux follows the warming.

Forward check on the laptop (GAMS controls, 20-substep reference vs the reference forward model `ch2_models.py`):

| Configuration | Peak T, optimiser | Peak T, reference | U_acc, optimiser | U_acc, reference |
|---|---|---|---|---|
| DICE-C | 2.5880 | 2.5873 | 12651.217 | 12651.212 |
| DICE-CP | 2.7683 | 2.7678 | 12639.549 | 12639.533 |

- 2 substeps match the 20-substep reference to 2e-6 °C.
- 1 substep is still unstable.

## 1. Copy files to the VM (run on the laptop)

```bash
cd ~/Desktop/Dissertation/Optimizing_DICE-CP
scp dice_pyomo_corrected.py gsavitsk@beckage-rs.uvm.edu:~/Optimizing_DICE-CP/
scp permafrost_revision_test/lenton_km145_iv_state2020.json gsavitsk@beckage-rs.uvm.edu:~/Optimizing_DICE-CP/permafrost_revision_test/
```

`perma_state2020.json` is unchanged and already on the VM.

## 2. Launch on the VM: two configurations × two starts

- Start 1: the GAMS controls.
- Start 2: last week's corrected variant (iii) optimum (`*_corr_ss5.json`, already on the VM).

New output names end in `_iv`, so nothing from last week is overwritten.

```bash
cd ~/Optimizing_DICE-CP && conda activate dice
ls controls.json dicec_corr_ss5.json dicecp_corr_ss5.json permafrost_revision_test/lenton_km145_iv_state2020.json
export DICE_HSLLIB=~/coinhsl_install/lib64/libcoinhsl.so DICE_SUBSTEPS=2 LENTON_VARIANT=iv
for cfg in dicec dicecp; do
  DICE_CONFIG=$cfg DICE_WARM_START=controls.json DICE_OUTPUT_FILE=${cfg}_iv_ss2_gams.json DICE_TRAJ_FILE=traj_${cfg}_iv_ss2_gams.json \
    nohup python -u dice_pyomo_corrected.py > ${cfg}_iv_ss2_gams.log 2>&1 &
  DICE_CONFIG=$cfg DICE_WARM_START=${cfg}_corr_ss5.json DICE_OUTPUT_FILE=${cfg}_iv_ss2_alt.json DICE_TRAJ_FILE=traj_${cfg}_iv_ss2_alt.json \
    nohup python -u dice_pyomo_corrected.py > ${cfg}_iv_ss2_alt.log 2>&1 &
done
disown -a
```

The first line of each log should read `[dicec, Lenton variant iv, K_M=145] chemistry fit ...`. If it says variant iii, the old file was used.

## 3. Check progress

```bash
cd ~/Optimizing_DICE-CP
grep -h "Lenton variant\|EXIT\|Objective (padded)\|peak T" *_iv_ss2_*.log
```

Paste the output back. Once both starts agree, run the 5-substep confirmation from the 2-substep solution:

```bash
export DICE_HSLLIB=~/coinhsl_install/lib64/libcoinhsl.so DICE_SUBSTEPS=5 LENTON_VARIANT=iv
for cfg in dicec dicecp; do
  DICE_CONFIG=$cfg DICE_WARM_START=${cfg}_iv_ss2_gams.json DICE_OUTPUT_FILE=${cfg}_iv_ss5.json DICE_TRAJ_FILE=traj_${cfg}_iv_ss5.json \
    nohup python -u dice_pyomo_corrected.py > ${cfg}_iv_ss5.log 2>&1 &
done
disown -a
```

## 4. Copy results back (on the laptop)

```bash
cd ~/Desktop/Dissertation/Optimizing_DICE-CP
scp 'gsavitsk@beckage-rs.uvm.edu:~/Optimizing_DICE-CP/*_iv_ss*.json' .
```

## Notes

- **TATM ≥ 0.5 °C bound.** As before, it may become active near 2300.
- **Chemistry fit range.** The fit covers −10% to +6% of 2020 dissolved carbon and −0.5 to +5 °C of ocean warming. Along the GAMS path the model stays within −5% to +3% of dissolved carbon. If a DIC bound becomes active, report it.

# VM runs: Chapter 3 re-optimisation with the corrected physics (2026-11-16; updated 2026-11-17)

> **2026-11-17:** the permafrost pools in 2020 now come from a 1800–2020 spin-up (`perma_state2020.json`). If `dice_pyomo_corrected.py` was already copied to the VM, copy it again together with this file.

Model: `dice_pyomo_corrected.py`, the same physics as the Chapter 2 Python model (`permafrost_revision_test/ch2_models.py`).
- `DICE_CONFIG=dicec`: Lenton carbon cycle.
- `DICE_CONFIG=dicecp`: Lenton + Hector permafrost + smooth abrupt thaw.
- DICE 2023 is unchanged, so no re-solve is needed for it.

**Annual RK4 (1 substep) is numerically unstable** for the corrected surface ocean, so use **`DICE_SUBSTEPS=2`**. At 2 substeps the forward simulation matches a 20-substep reference to 6e-7 °C.

## 1. Copy files to the VM (run on the laptop)

```bash
cd ~/Desktop/Dissertation/Optimizing_DICE-CP
scp dice_pyomo_corrected.py gsavitsk@beckage-rs.uvm.edu:~/Optimizing_DICE-CP/
ssh gsavitsk@beckage-rs.uvm.edu 'mkdir -p ~/Optimizing_DICE-CP/permafrost_revision_test'
scp permafrost_revision_test/lenton_km120_state2020.json permafrost_revision_test/perma_state2020.json gsavitsk@beckage-rs.uvm.edu:~/Optimizing_DICE-CP/permafrost_revision_test/
```

## 2. Launch on the VM: two configurations × two starts, in parallel

```bash
cd ~/Optimizing_DICE-CP && conda activate dice
export DICE_HSLLIB=~/coinhsl_install/lib64/libcoinhsl.so DICE_SUBSTEPS=2
for cfg in dicec dicecp; do
  DICE_CONFIG=$cfg DICE_WARM_START=controls.json DICE_OUTPUT_FILE=${cfg}_corr_ss2_gams.json DICE_TRAJ_FILE=traj_${cfg}_corr_ss2_gams.json \
    nohup python -u dice_pyomo_corrected.py > ${cfg}_corr_ss2_gams.log 2>&1 &
done
DICE_CONFIG=dicec DICE_WARM_START=dicec_fine_ss5_controls.json DICE_OUTPUT_FILE=dicec_corr_ss2_alt.json DICE_TRAJ_FILE=traj_dicec_corr_ss2_alt.json \
  nohup python -u dice_pyomo_corrected.py > dicec_corr_ss2_alt.log 2>&1 &
DICE_CONFIG=dicecp DICE_WARM_START=dicecp_nb_fine_ss5_controls.json DICE_OUTPUT_FILE=dicecp_corr_ss2_alt.json DICE_TRAJ_FILE=traj_dicecp_corr_ss2_alt.json \
  nohup python -u dice_pyomo_corrected.py > dicecp_corr_ss2_alt.log 2>&1 &
disown -a
```

- Start 1 uses the GAMS controls. Start 2 uses the previous fine-grid optimum. Check the alternative start files exist on the VM first: `ls dicec_fine_ss5_controls.json dicecp_nb_fine_ss5_controls.json`.
- Expected size per run: about 1.7 M Jacobian non-zeros, similar to the earlier 2-substep runs (157 s on the VM).

## 3. Check progress

```bash
cd ~/Optimizing_DICE-CP
tail -3 dicec_corr_ss2_gams.log dicecp_corr_ss2_gams.log dicec_corr_ss2_alt.log dicecp_corr_ss2_alt.log
grep -h "EXIT\|Objective (padded)\|peak T" *_corr_ss2_*.log
```

Paste the output back. Once both starts agree, run a 5-substep confirmation (`DICE_SUBSTEPS=5`, warm-started from the 2-substep solution).

## Notes

- **TATM lower bound.** The GAMS lower bound TATM ≥ 0.5 °C is kept. Under the GAMS controls, the corrected DICE-C/DICE-CP paths dip below 0.5 °C around 2286–2400, so the bound will be active near 2300 unless the optimal policy differs. Override with `DICE_TATM_LO` if you want to test it.
- **Carbonate chemistry fit.** Surface pCO2 uses a degree-10 fit of the exact carbonate chemistry, with a maximum error of 4e-6. The fit is valid for −10% to +6% of 2020 DIC, and the DIC states are bounded to that range. If a bound becomes active, report it: the fit range would need widening.

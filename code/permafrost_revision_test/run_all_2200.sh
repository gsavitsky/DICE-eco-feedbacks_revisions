#!/bin/bash
# Re-cut Chapter 2 at the 2200 horizon (variant iv). Shapley (already at 2200) is not rerun.
cd "$(dirname "$0")"; LOG="$(pwd)/logs_2200"; mkdir -p "$LOG"
R=/Library/Frameworks/R.framework/Resources/bin/Rscript
RP="$HOME/Desktop/DICE_eco_feedbacks/DICE2023 Rproj"
step() { echo "=== START $1 $(date +%H:%M:%S)"; shift; "$@" > "$LOG/$(basename "$2").log" 2>&1 && echo "=== OK $(date +%H:%M:%S)" || echo "=== FAILED (exit $?) $(date +%H:%M:%S)"; }
step "1 grid+sweep+ceiling" python3 ch2_final_runs_parallel.py
step "2 tables" python3 ch2_tables.py
step "3 warming over time" python3 ch2_warming_time.py
step "4 permafrost before/after export" python3 export_permafrost_before_after.py
cd "$RP"
step "5 R corrected figures" $R Visualizations_corrected.R
step "6 R hooks figures" $R Visualizations_permafrost_hooks.R
step "7 R before/after figures" $R Visualizations_before_after.R
echo "=== ALL DONE $(date +%H:%M:%S)"

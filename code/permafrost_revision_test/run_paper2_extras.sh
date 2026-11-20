#!/bin/bash
cd "$(dirname "$0")"; LOG="$(pwd)/logs_extras"; mkdir -p "$LOG"
step() { echo "=== START $1 $(date +%H:%M:%S)"; shift; env "$@" > "$LOG/$(echo "$*" | tr -c 'A-Za-z0-9_.' '_' | cut -c1-80).log" 2>&1 && echo "=== OK $(date +%H:%M:%S)" || echo "=== FAILED (exit $?) $(date +%H:%M:%S)"; }
step "1 SSP check (variant iv)" python3 ssp_check_v2.py
step "2 gradual-only" python3 ch2_gradual_only.py
step "3 warming over time + time series" python3 ch2_warming_time.py
step "4 Shapley K_M=100" SHAPLEY_KM=100 SHAPLEY_TAG=_KM100 SHAPLEY_SKIP_DFAIR=1 python3 ch2_scc_shapley.py
step "5 Shapley K_M=275" SHAPLEY_KM=275 SHAPLEY_TAG=_KM275 SHAPLEY_SKIP_DFAIR=1 python3 ch2_scc_shapley.py
echo "=== ALL DONE $(date +%H:%M:%S)"

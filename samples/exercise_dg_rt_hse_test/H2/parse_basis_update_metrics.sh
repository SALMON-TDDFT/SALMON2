#!/bin/bash
set -euo pipefail
CSV="${1:-basis_update_metrics.csv}"
awk -F, 'NR==1{next} {n++; if($4==1)u++} END{printf "rows=%d updates=%d\n", n, u}' "$CSV"

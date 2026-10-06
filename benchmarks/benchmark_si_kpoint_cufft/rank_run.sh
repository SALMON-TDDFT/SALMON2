#!/bin/bash
set -eu
rank=${OMPI_COMM_WORLD_RANK:?Use an Open MPI/HPC-X launcher}
: "${SALMON_EXE:?SALMON_EXE is required}"
sampler=''
cleanup() {
  if [[ -n "$sampler" ]]; then
    kill "$sampler" 2>/dev/null || true
    wait "$sampler" 2>/dev/null || true
  fi
}
trap cleanup EXIT
# Node-wide GPU memory, sampled rather than an exact process/device allocation peak.
if command -v nvidia-smi >/dev/null; then
  nvidia-smi --query-gpu=timestamp,index,memory.used --format=csv -l 1 \
    > "gpu-memory.rank-${rank}.csv" 2> "gpu-memory.rank-${rank}.err" &
  sampler=$!
fi
if [[ $(uname -s) == Linux ]]; then
  /usr/bin/time -v -o "host-memory.rank-${rank}.txt" "$SALMON_EXE" < inputfile
else
  /usr/bin/time -l "$SALMON_EXE" < inputfile 2> "host-memory.rank-${rank}.txt"
fi

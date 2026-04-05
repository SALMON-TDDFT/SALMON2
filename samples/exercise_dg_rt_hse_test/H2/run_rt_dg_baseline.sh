#!/bin/bash
set -euo pipefail

SALMON_EXE="/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/bin/salmon"
BASE_DIR="/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test/H2"
WORKDIR="${1:-/tmp/dg_rt_smoke_4frag}"
INPUT="${2:-inputfile_h2_periodic_20_dg_rt_4frag_nt40_smoke}"
NPROC="${3:-1}"
DCDATADIR="${4:-}"

mkdir -p "$WORKDIR"
cd "$WORKDIR"
cp "$BASE_DIR/$INPUT" inputfile
ln -sf "$BASE_DIR/H_rps.dat" H_rps.dat

if [ ! -e data_dcdft ] && [ -n "$DCDATADIR" ]; then
  ln -s "$DCDATADIR" data_dcdft
fi

if [ ! -f data_dcdft/fragments/000001/wavefunctions.bin ]; then
  echo "Missing DC-LCFO data: data_dcdft/fragments/000001/wavefunctions.bin" >&2
  echo "Pass 4th arg with DC data directory, e.g. .../data_dcdft" >&2
  exit 2
fi

# Keep MPI settings in input aligned with runtime process count for smoke checks.
perl -0pi -e "s/nproc_ob\s*=\s*\d+/nproc_ob = ${NPROC}/g" inputfile
perl -0pi -e "s/nproc_rgrid_tot\(1:3\)\s*=\s*\d+,\s*\d+,\s*\d+/nproc_rgrid_tot = ${NPROC},1,1/g" inputfile

OMP_NUM_THREADS=1 mpirun -np "$NPROC" "$SALMON_EXE" < inputfile > run.log 2>&1

rg -n "total calculation time," run.log
if rg -n "SIGSEGV|Program received signal|error\(s\) in input" run.log; then
  echo "fatal pattern detected in run.log" >&2
  exit 1
fi

echo "Baseline run succeeded: $WORKDIR/run.log"

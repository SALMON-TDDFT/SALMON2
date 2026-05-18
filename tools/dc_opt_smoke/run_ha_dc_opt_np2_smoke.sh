#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TESTDIR="/tmp/dc_opt_ha_np2_smoke_$(date +%Y%m%d_%H%M%S)_$$"
INPUT_SRC="$REPO_ROOT/tools/dc_opt_smoke/input_ha_dc_opt_np2_smoke.inp"
INPUT_DST="$TESTDIR/input_ha_dc_opt_np2_smoke.inp"
EXE="$TESTDIR/ha"
PSEUDO_SRC="$REPO_ROOT/samples/exercise_dg_rt_hse_test/H2/H_rps.dat"
PSEUDO_DST="$TESTDIR/H_rps.dat"
LOG="$TESTDIR/run_ha_dc_opt_np2.log"
ERR="$TESTDIR/run_ha_dc_opt_np2.err"

mkdir -p "$TESTDIR"
cp "$REPO_ROOT/build/salmon" "$EXE"
chmod +x "$EXE"
cp "$INPUT_SRC" "$INPUT_DST"
cp "$PSEUDO_SRC" "$PSEUDO_DST"

cd "$REPO_ROOT"
mpirun --bind-to none --wdir "$TESTDIR" -np 2 "$EXE" < "$INPUT_DST" > "$LOG" 2> "$ERR"

printf 'TESTDIR=%s\n' "$TESTDIR"
printf 'LOG=%s\n' "$LOG"
printf 'ERR=%s\n' "$ERR"

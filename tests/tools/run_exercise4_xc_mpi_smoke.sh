#!/bin/bash
# run_exercise4_xc_mpi_smoke.sh -- local MPI smoke for exercise_04 xc paths

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

SAMPLE_DIR="$REPO_ROOT/samples/exercise_04_bulkSi_gs"
SERIAL_BINARY="$REPO_ROOT/build-gcc15-env/salmon"
MPI_BINARY="$REPO_ROOT/build-mpi-gcc15/salmon"
RUNTIME_ROOT="$REPO_ROOT/build-mpi-gcc15/runtime-checks/exercise4_xc_mpi_smoke_$(date +%Y%m%d-%H%M%S)"
NP=4
OMP_THREADS=1
GRID=16
GAP_TOL_EV="1e-5"
SKIP_SERIAL=0
SKIP_STRESS=0

usage() {
  cat <<EOF
Usage:
  $0 [options]

Options:
  --serial-binary PATH   Serial binary for PZ/r2scan reference
  --mpi-binary PATH      MPI binary to use for 4-rank runs
  --runtime-root DIR     Runtime directory to create/use
  --np INT               Number of MPI ranks for k-point parallel run (default: 4)
  --omp INT              OMP_NUM_THREADS value (default: 1)
  --grid INT             num_rgrid value per axis (default: 16)
  --gap-tol-ev FLOAT     Allowed |serial-mpi| gap difference for PZ/r2scan (default: 1e-5)
  --skip-serial          Skip serial PZ/r2scan reference runs
  --skip-stress          Skip MPI PZ stress and r2scan stress-guard runs
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --serial-binary) SERIAL_BINARY="$2"; shift 2 ;;
    --mpi-binary) MPI_BINARY="$2"; shift 2 ;;
    --runtime-root) RUNTIME_ROOT="$2"; shift 2 ;;
    --np) NP="$2"; shift 2 ;;
    --omp) OMP_THREADS="$2"; shift 2 ;;
    --grid) GRID="$2"; shift 2 ;;
    --gap-tol-ev) GAP_TOL_EV="$2"; shift 2 ;;
    --skip-serial) SKIP_SERIAL=1; shift ;;
    --skip-stress) SKIP_STRESS=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

require_file() {
  local path="$1"
  [[ -f "$path" ]] || { echo "ERROR: file not found: $path" >&2; exit 1; }
}

require_exec() {
  local path="$1"
  [[ -x "$path" ]] || { echo "ERROR: executable not found: $path" >&2; exit 1; }
}

extract_gap_ev() {
  awk '/Fundamental gap\[eV\]/{val=$3} END{if(val==""){exit 1} print val}' "$1"
}

extract_pressure_gpa() {
  awk '/pressure =/{val=$3} END{if(val==""){exit 1} print val}' "$1"
}

assert_close() {
  local a="$1"
  local b="$2"
  local tol="$3"
  awk -v a="$a" -v b="$b" -v tol="$tol" '
    BEGIN {
      d = a - b
      if (d < 0) d = -d
      exit !(d <= tol)
    }
  '
}

format_abs_diff() {
  local a="$1"
  local b="$2"
  awk -v a="$a" -v b="$b" 'BEGIN { d = a - b; if (d < 0) d = -d; printf "%.12g\n", d }'
}

prepare_dir() {
  local dir="$1"
  mkdir -p "$dir"
  cp "$SAMPLE_DIR/Si_rps.dat" "$dir/"
}

write_input() {
  local path="$1"
  local sysname="$2"
  local xc="$3"
  local mode="$4"        # serial | mpi
  local init_wf="$5"     # empty | random
  local stress="$6"      # yes | no
  local nscf="$7"

  {
    echo "! exercise_04_bulkSi_gs/Si_gs.inp based smoke input"
    echo
    echo "&calculation"
    echo "  theory = 'dft'"
    echo "/"
    echo
    echo "&control"
    echo "  sysname = '${sysname}'"
    echo "/"
    echo
    echo "&parallel"
    if [[ "$mode" == "mpi" ]]; then
      echo "  nproc_k = ${NP}"
      echo "  nproc_ob = 1"
      echo "  nproc_rgrid(1) = 1"
      echo "  nproc_rgrid(2) = 1"
      echo "  nproc_rgrid(3) = 1"
    fi
    echo "  yn_gramschmidt_blas = 'n'"
    echo "/"
    echo
    echo "&units"
    echo "  unit_system = 'A_eV_fs'"
    echo "/"
    echo
    echo "&system"
    echo "  yn_periodic = 'y'"
    echo "  al(1:3) = 5.43d0, 5.43d0, 5.43d0"
    echo "  nelem  = 1"
    echo "  natom  = 8"
    echo "  nelec  = 32"
    echo "  nstate = 32"
    echo "/"
    echo
    echo "&pseudo"
    echo "  file_pseudo(1) = './Si_rps.dat'"
    echo "  izatom(1) = 14"
    echo "  lloc_ps(1) = 2"
    echo "/"
    echo
    echo "&functional"
    echo "  xc = '${xc}'"
    echo "/"
    echo
    echo "&rgrid"
    echo "  num_rgrid(1:3) = ${GRID}, ${GRID}, ${GRID}"
    echo "/"
    echo
    echo "&kgrid"
    echo "  num_kgrid(1:3) = 4, 4, 4"
    echo "/"
    echo
    echo "&scf"
    if [[ -n "$init_wf" ]]; then
      echo "  method_init_wf = '${init_wf}'"
    fi
    echo "  nscf      = ${nscf}"
    echo "  threshold = 1.0d-9"
    echo "/"
    if [[ "$stress" == "yes" ]]; then
      echo
      echo "&analysis"
      echo "  yn_out_stress = 'y'"
      echo "  yn_out_stress_decomp = 'y'"
      echo "  stress_fd_detail = 'low'"
      echo "  out_stress_step = 1"
      echo "/"
    fi
    echo
    echo "&atomic_red_coor"
    echo "  'Si' .0   .0   .0   1"
    echo "  'Si' .25  .25  .25  1"
    echo "  'Si' .5   .0   .5   1"
    echo "  'Si' .0   .5   .5   1"
    echo "  'Si' .5   .5   .0   1"
    echo "  'Si' .75  .25  .75  1"
    echo "  'Si' .25  .75  .75  1"
    echo "  'Si' .75  .75  .25  1"
    echo "/"
  } > "$path"
}

run_serial_case() {
  local dir="$1"
  (
    cd "$dir"
    OMP_NUM_THREADS="$OMP_THREADS" "$SERIAL_BINARY" < inputfile > outputfile 2>&1
  )
}

run_mpi_case() {
  local dir="$1"
  (
    cd "$dir"
    OMP_NUM_THREADS="$OMP_THREADS" mpiexec -n "$NP" "$MPI_BINARY" < inputfile > outputfile 2>&1
  )
}

require_file "$SAMPLE_DIR/Si_rps.dat"
require_exec "$MPI_BINARY"
if [[ "$SKIP_SERIAL" -eq 0 ]]; then
  require_exec "$SERIAL_BINARY"
fi

mkdir -p "$RUNTIME_ROOT"

SERIAL_PZ_DIR="$RUNTIME_ROOT/serial_pz"
SERIAL_R2SCAN_DIR="$RUNTIME_ROOT/serial_r2scan"
MPI_PZ_DIR="$RUNTIME_ROOT/mpi_pz"
MPI_R2SCAN_DIR="$RUNTIME_ROOT/mpi_r2scan"
MPI_TBMBJ_DIR="$RUNTIME_ROOT/mpi_tbmbj_random"
MPI_PZ_STRESS_DIR="$RUNTIME_ROOT/mpi_pz_stress"
MPI_R2SCAN_GUARD_DIR="$RUNTIME_ROOT/mpi_r2scan_stress_guard"

if [[ "$SKIP_SERIAL" -eq 0 ]]; then
  prepare_dir "$SERIAL_PZ_DIR"
  write_input "$SERIAL_PZ_DIR/inputfile" "Si_ex4_pz_serial" "PZ" "serial" "" "no" "300"
  echo "Running serial PZ..."
  run_serial_case "$SERIAL_PZ_DIR"

  prepare_dir "$SERIAL_R2SCAN_DIR"
  write_input "$SERIAL_R2SCAN_DIR/inputfile" "Si_ex4_r2scan_serial" "r2scan" "serial" "" "no" "300"
  echo "Running serial r2scan..."
  run_serial_case "$SERIAL_R2SCAN_DIR"
fi

prepare_dir "$MPI_PZ_DIR"
write_input "$MPI_PZ_DIR/inputfile" "Si_ex4_pz_k${NP}mpi" "PZ" "mpi" "" "no" "300"
echo "Running MPI PZ..."
run_mpi_case "$MPI_PZ_DIR"

prepare_dir "$MPI_R2SCAN_DIR"
write_input "$MPI_R2SCAN_DIR/inputfile" "Si_ex4_r2scan_k${NP}mpi" "r2scan" "mpi" "" "no" "300"
echo "Running MPI r2scan..."
run_mpi_case "$MPI_R2SCAN_DIR"

prepare_dir "$MPI_TBMBJ_DIR"
write_input "$MPI_TBMBJ_DIR/inputfile" "Si_ex4_tbmbj_k${NP}mpi_random" "TBmBJ" "mpi" "random" "no" "300"
echo "Running MPI TBmBJ (random init)..."
run_mpi_case "$MPI_TBMBJ_DIR"

if [[ "$SKIP_STRESS" -eq 0 ]]; then
  prepare_dir "$MPI_PZ_STRESS_DIR"
  write_input "$MPI_PZ_STRESS_DIR/inputfile" "Si_ex4_pz_stress_k${NP}mpi" "PZ" "mpi" "" "yes" "300"
  echo "Running MPI PZ stress..."
  run_mpi_case "$MPI_PZ_STRESS_DIR"

  prepare_dir "$MPI_R2SCAN_GUARD_DIR"
  write_input "$MPI_R2SCAN_GUARD_DIR/inputfile" "Si_ex4_r2scan_stress_guard_k${NP}mpi" "r2scan" "mpi" "" "yes" "1"
  echo "Running MPI r2scan stress guard..."
  set +e
  run_mpi_case "$MPI_R2SCAN_GUARD_DIR"
  rc=$?
  set -e
  if ! grep -q "STOP stress tensor does not support xc='r2scan'" "$MPI_R2SCAN_GUARD_DIR/outputfile"; then
    echo "ERROR: r2scan stress guard message not found" >&2
    exit 1
  fi
  echo "r2scan stress guard exited with rc=$rc (expected non-zero from MPI launcher)"
fi

MPI_PZ_GAP=$(extract_gap_ev "$MPI_PZ_DIR/outputfile")
MPI_R2SCAN_GAP=$(extract_gap_ev "$MPI_R2SCAN_DIR/outputfile")
MPI_TBMBJ_GAP=$(extract_gap_ev "$MPI_TBMBJ_DIR/outputfile")

echo
echo "=== exercise_04 xc mpi smoke ==="
echo "runtime root: $RUNTIME_ROOT"
echo "MPI layout   : ${NP} MPI x ${OMP_THREADS} OMP (nproc_k=${NP})"
echo "grid         : ${GRID}x${GRID}x${GRID}"

if [[ "$SKIP_SERIAL" -eq 0 ]]; then
  SERIAL_PZ_GAP=$(extract_gap_ev "$SERIAL_PZ_DIR/outputfile")
  SERIAL_R2SCAN_GAP=$(extract_gap_ev "$SERIAL_R2SCAN_DIR/outputfile")
  PZ_DIFF=$(format_abs_diff "$SERIAL_PZ_GAP" "$MPI_PZ_GAP")
  R2SCAN_DIFF=$(format_abs_diff "$SERIAL_R2SCAN_GAP" "$MPI_R2SCAN_GAP")

  assert_close "$SERIAL_PZ_GAP" "$MPI_PZ_GAP" "$GAP_TOL_EV" || {
    echo "ERROR: PZ gap mismatch exceeds tolerance: serial=$SERIAL_PZ_GAP mpi=$MPI_PZ_GAP tol=$GAP_TOL_EV" >&2
    exit 1
  }
  assert_close "$SERIAL_R2SCAN_GAP" "$MPI_R2SCAN_GAP" "$GAP_TOL_EV" || {
    echo "ERROR: r2scan gap mismatch exceeds tolerance: serial=$SERIAL_R2SCAN_GAP mpi=$MPI_R2SCAN_GAP tol=$GAP_TOL_EV" >&2
    exit 1
  }

  echo "PZ gap [eV]       serial=$SERIAL_PZ_GAP mpi=$MPI_PZ_GAP |diff|=$PZ_DIFF"
  echo "r2scan gap [eV]   serial=$SERIAL_R2SCAN_GAP mpi=$MPI_R2SCAN_GAP |diff|=$R2SCAN_DIFF"
else
  echo "PZ gap [eV]       mpi=$MPI_PZ_GAP"
  echo "r2scan gap [eV]   mpi=$MPI_R2SCAN_GAP"
fi

echo "TBmBJ gap [eV]    mpi(random)=$MPI_TBMBJ_GAP"

if [[ "$SKIP_STRESS" -eq 0 ]]; then
  MPI_PZ_PRESSURE=$(extract_pressure_gpa "$MPI_PZ_STRESS_DIR"/Si_ex4_pz_stress_k"${NP}"mpi_info.data)
  echo "PZ pressure [GPa] mpi=$MPI_PZ_PRESSURE"
  echo "r2scan stress     guard message confirmed"
fi

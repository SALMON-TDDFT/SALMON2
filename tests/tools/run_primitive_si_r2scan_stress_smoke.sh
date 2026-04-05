#!/bin/bash
# run_primitive_si_r2scan_stress_smoke.sh -- local MPI smoke for primitive-cell Si nonorthogonal r2scan stress

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

PSEUDO_FILE="$REPO_ROOT/testsuites/pseudo/Si_rps.dat"
MPI_BINARY="$REPO_ROOT/build-mpi-gcc15/salmon"
RUNTIME_ROOT="$REPO_ROOT/build-mpi-gcc15/runtime-checks/primitive_si_r2scan_stress_smoke_$(date +%Y%m%d-%H%M%S)"
NP=4
OMP_THREADS=1
GRID=16
KGRID=4
NSCF=200
RUN_PZ=1
RUN_R2SCAN_GS=1
RUN_R2SCAN_STRESS=1
STRESS_OUTPUT_LEVEL="high"

usage() {
  cat <<EOF
Usage:
  $0 [options]

Options:
  --mpi-binary PATH      MPI binary to use
  --runtime-root DIR     Runtime directory to create/use
  --np INT               Number of MPI ranks for k-point parallel run (default: 4)
  --omp INT              OMP_NUM_THREADS value (default: 1)
  --grid INT             num_rgrid value per axis (default: 16)
  --kgrid INT            num_kgrid value per axis (default: 4)
  --nscf INT             Number of SCF iterations (default: 200)
  --stress-output-level LEVEL   Stress output level: low, middle, or high (default: high)
  --skip-pz              Skip primitive PZ stress control run
  --skip-r2scan-gs       Skip primitive r2scan GS run
  --skip-r2scan-stress   Skip primitive r2scan stress run
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mpi-binary) MPI_BINARY="$2"; shift 2 ;;
    --runtime-root) RUNTIME_ROOT="$2"; shift 2 ;;
    --np) NP="$2"; shift 2 ;;
    --omp) OMP_THREADS="$2"; shift 2 ;;
    --grid) GRID="$2"; shift 2 ;;
    --kgrid) KGRID="$2"; shift 2 ;;
    --nscf) NSCF="$2"; shift 2 ;;
    --stress-output-level) STRESS_OUTPUT_LEVEL="$2"; shift 2 ;;
    --skip-pz) RUN_PZ=0; shift ;;
    --skip-r2scan-gs) RUN_R2SCAN_GS=0; shift ;;
    --skip-r2scan-stress) RUN_R2SCAN_STRESS=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

case "$STRESS_OUTPUT_LEVEL" in
  low|middle|high) ;;
  *) echo "ERROR: unsupported stress output level: $STRESS_OUTPUT_LEVEL" >&2; usage; exit 1 ;;
esac

require_file() {
  local path="$1"
  [[ -f "$path" ]] || { echo "ERROR: file not found: $path" >&2; exit 1; }
}

require_exec() {
  local path="$1"
  [[ -x "$path" ]] || { echo "ERROR: executable not found: $path" >&2; exit 1; }
}

prepare_dir() {
  local dir="$1"
  mkdir -p "$dir"
  cp "$PSEUDO_FILE" "$dir/Si_rps.dat"
}

write_input() {
  local path="$1"
  local sysname="$2"
  local xc="$3"
  local stress="$4"

  {
    echo "! primitive Si nonorthogonal lattice smoke"
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
    echo "  nproc_k = ${NP}"
    echo "  nproc_ob = 1"
    echo "  nproc_rgrid(1) = 1"
    echo "  nproc_rgrid(2) = 1"
    echo "  nproc_rgrid(3) = 1"
    echo "  yn_gramschmidt_blas = 'n'"
    echo "/"
    echo
    echo "&units"
    echo "  unit_system = 'A_eV_fs'"
    echo "/"
    echo
    echo "&system"
    echo "  yn_periodic = 'y'"
    echo "  al_vec1 = 2.715d0, 2.715d0, 0.000d0"
    echo "  al_vec2 = 0.000d0, 2.715d0, 2.715d0"
    echo "  al_vec3 = 2.715d0, 0.000d0, 2.715d0"
    echo "  nelem  = 1"
    echo "  natom  = 2"
    echo "  nelec  = 8"
    echo "  nstate = 8"
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
    echo "  num_kgrid(1:3) = ${KGRID}, ${KGRID}, ${KGRID}"
    echo "/"
    echo
    echo "&scf"
    echo "  nscf = ${NSCF}"
    echo "  threshold = 1.0d-9"
    echo "/"
    if [[ "$stress" == "yes" ]]; then
      echo
      echo "&analysis"
      echo "  yn_out_stress = 'y'"
      echo "  stress_output_level = '${STRESS_OUTPUT_LEVEL}'"
      echo "  out_stress_step = 1"
      echo "/"
    fi
    echo
    echo "&atomic_red_coor"
    echo "  'Si' 0.00 0.00 0.00 1"
    echo "  'Si' 0.25 0.25 0.25 1"
    echo "/"
  } > "$path"
}

run_mpi_case() {
  local dir="$1"
  (
    cd "$dir"
    OMP_NUM_THREADS="$OMP_THREADS" mpiexec -n "$NP" "$MPI_BINARY" < inputfile > outputfile 2>&1
  )
}

require_file "$PSEUDO_FILE"
require_exec "$MPI_BINARY"
mkdir -p "$RUNTIME_ROOT"

if [[ "$RUN_PZ" -eq 1 ]]; then
  PZ_DIR="$RUNTIME_ROOT/primitive_pz_stress"
  prepare_dir "$PZ_DIR"
  write_input "$PZ_DIR/inputfile" "Si_primitive_pz_stress" "PZ" "yes"
  echo "Running primitive PZ stress..."
  run_mpi_case "$PZ_DIR"
fi

if [[ "$RUN_R2SCAN_GS" -eq 1 ]]; then
  R2SCAN_GS_DIR="$RUNTIME_ROOT/primitive_r2scan_gs"
  prepare_dir "$R2SCAN_GS_DIR"
  write_input "$R2SCAN_GS_DIR/inputfile" "Si_primitive_r2scan_gs" "r2scan" "no"
  echo "Running primitive r2scan GS..."
  run_mpi_case "$R2SCAN_GS_DIR"
fi

if [[ "$RUN_R2SCAN_STRESS" -eq 1 ]]; then
  R2SCAN_STRESS_DIR="$RUNTIME_ROOT/primitive_r2scan_stress"
  prepare_dir "$R2SCAN_STRESS_DIR"
  write_input "$R2SCAN_STRESS_DIR/inputfile" "Si_primitive_r2scan_stress" "r2scan" "yes"
  echo "Running primitive r2scan stress..."
  run_mpi_case "$R2SCAN_STRESS_DIR"
fi

echo
echo "Primitive Si nonorthogonal smoke completed."
echo "Results: $RUNTIME_ROOT"

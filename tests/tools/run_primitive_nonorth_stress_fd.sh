#!/bin/bash
# primitive nonorthogonal stress finite-difference helper

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

PSEUDO_FILE="$REPO_ROOT/testsuites/pseudo/Si_rps.dat"
MPI_BINARY="$REPO_ROOT/build-mpi-gcc15/salmon"
RUNTIME_ROOT="$REPO_ROOT/build-mpi-gcc15/runtime-checks/primitive_nonorth_stress_fd_$(date +%Y%m%d-%H%M%S)"
MODE="scale"
XC="PZ"
NP=4
OMP_THREADS=1
GRID=10
KGRID=4
NSCF=80
DELTA_LIST=""

usage() {
  cat <<EOF
Usage:
  $0 --mode MODE [options]

Modes:
  scale                Isotropic scaling around a reference conventional lattice constant
  shear                Cartesian x<-x+delta*y shear of the primitive lattice vectors

Options:
  --mode MODE          One of: scale, shear
  --xc XC              XC functional to use (default: PZ)
  --mpi-binary PATH    MPI binary to use
  --runtime-root DIR   Runtime directory to create/use
  --np INT             Number of MPI ranks for k-point parallel run (default: 4)
  --omp INT            OMP_NUM_THREADS value (default: 1)
  --grid INT           num_rgrid value per axis (default: 10)
  --kgrid INT          num_kgrid value per axis (default: 4)
  --nscf INT           Number of SCF iterations (default: 80)
  --delta-list LIST    Comma-separated list such as -0.005,0.000,0.005
  -h, --help           Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --xc) XC="$2"; shift 2 ;;
    --mpi-binary) MPI_BINARY="$2"; shift 2 ;;
    --runtime-root) RUNTIME_ROOT="$2"; shift 2 ;;
    --np) NP="$2"; shift 2 ;;
    --omp) OMP_THREADS="$2"; shift 2 ;;
    --grid) GRID="$2"; shift 2 ;;
    --kgrid) KGRID="$2"; shift 2 ;;
    --nscf) NSCF="$2"; shift 2 ;;
    --delta-list) DELTA_LIST="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

case "$MODE" in
  scale|shear) ;;
  *) echo "ERROR: unsupported mode: $MODE" >&2; usage; exit 1 ;;
esac

require_file() {
  local path="$1"
  [[ -f "$path" ]] || { echo "ERROR: file not found: $path" >&2; exit 1; }
}

require_exec() {
  local path="$1"
  [[ -x "$path" ]] || { echo "ERROR: executable not found: $path" >&2; exit 1; }
}

default_delta_list() {
  case "$MODE" in
    scale) echo "-0.005,0.000,0.005" ;;
    shear) echo "-0.005,0.000,0.005" ;;
  esac
}

format_tag() {
  python3 - "$1" <<'PY'
import sys
delta = float(sys.argv[1])
tag = f"{delta:+0.3f}".replace("+", "p").replace("-", "m").replace(".", "")
print(tag)
PY
}

write_input() {
  local path="$1"
  local sysname="$2"
  local a11="$3"
  local a21="$4"
  local a31="$5"
  local a12="$6"
  local a22="$7"
  local a32="$8"
  local a13="$9"
  local a23="${10}"
  local a33="${11}"

  {
    echo "&calculation"
    echo "  theory = 'dft'"
    echo "/"
    echo "&control"
    echo "  sysname = '${sysname}'"
    echo "/"
    echo "&parallel"
    echo "  nproc_k = ${NP}"
    echo "  nproc_ob = 1"
    echo "  nproc_rgrid(1) = 1"
    echo "  nproc_rgrid(2) = 1"
    echo "  nproc_rgrid(3) = 1"
    echo "  yn_gramschmidt_blas = 'n'"
    echo "/"
    echo "&units"
    echo "  unit_system = 'A_eV_fs'"
    echo "/"
    echo "&system"
    echo "  yn_periodic = 'y'"
    echo "  yn_symmetry = 'n'"
    echo "  al_vec1 = ${a11}d0, ${a21}d0, ${a31}d0"
    echo "  al_vec2 = ${a12}d0, ${a22}d0, ${a32}d0"
    echo "  al_vec3 = ${a13}d0, ${a23}d0, ${a33}d0"
    echo "  nelem  = 1"
    echo "  natom  = 2"
    echo "  nelec  = 8"
    echo "  nstate = 8"
    echo "/"
    echo "&pseudo"
    echo "  file_pseudo(1) = './Si_rps.dat'"
    echo "  izatom(1) = 14"
    echo "  lloc_ps(1) = 2"
    echo "/"
    echo "&functional"
    echo "  xc = '${XC}'"
    echo "/"
    echo "&rgrid"
    echo "  num_rgrid(1:3) = ${GRID}, ${GRID}, ${GRID}"
    echo "/"
    echo "&kgrid"
    echo "  num_kgrid(1:3) = ${KGRID}, ${KGRID}, ${KGRID}"
    echo "/"
    echo "&scf"
    echo "  nscf = ${NSCF}"
    echo "  threshold = 1.0d-9"
    echo "/"
    echo "&analysis"
    echo "  yn_out_stress = 'y'"
    echo "  yn_out_stress_decomp = 'y'"
    echo "  stress_fd_detail = 'C'"
    echo "  out_stress_step = 1"
    echo "/"
    echo "&atomic_red_coor"
    echo "  'Si' 0.00 0.00 0.00 1"
    echo "  'Si' 0.25 0.25 0.25 1"
    echo "/"
  } > "$path"
}

prepare_case() {
  local delta="$1"
  local case_dir="$2"
  local sysname="$3"

  local coords
  coords=$(python3 - "$MODE" "$delta" <<'PY'
import sys
mode = sys.argv[1]
delta = float(sys.argv[2])
a = 5.43 / 2.0
a1 = [a, a, 0.0]
a2 = [0.0, a, a]
a3 = [a, 0.0, a]
if mode == "scale":
    s = 1.0 + delta
    a1 = [s * x for x in a1]
    a2 = [s * x for x in a2]
    a3 = [s * x for x in a3]
elif mode == "shear":
    def shear(vec):
        x, y, z = vec
        return [x + delta * y, y, z]
    a1 = shear(a1)
    a2 = shear(a2)
    a3 = shear(a3)
else:
    raise SystemExit("unsupported mode")
print(*(f"{v:.10f}" for vec in (a1, a2, a3) for v in vec))
PY
)

  mkdir -p "$case_dir"
  cp "$PSEUDO_FILE" "$case_dir/Si_rps.dat"
  # shellcheck disable=SC2086
  write_input "$case_dir/inputfile" "$sysname" $coords
}

run_case() {
  local case_dir="$1"
  (
    cd "$case_dir"
    OMP_NUM_THREADS="$OMP_THREADS" mpiexec -n "$NP" "$MPI_BINARY" < inputfile > outputfile 2>&1
  )
}

require_file "$PSEUDO_FILE"
require_exec "$MPI_BINARY"
mkdir -p "$RUNTIME_ROOT"

if [[ -z "$DELTA_LIST" ]]; then
  DELTA_LIST=$(default_delta_list)
fi

OLD_IFS="$IFS"
IFS=','
set -- $DELTA_LIST
IFS="$OLD_IFS"

for delta in "$@"; do
  tag=$(format_tag "$delta")
  case_dir="$RUNTIME_ROOT/$tag"
  sysname="Si_primitive_${MODE}_${tag}"
  echo "Preparing ${MODE} delta=${delta} -> ${case_dir}"
  prepare_case "$delta" "$case_dir" "$sysname"
  echo "Running ${sysname} ..."
  run_case "$case_dir"
done

echo
echo "Primitive nonorthogonal stress finite-difference runs completed."
echo "Results: $RUNTIME_ROOT"

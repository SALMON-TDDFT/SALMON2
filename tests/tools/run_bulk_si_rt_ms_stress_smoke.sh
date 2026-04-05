#!/bin/bash
# run_bulk_si_rt_ms_stress_smoke.sh -- bulk Si RT/MS stress smoke helper
#
# This helper follows the sample-style restart flow:
#   1. run GS and produce data_for_restart
#   2. expose that directory as restart/
#   3. run RT or MS without forcing yn_restart='y'
#
# Example emitted input lines:
#   theory = 'tddft_response'
#   theory = 'multi_scale_maxwell_tddft'
#   xc = 'PZ'
#   xc = 'r2scan'

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

PSEUDO_FILE="$REPO_ROOT/testsuites/pseudo/Si_rps.dat"
MPI_BINARY="$REPO_ROOT/build-mpi-gcc15/salmon"
RUNTIME_ROOT="$REPO_ROOT/build-mpi-gcc15/runtime-checks/bulk_si_rt_ms_stress_smoke_$(date +%Y%m%d-%H%M%S)"
NP=4
OMP_THREADS=1
GRID=16
KGRID=4
GS_NSCF=160
RT_NT=4
MS_NT=2
STRESS_OUTPUT_LEVEL="high"
RUN_PRIMITIVE=1
RUN_MS=1

usage() {
  cat <<EOF
Usage:
  $0 [options]

Options:
  --mpi-binary PATH      MPI binary to use
  --runtime-root DIR     Runtime directory to create/use
  --np INT               Number of MPI ranks for k-point-parallel runs (default: 4)
  --omp INT              OMP_NUM_THREADS value (default: 1)
  --grid INT             num_rgrid value per axis (default: 16)
  --kgrid INT            num_kgrid value per axis (default: 4)
  --gs-nscf INT          Number of GS SCF iterations (default: 160)
  --rt-nt INT            Number of RT time steps (default: 4)
  --ms-nt INT            Number of MS time steps (default: 2)
  --stress-output-level LEVEL   Stress output level: low, middle, or high (default: high)
  --skip-primitive       Skip primitive RT sanity runs
  --skip-ms              Skip conventional MS runs
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
    --gs-nscf) GS_NSCF="$2"; shift 2 ;;
    --rt-nt) RT_NT="$2"; shift 2 ;;
    --ms-nt) MS_NT="$2"; shift 2 ;;
    --stress-output-level) STRESS_OUTPUT_LEVEL="$2"; shift 2 ;;
    --skip-primitive) RUN_PRIMITIVE=0; shift ;;
    --skip-ms) RUN_MS=0; shift ;;
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

link_restart_dir() {
  local dir="$1"
  local src="$2"
  (
    cd "$dir"
    ln -s "$src" restart
  )
}

run_mpi_case() {
  local dir="$1"
  (
    cd "$dir"
    OMP_NUM_THREADS="$OMP_THREADS" mpiexec -n "$NP" "$MPI_BINARY" < inputfile > outputfile 2>&1
  )
}

write_common_prefix() {
  local sysname="$1"
  cat <<EOF
&control
  sysname = '${sysname}'
/

&parallel
  nproc_k = ${NP}
  nproc_ob = 1
  nproc_rgrid(1) = 1
  nproc_rgrid(2) = 1
  nproc_rgrid(3) = 1
  yn_gramschmidt_blas = 'n'
/

&units
  unit_system = 'A_eV_fs'
/

&pseudo
  file_pseudo(1) = './Si_rps.dat'
  izatom(1) = 14
  lloc_ps(1) = 2
/

&rgrid
  num_rgrid(1:3) = ${GRID}, ${GRID}, ${GRID}
/

&kgrid
  num_kgrid(1:3) = ${KGRID}, ${KGRID}, ${KGRID}
/
EOF
}

write_conventional_geometry() {
  cat <<'EOF'
&system
  yn_periodic = 'y'
  al(1:3) = 5.43d0, 5.43d0, 5.43d0
  nelem  = 1
  natom  = 8
  nelec  = 32
  nstate = 32
/

&atomic_red_coor
  'Si' .0   .0   .0   1
  'Si' .25  .25  .25  1
  'Si' .5   .0   .5   1
  'Si' .0   .5   .5   1
  'Si' .5   .5   .0   1
  'Si' .75  .25  .75  1
  'Si' .25  .75  .75  1
  'Si' .75  .75  .25  1
/
EOF
}

write_primitive_geometry() {
  cat <<'EOF'
&system
  yn_periodic = 'y'
  al_vec1 = 2.715d0, 2.715d0, 0.000d0
  al_vec2 = 0.000d0, 2.715d0, 2.715d0
  al_vec3 = 2.715d0, 0.000d0, 2.715d0
  nelem  = 1
  natom  = 2
  nelec  = 8
  nstate = 8
/

&atomic_red_coor
  'Si' 0.00 0.00 0.00 1
  'Si' 0.25 0.25 0.25 1
/
EOF
}

write_gs_input() {
  local path="$1"
  local sysname="$2"
  local xc_name="$3"
  local cell_kind="$4"  # conventional | primitive

  {
    echo "! bulk Si conventional/primitive GS smoke"
    echo
    echo "&calculation"
    echo "  theory = 'dft'"
    echo "/"
    echo
    write_common_prefix "$sysname"
    echo
    echo "&functional"
    echo "  xc = '${xc_name}'"
    echo "/"
    echo
    echo "&scf"
    echo "  nscf = ${GS_NSCF}"
    echo "  threshold = 1.0d-9"
    echo "/"
    echo
    if [[ "$cell_kind" == "conventional" ]]; then
      write_conventional_geometry
    else
      write_primitive_geometry
    fi
  } > "$path"
}

write_rt_input() {
  local path="$1"
  local sysname="$2"
  local xc_name="$3"
  local cell_kind="$4"  # conventional | primitive
  local stress_step="$5"

  {
    echo "! bulk Si conventional/primitive RT stress smoke"
    echo
    echo "&calculation"
    echo "  theory = 'tddft_response'"
    echo "/"
    echo
    write_common_prefix "$sysname"
    echo
    echo "&functional"
    echo "  xc = '${xc_name}'"
    echo "/"
    echo
    echo "&tgrid"
    echo "  dt = 0.002d0"
    echo "  nt = ${RT_NT}"
    echo "/"
    echo
    echo "&emfield"
    echo "  ae_shape1 = 'impulse'"
    echo "  epdir_re1(1:3) = 0.0d0, 0.0d0, 1.0d0"
    echo "/"
    echo
    echo "&analysis"
    echo "  yn_out_stress = 'y'"
    echo "  stress_output_level = '${STRESS_OUTPUT_LEVEL}'"
    echo "  stress_l_decomp = 'species'"
    echo "  out_stress_step = ${stress_step}"
    echo "  de = 0.05d0"
    echo "  nenergy = 40"
    echo "/"
    echo
    if [[ "$cell_kind" == "conventional" ]]; then
      write_conventional_geometry
    else
      write_primitive_geometry
    fi
  } > "$path"
}

write_ms_input() {
  local path="$1"
  local sysname="$2"
  local xc_name="$3"

  {
    echo "! bulk Si conventional MS stress smoke"
    echo
    echo "&calculation"
    echo "  theory = 'multi_scale_maxwell_tddft'"
    echo "/"
    echo
    write_common_prefix "$sysname"
    echo
    echo "&functional"
    echo "  xc = '${xc_name}'"
    echo "/"
    echo
    echo "&tgrid"
    echo "  dt = 0.002d0"
    echo "  nt = ${MS_NT}"
    echo "/"
    echo
    echo "&emfield"
    echo "  ae_shape1 = 'Acos2'"
    echo "  I_wcm2_1 = 1.0d10"
    echo "  tw1 = 0.4d0"
    echo "  omega1 = 1.55d0"
    echo "  epdir_re1(1:3) = 0.0d0, 0.0d0, 1.0d0"
    echo "/"
    echo
    echo "&analysis"
    echo "  yn_out_stress = 'y'"
    echo "  stress_output_level = '${STRESS_OUTPUT_LEVEL}'"
    echo "  stress_l_decomp = 'species'"
    echo "  out_stress_step = 1"
    echo "  out_ms_step = 1"
    echo "/"
    echo
    echo "&multiscale"
    echo "  nx_m = 1"
    echo "  ny_m = 1"
    echo "  nz_m = 1"
    echo "  hx_m = 20.0d0"
    echo "  hy_m = 20.0d0"
    echo "  hz_m = 20.0d0"
    echo "  nxvacl_m = 2"
    echo "  nxvacr_m = 2"
    echo "/"
    echo
    echo "&maxwell"
    echo "  boundary_em(1,1) = 'abc'"
    echo "  boundary_em(1,2) = 'abc'"
    echo "/"
    echo
    write_conventional_geometry
  } > "$path"
}

assert_exists() {
  local path="$1"
  [[ -e "$path" ]] || { echo "ERROR: expected artifact not found: $path" >&2; exit 1; }
}

run_gs_case() {
  local dir="$1"
  local sysname="$2"
  local xc_name="$3"
  local cell_kind="$4"

  prepare_dir "$dir"
  write_gs_input "$dir/inputfile" "$sysname" "$xc_name" "$cell_kind"
  echo "Running ${cell_kind} ${xc_name} GS..."
  run_mpi_case "$dir"
  assert_exists "$dir/data_for_restart"
}

run_rt_case() {
  local dir="$1"
  local sysname="$2"
  local xc_name="$3"
  local cell_kind="$4"
  local stress_step="$5"
  local restart_src="$6"

  prepare_dir "$dir"
  write_rt_input "$dir/inputfile" "$sysname" "$xc_name" "$cell_kind" "$stress_step"
  link_restart_dir "$dir" "$restart_src"
  echo "Running ${cell_kind} ${xc_name} RT stress..."
  run_mpi_case "$dir"
  assert_exists "$dir/${sysname}_stress.data"
}

run_ms_case() {
  local dir="$1"
  local sysname="$2"
  local xc_name="$3"
  local restart_src="$4"

  prepare_dir "$dir"
  write_ms_input "$dir/inputfile" "$sysname" "$xc_name"
  link_restart_dir "$dir" "$restart_src"
  echo "Running conventional ${xc_name} MS stress..."
  run_mpi_case "$dir"
  assert_exists "$dir/${sysname}_m/m000001/${sysname}_stress.data"
}

require_file "$PSEUDO_FILE"
require_exec "$MPI_BINARY"
mkdir -p "$RUNTIME_ROOT"

CONV_PZ_GS_DIR="$RUNTIME_ROOT/conventional_pz_gs"
CONV_R2SCAN_GS_DIR="$RUNTIME_ROOT/conventional_r2scan_gs"
PRIM_PZ_GS_DIR="$RUNTIME_ROOT/primitive_pz_gs"
PRIM_R2SCAN_GS_DIR="$RUNTIME_ROOT/primitive_r2scan_gs"

run_gs_case "$CONV_PZ_GS_DIR" "Si_conv_pz_gs" "PZ" "conventional"
run_gs_case "$CONV_R2SCAN_GS_DIR" "Si_conv_r2scan_gs" "r2scan" "conventional"

run_rt_case "$RUNTIME_ROOT/conventional_pz_rt" "Si_conv_pz_rt" "PZ" "conventional" 2 "$CONV_PZ_GS_DIR/data_for_restart"
run_rt_case "$RUNTIME_ROOT/conventional_r2scan_rt" "Si_conv_r2scan_rt" "r2scan" "conventional" 2 "$CONV_R2SCAN_GS_DIR/data_for_restart"

if [[ "$RUN_PRIMITIVE" -eq 1 ]]; then
  run_gs_case "$PRIM_PZ_GS_DIR" "Si_primitive_pz_gs" "PZ" "primitive"
  run_gs_case "$PRIM_R2SCAN_GS_DIR" "Si_primitive_r2scan_gs" "r2scan" "primitive"
  run_rt_case "$RUNTIME_ROOT/primitive_pz_rt" "Si_primitive_pz_rt" "PZ" "primitive" 2 "$PRIM_PZ_GS_DIR/data_for_restart"
  run_rt_case "$RUNTIME_ROOT/primitive_r2scan_rt" "Si_primitive_r2scan_rt" "r2scan" "primitive" 2 "$PRIM_R2SCAN_GS_DIR/data_for_restart"
fi

if [[ "$RUN_MS" -eq 1 ]]; then
  run_ms_case "$RUNTIME_ROOT/conventional_pz_ms" "Si_conv_pz_ms" "PZ" "$CONV_PZ_GS_DIR/data_for_restart"
  run_ms_case "$RUNTIME_ROOT/conventional_r2scan_ms" "Si_conv_r2scan_ms" "r2scan" "$CONV_R2SCAN_GS_DIR/data_for_restart"
fi

echo
echo "RT and MS stress smoke completed."
echo "Results: $RUNTIME_ROOT"

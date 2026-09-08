#!/bin/sh
set -eu

case_dir="${1:-}"
salmon_bin="${2:-}"
np="${3:-8}"
input_auto_n="${4:-inputfile_rt_auto_n_nt1}"
input_auto_y="${5:-inputfile_rt_auto_y_nt1}"

if [ -z "$case_dir" ] || [ -z "$salmon_bin" ]; then
  echo "usage: $0 CASE_DIR SALMON_BIN [NP] [INPUT_AUTO_N] [INPUT_AUTO_Y]" >&2
  exit 2
fi

if [ ! -d "$case_dir" ]; then
  echo "missing case directory: $case_dir" >&2
  exit 1
fi

if [ ! -x "$salmon_bin" ]; then
  echo "missing executable salmon binary: $salmon_bin" >&2
  exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
csv_tool="$script_dir/dg_bpw_auto_report_to_csv.py"

if [ ! -x "$csv_tool" ] && [ ! -f "$csv_tool" ]; then
  echo "missing CSV helper: $csv_tool" >&2
  exit 1
fi

cd "$case_dir"

if [ ! -f "$input_auto_n" ]; then
  echo "missing auto=n input: $case_dir/$input_auto_n" >&2
  exit 1
fi

if [ ! -f "$input_auto_y" ]; then
  echo "missing auto=y input: $case_dir/$input_auto_y" >&2
  exit 1
fi

sys_auto_n=$(sed -n "s/.*sysname *= *'\([^']*\)'.*/\1/p" "$input_auto_n" | head -n 1)
sys_auto_y=$(sed -n "s/.*sysname *= *'\([^']*\)'.*/\1/p" "$input_auto_y" | head -n 1)

if [ -z "$sys_auto_n" ] || [ -z "$sys_auto_y" ]; then
  echo "failed to read sysname from inputs" >&2
  exit 1
fi

common_env="-x OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} -x SALMON_DG_EXPDIAG_GLOBAL_FLUX=1 -x SALMON_DG_BPW_SELECT_MODE=shell_ecut -x SALMON_DG_BPW_ECUT=${SALMON_DG_BPW_ECUT:-2.0} -x SALMON_DG_MIXED_SPERP_TOL=${SALMON_DG_MIXED_SPERP_TOL:-1e-2} -x SALMON_DG_MIXED_Z=1"

# shellcheck disable=SC2086
mpirun -np "$np" $common_env "$salmon_bin" < "$input_auto_n" > check_dg_bpw_auto_n.log 2>&1
# shellcheck disable=SC2086
mpirun -np "$np" $common_env "$salmon_bin" < "$input_auto_y" > check_dg_bpw_auto_y.log 2>&1

pol_auto_n="${sys_auto_n}_dg_polarization.data"
pol_auto_y="${sys_auto_y}_dg_polarization.data"

if [ ! -f "$pol_auto_n" ] || [ ! -f "$pol_auto_y" ]; then
  echo "missing polarization output: $pol_auto_n or $pol_auto_y" >&2
  exit 1
fi

if ! cmp -s "$pol_auto_n" "$pol_auto_y"; then
  echo "auto=n/y polarization output differs" >&2
  exit 1
fi

if [ ! -f dg_bpw_auto_report.dat ]; then
  echo "auto=y did not generate dg_bpw_auto_report.dat" >&2
  exit 1
fi

if ! grep -q '^SUMMARY ' dg_bpw_auto_report.dat; then
  echo "dg_bpw_auto_report.dat has no SUMMARY line" >&2
  exit 1
fi

python3 "$csv_tool" dg_bpw_auto_report.dat --output dg_bpw_auto_check.csv

if ! head -n 1 dg_bpw_auto_check.csv | grep -q 'recommendation'; then
  echo "CSV output has no recommendation column" >&2
  exit 1
fi

echo "ok: dg_bpw_auto report-only mode leaves polarization unchanged and writes report/CSV"

#!/usr/bin/env bash
set -uo pipefail
cd /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac || exit 99
export OMP_NUM_THREADS=2
salmon_bin=../../../build-mpi-eigenexa-wannier-lib/salmon
input=inputfile_stage2c_laser_full_tddft_lcfo_dt002_nt22500
runlog=run_stage2c_laser_full_tddft_lcfo_dt002_nt22500.log
status=stage2c_laser_full_tddft_lcfo_dt002_nt22500.status
{
  echo "START $(date '+%Y-%m-%d %H:%M:%S')"
  echo "PWD=$PWD"
  echo "CMD=OMP_NUM_THREADS=$OMP_NUM_THREADS mpirun -np 8 $salmon_bin < $input"
} > "$status"
mpirun -np 8 "$salmon_bin" < "$input" > "$runlog" 2>&1
exit_code=$?
echo "SALMON_EXIT_CODE=$exit_code" >> "$status"
echo "SALMON_END $(date '+%Y-%m-%d %H:%M:%S')" >> "$status"
if [[ $exit_code -ne 0 ]]; then
  exit "$exit_code"
fi
if ! grep -q 'end SALMON' "$runlog"; then
  echo "COMPARE_SKIPPED=no_end_salmon" >> "$status"
  exit 2
fi
cd /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG || exit 99
base=samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac
python3 tools/compare_global_wfpw_hhg_sweep.py \
  --reference "$base/stage2c_laser_full_tddft_lcfo_dt002_nt22500_rt.data" \
  --reference-source rt_current_integral \
  --manifest "$base/stage2c_laser_global_wfpw_manifest.tsv" \
  --axis z --tmax-au 450.0 --detrend linear --damping 0.0 --pad-factor 4 \
  --emin 0.0 --emax 120.0 --hhg-emin 1.55 --hhg-emax 120.0 \
  --fundamental-ev 1.55 --harmonic-orders 3,5,7,9,11,13,15,17,19,21 \
  > "$base/stage2c_laser_full_vs_wfpw_summary.tsv"
compare_exit=$?
echo "COMPARE_PLUS_EXIT_CODE=$compare_exit" >> "$base/stage2c_laser_full_tddft_lcfo_dt002_nt22500.status"
python3 tools/compare_global_wfpw_hhg_sweep.py \
  --reference "$base/stage2c_laser_full_tddft_lcfo_dt002_nt22500_rt.data" \
  --reference-source rt_current_integral_minus \
  --manifest "$base/stage2c_laser_global_wfpw_manifest.tsv" \
  --axis z --tmax-au 450.0 --detrend linear --damping 0.0 --pad-factor 4 \
  --emin 0.0 --emax 120.0 --hhg-emin 1.55 --hhg-emax 120.0 \
  --fundamental-ev 1.55 --harmonic-orders 3,5,7,9,11,13,15,17,19,21 \
  > "$base/stage2c_laser_full_vs_wfpw_summary_minusJ.tsv"
compare_minus_exit=$?
echo "COMPARE_MINUS_EXIT_CODE=$compare_minus_exit" >> "$base/stage2c_laser_full_tddft_lcfo_dt002_nt22500.status"
python3 tools/plot_global_wfpw_hhg_sweep.py \
  --summary "$base/stage2c_laser_full_vs_wfpw_summary.tsv" \
  --plot-dir "$base/stage2c_laser_full_vs_wfpw_plots" \
  --x-axis PW_count \
  --metric cutoff_shift_eV --metric hhg_window_rms_error \
  --metric H03_rel_intensity_error --metric H05_rel_intensity_error --metric H07_rel_intensity_error \
  > "$base/stage2c_laser_full_vs_wfpw_plot.log" 2>&1
plot_exit=$?
echo "PLOT_EXIT_CODE=$plot_exit" >> "$base/stage2c_laser_full_tddft_lcfo_dt002_nt22500.status"
echo "COMPARE_DONE $(date '+%Y-%m-%d %H:%M:%S')" >> "$base/stage2c_laser_full_tddft_lcfo_dt002_nt22500.status"
exit $(( compare_exit || compare_minus_exit || plot_exit ))

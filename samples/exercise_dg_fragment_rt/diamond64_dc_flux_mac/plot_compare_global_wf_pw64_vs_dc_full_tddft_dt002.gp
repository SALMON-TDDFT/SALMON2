set terminal pngcairo size 1500,1200 enhanced font "Arial,18"
set output "compare_global_wf_pw64_vs_dc_full_tddft_dt002_hhg.png"
set multiplot layout 2,1 title "C64 laser response: global WF+PW64 vs DC + Full TDDFT"

set grid
set key outside top center horizontal
set xlabel "Time (a.u.)"
set ylabel "Polarization change (a.u.)"
plot "compare_global_wf_pw64_vs_dc_full_tddft_dt002_waveform.tsv" using 1:2 with lines lw 2 title "DC + Full TDDFT: integral Jz", \
     "" using 1:5 with lines lw 2 title "Global WF+PW64: Pz", \
     "" using 1:3 with lines lw 1 title "Global WF+PW64: Px", \
     "" using 1:4 with lines lw 1 title "Global WF+PW64: Py"

set logscale y
set xrange [0.5:15.5]
set yrange [1e-12:5]
set xtics 1
set xlabel "Harmonic order"
set ylabel "HHG intensity / H1 integral"
plot "compare_global_wf_pw64_vs_dc_full_tddft_dt002_spectrum.tsv" using 1:2 with lines lw 2 title "DC + Full TDDFT (dt=0.02)", \
     "" using 1:3 with lines lw 2 title "Global WF+PW64 (Exp, dt=0.1)"

unset multiplot

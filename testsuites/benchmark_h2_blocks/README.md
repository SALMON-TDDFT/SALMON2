# H2 block scaling

Each fragment core contains 2x2x2 H2 on the previous 8-bohr molecular lattice:
16-cubed-bohr cores, .5-bohr mesh, 1.4-bohr bonds. Buffered periodic DC-SCF
uses canonical exchange (no MLWF), 300 K common chemical potential, and 4-bohr
buffers on split axes only. Four extra fragment states are included beyond the
padded-cell occupied count. Random initial orbitals avoid the Gaussian pilot's
initial core-weighted capacity failure. PBE pre-SCF precedes PBEh40, with
Coulomb cutoff 4 bohr and no rVV10. This cutoff is held fixed for performance
comparison, not claimed physically converged.

ScaLAPACK LCFO output is reconstructed onto the real-space mesh. RT propagates
the occupied states, with no assigned electronic temperature and no retained
LCFO basis. Compare full support (fraction=1) and automatic .999 support on
the same seeds: impulse 1e-4, dt=.02, 16 steps, fixed ions, MLWF 5/1000/1e-6.
Positive fixed R is currently static-only and is not used here.

```sh
python3 testsuites/benchmark_h2_blocks/run.py --binary /absolute/salmon \
  --output /absolute/new-pilot --pilot --repeat 1
python3 testsuites/benchmark_h2_blocks/run.py --binary /absolute/salmon \
  --output /absolute/new-results --repeat 3
MPLCONFIGDIR=/tmp/salmon-h2-matplotlib python3 \
  testsuites/benchmark_h2_blocks/analyze.py /absolute/new-results/results.json \
  /absolute/report
```

Weak arrays and strong MPI counts match the earlier experiment, but each array
entry is now a core of eight molecules. Weak GS uses one rank per fragment.
Strong measurements refer to native global RT, with the same 128-H2 seed on
MPI1/2/4/8/16. Strong GS is not measured. A weak change in the number of split
axes also changes the padded fragment work; compare the 1-D family separately.
The one-fragment baseline has no buffer by input contract.

RT primary timing is the maximum-rank internal 16-step loop timer; report the
minimum of three sequential runs, median and range. Mode order alternates.
Memory is the median of each job's maximum-rank lifetime peak SALMON RSS,
including seed loading/reconstruction, not a phase-specific or summed node
peak. GS launcher time/RSS is measured once per shape separately from RT.
Other application contention (especially MPI16) remains possible.

Completion requires converged GS, correct charge, LCFO export, all rank exits,
16 finite RT rows and correct impulse/time grid. Energy width and differences
are recorded without relaxing or claiming accuracy certification. Local FFT
counters apply only to active-support updates; they exclude initialization,
communication, and global fallback work. Seed hashes are checked after reuse.
Use --resume only with identical binary/scripts/configuration; incomplete job
folders are preserved and rejected, never silently overwritten.

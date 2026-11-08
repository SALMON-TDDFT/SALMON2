# Internal Fortran W: small Si integration test

This is an integration fixture, not a converged Si spectrum calculation. It uses
an 8³ real-space grid, 2³ k/q meshes, 2³ auxiliary boxes, 1 Ha cutoff, eight GS /
projection states and three RT steps. Four occupied orbitals are propagated.
All SALMON runs use **4 MPI × 12 OpenMP**, sequential nested BLAS, and the
same audited executable SHA256. Supply the KY Si `Si_rps.dat` pseudopotential.

Python only prepares fixtures, independently calculates comparison W, and checks
results. `run-stage.sh` launches SALMON directly; no Python MPI process surrounds
SALMON. Each MPI rank explicitly changes to the stage directory and opens its
own input/output files. `run.log` is the rank-zero transcript; `launcher.log`,
`launcher.stderr`, `stdout.rank-*` and `stderr.rank-*` retain launch evidence.
`FLIB_FASTOMP=FALSE` prevents the Fujitsu fast OpenMP mode from overriding
`OMP_DYNAMIC=FALSE`. Do not run the launcher on a Fugaku login node: use a node allocation.
Include `#PJM -x PJM_LLIO_GFSCACHE=/vol0004` if Python reference libraries use Spack.

For local CTest runs, keep the configured Python interpreter and its dependency
environment available at test time as well as CMake configuration time. The
reference tests need NumPy/SciPy, and the original-H oracle additionally needs
mpmath. If dependencies were installed into a separate directory, preserve its
`PYTHONPATH` when invoking CTest. MPI tests also need permission to open local
communication sockets; a launcher failure before the test starts is not a
numerical comparison result.

From the source tree:

```sh
python3 testsuites/hsex/test_internal_si.py prepare /scratch/internal-si --pseudo /path/Si_rps.dat
bash /scratch/internal-si/run-stage.sh /scratch/internal-si /path/salmon AUDITED_SHA256 gs
python3 testsuites/hsex/test_internal_si.py reference /scratch/internal-si gs
bash /scratch/internal-si/run-stage.sh /scratch/internal-si /path/salmon AUDITED_SHA256 updated
bash /scratch/internal-si/run-stage.sh /scratch/internal-si /path/salmon AUDITED_SHA256 snapshot
python3 testsuites/hsex/test_internal_si.py reference /scratch/internal-si pump
bash /scratch/internal-si/run-stage.sh /scratch/internal-si /path/salmon AUDITED_SHA256 restart
bash /scratch/internal-si/run-stage.sh /scratch/internal-si /path/salmon AUDITED_SHA256 external
python3 testsuites/hsex/test_internal_si.py check /scratch/internal-si
```

The optional `run` CLI wraps this same shell launcher for local convenience.
Each directory/run log must be new; failed evidence is never overwritten.
GS data and checkpoint payloads stay on scratch. Collect only
`integration-check.json`, logs and selected small numerical comparisons into the
canonical calculation-results directory.

The check requires `fixture.json` to identify the prescribed 4 MPI × 12 OMP
case, every executed stage to contain `end SALMON`, and the updated-W log to
report an actual OpenMP team of 12. Updated, snapshot and external trajectories
must contain exactly steps 1–3; restart must contain exactly steps 2–3.

The checks require GS and pump W to agree within `atol=2e-12, rtol=2e-11`.
Pump eigenvalues/populations are checked separately (`atol=1e-7`); the update
itself also validates them. Restart and external-Python-W current comparisons use
`atol=1e-12` atomic units, with no relative tolerance. Bitwise equality is not
required for different OpenMP reduction orders.

The pump test defaults to numerical broadening input η=0.001 Ha (the accepted
W uses η/2=0.0005 Ha). `prepare --eta 0.01` and `check --eta 0.01` reproduce the
previous fixture; the reference reads eta from the fixture manifest. The
production η-convergence threshold 0.02 and head threshold 0.05 are unchanged.
The earlier approximately 1e-10 W discrepancy at small eta is addressed by
original-H iterative refinement, high+low solution storage, and compensated
Hermitian response accumulation; see `docs/hsex-internal-w.md`.

On an oversubscribed local Mac validation, `OMP_WAIT_POLICY=PASSIVE` can avoid
idle OpenMP workers consuming the cores needed by the guarded BLAS region.
This is a local testing setting, not a Fugaku performance recommendation.

To check the already executed local variants without copying GS payloads, `check`
accepts `--updated`, `--updated-log`, `--snapshot`, `--restart`, `--external`, `--gs-reference`, and
`--pump-reference`. The work directory still needs the correct fixture manifest. `--skip-external` explicitly records that comparison as
`NOT RUN`; it is not the complete acceptance procedure above.

For standalone native pipeline/metal tests on Fugaku, set `TMPDIR` to a shared
directory (for example `export TMPDIR="$PWD/tmp-native"; mkdir -p "$TMPDIR"`).
The tests explicitly enter that shared fixture directory on each MPI rank and
write rank-specific logs. Unexpected failures preserve their entire fixture
and print `MPI_TEST_FAILURE_FIXTURE_PRESERVED`. Successful fixture cleanup can
be disabled with `HSEX_KEEP_TEST_FIXTURES=1`.

The Fugaku launcher disables the hardware barrier with the existing SALMON
`FLIB_NOHARDBARRIER=TRUE, FLIB_CNTL_BARRIER_ERR=FALSE` policy for four ranks
per node. This avoids a vendor warning exit after an otherwise completed run;
process failure codes are still checked and are never ignored.

The default fixture also adds a delayed `Acos2` verification probe centered at
0.02 a.u., with width 0.02 a.u., amplitude 0.01 a.u., frequency 1 a.u., and
CEP 0.25. Its vector potential vanishes through the step-1 W update. The
initial impulse remains 1e-5 a.u.; the probe makes the updated-W action
resolvable over the three-step trajectory without changing eta or numerical
acceptance tolerances. This deliberately short probe tests integration only,
not a physical Si spectrum. A weak-kick-only fixture gives an unresolved
current difference at eta=0.001 even over 256 steps; it is not evidence that
W was ignored.

### Small-eta validation on Fugaku (2026-11-05)

The original-H revision passed the 4 MPI x 12 OpenMP native oracle and
fractional-metal pipeline (job 53659189). The final Si comparison completed
in job 53660764, reusing the completed SALMON stages from job 53659191.
The independent compensated dense-LU Python reference used `--mpi` with four
ranks and 12 BLAS threads per rank. Serial and multiprocessing reference jobs
exceeded their respective 10- and 15-minute limits; their evidence is retained.

The Si W maximum absolute difference was 9.87384e-13, with maximum mixed
tolerance ratio 0.327881. External-W and restart current errors were
5.39933e-17 and 1.62604e-19 a.u. Both final jobs exited with code zero.
The fractional-metal W error was 7.36327e-17. No tolerance was relaxed.

For this small Si case, Fortran response preparation took 12.3612 seconds;
the independent dense-LU Python reference took 355.869 seconds. These timings
include different numerical algorithms and do not measure OpenMP scaling
in isolation. The 64-band fractional-metal stress case took 481.1455 seconds
in the Fortran response stage. Larger production meshes require separate
performance and physical convergence measurements.

# Historical Si64 DG fixture adapted to fragment-periodic HSE

The coordinates and pseudopotential come from the actual Si64 divided-SCF fixture
identified through the chat **Task 1を実装計画どおり実行**. `provenance.json` records
the input, coordinate and pseudopotential hashes. The similarly named older
`si64_overlapping_wannier/inputfile_reference.in` actually selects carbon and was
not used.

Preserved: 64 Si atoms, 20.52 bohr cubic cell, 32^3 total grid, 2x2x2 fragments,
6-grid buffers, 400 retained states per fragment, 300 K, Gamma, 8 MPI ranks.
The core grid is 16^3 and the buffered fragment grid is 28^3. The spacing is
0.64125 bohr. All three axes are split, so the current DC implementation requires
Gamma for this geometry.

Changed: PZ to HSE06; unsupported DG controls removed; ScaLAPACK disabled; the
exchange-focused fixture keeps LCFO disabled as in its source. This is an input
for establishing a new full-HSE baseline, not a previously converged HSE result.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 8 /path/to/salmon < inputfile > outputfile
```

Cost matters: one complex orbital array per fragment is
28^3 * 400 * 16 = 140,492,800 bytes, about 134 MiB. There are several simultaneous
arrays, and eight fragments. Each complete ordered pair loop contains 160,000
pairs per fragment before screening. These are array sizes/counts, not measured
peak RAM or timings. Do not treat a one-iteration pilot as a converged baseline.

## Screening study

1. Converge this full-HSE reference and verify localization separately. Preserve
   occupation numbers and distinguish the 400-state DG retained space from the
   density's effective occupied rank.
2. Export localized density factors Q=Psi sqrt(f/2) U for each Gamma fragment.
   Save NPZ arrays `q` with shape (factor,nx,ny,nz), `spacing` in bohr, and `omega`
   in inverse bohr. Do not normalize fractional factors to unit occupations.
3. Run `../pair_screening.py fixture.npz --budgets 0 1e-6 1e-5 1e-4 --output result.json`.
   Budgets are **per-fragment periodic exchange** Ha, not the assembled DC core
   energy and not eV/atom. For each budget report the actual energy/action error,
   conservative omitted-energy bound, and retained ordered/unique pair counts.
4. Check the construction metric's anti-Hermitian part and minimum eigenvalue.
   The diagnostic intentionally does not certify an arbitrary screened-column
   construction as a common linear Fock operator or a variational SCF functional.
   Exact machine precision is not an accuracy target for the approximate physical
   model; these structural checks diagnose whether a proposed ACE construction is
   mathematically applicable.
5. Validate SCF observables and buffer convergence after implementing a consistent
   screened operator. Measure exchange, ACE build/apply, wall time and peak memory
   separately. Reuse the same initial conditions across comparisons.

Eight atoms per core make this a useful correctness/regression case. Because
28^3 is close to the total 32^3 grid, larger supercells are also needed to establish
asymptotic scaling; this fixture alone cannot prove it.

## Initial bounded cost pilot (2026-09-26)

The previous 3f81de4 MPI executable accepted this case with nscf=1 and
hse_mlwf_maxiter=1 (explicitly not a localization/convergence test). All eight
fragments reached their first HSE_WANNIER exchange refresh. The process group was
stopped at the 120-second limit before SCF completion. No Si64 screening rate,
converged HSE energy or speedup is claimed. Before extending the run, assess
retained-state convergence and localization of the occupation-weighted factors;
400 retained DG states need not imply 400 independent occupied exchange sources.

## Frozen-density pilot and exact source reduction (2026-09-26)

`inputfile.n128-pilot` records a bounded 128-state Broyden pilot. It inherited
`alpha_mb=0.75`; `mixrate=0.1` does not control Broyden. It shut down normally at
54 iterations/605 seconds, with residual 1.0376e-4 versus a 1e-7 threshold.
It is **not a converged HSE reference**. The highest retained state had zero
occupation, which alone does not establish retained-state convergence.

Fragment 1 had 89 states with strictly positive occupation. Omitting the other
39 from the exchange source is exact: source-target pairs fall from 128*128 to
89*128 (30.47% fewer) while all 128 action targets and ACE construction states
remain. This is a count reduction, not a measured runtime speedup. The code
performs this reduction dynamically using the common positive-occupation union
over k and resets the gauge when that set changes.

Offline localization with zero occupation cutoff left the density unchanged.
The gradient reached 1.2047e-5 bohr^2, then the line search stopped at iteration
557 before the 3000-iteration cap. It missed the requested 1e-6 criterion; these
are not certified MLWFs. Periodic exchange changed only from
-8.347162802248105 to -8.347162802248437 Ha after the gauge change.

`preliminary-screening.json` stores the frozen-density results and hashes. For
89 factors (7921 ordered pairs):

| Per-fragment error budget (Ha) | Removed ordered pairs | Actual exchange error (Ha) |
|---|---:|---:|
| 0 through 1e-4 | 0 | 0 |
| 1e-3 | 6 (0.0757%) | 4.1065e-5 |
| 1e-2 | 50 (0.6312%) | 4.1428e-4 |

This conservative bound gives little pruning on this particular frozen state.
The screened column metric loses Hermiticity (relative defects 8.29e-4 and
1.92e-3 at the two nonzero pruning levels), so these columns must not be fed
directly into production ACE. Even the unpruned Q-basis metric fails the strict
positive-definiteness diagnostic because tiny finite-temperature occupations
make Q ill conditioned. This does not invalidate the exact Fock operator or
the production ACE metric built on the orthonormal retained targets.

SCF convergence, state-count/buffer convergence, a consistent Hermitian screened
operator and actual timing/scaling remain outstanding. A DC wavefunction restart
attempt was rejected by the existing `yn_restart` guard; it is not a supported
continuation path. A subsequent pilot uses the supported `read_dns_cube` total
density initialization, explicit `alpha_mb=0.1`, and `ncg=20`; it reinitializes
orbitals and must not be described as a wavefunction restart.

The density-seeded pilot ended normally after 43 iterations in
610.3 seconds, still unconverged (residual
1.13426719e-03 versus 1e-7). Its input and summary are stored as
`inputfile.n128-densityseed` and `density-seed-pilot.json`. The seed density cube
is an intermediate checkpoint-derived file; its hash is recorded. This bounded
attempt did not improve convergence and is not the recommended production input.

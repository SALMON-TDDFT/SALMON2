# Embedded MPI TDDFT and E/H coupling validation

Implemented an experimental one-material/one-macro-point-per-group path. Python
owns MPI; the native communication layer duplicates the supplied solver group,
tracks/frees native descendants, and leaves MPI initialized after state close.
TDDFT state arrays are explicitly owned by the embedded context. Existing
initialization_rt/time_evolution_step and odd/even orbital swapping are reused.
The ordinary main_tddft path remains intact. Residual module configuration means
one native state lifetime per process; this is not general multi-state reentrancy.

Native E/H exposes a one-iteration entry through its existing eh_calc. Optional
physical current is injected at the existing current stage before E exchange.
Without the optional current, original numerical ordering is unchanged. LD
histories remain local to field ranks; electronic sites use vacuum constitutive
neighborhoods. Sampling uses native two-edge collocation; transpose deposition
adds only to rank-owned edges. Every macro point represents one field-cell volume.

The explicit clock contract is salmonpy.eh-centered-current/1: native E and saved
centered H are at n*dt; new A/J and raw H are at (n+1/2)*dt. Electronic local time
zero corresponds to field dt/2. Initial feedback J is zero. This is not a claim
of bitwise equivalence to native Weyl multiscale propagation.

## Verification

Local environment: macOS arm64, GNU Fortran, OpenMPI, OpenBLAS; one OpenMP/BLAS
thread per MPI rank. Source and all application-object eigensolver audits pass.
Pinned-source preparation and shared LD/SBE kernel byte-equality checks pass.

- Fresh source distribution -> repaired wheel -> isolated install: 107 Python
  tests pass. New modules, guides, sample and explicit build tool are included.
- Installed-package MPI regressions: three TDDFT internal-decomposition cases,
  vacuum fields on one/two ranks, LD fields, current injection across x/y/z rank
  boundaries, TDDFT coupling on one/two field ranks, and the shipped six-rank
  sample pass (12 cases). An additional two-LD-material + two-TDDFT-point mixed
  calculation passes, for 13 cases total.
- TDDFT with two independent groups: orbital splitting and unequal group sizes
  (one/two ranks) match the serial worker exactly; real-space splitting differs
  by at most 6.938893903907228e-18 a.u. in the short PZ fixture.
- TB-mBJ saved eight-atom Si state, 12^3 real-space and 2^3 k grid: two-rank k-point
  splitting differs from the serial worker by at most 1.6865367127543526e-19 a.u.
  This is an independent four-step inner-MPI comparison, not a full TB-mBJ E/H
  convergence study. Libxc inner-MPI combinations were not tested here.
- E/H one-iteration calls with zero external current reproduce original batch
  E/H observations exactly, both vacuum and background LD.
- Manufactured three-component current verifies -4*pi*dt, both split factors,
  overlapping adjacent-cell deposition and rank-boundary ownership in x/y/z.
- Coupled PZ fixture: 40^3 field grid, dt=.003 a.u., 64 steps, two electronic
  points. Field ranks=1 versus 2 give exactly equal A/E/J. Every electronic current
  matches serial replay of recorded A exactly. A differs from zero-feedback
  propagation by 8.454046411300287e-5 a.u.; the test contains real feedback.
- Mixed calculation: two LD models with one/two poles plus two internally parallel
  TDDFT points, six ranks; serial replay current matches exactly.
- Discrete A/E relation, previous-current lag and explicit half-step timestamps
  are checked at every coupled step. MPI remains usable after native close.
- Stale A and second native initialization are rejected. Review regressions cover
  missing ABI symbols before attachment, collective cleanup, source overlap,
  out-of-band advancement and invalidation after nonfinite native output.
- Ordinary CLI regression: all 19 FDTD material/interface cases pass. Ordinary
  TDDFT before/after output is byte-identical; worker batch/interleaving/rejection/
  recovery tests pass against the prior serial runtime.

The small PZ fixture has synthetic normalized orbitals, not a converged material
GS. Field/material convergence, energy-conservation benchmarking, production
scaling and multi-node shared-filesystem operation remain unverified. This path
rejects restart, nonzero initial fields/impulse, rotations, unequal macro volumes,
subcycling and TTM. A fresh MPI job is required for another native lifetime.

## Reproduction and artifacts

Use COUPLED.md / COUPLED.en.md for user setup. Explicit build tool:
`tools/build_embedded_tddft.py`; bounded regression runner:
`tools/run_coupled_mpi_tests.py`. Linux CI includes these steps but has not been
pushed or executed for the current changes; local results are not CI evidence.

Native library SHA256:
`1e9abfa895dc216399377c164c0cbefa606c436c7e7eee920f4aa8901fbdae2a`.
Runtime patch SHA256:
`14b82f79f73d589f9be92327e70c183d0a72a63f9a8a06812a6e0a66b6e744f7`.

Historical local evidence, not an installation default:
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/native-eh-tddft-mpi`.
Contains run logs, native build/source hashes, current comparisons, ordinary-CLI
reports, sample outputs and the wheel checksum.

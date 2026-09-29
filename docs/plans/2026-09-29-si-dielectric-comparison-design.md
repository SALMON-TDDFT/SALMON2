# Si dielectric comparison: approved design and execution plan

> Execute inline in the current feature checkout. User approved the comparison on
> 2026-09-29 and explicitly removed all zero-field runs and current subtraction.

**Goal:** Compare PBE, HSE06, PBE0 and PBEh(40) dielectric response and cost in
ordinary Si Nk=4^3, then 4x1x1 DC-initialized Gamma Si and a matching ordinary
Gamma reference. rVV10 is off.

**Architecture:** Reuse SALMON GS/native RT and DC/LCFO/native RT. Generate separate
immutable run folders and retain inputs, executable/input hashes, exit codes,
convergence evidence, timings and rank RSS. Analyze the impulse current directly;
never synthesize a zero-field subtraction. Keep raw DC response clearly labeled.

**Tech stack:** Existing GNU MPI/HSE/Libxc/ScaLAPACK build; Python/NumPy/Matplotlib.

## Physical conditions

- Diamond conventional cell a=10.26 bohr, eight atoms; same Si_rps.dat throughout.
- Start with 16^3 real-space points per conventional cell (h=.64125 bohr).
- Stage 1: full uniform 4^3 k mesh, no symmetry reduction, fixed nuclei.
- Stage 2: 4x1x1 conventional cells, 32 atoms, Gamma. Four DC cores with x buffers;
  y/z are unsplit periodic directions. Start x buffer at half a conventional cell
  and test a larger buffer/retained-state space before interpreting DC accuracy.
- Matching ordinary 32-atom Gamma GS/RT isolates DC effects. It is not equivalent
  to stage 1, which samples 4^3 k points rather than 4x1x1.
- Hybrid exchange uses full wavefunction support for baseline spectra. Compare
  99.9% support against full support in short Gamma tests before adopting it.
- PBE0/PBEh use the existing spherical Coulomb kernel. Record its effective
  radius: default half shortest Born-von-Karman supercell length (20.52 bohr in
  stage 1, 5.13 bohr in Gamma). Radius dependence is separate from WF truncation;
  cross-stage differences must not be attributed solely to DC or k sampling.
- DC electronic temperature 300 K, common chemical potential, no localization
  during fragment SCF. RT nuclei fixed, occupations fixed by the supported route.
- x impulse 1e-4 au, transverse, no zero-field run/subtraction. T=350 au
  (8.466 fs), nominal Fourier spacing .4885 eV. Initial dt=.08 au; compare to .04
  over equal short duration before launching 4375-step production. Refine if needed.
- No broadening disguised as physical lifetime; cubic finite-time window.
- Report first peaks and spectrum shape, not a KS gap as an optical excitation.

## Tasks and acceptance

1. Generate and inspect stage-1 input files and a stable executable snapshot.
   Test atom counts/cell/k grids, no zero-field branches, GS/RT functional and
   Coulomb settings match. Never overwrite completed runs.
2. Run fresh PBE and hybrid GS, require SALMON end marker and explicit convergence.
   Use density-converged PBE pre-SCF for hybrids. Failure stops the queue.
3. Run dt=.08/16 steps and .04/32 steps, compare current on common times, energy
   and electron norm. Require finite values and <=1% relative current discrepancy
   (absolute floor 1e-9 au); no assumption about stability of the whole run.
4. Run stage-1 impulse-only production and analyze eps_xx with the common window.
   Validate row counts/times and completion before calling results final.
5. Prepare and verify stage-2 DC/LCFO export and ordinary Gamma control. Check the
   existing HSE ordinary Gamma restart path before using it; do not silently route
   a nominal ordinary reference through DC. Validate DC buffer and local support.
6. Compare spectra and separately GS iterations/time, RT loop seconds/step,
   rank maximum RSS and sum of rank peaks (not simultaneous aggregate RSS), EXX
   refresh timings when available. Record actual parallel layouts and algorithm
   differences; multi-k global hybrid refresh currently runs on k-root.
7. Save Japanese results note with completed/pending cases and convergence limits.

## Scheduling and verification

Single measurement job at a time. MPI4/OMP1 for initial k-parallel trials,
BLAS1. Wait for the existing H2 PBE0 job before uncontended timings; never cancel
or rebuild its executable. Resource limits are not physical approximation knobs.
No MPI16 runs. No automatic future notification/automation has been requested.

## Execution ledger

- Approved: four functionals, three geometries/routes; user removed zero-field
  runs and subtraction. Existing H2 jobs are preserved.
- Source baseline 801fd01c; GNU build already completed. NVHPC is not used for
  these local CPU measurements.
- Code inspection: global-hybrid multi-k refresh gathers to k-root and uses a
  BvK-supercell Coulomb cutoff; timings are not comparable algorithm-for-algorithm
  with the distributed HSE k kernel. Preserve this caveat in the final note.

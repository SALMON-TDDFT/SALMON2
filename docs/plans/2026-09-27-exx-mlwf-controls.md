# EXX MLWF controls and local support implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Share MLWF controls and a spherical integration-support radius across HSE and PBEh.

**Architecture:** Canonical `exx_mlwf_interval/maxiter/tolerance` retain the
existing defaults (10/200/1e-6); old `hse_mlwf_*` inputs are compatibility aliases.
Explicit conflicting new/old values are rejected. `exx_mlwf_radius` is an input
length, converted to bohr; 0 retains full support. Reuse `lcfo_wf_plan_init`'s
periodic sphere/protected-center mask on occupation-weighted Wannier sources.
Pair formation and accumulation then use the existing compact support path.
This does not shrink the FFT domain or change the Coulomb kernel radius.

**Tech Stack:** Fortran, MPI, FFTW, Python unittest.

User approved this design, including old-name compatibility and withholding
finite-radius MD until force consistency is established. A finite source mask
is an explicit, gauge-dependent approximation; it does not establish a
variational SCF energy functional or correct ionic forces. Permit it only in
static DFT, and disallow checkpoint restarts and legacy snapshot export with
positive radius. Keep the separate legacy LCFO-RT radius unchanged (bohr).

### 1. Input compatibility
- Add compiled input tests in `developer_tests/653_functional/test_exx_inputs.py`:
  new/old equivalence, matching/conflicting dual assignments, validation and
  length-unit conversion, finite-radius MD/restart rejection.
- Replace internal MLWF control references in `src/io/salmon_global.f90`,
  `src/io/inputoutput.f90`, `src/xc/hse_native.f90`, `src/xc/lcfo_rt_wannier.f90`.
- Resolve sentinels after broadcast, before input validation; log canonical values.

### 2. Reuse spherical support
- Add a compiled probe for radius-zero/full limit, periodic-boundary masks,
  protected delocalized factors, discarded norm, Hermiticity, compact/dense
  equivalence for HSE and PBEh kernels.
- Use `lcfo_wf_support` in `src/xc/hse_wannier.f90`; compute periodic moments
  and protect any factor with center reliability below 0.1 on any axis.
- Mask both appearances of the source through the existing `op%source` action;
  do not renormalize or mask target states. Report removed norm and protected
  count at native refresh. Large enough radius must reproduce full support.
- Update standalone test compile dependencies for the shared support module.

### 3. Integration and documentation
- Build with/without HSE. Run new probes, MLWF and PBEh/DC regressions.
- Check finite-radius native SCF on a small periodic fixture and radius-zero
  parity. Update samples to canonical names and document radius convergence.
- Obtain independent review, address findings, and commit on the existing branch.

Ruling: finite-radius initialization must be independent of k-rank layout.
The new two/four-rank regression found a 0.112 eV difference. Existing Gaussian
initialization visits every global k but seeds by the first locally owned k;
that changes the initial occupied subspace and localization trajectory when
ranks change. For positive EXX radius only, seed the global Gaussian traversal
identically, and reseed random wavefunctions by physical k. Radius-zero behavior
is preserved. This touches `src/gs/init_gs.f90` beyond the original input/mask
files and is required for a reproducible finite-radius operator trajectory.


Ruling: reject a finite-support SCF result if the last MLWF minimization failed.
After fixing the initial seed, the DC fixture retained a 2.9e-4 eV difference
between layouts while both MLWF gradients remained far above tolerance. The
ordinary density convergence test cannot certify this gauge-dependent model.
Keep the requested mask during iteration, track the last spread-minimization
status through transport-only refreshes, and check it collectively at the SCF
exit before accepting a result. Radius-zero and masks that discard no norm are
exempt. The MPI regression now checks the same initial finite-mask update for
Gaussian/random seeds and rejects the known unconverged localization, rather
than treating those final energies as valid. A separate conventional test checks
converged finite-radius HSE/PBEh SCF. This is not a converged finite-radius DC
validation; its gauge/radius convergence remains to be established.

## Verification and review ledger

- Canonical and legacy inputs produce identical energies; matching dual values
  are accepted and conflicting values rejected. Input-length conversion, invalid
  controls, finite-radius MD/restart/snapshot guards are covered.
- Compiled sphere probe checks periodic boundary wrapping, ambiguous-center
  protection, norm diagnostics, compact/dense equivalence and Hermiticity for
  HSE/PBEh and one/two k points. Zero/large support radii reproduce full support.
- Conventional finite-radius HSE and PBEh+rVV10 converge in the small hydrogen
  fixture, with the new final-localization guard enabled. No finite-radius force
  or production-water convergence claim is made.
- MPI regression checks Gaussian and random initialization at the same first
  truncated DC update on two/four ranks, and rejects unconverged localization.
  The existing full-support one-fragment and multi-rank DC checks remain intact.
- HSE-enabled and HSE-disabled builds pass. Legacy exact-pair, LCFO transport,
  LCFO sphere and sphere-norm probes pass after updating the internal test stub
  names. Old input fixtures remain to exercise compatibility aliases.
- Independent review caught the LCFO test stub rename and stale localization
  status after failed polar transport. Both were fixed. A compiled regression
  first reproduced the stale status, then passed after invalidating the status
  when a new gauge is seeded or the overlap is lost.

# Stable Libxc HSE small-gradient evaluation (2026-11-18)

Scope: fast-PR HSE only, CPU scalar periodic full-q, symmetry opt-in with one
internal k rank. Local Apple Silicon / GNU Fortran15 / Libxc5.2.3. Synthetic8³
states, equivalent full8/reduced4 k points, C2z. No physical-convergence claim.
No GitHub Actions used.

## Cause and change

The original phase0 failure is retained. Libxc428 combines WPBEh(0),
-0.25 WPBEh(omega), and PBE correlation. WPBEh(0) differentiated separate
E1/log terms with cancelling1/sigma contributions. This amplified a density
roundoff difference into a local potential difference about0.016 and current
error about5.7e-10; refining dt did not remove it.

The complete hybrid_semilocal source override combines the cancellation before
energy differentiation, using u=s² for s<=1e-4. The E1 remainder is below double
precision in that region. The screened derivative uses its analytic quadratic
limit below roundoff in u; correlation remains Libxc. Spin scaling, fixed-core
NLCC, and Libxc's existing low-density screening are retained. Points outside
the branch retain the original Libxc output bytes. Pinned sources and old
binaries are unchanged. See runtime/NOTICE.md for MPL-2.0 provenance.

Both the full-k serial worker and host-MPI library were explicitly rebuilt in
new external directories. MPI library SHA256:
`9214851441d9be9ccc11783bbcd292defa44f04c0975e95140f58e8656cb789f`.
Symmetry protocol1 alone does not identify the new numerical evaluation.

## Measured checks

- Native regression fails with STOP2 on the old library; the MPI comparison also
  fails on the old library, without relaxing its2e-12 absolute current tolerance.
-240 density/screening/gradient cases:126 compare against independent80-digit
  energy differences;72 ordinary-gradient and42 low-density-screening cases
  retain raw Libxc output bytes. Max relative errors (energy,vrho,vsigma):
  `5.46e-16,6.83e-16,2.48e-15`. Absolute vsigma error64 occurs only where the
  low-density coefficient is of order1e17; report both relative and absolute errors.
- Polarized density and gradient energy differences, nonzero fixed-core
  density/gradient equivalence and single-spin threshold crossings pass.
- Original phase0 recorded local FDTD drive, dt=.003,32 steps:
  OMP1 max full/reduced current difference `4.2012834183813297e-19`, energy
  difference `8.881784197001252e-16`;
  OMP2 current `4.336808689942018e-19`, energy `4.440892098500626e-16`.
  Each host rank owns an independent COMM_SELF solver group, not internal-k MPI2.
- Original phase0 mixed field + two local LD materials + SBE + HSE passes the
  normal32-step FDTD run. Serial full-k replay max current difference8.95e-19.
  Injected late/deep/symmetry errors preserve originals and produce exactly the
  same final A/E/j arrays after retry.
- Four-rank protected restore fails at18 after saving16, preserves all original
  field/LD/electronic bytes and matches three future retry steps exactly;
  restored steps replayed:0. Fresh-process16→32 direct resume matches all
  continuous A/E/j and reference arrays exactly.
- Ordinary five/seven-rank PBE, HSE Taylor, HSE PT-CN and combined regressions
  pass with unchanged tolerances. PBE arrays remain byte equal. Max HSE/mixed
  current difference from prior binaries7.94e-17, A/c4.47e-19, E9.11e-18.
  HSE results are roundoff-equivalent, not claimed bit identical.

Source/installed Python suites and source-distribution-to-wheel/native-hash
checks are recorded in the external evidence bundle. Native numerical tests
must be rerun on other platforms; no cross-platform execution is claimed.

## Reproduction

Prepare a new complete source with `tools/prepare_hse_rt_state.py NEW_SOURCE
--explicit-routing`, then explicitly build matching MPI library/modules and
serial worker. Keep the former binaries for the failing regression. Compile
`tests/native/test_hse_low_gradient.f90` against the new hybrid_semilocal module
and library, and run:

```sh
python tools/test_hse_low_gradient.py BUILD_RECORD NATIVE_PROBE NEW_OUTPUT
mpiexec -n 2 python tools/test_fast_hse_low_gradient_mpi.py \
  BUILD_RECORD FULL_K_BUNDLE REDUCED_K_BUNDLE NEW_OUTPUT \
  --drive ORIGINAL_ARRAYS_NPZ --point-index 1
```

The default fixture phase remains0; no GS perturbation or field projection is
used. Raw `probe_hse_low_gradient.py` still detects sensitivity inside Libxc;
this diagnostic is distinct from the corrected wrapper verification.

Evidence directory (local execution, not a portable distribution URL):
`/Users/otobetoshihito/.codex/outputs/salmonpy-hse-low-gradient-20261118/`.
It retains source archives, build receipts/binaries, original drive, failed and
successful logs, selected rollback witnesses, direct checkpoint and comparison
arrays. Reduced internal-k MPI, broader initialization and converged physical
benchmarks remain separate work.

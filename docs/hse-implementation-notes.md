# HSE/SOC implementation map

Native hybrid orchestration is in `src/xc/exx_native.f90`; the `exx_*` modules
separate exchange kernels, ACE/history, symmetry reconstruction, orbital/data
layouts and worker communication. `hybrid_semilocal.f90` supplies the compatible
semilocal remainder. See [module names](exx-module-names-ja.md).

The shared scalar/SOC worker consumes symmetry-reduced Bloch orbitals. Automatic
FFT layout selection agrees trial errors, mode and worst-rank timings before
creating plans. Uneven strips retain exact coverage, quotient/remainder balancing
within each row, and bounded peak density storage.

KY pseudopotential extensions and noncollinear density/current treatment are
separate from hybrid exchange. The PP conversion utility is under
`utility/crps_soc_nlcc/`; its format documentation and generation inputs are
retained, while material-specific measurements and development probes are not
part of the proposed regression suite.

Real-space RT and checkpoint metadata retain their own compatibility and layout
checks. Fixed-LCFO-subspace RT has been removed; see [LCFO state preparation](inputs/lcfo-rt-development.md).
For accepted input controls and restrictions use [HSE inputs](inputs/hse.md) and
[feature selection](inputs/hse-feature-selection.md).

This map describes the submitted code, not a claim that every optional route is
validated. See [review scope](PR_REVIEW.md) and [validation](PR_VALIDATION.md).

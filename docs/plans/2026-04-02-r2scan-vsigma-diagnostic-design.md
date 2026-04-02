# r2SCAN vsigma Diagnostic Design

**Goal:** Add a diagnostic-only alternative representation for the r2SCAN density-gradient stress term so SALMON can compare the current `rdedd` route against a `vsigma` route on the same grid and with the same SCF state.

**Context:** Current `output_level='C'` analysis shows the residual XC-stress anisotropy is almost entirely in the `XC-grad` sector, while `XC-local` is isotropic and `XC-tau` is negligible. The current implementation stores `rdedd = -vgrad * grad(rho) / |grad(rho)|`, which introduces an explicit division by `|grad(rho)|`.

**Approach:**
- Keep the existing `rdedd` payload and stress assembly unchanged for the main path.
- Add a diagnostic scalar payload `vsigma` for builtin `r2scan` only.
- In `stress_fd_detail='high'`, compute an alternate gradient stress tensor from `vsigma` and the already available `grho_local`, and write it alongside the existing `XC-grad` decomposition.

**Why `vsigma` first instead of `h` and `vsigma` together:**
- `h` and `vsigma` are analytically equivalent when evaluated on the same `grho`.
- The real comparison we need first is “division by `|grad rho|`” vs “no division by `|grad rho|`”.
- `vsigma` is the minimal addition because it is a scalar field payload and reuses the current stress-side `grho_local`.

**Non-goals:**
- Do not change the total stress used by the code path yet.
- Do not remove the current `rdedd` route.
- Do not attempt a physics fix in the same step.

**Success criteria:**
- `output_level='C'` prints both `XC-grad-rdedd` and `XC-grad-vsigma`.
- Existing source tests still pass.
- Re-running the `18/22/26` MPI cases gives directly comparable anisotropy norms for the two gradient representations.

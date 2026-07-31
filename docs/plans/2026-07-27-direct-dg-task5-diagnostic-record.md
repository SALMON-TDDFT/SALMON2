# Direct-DG Task 5 diagnostic record

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Status

The direct real-space DG Task 5 gate did **not** pass. No accepted Task 5
checkpoint was published. The work described here is diagnostic evidence
only and is not an implementation prerequisite for the overlapping-Wannier
route.

The inspected worktree was at plan HEAD
`c323411361f6d6c8d8ca72cb1135d9048f3c63ab`. Its complete tracked
uncommitted diff was quarantined as
`/tmp/wpw_task5_uncommitted.patch` (SHA-256
`9ac4d1532c12da281a119bb1ddc6660251f977edcb542b795f915c832afe9dac`).
The untracked direct-DG gate checker was quarantined as
`/tmp/check_si64_dg_gate.task5-quarantine.py` (SHA-256
`de6de208780a3c0b013fa67cf8b8ec622fff2e0d283b391ae744d3e0c88c9927`).

## Inspected changes and disposition

The diff changed direct-SIPG face subtraction, projector acceptance,
direct-DG CG and Rayleigh-Ritz diagnostics, continuation ordering,
local-basis diagnostics, two direct-DG tests, and two direct-DG input
tolerances. None defines center-owned overlapping-Wannier metadata,
unique-core quadrature, or a nonorthogonal coefficient solve. Consequently
two minimal prerequisites were retained in production source:

- the `dcdft.f90` fragment state-count correction removes a fixed per-atom
  overwrite so that the namelist `nstate_frag` remains the authoritative
  candidate-window size;
- the production-state type retains `raw_basis_count` and
  `hamiltonian_operator_fingerprint`, which are already consumed by the
  committed checkpoint module, and the coefficient solver retains its
  optional iteration/residual outputs already used by committed
  `main_dft.f90`; these are required for a clean branch build.

No direct-SIPG diagnostic behavior was retained. The new route otherwise
starts from the committed direct-DG history and uses a distinct default-off
flag.

In particular:

- `dg_dc_direct_sipg.f90` subtracted a face action evaluated from the local
  buffer reference. Its focused test showed that a continuous local-buffer
  continuation should receive zero added face correction. This diagnoses
  the original full-face implementation as double counting the already
  present local-buffer kinetic continuation.
- `conjugate_gradient.f90`, `scf_iteration.f90`, and
  `subspace_diagonalization.f90` measured separate \(H_\mathrm{DC}\) and
  \(\lambda H_\mathrm{DG}\) norms, measured reduced-Hamiltonian
  Hermiticity, and applied the frozen direct-DG Hamiltonian during
  Rayleigh-Ritz.
- Projector changes separated residual and escape tolerances and reduced
  the required rank to the occupied fragment rank. They exposed, but did
  not eliminate, residual and tail-escape failures.
- Apart from the minimal prerequisites described above, the
  remaining edits were direct-DG diagnostic plumbing or tests. They are
  preserved only in the quarantined artifacts above.

## Raw evidence

The logs are under:

`/tmp/si64-task5-dg-first/samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac`

Key files:

| Log | SHA-256 | Evidence |
| --- | --- | --- |
| `run_dg_hnorm_fast.log` | `6c5e5b57a1f118cf6e8a74106855153fde81557a3f15062e02d8a64022e1232c` | The original full-face action was much larger than the DC action: at \(\lambda=0.125\), \(\|H_\mathrm{DC}\psi\|=3.0985\times10^2\) and \(\|\lambda H_\mathrm{DG}\psi\|=1.5566\times10^3\). The energy jumped from about \(-9.943\times10^3\) to \(-1.779\times10^4\). |
| `run_dg_lambda1e4_fullh.log` | `668fcb97339fb208eec2c069c8116716f33050140145ae5268e664599c7bcaf2` | At \(\lambda=10^{-4}\), the direct term norm was already about \(1.238\) against a DC norm of about \(309.8\); reduced Hermiticity defects varied and reached \(4.67\times10^{-5}\). |
| `run_dg_lambda005_fullh.log` | `4f1600239912cb97ef75d6c323de4cf91afa125d9af2fc97a60275cffc20c676` | At \(\lambda=5\times10^{-3}\), the direct term norm was \(6.2264\times10^1\), the reduced Hermiticity defect was \(3.8536\times10^{-2}\), and the SCF energy was displaced. |
| `run_dg_lambda1e4_handoff1e2_delta_tolsep.log` | `cb9a59551849c8ee9546b5e7da4953c9d7f9a296027eb0fb31c01c20028bf50e` | Local-buffer reference subtraction allowed stable continuation through \(\lambda=2.5749\times10^{-2}\). Across generations 1–12 the reduced Hermiticity defect increased from \(1.03\times10^{-9}\) to \(1.21\times10^{-4}\). The subsequent projector gate failed with forward/reverse residuals about \(1.15\times10^{-2}\) and \(9.92\times10^{-3}\), and escape about \(6.03\times10^{-2}\). |
| `run_dg_first_numeric.log` | `e7057dcc5c6a0a2291062308a01f97fd59c732f1157cf719df6fa346191b410f` | The initial projector gate retained rank 361 instead of required rank 400; residuals were \(9.81\times10^{-3}\) and \(9.62\times10^{-3}\), with escape \(5.16\times10^{-2}\). |
| `run_dg_projector_rank1e14.log` | `f7f9f8664325bc9651a541b0f7898b30e7dade14f1982a6bf2cccc2d353ecf7b` | Forcing rank 400 did not cure the physical gate: residuals remained about \(8.28\times10^{-3}\), escape \(4.35\times10^{-2}\), and the configured residual limit was exceeded. |

These results establish failure evidence, not acceptance. No log reached
full \(\lambda=1\), a verified unmixed fixed point, all projector gates,
and successful route-specific checkpoint publication.

## Task 1 verification evidence

The clean-first overlay root was `/tmp/ow-task1-overlay.ivUtqd`. Its base
was the committed plan HEAD, followed by only the Task 1 files and the
enumerated compile-only parent prerequisites recorded with hashes in
`parent-prerequisite-sha256.txt`. The prerequisites supply the parent WPW
API expected by the existing `lcfo_flux.f90`; none was copied into this
branch. The context overlay is the parent real-space metric context plus
the branch's existing H-epsilon-S/global-correction fields.

Configuration:

```text
cmake -S /tmp/ow-task1-overlay.ivUtqd/source \
  -B /tmp/ow-task1-overlay.ivUtqd/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
  -DCMAKE_Fortran_FLAGS=-fallow-invalid-boz \
  -DUSE_MPI=ON -DUSE_SCALAPACK=ON -DUSE_EIGENEXA=ON \
  -DUSE_WANNIER90=OFF
```

EigenExa 2.4b first required its dependency target to be serialized because
its Makefile has a module-order race under parallel build. The required
full command then completed:

```text
cmake --build /tmp/ow-task1-overlay.ivUtqd/build --target salmon -j4
[100%] Built target salmon
```

A fresh two-rank H2 normal-DC fixture, with the new flag absent, produced
both `start DC-LCFO`/`end DC-LCFO` and EigenExa initialization/diagonalize
markers. It contained no overlapping-Wannier marker. The raw log is
`/tmp/ow-task1-normal-dc.DMhjAQ/run.log`, SHA-256
`ac5d68f2532f2cb3a98d2fbdf156827afdcc675d3813fe42a229badd38dfeb05`.

The matching valid PZ/Gamma/non-SOI/DC fixture with the new flag enabled,
LCFO and EigenExa disabled, and normal restart publication disabled stopped
collectively with exactly `overlapping Wannier route: construction not
implemented`. Its log contains no LCFO, EigenExa, WPW, or checkpoint marker:
`/tmp/ow-task1-normal-dc.DMhjAQ/route_stop.log`, SHA-256
`7c7bd2935172684b8ae1eac662110d5b634a7c158554422c5fbaf010ee0e5af8`.

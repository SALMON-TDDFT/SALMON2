# Native HSE06 NLCC and SOC implementation plan

Approved scope: ordinary periodic fixed-ion non-DC HSE06 ground state and linear-response RT. Keep production Fugaku binary unchanged. Base 0f1bddf3c; branch codex/hse-nlcc-soc.

1. Establish an isolated build and baseline semilocal/exchange tests.
2. Add failing compiled NLCC tests: total-density Libxc agreement and finite-difference derivatives with fixed core. Implement total density and total gradient in hybrid semilocal energy/potential; preserve valence-only Fock sources. Verify NLCC GS and RT.
3. Add failing spinor exchange and spin semilocal tests: explicit Fock reference, off-diagonal spin density, SU(2) covariance, closed-shell limit and ACE interpolation. Implement full-spinor exchange using existing distributed density tiles and ACE over both spin components. Reuse existing noncollinear XC rotation and SOC Hamiltonian/current.
4. Enable supported noncollinear HSE GS and Taylor4 RT only. Preserve guards for DC, spatial/orbital decomposition, learned ACE, symmetry reduction and unsupported occupations.
5. Validate combined SOC+NLCC GS convergence, current and finite RT propagation, direct-vs-ACE, and original scalar regression. Record exact numerical limits. Review whole change before delivery.

Ledger:
- Existing SALMON2 checkout contains unrelated uncommitted changes; left untouched. Native worktree tool failed because parent repository metadata is unavailable; created worktree from clean pinned local source using git. No production binary altered.

- NLCC RED: compiled probe rejected absent total-density API. GREEN: Libxc equality and fixed-core density/gradient finite differences pass.
- SOC RED: absent spinor exchange module; polarized NLCC evaluator rejected spin dimension two. GREEN: explicit off-diagonal Fock, SU(2), ACE and two-rank k-distribution tests pass.
- Closed-shell integration RED: legacy spin rotation formed 0/0. Fixed with ordered Hermitian density eigenbasis via atan2 and an identity at zero magnetization.
- SOC energy correction: scalar Hpsi-inferred nonlocal energy includes Fock; direct pseudo_so energy does not. Removed the scalar subtraction from SOC only.
- Review found midpoint XC matrix was missing: time-order fixture RED ratio 2.0634; GREEN after snapshot/average ratio 5.4616. Restored scalar packed-block action index and restricted SOC q reduction.
- No unsupported spinor checkpoints, approximate support paths or hybrid NLCC ionic motion are enabled.
- Ruling: keep original scalar and j-dependent pseudopotential operators; do not assert their equivalence from an ad hoc SO-coefficient-zero fixture. Unit spinor exchange and equal-spin semilocal limits are checked independently. Cost if a legacy projector discrepancy matters: a separate pseudopotential-operator investigation is needed before quantitative scalar/SOC total-energy comparisons.

- Final audit: include core gradients during PBE pre-SCF even while exx_enabled() is false; select the hybrid functional type for the gradient path, leaving Fock disabled until the target-HSE switch.

- Compatibility ruling: consistently ordered spin eigenvalues are selected for HSE06 GGA only. Preserve the original nonhybrid spin chart/arithmetic away from singular axes; zero cases remain defined. This avoids propagating a new density-channel labeling into already-supported nonhybrid SOC.

- Final OpenMP audit: make new density-rotation temporaries tmp/tmp1 private in both CPU/OpenACC loop clauses. Shared temporaries caused schedule-dependent legacy Au reference drift; check multi-thread GS/RT after correction.

- Final verification: 8 compiled/unit/MPI tests PASS; threaded HSE Si GS/RT suite PASS; legacy Au GS and 1000-step RT reference verification PASS with MPI4 x OMP2. Preserve isolated branch/worktree; no merge, deployment, production binary replacement or remote job action.

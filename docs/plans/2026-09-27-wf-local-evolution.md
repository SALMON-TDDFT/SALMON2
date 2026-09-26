# Local WF construction and evolution implementation plan

**Goal:** Implement the user-approved sequence: exact spherical reconstruction, independent U transport cadence, then localized WF evolution with electromagnetic phase handling.

**Architecture:** Keep the validated dense path as a reference. Stage 1 groups grid rows by identical retained-WF masks and removes exactly zero basis columns in each group; dense BLAS operates only on retained products. Cache immutable basis blocks. Stage 2 schedules U transport by physical step independently of ACE, preserving predictor rollback and mandatory impulse handling. Stage 3 requires a separate propagation radius, overlap checks and phase-aware periodic transport; a compact exchange mask must never silently become a propagation truncation.

**Tech Stack:** Fortran, BLAS, MPI, existing native Diamond GS and regression probes.

## Execution order / ledger

- [x] Stage 1 RED: extend active-WF probe to require cached sparse work counts and agreement with dense mask for complex nonorthogonal, exact-zero, tiny nonzero, empty, full and protected cases.
- [x] Stage 1 GREEN: implement reusable row groups and cached nonzero basis blocks in lcfo_wf_support.f90. Group equal WF-mask and exact basis-support patterns; unmasked dense-basis data naturally form a single BLAS block. No magnitude threshold. Integrate explicit one-time preparation in lcfo_rt_wannier.f90; leave existing legacy API available for reference tests.
- [x] Stage 1 verify: active WF, transport, sphere and native regressions; clean release build. One numerical job at a time. Compare C64/C128 against frozen19b9333c and record raw timings plus numerical differences.
- [x] Stage 2: add independent physical-step U cadence with default1, master validation/broadcast, predictor/rollback tests and impulse first-step refresh. Validate intervals2/4 against1; record errors separately from exact stage1.
- [ ] Stage 3: design the full localized coefficient propagation against actual LCFO/ACE action interfaces, separating R_prop from R_int. Test phase-only gauge covariance and periodic boundary crossings, then dense-limit propagation, norm/orthogonality and finite-radius convergence before enabling any approximation.

Ruling: Existing dedicated branch dc-hse-mlwf-ace is clean and reused. No additional checkout is needed. User explicitly approved ordered implementation; proceed without repeated approval questions. Complete each stage's validation before changing the numerical approximation in the next. The precise stage3 integrator must be resolved from source, not assumed from the earlier architectural discussion.

Stage1 ledger: missing-kernel API RED observed; exact reconstruction/transport/sphere and native MPI regression GREEN. Reviewer found no production blocker; isolated tiny-term relative test and all-zero row added. Diamond timings in progress; block-count/performance caveat will be evaluated before making performance claims.

Stage2 ledger: initial held-frame reference caused secular spreading; reviewer identified it, 101-step stationary-density test reproduced it, and separate last-transported anchor fixes it. Added nontrivial cached-U acceptance assertion after rollback. Full/finite support OMP1/2/4, existing transport and native MPI regression all pass. Six sequential Diamond runs finished: U1/2/4 polar counts34/18/10; C64 source .8276/.7332/.6889s, C1281.6949/1.3483/1.0505s. E_inf relative to U1 is ~.00974% at2 and .02863% at4 for both sizes (16steps only). Default remains1; one timing sample per case, no definitive performance claim. Raw results and notebook stored under diamond-u-cadence-anchor and u-cadence-results.json.

## Stage3 interface investigation (not implemented)

Native LCFO RT still propagates real-space s_orbital through time_evolution_step/taylor and projects its Hamiltonian action through lcfo_hse_add_action. Initialization reconstructs from DC coefficients in init_conventional_from_dcdft_complex. Merely rotating the initial orbitals to WFs and bypassing C U is insufficient: stationary occupied energy differences can spread directly propagated WFs, the same phenomenon caught in stage2. The next implementation must include the occupied-gauge evolution (e.g. a verified parallel-transport equation) as well as sparse coefficient storage; renaming native orbitals is not completion of stage3.

A spatial electromagnetic phase is also not just an occupied U rotation. With finite fixed LCFO basis B, generally exp(-i A.r) B lies outside span(B). A projected phase B† exp(-i A.r) B is not generally unitary. Before using it as an acceleration, test projection leakage and metric covariance; do not replace it by an assumed unitary diagonal phase in coefficient space. A field-dressed moving basis needs its time-derivative connection and consistent periodic hopping/nonlocal/ACE action. Keep R_prop separate from R_int and test the dense/untruncated limit before a finite propagation radius.

Remaining stage3 gates: choose and verify direct-F gauge evolution against dense native density/current; implement sparse row/column ownership and Hamiltonian/ACE action; preserve overlap (or consistently use its inverse for nonorthogonal WFs); handle periodic phase and boundary crossing; then run R_prop convergence and longer optical-response tests. No direct-F production path or exp(i A.r) compensation has been enabled by stages1/2.

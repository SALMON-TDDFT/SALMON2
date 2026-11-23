# HSE occupation input validation — 2026-11-16

## Change and limits

Nine exx_native routines resolve nelec, temperature and nstate from the initialized system occupation snapshot. Arithmetic, branch order and update order are retained. An uninitialized system uses the legacy getter fallback. No new public Python API or live setter was added.

Independent native owners remain limited to periodic, scalar unpolarized, fixed integer occupation, serial CPU full-support fast HSE. Spin/DC/Wannier/spatial and other shared controls are not generalized. The outer TDDFT settings_match guard remains because non-HSE consumers still use shared values. Different electron counts in this native HSE probe do not establish arbitrary different-material full TDDFT creation. The nstate-dependent Wannier branch is routed but not behaviorally covered by the two-owner test (yn_hse_wannier=n).

## Fresh local evidence

Apple Silicon, GNU Fortran, OpenBLAS, Libxc 5.2.3; build workers 2, numerical tests OMP 1/2 with BLAS 1. No GitHub Actions, push or physical-convergence campaign.

- Final owner regression fails at binding (STOP 1) against the recovered pre-change build, and passes against the changed build, with ACE and full-action modes at OMP 1 and 2.
- A stores nelec=4/temp=-1/nstate=2; B stores nelec=2/temp=-2/nstate=3 and an empty second band. Globals are poisoned to 99/.5/0 before binding. Actions and exchange energies are byte-identical between independent subprocesses and interleaved runs; fresh reference owners, freeze flags, close/reopen and surviving-state invariance also pass.
- The former disjoint delta test orbitals produced a singular ACE metric with only one occupied source. Dense orthogonal plane waves make both target overlaps nonzero. This is a test-input correction, not a solver tolerance or algorithm change.
- Seven complete GS comparisons against the pinned original solver pass: PZ (3,-1), (4,-1), (4,0), (3,.05); HSE (2,-1), (4,-1), (4,.05). Electronic output/restart SHA256 and SCF residual text match exactly.
- Eight complete TDDFT configurations pass: k=2/3, exact/learned ACE Taylor, full Taylor, learned PT-CN. Solo/interleaved/reverse/close-reopen states and traces match exactly. Native original-driver current and energy use rtol=1e-12, atol=1e-14. A/B have distinct fixed dt, OpenMP route and step capacity, but common electron/temperature configuration as required by the outer guard.
- Recreated origin-atom synthetic input failed PT-CN with identical residual 3.1297286689933190e-10 before and after the HSE change. Restoring the existing multik fixture atom coordinates (.3,.4,.45) passes all eight configurations. Failure logs are retained; tolerances were not relaxed. These synthetic orbitals are not a converged physical GS.
- Source and installed-wheel Python suites each run 435 tests, 1 skipped, no failures. Source archive to wheel build passes.

All new outputs are persistent under `/Users/otobetoshihito/.codex/outputs/salmonpy-state-inputs-20261116`. This local path is evidence provenance, not a required installation path. Failed attempts are retained alongside final results. Prior temporary artifacts were lost; no claims rely on their continued availability.

## Evidence hashes

- `hse-input-build.json`: `d6d46e809cf5eabf32d81c92b4dcdaa75e4d4b48a451159f4f9f7fb03137cebc`
- `reference-build.json`: `b4bec1b51bec5d0a26a897e6b17db9a1397dc718029eb91227bb91ccc6994796`
- `hse-input-green-final/result.json`: `ac8f5facf35b0d9e699d9f3199be75fa0702a916db73b0e674d024f59713d9ea`
- `hse-input-green-omp2-final/result.json`: `750d8a00c8d663828a30c68f78fd9f8590c6bbd01117a3cf3487a3c84bec8051`
- `hse-input-gs-final/result.json`: `5ebf57126462eae34f7178fd1b9cd02d24659261e60106741d7c2e2e504aa810`
- `hse-input-tddft-standard/result.json`: `0d738db14e1e8ab95ec4938f4aa7517264fefc093a9e65e9b7007b21430a1489`
- `hse-installed-tests.log`: `ce8929b6b7be74e377acb97441e78263c494f50fa460dfe2d03db274dc6d0165`
- `hse-source-tests.log`: `8579fd109473afe4dc8c0f91cda56d45401532beb91e763f2bf6fba36e47fa57`

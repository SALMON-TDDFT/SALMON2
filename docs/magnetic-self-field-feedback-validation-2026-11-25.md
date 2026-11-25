# Orbital magnetic own-vacuum feedback — 2026-11-25

Implemented spatial E/H feedback in the existing total/own-vacuum validation driver using `--orbital-magnetic`. Each field group receives the identical sum of electric polarization current and interval orbital magnetization current exactly once. Fortran subtracts total/own E/H before molecular propagation. The policy removes the vacuum self-field while retaining environment-scattered return, and does not include vacuum radiation reaction or spin response.

## Field time alignment and persistence

Added the read-only `component_fields(indices)` query and `component_field_geometry()`. The owning field group performs one reduction of `[E,H]`; component staggering is preserved. H at field step n is the mean of the stored preceding raw half step and the current raw half step. A dedicated owned `coupling_h_left` array is used because observation arrays may already be time/spatially centered. It is captured before each field update and included in transaction copies and direct checkpoints. Current checkpoint signatures use revision2, without a legacy-format fallback. Added storage is three binary64 values per local cell including halo.

## Measured checks

- The previous validated changes were committed/pushed as `f9a495ce` to the independent SALMONpy branch before this implementation.
- Git-only source preparation/build and source/object eigensolver audits passed. Source inventory SHA256: `b0bc9897b12f79027e4588b04544b171f53da8dfad6118fc940949263d5eae90`.
- Field library SHA256: `5366b7005297fa9b5ed9fc1141f313a423b9a7ed9c2812f8fa28e3e007ccb600`. Reused validated molecular worker SHA256: `14557305ee88f46a99d3e15352f11ffb88c7c1c80626b93b21b2659a1b035897`; its matching grid-stream/density/orbital evidence was exercised without a new GS calculation.
- Spinless PZ C2H2 at off-grid position(-2.54,1.62,1.94) bohr; patch24^3, spacing2bohr; dt0.003au. Six outer MPI ranks (two ranks per FDTD group plus bookkeeping), one separate molecular worker.192 molecular updates end at field193. Nonzero interval magnetic current peak1.1991536694667786e-4au; selected predicted/accepted E/H maximum disagreement2.7755575615628914e-17au.
- Environment-scattered self-field relative to illumination-only: maximum absolute E difference3.8398232222666695e-7au and H difference4.1657759615627843e-7au. This is not a convergence or energy-conservation claim.
- Joint direct save molecular96/field97; new MPI launch imports both fields/LD plus quantum state and continues to192. Zero restoration evolution; all field histories and molecular observables match continuous execution bit-exactly. All six ordinary field/molecule output files match byte-for-byte.
- One field rank per group, explicit rejection/retry at molecular64: same histories, observables and six ordinary output files as the two-field-rank continuous control.
- Existing electric-only eight-step feedback retains the previous arrays bit-exactly.
- Separate direct E/H query tests cover native numerical overflow rejection, rollback/retry and process-exit field restart while retaining raw H brackets. E readout agrees with the existing raw E query.
- Source targeted106tests passed; sdist->wheel and fresh installed106tests passed. Tested wheel SHA256: `2665c933e84f1eab3979103cce865cac04ba1767ba2e6a4a39d9915e6a9ed0b4` (CPython3.14, Apple Silicon/macOS26).
- Installed-wheel six-rank continuation from molecular96/field97, with rejection/retry at128, matches continuous histories/observables and six ordinary files bit-exactly, again with zero restoration evolution.

## Limits and records

Fixed validation geometry/source, PZ and a single molecular worker. This is the connected E/H own-vacuum calculation, not an arbitrary-layout or universal failure-atomic production coordinator. Unexpected errors terminate that test MPI job. Spin magnetization, vacuum radiation reaction, physical convergence, full composed charge-continuity/energy validation and large-scale performance remain separate work. No extra heavy CI was manually dispatched and no physical calculations were added to required CI.

Authorized existing C2H2 GS reused. Canonical campaign: `calc2026/salmonpy-c2h2-nearfield-20261125`; machine-readable record `magnetic-self-field-feedback-validation.json`. Cases are named `magnetic-feedback-ld-long`, `magnetic-feedback-ld-restored`, `magnetic-feedback-ld-retry`, `magnetic-feedback-wheel-restored`, `magnetic-feedback-ld-external`, `magnetic-feedback-electric-regression`, `magnetic-fields-rollback`, `magnetic-fields-restored`. Local OneDrive storage does not establish cloud synchronization.

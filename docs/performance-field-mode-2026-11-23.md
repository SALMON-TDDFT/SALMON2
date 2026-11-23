# Field performance mode — 2026-11-23

User-directed optional unchecked stepping is implemented for NativeEHField and its embedded LD update. Select execution_mode="performance" at construction on every field-group rank. Default checked is unchanged; the property has no setter. No process-global mode flag is added. This is not yet a mode for standalone LD/SBE/TDDFT/HSE or for the outer coupling driver.

## Contract

Setup validates mode, inputs, field layout, ABI and recorded library SHA collectively once. The performance step packs current buffers and checks their local shape/lifetime, then enters the cached process-local I/O directory scope without a validation collective. The Fortran sp_eh_step_performance endpoint deposits currents with the same transpose map, executes the same E→LD→H update, reduces sampled E and integrates A with a binary64 volatile intermediate to preserve the existing separate multiply/subtract result. Compact owned-current rows are supported, including empty rows on ranks that own no mapped current.

It does not create the normal full-state before copy, scan all E/H/LD arrays for finite values, compare rank copies of A/J, call validation-only agreements, sample H solely for validation or perform automatic rollback. Mandatory build-time eigensolver audits and load-time SHA checks remain. Local ABI bounds, physical halo/sample communication, solver convergence/breakdown logic and requested ordinary output remain. The caller owns valid values/units/current ownership and matching collective order. Invalid values may propagate; a solver failure may terminate the job. There is no automatic state-preservation/retry guarantee.

Explicit protected-step targets, energy observations and checkpoint operations retain their own checks/copies. They are not silently disabled. Checkpoint format/physics are unchanged. The current mode is recorded in the field manifest; checked-mode step snapshots remain enabled by default. Other matter backends and outer response/coupling agreements remain checked.

C ABI: sp_eh_step_performance(npoints,nlocal,owned,J,A_old,A_next,E). Counts and owned0/1 are value int32. J is C(nlocal,3), old/new A and E are C(npoints,3), all float64; owned rows follow the existing ascending current-mask order. Caller-owned output buffers are consumed only on success. A nonzero status does not promise unchanged numerical state.

## Verification

- A unit fixture forbids validation collectives and supplies NaN: performance dispatch does not scan/reject it. Local buffer shape protection still rejects the wrong shape. This deliberately demonstrates unchecked behavior, not acceptable scientific input. Existing checked failure/retry tests remain enabled.
- Actual64-step reciprocal-tensor fields with a manufactured nonzero current after the initial step: checked/performance A/E history and complete per-rank checkpoint bytes match exactly.
- Spatial-owned current routing and step32 direct checkpoint/new-job resume match continuous execution exactly. Field MPI1/2/4 mapped A/E agree at rtol2e-10, atol2e-12.
- Conductive+LD and reference-frequency subpixel cases: checked/performance histories and final state bytes match exactly.
- Explicit energy observation with a caller-created backup works in performance mode and matches checked mode exactly; this includes the explicit observer/backup overhead.
- All730 solver inventory assets are tracked and byte-exact. Git-only export preparation/build and source/application-object prohibition audits pass. Final Git-built runtime and installed wheel match source controls, including direct restart.
- Source/wheel697 tests each, one optional skip each. Initial four failures were test-only __new__ fixtures missing the new mode field; fixtures now explicitly select checked, with no production fallback. Existing GS-free C/Fortran protocol queries and NLCC XC/core-force MPI1/2 checks pass. No fresh GS/SCF or remote CI dispatch.

## Short local timing control

72³ grid, four field ranks plus one bookkeeping rank,64 updates, one thread/rank. No explicit energy observation/outer backup; ordinary point-output settings are identical. Two repeats: checked loop6.680/6.740s, performance5.229/5.530s. Medians6.710 versus5.379s. Summed individual process peak RSS1567.8 versus1050.8MiB. All histories and final checkpoint bytes agree.

This is a short single-Mac control, not a universal speedup, multi-node scaling claim or simultaneous aggregate-memory bound. The timing library and final Git-built library have identical source hashes; both binary records are retained.

Canonical evidence: calc2026/salmonpy-performance-field-20261123/summary.json, timings.json and verified-artifacts/. Final wheel and both runtime binaries/records are retained. Local OneDrive storage is confirmed; cloud synchronization was not checked.

## Remaining rollout

The next stage is the same explicit per-instance choice for standalone LD/SBE, then TDDFT/fast-PR HSE and outer mixed coupling. Do not infer those checks are disabled by selecting the field mode. Keep checked as the default and require explicit performance opt-in.

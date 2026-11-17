# Protected restore in the ordinary mixed sample — 2026-11-17

## Change

The ordinary `mpi_mixed_native --direct-resume` now defaults to protected restore.
The same internal coordinator accepts a paused advanced session. It authenticates
the unchanged native-direct/2 inventory/input/topology, collectively agrees the
restore/protected policy, and validates/preloads history before any LOAD.

Native field candidates include PML and every field-local LD P/J/history array.
Standard TDDFT and fast-PR HSE use the previously verified numerical-candidate
ABI. Endpoint SBE snapshots load into independent handles and wrapper/adapter
candidates. Coupling arrays, clocks and owner lists are prepared before the common
MPI ready vote. Native candidates publish, then prepared Python bindings replace
the old bindings. Old snapshot handles retire only after success; release errors
are separate `checkpoint_cleanup_errors` diagnostics.

Handled preparation failures do not invalidate the session or replace diagnostics,
clocks, object references, field/LD/electronic data or owner lists. No propagation
or SCF occurs during restore. Older native ABIs, isolated HSE workers and custom
direct-only targets are rejected without fallback. Explicit legacy behavior is
`--unprotected-direct-restore` / `protected=False` and retains its old restrictions.

## Evidence

`tools/test_mixed_protected_restore_mpi.py` runs the ordinary sample's actual
configuration, session, stepper and checkpoint coordinator. After saving, it
advances two steps and injects one rank's exception after a valid later electronic
candidate has loaded. Every native shard saved before/after that attempt is
byte-identical, including field-local LD data. SBE snapshot bytes/integrity,
Python clocks, currents, adapter references and model-owner membership also remain
unchanged. Originals then advance once more. A successful retry restores the
saved step and reproduces the following three complete frames exactly, without
replaying time steps inside restore.

Four-rank configurations verified:

| Field | SBE | Electronic | Saved / failed step | Environment |
|---|---|---|---|---|
| 1 | 1 | Libxc-PBE TDDFT, orbital2 | 64 / 66 | Source and installed wheel |
| 2 | 1 | Libxc-PBE TDDFT, orbital1 | 64 / 66 | Source |
| 1 | 1 | Fast-PR HSE Taylor, k2 | 128 / 130 | Source |
| 1 | 1 | Fast-PR HSE PT-CN, k2 | 128 / 130 | Installed wheel |

Every configuration contains two field-local LD materials. MPI tests also reject
different protected policies before staging. Installed-wheel tests reject a
bad-but-rehashed history file before staging, proving the history check is not
merely a payload checksum failure.

A separate fresh four-rank job invokes the ordinary installed module CLI with
`--direct-resume` and no legacy flag. It resumes the field2 checkpoint at64,
advances only to68, and saves68. All A/E/j, previous-current and clock values at
65–67 match the independent continuous frames from the live-session test exactly.
The restored frame prefix is retained, giving68 frames in the new archive.

Source and installed-wheel suites each pass502 tests, one optional skip. New unit
tests reproduce collective policy rejection and a later invalid SBE candidate,
checking original ownership and exact future currents after retry. Existing
three-rank preflight tests still pass eight collective rejections with zero LOAD
calls. The source-distribution-to-wheel build/install and paired/AI guide sync pass.
Independent read-only review found a policy-branch mismatch. Its unit regression
failed before the fix and passes afterward; real MPI rejection and history
admission ordering were also checked.

Persistent artifacts, logs, continuous reference frames and wheel:
`/Users/otobetoshihito/.codex/outputs/salmonpy-mixed-protected-20261117`.
Existing audited native-candidates-20261117 v4 libraries were reused unchanged.
No native kernel, checkpoint encoding or numerical update formula changed. Builds
use2 workers; jobs run sequentially with at most4 ranks and OpenMP/BLAS1. No remote
Actions run. The campaign does not establish a hard RSS cap or converged optical
properties; these are connection/restart regressions on the local Mac.

## Limits

The host must pause evolution and keep healthy MPI participants. Native
publication still includes allocating assignments: OOM/process failure, Fortran
STOP, failed communicators or violating the prepared-publication contract cannot
be rolled back by this coordinator. An unexpected publication error invalidates
the session and does not claim original preservation. This is one native lifetime
per process, exact binary/input/topology only. Stable public native ABI, arbitrary
backend cloning and topology redistribution remain separate work. After success
retrieve replaced snapshot states through `session.states`; old aliases are retired.
Live rewind does not truncate already-written native output streams. The tests
compare scientific frame arrays, not append-only output files. Fresh CLI resumes
use a new output directory and retain the saved frame prefix.

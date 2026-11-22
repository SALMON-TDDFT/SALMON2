# Isolated worker live savepoint validation — 2026-11-22

## Implemented scope

PZ, isolated fixed-ion middlepoint worker, one MPI rank and one pending in-process
savepoint. Numerical ownership is local to main_tddft. Both orbital buffers,
temporary orbital, density/potential histories, Poisson state, full active RT
histories, energy, XC payload/coefficient, field vectors and accepted index copy
in Fortran. No new module-global selector or cache was added.

CAPTURE_E flushes the audited streams; RESTORE_E closes them and copies numerical
state; the validation driver truncates its captured output journal; RESUME_E
reconnects streams. RELEASE_E retires the backup. Partial stream-close failure
leaves writes suspended until explicit reconnection. These are development
commands, not public Python methods or durable checkpoints. Unsupported output
modes reject capture only. This is not a blanket allocation/OS/MPI failure claim.

A monotonic counter increments only after a real time_evolution_step call and is
never restored. PROBE_E exposes it without propagation. This distinguishes direct
copy restoration from replay even when resulting observations coincide.

## Evidence

Existing authorized C2H2 GS reused; no new SCF. dt=.003 au, capacity20, fixed ions,
three explicitly supplied E vectors cycling through nonzero x/y/z and sign changes.
Savepoints at steps0,7,8,20 restore in the same PID. Native update count changes
by zero during restore/resume. Continuing twice from each savepoint reproduces
all nine E3 observation values exactly, including norm and energy freshness.
The two principal output streams match continuous execution byte for byte.
A second live worker is unchanged; duplicate capture and absent restore return14
without altering the state. The validation driver rejects an old savepoint after
release/reuse before sending any restore command.

Ordinary E3 response/reference maximum dipole difference remains
3.85523175341039e-12 e*bohr and energy difference2.1316282072803006e-14 Hartree.
Independent workers and repeated runs are unchanged. The public loader, replay
snapshots, disk replay and killed-worker recovery also pass. Existing molecular
LD/FDTD MPI6 reject/retry with the new worker exactly matches all previous saved
arrays and every six-face energy record. That public coupling still uses replay
for electron cancellation; it does not yet call these direct commands.

The maintained solver was exported from Git-tracked files only, prepared and
built cleanly, retaining source/application-object prohibition audits. NLCC XC
and force local checks on MPI1/2 pass before advancing the source capability pin.
Source and installed-wheel suites each run689 tests, one optional skip, no failures.
Obsolete literal argument-list checks were removed; byte-exact source inventories
and real response/restore comparisons remain. Installed manuals, source identity
and worker source bytes match the working tree. No physical CI job was dispatched.

Executable SHA256:08b045275781442d990fc212ea91c21268551f3bde4dfe62f1fe9046b4f356d2.
Source inventory SHA256:a83d21f5ea6192ea850296f0b59374540ea4e87fad68dbd1f451c68f4d412a61.
Canonical evidence is under calc2026/salmonpy-c2h2-feedback-20261122:
native-savepoint-final-summary.json, native-savepoint-final-artifacts/,
native-savepoint-final0/7/8/20, native-savepoint-final-guard7 and
native-savepoint-final-feedback-retry. Earlier failed tests and binaries/records
were preserved separately. Local OneDrive storage is not confirmed cloud sync.

## Remaining

Public Python snapshot and FDTD cancellation integration, direct durable disk
persistence, and broader XC/internal-MPI/physical-convergence validation remain
separate. Current public IsolatedTDDFTState snapshot/restore and disk checkpoints
still replay accepted history. No claim of process-independent native snapshots.

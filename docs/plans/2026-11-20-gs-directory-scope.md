# Group5: retained GS calls without host chdir

Use an explicit native file-resolution scope for each retained-GS operation.
Keep the owner directory in existing Python state/session metadata. An additive
native enter/leave protocol sets a temporary resolver prefix; relative file
operations resolve under it, absolute paths keep their meaning. Never change the
process CWD. Pop the scope in finally, including rank-local admission failures.
Sequential native calls remain required; this is not threaded solver reentrancy.
Old binaries keep their current scoped-chdir implementation.

Route Fortran OPEN/INQUIRE and central filesystem/POSIX/MPI-file operations for
GS initialization/continuation/restart dependencies through the same resolver.
Keep logical namelist paths and numerical serialization unchanged; do not embed
host paths in physical inputs merely to work around CWD. Scope-free standalone
calls retain legacy behavior. No changes to immutable prepared source snapshots.

Verify a failing Python regression for host CWD changes; native synthetic scopes,
absolute/relative paths, missing/nested scopes and cleanup; actual ordinary and
full/reduced HSE owners, protected restore and previous-runtime numerical parity.
Use a host CWD sentinel/observer and check files remain in their owner folders.
No GS data collection for this campaign. Source/wheel checks, JP/EN/AI/list sync;
no routine CI. Other TDDFT/FDTD Python chdir wrappers remain separate work.

## Accepted checks

The pre-change Python path fails the explicit no-chdir regression. New source
and wheel each637 tests:636pass/1skip. Native MPI2 scopes preserve host raw file
lookup and relative/absolute paths, reject missing/nested scopes, and separate
siblings. Real MPI2 one-rank admission failure skips the callback everywhere,
cleans entered peers, and permits retry; callback exceptions also clean scopes.
Ordinary MPI2 and full/reduced-HSE MPI4 (three groups): no Python chdir on the new
runtime, numerical arrays bit exact to audited previous libraries, protected late
failure/retry exact and restore SCF work0. Installed ordinary/reduced-HSE controls
and public ordinary/HSE protected samples pass. Numerical stream writer bytes
and fields unchanged. Source/object audits and source snapshot/wheel/JP/EN/AI
parity pass. No GS datasets collected; scalar/synthetic evidence is stored at
calc2026/salmonpy-gs-directory-record-20261120. No CI or confirmed cloud sync.

# Direct molecular/FDTD bundle restart

Save only after acceptance: total FDTD step n, own-vacuum FDTD step n,
molecule step n-1. Check live Fortran field indices and the worker's real
accepted index. A read-only field ABI reports both live and pending indices;
normal timesteps gain no new communications.

The existing field binary stores E/H, LD, PML and current state per field rank.
The isolated direct disk API stores the molecule. A complete manifest binds
per-rank configuration/MPI layout, hashes, sampled A/E, selected prior E and
optional caller-owned JSON frame history. Publish only a completed new directory.
No automatic SCF, builds, elapsed-time replay or layout conversion.

Restore validates metadata, distributes hash verification, then stages both field
groups through NativeRestoreCandidate and the molecule in an independent worker.
Verify both staged Fortran steps n and molecular step n-1 with updates=0. Stage
field output streams with originals retained; complete I/O before publishing.
Collectively discard on preparation failure, preserving existing numerical states
and step numbers. OS/MPI failure or native allocation failure during publication
is outside this bounded guarantee and must not be reported as preserved.

Rebuild/audit a Git-only source export; reuse the already validated isolated
executable and GS. Compare continuous, stop/save/exit and new-MPI-job resume at
field step30; all returned frames and energy records must agree exactly. Inject
candidate preparation failure on a field leader and on the molecule owner, then
retry. Test metadata index mismatch before publication. Validate source/wheel
suites, installed example and unchanged baseline samples. Synchronize paired
manuals, AI JSON and current checklist. No extra physical CI.

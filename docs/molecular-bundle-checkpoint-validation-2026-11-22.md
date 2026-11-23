# Direct molecular/FDTD bundle restart — 2026-11-22

## Result

MolecularFeedback now saves/restores the total FDTD group, its own-vacuum reference
FDTD group and the isolated molecule together. Accepted indices are explicit:
field n / reference n / molecule n-1. Fortran sp_eh_step_indices reports both live
and pending field indices. The worker reports its native step and actual evolution
count. No normal-update collective or numerical integrator was changed.

A complete new-directory manifest binds field-rank binaries (E/H, LD, PML,
current), molecular direct state, configuration/runtime hashes, sampled A/E,
saved selected E, point-output streams and optional caller JSON history. Field
payload hashes are distributed over ranks, not redundantly computed by every
rank. A shared filesystem and identical input/capacity/runtime/MPI layout are
required. No old-format fallback or implicit replay is added.

Restore prepares detached native field candidates and an independent molecular
worker; validates all saved/native step indices before publication; stages output
streams with retained originals; then publishes after collective readiness.
Candidate/preparation failure discards staged payloads. Native fatal allocation
failure during publication, MPI loss and arbitrary OS failure are outside this
bounded guarantee. Candidate memory and output staging space are additional.

## Actual MPI evidence

Existing authorized C2H2 GS reused; no SCF. PZ, fixed ions, one molecular worker
rank/thread; MPI6 with two ranks in each FDTD group and bookkeeping ranks.
Scalar LD environment; dt=.003 au, capacity64. Save at field30/reference30/
molecule29, exit the entire original MPI job, then initialize a new job and
restore directly. All owning field ranks report30, the molecule reports29 and
updates=0. Only34 remaining electronic updates reach field64/reference64/
molecule63. The sample does not call start on resume.

Every saved frame (64) and six-face energy record matches uninterrupted
execution exactly. Normal execution also matches the previous field library.
All four field-rank final checkpoint binaries match exactly, covering E/H, LD,
PML, current and serialized field state. Both groups' point-output files match.
Installed-wheel save and a separate installed-wheel MPI job resume reproduce the
same results. Restored indices are included in the sample output.

Fault tests first advance the target to field5/molecule4, outside restoration.
Then they inject reference-stream preparation failure, molecular binding failure
and a complete-manifest/native-candidate index mismatch. For every failed case,
all original field numerical payload bytes, original molecular numerical payload
bytes, cached fields, PIDs and indices are preserved; candidates are discarded.
Retry succeeds and all64 returned frames/energy records equal the control.
Public field evolution and molecular step calls are forbidden during restore by
the test harness, so coincident values cannot hide an elapsed propagation replay.
Source and final installed-wheel variants both pass. The reusable development
entry is tools/test_molecular_checkpoint_mpi.py with the sample arguments/--resume.

## Molecular byte comparison: active values versus unused storage

Independent-job molecular checkpoint files are not wholly byte-identical.
A typed stream inspection consumed all12552252 bytes and attributed every
byte difference to tpsi.update_zwf_overlap, rt.Qp_e and poisson.wkbound.
Wavefunction arrays, density/potential, active RT histories, energy, fields and
other serialized values match exactly. The current fixed-ion worker never uses
the temporary orbital's overlap flag; restore invalidates orbital halos explicitly.
Qp_e is initialized at index0 but otherwise unused by this worker. Poisson writes/
reads only its direction-dependent boundary prefix; all differences are beyond
the maximum used6144-element prefix (24x24x32 grid, ndh4). That active prefix also
matches. No numerical buffers were rewritten merely to force an archive hash.
This scope-specific conclusion is not a generic permission to ignore arrays.
The complete field checkpoint binaries, unlike these unused molecular areas,
are byte-identical. See molecular-bundle-active-payload.log and its typed-IO
inspection script in the canonical artifacts.

## Build, tests and documentation

Tracked-only export preparation/build passed. Source/application-object
prohibition audits passed before newly linked artifacts were used. The only
Fortran addition is the read-only index query. Local prescribed-density NLCC XC
and core-force tests passed at MPI1/2; GS-free C/Fortran ABI queries passed.
The isolated worker is reused unchanged with its recorded SHA256, not rebuilt.

Source and installed wheel each ran692 tests with one optional skip. Final code,
paired guides, AI JSON, sample, header and example-index bytes match the installed
wheel. Documentation uses explicit accepted step numbers, not an undefined clock
or time-position shorthand. No extra physical CI was manually dispatched.

Canonical local OneDrive campaign:
calc2026/salmonpy-c2h2-feedback-20261122/molecular-bundle-summary.json and
molecular-bundle-artifacts. Important runs: bundle-source-save,
bundle-source-resume, bundle-source-live-fault-retry, bundle-wheel-save,
bundle-wheel-resume, bundle-wheel-live-fault-retry. Final-state controls:
bundle-continuous-final-state, bundle-resumed-final-state, bundle-wheel-final-state.
Local storage does not establish completed cloud synchronization.

Maintained source inventory SHA256:
20ab9af5c08d22046938d80e92e8677b54548b4ec6b7f02bb156314f1c6e5fa0

Field library SHA256:
7f3a80efe710820143306259ce78d96887a820f53f41e186c88f8c8807bf8442

Reused molecular executable SHA256:
473f3f8f91cc44e61995e55457e484b6fb8c1120f695e39ce28c83cd3bb750e8

Wheel SHA256:
358e82163ba5a4456ed50923adaa37b712f4e08b40208519d4689aa059b48989

Tensor/subpixel/arbitrary-background energy observation, physical convergence,
other XC/inner-MPI molecular profiles, cross-layout/platform restart and facility
scale remain separate work. The short C2H2 campaign is wiring/regression evidence.

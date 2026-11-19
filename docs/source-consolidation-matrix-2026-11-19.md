# Effective solver migration matrix

The raw fast-HSE reference is not the current fast-HSE engine. Resolve RT-state
and explicit-routing complete-file replacements before comparing sources.
`python tools/runtime_sources.py` prints verified path/hash/origin inventories
without creating work directories, executing code or changing a source tree.
Its results exactly match separately prepared standard670/fast-HSE725 files.
There are573 identical,77 differing,20 standard-only and75 fast-HSE-only paths.
These are byte comparisons, not a certificate of behavior/physics equivalence.

| Area | Current migration requirement |
|---|---|
| Retained GS / in-memory GS→TDDFT | Standard-only gs_model/storage/seed/embedded_gs and their initialization/input/communication contracts must move; raw fast-HSE cannot replace these APIs |
| Standard embedded TDDFT / failure handling | Retain embedded_tddft, worker and model-input interfaces; reconcile fast-HSE state/solver structures rather than overwrite either version |
| FDTD / LD / SBE | Retain embedded E/H, kernel, shared LD/SBE files; enforce native/solver shared-kernel equality |
| Fast PR HSE / symmetry / PT-CN | Keep current explicit-routing sources and host native_adapters; these are not replaced by legacy standard HSE |
| NLCC / XC / forces | Compare effective force/salmon_xc/semilocal paths and existing validation before porting; file-name or symbol presence alone does not establish that a fix is preserved |
| Checkpoint / MPI / model inputs | Reconcile differing implementations and model/partition identities; keep independent candidate restore and numerical rollback evidence |
| Build / source provenance | One maintained current tree with no-HSE and fast-HSE profiles; preserve origin commits/licenses/inventories and separately pinned reference evidence |

The first Python cleanup merges only GS archive internals. It does not consolidate
these native implementations. Three private modules become one, while public
gs_archive/gs_native/native_gs entry points and _gs_tddft lifecycle remain.
Original function ASTs are compared after only private export renaming and removal
of a now-unneeded self import. No compatibility forwarding modules are retained.

Current source trees are preserved pending native migration and parity checks.
Do not build or execute historical snapshots just because they have a recorded
name/hash. Keep the mandatory source/object prohibition checks on every new build.

The user-approved convergence workflow is separate: low-k PP-radius/3 to about
0.1 angstrom grid scan, then impulse dt scan from about0.2au; numerical failure
stops that trial only and nonbreaking dt differs from value-converged dt.
No real GS/SCF or convergence campaign is started in this cleanup batch.

## SBE sharing stage (current)

Five SBE paths in fast explicit routing now use existing standard source files:
CMakeLists, sbe_types, sbe_velocity_kernel, gs_info_ssbe and bloch_solver_ssbe.
Prepared SBE kernels are checked against native/ before approval. No third source
copy is introduced. The previous raw fast reference remains unchanged.

Current inventory:standard670/fast727;578 identical,74 differing,18 standard-only
and75 HSE-only paths. Both HSE-on/off builds pass with the same prepared source.
The HSE-off build exposed unguarded EXX recovery/history references in Taylor and
time_evolution_step; only these dedicated calls/imports are now guarded. Existing
wave/density/current finite checks remain. HSE preprocessing has identical
computational statements (NVTX line labels excluded). This is not complete native
source unification or standard embedded GS/TDDFT/EH migration.

## FDTD/LD/EH sharing stage

Fast routing also references six existing standard Maxwell files. The new host
library exposes the exact embedded E/H implementation; its record declares
field_abi protocol1 after a compiled-protocol probe. Legacy field records remain
accepted, and a separate executable-build.json serves batch/MPI FDTD. Host MPI
process separation and existing public field methods are unchanged.

Counts:standard670/fast730,584 identical,71 differing,15 standard-only,75 fast-only.
Five standard-only alternative HSE files are not porting targets under the user's
fast-PR-only policy. The bridge is already a fast host build input; output/guard
functions require organization, not blind duplicate copying. Primary unfinished
functionality is retained GS/model-input/ordinary TDDFT integration, including the
necessary differences in shared structures/input/initialization.

3D scattering can produce A components outside a reduced model's little group.
The ordinary mixed test therefore uses full-k HSE models, preserving all vector
components rather than projecting them. This is an implementation regression,
not convergence validation or evidence that arbitrary reduced meshes are valid
for a chosen structure/drive.

## IO initialization prerequisites

Shared IO CMake/model_inputs/step_output now build in the fast tree. The duplicate
native adapter output file is removed; host/serial builders reuse the solver
object and use the common source when linking older builds lacking it. Effective
inventory:587 identical,70 different,13 standard-only,75 fast-only (fast732).
The standard value-only input snapshot does not cover fast HSE's additional
settings: its existing full HSE snapshot remains authoritative. Parser workspace,
GS owners, ordinary TDDFT registration and initialization differences remain open.

Serial compilation exposed an unconditional MPI-only comm_is_attached reference
in symmetry initialization. The attached-host test is now conditional on USE_MPI;
standalone serial code retains the ordinary field-validation branch. This does
not weaken the arbitrary-polarization or little-group checks in the MPI path.

## Fast parser workspace migration

Fast inputoutput now adopts model-owned scratch/diagnostic/default/common-value
storage and broadcasts its metadata in the isolated path. HSE model registration
uses an isolated parser; complete HSE settings remain captured by the existing
HSE input owner. Fast namelists/defaults/conversions/validation are preserved.
Fresh legacy input logs are compared against the unchanged prior fast library,
not merely the refactored library's own legacy path. Ordinary parser tests still
pass. No additional solver path/Python module is introduced by this stage.

Parser workspace separation is implemented; common parser values are not a full
HSE model record. Retained GS/storage/ordinary TDDFT lifecycle migration and the
remaining source consolidation are still open. Effective path counts are unchanged
from the IO stage:587 identical,70 differing,13 standard-only,75 fast-only.

## Fast retained GS import / common numerical seed

Fast GS retains imported restart states without SCF and exports detached common
seed arrays/ownership metadata. The standard gs_seed source is referenced in the
fast inventory. Per-context child/cache/XC/HSE resources are closed independently.
Normal SCF iteration/output statements remain byte-exact; no new SCF run is used
for this migration check. Python GS handle/ordinary TDDFT initialization migration
is still open, and a common seed is not a full fast-HSE model or checkpoint.

Current inventory:588 identical,70 differing,12 standard-only,75 fast-only paths;
fast733/standard670. Do not remove the standard GS/TDDFT source tree yet.

## Python fast GS handle connection

NativeGSBackend.load now accepts the fast gs_abi/1 record, retaining the standard
Python client API. Fast C handles own stable context allocations and generation
identity; reads return detached density/occupation/orbital arrays. Full/reduced
model symmetry scratch is isolated. Metadata for uncomputed energy/residual/
tolerance/convergence is None; computed-only methods/arrays remain unavailable.
Standard SCF/export/GS→RT behavior is preserved. Existing sample/test files gain
an import-only mode; no additional Python module or sample file.

Fast GS SCF creation and GS→TDDFT initialization remain to port; native standalone
GS support does not yet make those Python fast-profile methods available. Current
source counts remain588/70/12/75. Both source profiles stay until migration passes.

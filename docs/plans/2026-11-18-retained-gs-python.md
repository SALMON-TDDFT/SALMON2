# Retained standard GS handles and Python connection

Goal: independently owned in-memory standard GS results for different DFT inputs,
readable from Python on a supplied MPI group, without changing SCF numerical order.

Architecture: move existing main_dft computational locals into stable s_gs_storage
owners; legacy wrappers still create temporary storage. Retaining GS runs keep
wavefunctions/density/occupation/eigen/force/XC and child MPI resources until close.
Generation-tagged C handles expose metadata and explicitly requested read-only
array copies, with local shapes/lower bounds. A Python NativeGSBackend session
owns one host MPI attachment; state close is independent, session close closes all.
Python DFT.run stays separate-process; retained GS→TDDFT handoff, GS checkpoint,
HSE/SOC/+U/DC and concurrent GS/TDDFT/threads remain separate work.

1. Red Python public entry/import tests and native protocol presence check.
2. Refactor standard GS storage/solver; preserve default arithmetic and exact SCF
   convergence flag. Extend internal model run with optional retained storage.
3. C protocol create/info/shape/read/close/resources with generation validation,
   collective rejection and explicit ownership cleanup. No automatic GS disk output.
4. Python backend/session/state/read-only array interface, immutable DFT input and
   material IDs, collective metadata admission and nonconvergence warnings.
5. Real MPI1/orbital2/spatial2: A/B retained memory unchanged across sibling solve,
   close/reuse/stale tokens; compare energies, density/orbitals/eigen/occupation/
   forces against independent process runs. Reuse audited controls, no GS bundle
   collection. Source/object audits for all new executables, no Actions.
6. Source/wheel verification; paired guides/sample/AI JSON/checklist and evidence.

## Measured implementation

Stable storage refactors existing computational locals into owned components; no
SCF arithmetic/iteration reordering. Native convergence is copied immediately from
the existing native report, without reinterpreting the final-iteration criterion.
The C owner keeps its Libxc and child MPI resources until close; request caches
retire before array lifetime ends. State metadata/arrays never consult another
model's scalar view. Python material ID/calculation and backend config are fixed.

Ten retained A/B vs independent MPI/OMP cases pass across PZ/Libxc-PZ/PBE,
real/complex orbitals, cell/grid/atom/electron/k changes and repeated slot reuse.
All density/orbital/occupation/eigen/force values agree exactly; no raw arrays are
collected. Legacy nonretained Fortran GS9 cases match old executable diagnostics.
TDDFT MPI6 currents/direct resume match old controls.

Darwin diagnosis: OpenMPI's wrapper injected flat_namespace, allowing NumPy's
Accelerate zdotc ABI to replace the configured OpenBLAS symbol. A two-level link
also needs commons,use_dylibs and libmpi ahead of tentative ignore-tkr sentinel
copies; otherwise MPI_IN_PLACE is a different address and nonzero admission
reductions fail. The permanent shared host/library mpi_link_command reproduces
the successful link exactly. Rejected binaries/logs are preserved as evidence.
Linux wrapper/flags are retained. No Actions launched in this chat; unrelated
LD-wheels email was verified as QEDFT-branch push0fd3b065, not d1b7b79d/05926ba0.

Final source and wheel:635 tests each, one skip. Installed MPI6 handles compare
all arrays to independent same-partition controls, reject rank-local invisible
work before native collectives, and preserve two different SCF caps/tolerances.
A separate lax-threshold diagnostic verifies native true convergence flags; the
normal four/six-iteration runs verify warning-only false flags. Installed MPI2
public sample passes. No physical-convergence certification. Selected evidence
at canonical calc2026/salmonpy-retained-gs-20261118, with no GS restart/wavefunction/
density collection. Rejected libraries remain explicitly NOT APPROVED. Local
OneDrive writes verified; cloud synchronization unverified.

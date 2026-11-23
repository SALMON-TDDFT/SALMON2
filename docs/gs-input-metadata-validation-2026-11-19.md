# Captured GS input descriptions — 2026-11-19

GSInputMetadata/describe_gs_input and NativeGSState.input_metadata provide a local,
immutable description of captured input text/assets. Explicit length units,
supported ordered species/coordinates, grid/k/XC/PP settings are described; missing
native defaults or unsupported/ambiguous data remain unknown. Source IDs identify
captured input/assets, not physical equivalence, binary authenticity or restore
compatibility. Imported provenance/convergence stay unknown and calculation None.
No uncomputed forces or imported-owner portable export are enabled.

Input text and the already obtained convergence flag are retained during existing
run/load. Description itself uses no source-file read, native call, MPI or solver.
It is valid after close; returned JSON copies are detached. Optional f90nml is
loaded on explicit description only via the metadata extra, not ordinary import/
run/load. NumPy remains the only required dependency. Standalone custom hosts may
pass explicit text and relative-name-to-bytes assets without a runtime.

Synthetic tests cover A_eV_fs/au units, nonorthogonal cells/external coordinates,
indexed species and channels, k tables, unknown units/defaults, mutation/closed
owners and missing optional dependency. A first MPI diagnostic run failed due to
an incorrect all(bool) assertion in the runner (not the solver); failed scratch
logs are preserved and the corrected runner passed in fresh directories.
Existing-GS/native controls: PZ2-k orbital2/spatial2, PBE8-k orbital2, TB-mBJ8-k
spatial2 and regular control, four total ranks. Describing inputs before/after
close leaves all16-step charge currents, independent-state and direct8+8 resume
results bit-exact with file controls, zero replay and zero SCF calls.

Source626/installed626 tests pass (one skip each). Installed native MPI4 control
and the public memory_gs --describe-input MPI2 example pass. Fresh installed import
loads neither f90nml nor mpi4py. Wheel built offline from sdist, bundled-kernel
source/all-object audits and native hashes pass. No electronic native rebuild,
new SCF, new reference GS or GS restart/orbital/density collection. Only existing
GS scratch working copies were used. Paired public/GS guides, optional sample,
API inventory and AI search fields/embedded guides synchronized.

This slice is declared configuration description, not canonical physical GS/native
model binding, a full SCF mixing checkpoint or changed-topology restoration.
It does not certify physical SCF/optical convergence. Measurement environment is
Apple Silicon/GNU/OpenMPI/Python3.14; other platforms remain unverified.

Selected logs/hash records/wheel (no GS arrays):
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-gs-input-metadata-20261119
Local writes verified; cloud synchronization unverified. No push/Actions this
slice; pending external-k fix and metadata work remain uncommitted afterd67ce847.

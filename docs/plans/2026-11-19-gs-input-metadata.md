# Declared GS input metadata

Continue approved GS input/ownership work without inventing legacy provenance.
Use a local optional description API over retained text/assets, rather than
making every native load depend on a new parser or trusting guessed defaults.
Reuse f90nml as an explicit extra and existing JSON freezing/atomic-unit constants.
No native ownership/propagation changes, implicit downloads or new SCF tests.

1. Detached immutable descriptor and JSON copies; source-byte identity distinct
   from user material/native IDs and physical/restore compatibility.
2. Explicit length/geometry normalization, declared grid/k/XC/potential references;
   unknown omitted/unsupported metadata noted instead of rejecting calculations.
3. Capture input and known convergence at existing state creation, expose local
   state.input_metadata before/after close without rereading original files.
4. Synthetic geometry/index/dependency/mutation tests; actual existing-GS MPI
   PZ/PBE/TB-mBJ current/independence/no-replay controls with description calls.
5. Installed wheel, public API, paired guides/AI JSON and optional sample.

This does not provide physical GS identity, imported DFT binding/portable export,
full SCF mixing snapshots or changed-topology native restoration.

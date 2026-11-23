# Heterogeneous TDDFT inputs in ordinary mixed sessions

Goal: make the approved per-material-ID design usable for different TDDFT inputs
in normal mixed calculations, while recording the remaining same-process work.
Architecture: reuse existing MaterialBinding/provider/backend, immutable input
bundles and separate MPI solver groups. Add an optional sample catalog mapping
material IDs to local bundle/input/decomposition/group-width. No new solver ABI,
no implicit downloads/builds, and no claim of different-input process reentrancy.
Defaults and existing archive signatures stay unchanged without the catalog.

1. Test catalog parsing, per-material settings, duplicate IDs, relative paths and
   invalid declarations without native initialization. Observe red before implementation.
2. Connect entries to normal native layer assembly, serial reference replay and
   input-hashed checkpoint identity; same-material sharing stays opt-in per ID.
3. Run PBE and TB-mBJ in separate internally parallel groups in one field/LD/SBE
   calculation; compare independent controls, shared-width changes and protected
   no-replay restart/failure-retry. Reuse existing GS fixtures/native binaries.
4. Update paired guides/AI/checklist with current support vs unfinished same-process
   ownership (inputoutput/read_input, initialization, common I/O/timers/communications,
   specialized resources); verify source/installed tests and sdist-wheel locally.

Decision: implement existing separate-group architecture first; a blanket global
variable capture/restore or arbitrary second read_input in a live model is unsafe.
Native model handles currently have no model ID and model_new rereads the original
GS under common input. Different-input native model handles remain later work.
No Actions; resource excess remains advisory; no GS data is collected.

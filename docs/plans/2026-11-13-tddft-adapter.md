# TDDFT A-to-J connection

Goal: preserve SALMON's existing middlepoint TDDFT update sequence and load a
reusable ground-state restart bundle into independent evolving contexts.

The full solver uses module globals and scratch storage. Reentrant shared-library
handles would require an audit/refactor across initialization, XC, communication,
PP and propagation. Instead the first adapter owns one persistent SALMON process
per state, retaining the existing initialization_rt/time_evolution_step calls.
A worker-only stdio protocol supplies endpoint A; ordinary execution is unchanged.
The executable is separately built/audited and hash-identified, not bundled in wheels.

Scope: serial MPI layout (OpenMP allowed), periodic fixed ions, transverse supplied
total Ac=A/c, middlepoint propagator, PZ (other functionals pending native validation), fresh GS rather than RT
checkpoints. Predictor-corrector retains its existing path. Reject unsupported
modes explicitly. The input bundle is snapshotted; each process receives a private
copy, including restart/pseudopotential files. dt/capacity are fixed per state.
J retains SALMON's electronic matter-current sign, matching the existing SBE.

Steps:
1. Tests for loader snapshots, protocol lifecycle, independent states and failure.
2. Fortran worker helpers + sentinel input + hooks around the unchanged timestep.
3. Python initial-state loader, owned process lifecycle and common batch cleanup.
4. Build/audit; compare ordinary before/after output, live step sequence and batch,
   interleave independent states using an authorized GS fixture.
5. Document supported modes, lifetime, errors, sign and reproducible validation.

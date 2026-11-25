# Molecular density and E-component current connections — 2026-11-25

Implemented explicit current-step electronic charge-density readout, including molecular coordinates and quadrature volumes, in Fortran. It does not propagate/reinitialize/replay a state. Signed charge is negative electron number density; fixed ionic charges are not included. A private worker binary transfer file is updated only on explicit requests; permanent density archives remain opt-in. An executable-bound optional validation receipt gates the public method.

Added checked rank-local E-component current injection, component-origin metadata and read-only raw E queries. Current lists may differ in length and only each rank's owned cells contribute. Macro-site averaging is not repeated. Existing E->LD->H evolution and checked recovery are reused; no application-wide MPI communication phase was added. Explicit raw E collection communicates within the field group.

## Executed checks

- C2H2 PZ fixed-ion existing GS,18,432 points. Density reads left all observations exact against an independent no-read control. Electronic step2 disk save, worker termination and a new worker's direct restart to4 matched density, coordinates, volumes and observations exactly. Restore added0evolutions.
- Actual charge changes were mapped in Fortran and injected into one/two field ranks; nonzero fields were observed. Two-rank lists contained900/520owned current entries. Raw injection matched the existing macro transpose exactly in a separate control.
- NaN current rejection retained state. Finite duplicated currents forced overflow inside the solver; both MPI field ranks recovered, and correct-input retry matched the normal trace. Direct live rollback/retry also matched.
- Prescribed radiation FDTD disk checkpoint at field3, process exit, direct restore in a new MPI job and continuation to5 matched all recorded macro/raw E arrays. Restore setup performed0field evolution calls. This check does not claim multi-file ordinary Fortran stream restoration.
- Updated installed-wheel density readout and the opt-in save example passed. Installed MPI source/restart tests matched source-array controls. The existing64-step total/own-vacuum/molecule calculation matched the previous executable/library's complete frames exactly.
- Source selected suite121passed; installed wheel158tests(156passed,2optional skips). Git-only source preparation/build, sdist-to-wheel, inventory tests, source/object eigensolver audits, localized guide/AI synchronization and whitespace checks passed. No CI or new SCF was run.

Executable SHA256: `a5bfe0e185fe10d036cffeceb74dcf4e6563c0025a9ce1e1b266370e2d86320c`.
Library SHA256: `43186f563eb297b5e901d15741211c94cbf68e08a438ba63e00306c11e1beaad`.
Wheel SHA256: `a6b6b80c7a27bfc479dcee3e08e20f7278b452cc5665ca3686d01cdf48f452bd`.

Detailed local records: OneDrive calc2026/salmonpy-c2h2-nearfield-20261125/live-electric-connections-validation.json. Cloud synchronization was not independently confirmed.

## Remaining scope

This connects real density-derived current to prescribed FDTD radiation. Live spatial-field streaming into the molecular Hamiltonian, a coordinated trial/own-vacuum feedback driver, magnetic response, discrete charge-continuity checks and physical convergence remain open. No complete spatial bidirectional-feedback claim is made.

## Execution incident

An unnecessary executable `--version` call with empty stdin entered the input reader and emptied four temporary inputs in the old SALMON-v.2.2.2 workspace: .namelist.tmp, .atomic_coor.tmp, .atomic_red_coor.tmp and .reenetrance.tmp. No computational input was supplied. The issue was reported to the user; subsequent executions used dedicated trial directories. A separate stalled read-only git-status query in that old directory was stopped; no calculation job was stopped or modified.

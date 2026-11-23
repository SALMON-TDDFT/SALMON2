# Fast model-owned parser migration

The fast parser now accepts the same optional s_input_parser workspace as the
standard parser. Each workspace owns scratch connections, input diagnostics,
coordinate flags and common default/value snapshots. Registered fast HSE models
use a local isolated workspace; complete HSE settings are still captured in the
existing HSE model owner. Legacy batch/serial input preserves named temporary
files, including those required by MD. Global activation remains serialized;
this does not establish threaded native reentrancy or complete the ordinary
retained GS/TDDFT migration.

The entire fast default block is moved, not replaced by standard defaults (only
trailing whitespace normalized). Deprecated local MLWF alias sentinels are
explicit initializer outputs. Unit initialization, KY resolution and input
validation are byte-identical to the previous fast source. Namelists, conversion
arithmetic and numeric calculation order are preserved. Isolated diagnostic and
coordinate metadata are broadcast; root-only file handles are not broadcast.
Long lines are continued for the portable free-form limit. No new IEEE inquiry.

## Verification on Apple Silicon

- MPI2 actual parser: A→B→A, reused workspace, interleaved input splitting,
  omitted defaults, allocated/unallocated atom storage, different coordinate
  sources/units and preexisting temporary files. Standard parser also passes.
- Every emitted variables.log entry equals fresh legacy parses from the
  unchanged preceding audited fast library. Explicit inverse-length HSE input
  and GPU-batch setting survive reading, and defaults are reset on the next
  input omitting them. No actual GPU execution is involved.
- Six ranks, internal MPI2: full/reduced heterogeneous model registry matches old
  independent current/energy/history, fresh state and direct restart exactly.
- Four ranks: LD/SBE/two full-k HSE models,128 steps, same result as prior runtime;
  direct save8/resume120 gives identical continuation with zero replay.
- Serial fast HSE: Taylor/PT-CN learned24-step current/energy and save8/resume16
  exactly match old controls. Same source builds host MPI HSE, serial HSE and
  MPI without HSE. Sixteen no-HSE SBE configurations are byte-identical to the
  preceding no-HSE executable and agree with Python currents.
- Source suite628 methods passes, one skipped. Distribution/installed checks are
  recorded separately with this document. No new test method or Python module.

All new solver, external adapter and parser-probe objects pass the mandatory
workspace eigensolver checker; checker SHA equals the repository checker SHA.
Historical/fixed builds are unchanged and verified by their recorded SHA.
No new SCF/GS collection, actual convergence campaign or GitHub Actions.
These are synthetic regression fixtures, not physical-convergence validation.

## Preserved development evidence

Initial actual-parser compilation failed before the workspace existed. The first
migration build found local deprecated aliases missing from the extracted
initializer; they are now explicit arguments. MPI testing found isolated metadata
had to be broadcast, as in the standard parser. A test initially specified the
fixed history length incorrectly, then used an incorrect inverse-length expected
conversion; those test conditions were corrected without relaxing engine checks.
Failed source/build/test logs are retained beside the passing evidence.

Evidence is stored locally at calc2026/salmonpy-hse-owned-parser-20261119 without
GS/checkpoint payloads. Cloud synchronization and other compiler/platform results
are not independently verified. Effective inventory counts remain unchanged:
587 identical,70 differing,13 standard-only,75 fast-only. Five standard HSE files
remain excluded under the fast-PR-only policy.

Final distribution check: clean installed wheel also passes628 methods (one
skip). Installed preparation reproduces every recorded solver source hash;
installed MPI2 parser ownership/old-log comparison and MPI6 heterogeneous HSE
registry controls pass. HSE-off executable has no FFTW dynamic dependency.
JP/EN guide snapshots match manual.json; diff whitespace check passes.

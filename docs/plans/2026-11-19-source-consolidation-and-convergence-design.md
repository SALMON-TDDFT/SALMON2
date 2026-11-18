# Single maintained solver source and SALMON-specific convergence workflow

Status: approved; first implementation batch underway. User direction 2026-11-19 supersedes
continuing SCF-checkpoint work as the next feature. Preserve the pending verified
SCF-progress/MPI-admission changes; do not delete sources or alter fixed binaries.

## Current evidence

Local file comparison excluding Python bytecode: standard670 files, fast-HSE
reference723;580 identical,69 different shared paths,21 standard-only and74
reference-only. These counts describe paths, not semantic equivalence. In
particular, standard-only files include retained GS owners, embedded TDDFT/EH,
model input capture and LD/SBE shared kernels. The current HSE maintained path
also applies complete stage-specific replacements; raw fast reference alone is
not the current functional engine. NLCC/XC fixes must be compared in effective
prepared source trees, not inferred from names or praise in an external review.

Top-level Python modules70;19 private modules excluding __init__. GS archive
helpers are145/212/158 lines plus an80-line GS-to-TDDFT adapter. Small files alone
are not a defect; organize around independent responsibilities/public boundaries.

## Source consolidation

Maintain one current complete solver tree, with fast-PR HSE as the HSE algorithm.
Use separate build profiles with HSE off/on, not separate maintained source trees.
PZ, TB-mBJ, Libxc-PZ/PBE and non-HSE execution remain available. Non-HSE users
must not acquire HSE-only dependency requirements. No conventional-HSE fallback.

Build a migration inventory of effective prepared sources and their inputs,
including NLCC/XC, symmetry, state ownership, direct checkpoints and recoverable
failure handling. Port standard-only functions and compatible fixes into the
single maintained tree. Retain independently pinned old binaries for comparisons;
never silently replace them. Retain licenses, origin commits and inventories.
Historical references remain in history/evidence, outside the normal current
source distribution once the migration is validated. No mass-rewrite of history.

Acceptance: no-HSE and fast-HSE complete builds pass prohibition audits; existing
LD/SBE/FDTD/GS/TDDFT APIs, mixed MPI/direct checkpoint/rollback controls remain
compatible; new runtime identities are explicit. Do not infer reference safety
from file names. No reference execution containing prohibited drivers.

Only then remove the duplicated current tree/preparation branches from normal
maintenance/distribution. Sample/pseudopotential copies follow the same inventory
and licensing policy; source-reference evidence is not casually deleted.

## Python organization

Preserve public imports and protocols. First consider combining archive loading,
export and MPI array assembly under one internal GS archive module, retaining
GS-to-TDDFT as a distinct lifecycle boundary. Verify the import/call graph before
moving code; do not create replacement forwarding modules indefinitely. Shared
runtime mechanisms may consolidate after the solver-tree migration simplifies
branches. No blanket merge based solely on line count or private names.

Curate three beginner examples (LD response, synthetic SBE response, Si current)
and separately index advanced/native/MPI examples. Mark the Si smoke example as
an execution check, not a physical-convergence certificate. Existing examples
remain callable during reorganization.

## SALMON-specific convergence workflow

One coherent public convergence entry/module, optional callback/backend support,
NumPy as the sole required scientific dependency. Runs and reports are separate:
users can inspect/refine a plan offline without starting native code or compiling.
Model/material IDs, PP hashes, XC, parallel settings, actual grids and runtime
identities follow each trial. Memory budgets remain advisory.

1. Spatial scan with very few fixed k points. Start from PP radius/3; decrease
   requested spacing to approximately0.1 angstrom by default. This is a default
   investigation limit, not a forbidden finer mesh. Record radius definition,
   channels/species, source units and any explicit override; do not mistake radial
   table extent or mesh step for a PP cutoff radius. Follow native KY/UPF radius
   semantics; missing/ambiguous radii require an explicit value, never a guess.
   In multi-species systems expose/record the smallest selected radius. Exact
   channel-reduction choice remains to be agreed before coding the parser.
   Map spacing to integer grid counts by ceiling through the existing geometry
   code, report actual per-axis spacing, and avoid duplicate effective grids.
   Start with GS energy per atom and optionally gap/user observable. Keep cell,
   atoms, XC, PP, bands and SCF accuracy fixed. SCF nonconvergence warns and does
   not halt the campaign; mark the affected convergence evidence provisional.
   k-point convergence is not established by this inexpensive spatial scan.
2. Freeze the selected grid and initial GS; scan dt from about0.2au downward
   (default halving). Use the same weak impulse/polarization and a few-fs physical
   duration, including x/y/z and arbitrary user polarization. In crystal velocity
   gauge the impulse is a jump of A/c with fixed integrated E, not a dt-dependent
   fixed E sample amplitude. Declare exact sample/time conventions and fixed
   excitation strength. Create independent RT states for every trial.
3. Prefer streaming step execution. Check returned observables and available
   native numerical status each step. NaN/Inf/numerical failure ends that trial
   immediately, preserves diagnostics and moves to smaller dt. Do not stop all
   trials, relabel a failure as zero, or use a partially failed trace for an error
   comparison. Generic I/O/worker/resource failures are recorded separately from
   numerical instability; resource measurement failure is not calculation failure.
4. Report tested nonbreaking dt separately from converged dt. A coarse dt that
   survives a few fs is evidence for this test window, not unconditional stability.
   Compare current traces on common physical times (nested dt preferred) with
   absolute and relative errors, and optionally spectra using identical analysis
   conventions. Need at least consecutive finer finite traces for comparison.
   Report all errors and user-visible tolerances; near-zero responses need an
   absolute tolerance floor. Do not invent a universal material-independent
   convergence tolerance or call the finest tested dt automatically converged.
   No converged selection if no trial meets the declared tolerance. Stopping at
   the default lower spacing/dt bound reports unresolved convergence honestly.
5. k-point convergence is a later independent stage. A low-k grid/dt diagnostic
   is not a complete optical-response validation, especially for other bands or
   excitations. Recheck the final conditions for production use.

Convergence driver acceptance tests use synthetic callbacks with known grid/dt
errors and controlled NaN/Inf onset, verifying trial-local stop, independent
states, physical-time alignment, no fake convergence, and warning-only SCF/resource
policy. Real-material campaigns need explicit selected material/observable/radius
and are not started by this design document. No new GS data collected now.

## Optional dependencies and priority

Do not automatically add the external review's list. Current process_memory uses
resource on Linux/macOS, not /proc. psutil could add current process/child RSS as
an optional extra; summed RSS is not a true aggregate high-water mark and can
count shared pages twice. Failure remains advisory. spglib reduction is separate
from k convergence and must respect the Hamiltonian/external-field symmetry;
never enable it silently for general polarization or fast HSE. Real zarr
interoperability, formatting/linting and property-based tests are separate scoped
changes. Existing unittest parametrization does not itself require pytest.

## Proposed implementation order

1. Verify effective-source parity/feature matrix and consolidate the solver base.
2. Consolidate GS archive internals and curate beginner entry points, preserving
   the public API and measured behavior.
3. Implement convergence planning/reporting and synthetic campaign tests.
4. Connect existing GS/RT backends and per-step numerical stop; verify a selected
   material within an agreed calculation budget. Complete full SCF continuation
   after this user-prioritized work unless redirected.

## First batch implementation

- tools/runtime_sources.py now resolves verified complete-file replacements in
  order, without copying/mutations; CLI emits effective-source origins and hashes.
  Standard670 and fast-HSE725 inventories exactly match independently prepared
  directories. Their effective comparison is573 identical,77 different,20
  standard-only and75 fast-HSE-only paths. Raw-reference counts above are a
  different comparison and remain historical evidence.
- GS archive load/MPI export/serial export now share _gs_archive.py; old three
  private modules removed, public entry points unchanged. Function AST comparison
  permits only the private export rename and removal of a self-import.
- Beginner index explicitly selects LD, synthetic SBE and the Si smoke example,
  with paired guides/AI snapshots.
- Full solver source consolidation is NOT completed by these steps. Port the
  verified standard-only state/bridge functionality and reconcile77 differing
  files before deleting either current tree. Host-side native_adapters are a
  separate build input and must also be retained in the migration inventory.

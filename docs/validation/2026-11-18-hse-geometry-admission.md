# Native fast-PR cell/atom/grid admission (2026-11-18)

Previous symmetry admission committed/pushed as6479a538 before this work. Python
only: no fixed Fortran binary, solver formula or shared kernel is modified.
Scope: existing native fast-PR periodic scalar CPU profile, one/two internal k
ranks, atomic-unit input and axis-aligned orthogonal cells. No general native
initialization recovery, restart authentication or physical convergence claim.

## Checks before native lifetime/work copy

- Focused numeric namelist reader retains assignment ordering, D literals,
  repeats, positive index/section forms, quoted path case and comments. Reject
  null/unsupported forms rather than guessing. Quoted numeric atom fields and
  invalid integer forms reject before native I/O.
- Match bundled resolve_cell priorities: al, diagonal vectors, or counts*dl.
  Match init_grid_whole_config rounding/tolerance and native integer limits,
  including transformed grid indices. Do not ceiling or normalize input.
- Validate finite positive effective cell lengths/volume and inverse; actual
  grid spacing/counts; one explicit inline/file Cartesian/reduced atom source,
  atom count, finite coordinates and species IDs1..nelem.
- With symmetry, check Cartesian metric, affine grid residual/index bounds over
  all eight corners and unique same-species atom permutations using native
  tolerances. No full-grid allocation is needed for admission.
- Preserve the original yn_symmetry spelling when validating canonical flags:
  the older generic scalar reader lowercases strings and is inappropriate here.

All create_state solver ranks collectively agree on admission errors. A rejected
call acquires no native lifetime or work copy. Repair input and construct a new
provider, keeping immutable-provider hash checks intact.

## Evidence

Local Apple Silicon/macOS, GNU Fortran15, OpenMPI5, Libxc5.2.3. Reuse the unchanged
explicitly built MPI library from the preceding validation. No GitHub Actions.

- Unit regressions cover valid scalar/index/repeat forms, cell fallback,
  malformed/singular/nonfinite cells, grid counts/spacing/overflow, large-grid
  accumulated mapping errors, transformed-index overflow, coordinate count/
  species/finite values, quoted numeric rejection, coordinate file case and
  read-only input/file preservation. Backend test forbids CDLL and work-copy
  creation on invalid geometry. Final source suite569 tests, one skip.
- MPI4, two independent groups of two internal k ranks: cell, cell-metric,
  cell-orientation, atoms, grid and grid-translation faults affect only the last
  rank of each group. All ranks reject before lifetime, repair and start in the
  same process. All normal current/energy traces match the prior32-step run
  bit-for-bit. Fresh save7→32 resume matches its continuous traces/final checkpoint
  bytes without accepted-step replay.
- Normal mixed FDTD + two local LD + SBE + HSE2 and former PBE/HSE Taylor/PT-CN/
  combined five/seven-rank regressions pass unchanged tolerances.
- Equivalent spacing-only grid input (dl rather than counts) and case-sensitive
  external Cartesian Atoms.dat input in au run with full/reduced MPI2 groups;
  both current and energy match the former input exactly.
- Source-distribution→wheel/native-hash checks and final installed Python/MPI
  verification are retained with the evidence bundle. New private module is
  explicitly installed. ASE/f90nml are not added as dependencies.

## Measured limitation and scope correction

Exploratory axis-permuted right/left-handed synthetic input variants reached
native ACE initialization STOP in the existing binary. This does not establish
support merely because the independent exchange core accepts general cells.
The admission profile therefore rejects off-axis/nonorthogonal vector inputs
before lifetime, without rotating/reordering a basis or reusing unmatched GS
states. These failed logs are retained; general-cell native initialization stays
open. A_eV_fs is already rejected by FastHSEMaterialProvider; retain its au-only
contract rather than extending units through this validation change.

K-star consistency, pseudopotential contents, GS identity/restart binary layout
and other initialization failures remain separate work. All ranks must enter the
collective call; provider construction itself remains a noncollective operation.
Cross-initialization checkpoint workspace bytes are not claimed invariant; direct
resume comparisons use the matching continuous trajectory.

Use the existing explicit MPI save probe with --group-size 2 and
--admission-fault cell|cell-metric|cell-orientation|atoms|grid|grid-translation.
Repair/retry is exercised before normal propagation; no compiler/download is
started by importing or running admission.

Local evidence:
`/Users/otobetoshihito/.codex/outputs/salmonpy-hse-geometry-admission-20261118/`.

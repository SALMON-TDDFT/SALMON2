# Standard cell/grid initialization ownership — 2026-11-17

Continuation of the approved remaining GS/TDDFT initialization work. Fast-PR HSE
already had explicit cell/grid resolution. This change ports its existing value
types/resolvers to the standard solver and optional RT-step structure profile;
no new numerical scheme or public Python/C API is introduced.

## Ownership scope

`s_cell_config`, `s_cell_geometry` and `s_rgrid_config` are value-only types.
Existing derived-type layouts are unchanged. Both profiles use capture-only
`cell_config_from_input`; explicit `resolve_cell`/`init_grid_whole_config` read
no global input. Valid legacy cell precedence, floating-point arithmetic,
periodic/isolated odd/even indices, epsilon thresholds, Maxwell consistency
lengths and downstream `al` normalization are preserved.

Rejected cell/grid configurations preserve the target geometry/index arrays.
The explicit grid helper builds whole-grid indices, not every coordinate,
parallel communicator, symmetry or FFT resource. Compatibility wrappers still
capture shared input. This does not complete different-input GS/TDDFT contexts
in one process or general native initialization-failure recovery. Molecular work
remains paused and memory budgets remain advisory.

## Verification

- Existing cell probe fails to compile against old standard modules (missing
  types/resolver). New standard and fast-PR cell/grid/ion probes pass OMP1/2;
  probes compile with bounds checking. Frozen original grid calculations cover
  periodic/isolated grids, odd/even placement, epsilon boundaries, Maxwell and
  explicit versus inferred spacing. Capture-then-poison-global checks confirm
  independent configuration values; invalid cells/grids leave originals unchanged.
- Standard and HSE cell/grid resolver bodies match; the port preserves existing
  type layouts. Standard inventory, RT-step before/after replacements/base pin,
  HSE replacement and upper manifest SHA256 links are refreshed. Strict fresh
  preparation passes for all three source profiles. Shared LD/SBE kernel bytes
  remain unchanged.
- Local Apple Silicon64GB/18-core host, at most seven MPI ranks, OMP/BLAS1,
  48GiB advisory budget. Fresh MPI2 GS warm-start three-iteration energy/info,
  eigenvalue, k-point and force files are byte-identical to the unchanged old
  library's executable. Comparison files are removed before each run. This is
  startup/numerical regression, not physical convergence.
- PBE TDDFT and fast-PR HSE Taylor/PT-CN five-rank mixed32-step late/deep failure,
  state/output rollback and valid retry pass. Seven-rank PBE+HSE mixture also
  passes. Final response arrays and every ordinary stream are identical to
  corresponding unchanged normal baselines.
- Seven-rank direct save16/late restore failure18 preserves original numerical
  bytes, allows original continuation and matches retry/three future frames
  without elapsed-history replay. Identity/admission negatives still pass.
- Source and fresh installed wheel suites545 tests, one skip each;
  source-distribution→wheel/native SHA256 checks and JP/EN/AI synchronization pass.
  Installed seven-rank mixed and five-rank1D two-TDDFT-group protected fault/retry
  cases also pass against unchanged baselines.

The optional RT-step profile also passes a full fresh four-worker build with
matching modules and all cell/grid/ion probes OMP1/2. This validates preparation,
build and explicit initialization helpers; it is not a full RT-profile physics
or restart campaign.

New standard/HSE libraries compile structures and initialization into NEW outputs,
reusing unchanged audited objects only where existing layouts/interfaces match.
The optional RT-step profile has a different s_rt layout: a mixed-profile partial
compile correctly rejected that mismatch. Full RT-profile build uses one matching
module tree. Early full configure attempts selected AppleClang through the default
MPI wrapper and could not find OpenMP_C; the existing GNU/OpenMPI build recipe
explicitly selects its GNU C compiler. These setup failures are retained separately
from successful scientific comparisons; they are not accepted solver results.

## Evidence

Persistent local evidence:
`/Users/otobetoshihito/.codex/outputs/salmonpy-cell-grid-initialization-20261117/`.
Contains prepared source archives, production receipts/libraries, actual build and
probe commands, independent baseline comparisons, output inventories, direct
restore evidence, installed tests and distributions. No remote Actions were used.
Local paths are provenance, not user installation requirements or portable runtime
release promises. Wider k/symmetry/communicator/initialization ownership and stable
native handles remain separate work.

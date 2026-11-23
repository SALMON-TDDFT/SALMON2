# Retained serial regular-k portable export — 2026-11-19

Preceding commit b66b97fd (and2340787d) pushed to
salmonpy/codex/salmonpy-state-inputs-recovery; remote SHA confirmed.
Commits use [skip ci]; no manual workflow dispatch.

NativeGSState.export_archive now covers immutable DFT-generated regular full
k meshes in serial subgroups. Native lattice.f90:init_kvector ordinary-case
order is x fastest, then y/z; coordinates (i+.5)/n-.5 for zero-based i,
with uniform1/nk weight. DFT.input_text does not select custom k files,
shifts, Gamma centering or symmetry reduction. No native solver/update was changed.
The Python export validates full native orbital/occupation/eigenvalue shapes.
Density is reconstructed in k/band order using weights and occupations.

Synthetic complex noncubic8×10×12 arrays and distinct bands/k points cover
2×1×1,3×1×1,2×3×1,2×2×2. Every orbital, occupation and eigenvalue survives
readback exactly. Coordinate/weight/electron contracts pass. Density agrees
with weighted abs(u)**2 to rtol1e-15, including rounding differences.
Closing/mutating source owners leaves archived arrays independent.
Both real/complex Gamma synthetic archives compared against the b66b97fd
implementation: every archive file's SHA256 is identical.
Source suite623 methods, one skip; focused GS/archive18 methods pass.

JP/EN guides, AI search fields/embedded guides, optional sample and checklist
are synchronized. The sample preserves its ordinary two-k-point models even
when exporting; --help now parses before MPI initialization. Help/sync/
whitespace checks pass. Optional SCF sample was not run.

No new SCF, real GS export/collection, native rebuild or wheel rebuild.
These tests do not establish physical native export/TDDFT parity. Retained
MPI/custom/reduced-k export, imported-owner physical identity, native GS
handle restoration and a full mixing-history checkpoint remain future work.
Native GS/2 import/conversion retains independent runtime/cell restrictions.

Selected logs/hashes only (no GS arrays):
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-retained-gs-multik-20261119
Local writes verified; cloud synchronization unverified.

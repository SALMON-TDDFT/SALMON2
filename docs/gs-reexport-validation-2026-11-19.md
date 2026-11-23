# GS/2-imported owner re-export — 2026-11-19

load_archive captures only small immutable source metadata/geometry/k/eigenvalue
arrays and potential objects. No source wavefunctions/live archive is cached.
A live imported owner reuses export_archive: native orbitals/occupations are read,
halos removed, density reconstructed with source sampling weights. Source geometry/
model/k data, eigenvalues and energy are carried and labelled rather than claiming
new native SCF computation. New retained-gs-reexport/1 provenance records the source
manifest SHA and inherited convergence separately; current convergence remains None.
calculation stays None and uncomputed forces are not enabled. Raw legacy owners
are not guessed from partial metadata. Old export provenance remains supported;
new provenance needs an updated reader.

Complex C-array reads previously assembled real+1j*imag, normalizing signed-zero
bits and potentially performing invalid FP arithmetic before collective rejection.
Component copies now preserve finite bits and reject NaN/Inf collectively even
with NumPy strict FP settings. No Fortran solver/update order changed. Synthetic
source-removal tests prove re-export works without original files; a native-vs-
source sentinel test confirms source eigen/energy values do not become fabricated
native measurements and native occupations take precedence. Original uncomputed
forces stay unavailable. Existing computed Gamma synthetic archives remain byte-
identical to b66b97fd (both real/complex fixtures); signed-zero getter fidelity is
an intentional correction where the former arithmetic lost those bits.

Seven actual MPI4 synthetic cases: PZ orbital/spatial/k, PBE orbital, complex Gamma
orbital/spatial, irregular k/nonuniform weights. Orbital payloads (including signed
zero) preserve bytes; occupations and carried model/eigen/k/energy data agree.
Reimport creates a second native model; its16-step current equals both the first
model and independent file controls. Independent state and direct8+8 checkpoint
resume agree with zero replay. Partial write failure at the third publication
while RT is active leaves GS resources/current unchanged, removes incomplete
output and permits exact retry. Later RT evolution still matches controls.
Comparisons are within the same MPI decomposition, not across decompositions.

Source628/installed628 tests pass (one skip each). Installed PBE/Gamma/nonuniform
actual MPI cases pass. Installed public sample calls twice with reexport_to in
one session, saves each material ID correctly and preserves the previous file-
current reference hash. JP/EN/AI guides/sample/checklist synchronized; offline
sdist/wheel bundled-source/object audits and native hashes pass. No new electronic
runtime, SCF or real GS data collection; the existing validated converter/runtime
were reused by recorded SHA. Measurement: Apple Silicon/GNU/OpenMPI/Python3.14.

This is labelled array reinitialization, not full SCF mixing/cache/iteration
restoration or migration of an existing RT checkpoint to a regenerated model.
Physical real-material convergence, raw legacy binding and other platforms remain
outside this evidence. Selected logs/hashes/wheel only (no GS arrays):
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-gs-reexport-20261119
Local files verified, cloud sync unverified. No push/Actions; pending work after
d67ce847, including prior import/metadata slices, remains uncommitted.

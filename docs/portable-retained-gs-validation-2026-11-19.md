# GS/2→retained native GS→TDDFT — 2026-11-19

NativeGSSession.load_archive explicitly connects a portable array archive to the
probed converter and existing native load. Root alone verifies/converts; members
agree selected arguments/environment, then initialize a new native GS owner.
Archive IDs and recorded convergence claims are retained separately from current
import convergence(None) in immutable archive_origin. Missing old compiler/source
attestation warns and proceeds; recorded incompatible markers/integer widths
reject before writing. New incomplete conversions clean up, existing targets
are preserved; completed bundles remain caller-owned. No implicit build/SCF.

A compiler-bound converter was built against the existing SHA-validated serial
runtime. Mandatory workspace source audit passed before extraction/build; the
builder checker is byte-identical to the workspace checker. All linked application
objects audited before linking/probes and again before evidence transfer. The
extracted native read/write routines remain unchanged; byte-parity and multi-k
probes pass. IEEE predicates in converter/exporter/stream sources are scalar.
No electronic SALMON solver or shared library rebuilt; none of its sources changed.

Actual native tests use synthetic orthonormal cell-periodic orbitals, two bands
and an8³ periodic cell. PZ two-k tests cover orbital2/spatial2/k2, PBE covers
orbital2, and complex Gamma covers orbital2/spatial2 (four total ranks per test).
All native local orbital tiles/replicated occupations equal source arrays exactly.
GS→memory TDDFT16-step currents equal independent converted-file controls, and
independent/sibling-state controls, checkpoint8+resume8 and failed-conversion/
retry checks pass with zero replay and zero SCF. Comparisons use the same MPI
partition on both sides; no cross-partition current-bit equality is claimed.

Gamma required an equivalent explicit k table: DFT num_kgrid=1 would allocate
real orbitals while the GS/2 writer preserves complex amplitudes. load_archive
stages that table and receipt hash before native import. Tests include nonzero
imaginary components and confirm exact preservation. Existing Gamma-to-TDDFT
converter behavior is unchanged. The earlier wheel candidate remains preserved
in private/tmp; only the later verified wheel includes this complete Gamma route.

Source627 and installed627 tests pass with one skip each. Installed PBE and
complex-Gamma native MPI cases pass. The installed public current_from_archive
sample called twice in one existing session matches the independent file-current
hash both times and releases owners. JP/EN/AI guides, sample, runner and package
manifests synchronized. Offline sdist/wheel source/object audits/native hashes pass.

This initializes orbital/occupation arrays and reconstructs native potentials/
density. It does not restore archived energy/eigenvalues, SCF mixing/cache/iteration
history or an arbitrary RT checkpoint on changed topology. Imported calculation
remains None and forces unavailable. General native STOP/MPI/process failure
recovery, physical real-material convergence and other platforms remain unverified.
No actual GS dataset was generated/collected; source arrays and native records
in tests are synthetic. Host: Apple Silicon/GNU/OpenMPI/Python3.14.

Selected logs/hashes/converter binaries/wheel only, no GS array/restart collection:
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-portable-retained-gs-20261119
Local storage verified; cloud sync unverified. Converter build.json is original
build evidence with original paths, not a relocated deployment catalog. No push
or Actions this slice; pending work afterd67ce847 remains uncommitted.

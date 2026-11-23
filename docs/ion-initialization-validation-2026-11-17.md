# Owned ion initialization inputs — 2026-11-17

First slice of the approved remaining GS/TDDFT initialization ownership work.
No Python API or stable C ABI is added. Standard and fast-PR full-source
`initialization_sub` now carry identical value-only `s_ion_initial_config`,
legacy capture `ion_input_config` and explicit `initialize_ions`.

## Contract and limits

Captured atom/element counts, Cartesian/reduced positions and species IDs are
independent allocatable values. The explicit initializer checks allocation,
shape, one-based lower bounds, species range and finite input/resolved positions
before changing the target. Reduced coordinates use that target's primitive
lattice. Numerical ion arrays are independently owned; no FFT/MPI/solver resource
owners are copied. Invalid input retains the original positions/species/counts,
mass allocation and velocity. Velocity is zero-initialized; mass/force retain the
legacy policy of being populated by their owning physical routines before use.

`init_dft_system` captures once and uses the same initializer. Its compatibility
wrapper retains reduced-coordinate normalization of the legacy global Rion for
remaining downstream consumers. The explicit helper never reads/writes shared
inputs. This does NOT remove cell/k/symmetry/communicator and other global
initialization dependencies or establish arbitrary full GS contexts in one
process. Allocation failure/global initialization STOP recovery remains open.

## Native tests and provenance

`tests/native/test_ion_initialization.f90` failed to compile against the old module
(missing type/helpers), then passed standard and fast-PR builds with OMP1/2.
Two different configurations survive poisoned shared counts/coordinates/species,
produce exact Cartesian/nonorthogonal reduced-coordinate references and can be
alternated repeatedly without changing the other system. Invalid shape/species/
NaN/unallocated/zero-based arrays reject without mutation. Read-only review found
the lower-bound case; explicit rejection and regressions were added before final
verification. IEEE predicates are called on scalars only.

Complete standard source inventory and fast-PR replacement hashes are refreshed,
including the standard overlay base pin. Fresh external source preparation
verifies these manifests; no compatibility patch or fixed binary was modified.
New libraries recompile only initialization_sub against unchanged audited core
objects and link to NEW output paths. Existing external builds/binaries remain
unchanged. Source/build receipts, actual commands and unchanged-core provenance
are retained. Shared native LD/SBE files remain byte-identical to the solver.

## Runtime comparisons

Local Apple Silicon64GB/18-core host, OMP/BLAS1, up to seven MPI ranks,48GiB
advisory budget, sequential MPI jobs. No remote Actions.

- Standard Libxc-PBE GS, internal MPI2, warm-start three-iteration regression:
  old/new energy/info, eigenvalues, k-point and force files byte-identical.
  Target output files are removed before running to establish fresh generation.
  This short test is not an SCF/physical-convergence certification.
- Five-rank FDTD2 + local LD + SBE1 + internally parallel electronic2:
  Libxc-PBE TDDFT and fast-PR HSE Taylor/PT-CN,32 accepted steps. First/warm/final
  late/deep failures preserve original native numerical bytes/Python states and
  ordinary output. Valid retries match prior unchanged binaries' normal arrays
  and every ordinary `.data` stream exactly.
- Seven-rank mixture includes separate PBE2 and HSE PT-CN2 groups plus SBE1 and
  field2/local LD. The same numerical/output rejection and retry comparisons pass.
- Seven-rank mixed direct checkpoint at16, injected late restore failure at18:
  all original state bytes preserved, original continues, retry plus three future
  steps exact without elapsed-history replay. Ten identity/admission negatives
  still reject before native load.
- Python source/fresh installed suites545 tests, one skip each.
  Source-distribution→wheel/native SHA256 validation passes. Installed seven-rank
  mixed and five-rank 1D two-TDDFT-group protected retry arrays/ordinary output
  also match unchanged normal baselines exactly.
  JP/EN API guides and AI JSON reflect the limited internal ownership scope.

## Reproduce and evidence

The native probe is shipped under `_runtime_source/tests/native/`; explicit
compiler/module/library commands are retained in the build receipt. Source
preparation uses `tools/prepare_salmon_runtime.py` and
`tools/prepare_hse_rt_state.py --explicit-routing` into new external directories.
Runtime cases use the existing native protected-step and mixed protected-restore
MPI tools with matching newly built libraries; nothing builds on import.

Persistent local evidence:
`/Users/otobetoshihito/.codex/outputs/salmonpy-ion-initialization-20261117/`.
This records local provenance, not portable installation requirements or a new
published runtime release. Wider initialization ownership/handles, platform
validation, generic archives and transferable resource models remain separate
work. Molecular work stays paused and budget excess only warns.

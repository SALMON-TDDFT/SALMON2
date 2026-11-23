# Distributed-field mixed-material validation — 2026-11-16

Extends the preceding native mixed fixture without solver changes. The sample adds --field-grid PX PY PZ and uses disjoint field/SBE/TDDFT rank ranges. It also supports --compare for another same-physics run, validates metadata,point order,array shapes and finite reference values,then compares full A/c,E,j histories with rtol2e-10/atol2e-12. These are fixed-sample checks, not a general archived-model compatibility API.

Apple Silicon/OpenMPI/GNU/OpenBLAS; reused audited native library and managed serial comparison runtime. Four MPI ranks,one OpenMP/BLAS thread each;128 steps,dt=.003,40³ field,synthetic LD/SBE and real Si starter GS for TDDFT. No builds/downloads on import or GitHub Actions.

- Field decomposition(2,1,1),(1,2,1),(1,1,2),with SBE1 and TDDFT1 ranks,all pass serial-material replay and native half-step-clock assertions.
- All A/c,E,j values match the preceding one-field-rank/one-TDDFT-rank baseline exactly (zero maximum difference),not merely within the allowed comparison tolerance.
- x splitting separates the two electronic point cells across field ranks. y/z splitting divides each9-cell LD material into3 and6 constitutive cells across the two blocks,exercising both field partitions with both LD materials.
- Source distribution -> wheel -> install -> four-rank y-split execution with --compare passes. Installed report explicitly records zero comparison differences.
- Independent read-only review found no rank overlap,point-order or comparison-target error.

Scope: field MPI and inner TDDFT MPI are separately verified at2 ranks. The simultaneous field2+SBE1+TDDFT2 five-rank configuration is not tested under the four-core budget. Large-cell scaling,aggregate RSS enforcement,whole-session restart,other platforms and physical convergence remain open. This tiny3D fixture does not change the recommendation of1D propagation for practical TDDFT.

Persistent evidence root: `/Users/otobetoshihito/.codex/outputs/salmonpy-state-inputs-20261116` (local provenance only).
- `native-mixed-one/result.json`: `90ffb11cf59413b196abfffeaa092f048c87f79c65683ee3cd47edec6824a31e`
- `native-mixed-field-x/result.json`: `97f9dfe51a52d48e125f502fcce975f4bc5799e1cbea20858a8921f6b1191f9f`
- `native-mixed-field-y/result.json`: `526b09a19397a6a922c71fb2170fc5b68cf07dfe337c3779a532786bb7bc6e99`
- `native-mixed-field-z/result.json`: `7b67560893a035d9005a37dd17512a37ff50d9e6a639e1dd33dcb3f22c7c6b0b`
- `native-mixed-field-installed/result.json`: `594491b88917b76a5fe060b632df9926a67cf8031f01af88d5f01122e963ace4`
- `field-mixed-dist.log`: `359babf10e4acd573724b7a27918f811d2c4900ac3456d216741d39f24aa632b`

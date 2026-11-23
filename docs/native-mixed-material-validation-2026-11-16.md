# Native field-local LD / SBE / internal-MPI TDDFT — 2026-11-16

A new explicitly built standard host-MPI library uses native E/H and TDDFT, PZ, no Libxc/HSE. Mac AppleClang behind mpicc initially lacked OpenMP; rebuilding in a new directory with OMPI_CC=/opt/homebrew/bin/gcc-15 succeeded. Toolchain paths are local evidence, not portable user defaults. Source/object forbidden-eigensolver audits passed; fixed binaries and pinned sources were not patched.

The installed-wheel mpi_mixed_native example runs128 steps at dt=.003 on a40³ field grid. One field rank owns native LD regions of9 constitutive cells each (Drude and two Lorentz poles). Separate SBE point10 occupies one rank; Si TDDFT point20 occupies a subgroup with one or two orbital ranks. The Si starter GS is the real54-iteration PZ GS from the preceding managed-worker test. This is a tiny3D connection fixture, not converged optics; production TDDFT still recommends1D propagation.

Acceptance:

- 4-rank run: field1 + SBE1 + TDDFT2. Both electronic currents are nonzero. SBE drive replay matches its original serial state exactly; native TDDFT replay against an independent managed serial worker differs by at most7.66e-19 a.u. (rtol2e-10,atol2e-12).
- 3-rank run: field1 + SBE1 + TDDFT1. Both serial replays match exactly. Full-history decomposition comparison with the4-rank installed-wheel run passes.
- Removing LD regions changes observed E and electronic j; local LD does participate rather than being a vacuous zero-response setup.
- Native LD-only regression matches the original batch-solver observation files exactly, preserving the existing constitutive update behavior.
- Feedback uses the preceding j, A advances by -dt E, and the electronic field clock is (n+1/2)dt. All128 steps are checked.
- Source distribution -> wheel -> installed four-rank execution passes. MPI failure cleanup was corrected after review: abort directly on error; only successful execution enters collective close operations. No recoverable-error guarantee is added.

Maximum1/2-TDDFT-rank differences: {"A_over_c": 1.1934753938724234e-20, "E": 4.336808689942018e-19, "j": 7.657177843178875e-19}.
Maximum LD-present/vacuum differences: {"A_over_c": 4.3502891187919523e-07, "E": 1.8713622656201022e-05, "j": 1.2088520445265873e-08}.

Scope limits: one field rank; small LD regions only; no many-field-rank scaling, whole-session snapshot/restart, GPU, other platforms or physical convergence. Each rank uses one OMP/BLAS thread; at most4 MPI ranks run at a time. Aggregate RSS was not measured/capped. No remote CI or push.

Persistent evidence under `/Users/otobetoshihito/.codex/outputs/salmonpy-state-inputs-20261116` (provenance, not a required installation path):
- `mixed-native-build-v2/build.json`: `d2b0569e6fad50f1deef1f6835b21ded28f72cb6e71f3b7db9a840fce50632ff`
- `native-mixed-one/result.json`: `90ffb11cf59413b196abfffeaa092f048c87f79c65683ee3cd47edec6824a31e`
- `native-mixed-installed/result.json`: `1fa5e4af375a047d52cf03c022f8c3480e369c54fc5767c0b511933a7b70bfc5`
- `native-mixed-vacuum/result.json`: `d6319e50332af0bdcb51dad3a13554a28b5d9b29f809b9ae4de4070334211f17`
- `native-ld-regression.log`: `2ed821c0dd95063c3a5c9c160dc00d0ff8acfd2ff41bfd4925fbc47808b9904b`
- `native-mixed-dist.log`: `1b682bafe1320a9b4f29d54aafac16d9239e804b8d02872087e98fe39ec755df`

# Computed fast-PR HSE GS and memory RT — 2026-11-19

Scope: GNU15.2/macOS arm64, OpenMPI, Libxc, OpenBLAS, CPU full-k unpolarized fixed-ion
HSE06. DFT(xc="hse06") now uses the declared fast-PR NativeGSBackend run capability.
Only the native entry is connected; batch HSE reports are not. New HSE descriptors
explicitly disable PBE pre-SCF so the requested iterations evaluate HSE immediately.
Unconverged SCF warns and preserves the last iterate. No conventional-HSE fallback.

A single prepared fast source was built into a fresh HSE-enabled runtime. Workspace
source and all-application-object prohibition audits passed. The previous fixed
executable was reused only after its recorded SHA256 matched. PZHEEV was never run.

Six host ranks form three two-rank k-only solver groups. Two differing materials
(Si1/8bohr/8^3/2bands/2x1x1 and Si2/9bohr/12^3/4bands/1x2x1) each run four HSE SCF
iterations. Retained A+B density/orbitals/occupations/eigenvalues/forces equal
independent A/B arrays exactly. Slot recycling, closing A and interleaving another
SCF preserve B and active RT owners. Native convergence metadata remains explicit.
Export assembly was inspected in RAM without writing GS archives. GS info and
eigenvalue reports are byte-identical to the fixed prior fast standalone executable
with the same inputs. Reports are only local controls; no GS numerical dataset is
copied into this evidence directory. This is not physical-convergence validation.

Source-GS-close RT owners use a three-component A/c drive. Sixteen-step currents
and HSE energies match independent execution exactly; a direct8→16 restart matches
without step replay. Transient ordinary/HSE-computed RT checkpoint controls are
removed after restore. Existing HSE learned full-k Taylor/reduced-k PT-CN controls
also retain exact fresh-state/direct-resume/history agreement with matched linking.

The first HSE release/debug native GS runs crashed in ZDOTC during orthogonalization.
Darwin mpifort -shared generated a flat namespace, resolving zdotc_ from NumPy's
Accelerate BLAS, whose complex return ABI differs. dladdr identified Accelerate;
object/disassembly evidence located the ZDOTC call. The builder now reuses ordinary
mpi_link_command for explicit two-level BLAS binding and MPI COMMON conventions.
Failed binaries/logs are preserved here. Fixed old binaries were not modified.

An unchanged-old-flat comparison failed strict equality with energy differences up
to4.4408921e-16 Ha in the first reported case. Its tolerance was not relaxed. A
separate new library linked from the prior audited objects using the same new
linker policy matches the new runtime exactly in the learned full/reduced controls.
This separates the link/BLAS effect from solver-source changes. The original
four-iteration PBE warmup control also stopped as unfinished; the log is retained.

Ordinary Libxc-PZ/PBE spatialMPI2 computed GS→RT remains exact. Source and installed
Python suites each run628 tests:627 pass,1 skip; no new test methods. Source sdist→
wheel→new environment install, installed HSE sampleMPI2, installed HSE independent/
resumeMPI6, JP/EN/AI/tool file parity, manual synchronization and diff whitespace
checks pass. Failed and final build records retain exact identities. No Actions,
commit or push in this slice. Local OneDrive storage does not certify cloud sync.

Remaining: batch GS output reports, one-session ordinary/HSE registry unification,
final source-tree consolidation and larger physical-convergence/platform campaigns.

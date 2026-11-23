# Fast-PR SBE source sharing and HSE-off build verification

Current fast explicit-routing preparation references five existing standard
files (SBE CMake/types/velocity kernel/GS wrapper/Bloch wrapper) in its complete-
file inventory. No new duplicate Fortran file or patch is introduced. Both
standard and fast preparation call a shared native/solver byte-equality check.
Reference-only snapshots are unchanged. The existing inventory test first failed
for the absent shared type path; it now validates native hashes and origins, and
its copy test covers mismatched kernels. No new Python test method/module added.

The same final prepared source builds with HSE ON and OFF. OFF is configured
without HSE FFTW flags and has no FFTW dynamic linkage. Its builds exposed
unguarded HSE-only imports/status calls and an endpoint-history call in Taylor/
time_evolution_step. These are guarded, retaining ordinary recovery_wave/scalar/
value/agreement checks. With HSE enabled, preprocessing gives the previous
computational statements, excluding only NVTX source-line labels. Source and
every linked application object, including host adapters outside the solver
build directory, pass prohibition checks. No new IEEE inquiry is introduced.

Numerical controls use the unchanged audited hse-layer-io-fixed executable and
library, with recorded SHA verification. Its original solver source hashes
match the pre-change effective fast tree; final differences are only five SBE
paths and two guarded RT files. Host adapter hashes are unchanged.

-16 SBE configurations (correction order0..3,nonlocal correction off/on,active
 bands2/3), each under HSE-OMP1/HSE-OMP2/no-HSE/HSE-MPI2:all pass. Single-thread
 files (current/energy/excitation) are byte-identical to old. OMP2 uses rtol3e-13
 and atol3e-15. Python current equals minus the native matter current within
 those tolerances. Inputs/orbitals are synthetic; no SCF is called.
- Learned Taylor and PT-CN, full and reduced k meshes, internal MPI2 (six total
 ranks) match independent old controls in current/energy/history. Fresh states,
 close ownership and direct checkpoint/resume also match with zero replay. This
 compares each full/reduced model to its own old control, not different-dt
 models to each other.
-628 source and628 installed-wheel methods pass,one skipped each. Installed
 preparation reproduces all compiled solver source hashes. Installed HSE MPI
 independent/current/energy/history/direct-resume and SBE example pass.
- Offline sdist→wheel/source/object hooks,installed native SHA and JP/EN/AI sync
 pass. Electronic runtime builds are local GNU/macOS tests; Fujitsu/Linux/Intel
 portability and real-material convergence are not established here.

Failed build logs for the progressively exposed HSE-off references are retained.
Initial HSE controls using old executed work-directory copies failed initialization;
old-only controls reproduced that failure. The regression oracle instead uses
the pinned original synthetic fixtures referenced in prior validation, which
pass. Those work copies are not treated as validated initial-state oracles; no
claim is made to fully diagnose them in this slice. A failed edit precondition
also left an unchanged intermediate preparation; it is retained, not the final
validated asset. No historical source/fixed binary was replaced or deleted.

Canonical logs/hashes/manifests/distribution and runtime evidence (no collected
GS datasets):
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-hse-shared-sbe-20261119

Current effective counts:standard670/fast727,578 same,74 different,18 standard-
only,75 fast-only. Embedded standard GS/TDDFT/EH migration remains; defaults and
old runtime receipts are unchanged. This is SBE consolidation, not completion of
a single maintained solver tree or the convergence campaign driver.

No CI dispatched or new SCF/real GS generated. Local record writes verified;
cloud synchronization unconfirmed.

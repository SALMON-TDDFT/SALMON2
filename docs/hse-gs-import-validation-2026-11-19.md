# Fast retained GS import and common seed validation

The internal fast GS context can import an existing compatible wavefunction
restart through the original initialization, return before SCF/warmup/GS-output,
and retain independent numerical/XC/HSE/MPI resources. Common input values are
captured per context. gs_export_seed returns detached orbitals, occupations and
partition/grid/input metadata through the existing standard gs_seed source.
That shared seed excludes native resource owners and is not a complete HSE model
or portable checkpoint. Python fast GS handles and ordinary embedded TDDFT seed
initialization remain unfinished; NativeGSBackend's supported runtime is unchanged.

Import-only resolves the unused no-mixing sentinel locally, disables irrelevant
SCF/pre-SCF controls, and retains mixing/buffer validation. Active states cannot
be overwritten. Closing a GS context retires request caches and owned child
communicators before clearing its allocation; repeated close is harmless. Import
probes omit energy/force/residual fields that are not established by importing
unconverged synthetic orbitals. Full convergence/energy/force certification is
not part of this test, nor are MD/DC/optimization import paths.

## Checks (Apple Silicon/GNU)

- Host-MPI fast HSE with two ranks and standalone serial fast HSE: four pairs
  HSE/PBE, HSE/HSE, PZ/PBE and TB-mBJ/Libxc-PZ. The two HSE states use different
  screening and coordinates. HSE-off MPI2 covers both non-HSE pairs.
- Finite complex orbitals, detached array/input copies, importing B preserves A,
  closing A preserves B, fresh A equals the original import byte-for-byte, closing
  B preserves A, repeated close and zero live child communicators after cleanup.
  Existing synthetic restart inputs remain unchanged. New SCF calls:0.
- Common seed transfer: analytic real/complex orbitals and empty-band/orbital
  redistribution under MPI2 pass. No GS solve or restart output.
- Existing HSE MPI2 registry current/energy/history and independent/fresh/direct
  restart controls match the preceding audited runtime. Four-rank LD/SBE/two-HSE
  mixed128 and direct8→128 resume match the prior results; zero replay.
- Normal SCF iteration/output body remains byte-identical. No new normal SCF run
  is performed under the current no-GS-collection instruction. Legacy GS probe
  compiles against the older serial build with the import extension disabled.
- Source628 test methods pass, one skipped; no additional Python module/test
  method created. Distribution/installed results are kept with this record.

Same prepared source builds HSE-on MPI/serial and HSE-off MPI. Mandatory workspace
checker source/object audits cover solver, adapters and every linked probe object;
repository/workspace checker bytes match. New IEEE inquiry tests use scalar
arguments only. Fujitsu/other platforms are not verified. Fixed/historical sources
and binaries remain unchanged; current build identities are recorded explicitly.

The new branch first failed compilation before export/shared-seed support existed.
Preserved intermediate failures include a test accessing private input fields,
legacy no-mixing sentinel rejecting an import, incompatible RT propagation in a
GS fixture, and a test checking Cartesian storage for reduced-coordinate input.
Tests use activation methods, GS-compatible fixture settings and reduced
coordinates; engine checks were not weakened to pass an invalid fixture.

Evidence:calc2026/salmonpy-hse-gs-import-20261119. Logs/results/hashes and audited
artifacts only; no GS/checkpoint payload collection or Actions. Local files are
verified, OneDrive cloud synchronization unconfirmed. Current effective inventory:
588 identical/70 differing/12 standard-only/75 fast-only (fast733). Five alternative
standard HSE modules remain outside the migration targets.

Final distribution checks: clean installed wheel passes628 methods (one skip);
installed preparation reproduces every recorded solver source hash and installed
MPI2 GS-import four-pair controls pass. HSE-off binary has no FFTW linkage.
The prior ordinary fast GS path's omitted no-mixing sentinel (-1) still needs
resolution audit with public GS lifecycle migration; the import-only path handles
its unused sentinel locally. No new ordinary SCF run has been authorized here.

# Python fast GS import handles

NativeGSBackend/NativeGSSession.load now reuse the existing retained-GS Python
API with a fast host library explicitly declaring gs_abi protocol1/import mode.
The builder probes the compiled protocol before recording that capability. Older
fast records without the declaration fail before loading native code; constructing
the backend still does not load/build a library or initialize MPI.

The fast GS module owns stable allocated contexts behind monotonically changing
opaque tokens. It supports import, info, array layout/read, close and resource
counts. Slot-table growth copies pointers rather than live owners. Matching token
queries are collective within the host subgroup. Closing/reusing slots invalidates
old tokens and does not close other owners. The same standard Python array/copy/
MPI/lifetime machinery is retained, rather than adding another Python module.

Each import parses into its own workspace, switches/restores symmetry and density-
symmetry scratch, preserves sym.dat, sets/restores GS calculation mode and follows
the original initialization. Unused RT factor-learning/pre-SCF controls are disabled
for import only; XC/screening/physical input remain. No SCF or GS restart output.

Available arrays are density, occupations and orbitals. Python returns None for
unestablished energy, residual, tolerance and convergence; native placeholder
values must not be interpreted as physical zeros. Computed eigenvalues/forces,
fast session.run, portable export and direct GS→TDDFT model construction are not
implemented. Standard-profile computed GS/export/GS→RT behavior is unchanged.
Fast load_archive is unverified and is not advertised as validated interoperability.

## Actual checks (GNU/Apple Silicon)

- MPI6, three groups of two: full/full and full/reduced GS density/occupation/
  orbital arrays match independently imported controls exactly. Distinct material
  IDs persist; closing A preserves B and fresh A matches before. Stale and mixed
  rank tokens reject without changing either live state. Live child counts return0.
- Different HSE/Libxc-PBE and TB-mBJ/Libxc-PZ inputs pass the same MPI6 checks.
  These are synthetic existing wavefunctions, not newly converged GS data.
- Existing retained_gs sample with --import-bundle/--second-bundle passes MPI2;
  JSON records material IDs and unknown metadata as null. No new sample file.
- Existing HSE registry current/energy/history/independent/fresh/direct-resume
  controls match prior runtime. LD/SBE/two-HSE MPI4 mixed128 and save8/resume120
  match previous results with zero replay.
- Same source builds HSE-on MPI, HSE-on serial and HSE-off MPI. The no-HSE binary
  has no FFTW linkage. Required source/object audits cover solver and adapters;
  repository/workspace eigensolver checkers have identical SHA. No new array IEEE
  inquiry. Fujitsu/other platform validation remains outstanding.
- Source and clean installed wheel:628 test methods each, one skipped. Existing
  NativeGS tests expanded; no new test method or Python module. Source→sdist→wheel
  passes. Installed preparation reproduces all recorded solver source hashes;
  installed MPI6 GS checks and MPI2 user sample pass.

Development evidence retains initial missing-backend admission, omitted symmetry
import, unnormalized RT-factor settings and missing GS calculation-mode failures.
The failures were corrected at the boundary, preserving the numerical solver.

Evidence:calc2026/salmonpy-python-fast-gs-20261119. Logs/results/hashes and audited
artifacts only; no GS/checkpoint datasets collected, new SCF or Actions. Local
storage verified; OneDrive synchronization not independently confirmed. Both source
profiles remain until fast SCF/export/ordinary TDDFT lifecycle migration completes.

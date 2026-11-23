# Built-in material configuration validation — 2026-11-15

Implemented explicit builtin_material_registry and material_from_catalog for
LD/SBE only. ConfiguredMaterial/ConfiguredState preserve material ID and definition
digest while delegating to existing models/states. No native kernel, numerical
update order, current sign or drive port was changed.

- Initial six tests failed because the module was absent; seven final focused
  tests pass using the audited LD/SBE handle library.
- Registration/resolution/model construction checked with native entry patched to
  fail if called: no native initialization occurs before state creation.
- LD and SBE configured states agree bit-for-bit with direct original API states
  on short synthetic drives. Separate states retain independent current/density;
  SBE snapshot restoration does not change the other state. Identity fields are
  immutable and retained.
- SBE file alteration is rejected against resolved hashes; an offline bundle
  can be relocated and the source originals deleted before model creation.
  Nondefault active_bands=1, correction_order=3 and nonlocal_correction=True survive.
- Read-only review found no important functional defect; SBE column-vector/bohr
  and ignored occupation-column contracts are explicit in paired guides.
- Full source suite: 344 tests, pass, one case-sensitive-filesystem skip, 8.377 s.
- Source distribution→wheel built locally with two jobs. Fresh venv isolated
  installed tests: 344 tests, pass, two skips (filesystem and optional plotting),
  8.424 s.
- Installed configured_ld sample ran under a socket audit guard: two synthetic
  material IDs, two independent states each, 20 steps. Input bundle integrity
  verified after results saved outside it. Installed module, example, guides and
  AI JSON equal source bytes.
- Manual sync and whitespace checks pass. No Actions and no full solver run.

Artifacts: `/tmp/salmonpy-builtin-source-tests-20261115.log`,
`/tmp/salmonpy-builtin-build-20261115.log`,
`/tmp/salmonpy-builtin-installed-tests-20261115.log`,
`/tmp/salmonpy-builtin-offline-smoke-20261115/validation.json`.

These are local Apple Silicon/macOS connection regressions, not physical material
validation or cross-platform/MPI evidence. Initial-state source/prepare recipes,
point/run IDs, MPI allocation and automatic output/restart identity propagation
remain unimplemented in this factory. DFT/TDDFT configuration adapters remain open.

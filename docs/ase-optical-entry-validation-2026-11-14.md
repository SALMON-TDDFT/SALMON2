# ASE optical entry validation — 2026-11-14

## Scope

Local macOS arm64, Python 3.14, ASE 3.29.0 / f90nml 1.5.0. This record covers
neutral electron inference, the new Atoms→GS→TDDFT entry, packaging, and fixed
Si numerical reproducibility. It does not establish optical convergence,
non-collinear/isolated-system support, force accuracy, or cross-platform native
results. Existing native kernels and propagation code were not modified.

## Tests and packaging

- Source Python suite: 155 tests passed.
- Built a source distribution, then a wheel from that source distribution using
  `python -m build --outdir NEW_EXTERNAL_DIRECTORY` (GNU Fortran 15.2.0).
- Native library hashes passed before wheel repair.
- New venv installed the wheel and binary ASE dependencies; installed-package
  suite: 155 tests passed. The test runner rejected any checkout import.
- Optional dependencies remain lazy: importing salmonpy.ase does not import ASE.
- New tests cover KY/UPF valence, invalid/fractional/missing charge rejection,
  explicit neutral-count mismatch, Calculator inference, preflight failure before
  native execution, existing GS/initial-loader composition, and AI guide sync.

## Real Si calculation

Ran the packaged example with an already installed, audited runtime:

```sh
python -I -m salmonpy.examples.ase_tddft RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY --verify-independence
```

Runtime executable SHA256:
`29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed`.
The runtime is the pinned macOS arm64 Si starter binary; no new native build was
needed for the integration test.

The example uses the existing starter potential, input, reference and their
provenance hashes: 8 Si atoms, 16³ grid, 2³ k mesh, 32 orbitals, PZ, dt=0.02 au,
20 steps. Electron count is omitted at the ASE entry and resolves to 32 from
8 × 4 pseudopotential valence electrons.

- SCF: 54 iterations; total energy -31.74653905752293 Hartree.
- Output: 21 A/current samples, charge-current convention j = -J_matter.
- Energy difference from starter reference: 0 Hartree.
- Maximum current difference from starter reference: 0 au.
- Comparison uses energy atol 1e-7 Hartree, current atol 1e-9 au / rtol 1e-5,
  not strict equality across platforms.
- An opposite drive produced a different current. Replaying the original drive
  from the same initial state matched the first response exactly on this binary.
  Only one worker was alive at a time.
- Reopening the GS after all consumers verified its file hashes remained intact.

The source run took 12.01 s inside the example. The installed-wheel run took
12.49 s inside the example and 12.68 s including the external observer/startup.
The source and wheel used different installed NumPy versions (2.4.6 and 2.5.4);
the native binary was identical.

## Resources and limits

User budget: total 4 GiB, 4 CPU cores. Runs used one native thread and one worker
at a time, with BLAS threads set to one. An external observer sampled the Python
process and all descendants using `ps` every 0.1 s. Maximum summed RSS was
263,061,504 bytes (about 250.88 MiB); shared pages may be counted more than once
and short peaks may be missed. A 3.2 GiB polling watchdog did not trigger; this
is not an OS hard memory limit.

The sample's report correctly says its own peak-memory measurement is
`not_measured`; the external observer stores a separate resource report. Do not
turn missing resource data into calculation failure or a claim of zero usage.

This remains the short Si connection regression. Its 0.4 au waveform is not a
physical multi-cycle laser example, and 16³ / 2³ is not a convergence guarantee.
Grid, k mesh and orbital count therefore remain explicit in the new ASE API.

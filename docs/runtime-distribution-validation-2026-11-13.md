# Runtime distribution validation — 2026-11-13

Implemented explicit status/install/verify/build commands and the same Python
interfaces. No download or compilation during import. Full runtimes use the
existing executable adapter, preserving solver source and A-driven boundaries.

Local macOS arm64 validation:

- 89 Python tests passed, including installer compatibility, hash failures,
  archive rejection, local CI archive, build diagnostics and existing APIs.
- Standalone sdist → wheel → repaired wheel → installed package: 89 tests passed.
- Installed-wheel preparation resources successfully prepared pinned upstream
  sources without a SALMONpy source-tree dependency.
- Public `runtime build` built/audited the pinned runtime with GNU 15, LP64
  OpenBLAS and Libxc 5.2.3. Missing LAPACK/Libxc now fails configuration rather
  than falling through to upstream dependency source downloads.
- Packager requires matching build.json provenance/hash, rejects a different
  executable architecture, bundles compiler/BLAS libraries and includes licenses.
- Resulting roughly 15 MB macOS arm64 archive was installed in a different
  directory. The ordinary TDDFT output matched the pre-relocation executable
  byte-for-byte; A→J, independent/interleaved states and recovery passed.
- The native fixture uses synthetic normalized orbitals, not a converged GS;
  this is relocation/adapter verification, not a scientific accuracy benchmark.

CI investigation:

- User-linked run 43354543306 failed installing NumPy 2.5.4 from source on an
  old manylinux environment (GCC 10.2.1; NumPy required >=10.3). It built a1,
  not this branch's a4. Latest pre-change Ubuntu job 131766181730 in run
  43389628991 passed with manylinux_2_28. Enforce binary NumPy/SciPy/Matplotlib
  dependencies in future wheel CI to prevent accidental source compilation.
- The same latest run had a Windows-only test failure: a macOS launch test
  assumed Unix permission 0644 on Windows. It now checks the actual mode is
  unchanged, retaining the mocked access checks and launch assertions.

Full runtime bundle CI targets Linux x86_64 and both macOS architectures with
PZ/TB-mBJ and no Libxc. Windows full runtime, MPI distribution and broad Libxc
functional support are not claimed. CI uploads only after relocated native tests.
Current branch assets are not upstream SALMON releases. The public catalog stays
empty until assets are actually published in the future independent repository.

## Completed branch CI

Implementation commit: `6fcfe2d5945ac90482cb3c72091803b183b55cdf`.

- [LD wheels run 43391897851](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43391897851):
  Ubuntu 22.04, Windows 2022, macOS 15 arm64, macOS 15 Intel, and sdist all passed.
  This includes installed-wheel tests with extras and a separate no-extras run.
- [Full runtime bundles run 43391897849](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43391897849):
  all three Linux/macOS profiles passed explicit builds, source/object audits,
  packaging, archive installation, band-path checks and relocated TDDFT tests
  (ordinary and predictor-corrector). The artifacts are approximately 10–14 MB.
- [TDDFT native worker run 43391897760](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43391897760):
  passed on Ubuntu.

Runtime artifacts are named `salmonpy-runtime-ubuntu-22.04`,
`salmonpy-runtime-macos-15`, and `salmonpy-runtime-macos-15-intel`. Download through
that run's Actions page; the outer ZIP contains `catalog.json` and the inner
runtime ZIP. Follow RUNTIME.md with `--catalog` and `--archive`. These are branch
artifacts with GitHub retention, not permanent release URLs. The three CI runtime
profiles do not include Libxc; the separate locally tested Libxc profile must not
be confused with them. Windows wheel success does not imply a Windows full solver.

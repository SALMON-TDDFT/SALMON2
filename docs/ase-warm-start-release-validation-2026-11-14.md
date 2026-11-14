# ASE warm start, XC forces and independent runtime distribution

Date: 2026-11-14. Local platform: macOS arm64, Python 3.14, ASE 3.29.
Budget: 4 GiB total, four cores; native trials used one thread and ran sequentially.
Runtime compilation used two jobs. No new memory enforcement or peak-memory
claim is made for the warm-start matrix.

## Native comparison

`tools/test_ase_warm_start.py` compared displaced eight-atom Si from fresh SCF
with opt-in wavefunction reuse. Cell 10.26 bohr, grid 16^3, k mesh 2^3,
32 states, displacement 0.05 bohr, SCF tolerance 1e-11. The tool checks the
written coordinates, retained source hashes, and native-to-ASE force conversion.

| XC / potential | Fresh iterations | Warm iterations | Energy difference (Ha) | Maximum force difference (a.u.) |
|---|---:|---:|---:|---:|
| PZ / KY | 62 | 43 | 4.0643e-11 | 2.7106e-11 |
| Libxc-PZ / KY | 63 | 52 | 2.8013e-11 | 1.7304e-11 |
| Libxc-PBE / NC-UPF | 71 | 56 | 7.3541e-13 | 3.7243e-11 |

Acceptance: energy absolute difference <=1e-7 Ha, force <=2e-6 a.u.
All passed; source state remained unchanged, and native force conversion was
exact. These are numerical connection tests, not grid-converged forces or a
validation of stress, molecular dynamics, isolated or noncollinear calculations.

The Libxc 5.2.3 runtime was explicitly prepared outside the checkout from the
pinned source and compatibility patch. Its executable SHA256 is
`d6afce46cc40db3e36a7e7b33e9901b9afbf63148a6f9cb8248b12973332c3fb`.
KY Si SHA256: `8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c`.
PBE NC-UPF SHA256: `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282`;
source/provenance is retained in the earlier Libxc validation record. This
external UPF and the local Libxc runtime are not included in the preview Release.

## Release provenance

The private [runtime-2026.11.14 Release](https://github.com/otobetomohito/SALMONpy/releases/tag/runtime-2026.11.14)
was published with nine assets (three unchanged outer CI ZIPs, three inner
runtime ZIPs, catalog, provenance and SHA256SUMS). All three jobs in
[own-repository CI 43554502676](https://github.com/otobetomohito/SALMONpy/actions/runs/43554502676)
passed at commit `8006261b9e48541b591e496a813b3f06ca5fc1b1`.
This CI predates the Python warm-start changes; it establishes the runtime
build/relocation/MPI regressions, not the new warm-start test across platforms.
The workflow now includes a separate PZ warm-start test for subsequent runs.

The downloaded outer ZIP hashes matched the CI-displayed digests. Installing
the Mac inner ZIP through the combined release catalog and running the ASE
Si-current example passed all reference checks and exact independent-state
replay. Its executable hash was
`29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed`.
This is the same binary hash as the original reference runtime; energy/current
differences were zero on this machine, but comparisons still use tolerances.

The original three Si artifact entries remain accepted for existing prepared
examples. Current guidance points to the new private Release, requiring an
authorized browser download and local `--archive` installation. Source acquisition
still points to SALMON upstream; no full solver sources were restored here.

## Package validation

Source suite: 158 tests passed. The source distribution was built into a wheel;
installed native library hashes passed. Installed suite: 158 tests passed. Running the installed wheel with `python -I`
prepared the newly pinned outer Release ZIP and completed Si current calculation
in 13.04 seconds at one thread. Calculation and reference validation passed;
energy/current differences were zero. Peak memory was not measured and was
recorded as null, not a calculation failure. The 4 GiB budget was recorded,
not enforced by the operating system.

Validated wheel SHA256: `93d2022219447d28ab63669a94a71a38eb137ec73929b963b018e79ce0953f51`.

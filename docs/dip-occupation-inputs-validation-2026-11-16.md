# TDDFT electron-count diagnostic input validation — 2026-11-16

`subdip` now obtains nelec/nelec_spin from get_occupation_input(system). The reduction, dipole outputs, strict lower/upper inequalities, spin-count precedence and legacy STOP are unchanged. No public API or relaxed outer TDDFT guard. Remaining I/O, ionic inputs and termination behavior are outside this change. The broad 0.1x–10x diagnostic does not certify precision conservation or physical convergence.

Local Apple Silicon/GNU/OpenBLAS/Libxc5.2.3 serial CPU checks:

- Final diagnostic probe fails to complete against the preceding HSE build after globals are poisoned (density is 4, legacy count 1000). The same probe passes against the changed build.
- Four original-subdip comparisons pass exactly for electron integral and electronic dipole histories: scalar count, two-spin positive sum, two-spin fallback to nelec, and uninitialized legacy fallback.
- Eight strict-bound cases (lower/upper for all four configurations) stop through the original diagnostic, including STOP's zero process exit code. The runner checks the diagnostic text and completion marker, rather than trusting exit status alone.
- Existing occupation solver probe still passes after its runner is generalized with --dip.
- Full TDDFT regression passes eight configurations: k2/k3, exact/learned ACE Taylor, direct Taylor, learned PT-CN; OMP2, distinct fixed dt/step capacities/OpenMP routes. Interleaved/reverse/recreated traces and currents/energies match solo output exactly; native original-driver comparison uses rtol1e-12/atol1e-14.
- Source-to-wheel build succeeds; source and installed Python suites each run435 tests,1skip, no failures. Installed wheel's own --dip native probe also passes, confirming the packaged Fortran fixture.

These are synthetic regression fixtures, not physical convergence, new spin-polarized full-TDDFT support, MPI validation or cross-platform evidence. No Actions or push. Persistent output root: `/Users/otobetoshihito/.codex/outputs/salmonpy-state-inputs-20261116` (local evidence only; not an installation requirement).

Evidence hashes:
- `dip-green/result.json`: `c8d42c7c5fb53cb692a9d99c08c079ae31d7afba67a0e555fde8a63d82f39585`
- `dip-occupation-regression/result.json`: `22779365582056d475acdf7e759a41d7783a69f701e43df1c85e5f5bbd2a103f`
- `dip-tddft/result.json`: `29eb79d1ddcf2a3c4c4db853e4f0175e4e56b541f542d0800cc62feb53defd32`
- `dip-installed-probe/result.json`: `5664d10638f2191950688c6ac9ab8215243f56bbbd98f5e1adaea29f8e9a0353`
- `dip-build.json`: `42780172ddc66b16891081fc2e9228a553b9e4f4ef20077fdc9266acf0241dbd`
- `dip-source-tests.log`: `efd15eeeed7e6ddce529cca444350a1583cd135e71049c34dfbb36ac15b57584`
- `dip-installed-tests.log`: `107a2d6c1d5dab95d9d9f22f88d206b8ba8eb8a15aea421a5dbafe4bb0ba8aba`
- `dip-dist.log`: `6124f1f0478db11bc6cb02da877f978955be731a804b251e060b5dcffd2560c7`

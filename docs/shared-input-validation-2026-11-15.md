# Shared-input extraction validation (2026-11-15)

Local Apple Silicon/GNU Fortran/OpenBLAS/Libxc/FFTW validation. No GitHub Actions or network builds were launched. Build concurrency 2; numerical probes use one thread, full-context interleaving uses two. All fixtures/build/output directories are external; no collected production GS data.

## Change and limits

Value configs now own k generation/symmetry, radial pseudopotential loading and XC identity/parameters. HSE consumes owned mesh/symmetry and functional selection/fraction/screening/radius; XC and total-energy selection agree. Compatibility adapters capture old inputs rather than swapping active global configurations. Scalar numerical statements and update order are retained.

Full native TDDFT still guards common grid/k/XC and other physical inputs. Grid PP construction/atom mapping, SCF phase, algorithm controls, communicators, fixed file units, and some unsupported paths retain globals. No heterogeneous full-TDDFT, concurrent-call, MPI-isolation, general atomic-failure or physical-convergence claim follows. Explicit lower-level input assumes validated configuration. Resource-bearing Libxc/HSE objects are not copied.

## Tests

- K config: seven modes (shift, even/odd gamma centering, both unfolding modes, normalized file weights, symmetry), alternation and poisoned globals, compared with frozen attributed original generation and symmetry code.
- PP config: scalar unmasked licensed Si KY input, legacy capture vs explicit poisoned inputs, alternate local channel/mass, copied file config and A/B/A arrays. This is loader-ownership evidence. Masked, UPF, SOC and generated atomic-orbital branches were not exercised by this probe.
- XC config: alternating HSE06/PBEH40 semilocal evaluations, poisoned global XC/omega/phase, explicit pre-SCF, distinct omega, invalid-input preservation. It does not validate all Libxc functionals.
- Native optimized HSE: independent FFT/ACE owners with distinct cell/screening, global mesh/symmetry/XC/radius poisoning, isolated process comparison, close/recreate.
- Full contexts: k2/k3 × Taylor ACE off/learned, full Taylor off, PT-CN learned (8 cases), differing fixed dt, OMP layout and capacity, forward/reverse interleave and recreation. Raw state traces compared against separate processes; native Acos2 compared to original fast PR.
- Original full solver: k2/k3/k6/k8 × dt 0.001/0.002 × four propagation/factor modes (32 cases). Current and energy output comparison against immutable fast PR c9e8ce431f43eefd4a88998134adb4796f8ec837.
- Python suite and source distribution → wheel → fresh environment tests; installed runtime-source preparation and byte comparison. Source and all object audits run before native link.

## Final results

All checks passed on the final source:

| Check | Result |
| --- | --- |
| k/PP/XC probes | All passed |
| Native optimized HSE owner/input independence | Passed, isolated/interleaved results identical |
| Full TDDFT contexts | 8 cases, raw traces bit-identical to separate processes |
| Original fast-PR driver | 32 cases, current/energy maximum output difference 0 |
| Source Python suite | 305 passed |
| sdist → wheel | Built successfully |
| Fresh wheel environment | 305 tests run, 1 optional plotting test skipped, no failures |
| Installed runtime preparation | All 23 replacement files byte-identical to tested source |
| Source/all-object eigensolver audits | Passed |
| Changed executable Fortran line lengths | At most 132 columns |
| AI manual/guide synchronization | Passed |

Final executable SHA256: `e2b7c2ae1cbd18637b92f07a7e84701a88463c9e1b8117c98bf9f8b8821edbdf`.
Replacement manifest SHA256: `29c6b5fcc8b9269f11c5c9a906219135e9fb4d35d8a0ec3f8c6dd260978c94da`.

External local receipts: `/tmp/salmonpy-shared-input-release-build-20261115.json` and
`/tmp/salmonpy-shared-input-release-{kprobe,ppprobe,xcprobe,native,contexts,default}-20261115/result.json`.
Package: `/tmp/salmonpy-shared-input-release-dist-20261115`; installed source:
`/tmp/salmonpy-shared-input-installed-source-20261115`.
These are machine-local validation records, not installation requirements or published artifacts.

Numerical fixtures are synthetic connection/regression checks, not converged material benchmarks.

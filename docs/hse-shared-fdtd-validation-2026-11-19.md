# Fast-PR source: FDTD, LD and native E/H validation

The explicit-routing preparation now references the six existing standard
Maxwell files (CMake, Lorentz–Drude, E/H kernel, E/H driver, batch entry and
embedded E/H). These are complete-file inventory references, not new source
copies or patches. The prepared LD kernel must equal the Python/native kernel.
Historical snapshots are unchanged. The earlier SBE sharing remains in place.

The fast host library declares `field_abi` only after probing the linked E/H
checkpoint protocol. NativeEHField accepts that declaration as well as the
legacy embedded record. The batch executable has its own executable-build.json
with its own SHA256; library and executable fingerprints are never exchanged.
The same prepared source builds with HSE enabled and disabled. The disabled
executable has no FFTW dynamic dependency. Source and all application objects,
including host adapters, passed the required eigensolver audits.

## Numerical checks on Apple Silicon, 2026-11-19

- Interface matrix: 38 cases, 76 observation files byte-identical to the audited
  standard executable. Includes single LD/dielectric, LD/LD and LD/dielectric
  interfaces, magnetic dielectric, three axes, ID relabelling and shape order.
- MPI standalone FDTD: 160 steps, 48×40×40, three LD materials, one versus four
  ranks and two spatial layouts. HSE ON/OFF both bit-identical. Magnetic
  dielectric also matches across 1×2×2, 2×1×2 and 2×2×1 decompositions.
- Native E/H: two ranks, 32 steps, two local LD materials. Mixed-entry rejection,
  real integrated-A overflow and capacity rejection preserve state; retry
  traces A/E/H equal the standard embedded library exactly.
- Ordinary mixed sample: four ranks (one field, one SBE, two shared-model HSE),
  two different full-k HSE inputs (Taylor and PT-CN learned), 128 steps. Using
  the common fast library for field and material gives exactly the standard
  field control. Direct checkpoint at step 8 resumes with 120 new steps and
  zero replay, matching uninterrupted execution.
- Protected restore: five ranks (two field, one SBE, two HSE), two local LD
  materials, failure after save at step 8. All original state bytes survive;
  three retry future steps match; zero replay. Memory advice uses 4 GiB.
- Source and clean installed-wheel suites: 628 tests each, one skipped.
  Installed-wheel MPI E/H recovery also passes and matches the standard traces.
  Preparing sources from the wheel reproduces all 255 recorded source hashes.

These are regression/connection tests with existing synthetic electronic
fixtures, not converged-material validation. No new SCF or GS collection was
performed. Arbitrary 3D scattered fields need full-k inputs when they break a
retained reduced-k little group; the symmetry rejection was preserved.

An initial source test lacked the native kernel path, and an initial 3D mixed
case correctly rejected a symmetry-incompatible reduced-k input. Both logs are
retained alongside the successful runs. A restricted local MPI socket launch
also failed before calculation; the permitted local launch then passed.

Records, logs and distribution artifacts are stored under the canonical
calc2026/salmonpy-hse-shared-fdtd-20261119 directory. No GS/restart payloads are
collected. GitHub Actions was not used.

## Remaining source consolidation

Effective fast inventory: 584 identical, 71 differing, 15 standard-only,
75 fast-only paths. The five alternative standard HSE modules must not be
ported. Ordinary GS/TDDFT initialization, parser/model ownership, output and
remaining build integration still need consolidation. Do not remove the
standard source tree until those paths have migrated and passed regressions.
The SALMON-specific grid/dt convergence driver remains a separate next task.

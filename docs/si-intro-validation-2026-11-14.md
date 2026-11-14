# Si introductory wheel workflow validation

2026-11-14, macOS 26.7.2 arm64, Python 3.14.5, NumPy 2.5.4.
User budget: 4 GiB total and 4 CPU cores. Native calculations used one thread,
sequential DFT then TDDFT, with numerical-library thread counts limited to one.

## What was supplied and executed

The updated wheel includes the Si example, fixed input, upstream Si KY
pseudopotential, expected endpoint charge currents, provenance and license.
The full solver was downloaded from successful Actions run 43541732653,
artifact 14092143545 (`salmonpy-runtime-macos-15`). No local SALMON source/build
was used for either installed-wheel trial. The package itself was explicitly
built by the maintainer from sdist, then compiler runtime libraries were bundled.

- Runtime outer ZIP SHA256:
  `ec98ba8df91dc26fd4389232dc66273d649da0977e2cb76bf22b30f81da4888c`
- Installed solver SHA256:
  `29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed`
- Tested repaired wheel SHA256:
  `c0716130b458cc0ab60ac999beb35a8e5455bc506936aac5bb45acded24f96bf`
- Si pseudopotential SHA256:
  `8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c`

Two fresh venv trials passed: ordinary binary-only dependency installation and
an offline installation using only the supplied SALMONpy/NumPy wheels. Both used
isolated Python (`-I`), verified the package path was inside the venv, diagnosed
LD/SBE libraries as loadable, installed the pinned prebuilt solver, and executed
GS and TDDFT. The offline kit requires Apple Silicon, macOS 26+ and Python 3.14;
other-platform runtime availability is not an assertion that this wheel works
there. No public release or PyPI publication was performed.

Reproduction from a source checkout, using NEW_OUTPUT outside the checkout:

```sh
python tools/test_wheel_si_intro.py UPDATED_WHEEL ACTIONS_OUTER_ZIP NEW_OUTPUT \
  --memory-budget-gib 4 --threads 1 --wheelhouse DIRECTORY_OF_BINARY_WHEELS
```

The installed example itself requires no checkout. See `salmonpy/SI_INTRO.md`
and `.en.md` for preparation/run commands. It refuses different runtime archive
hashes and modified input/reference files, and never overwrites previous results.

## Numerical result

Eight Si atoms, 10.26 bohr cubic cell, 32 electrons/states, 16³ real-space grid,
2³ k points, PZ; SCF residual threshold 1e-10. SCF converged in 54 updates.
TDDFT used 20 middlepoint steps, dt=0.02 a.u., the fixed short z-polarized pulse.

- GS energy: -31.74653905752293 Hartree.
- Maximum absolute current component: 1.4727524857430097e-6 a.u.
- Time samples: 21, including the physical initial current.
- Both wheel trials: zero energy and current differences from the fixed reference.
- Acceptance tolerances: energy absolute 1e-7 Hartree, current atol1e-9/rtol1e-5.
- Quantity: physical charge-current density j = -SALMON electronic matter J.

This is a real SCF-converged starter, not synthetic orbitals. Spatial/k-point/time
convergence of optical observables is NOT certified. A one-step shift and current
sign reversal are rejected by the comparison tests; automatic comparison is not
an independent validation of SALMON's underlying physical model.

## Resources and failure classification

An external observer sampled summed parent/descendant RSS approximately every
0.1 seconds, with a 3.2 GiB sampled-stop threshold to leave headroom under the
4 GiB budget. Shared pages can be counted twice, transient unsampled peaks can
be missed, and this is not an OS allocation limit.

| Trial | Install + workflow wall time | Largest observed summed RSS |
|---|---:|---:|
| Binary-only online install | 15.43 s | 259293184 bytes (247.28 MiB) |
| Offline kit install | 16.19 s | 259309568 bytes (247.30 MiB) |

No observer errors or budget stops occurred in these trials. The sample's own
report honestly leaves peak memory `not_measured`; the external observer report
is the separate evidence. Initial reference generation had failed `ps` access
and therefore no memory observation; its successful calculation is not
misreported as a failed calculation or as measured budget compliance. The
subsequent fresh-wheel trials used working observation, with separate outputs.

## Checks and retained evidence

- Source suite: 149 tests, 6 optional native cases skipped, no failures.
- Installed minimal-wheel suite: 149 tests, 7 ASE/plot cases skipped, no failures.
- Raw/repaired native library hashes, sdist-to-wheel build, guide/JSON consistency.
- Changed assets/unknown runtime/ZIP member traversal/sign/time/tolerance errors,
  native failure, measurement failure and preserved successful calculation on
  reference mismatch covered by unit tests.

Local evidence (temporary, not paths to assume on another machine):
`/tmp/si-intro-wheel-trial`, `/tmp/si-intro-offline-trial`, their adjacent
`.resources.json` files, and `/tmp/salmonpy-si-intro-macos-arm64-py314.zip`.
The kit includes `validation.json` and checksums; it contains no saved GS state
from the developer, so users recompute their own GS from the bundled inputs.

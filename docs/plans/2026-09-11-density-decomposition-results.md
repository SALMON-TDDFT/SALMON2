# Frozen-potential density comparison checkpoint

## Scope

The opt-in `SALMON_DG_DENSITY_DIAGNOSTIC_PREFIX` export performs one conventional
DC-LCFO solve after constructing Vlocal from the immutable DC density, before
the local DG stage. It reconstructs density directly from LCFO's orthonormal
core basis and coefficients. It does not smooth or stitch buffered WF tails.
The existing retained-state upper-edge occupation check remains active; no
material-specific state count is introduced.

The first DG density is captured before total-energy evaluation can update the
potential. Exact local-potential equality with the conventional solve is checked.
The last terminal density is exported separately, whether converged or not.
No additional conventional solve occurs when the environment variable is absent.

Per-rank text exports contain physical point IDs, integration weights, DC,
conventional LCFO, first DG and final DG densities, and frozen Vlocal. Headers
include MPI size, rank/fragment identity, global grid size, DC publication,
ownership and immutable fingerprints, electron count and temperature. Existing
files are never overwritten. These are diagnostic exports, NOT restart caches.

## Observed Si8 run

Directory: `/tmp/si8-density-decomposition-20260911`

Command (run in that directory):

```sh
OMP_NUM_THREADS=1 OMPI_MCA_rmaps_base_oversubscribe=1 \
SALMON_DG_DENSITY_DIAGNOSTIC_PREFIX=/tmp/si8-density-decomposition-20260911/density \
mpirun -np 8 /tmp/salmon-v230-dg-integration-20260907/build/salmon < inputfile > run.log 2>&1
```

Exit 0; all eight ranks emitted `end SALMON`. The first sandbox attempt failed
before MPI startup and is retained in `mpi-sandbox-denied.log`.
Run-log SHA256: `8571384fadabda0d1f4c79578bd338c3fdb85da76543bf2a70312ca62804176c`.

- DC seed read-only hit: 6737953370837868388; SCF skipped.
- WF cache read-only hit: 4254267444278987270; no Wannier regeneration.
- MPI ownership fingerprint: 398352827599223578.
- Temperature: 300 K (0.0009500434690366825 Hartree).
- One conventional LCFO solve; four terminal DG solves, still unconverged.

All density norms below share the denominator `||rho_DC||`:

| Diagnostic | Relative weighted L2 |
|---|---:|
| Conventional LCFO minus DC | 0.13268802131279162 |
| First DG LCFO minus conventional LCFO | 0.875022490960162 |
| First DG LCFO minus DC | 0.9469827106345656 |
| Last DG LCFO minus first DG LCFO | 0.051158279543262145 |

Integrated electron counts: DC 32.0000000000481, conventional
31.999999999999986, first DG 32.000000000000014, final DG 32.0.

The prior terminal fixed-point residual 0.766352494 has a different denominator
and is not any of the above accuracy measures. Scalar norms are not subtracted.
The normalized cross term is 0.113505783554759 and satisfies the squared-norm
decomposition identity.

This is evidence of a substantial DG-route density difference relative to
conventional LCFO, not proof of a specific bug or absolute DC accuracy. Basis
truncation and DG operator effects remain combined. Energy-error separation has
not yet been implemented.

## Validation and remaining work

The route assertion and snapshot-reader tests failed for missing functionality
before implementation and passed afterward. The reader rejects incomplete
rank sets, duplicate points, wrong rank/fragment mapping, inconsistent source
fingerprints, nonfinite values and truncated files. It emits SHA256 hashes of
the analyzed bytes for evidence tracking; these are not an authenticated cache
publication. The existing terminal/local route test also passes.

Read-only analysis:

```sh
python3 -B tests/dg/density_error_decomposition.py /tmp/si8-density-decomposition-20260911/density
```

Reference caching (to avoid the opt-in conventional solve on subsequent runs),
energy diagnostics, and production accuracy-gate wiring remain pending. Until
then, use this preserved export for analysis without rerunning DC or Wannier90.
The run predates only the follow-up collective environment-prefix validation
and equivalent formatted-output repeat-count cleanup; the final source is
rebuilt separately. No full regression certification is claimed.
Final rebuild exited 0 (`build-final.log`). The intermediate missing
`MPI_CHARACTER` import in the added prefix-consistency check was resolved by
using SALMON's existing `comm_bcast` character overload. Both Python unit cases,
the density snapshot route assertion, the existing terminal/local route
assertion, and `git diff --check` passed after that correction.

## Fourier diagnosis of the preserved snapshots

No additional SALMON, DC, SCF or Wannier calculation was run for this analysis.
The input gives a cubic cell L=10.26 Bohr with a 24x24x24 grid and 2x2x2
fragments. Physical point IDs were assembled with x-fastest ordering. The
analyzer verifies complete coverage, uniform volume weights, and the snapshot
hashes before using the densities. NumPy 2.4.6 was used.

```sh
python3 -B tests/dg/analyze_density_fourier.py /tmp/si8-density-decomposition-20260911/density --grid 24 --length 10.26
```

Detailed output is preserved at
`/tmp/si8-density-decomposition-20260911/fourier-analysis.json`.

| Wavelength band | Conventional LCFO minus DC power | DG minus conventional LCFO power |
|---|---:|---:|
| Long: 5.13 to 10.26 Bohr | 41.40% | 72.75% |
| Middle: 2.565 to below 5.13 Bohr | 33.64% | 22.06% |
| Short: below 2.565 Bohr | 24.96% | 5.20% |

Percentages are fractions of each difference's squared L2 norm, NOT relative
density errors and NOT an error budget for the orbital basis. The zero mode is
negligible, consistent with equal electron counts. These band boundaries are
diagnostic choices based on the fragment width, not universal definitions.

The DG increment is strongest at integer reciprocal-index squared magnitude
4 (wavelength 5.13 Bohr, 54.61% of power), followed by magnitude 3 (5.9236 Bohr,
16.78%). The very longest shells 1 and 2 together carry only 1.36%. Thus the
dominant discrepancy is at approximately the fragment length scale, not mainly
the longest whole-cell scale and not mainly fine-grid oscillations.

The DG-difference Parseval relative defect is 2.13e-16. An analytic test combining
unit-amplitude mode 1 and amplitude-2 mode 5 verifies the 20%/80% long/short
power split, total mean square 2.5, and constant/zero special cases. It failed
before the analyzer existed, then passed. Existing density-reader tests and
`git diff --check` also pass.

This supports investigating fragment-scale smooth-density reproduction before
attributing the discrepancy to short-wavelength PW shortage. It does not yet
distinguish basis representability from DG operator/interface effects. Density
is quadratic in orbitals, so its Fourier wave numbers cannot be mapped directly
to an orbital PW cutoff or used alone to prove a missing basis mode. Next check:
project the low-wave-number orbital probes onto the retained WF+PW space and
measure residuals, independently of the Hamiltonian solve.

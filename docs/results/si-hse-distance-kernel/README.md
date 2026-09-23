# Common distance-kernel exchange and ACE

The common electron-separation cutoff passes the Hermiticity and variational
checks that the earlier orbital-specific support cutoff failed. A new blocked
density-matrix implementation also reproduces full-support HSE without any
cutoff. This is an optional backend of the independent Python HSE reference,
not native SALMON hybrid dynamics.

## Method

For the periodic supercell, take the inverse FFT of the existing sampled HSE
multiplier and optionally set entries with minimum-image electron separation
`|r-r'| > Rc` to zero. This retains the original discrete kernel normalization
and aliasing; it is not an analytic isolated-space erfc integral. No individual
WF is clipped, and no MLWF optimization is necessary for this backend.

For each primitive-cell row block, form the occupied density matrix at each k,
including the phases for the shifted mesh. Transform over the k mesh to cell
translations, multiply by the common real-space kernel, transform back, and
apply to the target orbitals. Targets may differ from sources. The real,
inversion-even kernel gives a common Hermitian exchange action, invariant under
unitary mixing of the occupied states at each k. The same action enters the HSE
energy and its derivative.

Primitive columns that cannot intersect the radius are skipped. Translation
entries outside the radius are zero, but the k FFTs still process their dense
blocks. Thus this implementation does not yet exploit all possible sparsity.

## Fixed-state Si benchmark

Si8, 32 electrons, 16 occupied spatial orbitals, shifted4³ k mesh, 12³ primitive
grid, spacing0.855 bohr; HSE06 mixing0.25, screening0.11 bohr^-1. Input is the
converged original full-support HSE ground state. No k convergence study.

| Method / Rc | Exchange action (s) | Relative action error | HSE exchange error (meV/atom) | Min kernel eigenvalue |
|---|---:|---:|---:|---:|
| Original full MLWF | 12.212 | reference | reference | positive |
| New blocked, no cutoff | 7.692 | 7.55e-15 | 0 | 0.310 |
| 6 bohr | 7.138 | 1.316% | 32.930 | -1.630 |
| 10 bohr | 7.213 | 0.1055% | 1.329 | 0.102 |
| 16 bohr | 7.222 | 0.003253% | 0.01602 | 0.317 |

One warm-up and two timed applications; one BLAS/FFTW thread. The blocked
method uses NumPy FFT. Original MLWF timing excludes localization and transforms;
blocked timing includes density construction, transforms and action. The earlier
full-support RT job overlapped part of the benchmark, so modest timing differences
between radii should not be interpreted as robust cutoff speedups.

The measured ~1.6x exchange speedup is predominantly an algorithm change, not a
benefit uniquely attributable to cutoff. For Rc10 and16 every primitive column
is still selected. Rc6 reduces the processed primitive-pair fraction only to
0.937. A radius16 kernel has ~25% nonzero translation entries, but the FFT costs
remain dense. The explicit selected-array workspace estimate is ~255 MB; this
is not peak RSS and excludes internal FFT/BLAS temporaries and some persistent
and indexing arrays.

## ACE and cutoff acceptance

Every snapshot in the table has occupied-metric anti-Hermiticity around1e-15 and
ACE interpolation error around2e-15. However, the Rc6 kernel has negative Fourier
eigenvalues: occupied-space ACE success alone is not sufficient for general
use. `HSEFunctional` calls the backend's positivity validation and rejects such
a kernel before integration. Diagnostic `DistanceExchange.apply` remains
available to measure rejected cutoffs. Rc10 and16 pass positivity on this grid;
the constructor validation applies again on any other grid.

For Rc16, ACE compression is0.00761s and median application0.00376s. With local
rebuild cost L, compression B and application a, `L+B+m*a < m*L` for
`m > (L+B)/(L-a)`. Here that boundary is1.0016, so two uses already amortize
compression at this fixed state. This is not a full trajectory speedup estimate.

## Short PT-CN pilot

Optional `exchange_backend` is supported by `HSEFunctional` and
`benchmark_ptcn.run`. The existing full-MLWF default is unchanged. A backend
independent of the MLWF gauge skips localization. Pair screening cannot be
combined silently with this backend. ACE and endpoint residuals use the same
modified kernel; they are not residual checks against the unmodified kernel.

Rc16 was tested for one step dt0.32 a.u. (0.00774043 fs), at zero field and
constant A_z=1e-4. Both passed the unchanged full-backend residual tolerance1e-10,
electron-number and orthogonality gates.

| Pilot | Step time (s) | Endpoint residual | Energy change (Ha) |
|---|---:|---:|---:|
| Zero | 16.938 | 4.21e-11 | 0 at printed precision |
| Impulse | 25.688 | 9.47e-11 | -7.11e-14 |

The impulse used three exchange rebuilds and39 inner applications. The previous
full-MLWF PT-CN pilot with identical dt took37.823s: about1.47x faster propagation,
or49.267s versus33.281s including bootstrap (1.48x). These are separate runs and
the new run changes the kernel slightly; not a matched-error production claim.

At the same final time, impulse density differs from the prior full-MLWF PT-CN
state by4.41e-9 relatively, orbitals by7.64e-7, and current by9.51e-11 relatively.
Current is compared with the matching applied A, including nonlocal terms. Such
an early-time current agreement is not a bound on later spectral peak errors.

The starting GS was optimized for full HSE, not Rc16. Its zero-field PT residual
with the cutoff is9.63e-6 Ha (versus the original GS residual7.93e-7 Ha), so this
switch introduces a small quench. A ground state consistent with the chosen
cutoff is required before interpreting long-time cutoff spectra. These pilots
are local integration checks, not an optical spectrum or long-time validation.

## Scaling and limits

With fixed primitive grid G, occupied count b and row block B, density formation
and application cost O(Nk G² b), and the k transforms O(Nk log Nk G²). Working
arrays scale O(Nk B G); no full-supercell dense exchange matrix is stored. This
changes the Nk² dependence of the prior pair-exchange implementation even without
a cutoff. It can be expensive for large primitive cells because of G²; this is
not a universal linear-scaling algorithm in atom count. Sparse local blocks are
the next opportunity if further cutoff speedup is needed. No larger-k performance
measurement was done, and physical k convergence remains outside the scope.

## Verification and other calculation status

67 Python tests pass, including energy finite differences, arbitrary-target
Hermiticity/linearity, occupied-unitary covariance, full and finite-cutoff direct
parity, onsite limit, spatial skipping, shifted/permuted mesh3, and rejection of
an indefinite kernel. Independent review found no blocker; its additional random
comparisons agreed to2.7e-16–1.1e-15.

The separate original full-support impulse run completed130 steps to1.00625588fs,
with energy change -6.27e-12Ha and maximum electron-number deviation8.64e-11.
It was not restarted or switched to the new backend. Its completion metadata is
included in `validation.json`; it is not a completed exciton-spectrum analysis.

## Reproduction

```
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/bin/python3 \
 samples/hse_mlwf_reference/benchmark_distance_exchange.py \
 calculations/si_hse_reference/export calculations/si_hse_reference/scf/state.npz \
 /private/tmp/distance-kernel.json --radii 6 10 16 --repeats 2
```

For pilots, construct `DistanceExchange(model.shape, model.h, model.k, radius=16.)`
and pass it as `exchange_backend` to `benchmark_ptcn.run` with steps1, dt0.32 and
amplitude0 or1e-4. Raw measurements are `ground_state.json`, `zero_pilot.json`,
`plus_pilot.json`, and `validation.json`. The pilot reports correct an inherited
inactive FFTW-thread label; no numerical outputs or timings were changed.

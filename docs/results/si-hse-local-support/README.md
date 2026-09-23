# Fixed local-support exchange: accuracy and ACE gate

This experiment implements spatially local pair convolutions, skips disjoint
supports, and reuses reciprocal pairs. It is separate from the running
full-support HSE/PT-CN reference. It does not enable localized HSE propagation.

## Definition and validity

For fixed real orbital projectors `P_i`, define `v_i = P_i w_i` and use the
closed-shell exchange energy `E_x[v]`. The returned gradient is
`G_i = P_i K[v] v_i`, including all periodic source translations. Its real
directional derivative is `dE_x = 4 Re sum_i <dw_i,G_i>`. The gauge and supports
are held fixed in this identity. No tail renormalization is performed.

Pair density and both gradient contributions are restricted to the overlap of
the translated source and target supports. The convolution preserves the
original sampled periodic screened kernel. The local FFT width is the smaller
of twice the support width and the supercell width; the latter avoids padding
beyond a full periodic FFT. FFTW plans are cached for an operator instance.

The orbital-specific projectors break occupied-unitary invariance. This gradient
is therefore not necessarily the action of a common Hermitian operator on the
occupied orbitals. That distinction matters for ACE and for propagation with
orthonormal orbitals. A small scalar energy error does not remove it.

## HSE ground-state results

Si8, 16 occupied spatial orbitals, 4³ k mesh, 12³ primitive grid, spacing
0.855 bohr, HSE06 mixing 0.25 and screening 0.11 bohr^-1. Input is the converged
full-support HSE state. This is not the earlier PZ local-box benchmark.

Each row uses one warm-up followed by two timed applications, with single-thread
FFTW and BLAS. Setup is excluded and recorded in the JSON. A separate full-support
RT job was running concurrently; timings are indicative, not scaling results or
matched-accuracy production speedups. The total exchange action norm is used for
the relative error, not the much smaller impulse-induced current difference.

| Support width (grid points) | Side (bohr) | Action time (s) | Relative action error | HSE exchange error (meV/atom) | Max metric anti-Hermiticity |
|---|---:|---:|---:|---:|---:|
| 12 | 10.26 | 0.130 | 9.62% | 164.47 | 6.23% |
| 16 | 13.68 | 0.975 | 6.05% | 67.11 | 4.04% |
| 20 | 17.10 | 3.548 | 3.89% | 28.56 | 2.10% |
| 24 | 20.52 | 9.522 | 2.66% | 14.41 | 1.74% |
| 32 | 27.36 | 13.657 | 1.39% | 3.64 | 0.78% |
| 48 (new implementation, no cutoff) | 41.04 | 24.287 | 1.75e-16 | 7.6e-13 | 1.85e-15 |
| Full reference | 41.04 | 11.75–12.35 | 0 | 0 | ~2e-15 |

Metric defect is `max_k ||M_k-M_k†||_F / ||M_k||_F`, where
`M_nm=<u_n,G_m>`. Every truncated row fails the unchanged ACE constructor's
Hermiticity gate (1e-10). The full reference passes, with ACE interpolation error
about 2e-15. No occupied-metric symmetrization was applied to conceal the failure.

The smallest support reduces reciprocal convolutions from 8,256 to 688, and the
FFT grid from 48³ to 24³, explaining its speed. At width32 all 8,256 pairs remain
and the FFT is again 48³, explaining the absence of a speed benefit.
The width48 control recovers the original operator to roundoff. It is slower
than the specialized full-grid implementation because local gather/scatter and
masking overhead remains; it is an equivalence test, not a recommended full-grid
execution path. The last two rows express errors as fractions, not percentages.

## Decision

Do not use this orbital-specific hard-cutoff functional as a drop-in exchange
for the existing ACE/PT-CN driver. The fixed-support energy-gradient test passes,
but the common-Hermitian-operator gate does not. Adaptive support updates would
add another source of energy/gauge changes and are consequently not enabled.
The current full-support impulse run is unchanged.

This is a rejection of this particular truncation prescription, not of spatial
locality. A subsequent design needs to preserve a common Hermitian exchange
operator and its energy derivative, or explicitly formulate and validate a
different orbital-dependent evolution. Merely loosening ACE's check would not
solve that problem.

## Verification and reproduction

Tests include fixed-support energy finite differences, wrapped supports, full
support, disjoint pairs, FFTW/NumPy parity, and mesh3 inverse translations.
Independent review compared random complex mesh3 cases with direct periodic
convolution and found relative differences 1.6e-16–4.7e-16. The complete Python
suite passes 58 tests. The benchmark hashes the same archive bytes that it loads,
so atomically replaced live checkpoints cannot mismatch the recorded state hash.

```
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/bin/python3 \
 samples/hse_mlwf_reference/benchmark_local_support.py \
 calculations/si_hse_reference/export calculations/si_hse_reference/scf/state.npz \
 /private/tmp/local-support.json --widths 12 16 20 24 --repeats 2
```

Raw results: `ground_state.json` (widths12–24, recorded before snapshot hash
handling was revised; its legacy fingerprint includes the filename) and
`large_support.json` (width32 and full-support control, archive SHA256).

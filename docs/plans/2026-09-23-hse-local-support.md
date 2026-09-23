# Variational local support implementation plan

> Execute inline using executing-plans and test-driven-development.

**Goal:** Validate fixed local spatial supports for MLWF exchange before enabling them in HSE dynamics.

**Architecture:** For fixed real support projectors P_i, evaluate E_x[{P_i w_i}] and its derivative P_i K[{P_j w_j}]P_i w_i. Pair densities and output are confined to intersections of translated supports. Retain the sampled periodic HSE kernel through zero-padded local convolution, and reuse reciprocal pairs. This is an orbital-dependent functional: energy-gradient consistency alone does not imply a common Hermitian operator or compatibility with ACE. Measure the occupied metric before any dynamical integration.

**Tech stack:** Python, NumPy, existing single-thread FFTW, unittest.

## Approved scope and decisions

The user approved local domains, consistent energy/action, fixed domains within nonlinear iterations, and comparison with full support. First validate fixed domains, then consider adaptive domains. Keep the running full-support reference unchanged. Work on the existing user-requested TDCDFT branch; add isolated sample modules, not changes to the active propagator.

## Task 1: Mathematical tests

Create `samples/hse_mlwf_reference/test_local_support.py`. Test full-support equality with independent full convolution, boundary-wrapped intersections, finite-difference energy gradient at fixed supports, vanishing disjoint pairs, and invalid domains. Run unittest and confirm the new API is missing before implementation.

## Task 2: Local exchange

Create `samples/hse_mlwf_reference/local_support.py`. Store per-orbital cubic supports. For each reciprocal pair, intersect supports in periodic coordinates, convolve only a containing local box, and scatter both gradient contributions. Cache FFTW plans by box dimensions. Validate translation group and inputs. Return the gradient and statistics; explicitly do not advertise a generic Hermitian operator. Run Task 1 tests with NumPy and FFTW.

## Task 3: HSE benchmark and integration gate

Create `samples/hse_mlwf_reference/benchmark_local_support.py`. Load the persistent HSE ground state, compare widths with full exchange using equal FFTW settings, record setup and repeated timings separately, energy/action error, discarded norm, and per-k occupied-metric anti-Hermiticity. Attempt the unmodified ACE constructor and record acceptance/rejection. No symmetrization to hide an incompatibility. Use fixed supports during tests. Repeat on an available completed impulse snapshot if appropriate.

## Task 4: Review and report

Run the full Python suite and review the local-support change. Save results and explain whether it is suitable for ACE/PT-CN. Adaptive updates and RT integration are conditional on the mathematical gate, not automatic consequences of a fast convolution.

## Progress ledger

- Initial inspection: old source-only local support is not a variational SCF gradient. Existing full-support HSE job remains active.
- Ruling: test ACE compatibility explicitly before integration; MLWF-specific projectors break occupied-unitary invariance and may prevent a Hermitian ACE interpolation.
- Task 1 complete: observed missing-API failures, then fixed-support gradient, wrapped-domain and full-support tests pass. Corrected a test's pair-count bound from 20 to 24 (mesh2 self-inverse translations).
- Task 2 complete: local FFT, disjoint-pair omission and reciprocal reuse implemented in a separate experimental module. Added mesh3 inverse-pair regression and an explicit ACE rejection example.
- Task 3 complete: HSE ground-state widths12–32 measured; all fail existing ACE Hermiticity gate. Width48 (no cutoff) matches the original operator at 1.75e-16 and passes ACE. Reference time 11.75–12.35s, width16 0.975s but 6.05% action error; width32 13.657s and 1.39% error.
- Review: independent random complex mesh3 direct-convolution checks agree at 1.6e-16–4.7e-16. Addressed snapshot fingerprint race by hashing the same byte snapshot that is loaded; new test observed red then green.
- Ruling: do not integrate orbital-specific projected gradients into ACE/PT-CN or add adaptive domains after the Hermiticity gate fails. That would change the propagation model; a correct fixed-domain gradient alone is insufficient. Report this measured limitation rather than silently symmetrizing the occupied metric.
- Task 4 complete: independent review found no mathematical blocker for the experimental scope; 58 Python tests pass. Results and rejected integration gate documented in docs/results/si-hse-local-support/README.md. No production propagation changes.

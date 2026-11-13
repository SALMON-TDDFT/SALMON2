# Velocity-gauge SBE adapter validation

The original implementation is SALMON commit
`3322d0561c0cad5180679136b3eddbe8b909487c`.
`tools/test_sbe_native_regression.py` extracts its unmodified propagation,
current and trace routines, builds a separate oracle, and compares against the
new shared kernels before writing `salmonpy/tests/sbe_reference.npz`.

Local macOS arm64 / GNU Fortran 15 / Accelerate results:

- 4 current-correction orders × 2 nonlocal-correction flags × 41 endpoints:
  J and every density-matrix entry bitwise identical to the original oracle.
- Full SALMON executable: 16 configurations including both 2- and 3-band
  propagation. Current, energy and excitation output files byte-identical
  before/after extraction, OMP_NUM_THREADS=1.
- Python loader + common A→J entry: currents agree with full SALMON output
  in all 16 configurations (`rtol=3e-13`, `atol=3e-15`).
- Python suite: 30 tests pass, including independent states, failure atomicity,
  immutable initial arrays, text indices/units, partial-band loading and
  portable oracle comparisons (rtol=2e-13, atol=2e-15).
- Full application build and forbidden-eigensolver source/object audits pass.

Audited before executable SHA256:
`36ddc1b5d6527c89f0718be0e81cec0ae9fb2fe86fc5632f0b0551e961086c79`

After executable SHA256:
`32f9123cd400aab16750e6ee4c16253c26f257a99824d25dd957379088f345f7`

Reproduce the integration test with previously audited executables:

```
PYTHONPATH=. python tools/test_sbe_entrypoint.py BEFORE AFTER NEW_DIRECTORY
```

Fixtures are synthetic, not a physical material prediction or new GS calculation.
No exact-unitarity, adaptive-error, MPI reproducibility or length-gauge claim is
made. SALMON's existing MPI reduction stays outside the extracted local-current
kernel. The native Python adapter is serial and uses independent work arrays.

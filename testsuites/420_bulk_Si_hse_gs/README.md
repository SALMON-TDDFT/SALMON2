# Si HSE GS → impulse regression (420/421)

Si8, 12^3 real-space grid, shifted full 4^3 k grid, 32 electrons, 16 occupied
states, standard HSE06. Case 420 starts without external checkpoint data and
requires the density residual below 1e-8 (up to 160 SCF iterations). Case 421
requires successful GS verification and checks 64 Taylor4+ACE steps, dt=0.16 au,
z impulse 1e-4 au, finite response output, bounded current, and post-kick energy.
The dielectric output is checked against the conductivity relation.

The 10.24-au observation time is deliberately short for CI; the resulting
Fourier output is not a resolved optical spectrum. No experimental peak
position or k-grid convergence is asserted.

```sh
OMP_NUM_THREADS=1 ctest --test-dir build \
  -R '(420_bulk_Si_hse_gs|421_bulk_Si_hse_rt)' --output-on-failure
```

The tests use four MPI ranks. In a Fugaku compute-node allocation, set
OMP_NUM_THREADS=12 for the 4-rank × 12-thread configuration and use the site's
MPI launcher. Retain the same setting in both producer and consumer. Python 3
is used only by verification. CTest selecting case 421 automatically adds the
GS producer and verification through fixtures. Failed GS convergence blocks RT.

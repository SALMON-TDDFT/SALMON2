# Native HSE numerical regression tests

These small independent Python/Fortran fixtures test native exchange, ACE,
MPI row distribution, symmetry expansion, semilocal Libxc evaluation and
occupied-subspace observables. They do not run reference time propagation.

```
SALMON_TEST_MPI=1 OMP_NUM_THREADS=1 python3 -m unittest discover -s samples/hse_mlwf_reference -p 'test_native*.py'
```

Requires NumPy, a Fortran compiler, FFTW, BLAS and Libxc. Override FC,
FFTW_ROOT, OPENBLAS_ROOT and LIBXC_ROOT for your installation. MPI tests
require mpifort/mpiexec. End-to-end namelist and restart tests are provided
in ../hse_native/check_input_smoke.py and accept an existing converged Si GS.

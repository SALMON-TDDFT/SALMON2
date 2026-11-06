# Native CPU HSE build

CPU builds enable `USE_HSE` by default; accelerator builds leave it disabled.
The native exchange path requires FFTW, Libxc, BLAS/LAPACK, and optionally MPI.
CMake checks installed libraries and can prepare compatible fallback dependencies.
`USE_HSE=OFF` retains the ordinary non-hybrid build.

Run the source guard before configuring or building; stop if it fails:

```sh
python3 tools/check_forbidden_eigensolvers.py --source-root src &&
cmake -S . -B build -DUSE_MPI=ON -DUSE_HSE=ON -DCMAKE_BUILD_TYPE=Release &&
cmake --build build --parallel 4
```

CMake retains source and application-object guards during incremental builds.
Before deploying a new executable, inspect every linked application object with
the same checker (`--object-tree`), including objects outside the main build if
linked manually. Vendor-library exports alone do not prove an application call.
`PZHEEV` is prohibited; scalar IEEE inquiries are required. See [coding rules](../CODING_RULES.md).

The sole proposed regression system is the small Si cell, prepared as a GS
producer and RT consumer using the ordinary testsuite interface:

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 ctest --test-dir build -R '42[01]_bulk_Si_hse' --output-on-failure -j1
```

It requires four MPI ranks. Fixture dependencies produce a converged GS before
the 64-step RT run. This is a bounded numerical regression, not an optical
spectrum or performance benchmark. See [validation scope](PR_VALIDATION.md).

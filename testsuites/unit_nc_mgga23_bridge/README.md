# Discrete noncollinear meta-GGA bridge validation

This test uses synthetic complex occupied spinors, unequal occupations, nonzero
Bloch k, a skew Cartesian gradient transform, and allocated zero NLCC core arrays.
It performs no ground-state calculation and reads no material calculation data.

It checks the occupied-state XC band expectation against the stored contraction,
Hermiticity of the complete frozen local/tau/current operator, a centered orbital
finite difference of the complete discretized XC energy against the Hamiltonian
action, and energy invariance under a constant SU(2) rotation. Separate subprocesses verify explicit rejection of
nonzero core density and nonzero core kinetic energy density. Periodic indices
exercise derivative wrapping. The test asserts that SALMON `hpsi` adds the complete operator exactly once and
that its optional kinetic-only output excludes XC. A synthetic self-consistent
midpoint evolution uses SALMON's sixth-order Taylor propagator and checks time-step
convergence, norm, global spin without SOC, and payload reconstruction from the
propagated state. This is a numerical integration test, not a material benchmark.

The test also calls `init_xc` and `exchange_correlation` with `nc_mgga23` and
checks the dispatcher energy and zero scalar potential. It executes the actual
`write_rt_spin` writer on a synthetic propagated state; the Python runner verifies
its diagnostic header, six finite samples, and integral local spin >= total spin.
All output is written inside the temporary runner directory and deleted.

The optional `--mpi` mode splits occupied orbitals between exactly two MPI ranks.
The `--rspace` mode runs a dedicated two-rank x-domain regression against an
independently evaluated full-grid reference. It uses actual scalar/orbital
`init_sendrecv_grid` communication, four-cell halos, all four derivative stencil
radii, nonzero k, and a skew gradient transform. It compares the complete XC action,
energy and band expectation; y/z use periodic index wrapping.

Run against a current serial build (compiler and library options may differ):

```sh
python3 testsuites/unit_nc_mgga23_bridge/run.py \
  --build-dir /private/tmp/salmon-ncmgga23-serial \
  --checker /absolute/workspace/tools/check_forbidden_eigensolvers.py \
  --fc /opt/homebrew/bin/gfortran \
  --libs /opt/homebrew/opt/openblas/lib/libopenblas.dylib
```

The runner audits the current source, all linked application objects, and its new
test object before linking. Temporary artifacts are removed after execution.

NLCC restriction: nonzero core density or kinetic energy density is explicitly
rejected. A density-only frozen core lacks the consistent core tau/current gauge
response required by this functional. Allocated zero core arrays are allowed.
A physically consistent meta-GGA core construction remains future work; this
functional cannot yet be used with the Ni SOC+NLCC pseudopotential.

For the MPI build use `--fc mpifort --mpi --launcher mpiexec` and its build path.
Local MPI launch may require permission to bind local sockets.

For the dedicated halo regression, replace `--mpi` with `--rspace`. This mode uses
16³ synthetic grids and does not run the diagnostic writer or RT checks, which
remain in the serial/orbital-decomposition driver.

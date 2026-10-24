# KY V2 auxiliary-Pauli core reader and grid regression

`test_reader.py` extracts the production reader into the established minimal
reader fixture. It preserves V1 regression tests and rejects invalid V2 model,
units, angular partners, nonfinite scalar data and inconsistent primitive
invariants. Its V2 fixture contains paired regular s/p auxiliary primitives.

`run.py` links the full application objects from a current audited SALMON build.
It exercises the actual reader, core-table MPI broadcasts, and `calc_nlcc`.
Checks include retained V1 density/projector equality, full V2 core density tails,
core tau interpolation, Cartesian spin-current orientation and adaptive periodic
image support. A small cell makes contributions outside the old +/-2 images
numerically detectable. V1 retains both historical density truncations; V2
retains complete joint rho/tau/current support.

For V2 data it also constructs bounded smooth synthetic spinors on a 12-point
mesh and calls the actual noncollinear meta-GGA refresh and Hamiltonian action
using the loaded frozen core. This is a reader-to-operator coupling test, **not a
material trajectory, ionic SOC calculation, or ground-state benchmark**.

```sh
python3 testsuites/unit_ky_core23/run.py \
  --build-dir /private/tmp/salmon-ncmgga23-serial \
  --checker /absolute/workspace/tools/check_forbidden_eigensolvers.py \
  --fc /opt/homebrew/bin/gfortran \
  --libs /opt/homebrew/opt/openblas/lib/libopenblas.dylib \
  --pseudo /absolute/path/to/Ni_soc_nlcc_rps.dat
```

For the MPI build add `--mpi --launcher mpiexec` and use `--fc mpifort`. Two ranks
must produce byte-identical tables, grids and coupling diagnostics. This MPI test
validates broadcasts; the separate bridge `--rspace` regression validates domain
halos.

The model block uses modern Bohr conversion 0.529177210903 Angstrom, while SALMON
retains 0.52917721067 for the original KY rows. Cross-unit radial checks allow
2e-9 relative radius and 2e-8 relative density error; pointwise primitive
rho/tau/current closure retains 2e-9 tolerance. Stored primitive derivatives are
checked with fourth-order interior differences at 1e-3 of their maximum magnitude,
with second-order differences at the adjacent endpoints. The actual Ni model's
observed fourth-order residual is below 3.6e-4; second-order differences alone
have appreciable truncation error near the smooth core matching radius.

# Host finite-value validation (NVHPC 25.9)

NVHPC 25.9 was reported to reject the whole-array expression
`all(ieee_is_finite(SymMatA))` in `symmetry_validate_group` with
`Illegal call from host code to device subprogram __pgi_ieee_is_finite_dev_r8`.
The workaround copies each operation element into a host scalar before the IEEE
check. It retains checks of both rotations and translations.

The test extracts the production validator, retaining the module-scope allocatable
array, and checks a valid identity group and NaN / positive infinity / negative
infinity in rotation and translation entries. It has no SALMON or MPI dependency.

Run from the repository root:

```sh
python3 testsuites/unit_symmetry_finite/test_finite.py -v
FC=nvfortran FFLAGS='-O3 -gpu=cc90' python3 testsuites/unit_symmetry_finite/test_finite.py -v
```

For diagnosis, set FFLAGS to the flags used to compile symmetry.f90. The test can
also be copied to an unpatched checkout to compare the original failure.
A passing GNU test verifies preserved validation behavior, not resolution of the
NVHPC compiler error. NVHPC 25.9 and the full MIYABI build still need verification.

# ASE Calculator and egg-box diagnostics — 2026-11-13

Native changes only report existing forces after calc_force at full precision.
The Calculator exposes native energy/forces with unit conversion and isolated fresh
GS runs. No native force formula, potential, orbital update or force replacement
was changed. Initial Calculator scope is PZ/scalar KY, single rank, no stress.

Mac arm64 real test: Si cubic 8 atoms, lattice 10.26 bohr, atom 1 displaced by
0.12 bohr along x, 32 electrons/states, 2³ k mesh, rho_dne < 1e-12. Scalar KY Si
potential from exercise_06. All calculations start fresh. This is a discretization
study; pseudopotential generating-XC provenance has not been independently verified.

| Grid | Native Fx (Ha/bohr) | Central energy derivative (Ha/bohr), h=.01 | Relative difference |
|---|---:|---:|---:|
| 16³ | -0.0153748097599 | -0.00907322568331 | 69.45% |
| 24³ | -0.0161205154915 | -0.0184855829312 | 12.79% |
| 32³ | -0.0160818338746 | -0.0173110069319 | 7.10% |

Native force changes by about 0.24% from 24³ to 32³. This is not a demonstrated
fully converged force. At 16³ h=.02 and .01 energy differences agree closely;
at 24³ h=.02 gives -0.0183506065987 Ha/bohr, showing additional step dependence.
Only h=.01 was tested at 32³. Numerical derivatives are not treated as exact truth.

Rigid translation preserves all interatomic distances. Half a grid step changes
energy by -0.2684471031 eV at 16³ and +0.0441188250 eV at 24³. A full grid step
returns within 3e-11 eV at both grids. This establishes an egg-box contribution,
not that every remaining force discrepancy has been uniquely explained by it.
The rigid x energy slope at 16³ is about 0.00736097 Ha/bohr while the native total
Fx is about 0.000808846 Ha/bohr; translational energy ripple contaminates derivatives.

ASE BFGS smoke: 24³, atoms 2–8 fixed, maximum step .01 Å, exactly one optimizer
step. Energy decreased from -864.665562136653 to -864.6743497788104 eV. Maximum
force decreased from about .828950 to .699121 eV/Å. Optimizer convergence is false;
this validates connection and constraints, not a converged geometry or force accuracy.

Unit tests cover cache reuse, position/cell invalidation, unique result directories,
unsupported stress/settings, failure clearing stale values, copied returned arrays,
force order/shape/units and nonfinite values. ASE remains an optional separate module.

Reproduce from source using active Python and optional ASE dependencies:

```
python tools/check_force_grid.py SALMON_EXECUTABLE SI_RPS_DAT NEW_DIRECTORY --grid 24
python tools/test_ase_calculator.py SALMON_EXECUTABLE SI_RPS_DAT OTHER_NEW_DIRECTORY
```

The diagnostic evaluates fresh SCF calculations; account for memory, parallel budget
and runtime first. Scientific tests were local Mac only, not cross-platform native CI.

Scientific context: https://pubs.rsc.org/en/content/articlehtml/2015/cp/c5cp00351b

## Packaged Python verification

The repaired macOS arm64 wheel passed all 64 tests with ASE installed.
Without ASE, 58 tests passed and 6 optional ASE tests were skipped.
This includes the beginner resource-intake policy returned by manual search.

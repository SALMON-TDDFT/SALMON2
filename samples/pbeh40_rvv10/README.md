# PBEh(40)+rVV10 small water fixtures

Use an HSE-enabled CPU build. Copy `testsuites/pseudo/H_rps.dat` and `O_rps.dat` into a fresh run directory. Run `salmon < water.inp` for a static coarse-grid check, or `salmon < water_md.inp` for eight 0.05 fs fixed-cell BOMD steps.

Both use the inherited MLWF+ACE exchange path with 40% Coulomb exchange. The examples are a single periodic water molecule, not liquid-water production inputs. They deliberately retain finite grid and cell errors.

[Model, parameters, cutoff convention, limitations and validation](../../docs/inputs/pbeh40-rvv10.md).

`validate.py` repeats SCF at small H/O displacements and checks force against the energy derivative, or checks short NVE energy drift at two timesteps. Set `--atom O` for the oxygen displacement; `--xc pbeh40` omits rVV10. It refuses existing case directories and fails on unconverged SCF or failed numerical checks.

# PBEh(40)+rVV10 small water fixtures

Use an HSE-enabled CPU build. Copy `testsuites/pseudo/H_rps.dat` and `O_rps.dat` into a fresh run directory. Run `salmon < water.inp` for a static coarse-grid check, or `salmon < water_md.inp` for eight 0.05 fs fixed-cell BOMD steps.

Both use the inherited MLWF+ACE exchange path with 40% Coulomb exchange. The examples are a single periodic water molecule, not liquid-water production inputs. They deliberately retain finite grid and cell errors.

[Model, parameters, cutoff convention, limitations and validation](../../docs/inputs/pbeh40-rvv10.md).

`validate.py` repeats SCF at small H/O displacements and checks force against the energy derivative, or checks short NVE energy drift at two timesteps. Set `--atom O` for the oxygen displacement; `--xc pbeh40` omits rVV10. It refuses existing case directories and fails on unconverged SCF or failed numerical checks.

## Real-space Ehrenfest from a DC initial state

Run `salmon < water_ehrenfest_gs.inp` first, then
`salmon < water_ehrenfest_rt.inp` in the same fresh directory, containing the
H/O pseudopotentials and `water_ehrenfest_velocity.dat`. The GS stage writes
`data_dcdft`; the RT reader verifies its functional/run metadata and reconstructs
occupied mesh wavefunctions once. RT uses `yn_hse_lcfo_rt='n'`: no LCFO projection
or basis propagation occurs. MLWF/ACE accelerates exchange only.

This is a one-water, one-rank, 24³-grid smoke fixture: four 0.0005 fs steps after
an electronic impulse with NVE ions. It does not validate liquid-water dynamics
with a finite-duration laser. The Gamma native RT exchange now supports y/z spatial decomposition.
Orbital decomposition and giant-system performance validation remain outstanding.


For a finite pulse, use theory='tddft_pulse' and ae_shape1='Acos2', with positive
omega1/tw1, nonnegative t1_start, and linear epdir_re1. Frequency, duration and
amplitude follow the selected unit system. The executable H4 regression
`testsuites/unit_pbeh_rvv10/test_ehrenfest.py` contains the validated pulse input
and independent external-work integration; the water file above remains an
impulse smoke test.

The saved H4 inputs are `h4_ehrenfest_gs.inp` (two MPI ranks), followed by
`h4_ehrenfest_pulse.inp` (one rank) in the same directory with H_rps.dat and
`h4_ehrenfest_velocity.dat`. They reproduce the finest pulse test in atomic
units: dt=.02, 480 steps, width6.4, omega1=1.9634954084936207 and amplitude.03.


`h4_ehrenfest_pulse_spatial.inp` uses four MPI ranks with `nproc_rgrid=1,2,2`.
Prepare the same `data_dcdft` with `h4_ehrenfest_gs.inp` (two ranks), then run
`mpiexec -n 4 salmon < h4_ehrenfest_pulse_spatial.inp`. The RT stage stores local
wavefunction, MLWF-source and ACE rows; it uses distributed FFTW exchange.
When comparing serial and spatial outputs, use separate run directories with
the same GS data and pseudopotential/velocity files to preserve both outputs.
This full-support route does not yet support finite-radius MD or orbital MPI.

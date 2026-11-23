# Si starter: ground state to current

[日本語](SI_INTRO.md)

Compute the fixed Si SCF input and compare 21 samples of A/c-driven charge current against the reference.

## Required distributions

Obtain an updated SALMONpy wheel for your OS, CPU and Python, plus an outer runtime ZIP listed in the current `runtime_artifacts` of `provenance.json`. NumPy is the only Python dependency of this example.

The verified distribution is for Apple Silicon Macs. It is not available as a public Release; obtain the matching wheel and ZIP from the distributor. Check the filenames, hashes and platform conditions in `runtime_artifacts`. If no matching ZIP is available, use the wheel-only [LD/SBE test](README.en.md).

Intel Mac and Linux use requires a matching validated distribution for this example. No Windows DFT/TDDFT distribution is provided.

## Install, prepare, run

```sh
python -m pip install --only-binary=:all: ./SALMONPY_WHEEL.whl
python -m salmonpy.quickstart
python -m salmonpy.examples.si_intro prepare new-si --runtime-artifact RUNTIME_ARTIFACT.zip
python -m salmonpy.examples.si_intro run new-si --threads 1 --memory-budget-gib 4
```

Replace the wheel and ZIP placeholders with actual supplied files. Pass the outer ZIP without unpacking it. `prepare` verifies and places files; it never starts a calculation. Prepare on a management machine and invoke `run` within allocated compute resources. There is no implicit download/build or overwrite of existing directories.

Threads may be 1–4. The 4 GiB budget describes the whole job and is not an OS-enforced cap. Unavailable resource measurements are separate from calculation failure and do not certify compliance with the budget.

## Inspect results

Both `calculation` and `validation` must be `passed` in `new-si/result/report.json`. Resource measurements appear separately under `resources`.

- GS and restart: `result/gs/`; GS log: `result/gs/native/stdout.log`.
- Current: `result/current.csv`, `result/current.npz`; response archive: `result/rt/`.
- CSV columns: `time_au,Ac_x_au,Ac_y_au,Ac_z_au,j_x_au,j_y_au,j_z_au`.

Output is atomic-unit physical charge-current density **j**, with solver matter current `J = -j`. Do not negate the returned current again.

Expected GS energy is approximately −31.7465390575 Hartree. Acceptance uses energy absolute error 1e−7 Hartree and component-wise current `atol=1e-9, rtol=1e-5`; time differences are roundoff-only. Iteration count is diagnostic. The source run on 2026-11-23 took 54 iterations and matched the expected energy and all 21 current samples exactly.

## Scope and next calculations

The case is cubic Si with 8 atoms, cell length 10.26 bohr, 32 electrons/states, a 16³ grid, 2³ k points, PZ, and SCF threshold 1e−10. It propagates only 20 steps at dt=0.02 a.u. (0.4 a.u.). This checks installation and interfaces, not a realistic laser pulse or converged optical properties.

For research, extend the pulse to several periods and the observation time to the target phenomenon. Use [convergence diagnostics](CONVERGENCE.en.md) to assess grid, k points and dt. Grid screening compares current waveforms at converged dt, using 5% relative RMS as the default; dt screening uses measured relative electron-number decrease of 1e−5. These do not certify dielectric or harmonic spectra.

The Si potential is already bundled in the wheel. For separate acquisition and other referenced files, see [pseudopotential downloads](PSEUDOPOTENTIALS.en.md).

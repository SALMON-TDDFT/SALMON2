# Short Ni SOC test

[日本語](NICKEL_SOC.md)

The small SOC CI case uses a one-atom fcc primitive cell and SALMON's standard 2×2×2 grid (8 k points). Python does not generate an alternative k-point set.

## Potential and conditions

- PseudoDojo PBE v0.4 fully relativistic NC-UPF `Ni_std.upf`, with SOC and NLCC.
- Exact catalog ID: `pseudodojo-pbe-0.4-fr-standard-ni`; commit and SHA256 are pinned.
- **18 valence electrons including 3s/3p semicore**. The atomic reference configuration is 3s² 3p⁶ 3d⁸ 4s², not a constraint on solid-state band occupations.
- 36 spinor states, 16³ grid, Libxc-PBE, electronic temperature 300 K.
- Seeded random spinors, seed0; no prescribed magnetization direction or moment.
- `simple_dm` mixing with rate0.2; exactly three SCF iterations.

Grid, k mesh, temperature and lattice constant serve a short execution test. They do not establish converged Ni magnetism, anisotropy or optical properties.

## Explicit execution

Use a `build.json` from a Libxc-enabled build with the required source/object eigensolver audits:

```sh
python tools/test_nickel_soc.py BUILD_JSON NEW_OUTPUT --potential NI_UPF
```

Replace the placeholders with actual paths. **This command starts a short calculation**, unlike input preparation. Run within allocated compute resources or a local test environment, not on a management node. Existing output directories are rejected.

Replace `--potential NI_UPF` with `--allow-download` only when fetching the pinned file is explicitly authorized. For offline compute nodes, acquire the [UPF beforehand](PSEUDOPOTENTIALS.en.md) on a connected preparation machine.

The wheel includes these tools under `salmonpy/_runtime_source/tools/` and `tests/runtime/`. Importing SALMONpy or running the initial LD/SBE examples does not invoke this compilation/calculation.

## Assertions and output

Fortran checks finite density/XC/mixing arrays, Hermitian density matrices, electron count18, independent noncollinear-array copies and exact restoration. Cases use one and two threads.

Nonconvergence after three iterations is recorded as expected for this smoke test; it never triggers automatic RT continuation. No wavefunction, density or GS restart is exported. `result.json` and `execution/test-record.json` contain ordinary logs, residuals and outcomes.

**Restoration here concerns noncollinear work arrays, not the entire SCF/TDDFT state or a disk restart.** End-to-end SOC checkpoint support remains separate work.

## CI scope

The existing manual `NLCC native consistency` workflow reuses its Libxc build for this test. There is no extra solver build, automatic physical-property campaign or new workflow. Local Apple Silicon execution passed on 2026-11-23. Remote GitHub Actions execution has not been performed.

## Independent SCF states and memory snapshots

`tools/test_nickel_scf_mpi.py` compares interleaved states with different cells and mixing settings against independent runs.
It uses 18 electrons, 36 states, 300 K and SALMON's standard 2×2×2 k mesh.

```sh
mpiexec -n 3 python tools/test_nickel_scf_mpi.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 1
mpiexec -n 6 python tools/test_nickel_scf_mpi.py BUILD_JSON NEW_OUTPUT_MPI --potential NI_UPF --size 2
```

Mac tests on 2026-11-23 matched density, occupations, eigenvalues, orbitals and residual history exactly with one or two internal MPI ranks.
Restoring performed zero SCF iterations. Projector construction scratch is call-local; mutable noncollinear arrays belong to each SCF state.
No orbital, density or GS restart files are exported. This short unconverged SCF test does not validate SOC disk restart, RT, forces or fast HSE.

## SCF restart after process exit

```sh
python tools/test_nickel_scf_disk.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 2
```

This manual test saves electronic states, mixing history and comparison arrays. Confirm the destination and permission to retain data first.
It saves after two SCF iterations, restores in a fresh MPI process and runs two more iterations against a continuous four-iteration control.
Mac tests on 2026-11-23 matched density, occupations, eigenvalues, orbitals and residual history exactly with one or two internal MPI ranks.
Restoration executed zero SCF iterations. Stream format2 includes noncollinear arrays; old streams are not accepted.
This is a short unconverged SCF test with the same binary and MPI layout. SOC rank-count changes, RT, forces and fast HSE remain unvalidated.

## SOC A-to-current response and direct checkpoints

```sh
mpiexec -n 3 python tools/test_nickel_tddft_mpi.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 1
mpiexec -n 6 python tools/test_nickel_tddft_mpi.py BUILD_JSON NEW_OUTPUT_MPI --potential NI_UPF --size 2
```

This manual Ni integration probe uses 300 K,18 electrons,36 states and SALMON's original2×2×2 k mesh.
Its four-iteration unconverged GS stays in memory and is only a software-comparison seed, not a physical result.
Two independent states receive different three-component A drives for eight steps and match independent controls exactly.
Electronic states saved at step4 are directly restored into fresh states; currents at steps5–8 match exactly.
Mac tests pass with one/two internal MPI ranks. Scalar Si currents also match the preceding library exactly.

This test retains electronic checkpoints; confirm storage permission and destination. Use `--without-checkpoint` to compare live states only.
The validated case is fixed-ion,full-k,CPU Libxc-PBE Ni. Noncollinear arrays belong to each state.
Crystal input is A/c and output is charge current. Scalar equations and update order are unchanged.
The TDDFT stream is format3, without implicit old-stream conversion.
Fresh-process restart is also validated below. Changed MPI layouts,
deep numerical failure recovery,fast HSE and physical convergence remain unvalidated. WSe₂ validation is separate.

## TDDFT restart after process exit

```sh
python tools/test_nickel_tddft_disk.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 1
python tools/test_nickel_tddft_disk.py BUILD_JSON NEW_OUTPUT_MPI --potential NI_UPF --size 2
```

This manual test retains GS and RT electronic data; confirm storage permission and destination first.
For each MPI layout, the first process saves a four-iteration software-test GS once and computes an eight-step continuous control.
The next process loads that GS, runs RT steps1–4, saves and exits.
The final process loads GS and RT data and executes only steps5–8.

Mac tests pass with one/two internal MPI ranks. Both reload processes execute zero SCF iterations;
RT loads step4 and executes four additional steps. Charge currents and measured electron numbers match the continuous control exactly.
The binary and MPI layout are unchanged. Work arrays, communicators and potentials are initialized, but SCF iterations and earlier RT steps are not replayed.
The short unconverged-GS probe does not establish physically converged Ni properties.


## Rollback and retry including output streams

```sh
mpiexec -n 2 python tools/test_nickel_tddft_retry.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 2
```

The unconverged comparison GS remains in RAM; no GS/electronic checkpoint is exported.
The probe requests per-step spin output and raises a handled host exception after a step.
With one/two internal MPI ranks, current/energy/spin files return to their original bytes,
retry currents/magnetization/spin currents match continuous controls, and sibling outputs are preserved.
This does not validate deep numerical failures, MPI failure or process-crash recovery.
See [reading and saving examples](OUTPUTS.en.md).

## Checked and performance modes for SOC

```python
# model = gs.tddft_model(dt=..., max_steps=...)
audited = model.new_state('audit', execution_mode='checked')
fast = model.new_state('fast', execution_mode='performance')
```

| Mode | Per-step behavior | Failure contract |
|---|---|---|
| `checked` (default) | Input/MPI agreement, finite checks for wavefunctions, density, Hartree, XC/SOC arrays and observables; numerical backup | Supported numerical failures automatically restore the accepted state and preserve current/energy/spin files. `NumericalError` permits retry after correcting the cause |
| `performance` | Omits numerical diagnostic checks, diagnostic agreement and automatic state copies inside the update | No automatic recovery guarantee. Explicitly use `guarded_step` or checkpoints where protection is needed |

Checked mode does not certify physical convergence. Source/object forbidden-eigensolver audits remain mandatory at build time in both modes.
Reuse of an unchanged calculation library verifies its recorded SHA256; source audits are not repeated at job startup or every step.
SOC-array failures join the existing XC-stage agreement; magnetization/spin-current checks join the existing observable-stage agreement.
Observables are computed before output and reused by the writer, without duplicate spin-observable computation.
Performance updates do not invoke these checks or the automatic backup.

For a normal production library:

```sh
mpiexec -n 2 python tools/test_nickel_tddft_retry.py BUILD_JSON NEW_CHECKED --potential NI_UPF --size 2 --execution-mode checked
mpiexec -n 2 python tools/test_nickel_tddft_retry.py BUILD_JSON NEW_FAST --potential NI_UPF --size 2 --execution-mode performance
```

Fault injection requires a separate test-only build with `SALMONPY_RECOVERY_TESTING`; add `--deep-faults` to the checked-mode probe.
Production libraries do not export this hook. Mac Ni Libxc-PBE tests cover real overflow and wave/density/Hartree/XC/observable faults
with one/two internal MPI ranks. Accepted state and all three stream files are preserved; retry currents, magnetization and spin currents
match controls exactly, with sibling isolation. Performance results match checked results and explicit protection restores handled host exceptions.
The validated profile is fixed-ion, full-k, CPU k-point MPI. MPI/process-crash recovery is excluded.
The short unconverged-GS software probe is not a physical convergence test.

For auxiliary-field SCF stages and MPI1/2 restart tests, see [SOC auxiliary spin fields](NATIVE_GS.en.md#stabilize-soc-scf-with-an-auxiliary-spin-field). Long fixed-input trials show alignment, but final zero-field SCF convergence remains unverified.

## Manual SCF comparison with1000-iteration stage budgets

`tools/test_nickel_scf_mpi.py --alignment-study` compares a zero-field control
with staged removal of an auxiliary field, separately from the short state tests.
This command **runs calculations**. It is not a management-node input checker;
run it explicitly in a calculation environment or allocated job.

```sh
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 mpiexec -n 2 python tools/test_nickel_scf_mpi.py \
  BUILD_JSON NEW_OUTPUT --potential PINNED_NI_UPF --size 2 \
  --alignment-study --stage-limit 1000 --final-stage-limit 2000
```

The fixed case uses18 electrons,36 spinors,300 K,16^3 points and SALMON's2x2x2
k mesh. It explicitly selects `simple_dm`, coefficient0.3 and downward adaptation.
Seed0 initial arrays are compared with recorded tolerances/differences; threading
can introduce roundoff. The auxiliary schedule is Bz=0.002,0.0005,0.0001,0 (au),
up to1000 iterations per stage; the zero-field control has a1000-iteration budget.
The current magnetic residual target is1e-5, configurable with `--charge-tolerance`.
Final zero-field observation changes use1e-4 mu_B for moment components and1e-5 Ha for energy,
configurable with `--moment-tolerance` and `--energy-tolerance`. Preliminary nonzero-field stages
use residual1e-5 and moment changes1e-3 mu_B (`--preliminary-residual`, `--preliminary-moment-change`).
Final mixing is0.15 (`--final-mixing-rate`); `--final-stage-limit` sets its iteration budget.
Sampling switches from25 iterations to every iteration below ten times the stage residual target.
A convergence-confirmation call with no accepted iteration retains the last genuine changes.
A stage that misses its targets does not advance. `--control-record PATH_TO_ZERO_CONTROL_STAGES_JSON`
explicitly reuses a matching control after checking its input, potential and configuration; provenance is recorded.

Outputs include `history.csv`, `progress.json`, stage JSON, residual CSV/optional
plots, `scf-alignment.png` and `result.json`. Wavefunctions, density and restart
state remain only in RAM and are not retained for another process.
This studies SCF stability on a fixed coarse input, not grid/k convergence or a
certified Ni property. The long calculation is excluded from required CI.

### Historical observations with the strict1e-9 target (2026-11-23)

With MPI2 and4 threads per rank, the zero-field control had residual about3.0e-2
at iteration1001, versus about4.0e-6 for fixed auxiliary Bz=0.002 at iteration1001.
Using separate preliminary targets allowed Bz=0.002,0.0005,0.0001,0 stages.
After1000 zero-field iterations, at total iteration2801, the residual was about5.9e-6,
not the final1e-9 target. The spin moment was approximately(0.2844,0.2877,0.2893) mu_B;
its maximum component change over the last25 iterations was about1.5e-4 mu_B.
Alignment and staged removal were observed; final zero-field SCF convergence was not.
Downward adaptation reached0.01 in the zero-field control and0.15 after the ramp,
so this compares the complete strategies, not an isolated field-only effect.
Residual stagnation needs diagnosis rather than declaring completion by raising the iteration budget.

### Magnetic-material SCF threshold

Explicitly select `SCFSettings.for_magnetic(max_iterations=4000)` to start with a residual tolerance of `1e-5`. Override it with `tolerance=`. Ordinary `SCFSettings()` retains `1e-9`; SOC does not automatically imply magnetism.

Inspect energy and magnetic-moment changes as well as the residual. Reduce mixing progressively when oscillations persist. The Ni diagnostic uses mixing 0.15 for the final zero-field stage, energy change `1e-5 Ha`, and maximum moment-component change `1e-4 μB` between observations. These configurable diagnostics depend on the input and sampling interval and do not certify physical accuracy. Native SCF stops at its charge threshold; unsettled observables are reported for review rather than labeled converged.

2026-11-23: MPI2 × 4 threads reached all three auxiliary-field stage targets. After 1000 final zero-field iterations, residual was 1.97e-5, energy change 1.35e-5 Ha and maximum moment-component change 2.19e-5 μB over the last iteration. Residual and energy targets remained unmet; review is required. No electronic datasets were saved and no TDDFT was started.

2026-11-24: Extending the final zero-field limit to 2000 reached convergence after 1437 accepted zero-field iterations (3181 total), retaining mixing 0.15. Residual was 9.929e-6; the last actual iteration changed energy by 5.239e-6 Ha and moment components by at most 1.095e-5 μB. The auxiliary-field-free control remained unconverged.

### TDDFT from the converged GS

Add `--rt-steps 8` to pass the successful final zero-field GS directly from memory to TDDFT.
Read observation streams after closing states and the model. Use a rebuilt calculation library containing the normal SOC spin-output close fix. Charge current covers all eight updates; default spin observations cover steps0,1,8.
Checked and performance states receive three-component A/c; charge current, electron number,
magnetization and Pauli spin current are compared. Capture a memory backup midway, run the second half,
roll back and compare repeated values with `rtol=2e-12, atol=2e-14`; byte equality is recorded separately. No SCF replay or electronic disk checkpoint is used.
TDDFT is not started if final SCF targets are unmet. This short connection/recovery check does not
certify spectra or timestep accuracy.

Converged-GS trial on2026-11-24: internal MPI2 × 4 threads passed mode and step4 direct-memory rollback comparisons. Maximum current/electron-number difference was3.55e-15; flushed observation files, including magnetization and spin current, also agreed within the declared tolerance. This is not bitwise identity. GS counters were unchanged and no electronic disk data were saved. The auxiliary-field-free control remained unconverged; its overall review status is separate from the successful field-ramp/TDDFT route.

### Normal SOC output closure

Run `tools/test_nickel_tddft_retry.py BUILD_JSON NEW_OUTPUT --potential PINNED_NI_UPF --size 2 --close-outputs` explicitly with MPI2 ranks to verify that normal closure alone makes all five spin rows for steps0..4 readable and preserves closed sibling streams. This short software probe uses an unconverged GS, without an additional flush or protective backup. Internal MPI1/2 passed on2026-11-24; pre/post-fix current rows matched exactly.

## SOC orbital and mixed MPI validation (2026-11-24)

The fixed Ni input uses18 electrons,36 states,300 K and SALMON's2×2×2 k mesh. Internal `(k, orbital, x, y, z)` layouts tested are `(1,2,1,1,1)` with2 ranks and `(2,2,1,1,1)` with4 ranks. Checked mode verifies recovery from an actual overflow caused by a large finite A/c: accepted A, current, electron number and output remain intact, and a normal retry matches an independent control. Performance mode does not automatically detect numerical faults; it verifies explicit rollback after a handled host exception. Charge current, magnetization, spin current and sibling isolation are checked. Five spin samples are compared. These are software checks using a short unconverged GS, not physical convergence tests; electronic state files are not saved.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python tools/test_nickel_tddft_retry.py \
  BUILD_JSON NEW_OUTPUT --potential PINNED_NI_UPF --size 4 \
  --decomposition 2 2 1 1 1 --numerical-overflow
```

For performance mode, replace `--numerical-overflow` with `--execution-mode performance`. The test tool no longer caps rank counts at2; the decomposition product must match the launched rank count.

### CPU spatial SOC with orthogonal cells

The CPU calculation library aggregates distributed nonlocal projector overlaps and computes spatial spin current.
Python requires `sp_tddft_soc_spatial_protocol()` to return1 before enabling this layout; unsupported libraries are rejected before initialization.
GPU builds return0 and are outside this validation. SALMON does not support spatial splitting of the nonorthogonal primitive Ni cell; use cubic fcc Ni.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 2 python tools/test_nickel_tddft_retry.py \
  BUILD_JSON NEW_OUTPUT --potential PINNED_NI_UPF --size 2 \
  --decomposition 1 1 2 1 1 --conventional-cell --numerical-overflow
```

The cubic test has4 atoms,72 electrons,144 states,300 K,2×2×2 k points and a16³ grid.
Two-way x and y spatial splits verify Pauli0 versus matter current Jm, actual overflow recovery, output preservation and normal retry.
Performance mode verifies explicit protected rollback. Use `--decomposition 1 1 1 2 1` for y splitting.

### Checking outputs

The Pauli identity channel matches native matter current Jm at matching times and units. Python charge current is j=-Jm.
Read magnetization and spin current with their own time arrays and units.
For research calculations, check grid, k-point, timestep and simulation-duration convergence for the observable of interest.

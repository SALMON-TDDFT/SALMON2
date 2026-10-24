# Corrected noncollinear MGGA versus LSDA

`xc='nc_mgga23'` selects the adiabatic gKS implementation of Tancogne-Dejean,
Rubio and Ullrich, PRB 107,165111 (2023), corrected by PRB 110,199902 (2024).
Exchange is the noncollinear BR model with gamma=0.8; correlation uses the
corrected CS expression. It is not a spin-rotated collinear meta-GGA and does
not require USE_HSE. HSE06 and the existing LSDA selector remain available. For noncollinear
predictor/corrector propagation, LSDA now also averages its full XC matrix
between endpoints, matching the midpoint treatment of MGGA.

## Supported scope

Conventional CPU fixed-ion calculations with complex noncollinear spinors and
SOC representation; serial, orbital-decomposed MPI and real-space-decomposed MPI tested. RT requires
`propagator='middlepoint',yn_predictor_corrector='y'` and self-consistent XC.
Microscopic and band-decomposed charge/spin-current outputs are rejected because
the added frozen-core current has no assigned band decomposition.
Accelerators, DC, MD/optimization, symmetry reduction, frozen XC and other
propagators are rejected. Absorbing/non-Hermitian boundaries have not been validated. There is no OEP/KLI approximation.

## Consistent frozen-core NLCC

V2 KY files carry `CORE_AUX_PAULI_V1` metadata, rho, kinetic density (without
one-half) and the closed-shell spin-current tensor derived from the same
regularized auxiliary Pauli orbitals. The scalar charge current vanishes in
the ionic rest frame. In velocity gauge the code applies the same phase shift
to core kinetic density and current as to valence spinors, and includes the
explicit core field derivative in the macroscopic current. All core fields
enter XC together; band subtraction contracts valence fields only. Nonzero
legacy V1 NLCC is rejected for MGGA because its kinetic/current data are absent.
LSDA accepts the same V2 file and uses its rho; use V2 for **both** runs.

This is a declared frozen auxiliary-Pauli model, not an exact reduction of
Dirac core orbitals. It does not certify pseudopotential transferability or
remove the LDA-reference pseudopotential approximation in MGGA.
See [the model and validation plan](../../docs/plans/2026-10-24-nc-mgga23-nlcc-design.md).

The generated Ni candidate is in local OneDrive storage:
`/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/ni-soc-nlcc-20261024/v2-aux-pauli-validated/Ni/Ni_soc_nlcc_rps.dat`.
Its 16 valence electrons, s-local channel and s/p/d radii 3.2/2.8/3.2 bohr
match the selected V1 potential. All existing numerical radial rows are
unchanged. Provenance and radial screening are described in
[NI_VALIDATION.md](../../utility/crps_soc_nlcc/NI_VALIDATION.md).

## Matched RT inputs

The `functional_*.inp` files are fragments. `Ni_*_rt.inp.template` provides
matched four-atom fcc Ni RT inputs using the generated V2 potential. Replace
the checkpoint placeholder with an existing compatible state and match every
structural/discretization/occupation parameter to that state before running.
The supplied mesh, bands, k sampling, timestep and pulse are starting values,
not converged material settings. Run each in its own directory. First use
zero pulse intensity for the reference and halve dt for the convergence test.
The templates have deliberately unresolved checkpoint paths; no GS calculation
or GS data collection is launched by this workflow.
For each functional use its own consistent initial state. Keep pseudopotential,
geometry, k grid, real-space mesh, pulse, occupations and propagation settings
matched. Reusing an LSDA state with MGGA is a functional quench, not an
independently equilibrated comparison; label it explicitly if intended.
No material initial states or GS data are supplied/collected by this feature.

Use a short timestep-convergence study and a no-pulse reference before
interpreting laser-induced changes. The trajectory tests use synthetic
spinors without a physical SOC potential. Actual Ni V2 reader/grid tests check
the core data path, but do not establish material accuracy.

## Output and analysis

Existing `SYSname_rt_spin.data` retains its established columns and spin units.
It now includes the initial/restart sample before propagation. The comparison
records the first sample time explicitly as `reference_time`; restarted or
older files starting after zero are normalized to that time, not to t=0.
The added `SYSname_rt_magnetization_norm.data` contains time, the magnitude of
total spin and the integral of local spin magnitude. This separates coherent
rotation, loss of local moments, and cancellation between local directions.
For MGGA it additionally reports the integrated LOCAL multiplicative
commutator in the standard Pauli frame; tau/current contributions are excluded,
so it is not the total XC torque. LSDA's local commutator is zero.
Current columns in `rt_spin.data` retain the legacy bare kinetic/ionic
diagnostic; they exclude XC/core field response, as marked in the MGGA header.
Use the standard macroscopic RT current stream for corrected charge current.
Native legacy spin-output axes/normalization are preserved; no mu_B conversion
is inferred. Antiferromagnetic order requires a defined sublattice partition,
not just total magnetization.

```sh
python3 tools/compare_spin.py \
  --run LSDA=/absolute/path/lsda/SYSname_rt_spin.data \
  --run MGGA23=/absolute/path/mgga/SYSname_rt_spin.data \
  --output /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/magnetic-comparison/example \
  --plot
```

The output directory must be new. JSON/CSV are always produced; PNG plotting
requires matplotlib. Time units and time grids must agree; there is no implicit
interpolation. A zero initial magnetization has no normalized magnitude curve.
Compare absorbed energy separately; equal incident pulses need not deposit
equal energy. Canonical calc2026 is local OneDrive storage; a file write does
not establish cloud synchronization.

## Verification

The kernel tests cover all 36 partial derivatives, corrected correlation,
hydrogenic exchange, BR zero curvature, scalar NaN/Infinity rejection, and
U(1)/SU(2) invariance. Bridge tests cover orbital energy differences,
Hermiticity, complete band contraction, nonzero Bloch k and skew coordinate
transform, and 4-point MPI halos against a full-grid reference. The Taylor test checks actual hpsi integration, no-SOC spin/norm
conservation, midpoint timestep convergence and rebuilding the operator from
saved synthetic orbitals. This is not a full checkpoint-file round-trip or a
physical spin-wave/demagnetization benchmark. GNU Fortran tested; Fujitsu not
available. See the `testsuites/unit_nc_mgga23*` READMEs for commands.

NLCC-specific verification additionally checks explicit core-orbital phase
transformation, all three field derivatives, valence-only band subtraction,
core-plus-valence gauge covariance, actual Ni V2 read/grid/operator integration,
and macroscopic core current added exactly once. These tests pass in serial
and, where applicable, two-rank MPI. Fujitsu compilation remains untested.

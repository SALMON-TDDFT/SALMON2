# Optical observables and material reuse — 2026-11-20

This records the new feature request against05827a7e, separately from the prior
bug/operations review. It is an implementation queue, not a claim that the new
optical APIs already exist. The user's CI clarification takes precedence over
the request's suggestion to reproduce real physical observables in CI.

## Existing baseline

The prior review fixes include explicit native_host parser selection, verified
runtime loading and FDTD examples, PZ fixture preparation from existing GS,
one-time compatibility warnings and transverse-periodic FDTD. Synthetic
periodic dielectric/LD slabs have a local C-reference Fresnel comparison.
The new feature request does not require repeating that prior confirmation.

## CI and validation policy

Required automatic checks: inventories, Git-only preparation/build, package
integrity/imports, synthetic kernel tests and host/API contracts. No real-material
observable reproduction, long SCF/TDDFT runs, optical convergence campaigns or
experimental-spectrum matching in required CI. Existing NLCC, TDDFT-worker and
full-runtime campaigns are manual-only. Local physical validation remains
separately planned, explicitly executed and recorded, with inputs, tolerances,
runtime identity and limitations. CI configuration is not evidence of a run.

Do not acquire or copy new GS datasets for these tasks without the applicable
permission. Reuse existing authorized references. Import/checking never launches
a solver, download, build or job submission. Numerical work stays in Fortran/C;
Python owns inputs, orchestration and records. NumPy remains the core dependency.

## Implementation order

1. **Public TB-mBJ GS input.** Inspect native cval/default handling and expose
   supported parameters through immutable DFT settings and GS-to-RT identity.
   Verify Libxc-on/off support rather than assuming it. Separate band/response
   use from energy/force and ASE optimization claims; accepting an XC string
   alone does not validate those quantities. No long-range excitonic kernel.
2. **Frequency-domain observable contract.** Define optical_conductivity and
   dielectric_function above the existing spectrum conventions. Require explicit
   charge-current sign, driving field/gauge, units, window, dt, duration and
   amplitude; preserve available grid/k metadata and provenance. Decide tensor
   versus single-polarization semantics and the meaning of externally supplied
   versus total fields before coding the ratio. Frequencies with zero/insufficient
   drive do not yield invented numerical observables. Invalid observable requests
   must not stop an already runnable time-domain solver. Document scope explicitly.
   Excited-electron projections need the required reference/occupation data;
   absorbed-energy, HHG form and R/T normalization follow as separate observables.
3. **Cheap diagnostic records.** Linearity compares supplied or explicitly
   requested amplitude runs; never launches extra calculations implicitly.
   Sum-rule diagnostics record frequency range, missing-band limitations,
   current-operator convention and the applicable expected normalization before
   classifying agreement. They are not a universal operator-error detector or
   a substitute for user-owned k convergence. Attach reports to optical results.
4. **Native LD epsilon_infinity, before fitting.** Add a default1 model parameter,
   propagate it through configuration, material registration, native coefficients
   and checkpoint identities; preserve old default outputs. Implement the intended
   D=epsilon_infinity E+4pi P consistently in field updates, energy/observables
   and standalone response. Never substitute a fictitious high-frequency pole.
   CFL must use the fastest cell in the whole domain, including vacuum/PML and
   magnetic coefficients: a slower material alone does not relax a vacuum CFL.
   Keep parameter validity and stability constraints explicit and minimal.
5. **Constrained material fitting and provenance.** Fit sampled epsilon(omega)
   from either TDDFT observables or user-provided n,k data through a common entry.
   State time-factor/absorption sign conventions. Native numerical optimization
   must respect nonnegative strength/damping/resonance and report failures without
   emitting unstable coefficients. FitReport retains usable band, residual,
   pole count, timestep advice/stability margin and source-file SHA256/reference.
   Compare pole-count choices explicitly; do not fit or download automatically.
   Out-of-band excitation warns rather than silently claiming calibrated accuracy.
6. **Reuse existing trajectories.** Compare durations/windows/damping on the same
   data without extra solver runs. Derive post-hoc k-contribution diagnostics
   without reapplying occupation/k weights already present in decomposed currents.
   Add measured/unknown cost and memory advice calibrated to hardware and MPI
   layout, not an unqualified formula from k count alone. No automatic k scan.
7. **Geometry improvement, separately designed.** Use the completed periodic
   planar reference first. Review subpixel/interface tensor treatment and
   dispersive-pole/background separation, state/checkpoint representation and
   numerical energy/passivity implications. A tensor response may need more than
   scalar coefficient initialization; do not claim otherwise without inspection.
   Preserve the existing staircase default until an explicit new option is
   validated. Sphere remeasurement and monotonic/<0.5% claims require actual local
   measurements; they are not CI requirements or predeclared outcomes.

## Explicit exclusions and deferred decisions

- No replacement of the time-domain Yee core with FEM/DG.
- No plan_nstate_scan or automatic plan_kpoint_scan. Users choose k sampling and
  remain responsible for its convergence; the lack of a planner is not a proof
  that k convergence is unnecessary.
- No bundled Palik/Johnson–Christy tables. Follow the pseudopotential policy:
  references and hashes; downloads require user authorization and verification.
- A future external FEM/DG field adapter, excitonic XC kernels and decisions to
  return NLCC fixes upstream remain the project lead's decisions, not tasks
  authorized by this feature note.
- No guarantee for subnanometer gaps or sharp-corner local-response physics.
  Do not conflate geometry/discretization error with material-model validity.

- [x] 2026-11-20: public serial TB-mBJ DFT/cval, ASE structure bridge and material JSON; coefficient preserved through GS/TDDFT/archive and restart checks. Synthetic native XC probes agree bitwise with Libxc ON/OFF. No new GS collected; real-material precision and native parallel TB-mBJ SCF are not claimed. Optical observable work remains next.

- [x] 2026-11-20: single-polarization optical_conductivity/dielectric_function, explicit total-E/charge-current/Gaussian-AU contract, masked weak-drive/DC bins; native Fortran ratios and finite-band/linearity diagnostics on supplied data. Synthetic phase, sign, amplitude scaling and finite-band integral tests; no material accuracy certification or extra solver runs.
- [x] Native LD epsilon_infinity completed; fitting/provenance work is recorded below. HPC remains after the optical queue.

- [x] 2026-11-20: LD epsilon_infinity (default1), native continuous permittivity/D/background-energy and FDTD coefficients; JSON/provider/config/checkpoint identity propagation. Preserve pole P/J and default update order. Non-unit-background legacy impulse spectra warn/omit without stopping evolution. Source synthetic and native MPI1/2 Fresnel/default parity/direct checkpoint checks passed; real-material certification is not claimed.
- [x] Constrained LD fitting/provenance completed below; remaining optical work is duration/window/k post-analysis, then subpixel design. HPC remains deferred.

- [x] 2026-11-20: native passive LD fitting, direct epsilon/nk/masked-optical inputs, explicit pole counts and deterministic starts, optional initial model/background, residual/optimizer separation, data/file SHA256 and literature reference, calibrated band and warning-only excitation checks. Synthetic tests cover fit accuracy, passive failures, source identity and input independence. Local optimization and finite-band calibration do not certify a global optimum or real materials.
- [ ] Next optical work: duration/window sensitivity and post-hoc k contributions, then subpixel/interface design. HPC work remains later per user priority.

- [x] 2026-11-20 user explicitly authorized acquiring/retaining legally redistributable experimental optical data and fitting a reusable list. This supersedes the earlier exclusion of bundled Johnson-Christy numerical data for the selected CC0 database subset only. Preserve original measurements, license, pinned commit, SHA256, DOI and conditions; no journal articles/figures. Five materials,22 retained candidates; choose defaults using residual and time-resolution screening, retain warnings and optimizer nonconvergence. All experimental refits stay explicit/local, not required CI.

- [x] CC0 subset data library:5 sources,22 stored fit trials,4 screened defaults. Experimental comparison exposed Ag local weighted errors up to38.6% despite4.34% global RMS; mark needs_refinement, retain warning-only explicit reuse. Al stiff low-pole fits stay available but are not defaults. Preserve raw/license bytes across platforms with attributes.
- [ ] Improve Ag local-shape fitting/optimizer termination before presenting it as a default; do not relax/hide screening to obtain a pass. All new experimental-data matching/refits remain explicit local tools, not required CI.

## SciPy optical fitting (2026-11-20)
- Implemented optional optimizer="scipy" and fit extra; native objective/seed/statistics, SciPy bounded TRF and derivative estimates. No import-time install or fallback.
- Local CC0 study:44 SciPy trials plus22 retained native trials; selected SciPy fits locally terminated, Ag remains needs_refinement. Smallest screened defaults Au6/Cu6/Al4/Si6. Raw data preserved.
- Source and sdist-to-wheel installed validation:649 tests passed (one optional skip) in each; final artifact/record checks follow. No experimental matching in required CI.

## Optical oxides and tensor direction (2026-11-20)
- Added pinned CC0 scalar n,k films SiO2/Al2O3/TiO2/Ta2O5/HfO2 and22 SciPy trials (88 retained total); preserved prior five-material fit records. No n-only zero-absorption substitution. TiO2 residual passes but optimization-budget warnings remain.
- Common tensor representation in ALL modes supersedes3D-only scope. Implemented immutable measured/model(...,3,3), native isotropic embedding, compressed LDTensorModelConfig, FDTD/standalone/MPI-provider acceptance. Unknown directional response uses one independent scalar on the diagonal.
- Source tests651 pass, one optional skip. Audited unchanged FDTD binary:scalar/tensor inputs and observations bitwise equal for xyz polarizations, MPI1/2 (six paired cases). Wheel validation follows.
- REMAINING:independent diagonal axes/off-diagonal native propagation and fitting/rotations across modes; do not claim general anisotropy from isotropic tensors.

## Reciprocal tensor material response (2026-11-20)
- Implemented native PSD/SPD matrix validation, tensor epsilon, D/background energy and persistent coupled P/J updates. Isotropic kernel path retained. TensorLorentzDrude supports separate diagonal axes and reciprocal off-diagonal couplings; no automatic rotation or tensor fitting.
- Direct process-local snapshots restore without replay; overflow rejection preserves accepted P/J/E and permits retry. Independent states do not interfere.
- Local source suite653 tests passed, one optional skip. Actual MPI3 with two different tensor points:40 endpoint steps bitwise equal to independent native controls; caller MPI survives.
- REMAINING: general anisotropic native FDTD field/coefficient update with Yee-grid interpolation and checkpoint/schema wiring across field modes. General tensors are rejected by existing native FDTD rather than projected away. Source-to-sdist-to-wheel installation verified:653 tests passed, one optional skip; wheel MPI3/two-point controls and direct snapshot continuation match. Previous isotropic kernel responses match bitwise in three500-step controls.

## Native diagonal E/H and waveplates (2026-11-21)
- Principal-axis diagonal tensors connected to native E/LD coefficients, preserving P/J/Q ordering and existing history arrays. Explicit ld-diagonal input asset, hash binding, matching-runtime check and both batch/host-MPI staging. No new process-global material arrays.
- First prototype:old/new scalar, diagonal isotropic limit and MPI1/2 outputs bitwise equal. Dispersive distinct-axis controls match separate scalar native FDTD to below3e-18 absolute field difference.
- Synthetic quarter/half-wave plates:normalized |S3|0.99669 and S2=-0.999697 respectively; Stokes-norm discrepancy vs Fresnel slab0.00305/0.00563 at spacing8. This is not real-material or complete convergence validation. Refinement verified:error0.000760 at spacing4 (4.01-fold reduction). Git-only739-file export/build and source/installed-wheel655-test suites passed (one optional skip each). Native MPI2 direct checkpoint A/E/H and RT outputs match continuous execution bitwise, invalid/missing assets return1110 collectively. Final record:salmonpy-diagonal-waveplate-record-20261121 in the canonical results store.
- Native carrier/analytic Stokes diagnostic and explicit --run sample added. Off-diagonal Yee field coupling, automatic rotations and all-mode propagation coverage remain.


## LD-only reference-frequency interface smoothing (2026-11-21)
- [x] Dedicated native coefficient module, spatial tensor hookup and versioned
  assets, LD-only gates, electronic-point exclusion, local warning/fallback and
  direct checkpoint identity. Default off retains the staircase result.
- [x] Fixed-dt axial-cylinder grid refinement, homogeneous/legacy parity,
  rotated slab polarization, MPI2 direct restore/failed-trial recovery and
  source/wheel verification; see2026-11-21-reference-subpixel-smoothing.md.
  Claims are local/synthetic; no general broadband convergence guarantee.
- [x] Observation-only Fortran A4 low-frequency correction with raw/corrected
  spectra and explicit insulating assumptions; see2026-11-21-optical-a4-correction.md.
- [ ] Next: reuse trajectories for duration/window sensitivity and post-hoc
  k-contribution diagnostics. Earlier experimental Ag local-shape screening and
  the queued HPC receipt/launch/scale work remain distinct pending tasks.


## Stored-trajectory diagnostics (2026-11-21)
- [x] Native explicit common-frequency dt-sum transforms for uniform stored
  trajectories, none/Hann/causal cubic windows, prefix duration comparisons,
  Nyquist/weak-drive masking and descriptive relative RMS differences. Native
  k/band/spin reduction consumes already weighted currents, converts declared
  matter sign once, and optionally compares simultaneous physical charge current.
- [x] Retained current by k/total, optional charge-current spectra, signed
  cancellation diagnostics, immutable arrays and input content SHA256. No
  missing-entry inference, hidden resampling, automatic k scan or convergence
  certificate. Transform wall time/work items recorded; large cost warns only.
- [x] Paired OPTICAL_TRAJECTORIES guides, AI JSON and stored-NPZ sample added;
  preview starts no solver or output directory, --run is postprocessing only.
- [ ] Measured peak-memory advice for very large trajectory collections remains
  with the queued HPC observation-I/O work; current timing is measured local cost,
  not a site-scale resource guarantee. Ag local-fit screening remains separate.

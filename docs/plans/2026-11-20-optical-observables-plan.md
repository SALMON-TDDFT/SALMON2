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
- [ ] Next optical task: native LD epsilon_infinity, then constrained fitting/provenance, time/window/k post-processing, and subpixel design. HPC receipt work follows this optical queue per user steering.

- [x] 2026-11-20: LD epsilon_infinity (default1), native continuous permittivity/D/background-energy and FDTD coefficients; JSON/provider/config/checkpoint identity propagation. Preserve pole P/J and default update order. Non-unit-background legacy impulse spectra warn/omit without stopping evolution. Source synthetic and native MPI1/2 Fresnel/default parity/direct checkpoint checks passed; real-material certification is not claimed.
- [ ] Next: constrained LD fitting/provenance, then time/window/k post-analysis, followed by subpixel design; HPC queue remains deferred until those optical items.

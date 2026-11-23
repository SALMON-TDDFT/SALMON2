# FDTD/multiscale energy diagnostics (2026-11-21)

User request: absorbed_energy_density per macro point, R(omega)/T(omega), and
an FDTD energy conservation check. Implemented as explicit observation-only
functions; no solver update or automatic MPI reduction changes.

- Native Fortran trapezoidal total E dot charge j, retaining all three components,
  per-point signs, and finite checks with scalar IEEE calls.
- Native complex Poynting spectral flux ratios and raw unwindowed time fluences.
  Compiled NumPy FFT is the existing spectrum backend. Threshold masks weak input.
- Native weighted work, area conversion and signed total-budget residual.
  Field-energy change is explicit; net work already includes material storage.
- First scope is collocated samples and equal-area uniform plane ports. Users
  separate incident/reflected fields and account for lateral/PML flux separately.
  General nonuniform surface integration is not claimed.
- Energy mismatch warns (default 1% relative); normal calculations are not stopped.
  No physical-accuracy claim is inferred from a closed budget.
- Stored-NPZ example computes reports only, never starts a solver or MPI.

Verification: two-material points with known work distribution, R=.25/T=.5,
unequal outgoing impedance, vacuum R=0/T=1, changed field/other-outflow budget,
invalid input and direct native NaN/Inf rejection; chunk work sums match a
continuous observation trace; input arrays remain unchanged/outputs immutable.
Full source suite 674 tests with one optional skip (before distribution check).
This tests observation algebra/contracts, not a new full FDTD physical campaign.

Distribution verification: sdist→wheel→fresh environment passed; both source
and installed suites: 674 tests, 1 optional skip. Stored-array CLI reproduced
R=.25/T=.5 and a closed budget; intentional missing energy warned with exit0.
Build source/object audits passed. Results/logs are retained in canonical
calc2026/salmonpy-fdtd-energy-20261121; no GS or solver runs, no cloud-sync claim.

## Actual native coupling follow-up

Commit 58edd163 pushed; lightweight workflow44864514394 succeeded. Existing
audited executable SHA256 e6f75bd582a034bcaf9246f996477ca8f3184fa55d52ef795435822319427b61
was reused by SHA verification; no production-startup source/object reaudit.
The new C observation reference is separately source/object audited when built.

Actual LD slab wp20,w025,gamma2,epsilon_inf1: 256x8x8 grid, dx=.5, dt=.0015,
4000 steps, 40 material sample points. Signed relative budget residual
-1.0825977455e-4. Native final checkpoint matches reconstruction of LD j from
observed E (rtol1e-10/atol1e-18). Vacuum control residual -9.6777657479e-5.
Already-running dx=.25/dt=.00075 check: 512x8x8, 8000 steps, 80 sample points,
residual -2.7004059753e-5. Saved LD/field checkpoint data are not GS data.

Actual Si TDDFT uses existing 16^3/2^3-k/32-state GS in place, PZ, fixed ions,
2 macro points, macro dx20, 8 MPI ranks, dt=.003, 5334 steps. Delayed incidence
t1_start3,omega.5,tw8 leaves initial fields empty. Signed relative residual
-7.1097246880e-6; material work2.9678358087e-7 Hartree vs incident2.5712211069e-5.
This is non-negligible energy transfer (~1.15%), unlike the retained high-frequency
control. No new GS, occupations or wavefunctions were collected. Source-GS byte
hashes are verified unchanged by the runner.

The separate prefilled-wave low-frequency case closes at -1.1943651374e-4
only when initial field storage and half endpoint volumes are included. The early
analysis with incorrect endpoint weights remains as an error log, superseded by
corrected results; it is not accepted evidence. Finite-time truncated incident
traces explain apparent T>1 there, not gain. No clipping or threshold relaxation.

This is a local integration/conservation check, not real-material optical
convergence or multi-node scaling. User accepted the precision; no further
precision search. Reusable observation-only analyzer: tools/analyze_fdtd_energy.py,
with tests/native/test_energy_budget.c. Kernel calculations are Fortran/C; Python
reads, dispatches and records. Real physical calculations stay outside required CI.

Canonical records: calc2026/salmonpy-real-energy-20261121. No assertion of cloud sync.

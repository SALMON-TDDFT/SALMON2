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

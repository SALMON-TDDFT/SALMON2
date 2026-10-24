# Joint auxiliary-core NLCC for noncollinear MGGA23

User authorization: retain NLCC, implement and validate core kinetic density
and current consistently, and enable Ni SOC+NLCC LSDA/MGGA spin comparisons.
This extends the approved MGGA implementation; old files/results remain evidence.

## Chosen approximation

Use one explicit closed-shell auxiliary Pauli ensemble for all three frozen-core
fields N_c, T_c and J_c. The starting radial blocks are the generator's occupied
Dirac large/small components. Each block is treated as a Pauli angular shell,
with regular origin behavior and smooth matching at the core boundary, followed
by a common amplitude factor reproducing the existing KY partial-core density.
The retained ionic potential, valence projectors and partial-core density samples
must be numerically unchanged in the new V2 dataset.

This is an auxiliary-Pauli frozen-core approximation, NOT an exact Dirac/FW
reduction. Its positive-gradient T (no factor1/2) is the T in the 2023 functional;
it is not the Dirac kinetic-energy operator discussed in PRB114165112(2026).
Using tauW alone would eliminate multi-orbital core information; separately
smoothing arbitrary rho/tau tables would lose orbital representability. Neither
is used. A fully relativistic all-electron or PAW reconstruction is out of scope.

## Fields and variation

N=N_v+N_c, T=T_v+T_c, J_i=J_vi+J_ci; n gradients/Laplacians use the TOTAL density.
The core is fixed during orbital variation. The full XC band expectation is
contracted with VALENCE fields only, so the frozen background is not subtracted
as though it were an occupied propagated state. LSDA uses the same core density
and unchanged ionic potential; extra T/J metadata do not change LSDA's functional.

Reference closed-shell Pauli traces are N_c=(n_c,0,0,0), T_c=(tau_c,0,0,0),
J^0_ci=0, J^a_ci=C(r) epsilon_iab rhat_b. Include both relativistic radial blocks
with their own angular labels; partner spin currents need not cancel.

## Electromagnetic reference frame

Freeze the auxiliary core in the ionic rest/length-gauge frame. In the dipole
velocity-gauge representation use phase gradient q=-A/c:

J_ci(A)=J_ci(0)+q_i N_c,
T_c(A)=T_c(0)+2 sum_i q_i J_ci(0)+|q|^2 N_c.

Use the phase response at each predictor/corrector field and at endpoint
observables. The explicit core contribution dE_xc/d(A/c) is included ONCE in
the macroscopic current, in SALMON's positive velocity-current convention.
This term is not a separate occupation or arbitrary band-decomposed current.
Covariance tests co-transform BOTH core and valence; rotating valence alone
changes the relative state of a frozen background and is a different test.

## Data and implementation boundaries

- V1 remains supported for LSDA and preserved byte-for-byte as evidence.
- V2 adds a versioned joint auxiliary-core block with model/units/normalization,
  radial amplitudes and derivatives, and reconstructed n/tau/spin-current fields.
  Validate finite values elementwise, angular labels, density matching, and
  kinetic/current reconstruction. Missing/incomplete metadata fails MGGA+NLCC.
- Extend pp species tables and ppn grid fields explicitly; MPI-broadcast all
  metadata/tables. Reject mixed incomplete-core species, do not silently fallback.
- Interpolate core quantities on the same atom/image support; preserve the joint
  V2 tail rather than applying a scalar-density-only premature cutoff.
- Preserve current options and use a distinct model log; separate old validation
  records from this new core extension.

## Acceptance sequence

1. Radial/angular shell sum vs explicit spin spherical harmonics; regular origin,
   positive density/kinetic density, exact retained rho samples, convergence.
2. V1 regression; V2 valid/malformed/NaN/Infinity/mismatch reader tests; serial/MPI.
3. Core-on orbital energy finite differences, Hermiticity, valence-only band
   contraction, core-off regression; U(1)/SU(2) covariance with transformed core.
4. Finite-difference A derivative vs core current; midpoint/endpoint refresh;
   field-work consistency and no-field stability in synthetic propagation.
5. Read generated Ni V2 in SALMON, construct actual core/SOC projectors and
   matched LSDA/MGGA comparison inputs. Any nonconverged initial-state operator
   test is labeled as such; no physical demagnetization conclusion from it.
6. Real-space/k/time-step tests as feasible. No GS data collection under the
   current instruction. A physical comparison needs approved consistent initial
   states, pulse and material convergence; provide reproducible inputs without
   pretending these are already converged.

All new builds retain mandatory source/object PZHEEV audits. All Fortran IEEE
inquiries are scalar. Store results/records in canonical OneDrive calc2026,
without claiming confirmed cloud synchronization. No production job modified.

## Generator/consumer responsibility (user clarification)

The current generator defaults to the complete V2 core model. Selection belongs
to SALMON: LSDA consumes rho, while corrected noncollinear MGGA consumes the
coherent rho/tau/current set. Explicit V1 generation remains only a reproduction
compatibility mode. Historical source/data snapshots remain unchanged.

# DC force response gate

The approved BOMD design requires force/energy consistency before ionic motion.
The new explicit-force diagnostic passes one-fragment/full-buffer H4 limits but
fails that gate for genuinely truncated fragments. Two displacement sizes
confirm that the discrepancy is not the central-difference truncation error.
This is evidence against using the frozen-orbital force as the MD force.

## Potential and response variables

Collect independent real components of fragment orbitals, occupations, mixed
total density and chemical-potential variables into q. Include normalization
and gauge constraints in the converged DC equations C(q,R)=0. Choose the ionic
potential Phi explicitly: zero-temperature ground-state energy with fixed
integer occupations, or a defined finite-temperature free-energy functional.
The existing core-weighted TS diagnostic alone does not prove stationarity.

Along the self-consistent solution,

    dPhi/dR = Phi_R + Phi_q dq/dR,
    C_q dq/dR = -C_R.

An adjoint solve C_q^T lambda = Phi_q gives

    F = -Phi_R + lambda^T C_R.

The new diagnostic covers -Phi_R only for fixed orbitals/occupations and fixed
partition/support. Nuclear-position response of local/nonlocal projectors is
included; orbital/occupation response and moving partition membership are not.
The response operator must differentiate the actual fragment SCF and total-grid
rVV10 equations, not replace Phi_q by zero merely because fragment eigenstates
have converged. Gauge nullspaces and electron-count constraints must be handled.

## Measured discrimination

H4 at300 K,2 MPI ranks,density threshold1e-10,delta0.002/0.001 bohr:
full-cell buffers reproduce the energy derivative to about8e-6 eV/angstrom at
the smaller displacement. Buffer2 leaves~0.0504 eV/angstrom against E, and
~0.0496 eV/angstrom against the tested E-TS. Buffer1 leaves~0.0792 and
~0.0780 eV/angstrom respectively. These are combined electronic-response
residuals; they are not uniquely orbital response and are not a water benchmark.

## Pending choice and implementation

The user has been asked whether to prioritize fixed finite electronic temperature
with a free-energy formulation or strict zero-temperature integer occupations.
This selects Phi and the occupation constraints in C. No further permission is
needed to continue the approved force validation. Implementing the dependent
occupation/response model awaits that scientific choice.

Moving atom lists, Verlet integration, boundary-crossing tests and production
DC-MD remain unimplemented. The native DC-MD guard remains intentional.

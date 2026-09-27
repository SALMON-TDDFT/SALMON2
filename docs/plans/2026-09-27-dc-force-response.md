# DC force response gate

Scope update, 2026-09-28: this is the equilibrium SCF/BOMD response analysis.
The primary laser-excited route is now Ehrenfest dynamics, described in
2026-09-28-dc-ehrenfest-design.md. Its instantaneous-state force is not obtained
by automatically adding the BO adjoint correction derived here.

The approved BOMD design requires force/energy consistency before ionic motion.
The new explicit-force diagnostic passes one-fragment/full-buffer H4 limits but
fails that gate for genuinely truncated fragments. Two displacement sizes
confirm that the discrepancy is not the central-difference truncation error.
This is evidence against using the frozen-orbital force as the MD force.

## Potential and response variables

Collect independent real components of fragment orbitals, occupations, mixed
total density and chemical-potential variables into q. Include normalization
and gauge constraints in the converged DC equations C(q,R)=0. The selected ionic
potential is Phi=E-T*S at fixed finite electronic temperature and electron count.
Here T denotes kBT in Hartree; S is the dimensionless core-weighted Fermi entropy.
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

## Accepted choice and implementation

On 2026-09-28 the user selected fixed finite electronic temperature and free
energy. Occupations obey a common chemical potential and weighted total charge.
Both occupations and core weights contribute to the entropy derivative.
The occupation block is implemented first; the coupled orbital/density adjoint
and moving partition response remain required before MD can be certified.

Moving atom lists, Verlet integration, boundary-crossing tests and production
DC-MD remain unimplemented. The native DC-MD guard remains intentional.

## Fixed-temperature occupation block

`dc_thermal` uses f_i=1/(1+exp((epsilon_i-mu)/T)), weights
w_i=kweight_i*integral_core(|psi_i|^2), spin degeneracy g and
N=g*sum(w_i*f_i). Phi uses TS=-g*T*sum(w_i*[f_i log(f_i)+(1-f_i)log(1-f_i)]).
At fixed T,N, write b_i=f_i*(1-f_i). The implemented directional response is

    dmu = [sum(w*b*depsilon) - T*sum(f*dw)] / sum(w*b)
    df = b*(dmu-depsilon)/T
    dTS = g*sum(T*s(f)*dw + (epsilon-mu)*w*df).

Both terms in dTS are required. Tests separately perturb eigenvalues and core
weights, then both, and compare against independently reconverged occupations.
The response routine is not yet called by a nuclear-force calculation: obtaining
depsilon and dw requires the coupled SCF response described above.

Positive-temperature PBEh DC now uses a bounded, charge-checked chemical solve.
Other functionals and the zero-temperature legacy path are unchanged. A capacity
overshoot within 1e-10+64*epsilon*capacity is accepted only as normalization
roundoff; interior targets are solved normally. Numerically saturated states can
pass the charge solve, but response returns an explicit singular status when
sum(w*b)<=128*epsilon*sum(w). Occupied-only fixtures therefore do not certify
finite-temperature response for MD. No arbitrary zero response is substituted.

The stable entropy calculation retains minority tails even if f rounds to one.
The existing static diagnostic computes the same Fermi entropy from stored
occupations; it can lose sub-roundoff minority tails. This does not establish
stationarity or replace the full force gate.

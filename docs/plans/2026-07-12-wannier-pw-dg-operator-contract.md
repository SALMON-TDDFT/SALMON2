# Wannier+PW DG mathematical operator contract

This document freezes the operator used by the first global-ownership
Wannier+PW DG-DC and real-time milestone.  It is restricted to **gapped,
integer-occupation LDA** calculations.  It defines the physical discrete
space independently of its later MPI distribution.  The supported scope is gapped, integer-occupation LDA.

## Fragment domains, basis support, and ownership

The periodic simulation cell is the disjoint union of fragment core domains
`Omega = union_K Omega_K`.  Every real-space grid point belongs to exactly one
core.  This **unique grid ownership** is used for density, local-potential, and
energy quadrature; buffer and halo copies never contribute a second volume
weight.

A **fragment-restricted Wannier** basis function is the broken function

`W_(K mu)(r) = w_(K mu)(r)` for `r in Omega_K`, and zero on every other
element.  Its buffer samples are stencil and face-trace data only.  They do not
extend its volume support and do not create a second owner.

The plane-wave enrichment support is the periodic, windowed function

`P_(K G)(r) = chi_K(r) exp(i G.r) / sqrt(Omega_cell)`.

`chi_K` is supported on the fragment core plus its configured transition
buffer and is evaluated as one global periodic function, including on a
neighbor's core.  The normalized windows obey `sum_K chi_K(r)^2 = 1` at every
owned grid point.  This is the **plane-wave enrichment normalization**: the
underlying plane wave has cell normalization `1/sqrt(Omega_cell)`, whereas an
individual windowed column is not rescaled to unit norm.  Its actual norm and
all linear dependence are represented by the assembled overlap matrix `S`.
Renormalizing each windowed column separately is rejected because it destroys
the partition identity and changes WP/PP couplings.

Wannier gauges generated independently on adjacent fragments must be stitched
before any block below is assembled.  The DG buffer is part of trace
evaluation, not a physical overlap convention.

## Selected weak form

For each unordered face `F`, choose one canonical face normal `n_F` pointing
from the lexicographically smaller fragment id, denoted `K-`, to `K+`.  For a
scalar broken function,

`[[u]] = u^- - u^+`,

`{{grad u}}.n_F = 0.5 (grad u^- + grad u^+).n_F`.

Changing the canonical orientation changes both jump and normal signs and
therefore leaves the assembled bilinear form invariant.  Periodic neighbor
faces use the same rule after wrapping fragment ids; a face is assembled once.

The accepted zero-field Hamiltonian is the symmetric interior-penalty form

```text
a_DG(v,u) = sum_K [ 1/2 int_K grad(v)*.grad(u)
                    + int_K v* V_eff u
                    + <v|V_nl|u>_K ]
              - 1/2 sum_F int_F [[v*]] {{grad u}}.n_F
              - 1/2 sum_F int_F {{grad v*}}.n_F [[u]]
              +     sum_F int_F sigma_F [[v*]] [[u]].
```

The first implementation uses `sigma_F = 10 / h_F`, where `h_F` is the grid
spacing normal to the face.  This reproduces the existing self-face factors
`-1/4, -1/4, +sigma_F` and cross-face factors
`-1/4, +1/4, -sigma_F`.  A later polynomial-order-dependent penalty may be
introduced only by a new contract and convergence test; it is not an
environment-variable choice.

Local, Hartree, and LDA exchange-correlation potentials use unique core
quadrature.  The nonlocal pseudopotential is included consistently in the
volume operator and in the velocity commutator.  The fixed-basis DG-DC SCF
updates density-dependent potential matrices, but not `S`, kinetic, ionic, or
face matrices.

## Matrix block contract

Let `B = {W_(K mu)} union {P_(K G)}`.  Covariant matrices are
`S_ab=<B_a|B_b>` and `H_ab=a_DG(B_a,B_b)`.  The following blocks exist:

- `H_WW`: same-fragment volume WW blocks plus WW blocks on face neighbors from
  SIPG traces.  There is no dense all-fragment WW volume coupling.
- `H_WP` and `H_PW=H_WP^H`: every volume block whose Wannier core intersects a
  windowed PW support, plus the SIPG face terms obtained by inserting W and P
  into the same weak form.  Neighbor WP terms are physical support/face
  couplings, not communication artifacts.
- `H_PP`: every overlapping-support PP volume block.  Because each `P_(K G)`
  is a single periodic `H^1` function, its two traces agree and `[[P]]=0`;
  consequently PP jump/penalty contributions vanish.  No PP face penalty may
  be invented to compensate for distribution.
- `H_face`: the sparse union of all nonzero face-neighbor WW, WP, and PW
  contributions from the formula above.  It is a component of `H`, not a
  separate Hamiltonian and not an MPI halo term.

The support qualifications above apply to local volume and face terms, not to
the separable nonlocal pseudopotential.  For
`V_nl=sum_p |beta_p> D_p <beta_p|`, the Hamiltonian retains every nonzero
projector-mediated block

`<B_a|V_nl|B_b> = sum_p <B_a|beta_p> D_p <beta_p|B_b>`.

This rule is independent of direct support intersection and includes WW, WP,
PW, and PP pairs that are not face neighbors when both basis functions overlap
the same projector.  Nonlocal sparsity may be inferred only from vanishing
projector overlaps, never merely from basis-support or fragment adjacency.

`S` contains every support-overlap block and no face integral.  Its sparsity is
set by the actual Wannier and window supports; it is not restricted to immediate
face neighbors when a configured window reaches farther.  All retained mixed
basis columns must pass the configured positive-definiteness/conditioning gate
before diagonalization or propagation.

Global coefficient ownership in the first milestone means every participating
rank can address the same assembled coefficient vector.  Physical coupling is independent of MPI ownership.
A later distributed layout may change who
stores or communicates a block, but it must reproduce the identical `H` and
`S`; rank, halo, and reduction topology may not create or remove a physical
WW/WP/PW/face block.  Rank-local invalid traces, non-finite entries, or overlap
failures are diagnosed on the detecting rank before collective reduction.

## Position, periodic length gauge, and polarization

For a periodic cell, a globally single-valued linear `z` is not an observable.
The accepted periodic length-gauge operator starts from the covariant phase
matrix, evaluated by unique-core volume quadrature,

`U_S = <B|exp(i 2 pi z/L_z)|B>`.

Linear dependencies are removed before this construction by rebuilding the
retained mixed basis.  In that retained basis `S` must be positive definite.
X is square and full rank, with `X^H S X=I`; a rectangular orthogonalizer is
not admitted by this contract.  Form `U_tilde=X^H U_S X`, and replace numerical
non-unitarity caused by finite projection with its unitary polar factor `Q`:

`Q = U_tilde (U_tilde^H U_tilde)^(-1/2)`.

The minimum singular value of `U_tilde` must exceed the configured polar
tolerance.  A singular or ill-conditioned polar factor is fatal on the
detecting rank before collective reduction.

For a documented branch cut `b`, define the Hermitian projected position

`Z_tilde = (L_z / 2 pi) (-i Log_b Q)`,

`Z_S = S X Z_tilde X^H S`.

For the required square full-rank `X`, this is equivalently
`Z_S = X^-H Z_tilde X^-1`; the first form is normative because it does not
expose explicit inverses.  There is no discarded subspace after the retained
basis has been rebuilt.

Here `Log_b` is the spectral matrix logarithm whose eigenphases are followed
continuously from the zero-field reference.  The branch is admissible only
while the spectrum of `Q` remains separated from its cut by the configured
tolerance; loss of that separation is fatal rather than a silent phase jump.
The same S-orthogonalization convention must be used by initialization,
propagation, and observables.  This construction, rather than a pointwise
sawtooth matrix, is the sole `Z_S` used below.

The polarization branch is selected by continuity from the zero-field state
and unwrapped in time; `Pz(t)` is the primary observable and
`Jz(t)=dPz/dt`.  A raw sawtooth coordinate with an unreported cut, independent
per-fragment branch choices, or a direct elementwise matrix logarithm is
rejected.  Under a common origin shift that does not cross the tracked cut,
`U_S` gains the scalar phase `exp(i 2 pi z0/L_z)` and the continued logarithm
gives `Z_S -> Z_S + z0 S` modulo one polarization quantum.  The commutator and
unwrapped response are therefore origin invariant.

All WW, WP, PW, and overlapping-support entries of `U_S`, and hence `Z_S`,
come from unique core quadrature.  There is no independent interface position correction.
In particular, adding a face matrix merely because `H_DG` contains face terms
is rejected.

The field coupling is `H(t)=H_DG[n(t)] + E_z(t) Z_S` using the same branch and
S metric as initialization, propagation, and observables.  Scientific controls belong in the namelist;
environment variables may enable diagnostics only and
must not select the PW window, penalty, branch, operator blocks, or field
coupling.

## Discrete velocity relation and interface decision

In a non-orthogonal basis, a covariant operator matrix `O_S` acts on
coefficients as `S^-1 O_S`.  Define the covariant commutator velocity by

`(V_z)_S = i (H S^-1 Z_S - Z_S S^-1 H)`.

`V_z = i S^-1 (H S^-1 Z_S - Z_S S^-1 H)`.

Thus `v_z = i s^-1 (h s^-1 z_s - z_s s^-1 h)` is the required
coefficient-space formula.  Omitting either metric inverse is rejected.

To see why no separate face bilinear form belongs to position, first use an
ordinary local coordinate `z` on an internal face of the SIPG mesh.  It is
continuous there, so `[[z]]=0`, `[[z u]]=z[[u]]`, and
`grad(z u)=e_z u+z grad(u)`.  The product-rule pieces survive in the
Hamiltonian-position commutator as the DG-consistent face contribution to
velocity.  In the periodic operator, the multiplicative phase
`exp(i 2 pi z/L_z)` is itself continuous across the wrapped boundary and its
matrix is assembled entirely by volume quadrature.  Applying the global polar
factor and spectral logarithm does not create an element-face bilinear form.
DG face effects therefore enter through `H` in the S-metric commutator, not
through an independently postulated correction to `Z_S`.

The periodic `Z_S` is the phase/polar/log operator defined above; it is not
identified with the finite-space projection of a raw sawtooth coordinate.
Consequently, the numerical relation is checked directly in the observable
form.  A deterministic generalized model propagates the same initial state
with the complete `H_DG=H_volume+H_face`, obtains `d<Pz>/dt` by a centered time
difference, and compares it with the expectation of
`i(H_DG S^-1 Z_S-Z_S S^-1 H_DG)`.  The finite-difference dPz/dt residual must
converge at the expected time-difference order.

The same fixture also evaluates a velocity in which `H_face` is deliberately
omitted from the commutator while propagation retains it.  Its face-omitted velocity residual must be nonzero.
This demonstrates that DG face physics
enters velocity through `H_DG`; it does not justify or require a separate face
bilinear form in `Z_S`.  A direct-coordinate product-rule test is rejected for
this periodic contract because projection, polar unitarization, and the matrix
logarithm do not commute in a finite mixed basis.

The rejected alternative is a nonzero interface position correction added
without a weak-form derivation and a small-model identity test.  Such a term
would double count the product-rule face contribution and generally violate
Hermiticity/origin covariance.  A future nonzero interface position correction
would require changing the discretization itself and demonstrating both the
weak identity and the numerical commutator residual.

For every eigenpair `H c_n = S c_n epsilon_n`, the diagnostic identity is

`c_m^H (V_z)_S c_n = i (epsilon_m-epsilon_n) c_m^H Z_S c_n`.

The accepted initial occupied states are obtained from the self-consistent,
complete zero-field DG generalized eigenproblem `H_DG C = S C epsilon`, with
all volume, nonlocal, WP/PW, flux, and penalty blocks present.  Projected
conventional orbitals and eigenstates of a Hamiltonian with DG terms omitted
remain inadmissible.

## Scope and deferred implementation

This Task 1 contract and its textual regression test do not implement the
shared algebra module, DG-DC solver, exponential propagator, or new namelist
variables.  Those belong to Task 2 and later.  Any later implementation that
cannot realize this exact operator must stop and revise the contract before
changing production code.

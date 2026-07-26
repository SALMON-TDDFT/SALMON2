# Core-Centered Buffer-Supported Wannier Galerkin Design

## Purpose

Build a default-off Gamma, non-SOI, LDA route that constructs localized
Wannier functions in fragment-local periodic boxes, retains their tails
through the configured buffer, and solves the global ground state and
real-time dynamics in the resulting overlapping nonorthogonal basis.

The production route does not add a direct real-space SIPG face correction
to the existing buffer-extended DC kinetic stencil. The direct-DG work
through commit `05e63f7` and the later uncommitted diagnostics remain
investigation evidence only.

The route must not promote itself to LCFO, EigenExa, WPW, normal checkpoint
publication, or conventional RT. The established non-DG DC route, including
LCFO and EigenExa, remains unchanged.

## Mathematical space

Let the physical grid be partitioned into uniquely owned, nonoverlapping
fragment cores:

\[
\Omega = \bigcup_F \Omega_F,
\qquad
\Omega_F\cap\Omega_G=\varnothing\quad(F\ne G).
\]

Each Wannier function \(w_{\alpha F}\) has one center-owner fragment \(F\),
but its retained support extends through that fragment's buffer:

\[
\operatorname{supp}w_{\alpha F}
\subset \Omega_F^{\rm core+buffer}.
\]

The owner stores the only authoritative copy of the function and its
derivatives. Overlapping fragment boxes must not define competing copies of
the same Wannier function.

The finite basis is generally nonorthogonal under the global physical
metric. Its overlap matrix is

\[
S_{ij}
=
\sum_K\int_{\Omega_K}w_i^*(\mathbf r)w_j(\mathbf r)\,d\mathbf r.
\]

Every physical quadrature point is counted once, while every retained
Wannier tail covering that point participates. The basis is accepted only
when \(S\) is Hermitian positive definite on the retained subspace and its
condition number satisfies the configured gate.

## Candidate construction and localization

The conventional fragment DC solve supplies occupied and empty candidate
orbitals in each buffer-extended periodic fragment box. Candidate
construction is complete before localization and does not call LCFO,
EigenExa, or WPW.

Localization selects center-owned Wannier functions while preserving their
complete retained buffer tails. The candidate-window and target ranks are
recorded. The retained space must contain the requested occupied subspace
under the global core quadrature metric.

The buffer exterior is an auxiliary numerical boundary, not a physical DG
interface. The value and gradient norms of every retained Wannier function
at that boundary must be below their tolerances. Increasing the buffer must
converge the overlap, Hamiltonian, kinetic, velocity, and density matrices.
Abrupt truncation at a core or buffer boundary is forbidden.

## Unique-core quadrature

For a local or differential one-body operator \(\hat O\), matrix elements
are assembled as

\[
O_{ij}
=
\sum_K\int_{\Omega_K}
w_i^*(\mathbf r)\hat O w_j(\mathbf r)\,d\mathbf r.
\]

The integration domain is the unique core partition. Buffer data supplies
the values and stencil halo of any owner-centered Wannier tail that reaches
the owned point. Integrating every fragment buffer independently is
forbidden because it double counts physical volume.

The distributed assembly records the contributing Wannier owner, the
physical core-point owner, the source buffer fingerprint, and the operator
generation. Missing tails, duplicate point ownership, stale generations,
and incomplete owner pairs are collective failures.

## Kinetic energy

The production kinetic matrix is constructed as one symmetric weak
bilinear form:

\[
T_{ij}
=
\frac12\sum_K\int_{\Omega_K}
\nabla w_i^*(\mathbf r)\cdot\nabla w_j(\mathbf r)\,d\mathbf r.
\]

Core-boundary derivatives use retained buffer tails. This is not combined
with the direct real-space SIPG correction previously added to the
buffer-extended DC strong stencil.

If a discontinuity remains in the assembled overlapping basis, its
interface contribution must be derived as part of the same symmetric weak
form and accumulated once by a canonical face owner into both conjugate
blocks. An independent additive face Hamiltonian is forbidden.

The kinetic gate requires

\[
\|T-T^\dagger\|\le\epsilon_T
\]

and agreement between the weak kinetic energy and an independently
evaluated reference on accepted small systems.

## Local and nonlocal Hamiltonian

The density-dependent local Hamiltonian uses the existing DC ionic,
Hartree, XC, and local-potential implementations. At each outer Wannier SCF
iteration, these fields are evaluated on uniquely owned physical core
points and projected with the same global quadrature used for \(S\).

Nonlocal pseudopotential projectors are assigned a unique atom owner.
Wannier-projector overlaps include every retained tail that reaches the
projector support, but each atom/projector contribution is accumulated once.
The resulting nonlocal matrix must be Hermitian independently of the total
Hamiltonian.

The complete matrix is

\[
H[n]=T+V_{\rm ion}^{\rm local}+V_{\rm NL}
     +V_H[n]+V_{\rm xc}[n].
\]

It must satisfy \(H=H^\dagger\) before any eigensolve or time step.

## Density and electron count

For density matrix \(P\), the density on the core owned by \(K\) is

\[
\rho(\mathbf r)
=
\sum_{ij}w_i(\mathbf r)P_{ij}w_j^*(\mathbf r),
\qquad \mathbf r\in\Omega_K.
\]

The sum includes all Wannier tails covering the point, not only functions
centered in \(K\). No density is independently published from buffer
points.

The electron-count identities

\[
N_e=\operatorname{Tr}(PS)
\]

and

\[
N_e=\sum_K\int_{\Omega_K}\rho(\mathbf r)\,d\mathbf r
\]

must agree within the configured tolerance.

## Ground-state solve

The fixed-basis inner problem is the generalized Hermitian eigenproblem

\[
H[n]C=SC\varepsilon.
\]

It is solved without EigenExa. A local dense solver is permitted only for
focused fixtures and sizes covered by an explicit memory gate. Production
uses the planned distributed orthogonalized coefficient solver.

The outer loop:

1. assembles \(S\) once for a fixed Wannier generation;
2. assembles \(H[n]\);
3. solves the occupied and required empty coefficient subspace;
4. constructs the density on uniquely owned cores;
5. updates the existing DC Hartree, XC, and local potentials;
6. mixes density using preserved, seeded history;
7. repeats until the unmixed density and coefficient residuals converge.

Changing the candidate window, buffer, or Wannier functions starts a new
explicit outer-basis generation. It is forbidden inside a coefficient solve.

## Position, momentum, velocity, and transitions

Position matrix elements use unique-core quadrature and all retained tails.
For periodic observables, Berry-compatible or cell-periodic forms are used
instead of an unqualified unbounded position operator.

The canonical momentum diagnostic is

\[
p_{ij}
=
-i\sum_K\int_{\Omega_K}w_i^*\nabla w_j.
\]

The physical velocity matrix is derived consistently from the complete
discrete Hamiltonian:

\[
v=i[H,r].
\]

It includes the nonlocal-pseudopotential commutator. Optical transition and
current observables must not silently substitute canonical momentum for
velocity.

The operator gate checks Hermiticity or anti-Hermiticity as appropriate,
gauge covariance under retained-space rotations, origin/cell conventions,
and agreement with a direct small-system reference.

## Real-time propagation

After the ground-state and fixed-basis gates pass, real time evolves only
the coefficients:

\[
iS\dot C(t)=H(t)C(t).
\]

The Wannier functions, overlap matrix, ownership, and buffer fingerprints
remain fixed. The propagator must preserve \(C^\dagger SC\) in the field-free
case and reproduce field-free energy conservation within tolerance.

Tail escape or loss of fixed-space validity is an explicit failure. There is
no automatic basis update and no fallback to conventional RT.

## Gates

### Algebra and ownership

- unique physical core-point ownership;
- unique Wannier center ownership;
- complete tail coverage for every contributing owner pair;
- \(S=S^\dagger\), positive retained metric, acceptable condition number;
- invariance to sign, permutation, and retained-subspace rotation.

### Basis and buffer

- candidate and target rank policy;
- occupied-subspace inclusion;
- centers and spreads;
- buffer-boundary value and gradient norms;
- convergence of \(S,H,T,v,\rho\) across at least two buffer widths.

### Operators

- \(T=T^\dagger\);
- \(V_{\rm NL}=V_{\rm NL}^\dagger\);
- \(H=H^\dagger\);
- velocity/nonlocal-commutator consistency;
- small-system reference agreement.

### Ground state

- generalized eigen-residual tolerance;
- \(C^\dagger SC=I\);
- unmixed density fixed point;
- \(\operatorname{Tr}(PS)\), integrated charge, and target electron count;
- energy and density agreement across accepted windows and buffers.

### Real time

- metric norm conservation;
- field-free energy conservation;
- short-run and restart-split equivalence;
- explicit rejection of tail escape or stale basis generations.

### Si64 production matrix

The Gamma, non-SOI, PZ Si64 gate uses 8 MPI ranks and 1 OpenMP thread. It
compares at least:

- two equivalent fragment decompositions;
- two buffer widths;
- two candidate/target Wannier windows.

It records runtime, peak memory, ranks, metric spectrum, centers, spreads,
tail norms, operator defects, density and coefficient residuals, charge,
energy, and restart state. Wannier coefficient RT remains blocked unless
one full matrix configuration passes the ground-state gate.

## Route isolation

With the new route enabled:

- do not call the direct real-space DG continuation;
- do not call LCFO, EigenExa, WPW, normal checkpoint publication, or
  conventional RT;
- publish only the route-specific overlapping-Wannier checkpoint after the
  complete ground-state gate.

With the route disabled, established DC SCF, LCFO, EigenExa, checkpoint, and
RT behavior is unchanged.

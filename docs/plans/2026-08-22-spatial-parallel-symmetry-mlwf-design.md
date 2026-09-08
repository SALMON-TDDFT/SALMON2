# Spatially Parallel Symmetry-Constrained MLWF Design

## Scope and decision

This is the production localization design for Wannier-representable insulating and semiconducting DC+LCFO DG calculations from roughly 1,000 toward 100,000 atoms, primarily interfaces, surfaces, liquids, and other large Gamma-point supercells without an internal primitive translation lattice. Metals are out of scope because the present DC method does not provide their ground state. The first version does not classify topology; representation, decay, or localization feasibility failure, including a possible obstruction, fails the complete block closed to WPW.

Use spatially owned local Jacobi localization, with every update constrained to a complete verified finite-symmetry pair orbit. The result is a stationary local minimum of the stated Gamma-supercell periodic spread on the symmetry-constrained, fixed-occupation unitary manifold; symmetry constraints and finite initialization do not imply the global minimum. Retain well-localized complete symmetry blocks as orthogonal Wannier functions and replace rejected blocks with the existing windowed-plane-wave (WPW) complement. Wannier90 remains a small-system development/CI reference oracle, not a runtime, per-fragment, or production dependency. The production executable neither calls Wannier90 nor gathers distributed orbitals/DMN data for it.

The design combines the robust pair-Jacobi idea used in linear-scaling MLWF work and CP2K with ONETEP/BigDFT-style localized ownership, bounded halos, and sparse operator application. It does not adopt nonorthogonal NGWFs as the final RT basis.

## Established Si64 reference baseline (2026-08-24)

Before replacing Wannier90 with the spatially parallel backend, the corrected-unit Wannier90 plus WPW Hybrid-DG route was run through the complete Si64 ground-state handoff with 8 MPI ranks, `OMP_NUM_THREADS=1`, and no time cutoff. This is reference evidence for the handoff and DG-SCF layers; it is not evidence that the spatially parallel localization backend is implemented or weak-scaling compliant.

The verified affine group had order 1536, factored into a translation subgroup of order 32 and a point co-group of order 48, with five generators. The corrected-unit Wannier90 result reported gauge spread about 2205.6 and total spread about 3368.8 in its output convention. Its fixed-center covariance defect was `2.38e-13`; the retained Wannier center action was nonmonomial for at least operation 2. Therefore the reference demonstrates a symmetry-covariant subspace and a usable Hybrid-DG handoff, not a set of uniformly fragment-local, symmetry-permuted individual Wannier functions. This distinction remains a primary motivation for the spatially local backend and complete-block WPW fallback.

After correcting the pullback representation traversal described below, the Hybrid-DG SCF converged in 66 iterations. Final residuals were density `2.50e-9`, band-energy change `8.89e-8`, generalized eigensystem `1.75e-15`, electron-count error `4.32e-12`, and symmetry `2.28e-14`. The converged occupied checkpoint was published. During the converged rerun the eight ranks remained near full CPU utilization; observed resident memory was about 1.37--1.44 GiB per rank in the SCF phase, while the monitored peak over the earlier complete diagnostic run was about 1.80 GiB per rank. These are Si64 reference-route measurements, not asymptotic memory bounds.

The pre-handoff symmetry construction was materially more expensive than the final 66 DG-SCF iterations. In particular, the full affine order 1536 makes any design retaining an `N_basis x N_basis x |G|` representation array unacceptable. The production spatial backend must retain generators, factored subgroup/co-group metadata, sparse transporters, and bounded probes only.

## Required invariants

1. Rotations never mix occupation blocks with unequal occupations. Occupied and empty sectors are localized separately.
2. Localization updates are unitary before controlled sparsification. Fixed-accuracy truncation is followed by a bounded sparse metric/polar correction, and all induced defects are included in sweep acceptance. No global dense reorthogonalization repairs drift.
3. The independent optimization variable is a complete symmetry-connected generator block, not an arbitrary pair.
4. Every provisional unitary candidate decreases the explicitly defined Gamma-supercell periodic spread within tolerance. A sweep is accepted only when the complete objective after truncation and metric correction does not exceed its checkpoint value within tolerance, and final convergence bounds the gradient over retained and omitted admissible directions.
5. Wannier acceptance/rejection and WPW augmentation act on complete symmetry/irrep blocks.
6. No rank materializes all orbitals, the full field, a dense global gauge, or a global DMN matrix.
7. Canonical identifiers, coloring, and tie-breaking make invariant results MPI-layout independent within reduction tolerance.
8. Disjoint fragment cores contain the only authoritative real-space field values. Fragment buffers are generation-matched read-only caches used to evaluate operators on core output points; they are never independent physical copies or integration domains.

## Where symmetry enters

### Geometric products and pullback column actions

The product table stores the geometric left product. Real-space target maps implement the pullback

    (T_g psi)(r) = psi(g^-1 r),

and consequently compose in reversed table order when applied successively. If basis column actions are published by

    T_g B = B D_g,

then the implementation convention is an anti-representation of the geometric table,

    D_(g h) = D_h D_g.

This is not an optional gauge choice. A breadth-first traversal that stores `D(parent) D(generator)` must advance the geometric label on the left,

    operation = product_table(generator, parent),

so that `D(generator * parent)=D(parent)D(generator)`. Advancing with `product_table(parent,generator)` is wrong for these pullbacks. It passed earlier C2 tests because an abelian involution cannot expose multiplication reversal, but on Si64 it left post-average H and density covariance defects of order `1e-4`.

Every implementation must test the convention with a noncommuting group and at least one redundant generator whose table relation is order-sensitive. The reference fixture uses the right-regular anti-representation of `D3=C3 semidirect C2`, compares the factored result with a direct full-group Reynolds average, and separately rejects unitary generators that violate the supplied product table. A trivial or C2-only fixture is insufficient. Because the translation subgroup is normal, both left and right cosets cover the full group; agreement or disagreement of a factorized average alone does not identify the representation convention. The decisive checks are the multiword transporter relation and the direct full-group oracle on a noncommuting fixture.

For large groups, relation validation uses a small fixed set of dense deterministic row probes propagated over the generator tree. Its storage is `O(k N_basis |G|)` and its work is `O(k N_basis^2 |G| n_generator)` for bounded `k`; it must never materialize all full representation matrices. A failed relation check is fatal before Hamiltonian averaging and reports the left operation, generator operation, and defect.

### Initial basis

The coherent DC+LCFO target space is divided first by occupation, then into complete center orbits and small representation blocks. Only finite operations verified against the actual atomic structure, production grid, LCFO basis, H, S, and target projector are imposed. Supercell periodicity enters distances, moments, support images, and halo exchange; it does not create an artificial internal primitive translation orbit. Symmetry mixing unequal occupation blocks is a construction error.

Before localization, a separate distributed seed stage constructs local columns from the sparse generalized LCFO pair (H,S). It first verifies that S is Hermitian positive definite within tolerance, builds a sparse approximate orthogonalizer X by a scaled metric Newton-Schulz/purification iteration, and verifies ||X^dagger S X-I||. It then applies a bounded-degree Chebyshev filter to H_bar=X^dagger H X using matrix-free sparse H/S actions and verified spectral bounds. If Q contains the filtered, rank-revealed columns in the orthogonalized representation, the LCFO-basis columns are C=XQ and their S-orthogonal projector is P_C=C(C^dagger S C)^-1 C^dagger S. Coefficient-space seed idempotency, mutual occupied/empty orthogonality, subspace residuals, and principal angles use this S-metric definition and a stated induced norm. After stitching, localization uses the disjoint-core physical inner product; the later hybrid P_target comparison uses the production PoU metric only after its bounded-block equivalence to that core representation has passed. No dense eigenvector matrix is materialized. Empty response seeds use a smooth upper-energy filter followed by S-orthogonal projection out of the occupied filtered space. Symmetry-block pivoted rank revelation defines the retained target columns and P_target. Gap, metric condition, spectral bounds, filter error, fill, iterations, residual, projector idempotency, and memory are explicit feasibility gates; failure sends the complete requested response block to WPW rather than invoking dense diagonalization.

After rank revelation and S-orthonormalization, C^dagger S C=I within tolerance. With the Ritz matrix Lambda=C^dagger H C and generalized residual R_H=HC-SC Lambda, production gates include

    ||C^dagger S C-I||_2,
    ||P_C^2-P_C||_S,
    ||P_C^dagger S-S P_C||_2,
    ||R_H||_{S^-1}/max(1,||HC||_{S^-1}),
    ||P_occ P_empty||_S.

Let A be the linear LCFO-to-authoritative-core stitching map, including the existing partition-of-unity composition but excluding read-only buffer duplication. The coefficient metric used by localization must be the pullback of the physical core inner product on every retained target block. Before any spread evaluation or update, require

    W = A C,
    S_map = A^dagger A,
    epsilon_map = ||C^dagger(S-S_map)C||_2 <= tolerance,
    ||A(CU)-(AC)U||_(core,F) <= tolerance

for dense fixture unitaries at setup. Also require symmetry equivariance ||A(B_g C)-T_g A(C)||_(core,F) within tolerance. Production evaluates epsilon_map from the two small retained-block Gram matrices C^dagger S C and W_core^dagger W_core; it never forms global A or S_map. After these gates pass, C and the map fingerprint become immutable seed provenance and W=AC becomes the sole mutable localization state. The map is not reapplied during a sweep. These gates justify initial physical orthonormality and the subsequent real-space identity Z'_I=U^dagger Z_I U. A failed map/metric/equivariance gate is a construction failure; localization cannot repair it by polar correction.

Here ||x||_S^2=x^dagger Sx, ||Y||_(S,F)^2=tr(Y^dagger S Y), ||y||_{S^-1}^2=y^dagger S^-1 y is evaluated by a certified sparse metric solve rather than a formed inverse, and ||A||_S=||S^(1/2) A S^(-1/2)||_2. These coefficient-space definitions are used only by the seed, generalized projector, Ritz, and their regression gates. After stitching, localization uses the core physical inner product and hybrid-span validation uses the separately defined production PoU physical-field inner product; neither applies coefficient-space S to a physical field. Dense seed references additionally compare principal angles after mapping both coefficient spaces through a verified S^(1/2) action.

The empty sector is finite and input-defined by an energy/rank density target plus a smooth filter width. Degenerate shells and symmetry irreps are never partially selected. The resulting number of response states per atom and filter degree must remain bounded. The constructed subspaces satisfy P_occ P_empty = 0 within tolerance, and occupied/empty rotations remain separate. A sharp empty-window projector is permitted only when an actual upper spectral gap is verified.

### Gamma-supercell localization functional

Let {G_I,w_I}, I=1,...,N_G, be a fixed inversion-closed shell of nonzero supercell reciprocal vectors and positive weights. They are accepted only when they reproduce the Cartesian metric,

    A_G = sum_I w_I G_I G_I^T,
    ||A_G-I_3||_2 <= epsilon_metric,

with both +G and -G carrying half of a pair's total weight when both are stored. For an S-orthonormal real-space Wannier field w_n, define

    Z_I,mn = <w_m|exp(i G_I.r)|w_n>,
    z_I,n = Z_I,nn,
    Omega_Gamma(W) = sum_n sum_I w_I (1-|z_I,n|^2).

Write the supercell lattice matrix as L and its reciprocal matrix as B=2 pi L^-T. Each shell vector is stored canonically by an integer Miller column q_I in Z^3 and evaluated as G_I=B q_I. For every verified finite spatial operation g=(R_g|t_g), first obtain its exact integer unimodular fractional-lattice action A_g from R_g L=L A_g using the same crystallographic integer data as the grid permutation. Construct the reciprocal integer automorphism as

    M_g = A_g^T,
    B M_g = R_g^T B.

M_g and its products are evaluated by integer arithmetic. The Cartesian equality is a residual check against L/R_g, never the source of a rounded Miller map. The shell and weights must be closed under this action through a canonical permutation pi_g satisfying

    q_pi_g(I) = M_g q_I,
    R_g^T G_I = G_pi_g(I),
    w_pi_g(I) = w_I,
    M_gh = M_h M_g,
    pi_h o pi_g = pi_gh  when T_g T_h=T_gh.

The field action convention is fixed as

    (T_g psi)(r) = psi(R_g^-1(r-t_g)).

With this convention and T_g W=W D_g, change of variables gives the moment covariance

    D_g^dagger Z_I D_g
      = exp(i G_I.t_g) Z_pi_g(I),

including the stored periodic-image shift in t_g. The sign of the phase is tested directly against the discrete grid permutation; changing the field-action convention requires changing both the formula and its tests. At Gamma a full supercell lattice shift has unit phase, but it remains part of the integer-grid transporter and its fingerprint. This covariance implies Omega_Gamma(T_g W)=Omega_Gamma(W). Metric completeness without this group closure is insufficient and fails the requested symmetry route.

The weights have dimensions of length squared, so Omega_Gamma has the dimensions and large-cell limit of the quadratic spread. This is the sole production objective. The alternative logarithmic localization tensor -sum_I w_I log|z_I,n|^2 is reported only as a diagnostic and is never mixed into acceptance. A shell with singular/ill-conditioned A_G or disagreement with the dense complete-field evaluator fails closed. A moment |z_I,n| below the configured reliability floor leaves that center component undefined and forces anchor-based ownership during iteration; it does not invalidate the polynomial objective, but it must be resolved before final Wannier acceptance.

### Units and discrete quadrature contract

All lengths are converted once to the production internal length unit before L, B, G_I, grid coordinates, or Delta V are formed. For a uniform periodic N_1 by N_2 by N_3 supercell grid with no duplicated endpoint,

    Delta V = |det L|/(N_1 N_2 N_3),
    sum_{all unique grid points} Delta V = |det L|.

If the production grid supplies an equivalent canonical quadrature weight, that single value is used instead; no caller multiplies by the cell volume or Delta V again. Input LCFO/grid arrays may use an implementation-specific scaled storage convention, so the stitching map A contains the sole conversion into authoritative physical-sample fields W. Setup fingerprints the source convention and proves that A maps both an analytic constant and a normalized fixture to physical samples. After this boundary, core W and every buffer copy store the same physical sample with dimensions length^-3/2; no sqrt(Delta V)-scaled alternative is permitted inside localization. This makes the dimensional ledger unambiguous:

    r,L                         length
    G_I                         length^-1
    Delta V                     length^3
    w_n(r), p_a,mu(r)           length^-3/2
    <w_m|w_n>, Z_I,mn, q_a,mu   dimensionless
    w_I                         length^2
    Omega_Gamma, raw A          length^2
    N_W90(Pi_G(A)), K, U        dimensionless.

Here projector normalization is understood in the convention used by the pseudopotential implementation; if a stored projector carries a different internal prefactor, the q dimension and coefficient path are recorded explicitly rather than forced into this shorthand. H matrix elements have the production energy unit and S/Gram matrices are dimensionless. "Delta V exactly once" refers to the mathematical quadrature represented by each integrated quantity. Because localization W is explicitly physical-sample storage, that quadrature is an explicit single Delta V multiplication there. A source array that already absorbs sqrt(Delta V) is converted once by A and is never multiplied under its old convention inside localization. A finite-difference operator application contains its own derivative coefficients but receives no extra Delta V; Delta V enters only when its output is integrated into a matrix element.

Quadrature weights, PoU weights, reciprocal-shell weights, and symmetry multiplicities are distinct quantities and use distinct variables. PoU weights compose overlapping fragment contributions before the unique core field exists; they are never substituted for Delta V or reapplied to a stitched field. Buffer copies and periodic endpoint images contribute zero additional integration weight. If both +G and -G are stored, each receives half the intended pair weight exactly once. Setup gates include constant-field volume, normalized-orbital norm, identity overlap, direct-versus-core moment, nonlocal projector overlap, +G/-G pair weight, and unit-rescaling fixtures. Each fixture deliberately injects one missing or doubled Delta V/PoU/shell factor and must fail.

Wannier centers are obtained from phi_I,n=Arg(z_I,n) by solving the wrapped weighted system

    r_n = argmin_r,min_{m_I in Z} sum_I w_I [G_I.r-phi_I,n-2 pi m_I]^2.

The integer branch is chosen deterministically: first by the symmetry-predicted center/seed anchor, thereafter by the periodic image nearest the previously committed center. The wrapped residual and branch fingerprint are color-candidate metadata. Center comparisons use the minimum-image displacement only after orbit-consistent branch alignment.

If W'=WU for a unitary block update, the moment matrices transform exactly as

    Z'_I = U^dagger Z_I U,
    delta Omega_Gamma = -sum_I w_I sum_n (|Z'_I,nn|^2-|Z_I,nn|^2).

For an infinitesimal U(t)=exp(tK), K^dagger=-K,

    dZ_I/dt at t=0 = [Z_I,K],
    dOmega_Gamma/dt at t=0
      = -2 sum_I w_I Re sum_n z_I,n^* [Z_I,K]_nn.

Define the raw anti-Hermitian gradient A by dOmega_Gamma=Re tr(A^dagger K). Because Omega_Gamma has dimensions of length squared, the numerical values 0.5 and 0.2 are not applied directly to A. Inspection of the pinned official Wannier90 v3.1.0 `wannierise.F90` establishes two separate facts: `num_cg_steps=0` sets the Fletcher-Reeves coefficient to zero every iteration, and the trial generator is the current steepest-descent array multiplied by `trial_step/(4*wbtot)`. No gradient-norm normalization enters that trial generator. Define

    W_G = sum_I w_I                 (same stored-shell convention),
    G_c = N_W90(Pi_G(A_c)) = Pi_G(A_c)/(4 W_G),
    G_c^dagger = -G_c,
    K_c^trial = -alpha_l G_c.

W_G is a bounded replicated scalar from the fixed reciprocal shell, not a gradient norm or system-size reduction, and W_G>0 is required. This mapping preserves anti-Hermiticity and finite-symmetry conjugacy and introduces no cross-block norm. The values 0.5 and 0.2 are trial scales, not claims about Wannier90's final accepted rotation: Wannier90 subsequently performs a parabolic line search, whereas this design uses its own bounded monotonic backtracking from K_c^trial. A development/CI fixture compares the synthetic raw-gradient-to-trial-K mapping, including the factor 4, W_G convention, units, and sign, against pinned v3.1.0 at `trial_step=0.5` and 0.2. Trial U, transformed Z_I, and trial delta Omega_Gamma are compared only in a fixture whose moment shell, weights, normalization, and objective are proven identical. The final accepted U is not required to equal Wannier90 because the line searches differ. The reviewed formula, Wannier90 revision, source-line fingerprint, objective fingerprint, normalization fingerprint, length unit, Delta V convention, and shell-weight convention are frozen into test evidence. Production executes the internal formula without linking, launching, or reading Wannier90. Failure of the development/CI gate blocks release of the backend; it is not a runtime distributed-calculation dependency.

The local pair/small-block oracle evaluates these same equations on the union of untruncated active supports and attaches a certified tail interval |delta Omega_Gamma-delta Omega_hat|<=epsilon_tail. A local update is monotone-certified only when delta Omega_hat+epsilon_tail<=-epsilon_accept; an inconclusive interval triggers a complete distributed evaluation or rejection. Its delta interval and directional derivative are compared with a complete distributed evaluation; no surrogate center-distance or nonperiodic r^2 objective may decide acceptance.

All real-space global moments and norms use each physical grid point exactly once. If {C_f} is the disjoint periodic partition of the supercell into fragment cores, production moments are accumulated as

    Z_I,mn = sum_f sum_{r in C_f} Delta V
               w_m(r)^* exp(i G_I.r) w_n(r).

Buffer copies never enter this sum. The existing smooth partition-of-unity weights are used when composing overlapping LCFO fragment contributions and constructing WPW windows; after the unique global field has been stitched onto core ownership, localization does not multiply it by those weights again. Local H/S/gradient actions may read buffer data but assign every output to one core point or canonical coefficient row. The established PoU-weighted nonlocal path is retained at the hybrid handoff and must pass the core-reference equivalence gates below. This prevents both double counting and artificial attenuation at fragment boundaries without introducing a second nonlocal implementation merely for localization.

### Search direction

For every verified finite operation g, let B_g be its action on LCFO coefficients/fields and D_g the desired induced representation on the ordered Wannier center-orbit and irrep columns. They satisfy

    B_g^dagger S B_g = S,
    B_g^dagger H B_g = H,
    B_g C = C D_g,
    D_g^dagger D_g=I,
    B_g B_h = B_gh,
    D_g D_h = D_gh,

after periodic-image shifts and nonsymmorphic Gamma phases have been included in B_g and D_g. If a future spinor/projective representation is enabled, both sides must use the same verified cocycle; it is not part of this spinless first version. The normalized intertwining residual is

    eta_g = ||B_g C-C D_g||_(S,F)/max(1,||C||_(S,F)).

The initial seed must pass max_g eta_g before localization. If a filtered trial block Y has the requested complete representation D_g but is not yet an intertwiner, construct the sole permitted repair by the Reynolds intertwiner

    P_D(Y) = (1/|G|) sum_g B_g Y D_g^dagger.

Then perform S-metric rank revelation only in complete irrep blocks, S-orthonormalize the retained columns, and re-evaluate rank, projector, Ritz, occupation, and max_g eta_g gates. Failure of rank or any gate rejects the complete requested block. Averaging a later generator is not a repair for a nonintertwining seed, and no localization mutation begins before this step passes.

Once C is an intertwiner and W=AC has passed map equivariance, T_g W=W D_g. A real-space update W'=WU preserves the specified centers and irreps exactly only if

    D_g U = U D_g  for all g,

or infinitesimally [D_g,K]=0. For a representative anti-Hermitian trial generator K_ij, the orthogonal projection onto this commutant is

    Pi_G(K_ij) = (1/|G|) sum_g D_g^dagger K_ij D_g,
    K_sym = (Pi_G(K_ij)-Pi_G(K_ij)^dagger)/2.

This average is valid only after the B_g C=C D_g contract and the representation multiplication table have passed. The implementation verifies ||[K_sym,D_g]||, the predicted induced center permutation, irrep action, and exp(K_sym) commutators before field mutation.

Finite symmetry images form connected components induced by shared Wanniers; each bounded component assembles a small anti-Hermitian K_sym^(c) and applies exp(K_sym^(c)) once. Independent disjoint images may retain the 2-by-2 fast path only after commutativity is proven. Degenerate irreps use the same small-block rule. With a trivial instantaneous symmetry group, as is usual for a liquid snapshot, this reduces to ordinary local pair/block Jacobi updates while retaining supercell-periodic support handling and MPI-layout invariance.

The first production version has no primitive-translation FFT update. A block-circulant generator is invalid for the intended interface, surface, and liquid supercells unless an internal translation subgroup is independently verified against the structure, grid, LCFO basis, H, S, occupations, and target projector. Even when such a subgroup exists, exp(K) of a finite-range generator is generally not finite range, so it cannot be treated as a bounded local update. Any future translation-sector accelerator must therefore be a separately designed collective algorithm with its own data redistribution, communication, truncation/correction error, checkpoint, and scaling contract. Colored noncommuting translation layers, polynomial/Krylov/Cayley approximations, truncated-FFT updates, and artificial primitive translations are not permitted in this version.

### Acceptance

A symmetry-connected unitary candidate is provisionally accepted only when a complete evaluation or certified interval for its exact-objective delta Omega_Gamma proves monotonicity and its by-construction intertwining/equivariance, orthogonality, branch, and occupied-projector residuals pass on every affected rank. The truncated shadow is separately checked with moments evaluated from the truncated fields plus certified tails; its result is not inferred from U^dagger Z_I U. A bounded backtracking search may scale the current block generator; failure rejects the candidate without an unconstrained fallback. Provisional color commits remain protected by the sweep checkpoint until post-correction sweep acceptance.

The acceptance scope is the complete verified finite-group orbit, even when its images are physically distant and split across non-neighbor MPI ranks. Spatially disconnected image blocks share one canonical representative generator, step scale, and epoch. A failure at any image rejects and rolls back every image in that orbit.

### Long, history-free convergence schedule

Wannier90 experience for these large Gamma-point localization problems shows that convergence can require several thousand iterations, progressively smaller trial steps, and `num_cg_steps=0`. The production default therefore uses the corresponding history-free steepest/Jacobi direction: every sweep constructs its generator solely from the current symmetry-projected gradient,

    K_s^trial = -alpha_l_eff N_W90(Pi_G(A[W_s])),
    alpha_0 = alpha_search = 0.2,
    alpha_l = max(alpha_min, alpha_finish,1 2^(-(l-1)))  for l >= 1,
    alpha_finish,1 = 0.1,
    l_max = min {l >= 1 : alpha_finish,1 2^(-(l-1)) <= alpha_min},
    n_hold = ceil(f_hold n_acc,max),   f_hold = 0.70,
    l_eff = min(l_base + r_retry, l_max).

Here n_acc is the number of accepted complete sweeps; a rejected attempt and its retry do not increment n_acc, and n_acc,max is the configured maximum. The integer n_hold is computed once at sector entry. Let l_base be the committed/checkpointed scale level and r_retry>=0 be an attempt-controller failure count outside the checkpoint. Level zero is the search scale 0.2; levels one and above are the finishing sequence 0.1, 0.05, 0.025, ... down to alpha_min. There is no previous-gradient vector, conjugate direction, Polak-Ribiere/Fletcher-Reeves coefficient, or localization-CG restart state. All members of one finite symmetry orbit use the same alpha_l_eff. Initially n_acc=0, l_base=0, and r_retry=0. A restored full-sweep monotonicity failure may temporarily select l_eff>l_base for that retry. A REBUILD preserves the current failure count so the restarted attempt uses the same effective scale; it neither advances nor increases the scale. Overrides are recorded and must satisfy

    0 < alpha_min < alpha_finish,1 < alpha_search,   0 < f_hold < 1.

Strict ordering keeps every controller-level advance equal to one real decrease of the trial scale. Duplicate scales are rejected at input rather than represented by fictitious levels.

The search phase deliberately keeps alpha_search=0.2 while n_acc<n_hold: a plateau is logged but cannot commit a smaller level. Convergence gates may nevertheless stop the calculation early. This rule follows the reported convergence trace: fixed 0.2 gave 917.34 Angstrom^2 after 10,000 iterations (best 916.80 Angstrom^2), whereas unconditional staged decay gave 965.49 Angstrom^2; already at iteration 500 the fixed-0.2 result was about 10 Angstrom^2 lower. Thus 0.2 is part of the exploration, not merely a coarse transient step. These values are controller evidence for the default, not a transferable physical-convergence certificate; the user remains responsible for system-specific convergence studies.

For each committed level, retain a bounded rolling queue of exactly the last at most m_stall+1 accepted complete-sweep objectives at that unchanged level. Update this queue after every accepted sweep. At each configured convergence-check sweep, once the queue is full, define Omega_ref as the oldest value, exactly m_stall accepted sweeps before Omega_now, and test a plateau from

    Omega_ref - Omega_now
      <= max(epsilon_stall,abs,
             epsilon_stall,rel |Omega_ref|),

and at least one final retained-plus-omitted gradient gate evaluated at that same convergence check is still above tolerance. If every final gate passes, terminate instead of reducing the scale. At or after n_hold, a plateau at level zero commits level one (0.1); later plateaus halve to 0.05, 0.025, ... . At l_max, record `stalled_at_alpha_min` without inventing another level or clearing the queue. An accepted smaller retry before n_hold is safety-only: reset r_retry, keep l_base=0, return to 0.2 on the next sweep, and reseed the level-zero queue because the accepted objective came from another scale. At or after n_hold, an accepted smaller retry commits l_base=l_eff and seeds that level's queue. This prevents a numerical rejection from permanently destroying early exploration while retaining bounded safe retries. Objective rollback restores the complete queue exactly; REBUILD carries that restored queue into the replacement checkpoint. The stall window, tolerances, scales, f_hold/n_hold, convergence-check interval, counters, and maximum accepted-sweep count are explicit inputs or log fields. Backtracking within a candidate remains bounded and uses only the current K_s^trial.

Convergence is decided only by the spread-change, generator, and complete retained-plus-omitted gradient gates. Reaching the configured maximum sweeps reports nonconvergence with the final residuals; it is not treated as evidence that the method is unsuitable merely because a short iteration budget was exhausted. The optional future use of localization CG requires a separate design and comparison gate and is outside the first production route. This restriction does not set the iteration count of the later sparse metric, block-CG, or hybrid-SCF linear solvers, which solve different problems and retain their existing bounded controls.

### Wannier/WPW selection

Localization and tail gates reduce over a complete symmetry block. A block is retained wholly or rejected wholly. Rejected blocks produce symmetry-complete WPWs, including complete reciprocal stars. Acceptance also requires the hybrid span residual ||(1-P_hybrid)P_target||, or equivalently the worst principal angle, to pass. Occupied-block replacement must reproduce the occupied projector explicitly; rank equality alone is insufficient.

## Distributed representation

Each Gamma-supercell Wannier has an immutable canonical anchor fragment chosen from the symmetry-adapted seed, and its control owner is the rank that owns that anchor. Independently, every physical grid point belongs to exactly one disjoint fragment core. The core owner holds the authoritative values of every active orbital at that point; a fragment buffer contains only cached copies imported from the corresponding core owners. The control owner stores canonical IDs, occupation/representation block IDs, active support metadata, adjacent distributed tails, periodic moments, current wrapped center, overlap-graph neighbors, core generation, and buffer generation. A Wannier center may cross a fragment boundary without changing its canonical ID, control owner, orbit representative, or coordinator. Its field remains distributed over the authoritative cores intersecting the active support; anchor ownership does not mean gathering an orbital onto its owner. If load imbalance from immutable anchors violates a production bound, the route fails for a separately designed collective repartition rather than migrating owners inside localization.

The fixed fragment buffer and the movable Wannier support envelope are different objects. The support envelope may cross any number of fragments and determines which core owners participate; it is not required to fit inside one fragment buffer and never causes that buffer to become the authoritative storage for a whole orbital. Operator application has exactly two permitted data paths.

For a local stencil operator O evaluated on core C_f, the fragment buffer supplies the exact input neighborhood N_O(C_f), and setup proves

    N_O(C_f) subseteq B_f

for every locally evaluated finite-difference Laplacian/gradient, local retained H/S stencil, and partition-of-unity transition derivative used by seed or WPW construction. The bound includes the largest finite-difference half-width and local transition/stencil reach.

A separable nonlocal pseudopotential never requires its full projector support to fit in one fragment buffer. The production implementation keeps the existing PoU-weighted projector assembly and its established ownership/communication path. It must nevertheless be equivalent, on bounded verification blocks, to the disjoint-core reference

    q_a,mu(psi) = sum_f sum_{r in C_f} Delta V
                    p_a,mu(r)^* psi(r),
    (V_NL psi)(r in C_f) = sum_a,mu p_a,mu(r) q_a,mu(psi).

The verification proves PoU completeness, projector single counting, Hermiticity, and equality of small-block nonlocal matrix elements and actions against this reference. A focused failing test may authorize a minimal correction to the existing path; the first implementation does not create a new canonical-projector owner or reduction schedule. Failure of local-stencil closure, PoU completeness, nonlocal single counting, or the bounded reference comparison is FATAL. Values at the fragment-box edge are never treated as physical zeros or as a periodic boundary of the fragment.

Only authoritative core field values, not buffer copies, are checkpointed, shadowed, and committed. A sweep checkpoint stores those fields and only the minimal replay state: previous wrapped centers/branches, active and truncation masks, n_acc, l_base, the bounded same-level rolling queue of at most m_stall+1 accepted objectives, support/graph/orbit epochs, and core generations. The attempt-controller r_retry is deliberately outside the restored payload: it increments after an objective restore, survives a REBUILD of that same attempt, and resets only after acceptance or terminal failure. The verified N_W90/objective/unit/quadrature fingerprints and input scales are immutable setup provenance rather than checkpoint payload. Moments, spreads, defects, correction-cluster membership, participant-local cache state, and buffer generations are derived data; after restore they are discarded and recomputed deterministically from the restored cores and versioned graph. Buffers carry (source_core_generation,graph_epoch) and are valid for candidate evaluation only when their generation equals the committed generation of every dirty source. A stale buffer before refresh is normal and is refilled; a mismatch that remains after the scheduled refresh is FATAL. The color order is fixed: refresh required buffers from committed cores; compute core shadows; reduce OK/REJECT/REBUILD/FATAL; swap or discard core shadows and authoritative metadata; synchronize; refresh dirty buffers; then enter the next color. REJECT never publishes its buffers, and checkpoint restore invalidates all newer derived data and refills buffers from restored cores. Remote T_g audits compare authoritative source and target core values; buffer equality is only a cache-consistency diagnostic.

Setup constructs a distributed symmetry-orbit directory using the existing affine fragment/point maps and ownership exchange machinery. For each canonical orbit it records the representative ID, canonical anchor, symmetry operation and periodic image, image fragment/Wannier ID, control owner rank, local representation action, and a bounded owner list. It also stores one canonical transporter from the representative to every image: exact global-grid permutation, source/target field-owner maps, periodic integer shift, nonsymmorphic Gamma phase, and irrep matrix. The directory is sharded by canonical orbit ID; no rank stores all system orbits.

Remote transporters must be path independent. For every g,h and every canonical row/column ID, direct transport by gh and composed transport by h then g must produce the same target ID, grid permutation, periodic shift modulo the supercell, phase, field-owner endpoint, and representation action. With T_g T_h=T_gh, the affine data obey R_gh=R_g R_h and t_gh=R_g t_h+t_g, while the reciprocal-shell permutation is contravariant, pi_gh=pi_h o pi_g. The directory verifies the global representation laws

    T_g T_h = T_gh,
    B_g B_h = B_gh,
    D_g D_h = D_gh,

For a cluster-dependent column transporter defined below, the corresponding basepoint-dependent cocycle is instead

    D_gh,c = D_g,hc D_h,c,

because T_h maps c to hc before T_g acts. The equality includes column order, periodic-image phase, and irrep action. Both global laws and this cluster cocycle use fingerprints that include source/target canonical IDs, integer shifts, phases, graph/orbit version, and participant ownership. A mismatch is a construction failure; no arbitrary shortest group word is accepted as authoritative.

T_g must map the production grid by an exact integer permutation. Symmetry enforcement never interpolates between grid points. If a required operation cannot be represented exactly on the chosen cell/grid, the complete localization/hybrid route stops; WPW replacement cannot repair a grid that lacks the required action. A reduced verified group is permitted only by an explicit input policy and requires discarding the old seed/checkpoint and rebuilding the seed, orbit directory, induced representations, graph, and all gates from the beginning. No operation is silently dropped during a run.

Each finite point-block update of W is unitary in the core physical inner product before sparsification. Update colors create symmetry-related sparse W candidates in shadow buffers but do not perform polar correction. At the color/sweep boundary, sparse physical-overlap edges define bounded correction clusters; every cluster is closed under all finite symmetry orbits it touches and includes every coupled core/buffer rank. A separately colored correction phase also uses shadow buffers. No post-stitching truncation is projected back to LCFO coefficient space.

For each cluster orbit, one canonical representative fixes the truncation-mask orbit and cluster membership. Let H_c={g:T_g c=c} be the cluster stabilizer. The code always assembles the bounded complete core overlap M_c=W_c^dagger W_c, with every core point counted once, verifies its positive spectrum and H_c covariance, and applies one simultaneous block-polar correction

    Q_c = M_c^-1/2,
    W_c' = W_c Q_c,
    D_h,c^dagger M_c D_h,c=M_c,
    [Q_c,D_h,c]=0  for all h in H_c.

The inverse square root is evaluated as one bounded small-cluster matrix function with verified residuals; it is not a global dense orthogonalization. For the fixed ordered columns of c and its image gc, define the verified column-space unitary D_g,c, including column permutation, phase, and irrep action, by

    T_g W_c = W_gc D_g,c.

The remote correction is therefore

    Q_gc = D_g,c Q_c D_g,c^dagger,
    T_g(W_c Q_c) = (T_g W_c)Q_c
                  = W_gc Q_gc D_g,c.

The code compares this Q_gc with the inverse square root independently formed from M_gc and verifies D_gh,c=D_g,hc D_h,c along direct and multiword remote paths. For a stabilizer h in H_c, every covariance/commutator uses the basepointed D_h,c. T_g is never applied to or composed with the column-space matrix Q_c. There is no ordered multiplicative-Schwarz alternative.

Before correction, the cluster must also certify its exterior physical-overlap coupling. With W_bar_c denoting all columns outside the cluster, retained core boundary blocks and omitted-tail bounds provide

    eta_c^ext >= ||W_bar_c^dagger W_c||_2.

Correction uses one read-only preparation pass before any correction color commits. For every candidate cluster c, form M_c and Q_c from the same post-update committed W and obtain a certified q_c>=||Q_c||_2 from the verified eigenvalue interval of M_c. For every exterior coupling between prepared clusters c and d, require the two-sided bound

    ||Q_d^dagger W_d^dagger W_c Q_c||_2
      <= q_d eta_dc^ext q_c <= epsilon_ext.

If it fails, merge/add the offending physical-overlap neighbor and its complete finite-symmetry orbit, discard the entire prepared correction plan, and restart the read-only preparation pass. No correction shadow has yet been committed. Continue until every two-sided bound passes or the configured cluster/envelope bound is reached. Then freeze that prepared cluster/Q/q plan for the correction transaction and apply its colors. A color rejection discards its shadow without changing the prepared plan. After all correction colors, evaluate one complete retained-plus-omitted global post-correction exterior bound before sweep acceptance. Failure restores the whole sweep checkpoint; it never expands or compounds a previously corrected shadow. Envelope escape returns REBUILD after restore. Size, conditioning, fill, iteration, or communication exhaustion also restores the whole sweep before the complete affected-block fail-closed/WPW decision. The global sparse orthogonality defect, truncation, correction, spread, projector, and correction-intertwining defects gate the sweep checkpoint. Cluster radius, two-sided prepared and final global exterior bounds, q_c, size, fill, inverse-square-root iterations, overlap, stabilizer, and conjugacy are feasibility gates.

A final support compression uses the same error contract and additionally checks energy and local H/Z defects. Failure enlarges the local correction/support radius or sends the complete block to WPW; it never invokes dense orthogonalization.

The sparse pair graph has a vertex per Wannier center in the Gamma-point supercell and an edge for a periodic-image overlap that can occur anywhere inside a configured maximum active/correction envelope around the immutable seed anchor a_n. Define its conservative geometric radius by the minimum-image Minkowski sum

    R_env,n = R_center,max,n + R_active,max,n
              + R_correction,max + R_overlap,

where R_center,max,n bounds ||r_n-a_n||, R_active,max,n bounds field support about r_n, R_correction,max covers cluster closure, and R_overlap covers only retained physical-overlap reach plus its certified tail margin. Finite-difference, H, coefficient-space S, and nonlocal-projector reach belong to independent fragment-operator schedules and do not enlarge the Wannier conflict graph. Setup constructs the core points, source buffers needed by localization, rank sets, halo endpoints, periodic images, update/correction participants, both conflict graphs, and deterministic colors intersecting this whole region and all of its finite-symmetry images. The correction conflict graph treats two canonical clusters as conflicting whenever their maximum possible symmetry-closed correction envelopes overlap; therefore runtime expansion inside those envelopes cannot introduce a same-color conflict or require recoloring. The current center and nonzero support may move inside these bounds without discovery, rehashing, rebuilding, or recoloring; communication is still restricted to dirty active data. A center, field support, physical-overlap neighbor, or correction closure attempting to exceed any component bound returns REBUILD. Then no later mutation in that sweep is allowed: all ranks restore the sweep checkpoint, enlarge and rebuild the envelope and schedules collectively, advance the graph epoch, replace the restored old-epoch checkpoint with a new checkpoint of the same W/counters under the new epoch, and only then restart. Repeated envelope growth beyond configured radius, degree, memory, or restart bounds fails closed.

Graph sparsity is not itself a stationarity proof. Let A be the full anti-Hermitian gradient defined by dOmega_Gamma=Re tr(A^dagger K). Retained edges provide exact blocks A_E. For every omitted pair (i,j), support/tail Cauchy bounds provide beta_I,ij >= |Z_I,ij| and hence a certified gradient bound gamma_ij >= |A_ij|; a valid conservative choice is gamma_ij=4 sum_I w_I beta_I,ij. Accumulate

    r_i^omit = sum_{j:(i,j) notin E} gamma_ij,
    ||A_omit||_2 <= sqrt(max_i r_i^omit max_j c_j^omit),
    ||A||_2 <= ||A_E||_2,bound + ||A_omit||_2,bound.

The retained contribution is also evaluated by block row/column sums, not only the largest pair. Convergence requires the last bound, the maximum symmetry-projected block directional derivative, and delta Omega_Gamma to pass. If omitted row sums grow with system size, support must expand or the complete block fails to WPW; an entrywise threshold alone is invalid.

Steady-state communication is limited to dirty neighbor-tail exchanges, small connected-block coordination messages, two scalar synchronizations per executed color, and infrequent convergence reductions of spread and maximum residual. There is no field all-gather or N_W by N_W all-reduce. The color synchronizations are the deliberate cost of the simpler atomicity rule and must be measured rather than hidden.

Communication has two distinct planes:

- the field plane exchanges values, tails, and sparse-correction halos only between spatial support neighbors;
- the symmetry control plane exchanges representative generator parameters, step scale, truncation-mask descriptors, correction status, and color validity among the prebuilt bounded participant set of a verified finite-group orbit, regardless of MPI-topology distance.

A third, low-frequency symmetry audit plane streams bounded chunks of T_g-mapped field/support data directly between symmetry-related field ranks. Before a sweep begins, committed fields and metadata receive a recoverable sweep checkpoint. The audit runs after update and correction phases but before the checkpoint is released. Audit failure restores the complete sweep checkpoint; it is not attributed to any one previously committed update color.

The orbit coordinator is selected deterministically from the canonical orbit ID. Setup derives each update participant set from the maximum envelope and each correction participant set from its finite-symmetry closure. These immutable per-graph-epoch sets are checked once by a two-sided exchange; they are not rediscovered in the sweep.

Atomicity is enforced at color granularity with one synchronous procedure for both update and correction colors. Each rank returns one ordered status

    OK < REJECT < REBUILD < FATAL,

and a blocking MPI_MAX reduction selects the common action:

1. every participating rank first verifies/refills generation-matched read-only buffers, then builds core-only local candidates in shadow storage and validates the common graph epoch, transporter, objective/correction residuals, and bounds;
2. the reduction yields OK when all candidates pass, REJECT for an ordinary objective/bounded-backtracking rejection, REBUILD for an envelope escape, or FATAL for a non-finite value, representation/version inconsistency, or other construction failure;
3. on OK all ranks swap the color's shadows; on REJECT all discard and continue; on REBUILD all discard, restore the sweep checkpoint, and rebuild before any further mutation; on FATAL all discard and terminate the route; and
4. after OK or REJECT, a blocking post-action synchronization completes, dirty committed cores refill their dependent buffers, and only then may any rank start the next color.

Thus no rank can enter a later color while another is deciding or applying the current one, and no per-pair multistage state machine is required. MPI process-failure recovery is out of scope. Per-update acceptance checks parameter/generator/mask equivariance by construction. The post-correction audit evaluates ||W_g-T_g W_rep|| and conjugacy residuals of metric, moments, truncation, and correction. Audit or envelope failure is handled only at a synchronized color or sweep boundary and collectively restores the sweep checkpoint; the audit detects faults but is not the mechanism that constructs remote covariance.

When the verified group contains only the identity, the implementation bypasses the orbit directory, representation tables, remote symmetry control, and remote audit. It uses the same conservative spatial graph, color-level shadow commit, block-polar correction, checkpoint, and numerical gates. This is the normal liquid-snapshot fast path, not a relaxation of supercell periodicity or MPI-invariance requirements.

## Localization flow

1. Define occupied and finite empty windows, closing their boundaries over degeneracies and irreps.
2. Construct symmetry-complete local seeds from distributed sparse spectral projectors and pass seed/projector feasibility gates.
3. Assign immutable canonical anchor owners while leaving field values spatially distributed. For a nontrivial verified group, build and composition-check the sharded canonical-transporter directory and remote audit schedule; for the identity group, bypass them.
4. Build the maximum admitted spatial envelope, its fixed participant sets, periodic-image overlap graph, finite-symmetry image blocks when present, separate update/correction conflict graphs, and deterministic colors.
5. Localize one occupation sector at a time: the complete occupied sector, then the requested finite empty sector. Compute n_hold=ceil(f_hold n_acc,max) once, keep committed level zero through the search phase n_acc<n_hold, and select the deterministic history-free retry level l_eff=min(l_base+r_retry,l_max). Save a sweep checkpoint. For each update color, construct K^trial from the current symmetry-projected gradient and the verified 1/(4 W_G) factor, apply bounded monotonic backtracking, build shadows, reduce the four-state action, collectively swap or discard, and synchronize. Prepare the complete symmetry-closed correction plan read-only, merging clusters and restarting preparation until every two-sided bound passes; freeze it before applying any correction color. Apply the prepared correction colors, then evaluate the global post-correction exterior bound and complete Omega_Gamma and classify one collective sweep outcome. ACCEPT requires monotonicity and every finite numerical validity gate. SCALE_REJECT is restricted to a finite, scale-responsive objective or truncation/correction residual failure for which the configured envelope and all construction invariants still pass; it restores the checkpoint, leaves n_acc/l_base/queue unchanged, increments r_retry once, and retries up to the fixed limit. REBUILD is restricted to envelope escape; it restores while preserving r_retry and rebuilds as specified below. FATAL covers non-finite data, post-refresh generation mismatch, representation/transporter/version inconsistency, PoU/core reference failure, or any other construction invariant failure; it restores if mutation occurred and terminates without changing the scale controller. Any bounded block-local feasibility failure first restores the complete sweep checkpoint and invalidates shadows, buffers, and derived state. Only afterward may the collective route rebuild the W/P partition, graph, orbit directory, participant sets, and epochs for the explicitly identified complete block or fail closed; it is never relabeled SCALE_REJECT. The remote audit runs only for nontrivial symmetry; its failure restores and terminates because repeating identical data cannot repair covariance.

On ACCEPT, release the checkpoint and increment n_acc once. If a smaller retry succeeds before n_hold, reset r_retry, retain l_base=0, reseed the level-zero queue, and return to alpha_search on the next sweep. At or after n_hold, commit a successful smaller retry into l_base=l_eff and seed its queue. Otherwise append the accepted Omega, retaining at most m_stall+1 values at the unchanged level. On a configured convergence-check sweep, stop if all fresh gates pass. If a gate remains open and the full queue passes the plateau inequality, log but ignore the plateau before n_hold; at or after n_hold advance exactly one level and reseed the queue, unless already at l_max, where the controller retains both and records `stalled_at_alpha_min`. SCALE_REJECT retry exhaustion reports the entire currently localized occupation sector as nonconverged with residuals; it never infers offending blocks from the aggregate sweep result. A block-local handoff is permitted only after the full restore specified above. On REBUILD, restore without changing n_acc, l_base, r_retry, or the rolling queue; enlarge and rebuild the envelope/participants/graphs/colors; advance the graph epoch; discard the old checkpoint; save a new checkpoint of the same restored W/counters/queue under the new epoch; and only then restart the attempt at the same effective scale. Fail closed if configured rebuild bounds are exceeded.
6. Reduce scalar convergence metrics only at the configured sweep interval.
7. Stop only when spread change, maximum generator, symmetry-projected block derivative, retained-gradient row/column bound, and omitted-gradient accumulation all pass. This establishes a fixed-accuracy symmetry-constrained stationary local minimum, not a proof of the global minimum. Several thousand sweeps are permitted. Maximum-sweep exhaustion reports nonconvergence and residuals rather than silently accepting or switching to localization CG.
8. Emit the fingerprints and invariant diagnostics needed for user-run fragment-buffer, active-support, and compression convergence studies; do not launch comparison runs inside the localization route.

The inner loop performs no allocation, dense group average, avoidable division, or complete spread reevaluation. Phase factors, reciprocal lengths, symmetry maps, schedules, and connected-block metadata are cached before sweeps.

## Fragment-buffer validity and user convergence

The code proves discrete buffer validity for the selected run: unique core ownership and coverage, partition-of-unity sum/gradient identities where used, operator-neighborhood closure, nonlocal-projector single counting, generation consistency, and equality of buffer-assisted versus direct core-reference H/S/gradient/Z actions on bounded test cases. These are correctness checks, not a claim that the user-selected physical buffer or localization radius is converged.

Production convergence across spatial sizes is a user-run study. The program reports invariant diagnostics needed to compare separate runs but does not launch a second calculation or automatically declare physical convergence. The user varies one quantity at a time while keeping structure, grid, seed fingerprint, tolerances, MPI-invariant IDs, and all other radii fixed:

1. fragment-buffer study: increase the fragment buffer while keeping core partition and active support fixed; compare core H/S/gradient actions, Omega_Gamma, center/projector orbits, energy, and transitions;
2. active-support study: rerun localization from the same canonical seed at larger R_active,max and compare aligned center/projector orbits, spread, tail, covariance, local H/Z, energy, and transitions; and
3. compression study: apply two final compression radii to the same converged orbitals and compare the same invariant quantities.

The user records the chosen widths and differences in the production-results document. The code only rejects an internally invalid buffer, failed operator closure, or a result that violates the configured single-run numerical gates.

## Hybrid handoff and initial state

Accepted real-space blocks first form a fixed core-field W and its physical projector

    P_W^core = W (W^dagger W)^-1 W^dagger,

where every grid point is counted once. For each rejected complete target block, generate partition-of-unity WPW candidates whose windows, centers, reciprocal stars, periodic phases, and irrep action form complete induced space-group orbits. Project each candidate out of W with P_W^core, then perform rank revelation and orthonormalization only by complete symmetry blocks. Bounded window expansion is the only retry; exhaustion of its support/rank/conditioning bound fails closed. The resulting P columns are frozen with W to form B=[W,P].

This is the explicit bridge between core-only localization and the existing PoU hybrid assembly. On bounded blocks, assemble S_B^PoU and H_B^PoU through the production path and independently assemble the disjoint-core references S_B^core=B^dagger B and H_B^core. Require bounded ||S_B^PoU-S_B^core|| and ||H_B^PoU-H_B^core|| together with Hermiticity, symmetry, nonlocal single-counting, and energy/action residuals. The production route keeps the PoU matrices after these gates pass; a mismatch is a representation-handoff failure, not a reason to reinterpret W in the LCFO coefficient metric.

Hybrid-span validation uses the production PoU inner product as its authoritative metric; coefficient-space LCFO S is never applied directly to the post-stitching fields B, W, P, or a physical target field. Map each retained coefficient-space target block once through the already verified stitching/target-field construction. The resulting symmetry-complete bounded-support physical column blocks X span range(P_target); the projector operator itself is never treated as a column block. For each block assemble through the production path

    S_B^PoU = <B|B>_PoU,
    C_BX^PoU = <B|X>_PoU,
    S_X^PoU = <X|X>_PoU,
    S_B^PoU Y = C_BX^PoU,
    P_B X = B Y.

The equivalent core-reference Grams S_B^core, C_BX^core, and S_X^core are assembled independently and must agree on bounded verification blocks before the PoU result is used. The residual Gram G=<R|R>_PoU is formed entirely in the PoU physical-field metric, either explicitly for R=X-BY or by the algebraically equivalent Hermitian Schur expression with a verified solve residual. Because X need not be PoU-orthonormal, the invariant squared span error is the largest generalized eigenvalue

    epsilon_span^2 = lambda_max(G,S_X^PoU)
                   = ||(S_X^PoU)^(-1/2) G (S_X^PoU)^(-1/2)||_2.

Production uses one conservative bound and no generalized eigensolver. Partition the canonical target columns into the fixed bounded blocks used by the sparse Gram representation. With S_X,ij^PoU denoting those blocks and epsilon_S,omit a certified upper bound on the omitted block-row operator contribution, define the block-Gershgorin lower bound

    s_X,min = min_i [lambda_min(S_X,ii^PoU)
                     - sum_(j != i, retained) ||S_X,ij^PoU||_2]
              - epsilon_S,omit.

Require s_X,min>0; the omitted budget is always subtracted from the retained lower bound. Separately accumulate retained and omitted residual-Gram blocks into verified one- and infinity-norm bounds G_1_hat and G_inf_hat. Then

    G_2_hat = sqrt(G_1_hat G_inf_hat) >= ||G||_2,
    epsilon_span^2 = lambda_max(G,S_X^PoU)
                   <= G_2_hat/s_X,min,
    epsilon_span <= sqrt(G_2_hat/s_X,min).

This sufficient bound is evaluated for the single deterministic canonical target-column normalization/order frozen by the seed fingerprint. The exact generalized eigenvalue is basis invariant, but the conservativeness of this bound need not be; arbitrary runtime target-basis changes are forbidden. Dense CI fixtures compare the bound with the exact principal angle and exercise rescaled/ill-conditioned targets only to prove the bound remains safe or fails closed, not numerically invariant. Failure of Hermiticity, positive s_X,min, conditioning, fill, omitted budgets, or the final span tolerance fails closed. No inverse square root, target-metric solve, Lanczos iteration, dense complement projector, or dense global target Gram is formed. O(N) is claimed only where target support/degree, G/S_X fill, row/column accumulation, omitted budgets, workspace, and communication per rank remain bounded in the weak-scaling gate.

Initial hybrid coefficients use the established bounded block-CG and safeguarded hybrid mixing. Inner CG counts remain capped because excessive inner iteration can oscillate. This does not add Pulay mixing to the DC solve.

The 2026-08-24 Si64 reference establishes the following minimum regression for this handoff. Starting from the converged DC+LCFO density, assemble kinetic, local, and nonlocal actions in the WF+WPW basis, apply the anti-representation-consistent full-group pencil average, solve the distributed generalized eigensystem, rebuild density, and repeat the DG Hamiltonian SCF. On the accepted run the first post-average residuals were approximately H `2.16e-14`, S `2.19e-14`, and density `1.03e-15`; these stayed at roundoff through convergence. The raw local component was not separately symmetric before averaging and decreased with the SCF density residual, whereas kinetic and nonlocal components were already covariant to roundoff. Production tests must preserve this distinction rather than require every density-dependent raw local component to commute before group projection.

## Scaling model

With N_W=O(N), fixed requested accuracy, bounded core and fragment-buffer points per rank, active support n_s, graph degree z, point-block size b, correction-cluster radius/overlap/iterations, WPW/target-block counts and supports per atom, metric and hybrid-projector condition numbers, and solver/sweep iteration counts all bounded:

- authoritative field storage is O(N_W n_s)=O(N), while read-only fragment-buffer duplication is bounded per rank and is not multiplied by the number of remote symmetry images;
- local point/correction work per sweep is O(N);
- field-halo traffic is bounded per rank in weak scaling;
- symmetry-control traffic is proportional to bounded color-participant fanout and absent on the identity-group fast path;
- low-frequency audit traffic is O(N) globally and bounded per rank at fixed local support;
- color validity/post-swap synchronizations add O(n_color log P) latency per sweep, with n_color required to remain bounded at fixed accuracy;
- convergence reductions add O(log P) latency at reduced cadence; and
- bounded-window WPW sparse application remains O(N).

These are conditional fixed-accuracy bounds. The localization sweep count may be a large size-independent constant, including several thousand, and is measured rather than assumed small. Weak-scaling gates measure authoritative core bytes, read-only buffer duplication/refresh bytes and generation stalls, n_s, z, b, update/correction color counts and collective latency, correction-plan restart count, cluster reach/overlap/merges, q_c, two-sided/final exterior bounds and inverse-square-root iterations, stabilizer-block size/conditioning/frequency, color participant fanout and control bytes, immutable-anchor load imbalance, envelope size and rebuild frequency/cost, checkpoint/audit bytes, n_pw, target-block support/RHS counts, G/S_X fill, s_X,min, row/column work and omitted-operator budgets, metric/filter/projector iterations and conditioning, localization sweeps, accepted sweeps per scale, plateau detections, block-CG iterations, and SCF iterations.

## Fail-closed conditions

- incomplete group orbit or inconsistent multiplication table;
- a geometric-product/pullback convention mismatch, a path-dependent generator word, or a generator relation defect under the anti-representation convention;
- failure of the requested block's representation/localization feasibility gates, including possible topological obstruction;
- sparse seed projector or metric purification exceeding fill, residual, or memory bounds;
- inconsistent internal-unit conversion, cell/grid volume, Delta V multiplicity, PoU/quadrature separation, projector convention, or reciprocal-shell weight multiplicity;
- N_W90 step mapping mismatch, or an attempted one-step U/Z/Omega comparison between unequal objective fingerprints;
- LCFO-to-core stitching map failing metric pullback, linearity, or symmetry-equivariance gates;
- missing/invalid generalized H/S spectral bounds, gap, or metric conditioning;
- empty-window boundary that cannot be closed over a bounded symmetry-complete shell;
- symmetry mixing unequal occupations;
- missing bounded core-owner exchange/reduction needed to supply an operator neighborhood or connected block;
- incomplete or inconsistent remote symmetry-orbit owner directory;
- reciprocal moment shell/weights not closed under the verified finite group;
- noninteger/nonunimodular reciprocal lattice action or inconsistent field-action phase convention;
- path-dependent remote transporter, inconsistent periodic shift/phase, or mismatched source/target owner endpoint;
- non-permuting production-grid action or a symmetry operation requiring interpolation;
- required symmetry that is not representable on the grid (route failure, not WPW fallback), or an attempted mid-run reduction of the verified group;
- an internal primitive translation or block-circulant structure inferred without verification against the structure, grid, LCFO basis, H, S, occupations, and target projector;
- remote field/metric/moment/truncation/correction covariance residual beyond tolerance;
- remote orbit-owner fanout exceeding the finite-group bound;
- stale symmetry epoch or inconsistent representative generator/step scale;
- graph-epoch or fixed participant-set disagreement at a color boundary;
- non-finite moment, angle, spread, or residual;
- owner disagreement during connected-block commit;
- monotonicity failure after bounded backtracking;
- orthogonality/projector/covariance violation;
- truncation/correction failure or support growth beyond the fixed-accuracy bound;
- nonconjugate remote block-polar correction or truncation mask;
- a stabilizer-closed block-polar cluster exceeding size, conditioning, fill, or communication bounds;
- attempted support/correction outside the admitted envelope;
- correction-cluster or envelope-rebuild growth beyond the fixed-accuracy bound;
- hybrid matrix-free projector solve failing residual/condition/iteration bounds;
- invalid core ownership, buffer generation, operator closure, PoU completeness, nonlocal-projector single counting, or PoU/core handoff contract; or
- hybrid rank loss or unacceptable sparse-metric conditioning.

No failure silently invokes global dense localization or accepts a partial orbit.

## Validation

Small exact tests compare generalized sparse-filter seeds, Ritz/generalized residuals, and S-projectors with dense generalized eigensystem references; reciprocal-shell metric completeness, wrapped centers, analytic directional derivatives, and local delta Omega_Gamma with complete-field finite differences; and finite point-block updates against explicit dense K_sym exponentiation. They include a trivial-symmetry liquid-like fixture, a noncommuting shared-Wannier orbit, periodic-boundary support pairs, branch crossings, an accumulated omitted-gradient counterexample, rejection of unverified internal translations, a nonintertwining seed, nonsymmorphic operations, multidimensional irreps, controlled truncation/block-polar correction, projector preservation, monotonicity, and color shadow-buffer rollback.

Remote-symmetry tests deliberately place rotated/nonsymmorphic images of one finite-group orbit on physically distant non-neighbor MPI ranks. They verify exact integer/unimodular reciprocal actions, group-closed shells, field-action phase sign, and moment covariance; exact global-grid permutations T_g; equality of direct and multiword transporter paths including target IDs, periodic shifts, phases, core-owner endpoints, and D_g actions; setup-time equality of conservative participant sets; conjugated generators and block-polar corrections; color-wide discard after one remote rejection; field/metric/moment/truncation/correction covariance within reduction tolerance; center and active-support motion inside the envelope without owner migration or graph rebuild; collective checkpoint restore and rebuild on envelope escape; core-only remote audits and buffer-generation diagnostics; rank reordering invariance; rejection of interpolation-only maps; and that remote symmetry enforcement does not enlarge the fragment buffer. A separate identity-group test proves that liquid-like inputs allocate no orbit directory and send no remote symmetry-control or audit traffic.

Small periodic insulating fixtures compare spread, center orbits, invariant projectors, and local H/Z matrices with the full-cell/Wannier90 route. Orbital sign, phase, and rotations within degenerate subspaces are not comparison targets. The symmetry-pencils fixture includes noncommuting D3 with a redundant order-sensitive generator, a direct full-group average, and a deliberately unitary but relation-inconsistent generator set.

MPI tests run with 1, 2, 4, and 8 ranks, checking identical block selection/coloring, invariant observables, bounded seed/localization workspace, and absence of N_W-squared traffic. Material progression is: dense-versus-sparse generalized seed feasibility, exact synthetic localization fixtures, a physically meaningful buffered cell, Si64, a surface or interface, an ionic oxide, liquid water, then replicated disordered/surface weak scaling toward 1,000, 10,000, and 100,000 atoms.

The completed 2026-08-24 Si64 Wannier90 plus Hybrid-DG run is the reference baseline summarized above. It validates corrected units, the pullback convention fix, the WF+WPW Hamiltonian handoff, full-group pencil averaging, generalized diagonalization, density reconstruction, and DG-SCF convergence. It is not an implementation prerequisite for each spatial-localization test and does not satisfy the spatial backend's production or weak-scaling acceptance by itself.

## Primary-method references

- N. Marzari and D. Vanderbilt, [Phys. Rev. B 56, 12847 (1997)](https://doi.org/10.1103/PhysRevB.56.12847), defines MLWFs by unitary minimization of total spread on a fixed band subspace.
- P. L. Silvestrelli, [Phys. Rev. B 59, 9703 (1999)](https://doi.org/10.1103/PhysRevB.59.9703), and G. Berghold et al., [Phys. Rev. B 61, 10040 (2000)](https://doi.org/10.1103/PhysRevB.61.10040), give Gamma-point periodic-supercell localization functionals and arbitrary-cell metric weighting.
- R. Sakuma, [Phys. Rev. B 87, 235109 (2013)](https://doi.org/10.1103/PhysRevB.87.235109), gives the site-symmetry/induced-representation constraint on the Wannier transformation and notes that a symmetry-constrained solution need not be the unconstrained global minimum.
- H. J. Xiang et al., [J. Chem. Phys. 124, 234108 (2006)](https://doi.org/10.1063/1.2207622), demonstrates local-orbital initialization and Jacobi-sweep linear-scaling MLWF localization under locality assumptions.

## Non-goals

- metals or entangled metallic Fermi surfaces;
- many independent Wannier90 processes;
- full-system field/gauge gathering;
- dynamic basis adaptation during RT;
- post-hoc unconstrained symmetry repair;
- primitive-translation FFT or truncated matrix-exponential acceleration;
- automatic irregular repartitioning; and
- orbital-by-orbital equality with Wannier90 in degenerate subspaces.

## Production acceptance

The route is ready only after exact and MPI invariant tests pass, including group-closed reciprocal moments, multiword remote-transporter equality, synchronous color commit/discard, conjugated block-polar correction, immutable-anchor center crossing, non-neighbor remote orbits, the identity-group bypass, periodic boundaries, and bounded envelope rebuild; Si64 plus at least one surface/interface and one liquid fixture reach the hybrid initial state without global replicated storage; the W/PW basis passes S-projector/symmetry/rank/metric/energy/transition gates; maximum remote orbit-owner and color-participant fanout is bounded by the verified finite action and admitted envelope; memory per rank remains bounded under weak scaling; and measured sweep work plus field/control traffic are linear in local problem size over the tested range.

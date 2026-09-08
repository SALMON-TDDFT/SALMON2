# Spatially Parallel Symmetry-Constrained MLWF Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a spatially owned, symmetry-constrained local Jacobi MLWF backend for large Gamma-point interface, surface, and liquid supercells, and hand complete failed localization blocks to the existing WPW complement.

**Architecture:** Construct S-orthogonal local seeds from an explicit sparse generalized (H,S) filter. Apply only finite symmetry operations verified against the actual Gamma-point structure, grid, LCFO operators, and target projector; otherwise use local pair/block Jacobi updates with periodic-image support. Prebuild one conservative maximum-support envelope, commit shadow updates synchronously by color, and always use bounded symmetry-closed block-polar correction. Bypass remote-symmetry machinery for the identity group. Keep Wannier90/full-cell as a small reference oracle and reuse existing hybrid sparse solvers and RT handoff.

**Tech Stack:** Fortran 2008, MPI, SALMON DC+LCFO/DG modules, Python regression drivers, CMake, and Wannier90 reference fixtures.

---

## Execution rules

- Work test-first: create one failing contract, observe the expected failure, add minimal implementation, rerun focused tests.
- Do not alter corrected SALMON/Wannier90 units or existing symmetry-file semantics.
- Never add a global field all-gather, dense global gauge/DMN, per-fragment Wannier90 call, primitive-translation FFT/truncated-exponential path, artificial internal translation orbit, or partial symmetry fallback.
- Keep occupied and empty sectors separate.
- After every route change run `python3 tests/dg/check_dg_overlapping_wannier_route.py`.

## Accepted reference baseline and inherited regression gates

The corrected-unit Wannier90 plus WPW Hybrid-DG Si64 reference completed on 2026-08-24 with 8 MPI ranks, `OMP_NUM_THREADS=1`, and no time cutoff. It converged the DG Hamiltonian SCF in 66 iterations with density residual `2.50e-9`, band-energy change `8.89e-8`, eigensystem residual `1.75e-15`, electron-count error `4.32e-12`, and symmetry residual `2.28e-14`, and published `overlapping_wannier_occupied.chk`. The affine group order was 1536 = 32 translations x 48 point-co-group representatives with five generators. This is an accepted handoff/SCF reference, not completion evidence for Tasks 1--9 or for weak scaling.

The reference exposed and fixed a geometric-product/pullback mismatch in `symmetrize_dg_distributed_pencil_rows`. Grid pullbacks publish column actions satisfying `D(g h)=D(h)D(g)`. Therefore a traversal that appends `D(generator)` on the right advances labels with `product_table(generator,parent)`, not `product_table(parent,generator)`. Preserve the following regression gates throughout this plan:

1. a noncommuting D3 right-regular anti-representation with a redundant order-sensitive generator must match a direct full-group average at 1, 2, 4, and 8 MPI ranks;
2. a unitary generator set that violates the supplied group relation must fail before operator averaging;
3. relation checking must use bounded deterministic probes and must not retain `N_basis x N_basis x |G|` matrices;
4. the Si64 pencil average must reduce H, S, and density covariance to roundoff and advance beyond the first Hamiltonian assembly; and
5. the converged reference residuals above must remain available as comparison evidence, while spatial-backend acceptance still requires the later material and scaling tasks.

Focused commands already passing for the convention fix are:

```bash
OMP_NUM_THREADS=1 python3 tests/dg/run_dg_overlapping_wannier_fragment_symmetry_mpi.py
OMP_NUM_THREADS=1 python3 tests/dg/check_dg_overlapping_wannier_route.py
```

### Task 1: Prove sparse seed and support feasibility

**Files:**

- Create: `src/gs/dc/dg_spatial_mlwf_seed.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_spatial_mlwf_seed_mpi.f90`
- Create: `tests/dg/run_dg_spatial_mlwf_seed_mpi.py`

**Steps:**

1. Write dense-versus-sparse tests for an explicit generalized (H,S) occupied filter and smooth empty-response filter. Require S Hermitian positive definite, X^dagger S X=I within tolerance, C=XQ back-transformation, and P_C=C(C^dagger S C)^-1 C^dagger S. For a fixture stitching map A require C^dagger S C=(AC)^dagger_core(AC), A(CU)=(AC)U, and A(B_gC)=T_gA(C) within tolerance. Include spectral-bound errors, gap closure, ill-conditioned/indefinite S, map-metric/equivariance failure, degeneracy/irrep closure, S-projector idempotency and S-self-adjointness, P_occ P_empty residual, generalized Ritz residual R_H=HC-SC(C^dagger HC), symmetry-intertwining residual, target-rank revelation, fill growth, and 1/2/4-rank invariance in the design's S-induced norms. Through a pure seed API taking already verified B_g,D_g tables, include a nonintertwining trial Y and verify P_D(Y)=|G|^-1 sum_g B_g Y D_g^dagger followed by complete-irrep S-rank revelation, or complete-block rejection on rank loss.
2. Run `python3 tests/dg/run_dg_spatial_mlwf_seed_mpi.py --ranks 1,2,4`; confirm the missing-module failure.
3. Implement scaled sparse X approximating S^-1/2 with Newton-Schulz/purification and gate positive definiteness, ||X^dagger S X-I||, fill, conditioning, and iterations. Estimate verified spectral bounds for H_bar=X^dagger H X; apply Chebyshev occupied/smooth-empty filters matrix-free. When verified B_g,D_g tables are supplied, construct P_D(Y), rank-reveal and S-orthonormalize only complete verified irrep blocks, and return C; the seed module does not discover D_g. Define P_target with the S-orthogonal projector formula; evaluate R_H in the S^-1 norm through a certified sparse metric solve; and recheck intertwining, rank, occupation, and projector gates before localization. Do not materialize dense eigenvectors, S^-1, or an N_W by N_W transform.
4. Compare with dense reference projectors and measure fill/support growth on increasing replicated fixtures.
5. Stop before production input if seed memory/work is superlinear or filter/projector residual and fill cannot be bounded. Localization-support feasibility is evaluated later after the minimal sweep exists. Otherwise commit with `feat: construct sparse spatial MLWF seeds`.

### Task 1B: Define input and data contracts

**Files:**

- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Create: `tests/dg/check_dg_spatial_mlwf_input.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Steps:**

1. After Task 1 passes, add a failing source/input test for backend, occupied projector and smooth empty-filter settings, reciprocal-shell/metric/group-closure and moment-reliability settings, finite-group reduction policy, a maximum accepted-sweep count admitting several thousand sweeps, convergence-check interval, history-free localization mode (`num_cg_steps=0` semantics), search scale alpha_search=0.2, first finishing scale alpha_finish,1=0.1, hold fraction f_hold=0.70, minimum scale alpha_min, stall window m_stall, absolute/relative stall tolerances, sweep retry count, localization tail-interval and retained/omitted-gradient tolerances, truncation/correction/exterior-coupling tolerances, maximum support/correction envelope and fill limits, active-radius tolerance, envelope-rebuild limits, immutable-anchor imbalance limits, and omitted-Gram error budgets. Compute n_hold=ceil(f_hold n_acc,max) once per sector. Require 0<alpha_min<alpha_finish,1<alpha_search and 0<f_hold<1, reject duplicate scales, require m_stall>=1, nonnegative stall tolerances, and at least one positive stall tolerance. Keep the reference backend as default.
2. Run the input test and confirm failure.
3. Add declarations, validation, broadcast, logging, and compact canonical/block/statistic types.
4. Run the input and route tests.
5. Commit with `feat: define spatial MLWF contracts`.

### Task 2: Build bounded ownership and overlap graph

**Files:**

- Create: `src/gs/dc/dg_spatial_mlwf_graph.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_spatial_mlwf_graph_mpi.f90`
- Create: `tests/dg/run_dg_spatial_mlwf_graph_mpi.py`

**Steps:**

1. Write MPI tests for exactly-once fragment-core ownership, read-only duplicate buffer IDs, immutable seed-anchor control ownership, and R_env=R_center,max+R_active,max+R_correction,max+R_overlap under minimum-image geometry. For local finite-difference, retained local H/S, and partition-transition operations, prove N_O(C_f) is contained in the fragment buffer. Separately compare the existing PoU-weighted nonlocal assembly for a projector spanning several cores with a direct disjoint-core reference; require PoU completeness, single counting, Hermiticity, and equal bounded-block actions without requiring full projector support in one buffer. Move centers and active supports within each envelope component bound without ID/owner migration or graph rebuild; make each component escape separately and require REBUILD. Test all conservative core/support/correction ranks, periodic edges, stable IDs, and equal global edge sets at 1/2/4 ranks. For nontrivial groups build a sharded orbit directory with one canonical representative-to-image transporter and exact T_g actions. For every g,h compare direct gh transport with composed h-then-g transport for canonical IDs, grid permutation, integer periodic shift, phase, source/target core-owner endpoint, B_g, and D_g. For the identity group require no orbit directory or remote symmetry schedule. Assert local storage depends on owned core vertices plus bounded buffer/envelope/orbit metadata.
2. Run `python3 tests/dg/run_dg_spatial_mlwf_graph_mpi.py --ranks 1,2,4`; confirm the missing-module failure.
3. Implement compact vertex/edge arrays, unique core ownership, read-only buffer-source maps, immutable anchor ownership, canonical ordering, adjacency, setup-time two-sided validation of fixed conservative participant sets, and reuse of existing affine maps. Preserve the existing PoU nonlocal ownership and communication path; modify it only if the focused core-reference test identifies a concrete defect. Store two typed records in each orbit-directory entry: the real-space field transporter T_g contains source/target canonical IDs and core owners, exact integer grid permutation, periodic shift, and nonsymmorphic Gamma field phase; the separate column transporter D_g or D_g,c contains column permutation, column phase, and irrep action. Verify their linked composition laws but never place the irrep matrix inside T_g; reject interpolation or path-dependent composition. Build and version the Wannier graph, core/source-buffer rank sets, halo endpoints, periodic images, and update/correction participants from R_center,max, R_active,max, R_correction,max, and R_overlap only. Keep finite-difference/H/nonlocal operator schedules separate. Record maximum possible symmetry-closed correction envelopes so Task 6 can conflict-color their overlaps. Add the identity-group bypass. Do not color yet.
4. Rerun the graph and route tests.
5. Commit with `feat: add bounded spatial MLWF graph`.

### Task 3: Add an exact local pair-spread oracle

**Files:**

- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_localization_mpi.py`

**Steps:**

1. Store each reciprocal vector as an integer Miller column q_I with G_I=Bq_I and test sum_I w_I G_I G_I^T=I for orthorhombic and triclinic cells. Add an explicit dimensional/quadrature test ledger: length conversion occurs once; Delta V=|det L|/(N_1 N_2 N_3) on the uniform periodic grid without a duplicated endpoint; sum Delta V=|det L|. Fingerprint whether each source LCFO/grid array is a physical sample or has sqrt(Delta V) absorbed; require the stitching map A to convert once into physical-sample core W and make buffers exact copies in that convention. Accumulate Z_I, norms, projector overlaps, and complete-field matrix elements only over disjoint authoritative cores with one represented quadrature factor. Prove PoU, Delta V, reciprocal-shell weights, and symmetry multiplicities use distinct variables. Inject duplicate buffers, a missing/doubled represented Delta V, a reapplied PoU factor, a duplicated periodic endpoint, and doubled +G/-G pair weights and require failure. Add unit-rescaling fixtures that preserve dimensionless Z/S and convert Omega/H consistently. For every verified g obtain the exact integer unimodular fractional action A_g from R_gL=LA_g, construct M_g=A_g^T and pi_g from q_pi=M_gq, and use BM_g=R_g^TB only as a Cartesian residual check; reject rounded Cartesian matching. Under T_g T_h=T_gh require M_gh=M_hM_g, pi_gh=pi_h o pi_g, R_gh=R_gR_h, and t_gh=R_g t_h+t_g. Fix (T_g psi)(r)=psi(R_g^-1(r-t_g)) and verify its discrete core-owner action and D_g^dagger Z_I D_g=exp(i G_I.t_g)Z_pi_g(I), including phase sign, across remote fragments. Add real and complex two-orbital fast-path tests comparing z_I,n, Omega_Gamma=sum_nI w_I(1-|z_I,n|^2), wrapped centers, certified local-delta intervals, and analytic directional derivatives with complete-field finite differences; include periodic branch crossing, an initially unreliable |z| that later becomes localized, zero-overlap, metric-complete but group-nonclosed shells, noninteger reciprocal actions, and deliberately invalid shell cases. From pinned official Wannier90 v3.1.0 verify that `num_cg_steps=0` resets to steepest descent and that its trial generator uses `trial_step/(4*wbtot)`. Compare the internal factor 1/(4 W_G), direction, sign, and trial K at 0.5 and 0.2. Compare trial U/Z/delta Omega only for an objective-identical fixture and never require final accepted U equality because Wannier90 uses a parabolic line search while this backend uses bounded monotonic backtracking. Freeze source/objective/normalization/unit/quadrature/shell fingerprints as CI evidence. Add a route test proving the production binary neither links/calls Wannier90 nor gathers orbitals/DMN; no Wannier90 installation is required at runtime. Add a small general-block objective evaluator for later noncommuting symmetry tests.
2. Run the localization driver and confirm failure due to the absent oracle.
3. Implement allocation-free pair evaluation using Z'_I=U^dagger Z_I U on the union of both untruncated active supports and a certified tail interval. Accept locally only when the interval's upper endpoint proves decrease; otherwise run the complete distributed evaluator or reject. Implement deterministic anchor/previous-center branch selection and a bounded small-block objective evaluator. Keep the complete distributed Omega_Gamma, center, and finite-difference routines as reference code; never decide acceptance from a different spread surrogate.
4. Rerun the driver and require numerical agreement.
5. Commit with `feat: add local Jacobi spread oracle`.

### Task 4: Build symmetry-connected generator blocks

**Files:**

- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Create: `tests/dg/test_dg_spatial_mlwf_symmetry_mpi.f90`
- Create: `tests/dg/run_dg_spatial_mlwf_symmetry_mpi.py`

**Steps:**

1. Test inversion, a nonsymmorphic fractional translation attached to a finite point operation, a fourfold center orbit, a two-dimensional site-symmetry irrep, a shared-Wannier noncommuting orbit, and a trivial-symmetry liquid-like fixture. Construct LCFO actions B_g and induced Wannier actions D_g; require B_g^dagger S B_g=S, B_g^dagger H B_g=H, matching multiplication tables after periodic-image phases are included, B_g C=C D_g, and [K,D_g]=0. Explicitly distinguish the geometric product from the pullback anti-representation `D(g h)=D(h)D(g)`. Include the accepted D3 direct-average fixture with a redundant order-sensitive generator, plus a unitary but relation-inconsistent negative fixture; C2 alone is not sufficient. Include a seed that is group-closed but not an intertwiner and require rejection/reconstruction. Force images and field holders onto non-neighbor ranks and verify T_g. Test periodic-boundary support images without turning them into internal translation orbits. Reject incomplete, unequal-occupation, interpolation-only, or unverified internal-translation constructions.
2. Run the new driver at 1/2/4 ranks and confirm failure.
3. Build site-symmetry center orbits and induced D_g, verify every B_g against the atomic structure, integer production-grid permutation, LCFO basis action, H, S, occupations, and P_target, then pass these tables to the Task 1 seed API for P_D(Y), complete-block rank revelation, and B_g C=C D_g verification. Build multiword paths by advancing geometric labels on the side required by the pullback convention and validate redundant-generator relations with bounded probes before any group projection. Only after that construct finite verified-symmetry K_sym components; do not mutate fields earlier. Verify exp(K_sym) commutes with every D_g. Treat supercell wrapping only as periodic support metadata. Never use generator averaging to repair a nonintertwining seed or create primitive-translation FFT metadata, colored translation layers, polynomial/Krylov/Cayley/truncated exponentials, full dense translation paths, or a full `N_basis x N_basis x |G|` representation cache.
4. Rerun the new driver plus `run_dg_overlapping_wannier_symmetry_projection_mpi.py` and `run_dg_overlapping_wannier_fragment_symmetry_mpi.py`.
5. Commit with `feat: build symmetry connected generators`.

### Task 5: Implement atomic orbit updates

**Files:**

- Create: `src/gs/dc/dg_spatial_mlwf_jacobi.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_spatial_mlwf_jacobi_mpi.f90`
- Create: `tests/dg/run_dg_spatial_mlwf_jacobi_mpi.py`

**Steps:**

1. Test synchronous update-color commits across immutable anchor owners and all conservative core ranks, plus separate correction colors over the symmetry closure of touched core/buffer ranks. Exercise the ordered OK/REJECT/REBUILD/FATAL MPI_MAX result: ordinary objective rejection discards the color and continues, envelope escape restores and rebuilds, and non-finite, post-refresh buffer-generation mismatch, or representation/version inconsistency terminates the route. A stale pre-refresh buffer must be refilled normally. Test that checkpoints contain core W plus only previous centers/branches, masks, n_acc, l_base, the bounded same-level rolling queue of at most m_stall+1 accepted objectives, epochs, and core generations; immutable objective/normalization/unit/quadrature fingerprints remain setup provenance, while r_retry remains outside restored payload. Consecutive failures must produce effective levels min(l_base+1,l_max), min(l_base+2,l_max), ...; after l_max only the retry count may advance. REBUILD preserves r_retry and therefore the current effective level, but does not itself commit that level. On later acceptance, scale commit/reset follows the common n_hold rule: before n_hold retain l_base=0, reset r_retry, and return to alpha_search on the next sweep; at or after n_hold commit l_eff into l_base and reset r_retry. A real scale change must clear the old-level queue and seed the new queue from the current accepted objective; a plateau at l_max changes neither level nor queue. No previous gradient, conjugate direction, or localization-CG restart state may exist. After restore, moments, spreads, defects, correction membership, participants' cache state, and buffers must be discarded and deterministically recomputed. REJECT publishes neither field nor metadata; OK swaps both, synchronizes, and refreshes dirty buffers before the next color. Test provisional unitary monotonicity, truncated-shadow moment/tail evaluation, post-correction sweep monotonicity and rollback, block-polar residuals/covariance, and the identity-group bypass.
2. Run at 2/4 ranks and confirm failure.
3. Implement finite point-block or ordinary local pair updates W'=WU into core-only field shadows plus minimal authoritative metadata shadows. At every sweep form K_trial=-alpha_l_eff Pi_G(A[W_s])/(4 W_G) from the current field only. Define alpha_0=alpha_search=0.2 and alpha_l=max(alpha_min,alpha_finish,1 2^(-(l-1))) for l>=1 with alpha_finish,1=0.1; use l_eff=min(l_base+r_retry,l_max) and bounded monotonic backtracking. Keep l_base=0 while n_acc<n_hold; a plateau cannot commit a finishing level in that phase. Do not implement Wannier90's parabolic line search, call Wannier90, or store localization-CG history. Refill stale read-only source buffers before evaluation and require generation agreement afterward; check unitary-candidate monotonicity, then evaluate each truncated shadow directly with certified tails and core-only global reductions. Use one blocking communicator-wide MPI_MAX over OK/REJECT/REBUILD/FATAL, take the specified common action, synchronize after OK/REJECT, and refresh dirty buffers only after an OK core/metadata swap. Implement correction in three steps: read-only preparation, colored application, and one final global exterior check. During preparation form every symmetry-closed M_c/Q_c/q_c plan from the same committed post-update W and require q_d eta_dc^ext q_c<=epsilon_ext for every exterior cluster pair; on failure merge the complete offending orbit, discard the whole plan, and restart preparation before any correction commit. Freeze the passing plan, apply correction colors with shadows, then evaluate the retained-plus-omitted global post-correction exterior bound. Any final failure restores the whole sweep checkpoint. Verify [Q_c,D_h,c]=0. For every image define T_g W_c=W_gc D_g,c, Q_gc=D_g,c Q_c D_g,c^dagger, and D_gh,c=D_g,hc D_h,c including column order and periodic phase; verify T_g(W_cQ_c)=W_gc Q_gc D_g,c and direct/composed paths. Never conjugate Q_c by the real-space T_g. Save/restore core W and minimal replay metadata around both phases and the nontrivial-group remote audit; restore invalidates all derived quantities and newer buffers, then recomputes/refills them from restored cores.
4. Run at 1/2/4 ranks.
5. Commit with `feat: apply atomic symmetry Jacobi updates`.

### Task 6: Add conflict-free sweeps and halo exchange

**Files:**

- Modify: `src/gs/dc/dg_spatial_mlwf_graph.f90`
- Modify: `src/gs/dc/dg_spatial_mlwf_jacobi.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Modify: `tests/dg/test_dg_spatial_mlwf_jacobi_mpi.f90`
- Modify: `tests/dg/run_dg_spatial_mlwf_jacobi_mpi.py`

**Steps:**

1. Build separate update and correction-cluster conflict graphs for the maximum admitted envelope; two correction vertices conflict whenever their maximum possible symmetry-closed correction envelopes overlap. Add tests for coloring, runtime cluster expansion without a new same-color conflict, stale graph epochs, periodic-boundary pairs, centers and support membership moving within the envelope under immutable anchor ownership without rebuild, global sparse orthogonality defect, symmetry-conjugated block-polar correction, exterior cluster coupling that forces symmetry-closed expansion, symmetric truncation, color discard and sweep restore, rank invariance, and absence of N_W-squared storage/messages. Exercise the complete correction transaction: a read-only preparation failure discards the whole plan before commit; a passing plan satisfies q_d eta_dc^ext q_c for every cluster pair; later colors cannot invalidate it silently; the final global exterior check catches accumulated coupling and restores the whole sweep on failure. Put remote images on non-neighbor ranks and verify Q_gc=D_g,c Q_c D_g,c^dagger and D_gh,c=D_g,hc D_h,c including column order/phase; reject an attempted real-space T_g conjugation of Q_c. For a 10,000-sweep budget require n_hold=7,000, alpha=0.2 throughout accepted sweeps 1--7,000 despite plateaus, and then plateau-driven 0.1, 0.05, 0.025, ... . Permit early termination when all convergence gates pass. Test that a safety retry before n_hold may use 0.1 or below, but after its acceptance the next sweep returns to 0.2, l_base remains zero, and the same-level queue is reseeded; after n_hold the same accepted retry commits the smaller level. Require a rolling queue of m_stall+1 accepted objectives at one unchanged scale, update it after every accepted sweep, compare its oldest and newest entries only on a configured convergence-check sweep using freshly evaluated convergence gates, and restore the exact queue after rollback/REBUILD. Include the reported regression trace (fixed 0.2: 917.34 Angstrom^2 at 10,000, best 916.80 Angstrom^2; early staged decay: 965.49 Angstrom^2; fixed 0.2 about 10 Angstrom^2 lower at iteration 500) as controller evidence, not as a universal convergence threshold. Reset the queue after a real scale change or a return from a differently scaled search-phase retry, do not change it at an alpha_min plateau, and prevent double reduction. Require `stalled_at_alpha_min` without a fictitious new scale. Explicitly count n_acc only after a complete accepted sweep. Test the collective outcome partition separately: finite scale-responsive objective/correction SCALE_REJECT restores and advances r_retry; REBUILD restores and preserves r_retry; FATAL construction failures restore if needed and terminate without scale change; audit failure restores and terminates. Any independently identified bounded block-local failure first restores the complete sweep, invalidates shadows/buffers/derived data, and only then collectively rebuilds the W/P partition, graph, orbit directory, participant sets, and epochs or fails closed. Add a global-sweep SCALE_REJECT fixture with several previously committed colors: retry exhaustion must report the entire current occupied or finite-empty sector as nonconverged and must not select an alleged offending pair or block. Require REBUILD to replace the old-epoch checkpoint before restart and fail a fixture that restores it again. Add a long convergence fixture whose stable solution needs more than 1,000 accepted sweeps, identical results after checkpoint replay and at 1/2/4 ranks, and successful convergence with history-free directions. Add a guard fixture that fails if previous-gradient or conjugate-direction state affects an update. Build retained-gradient block row/column bounds and omitted-pair gamma_ij bounds from certified support tails; include a case where every omitted entry passes an entrywise threshold but their row sum fails. Force an envelope escape during a sweep and require immediate checkpoint restore before collective rebuild or failure.
2. Run `python3 tests/dg/run_dg_spatial_mlwf_jacobi_mpi.py --ranks 1,2,4,8 --check-traffic`; confirm the new assertions fail.
3. Localize the complete occupied sector and requested finite empty sector separately. Execute update colors with K_trial=-alpha_l_eff Pi_G(A)/(4 W_G), l_eff=min(l_base+r_retry,l_max), and bounded monotonic backtracking. Prepare all correction clusters/Q/q read-only from that same committed post-update W, merge complete offending orbits and restart the whole preparation until every two-sided bound passes, freeze the plan, then apply its correction colors using the four-state reduction, common action, and post-action synchronization. Evaluate the final global exterior bound and complete post-correction Omega_Gamma and classify one collective outcome. ACCEPT requires monotonicity plus finite orthogonality, truncation-tail, correction, projector, branch, representation, covariance, graph/version, and buffer-generation validity gates; final convergence gates are not required for sweep acceptance. SCALE_REJECT is limited to finite scale-responsive objective or numerical truncation/correction rejection with every construction invariant intact: restore n_acc/l_base/queue, increment r_retry, and retry to the fixed limit; exhaustion reports the entire sector as nonconverged without attributing a block. REBUILD restores n_acc/l_base/queue, preserves r_retry, enlarges/rebuilds the envelope/participants/graphs/colors, advances the graph epoch, replaces the old checkpoint, and retries at the same effective scale. FATAL includes non-finite data, post-refresh generation mismatch, representation/transporter/version inconsistency, PoU/core reference mismatch, and other construction failures; restore if mutation occurred and terminate without changing r_retry. Any bounded block-local feasibility failure restores the complete sweep and invalidates shadows/buffers/derived data before the collective W/P partition, graph, orbit directory, participant sets, and epochs are rebuilt for the explicitly identified complete block or the route fails closed. For a nontrivial group, audit failure restores and terminates.

   On ACCEPT increment n_acc once. Before n_hold, a successful smaller retry resets r_retry but does not commit l_eff: retain l_base=0, reseed its queue, and return to alpha_search next sweep. At or after n_hold, commit such a retry into l_base=l_eff and seed its queue. Otherwise append the accepted Omega and retain at most m_stall+1 entries at the unchanged level. On configured convergence-check sweeps, terminate if all fresh gates pass. If a gate remains open and the full queue meets the plateau test, log without reduction before n_hold; at or after n_hold increment l_base once and reseed only when l_base<l_max. At l_max retain both level and queue and log `stalled_at_alpha_min`. Fail if configured rebuild bounds are exceeded.
4. Run replicated long-sweep tests without a short timeout. Require bounded support/fill, retained plus omitted full-gradient bound, symmetry-projected block derivative, complete per-scale objective histories and plateau decisions, update/correction color counts and collective latency, correction cluster radius/overlap/inverse-square-root iterations, color participant fanout/control bytes, immutable-anchor imbalance, envelope size/rebuild frequency/cost, and checkpoint/audit bytes. Permit several thousand sweeps and report maximum-sweep exhaustion as nonconvergence with residuals, never as success or an automatic switch to localization CG. Report convergence as a symmetry-constrained stationary local minimum. This is the localization feasibility go/no-go gate.
5. Rerun the instrumented driver and commit with `feat: distribute spatial MLWF sweeps` only after the gate passes.

### Task 7: Wire occupied/empty localization into construction

**Files:**

- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `src/gs/main_dft.f90`
- Create: `tests/dg/check_dg_spatial_mlwf_route.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Steps:**

1. Add a failing route test requiring `spatial_jacobi` to consume Task 1 sparse projector/intertwiner seeds, stitch them through the production LCFO-to-core map A, and reject metric-pullback, small-block linearity, or symmetry-equivariance residuals before any spread evaluation. Then require separate occupied and symmetry-closed finite empty blocks, connected commutant generators, and the existing basis/checkpoint interface without spawning Wannier90 or gathering coherent eigenvectors. A required grid symmetry that is not an integer permutation must fail the complete route; an explicitly allowed reduced group must restart from seed construction rather than fall back blockwise to WPW.
2. Run the new test and confirm failure.
3. Add one backend dispatch point. Reuse existing occupation checks; at setup compute C^dagger S C and W_core^dagger W_core as small block Grams without forming A^dagger A, verify the fixture identities A(CU)=(AC)U and A(B_gC)=T_gA(C), then freeze C and the map fingerprint. Hand W=AC to the localization backend as its sole mutable state and record the setup residuals with seed projector/fill diagnostics. Preserve the reference branch outside dispatch; do not reapply A or recheck A(CU) for every update.
4. Run both route tests and `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py`.
5. Commit with `feat: route DC LCFO states through spatial MLWF`.

### Task 8: Expose buffer/support convergence diagnostics

**Files:**

- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Create: `tests/dg/test_dg_spatial_mlwf_buffer_mpi.f90`
- Create: `tests/dg/run_dg_spatial_mlwf_buffer_mpi.py`

**Steps:**

1. Add correctness tests for unique core ownership, pre-refresh stale versus post-refresh mismatch behavior, buffer source IDs/generations, local-stencil neighborhood closure, partition-of-unity identities where used, existing PoU nonlocal assembly versus direct disjoint-core reference, nonlocal single counting, core-only moment/norm reductions, setup-only LCFO-to-core map/metric/equivariance residuals, and buffer-assisted versus direct core-reference local H/S/gradient/Z actions. Check Delta V appears exactly once in every integrated reference, never in a raw stencil application, and never aliases a PoU or reciprocal-shell weight. Also test diagnostic output for the three distinct user studies: fragment-buffer width, active-support radius, and final compression radius.
2. Run at 2/4 ranks and confirm failure.
3. Implement only the single-run validity checks and invariant diagnostic output. Reuse the in-sweep controlled truncation contract for final compression. Do not launch comparison calculations or automatically declare buffer/support convergence; provide seed, core partition, radius, tolerance, center/projector, spread, tail, covariance, local H/Z, energy, and transition fingerprints so the user can compare separately submitted runs.
4. Run the correctness/diagnostic tests at 1/2/4 ranks.
5. Commit with `feat: expose spatial MLWF buffer diagnostics`.

### Task 9: Connect complete rejection to WPW

**Files:**

- Modify: `src/gs/dc/dg_hybrid_wannier_selection.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_projection.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_hybrid_wannier_selection_mpi.f90`
- Modify: `tests/dg/test_dg_hybrid_wannier_complement_mpi.f90`
- Modify: `tests/dg/test_dg_hybrid_windowed_pw_basis_mpi.f90`

**Steps:**

1. Create an orbit with one delocalized member and require whole-block rejection. Freeze accepted core-field W, generate complete induced-symmetry WPW candidates, apply (I-P_W^core) with P_W^core=W(W^dagger W)^-1W^dagger, and rank-reveal only complete blocks before forming B=[W,P]. Require bounded-block agreement of the production PoU overlap/Hamiltonian matrices with disjoint-core references before using the hybrid representation. WPW windows, centers, reciprocal stars, phases, and irrep matrices must remain complete induced space-group orbits. Map retained coefficient-space target blocks once through the verified target-field construction and decompose the resulting physical target span into symmetry-complete bounded-support column blocks X with deterministic canonical normalization/order. Include W/P overlap, rank loss after projection, PoU/core Gram and cross-Gram handoff mismatch, an attempted coefficient-space S action on post-stitching B, singular/ill-conditioned S_X^PoU, forbidden runtime changes of X basis, cross-block residual coupling, bounded window expansion, wrong span, unbounded support, and many individually small omitted G or S_X blocks whose accumulated row/column budget fails. Dense rescaling/change-of-basis fixtures compare with the exact principal angle and require the conservative bound to remain safe or fail closed; they do not require equal numerical bounds for different X bases.
2. Run the three existing MPI drivers and confirm the new scenario fails.
3. Implement the stated construction order: freeze W; generate symmetry-complete P_raw; compute P=(I-P_W^core)P_raw in the core physical inner product; perform complete-block rank revelation/orthonormalization; expand windows within bounds on rank failure; freeze B=[W,P]; assemble S_B^PoU and H_B^PoU through the existing hybrid path; compare bounded blocks with S_B^core and H_B^core; then validate the span. For batched bounded-support physical target column blocks X spanning range(P_target), assemble S_B^PoU=<B|B>_PoU, C_BX^PoU=<B|X>_PoU, and the full sparse block Gram S_X^PoU=<X|X>_PoU including cross-block terms; solve S_B^PoU Y=C_BX^PoU and form G=<X-BY|X-BY>_PoU or its verified Hermitian Schur equivalent. Never apply LCFO coefficient-space S directly to B, W, P, or X. For the deterministic canonical X normalization/order, obtain s_X,min>0 from retained block Gershgorin/row sums minus the certified omitted S_X operator budget. Accumulate retained plus omitted G blocks into G_1_hat and G_inf_hat, set G_2_hat=sqrt(G_1_hat G_inf_hat), and accept only when sqrt(G_2_hat/s_X,min) passes the span tolerance. Gate Hermiticity, rank, s_X,min, conditioning, support/degree/fill, RHS count, both S_B solves' residuals/iterations, row/column work, workspace, traffic, and omitted budgets. Dense fixtures compare this sufficient bound with the exact principal angle and prove rescaled/ill-conditioned targets remain safe or fail closed. Never run a target generalized eigensolver/Lanczos, form S_X^-1/2, solve one global RHS independently for every target column, or form a dense global target Gram.
4. Rerun those drivers plus `check_dg_hybrid_windowed_pw_route.py`.
5. Commit with `feat: hand rejected MLWF blocks to windowed PW`.

### Task 10: Verify hybrid Hamiltonian and initial state

**Files:**

- Create: `tests/dg/test_dg_spatial_mlwf_hybrid_state_mpi.f90`
- Create: `tests/dg/run_dg_spatial_mlwf_hybrid_state_mpi.py`
- Modify only if a focused failure proves necessary: `src/gs/dc/dg_hybrid_scf.f90`
- Modify only if a focused failure proves necessary: `src/gs/dc/dg_hybrid_block_cg.f90`
- Modify only if a focused failure proves necessary: `src/gs/dc/dg_hybrid_generalized_eigensystem.f90`
- Modify only if the PoU/core reference test proves a defect: `src/gs/dc/dg_overlapping_wannier_nonlocal.f90`
- Modify only if the basis handoff itself proves a defect: `src/gs/dc/dg_hybrid_full_cell_operator_adapter.f90`

**Steps:**

1. Test bounded-block PoU/core agreement for WW/WP/PP overlap and local/nonlocal H actions, Hermitian WW/WP/PP H, positive retained S, symmetry covariance, occupied projector, energy, transition/position matrices, bounded block-CG iterations, and existing safeguarded hybrid mixing. Include an excessive-inner-CG oscillation fixture. Add the accepted Si64 reference assertions: full-group pencil averaging reaches roundoff covariance, the first Hamiltonian assembly completes, density and band energy decrease without symmetry drift, and the occupied checkpoint is published on convergence.
2. Run the new driver plus `run_dg_hybrid_block_cg_mpi.py` and `run_dg_hybrid_sparse_operators_mpi.py`.
3. Wire the spatial basis into existing interfaces and retain the established PoU nonlocal path. Change the adapter, nonlocal path, or solver only in response to its focused reference failure; do not introduce a parallel projector-reduction design or replace the established mixing policy.
4. Also run `run_dg_hybrid_scf_mpi.py`, `run_rt_dg_hybrid_metric_solver_mpi.py`, and `run_rt_dg_hybrid_length_gauge_mpi.py`. On Si64 compare against the accepted 66-iteration reference numerically rather than requiring identical wall time or orbital gauges.
5. Commit with `test: validate spatial MLWF hybrid initial state`.

### Task 11: Compare invariant results with Wannier90/full-cell references

**Files:**

- Create: `tests/dg/compare_spatial_mlwf_reference.py`
- Modify: `tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py`
- Modify: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`
- Create fixtures under: `tests/dg/data/spatial_mlwf_reference/`

**Steps:**

1. Add a comparison driver for the explicitly selected Gamma-supercell Omega_Gamma, wrapped center orbits, symmetry intertwining residuals, occupied/empty S-projectors, and local H/Z blocks. Compare invariant subspaces, not signs/phases/member order. When the reference tool reports a different Gamma-only functional, compare common moment matrices and centers and label spreads non-equivalent rather than treating their raw values as equal.
2. Run it for a small insulator at 1/2/4 ranks and confirm missing-reference failure.
3. Generate the smallest physically meaningful corrected-unit reference and record commands and Wannier90 version.
4. Rerun comparison and both existing reference drivers.
5. Commit with `test: compare spatial MLWF invariant subspaces`.

### Task 12: Validate representative Gamma-point supercells without a short timeout

**Files:**

- Modify: `tests/dg/run_si64_overlapping_wannier_gate.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_gate.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`
- Create: `docs/plans/2026-08-22-spatial-parallel-symmetry-mlwf-production-results.md`

**Steps:**

1. Require logs for core coverage/uniqueness, internal length/energy units, cell volume, grid dimensions, Delta V and summed unique-core quadrature volume, objective/normalization/unit/quadrature/shell-weight fingerprints, setup LCFO-to-core metric/linearity/equivariance residuals, fragment-buffer width and generation stalls, local-operator closure margins, PoU completeness, nonlocal single counting, PoU/core overlap and local/nonlocal H residuals, partition identities, reciprocal-shell metric defect and +G/-G weight sum, moment reliability/branch residual, peak seed/localization workspace per rank, graph degree, connected-block sizes/counts, correction-plan restarts/merges, q_c and two-sided/final exterior physical-overlap bounds, buffer-refresh bytes, envelope rebuilds, localization mode=`history_free`, Wannier90 revision and N_W90 normalization fingerprint, alpha_search/alpha_finish,1/alpha_min, f_hold/n_hold and hold progress, m_stall/absolute-relative stall tolerances, l_base/r_retry/effective level, attempt-local versus committed scale, return-to-search events, accepted attempts and objective improvement per level, plateau threshold/decision including ignored search-phase plateaus, stalled_at_alpha_min, n_acc, total attempts, scale-change reason, retries, pre- and post-correction Omega_Gamma with units, retained/omitted full-gradient bounds with units, S-projector/Ritz/intertwining residuals, PoU/core target Gram and cross-Gram residuals, S_X^PoU rank/conditioning/s_X,min, G_1_hat/G_inf_hat/G_2_hat and retained/omitted operator budgets, W/PW ranks, conservative hybrid span bound, energy, and transitions. First make the checker reject incomplete synthetic logs.
2. Run its self-test.
3. Run Si64 with 8 MPI ranks and `OMP_NUM_THREADS=1`, no time cutoff, and no concurrent production run: `python3 tests/dg/run_si64_overlapping_wannier_gate.py --repo "$PWD" --backend spatial_jacobi`. Compare its invariant handoff and SCF results with the accepted 2026-08-24 Wannier90+WPW reference: group order 1536, translation/point orders 32/48, five generators, roundoff post-average covariance, 66 reference SCF iterations, final density `2.50e-9`, band-energy change `8.89e-8`, eigensystem `1.75e-15`, electron error `4.32e-12`, and symmetry `2.28e-14`. Do not require equal individual Wannier gauges or spreads.
4. Run at least one insulating surface/interface fixture with reduced finite symmetry and one insulating liquid-like disordered fixture with trivial symmetry through localization, W/PW selection, hybrid initial-state construction, energy, and transition gates. Neither fixture may rely on an inferred internal primitive translation.
5. Record revision, compiler, MPI/OMP, inputs, phase timings, single-run numerical convergence, memory, and every gate for all three classes. Also record user-provided paired-run evidence for fragment-buffer width, active-support radius, and final compression radius; the implementation and checker report differences but do not submit those convergence runs or choose their widths. On failure, return to the first violated focused invariant instead of adding heuristic recovery. Commit passing evidence with `test: validate spatial MLWF Gamma supercells`.

### Task 13: Establish weak-scaling evidence

**Files:**

- Create: `tests/dg/run_spatial_mlwf_weak_scaling.py`
- Create: `tests/dg/check_spatial_mlwf_scaling.py`
- Create: `docs/plans/2026-08-22-spatial-parallel-symmetry-mlwf-scaling-results.md`

**Steps:**

1. Write a checker requiring near-constant active support/fill, maximum-envelope overhead, graph degree, point-block size, update/correction color counts and collective latency, correction-plan restarts/merges, cluster reach/overlap, q_c and two-sided/final exterior bounds, stabilizer-block size/conditioning/frequency, color participant fanout and control bytes, immutable-anchor load imbalance, envelope-rebuild frequency and cost, minimal-checkpoint bytes, existing PoU nonlocal traffic, PoU/core handoff residuals, localization sweeps, accepted sweeps per scale and plateau counts, WPW and target-block RHS/support per atom, G/S_X degree/fill, s_X,min, row/column work and omitted-operator budgets, workspace, and field/control/audit traffic per rank. Treat a large but size-independent localization count, including several thousand sweeps, as compatible with fixed-accuracy weak scaling; reject growth with system size rather than imposing an artificially small count. Require zero remote-symmetry control/audit traffic for the identity group, O(N) local work over the tested fixed-accuracy range, plus bounded setup map/metric residuals and metric/filter/hybrid-projector conditioning and iterations.
2. Self-test it with synthetic passing and failing data.
3. Run replicated and independently perturbed insulating interface/surface/liquid-like cells at sizes 1/2/4/8, then around 1,000 atoms. Do not use an exact internal primitive translation to manufacture favorable scaling. Advance toward 10,000 and 100,000 only after the prior stage passes.
4. Document setup, graph, sweep, hybrid, and operator timings, tested range, and uncertainty without unsupported extrapolation.
5. Commit with `perf: establish spatial MLWF weak scaling`.

### Task 14: Final review and regression

**Files:**

- Modify: `docs/plans/2026-08-22-spatial-parallel-symmetry-mlwf-design.md`
- Modify: `docs/plans/2026-08-22-spatial-parallel-symmetry-mlwf.md`
- Modify implementation only for demonstrated failures.

**Steps:**

1. Run all new graph, symmetry, Jacobi, buffer, hybrid, reference, route, and 1/2/4/8-rank tests.
2. Search source and inspect counters for forbidden field/eigenvector all-gathers, global dense gauge/DMN, per-pair allocation, per-fragment Wannier90, partial remote-orbit acceptance, untracked or non-covariant truncation, interpolation-based symmetry maps, and unbounded inner CG.
3. Review generalized H/S back-transformation, Ritz residual, and S-projector accuracy; the complete unit ledger and one-time internal-unit conversion; unique periodic grid points and sum Delta V=cell volume; exactly one Delta V in norms/moments/projector overlaps/matrix elements and none in raw stencil application; separation of Delta V, PoU, reciprocal-shell, and symmetry weights; +G/-G half weighting; setup-only LCFO-to-core metric pullback, linearity, and symmetry equivariance; immutable C/map and sole mutable core-field W; exactly-once core ownership and core-only localization integration; read-only generation-matched buffers; local finite-difference/H/S/partition-transition closure; the existing PoU nonlocal path against a disjoint-core reference; PoU/core overlap and Hamiltonian handoff; Gamma reciprocal-shell metric completeness, exact integer group action, field-action phase convention, geometric-product versus pullback anti-representation convention, redundant-generator relation closure, group closure/moment covariance, exact Omega_Gamma deltas, wrapped centers, and full-gradient stationarity; objective-independent N_W90 trial-step mapping and objective-equal-only trial U/Z/Omega comparison, with 0.5 retained only as a pinned Wannier90 fixture scale; B_g C=C D_g, direct-versus-composed real-space T_g and basepointed column-space D_g,c transporter equality, and commutant-constrained finite point exponentials; history-free localization equivalent to `num_cg_steps=0`, alpha_search=0.2 held through n_hold followed by plateau-triggered 0.1→0.05→0.025→... down to alpha_min, same-level stall histories, consecutive retry offsets, absence of conjugate-direction state, and long-sweep convergence; no post-sweep failure attribution to individual blocks; immutable anchor ownership and fixed conservative participant sets; synchronous color core/metadata shadow commit/discard and buffer refill; minimal checkpoint state and deterministic reconstruction of derived data; stabilizer covariance and bounded D_g,c-conjugated block-polar correction; identity-group bypass; rejection of artificial internal translations and nonrepresentable required grid symmetries. Verify that relation diagnostics use bounded probes and never retain full-group dense representation matrices. Verify the state paths separately: SCALE_REJECT restores and changes only r_retry; envelope REBUILD restores then rebuilds; final correction failure restores; FATAL restores after mutation then terminates; block-local handoff restores before collective W/P/graph/orbit/participant rebuild; audit failure restores then terminates without rebuild. Finally review graph validity, global orthogonality, availability of user-run fragment-buffer/active-support/compression evidence, conservative omitted-Gram-aware hybrid span bound with explicit s_X,min subtraction, complete W/PW induced symmetry, Hermitian matrices, energy, and transitions. Confirm separately that localization CG is disabled while sparse metric, hybrid block-CG, and SCF solver controls remain unchanged.
4. Update documentation only with measured limits and deviations. Keep metals and dynamic RT adaptation out of scope.
5. Commit documentation with `docs: finalize spatial MLWF production plan`.

# Internal microscopic screening and W updates

## Experimental projection-free full-q W update

Select `hsex_w_solver='frozen'`, `yn_hsex_update_w='y'`, a positive
`hsex_frozen_rank` and `hsex_screen_step`, and optionally a positive
`hsex_screen_interval`. Use TB-mBJ HSEX RT, `projection_option='no'`,
`nstate=nelec/2` and `yn_hsex_frozen_reference='n'`. The existing `auto`
single-q diagnostic remains separate and does not replace W.

All q sectors and all G/G' entries within `hsex_w_cutoff` are computed from
the current coherent occupied projector and the frozen endpoint propagation
Hamiltonian. No physical conduction-state expansion or instantaneous-state
projection is used. Both eta values share response spaces; every equation
must meet `hsex_frozen_tol`. Capacity exhaustion rejects the update.
The full-q path evaluates a conservative spectral-frame residual upper bound:
it retains the off-frame H defects and source-projection defects and uses
triangle/Frobenius bounds plus the previous roundoff allowance. It does not
discard coherences or change occupations. This avoids dense residual-factor
QR and source/shift-dependent large basis rotations. Reported residuals in
this mode are upper bounds, which can be less sharp than the single-q norm.

The default strict closure is a reciprocal, Hermitian **adiabatic static
snapshot** approximation. Before taking its static part, relative q/-q
reality and anti-Hermitian defects are printed as `HSEX_FROZEN_STATIC_DEFECT` and must not exceed
`hsex_w_eta_tol`. The same tolerance checks assembled W's eta dependence.
The old assembler's positive-response convention receives `-chi_raw/2`:
the new commutator already includes spin and both transition directions.
No second spin or +/-q transition factor is applied.

The experimental insulating mini-BZ head uses the average inverse head on
the shortest finite-q shell (`HSEX_FROZEN_HEAD_MIN_SHELL`). This is a
finite-mesh closure, not analytic optical head/wings or a metallic Drude
limit. Currently Nq must equal Nk, each k dimension must be a power of two,
and only scalar closed-shell CPU k-only MPI is supported. The response
coordinator is a validation implementation; production scaling is unproven.

### Stationary insulating optical limit (opt-in)

Set `yn_hsex_frozen_optical='y'` with the full `frozen` solver. This replaces
the **entire** Gamma mini-BZ cell by the existing analytic insulating
head/wing/body integral, retaining local fields and box form-factor phases.
The default `n` preserves the earlier shortest-finite-q-shell experiment.

Only the small occupied-space Hamiltonian is diagonalized. The projector is
first checked using `||[H,P]||F/||P||F`, in Ha, against
`max(1e-8,hsex_frozen_tol)`. This threshold is a stationarity acceptance
criterion, not a proof of the gap. The subsequent projected Sternheimer
solve must also have positive curvature and converge. The response uses
kinetic, nonlocal pseudopotential and frozen HSEX velocity contributions;
HSEX is included exactly once. Occupied rotations leave the result invariant.
No empty RT orbitals or instantaneous-state occupation projection is added.

At each finite q, responses at eta, eta/2 and eta/4 share H actions. Their
static Hermitian parts are extrapolated to zero eta with
`(chi(eta)-6*chi(eta/2)+8*chi(eta/4))/3`, allowing both linear and quadratic
terms. The estimated extrapolation correction and nested linear-fit
disagreement must satisfy `hsex_w_eta_tol`. Dielectric positivity is checked
again after extrapolation, and the change of assembled W from eta/4 is
also checked. Three samples provide a convergence estimate, not a general
existence proof. Further eta refinement is required for a convergence study.

`HSEX_FROZEN_OPTICAL_STATIONARITY`, `HSEX_FROZEN_EPS_LONGWAVE`,
`HSEX_FROZEN_ETA_CHI_LIMIT` and `HSEX_FROZEN_LIMIT_PASS` identify this path.
Any failed check aborts collectively without publishing or replacing W.
There is no automatic fallback to the minimum-shell closure.

This first limit mode is **stationary and insulating**. A coherent excited
projector or a metallic state must not be treated as an equilibrium insulator.
Its nonstationary optical derivative response and metallic static/Drude
limits remain to be implemented and validated. Static W takes eta to zero
at fixed nonzero q before the long-wavelength continuation; optical Drude
response takes q to zero at finite frequency and is a separate calculation.

After positivity, eta, finiteness and reciprocity checks, a new HXT1-format
`hsex-frozen-updated-NNNNNNNNN.bin` is atomically published with the unchanged
TB-mBJ generating metadata. It contains no fabricated eigenvalues or
projected diagonal occupations. Only then is the active W replaced and the
endpoint exchange action/current/energy refreshed. Failed checks preserve
the active W and abort. Propagated orbitals remain unchanged in number.

For a fixed/updated comparison, first generate a common equilibrium W with
this same response/closure so a change of initial screening method does not
masquerade as a pump effect. Keep zero-field controls when evaluating small
changes. This mode does not claim general two-time nonstationary RPA.

## Coherent electric-impulse optical diagnostic

### Unexcited finite-frequency comparison

`yn_hsex_sternheimer_optical='y'` is a separate **diagnostic-only** RT
invocation: it reads an unexcited GS, calculates the frequency-domain response,
and returns without propagating the physical RT orbitals. Run an ordinary
impulse separately with this flag off. It requires scalar CPU TB-mBJ,
`kli_mode='off'`, `theory='tddft_response'`, `ae_shape1='impulse'`,
`trans_longi='tr'`, uniform full-k weights, k-only MPI and occupied-only states.
Use `hsex_w_solver='auto'`, a positive `hsex_frozen_rank` and tolerance,
`hsex_screen_step=0`, no W preparation/update, and no eigenstate projection.
The existing `nenergy,de` define frequencies from zero through `nenergy*de`,
`hsex_w_eta` is positive damping, and `hsex_coherent_axis` selects the current.
Existing energy/time unit conversions apply. `nt,dt` define only the
finite-window companion table. `hsex_sternheimer_kernel` selects `both`
(default), `independent` (only frozen-potential optical response),
`selfconsistent` (only Hartree/TB-mBJ optical response),
`density` (only density probes),
or `reference` (only static current/contact diagnostics). `yn_fix_func='y'`
excludes the self-consistent kernel; it is incompatible with an explicit
`selfconsistent` selection.

The real-linear occupied tangent equation includes the same native Hartree
and TB-mBJ routines, including kinetic-density/current and explicit-A
dependence, via a central directional derivative. Two derivative steps are
compared. A distributed real rational Krylov representation solves `(s-K)x=b`,
`s=eta-i*omega`. Occupied-only static Sternheimer inverses precondition the
frequency residuals; no empty-state diagonalization is needed. The weighted
full-space defect is evaluated using a distributed real QR factorization.
After convergence, every requested frequency is checked by applying the
original tangent operator to its real and imaginary frequency quadratures.
The older unpreconditioned Arnoldi path remains available to native unit tests.
Native current state/A derivatives are independently checked.
No empty eigenstates are used; physical impulses still propagate on the grid.

`sternheimer-{independent,selfconsistent}.dat` contains raw optical current
responses, including the contact term; no DC correction is applied.
The `*-operator.bin` files store reduced operators, not orbital restarts.
Their `HIO2` header precedes the rank, source/defect norms, contact, real
projected matrix, observation vector, explicit validation residual and optional
real residual factor. The native replay utility also accepts the legacy format.
The `*-time.dat` tables are finite-time companions of the same linear operator;
they permit the same sampling, damping and SALMON cubic window as physical RT.
The spectral time reconstruction is checked against an independent shifted
linear solve before its output is accepted. A right-eigenvector linear solve
handles degenerate eigenspaces without relying on individual left/right pair
normalizations. Frequency and window outcomes have separate PASS/REJECT
markers: a rejected window reconstruction must not be used as a matched-window
comparison even when its frequency solve passed.

For stationary density chi0, a separate k-distributed solver streams only the
needed occupied +/-q projectors and solves rectangular one-electron
Sternheimer equations, one trial potential at a time to avoid a tripled
factorization basis. All three diagonal channels retain independent residual
checks. Its residual is evaluated without replacing the
weighted Hamiltonian defect by a bound on the entire response space.
`sternheimer-density-q0.dat` probes G=100,010,001 at zero momentum transfer;
`sternheimer-density-qfinite.dat` probes G=0,100,010 at q=(1/Nk_x,0,0).
Both +/-q terms and reciprocal folding are retained. These are diagnostic
diagonal elements, not a new W or a full local-field dielectric matrix.

`sternheimer-density-optical-chi0.dat` also evaluates the insulating optical
density head from the stationary velocity/dipole response. Its static
reference is obtained independently by occupied-only Sternheimer solves.
It is an operator conversion from current to density response, not an
arbitrary cancellation of a pumped Drude term. This option is confined to an
unexcited insulating reference and does not extend that conversion to metals.
The old nonstationary W-update gates and production kernels remain unchanged.

### Coherent snapshot diagnostic

Use `yn_hsex_coherent_optical='y'`, `hsex_w_solver='auto'`,
`yn_hsex_update_w='n'`, a positive `hsex_frozen_rank`, and
`hsex_screen_step` with `hsex_screen_interval=0`. This selects an additional
**diagnostic**, without replacing W. Use occupied-only scalar CPU TB-mBJ
HSEX RT, uniform full k weights, k-only MPI, `projection_option='no'` and
`yn_hsex_frozen_reference='n'`. Keep `yn_hsex_frozen_optical='n'`.
The momentum transfer is q=0; a Bloch k=Gamma sample is not required.

`hsex_coherent_axis` selects Cartesian x/y/z as 1/2/3 (default 3).
`hsex_coherent_omega` is an energy in the selected input units (default
0.125 Ha). The probe evaluates zero and this frequency, each at
`eta=hsex_w_eta` and eta/2. It writes `hsex-coherent-optical.dat` with the
complex conductivity and **independent-particle** dielectric response,
separate residuals for both resolvents, and induced-density coefficients
at G=0,(100),(010),(001). These are not a full screened head/wing/body matrix.

With L X=[H,X] and R=(z-L)^(-1), the electric-impulse response is

```
sigma(z) = i/volume Tr(v R [v,R rho] + a R rho)
v = dH/dA, a = d2H/dA2, Im(z)>0
```

The background evolves under the frozen endpoint H. Replacing the inner
`R rho` by `rho/z` would impose stationarity and is not done here. The
contact term includes kinetic, nonlocal pseudopotential and frozen HSEX
derivatives. There is no induced Hartree or induced HSEX vertex in this chi0
diagnostic and no finite-frequency Hermitization. G-shifted occupied vectors
seed the iterative spaces with **zero source coefficients**: rho, physical
occupations, and the residual tolerance remain unchanged.

Each occupied state's actual H action at A+/-1e-4 is compared with the first
and second derivative actions. Relative discrepancies must be below 1e-6
and 1e-4 respectively. The HSEX density kernel is held fixed during this
operator-derivative check, matching the frozen self-energy in chi0. The
uniform induced charge must be below max(1e-8,10*hsex_frozen_tol). All MPI
ranks reject failures collectively before a PASS marker is printed.

The same opt-in diagnostic reports `HSEX_GRID_LOCAL_TRANSLATION` for joint
integer-grid and half-grid translations of the actual local field and each
occupied periodic orbital, along all three lattice directions. The output
gives maxima over k/orbitals of `||V_delta*Tpsi-T(V*psi)||/||V*psi||` and the
absolute change of the local-potential expectation value, in Ha. Integer
covariance must hold within 1e-10. Half-grid defects measure local collocation
aliasing and are only reported; neither states nor response poles are removed.
The fractional shift interpolates a real potential with symmetric Nyquist
coefficients and orbitals with a unitary complex Fourier shift. This component
test does not resample ionic/nonlocal projectors or recompute a translated SCF
state, and does not establish convergence of the complete Hamiltonian.

`screen_q2_head` is a separately tested numerical building block for the next
Gamma closure: nested 26/52-point quadratic fits of `F(q)=q^2 W(q)` and an
equal-volume spherical average `3*A/qc^2+Tr(C)/3`. Here `C` is the quadratic
coefficient in physical Cartesian q, and the first output is A, not W(0) or
an inverse dielectric. In a metal A can vanish while C remains finite.
This is a zero-translation matrix average; a complete real-space kernel
connection must also account for the lattice-translation phase within the cell.
The routine retains full matrix coefficients, rejects inadequate samples,
failed/nested fits and invalid matrices, and preserves the existing regular-W
rectangular-cell fitter's default behavior. It is not yet enabled in the
nonstationary full-W driver; small-q/eta convergence and full-operator
translation tests are still required before publishing such a Gamma cell.

The auxiliary HSEX velocity previously missed the `1/hvol` used in
`hsex_add`. Both the frozen-W and projected-W update callbacks now use the
same tested derivative helper. Static equilibrium tests with zero induced
HSEX could not expose that defect. This fixes auxiliary response operators;
previous nonzero-HSEX screening data require rechecking before reuse.

This diagnostic does **not** define the static Gamma-cell W of a general
coherent/metallic state. In particular, a low-frequency pole already present
in the unpumped finite-grid/k reference must not be identified as a
photoinduced Drude signal. Keep zero-field controls and convergence tests;
do not force the pumped conductivity to vanish or insert the raw optical
dielectric into a static metallic mini-BZ formula. The existing stationary
insulating W-update gate remains intact.

When all k sectors pass a separate 1e-8 Ha stationarity check, an independent
occupied-space/Sternheimer calculation prints
`HSEX_COHERENT_STATIONARY_REFERENCE D_raw/epsilon_IP_static`. It evaluates
the contact plus zero-frequency paramagnetic coefficient and the regular
static IP polarizability without empty eigenstates. This is a diagnostic
reference only: no DC subtraction is applied to conductivity, RT or W.
If using an equilibrium reference subtraction in post-processing, hold that
same reference fixed for the pumped response so a pump-induced low-frequency
term is retained. A two-width fit alone does not establish the static limit.

## Frozen-H occupied-orbital response diagnostic

`hsex_frozen_rank>0` enables an experimental single-nonzero-q diagnostic.
It forms the commutator directly from the occupied real-space orbitals, without
instantaneous eigenstate projection or a conduction-state expansion. The
rank parameter limits each internal response space, not the number of RT
orbitals. The main RT state and active W are left unchanged.

| Variable in `&functional` | Default | Meaning |
|---|---|---|
| `hsex_frozen_rank` | `0` | Positive maximum dimension of each response space; 0 disables this path. |
| `hsex_frozen_qpoint` | `1,0,0` | Integer point on `hsex_w_qgrid`, nonzero modulo that mesh. |
| `hsex_frozen_tol` | `1d-5` | Relative equation residual with a roundoff safety term. |
| `yn_hsex_frozen_reference` | `'n'` | Opt-in full-H spectral validation, restricted to at most 512 grid points. |

Use `hsex_w_solver='auto'`. For GS, select non-SIC scalar preparation with
`yn_hsex_prepare_w='y'`; a converged GS and its usual auxiliary states are
still required. For TB-mBJ HSEX RT, set `projection_option='no'`,
`yn_hsex_update_w='n'`, `hsex_screen_interval=0`, and choose the absolute
`hsex_screen_step`. At this step, endpoint A, Vlocal and the propagation HSEX
operator are frozen for the probe. HSEX is included once per H action.
Closed-shell even electron counts, actual occupations 2/0, uniform k weights,
CPU and k-only MPI are checked.

Three trial potentials use periodic parts with G=(0,0,0), (1,0,0), (0,1,0),
including the k+q reciprocal-lattice fold. Vertices that reach the grid
Nyquist boundary are rejected. Trial potentials and both shifts i*eta and
i*eta/2 share the response spaces and cached H actions. Both directions of
the commutator and spin occupation are included once; no additional legacy
+/-q or spin factor is applied. `HSEX_FROZEN_Q_CHI` is the physical retarded
response divided by Nk*cell volume.

`hsex-frozen-q.dat` contains each k pair, trial and shift, sampled response,
optional spectral reference, sampled relative error, equation residual,
response dimensions, shared pair H-action count and seconds, and optional
relative L2 error of the entire induced real-space density for each k pair,
trial and shift. This is not a k-summed density error or a W error. Shared counts
and times repeat within a pair and must not be summed over rows.
`HSEX_FROZEN_Q_SUM_DENSITY_ERROR` separately reports full-grid relative L2
errors after k summation, for eta and eta/2 in each trial. The k-pair density
is folded by exp(-i G_fold.r) and weighted by Ng/(Nk*volume). A value of -1
means that the optional reference was disabled. `HSEX_FROZEN_Q_TRACE_CLOSURE`
checks that the canonical-G Fourier contraction of this density agrees with
the independently accumulated physical trace. These checks do not establish
an error bound for W or an optical spectrum.
`HSEX_FROZEN_BACKGROUND` reports the weighted occupied-projector commutator
rate and its ratio to eta/2. `HSEX_FROZEN_Q_DIAGNOSTIC_DONE` reports unfinished
sources explicitly: diagnostic completion is not numerical convergence.
`HSEX_FROZEN_Q_STATUS` separately reports diagnostic errors.
`HSEX_FROZEN_RSS_MIB` reports process high-water marks before and after the
diagnostic, maximized over MPI ranks. The latter includes any enabled spectral
reference and prior RT allocations; it is not the incremental solver memory.

The norm calculation uses untruncated Householder QR. Identical Q factors in
the residual are combined algebraically, and only the triangular R factors
are formed. This reduces work without truncating the residual; the original
conservative roundoff bound is retained. Ambiguous near-cancelled
sources reject instead of becoming zero responses. Capacity limits never
silently turn into accepted convergence. The test suite includes independent
small Liouvillian and orbital-time references, a coherent evolving-background
Hartree-RPA comparison, NaN/Infinity, near cancellation, resonant roundoff,
remote callback failure and root output rejection.

This coordinator streams occupied states and delegates each H action to its
k-point owner. It keeps the response spaces on the coordinator and is a
validation path, not a scalable production W driver. No W is published,
no optical head/wings are changed, and no real/reciprocal interaction is
silently constructed from a coherent nonstationary response. General
two-time RPA and single-frequency snapshot RPA remain distinct. See
`docs/plans/2026-11-22-frozen-screening-research.md` and
`samples/hsex_frozen_response/README.md`.

## Nonthermal projection capture

At `hsex_screen_step`, `projection_option='td'` freezes the actual RT endpoint
local potential and A, includes the propagation HSEX operator in the auxiliary
eigenproblem, and projects the occupied real-space RT orbitals. The raw diagonal
occupations feed W; no thermal refill, electron-count repair or rescaling is
performed. Inter-eigenstate coherences are discarded by this approximation.

For repeated internal updates, set `hsex_screen_interval` to a positive number
of RT steps. Events occur at `hsex_screen_step + m*hsex_screen_interval`
(`m=0,1,...`). The default interval 0 preserves the single-update behavior.
Convert a desired physical interval to steps using the chosen `dt`; for
example, `dt=0.0516767166707017` atomic units and interval 800 give 1 fs.
Periodic updates require `yn_hsex_update_w='y'` and `hsex_probe_w='none'`.
Normal restart (`yn_reset_step_restart='n'`) retains the absolute step cadence;
explicitly resetting the restart step to zero also reanchors this schedule.

`HSEX_SCREEN_EVENT` records the step, physical time in atomic units and cadence.
`HSEX_SCREEN_SECONDS` reports projection, W update, post-update exchange refresh
and total event wall seconds. `HSEX_W_UPDATE_SECONDS` splits the internal W
update into preparation, computation and save phases. Each time is the maximum
of local elapsed clock differences across the relevant MPI communicator;
individual phase maxima need not add to the reported total. W timing includes
the existing convergence checks and publication of the updated kernel.

`hsex_projection_capture_tol=1d-8` bounds the missing norm per occupied RT
orbital. `HSEX_PROJECTION_POP` also reports weighted propagated/projected counts,
lost electrons, excited electrons and holes. When capture fails, auxiliary
states double up to `hsex_projection_max_states` (0 selects the smaller of twice
the initial `nstate` and the real-space dimension). New orthogonal guesses are
formed from missing RT components, then deterministic grid vectors, and solved
with the same frozen Hamiltonian. Only the auxiliary guesses are normalized;
RT orbitals and occupations are untouched. The initial `nstate` must exist in
the GS restart; additional auxiliary states need no GS wavefunctions.

CPU scalar k-only MPI/LAPACK supports growth. GPU and distributed eigensolver
growth stop explicitly before modifying W; sufficient initial projection states
retain the existing path. The actual initialized grid bounds the state count,
including inputs specified using `hgs`. At the cap, reject collectively before
W publication. Projection halo/CG scratch, ownership and eigenvalues are rebuilt,
but RT ownership and propagated orbital count remain unchanged. Active-W
checkpoint layouts are unchanged; projection restarts from GS and may regrow
on a later screening step.

Norm capture and eigenpair convergence are separate requirements. A captured
RT subspace does not establish convergence of the highest auxiliary empty
state. Dense screening retains its full-H eigenpair and orthogonality tests;
`HSEX_DENSE_TARGET_REJECT` identifies the failing rank, k index, band or matrix
condition. A k index is an array index, not a Gamma-point label. In particular,
the default even 4x4x4 Si grid has fractional components +/-1/8 and +/-3/8.
Changing the auxiliary count for a weak-field timing case still requires all
retained eigenpairs, raw occupation capture and W acceptance tests to pass.
It does not validate that smaller auxiliary space for a stronger pump.

Production screening is now callable in SALMON. Python tools remain independent
reference/analysis programs and are not launched by the Fortran production path.
The input variables below belong to `&functional`.

| Variable | Default | Meaning |
|---|---|---|
| `yn_hsex_prepare_w` | `'n'` | Generate `hsex-w.bin` after a converged GS calculation. |
| `yn_hsex_head_only` | `'n'` | For finite-temperature PZ SOC GS screening, solve only the inversion-symmetric small-q stencil and apply the usual metallic-head acceptance test. |
| `hsex_w_boxes` | `6` | Auxiliary boxes per axis for GS preparation; must divide each grid dimension. RT inherits the existing W box layout. |
| `yn_hsex_update_w` | `'n'` | Recompute W internally after TD projection at scheduled screening steps. |
| `hsex_screen_interval` | `0` | Repeat internal updates at this interval in RT steps; 0 selects one update at `hsex_screen_step`. |
| `hsex_w_qgrid` | `0,0,0` | Screening q mesh; zero triplet selects electronic k mesh. Each dimension must divide k mesh. |
| `hsex_w_solver` | `'auto'` | Keep existing scalar dense / SOC Sternheimer defaults; `'sternheimer'` selects matrix-free scalar insulating GS screening, `'dense'` explicitly selects the scalar dense reference; `'hybrid'` keeps dense finite-q solves and uses real-space optical head/wings. |
| `hsex_w_cutoff` | `3 Ha` | Dielectric plane-wave energy cutoff, entered in the selected energy units. |
| `hsex_w_eta` | `0.001 Ha` | Positive numerical broadening, entered in selected energy units. RT computes eta and eta/2. |
| `hsex_w_head_tol` | `0.05` | Relative guarded metallic-head fit/nested-stencil threshold. |
| `hsex_w_eta_tol` | `0.02` | Relative eta convergence threshold, including charge and wing channels. |
| `hsex_screening_w` | `'none'` | Independent HXL1 LDA kernel for TB-mBJ HSEX RT; the ordinary TB reference file remains required. |
| `hsex_w_scale` | `1` | Dimensionless diagnostic multiplier of the induced screened exchange; range 0–1. Nondefault values require `kli_mode='hsex'`. |

For a strength diagnosis, keep the same GS and W and set `hsex_w_scale=0.5`
or `0.7` in the RT input. The scale multiplies `-W(P-P0)` once, so the
Hamiltonian, current and induced exchange energy use the same strength.
Microscopic W files, their dielectric diagnostics and the reference density
remain unscaled. Zero removes the induced exchange action while retaining
the configured Hartree and XC updates. New HXC3 RT checkpoints record the
scale and reject a changed scale on restart; legacy HXC2 is accepted at scale 1.

GS preparation supports static KLI, TB-mBJ and unpolarized PZ LDA, closed-shell periodic full uniform
k meshes and k-only MPI. It retains the real-space Hamiltonian, finite-grid box
factors and local-field dielectric matrix. It does not require a large conduction
manifold. Normal GS partial diagonalization still uses its configured empty states.

For a filled gapped manifold, the native path computes the optical head and
wings from occupied-state Sternheimer solves with the analytic real-space
`dH/dk` action. This includes the finite-difference kinetic operator and the
nonlocal pseudopotential, including spinors. An internal RT update also includes
the derivative of the staged projected HSEX operator. The dielectric tensor
uses the full local-field Schur complement. Its tensor and the shortest-q
dielectric estimate are printed separately as `HSEX_W_EPS_LONGWAVE` and
`HSEX_W_EPS_FINITE_Q`.

`HSEX_W_FINITE_Q` records each successfully assembled finite q once on the
communicator root. Columns are the one-based q index, three fractional reciprocal
coordinates, the minimum Cartesian `|q+G|` in bohr^-1, and the corresponding
diagonal element of the inverse dielectric matrix including local fields. The
radius uses the shortest reciprocal representative, so it need not equal the
norm of the printed unwrapped fractional q. Its inverse is a finite-q screening
diagnostic, not the optical macroscopic dielectric constant. Gamma is excluded;
use `HSEX_W_EPS_LONGWAVE` for the optical limit. Head-only metallic runs print
only the finite q points actually computed. For broadened partial occupations,
the record uses the final broadening channel, as does `HSEX_W_EPS_FINITE_Q`.

The Gamma cell uses the existing equal-volume sphere. The analytic head,
wings and constant finite-G body are projected onto the auxiliary boxes and
integrated with successively refined angular and radial quadratures. The
relative change of the projected cell average must be below 1e-7. This is a
long-wavelength cell approximation, whose accuracy still depends on the
electronic k mesh. The optical tensor also controls the interpolation tail.
Partially occupied spinor metals use the guarded metallic closure.

Static dielectric inversions check Hermiticity and use LAPACK `ZHETRF/ZHETRS`.
This avoids the spurious invalid/zero-divide exceptions reproduced in the
complex LU path on macOS, without disabling IEEE exception traps. The dedicated inverse
regression enables invalid, zero-divide and overflow traps around the production
inversion. Reference arithmetic stays outside this region because relaxed
complex division/ABS can itself raise flags on some compilers. The optical
regression also checks an indefinite two-by-two pivot and rejects singular,
non-Hermitian and nonfinite matrices. Hermiticity uses scaled real-component norms, retaining the same 1e-10
relative threshold while avoiding relaxed complex-ABS exceptions at zero.
Finite-q screening remains checked against an independent LU-based reference.

For finite-temperature PZ SOC ground states, `yn_hsex_prepare_w='y'` and
`yn_hsex_head_only='y'` run the static Sternheimer solver only at the finite-q
points in the nested head-fit stencils. Source orbitals remain k-distributed.
Each rank exchanges only the source orbitals for its current q batch.
The selected samples are stored as
`hsex-spinor-head-samples.bin`. A successful test prints
`HSEX_INTERNAL_GS_HEAD_PASS` after normal GS checkpoint writing.
It does not produce a full W kernel. The q mesh and local-field cutoff should
be recorded when comparing these diagnostics.

RT updating currently uses TB-mBJ HSEX, `projection_option='td'`, and enough
projection states to capture the population. `hsex_probe_w` must be `'none'` when
internal updating is enabled. The same staged HSEX action participates in the
projection and Hamiltonian evaluation. Only propagated occupied orbitals undergo
time evolution; projection bands do not become propagated states.

The update writes `hsex-updated-SSSSSSSSS.bin` in HXP1 format, atomically publishes
it after collective validation, then replaces the active pair-owned W. It refreshes
same-step energy/current caches. Output files are not silently overwritten.

Fractional populations or the absence of a demonstrated common filled gapped
manifold require q dimensions >=6. The full box-projected finite-q W is fitted
on two inversion-symmetric stencils; unresolved or nonpositive limits reject the
update. Eta is a numerical regulator, not a fitted Drude scattering rate. The
runtime residual, population, eta and head checks retain the Python thresholds.

RT checkpoints include `hsex-active.bin` (HXC2), containing active W, original
reference density, reference XC energy and GS reference identity. TB-mBJ HSEX
restart restores these before propagation, including with reset time counters.
Older HSEX checkpoints without this information are rejected. As in SALMON's
existing projection restart convention, provide original GS projection data at
`./restart/`; propagated orbitals come from `directory_read_data`.

The driver distributes destination k points and q/-q dielectric ownership, then
retains only owned box pairs for interpolation and transforms. Initial W reads
use owned file slices. One dense grid H is built at a time; no dense H files are
needed for internal updating. Source projection orbitals remain distributed across k ranks. The driver
exchanges the orbitals needed by each q batch and releases the exchange
buffers before solving that batch.

Use 4 MPI ranks and `OMP_NUM_THREADS=12`, `OMP_DYNAMIC=FALSE` for the requested
Fugaku test. BLAS work inside q parallel regions uses one thread per worker.
The mixed GNU/LLVM OpenMP OpenBLAS configuration on macOS uses the existing
vendor-safety policy to serialize q BLAS work; this guard does not apply to Fugaku.

Validation programs are in `testsuites/hsex/test_fortran_screen*.py` with audited
builder `build_fortran_screen.sh`. These Python programs validate native modules;
they are not part of a production SALMON run. Test fixtures distinguish numerical
port equivalence from physical k/q/cutoff convergence.

The complete Si integration test is documented in
`testsuites/hsex/internal_si/README.md`. It checks GS W, projected pump W,
restart and propagation with external Python-reference W at 4 MPI × 12 OMP.
Its small mesh and short trajectory test implementation equivalence, not Si
spectral convergence. The legacy Python oracle retains its finite-q Gamma
rule. When the analytic closure is active, this oracle is compared at every
nonzero q. The optical unit test independently checks Sternheimer against an
explicit conduction-state sum, the local-field finite-q limit, and the
isotropic mini-sphere normalization. The external replacement test reads the
native updated W to check propagation through the file-loading path.
The default integration fixture now uses eta=0.001 Ha;
`--eta 0.01` retains the earlier fixture. Neither the eta convergence check nor
the W/current comparison thresholds are relaxed.

### Small-eta numerical stability

Pump response solves reuse the Hermitian tridiagonal factors as a preconditioner,
then refine against the original grid operator. The operator is the same
lower-triangle Hermitian matrix selected by `ZHETRD('L')`: upper entries are
conjugates of lower entries and the diagonal is real. Existing finite,
Hermiticity and projection-residual checks still run before this normalization.
No eigenstates are anchored and no source occupations are changed.

The residual uses fixed-point matrix slices and BLAS with compensated summation;
the solution retains high and low components. Hermitian Gram entries are
accumulated before rounding away the cancelling anti-Hermitian near-pole part.
Two to ten corrections must stabilize every response entry within
`1e-12*(1+abs(entry))`; otherwise the calculation rejects the update. The final
original-H residual must still satisfy `1e-8`. A large residual of the initial
transformed estimate alone is not a reason to reject a successfully refined
solution; nonfinite initial values and factorization failures remain errors.

Five matrix slices are used at eta/2 >= 0.0005 Ha. Smaller broadenings add slices
to compensate eta^-2 error amplification, with a maximum of twelve; requests
beyond that precision budget are rejected. This is a tested numerical algorithm,
not a guarantee for arbitrary conditioning or arbitrary tiny eta. In particular,
a nonconvergent response is rejected rather than silently accepted after a fixed
number of iterations.

This improvement adds computation and memory: the original H and 5–12 shared
H slices per MPI rank, plus sliced RHS workspaces per q worker. It does not
reduce the cost of the dense response workspaces. Measure memory
and execution time before large Nk production runs. The independent Python
reference uses dense LU and its own compensated matrix-product implementation;
production SALMON invokes neither Python nor that reference implementation.

`test_original_h_refinement.py` compares complex diagonal and off-diagonal
responses with a 90-digit original-H oracle and tests nonfinite/unconverged
rejection. `test_refinement_large.py` uses a dense 2048-point exactly dyadic
Hamiltonian with an analytic spectral oracle. `test_pump_lower_triangle.py`
checks the preserved lower-triangle convention. The audited native `refine`
probe reads a synthetic `oracle.bin` exported by setting `HSEX_ORACLE_FIXTURE`
when running `test_original_h_refinement.py`, allowing the same oracle to be
checked on four MPI ranks and twelve OpenMP threads without mpmath on Fugaku.

Near a diagonal resonance, the residual evaluates `E-H_ii` as an error-free
high/low difference before multiplication. Both parts multiply both parts of
the solution. The real diagonal is excluded before slicing the remaining
matrix, so it cannot dominate the row scale; any imaginary diagonal remains
in that matrix. Only identically zero trailing matrix slices are skipped.
There is no small-element cutoff.

With GNU Fortran, SciPy and mpmath available, CTest registers
`hsex_original_h_refinement`, which builds an audited native library and runs
the 90-digit, diagonal-resonance and dense 2048-dimensional regressions.
`testsuites/hsex/test_original_h_diagonal.py` includes adjacent binary64
energies and a subtraction-low-part check. Full Si and fractional-metal MPI
validation remain explicit integration runs described in the test README.

The small-eta original-H revision has passed Fugaku 4 MPI x 12 OpenMP
validation with the production compiler flags. The Si updated-W difference
from the independent Python reference was 9.87e-13 within the unchanged
mixed tolerance; the fractional-metal full-W difference was 7.36e-17.
Details, measured costs, and the verification-only delayed probe are in
`testsuites/hsex/internal_si/README.md`. These are integration and numerical
checks, not production spectrum convergence results.

### Independent LDA screening for a TB-mBJ reference

For a controlled starting-point comparison, unpolarized `xc='pz'`,
`kli_mode='off'` GS runs support `yn_hsex_prepare_w='y'`. Their `HXL1`
file records the LDA generating XC potential (the scalar metadata field is
zero), rather than pretending to be a TB-mBJ or KLI reference.

In TB-mBJ HSEX `tddft_pulse` and `tddft_response` runs, `hsex_screening_w='path/to/lda/hsex-w.bin'`
selects this independent initial kernel. Default `none` preserves existing
behavior. The ordinary `hsex-w.bin` must still be the TB-mBJ reference W file:
its potential is checked against the actual TB-mBJ GS before the independent
kernel is loaded. Crystal, full k mesh, real-space grid, occupied count,
box layout and translations must agree. Generating metadata and all local
kernel entries must be finite. The TB-mBJ initial density and reference
potential remain the subtraction reference; neither file is rewritten.
The active kernel is included in the normal RT checkpoint. This initial
comparison mode rejects simultaneous pump W updates or external probe W.

Use identical pseudopotentials and screening parameters and retain a manifest
of both generating calculations. The legacy W header does not contain a
pseudopotential fingerprint or q/cutoff provenance, so matching those remains
the caller's responsibility. A small k mesh comparison tests the mechanism;
it does not establish a converged optical peak or a preferred starting point.

The opt-in `testsuites/hsex/test_lda_screening_rt.py` regression uses existing
matched TB-mBJ GS and TB/LDA W files. It checks both RT theories, identical
three-step currents, and the native response spectrum output. All files stay
in a caller-selected fresh scratch directory; W generation is not controlled
by the Python test harness.

## Scalar matrix-free GS screening

Set `hsex_w_solver='sternheimer'` with `yn_hsex_prepare_w='y'` to solve the
scalar insulating GS response using SALMON's real-space Hamiltonian action.
The default `auto` preserves existing results and solver choices. The CG route
requires integer occupations, a common positive indirect gap, converged occupied
eigenpairs and an empty level for the positive-gap check, and the existing residual
checks. Scalar fractional/pumped updates are deliberately not routed through
this new GS-only path. SOC behavior is unchanged. Spin degeneracy, local
fields, optical head/wings, q interpolation and W-file formats are unchanged.

`HSEX_W_PHASE_SECONDS` reports rank-maximum accumulated seconds for dense H
construction, target factorization/checks, optical response, source exchange,
vertex construction, response solves, chi reductions, dielectric/Gamma
assembly, box redistribution/interpolation, Fourier transform and reciprocity
validation. These phase maxima need not add to the overall wall time, because
the maxima can come from different ranks and initialization is not itemized.
The established `HSEX_W_SECONDS` aggregate line remains unchanged.

For the measured small Si case, optical-only real-space action avoids the
expensive dense-H products in the optical CG solve. Use `hsex_w_solver='hybrid'`
to keep dense finite-q solves while selecting this optical action. This mode
still builds a dense H for finite q; only `sternheimer` avoids dense H entirely.
Neither option changes the cutoff, q mesh, tolerances or physical kernel.

The scalar CPU internal RT update also accepts `hsex_w_solver='hybrid'`.
It retains the existing dense finite-q treatment of projected nonthermal
occupations and eta checks; the optical action uses the same endpoint
Hamiltonian and HSEX velocity correction. The compressed `reduced` solver
remains restricted to equilibrium GS preparation. A weak-pulse insulating
timing test does not validate a metallic or strongly excited update: the
existing metallic-head q-mesh and acceptance requirements still apply.

Box form factors reuse the three one-dimensional finite-grid sums within each
call, before forming their Cartesian products. The cache is local to the call;
there is no stored state across geometries or threads. This preserves the
finite-grid definition, phase convention and original product/division order.
The optical head/wings and Gamma-cell quadrature nodes, refinement criterion
and tolerance are unchanged. `fortran_screen_forms.f90` compares the reused
sums with a direct reference for boxes 1, 3, 4, 6 and 12 on an anisotropic grid;
its optional `--benchmark` mode is a local performance check, not a CI timing
requirement.

The insulating zero-limit regression in `fortran_longwave.f90` uses an
independent gapped four-site Hamiltonian. Finite-q density vertices are built
from `exp(i*q.r)` and optical velocity vertices from `i*[H,r]`. It checks the
uniform-density interband response, the O(q^2) head and O(q) wings, quadratic
convergence of their normalized coefficients, and invariance to occupied
components of the vertices. A model current sum-rule check independently
compares the real-space contact term with the static paramagnetic response.
These are model/static-response checks: they do not replace a full Si RT DC
current, f-sum or nonlocal-pseudopotential velocity validation. No zero-DC
constraint is imposed on metallic or fractionally occupied response paths.

Static dense finite-q screening packs up to four nonzero-occupation source bands
into one BLAS projection/target-rotation panel. Additional wide buffers are
bounded to 64 MiB per OpenMP worker (with a one-band minimum). Each band's
tridiagonal solve, residual acceptance and occupation weight remain separate.
Pump/original-H refinement keeps the existing projected RHS arithmetic. The
C ABI and W-file format are unchanged. `q_batch` tests mixed G counts, skipped
and fractional occupations, partial panels, empty source sets, an independent
diagonal-resolvent oracle and 1/2 OpenMP workers. Its optional `--benchmark`
reports timing without a routine CI speed threshold. Local standalone tests use
Accelerate on macOS; SALMON retains its existing GNU/sequential-OpenBLAS
concurrent-call guard. The standalone kernel is not built with that SALMON guard.

Dense finite-q screening now keeps one OpenMP team alive across all bounded q
batches for a destination k. Each worker allocates its private scratch once per
target; original-H refinement slices are shared and prepared once per target.
The q charge vertices and response solves are distributed by the outer q loop.
Worker BLAS remains single-threaded, with the established backend policy restored
on the same thread at exit. The C ABI2 reference entry remains available.

MPI source exchange and collective rejection inside this team are executed only
by `omp master`, the MPI-initializing thread. Multiple-worker dense screening
requires at least `MPI_THREAD_FUNNELED` and the primary thread; insufficient MPI
support or a smaller-than-requested OMP team rejects collectively. Non-MPI
communication retains its dummy implementation. Bounded source buffers are
recycled by the master; workers see new buffers only after the top barrier, and
the q worksharing end barrier protects completed outputs before recycling.
`HSEX_W_Q_PARALLEL regions=... batches=...` records the actual per-target region
count versus the number of bounded q batches on the reporting rank.

Dense `vertices` and `solve` phase counters now report maximum accumulated worker
compute time per target, then the usual rank maximum. Worker setup, shared
refinement preparation and synchronization are included in aggregate response
wall time, not in those compute-only phase counters. Compare the aggregate
`HSEX_W_SECONDS` lines for old/new whole-response performance.

Dense H preparation and optical/grid callbacks remain outside the persistent
team. The matrix-free/SOC adapter shares mutable x/hx/target_info buffers and may
invoke Hamiltonian communication; its q/band loops are not made concurrent by
this change. It retains existing Hamiltonian/grid OpenMP parallelism. Moving
those callbacks outwards requires a separate reentrant Hamiltonian adapter.

Regressions include persistent C-ABI/worker numerical equivalence, 4MPI×12OMP,
uneven5MPI×2OMP, more ranks than k points, empty output ownership, nonfinite/vertex
rejection, missing MPI thread support, and a single rank receiving too few OMP
threads. Local same-backend benchmarks are distinct from Fugaku performance.

### Adaptive-projection regression campaign

`testsuites/hsex/check_nonthermal_projection.py WORK --output verification.json`
checks an explicitly launched scratch campaign. Use the `internal_si/rt.inp`
fixture (Nr8³, Nk=Nq2³, four occupied propagated orbitals, dt0.01au, three
steps, step-1 internal W update, eta0.01Ha), existing eight-state Si GS data,
and the same initial W/pseudopotential for each case. Run each on 4MPI×12OMP.
Keep the GS payload and checkpoints on scratch.

* `reference`: previously validated executable, initial nstate8.
* `captured`: new executable, nstate8, default capture controls.
* `initial5`: new executable, nstate5, default capture controls.
* `cap`: nstate5, capture tolerance1e-10, maximum5; expect failure before W save.
* `grown`: nstate5, capture tolerance1e-10, maximum20; expect5→10 expansion.
* `restart`: restart the grown case's step-1 RT checkpoint, with `./restart`
  pointing to the original GS directory, and retain all functional controls.

These deliberately tighter-than-default tolerances force growth with a weak
kick without conflating the test with metallic small-q closure. The checker
compares unchanged W, continued currents, endpoint A, raw occupations stored
in W against projection output, and occupied-only RT checkpoint counts. Its
checkpoint-count reader targets the GNU four-byte record-marker fixture. This
is a small integration regression; it does not establish a strong-pump screening
convergence result, nor validate Fujitsu/GPU/distributed-eigensolver growth.

## Finite-q PDEP experiment (equilibrium diagnostic only)

The empty-state-free reference is `hsex_w_solver='sternheimer'`, using the
occupied projector and real-space Hpsi. It is not an empty-band-sum reference.
Set `yn_hsex_prepare_w='y'`, `hsex_pdep_rank=-1` for the full-G reference,
`hsex_pdep_qpoint=1,0,0` (integer coordinates divided by `hsex_w_qgrid`), and
`hsex_pdep_reference='absolute/new/reference.bin'`. Positive rank selects the
PDEP diagnostic, while default0 preserves production screening. The diagnostic
returns without publishing/replacing a production W. Run reference and each
rank in a new scratch directory/process, sharing only the small reference file.
Only insulating unpolarized CPU GS with uniform full k weights and k-only MPI
is supported. GS must contain at least one empty state for the observable probe;
only occupied states enter screening, and no empty-state sum is performed.

The reference uses `sternheimer_band` for each G right-hand side. The PDEP
callback solves a single trial potential using the same projected CG and true
residual tolerances (relative1e-10, absolute1e-12 times vertex norm), includes
both q orientations and the same spin, volume and uniform k normalization as
`screen_assemble`. Negative-q G columns are relabeled to match `-q-G`; k+q
folding is explicit. For B=-sqrt(v)chi0sqrt(v), the approximate inverse is
I-U[lambda/(1+lambda)]U†. Construction uses m=min(NG,2*rank+8) fully
reorthogonalized Krylov vectors, with deterministic restarts after breakdown,
then Hermitian Rayleigh-Ritz. This is a fixed-subspace approximation; report the
Ritz residual and `PDEP_RITZ_CONVERGED` separately. Unconverged modes must not be
presented as pure rank truncation error. The original G-basis polarization matrix is not formed or loaded in the
PDEP branch. When the Krylov dimension reaches NG, the Rayleigh-Ritz matrix
is NG by NG and construction needs all NG actions; this is explicitly
reported and provides no RHS-count advantage. The reference file contains only geometry/G order,
six probe vectors and their W actions. It is not a stored dielectric matrix.

Two deterministic complex trial vectors (one is the Krylov seed, the other
is independent of it) and two electron/hole transition-charge pairs
probe the full local-field reciprocal-space W *before box projection*. Each
pair uses the first empty conduction band at k and k-q and a valence band pair;
its direct-kernel element is rho_e†W rho_h / volume (Ha). This is a representative
BSE direct term, not a BSE excitation or an E1 calculation. Compare W action
relative errors, complex direct terms and absolute errors; initial acceptance
is1% for both, with near-zero denominators reported explicitly. Gamma is rejected:
production optical head/wings, analytic velocity response, and metallic closure
are not changed or approximated by the experiment. No symmetry reduction is used.

`PDEP_SOLVES` counts all occupied-band RHS/CG iterations/Hpsi actions in the
measured response; `polarization_actions` is the reference G-column count or
the actual Krylov action count. `PDEP_SECONDS` reports setup, response/basis
construction, inverse/application and diagnostics, and their total, as rank-max
wall times. MPI reduction is included in the corresponding phase. The measured
total excludes SCF but includes basis construction and reference-file validation.
`PDEP_PHASE_SECONDS` additionally separates CG, k-sum/reduction (including
waiting), remaining vertex/response/epsilon/Ritz assembly, and inverse/application.
These are maxima over ranks and need not sum to the maximum wall time because
the maxima can occur on different ranks. The top-level setup/build/application
timings account for the total.
`PDEP_PEAK_RSS_MIB` is the maximum process high-water mark before response
construction and after the diagnostic, including the preceding GS restart/SCF; each variant requires a
fresh process. RSS is not a standalone estimate of dielectric storage, and the
before/after difference must not be interpreted as allocated workspace size.

The single-q experimental adapter currently gathers the occupied manifold plus
one conduction *probe* state once onto each rank. Both reference and PDEP use
this same setup, included in time/RSS. Hpsi and the k sum remain k-MPI-local,
and there is no dense grid Hamiltonian. This diagnostic does not replace the
production driver's distributed source exchange and is not an Nk32³ production
memory-scaling implementation. Interpret the measurements as a controlled
single-q algorithm experiment. If PDEP proves useful, reuse the existing
production source exchange before deploying a production PDEP path.


## Gamma-only production PDEP for insulating supercells

Use `yn_hsex_prepare_w='y'`, `hsex_w_solver='pdep'`, electronic/q meshes
1×1×1, `nproc_k=1`, and `nproc_rgrid=1,1,1`. Source occupied bands are
MPI-distributed through `nproc_ob`; the occupied projector and full real-space
grid remain replicated. Hpsi uses OpenMP. This CPU equilibrium implementation
supports scalar closed-shell insulating PZ LDA and TB-mBJ. It does not extend
SOC, metallic/nonthermal Gamma screening, GPU support or HSEX RT MPI layouts.
No grid Hamiltonian or grid-by-G phase matrix is constructed.

The augmented polarization coordinates contain three optical velocity-response
channels and the nonzero-G body. A time-reversal-adapted real basis combines
G and -G with cosine/sine coordinates and imaginary optical channels. Each
operator action requires one projected Sternheimer solve per occupied source
band. Source response densities are summed before one forward FFT per rank
and trial. The optical response and spin/time-reversal factors use the same
normalization as the independent full-G reference.

`pdep_modes` uses restarted residual expansion with double orthogonalization,
full-support initial probes and four guard modes. A warm subspace therefore
cannot silently pin the solver to an old nonleading invariant subspace.
`hsex_pdep_rank=0` starts at min(2*Noccupied, NG_body+3); a positive rank is
explicit. `hsex_pdep_tol=1e-6` and `hsex_pdep_maxiter=40` govern numerical
Ritz convergence, not truncation accuracy. Increase rank and compare relevant
observables; there is no built-in 1% physical acceptance rule.

The low-rank dielectric is reduced through M=I+F†F, where
F=U_body sqrt(lambda). Optical coupling H=U_head sqrt(lambda) gives the
macroscopic dielectric I+H M^-1 H†. The body inverse and screened wings
are contracted directly with box form factors using BLAS, then integrated
with the same radial/angular Gamma-cell quadrature as `longwave_cell`.
The normal W writer is used. For a full-G reference at Gamma select
`hsex_w_solver='sternheimer'`, `hsex_pdep_rank=0`.

`hsex_pdep_basis_read/write` default to `'none'`. The GPC1 cache contains only
initial vectors and compatibility metadata, never reusable eigenvalues or W.
Current response actions and convergence checks run after every read. Cell,
grid, G order, cutoff, occupied count and rank must match, while ionic positions
may change. Publication uses a new temporary file and exclusive final creation;
existing outputs and identical read/write paths are rejected collectively.

### Verification and example

`samples/hsex_gamma_pdep/water_gs.inp` is a one-water smoke test, not a
converged liquid model. Use 4 MPI×12 OMP, with one BLAS thread, for the default
orbital decomposition. `testsuites/hsex/check_gamma_water.py WORK --output
verification.json` checks separately launched `reference`, `full`, `warm`,
`reduced`, `warm-reduced`, `auto`, and `serial-full` stages. The latter uses
1 MPI×12 OMP. Keep GS/restart payloads on scratch.

Native `fortran_gamma_pdep_mpi.f90` / `test_gamma_pdep.py` compare the FFT
operator to independent real-space vertices and full optical quadrature,
including changed-H cache reuse, empty-source ranks, collective solver/cache
failure, NaN/Infinity rejection and refusal to overwrite a cache. Adaptive-mode
tests include an obsolete exact eigenvector that is no longer dominant.
`fortran_gamma_optical.f90` compares the low-rank optical Schur complement with
the original full-matrix result. Build these through the guarded native builder.
Large amorphous/liquid rank convergence, speedup and distributed-real-space
capacity remain to be measured; the small smoke test does not establish them.


## Fixed-target-k shared response-space experiment

The opt-in `hsex_reduced_target>0` GS diagnostic is independent of Gamma PDEP.
Use the Sternheimer solver, `hsex_pdep_rank=0`, scalar insulating non-SIC GS,
full uniform k mesh and k-only MPI on CPU. Defaults are target0 (disabled),
`hsex_reduced_qcount=4`, `hsex_reduced_maxrank=384`, and
`hsex_reduced_tol=1e-6`. The selected target is a one-based global k index.
The shortest nonzero q points are selected with production mesh enumeration,
metric, G cutoff, source lookup and k+q wrapping. No W file is published.

For fixed target k, let A=Q_k H_k Q_k on the unoccupied complement. Source
energies enter as exact shifts e_j in (A-e_j)x_j=b_j; right-hand sides differ
with source band, q and G. A shared orthonormal basis U gives
x_j≈U(T-e_j I)^-1 U†b_j, T=U†AU. Each added basis vector needs one Hpsi action,
reused across all source energies. Pivoted residual directions are twice
orthogonalized and added in blocks of at most16. The algorithm retains AU,
uses BLAS3 for projection/reconstruction, and diagonalizes only T. There is no
energy interpolation, empty-band sum, dense real-grid Hamiltonian, change to
the physical grid or change to the input orbitals.

The original real-space Hamiltonian is applied independently to every final
solution. Residual normalization is max(norm(Qb),0.01*norm(b),tiny), and the
maximum must satisfy the requested tolerance. Exhausting capacity is status2,
reported as `converged=0` even when the polarization observable is accurate.
For exact arithmetic, the same-energy Galerkin response error has the form
R†(A-e)^-1 R and is positive semidefinite when A-e is positive definite. Thus
observable errors can be much smaller than solution residuals. This identity
is explanatory, not an implemented rigorous stopping bound: no certified
coercivity bound is currently available.

For each q the diagnostic compares sum_n V_n†x_n, a single-target contribution
before global k normalization, q/-q assembly and dielectric inversion. These
are relative errors of target polarization blocks, not errors of complete W,
exciton energies or optical spectra. The two q orientations, if selected,
are independent contributions and are not added a second time. SVD of the
projected RHS reports relative singular-value ranks at1e-2,1e-4,1e-6,1e-8;
this SVD is only a diagnostic and its time is included in setup costs. The
solver itself uses pivoted residual expansion, not those singular vectors.

The full occupied k manifold is gathered on every rank, and selected RHS,
reference and approximate solutions are materialized. Only the owning k rank
runs Hpsi; the other ranks wait for a final collective status. This is a
fixed-target experiment, not a scalable production layout. Shared Q/AQ storage
is bounded by maxrank, but total memory also scales with selected RHS count.
The default production screening and all remote jobs/binaries are unchanged.

### Tests

`fortran_response_space.f90` compares exact shifted solutions for dependent
RHS, complex noncommuting projectors against independent dense positive-definite
solves, capacity exhaustion, the quadratic Galerkin observable-error identity, and
NaN/Infinity rejection.
`fortran_reduced_q_mpi.f90` / `test_reduced_q.py` verify an independent analytic
Fourier response including a nonzero wrapping vector, non-root target ownership,
serial/MPI execution, callback failure and empty cutoff-G collective rejection.
The guarded native builder produces `response_space` and `reduced_q`, and the
existing screening aggregate registers both. Use4MPI×12OMP with1 BLAS thread.


## Symmetry combined with a fixed-target common response space

Set `yn_hsex_reduced_symmetry='y'` within the opt-in fixed-k diagnostic;
default `n` preserves the previous path. Keep the GS full-k and
`yn_symmetry='nnn'`. The first supported subgroup is positive fractional-axis
permutations with zero fractional translation, using the exact grid convention
of the EXX symmetry implementation. Candidate operations preserve the cell
metric, real grid, k/q meshes, species-labelled ions and target k modulo a
reciprocal vector. General signed rotations, nonsymmorphic translations,
time reversal, SOC and global irreducible-target-k W assembly are not included.

For a candidate R, the periodic orbital transforms as
U_k u(r)=exp[i(Rk-k')·r]u(R^-1 r). Occupied covariance is tested on projectors
rather than individual orbital gauges, and C=Psi_destination† U Psi_source
must intertwine occupied energy matrices. Tolerances are1e-6 for normalized
projector residuals and1e-7Ha for energy intertwining. Each nonidentity target
operation also compares H Uv and U Hv for two deterministic complex probes
(relative tolerance1e-7). Sampling the operator is not an exact certificate;
independent all-q CG response restoration below is also mandatory. Rejected
operations/relations simply do not reduce that q; identity remains available.

For q_destination=R q_representative-L, the corresponding reciprocal column
is G_destination=R G_representative+L. The representative raw polarization
matrix is reconstructed using this bijective G map. No extra orbital phases,
star weights or orientation factors are applied to chi for these zero-translation
operations. Only one representative's RHS enters the shared response space;
all selected q remain present in the restored observable. Electronic source
k points related by the operation are handled through the source projector
and energy checks, while the target k itself is held fixed. This is not a
reduction of the complete electronic GS mesh.

The all-q independent CG calculation is retained as a diagnostic reference.
Symmetry-only restoration must match each target-chi block within1e-6 relative
Frobenius norm or the diagnostic fails collectively before any W publication.
Combined approximation errors are reported separately. The existing
`solution/true_residual/status` line applies to representative RHS only;
`TARGET_CHI_ERROR` covers every selected q. Rank-cap status2 is never renamed
convergence, even if observable errors are small.

Timing setup includes covariance discovery, RHS compaction and the diagnostic
SVD. Shared build includes the candidate contractions/restoration AND full-CG
reference chi contractions/error checks, a conservative diagnostic-inclusive
cost. Representative independent-H residual checking is separate. Sum these
three times for the candidate comparison, and include symmetry H_probes in
addition to shared_H for Hpsi accounting. The total-run timer also includes
the all-q independent CG reference. `SYMMETRY_CG` is only accumulated solve
time/counts of representative RHS inside that reference run; it is not a full
separately launched symmetry-only pipeline. Memory still includes all selected
RHS and the replicated occupied k manifold.

Native `fortran_response_symmetry.f90` covers arbitrary degenerate occupied
gauges, broken ionic symmetry, broken H, mismatched source spans, split source
energies, nonzero target reciprocal folds and q_delta/G shifts. The MPI fixture
`fortran_response_star_mpi.f90` / `test_response_star.py` restores an independently
computed all-q Fourier-model response and tests broken symmetries and owner-only
failure on1/4MPI×12OMP. The earlier off-mode fixed-target tests remain registered.


## Adopted full-W response-compression option and RT comparison

`hsex_w_solver='reduced'` selects the all-k/all-q bridge with
`hsex_reduced_target=0`. It uses the existing W assembly and writer. Each owned
target k gathers the occupied source manifold, establishes optional validated
little-group q stars, and batches representative q RHS. Defaults are rank384,
q-batch4 and residual tolerance1e-6. `hsex_reduced_batch` controls the batch;
`yn_hsex_reduced_symmetry` controls the internal q maps, independently of the
full-k GS setting. The source projector remains replicated across ranks;
large-k streaming and global target-k star reduction are not implemented.

Only scalar closed-shell equilibrium insulating CPU GS preparation is allowed.
The optical analytic-velocity head/wings are accumulated by the unchanged
longwave routine before the approximate finite-q branch. A scalar occupied
source contributes V†Q(QHQ-e)^-1QV: the 0.5 Sternheimer factor and scalar
spin2 cancel exactly once, leaving downstream2/(Nk*volume) and q/-q assembly
unchanged. Per-q G relabeling is applied before existing reductions, dielectric
inversion, box projection and reciprocal transform.

Status2 from the common-space solver means rank exhaustion or stagnation. It
is accepted only in this explicit approximate mode, counted as unfinished,
and its maximum true residual is reported; publication is not a claim of
orbital-response convergence. Other errors reject collectively. Per-target
Hpsi counts include basis construction, representative independent residual
checks and symmetry probes. Neither GS states nor RT orbital counts change.

`yn_hsex_reduced_reference='y'` first generates `hsex-reference.bin` with a
local hybrid override from exactly the same in-memory GS, then writes the
reduced `hsex-w.bin`. Input globals are not mutated, and the two files must
have identical metadata preceding the W kernel. The override is confined to
this equilibrium comparison. Repeated source/object audits belong to new
builds, while scratch launches reuse an audited binary by SHA256.

`test_fortran_scalar_screen.py` compares full-rank reduced W with dense and
matrix-free references, symmetry on/off, and verifies finite truncated-rank
output. `test_reduced_w_layout.py` covers replicated/distributed source APIs
with64k on5MPI and2k on4MPI (empty owners), each12OMP. Per-slot status
collectives are used only with globally padded distributed k slots; replicated
uneven ownership reaches a single common final rejection reduction.

`tools/hsex/analyze_reduced_rt.py` is offline analysis only. Its `replay` mode
checks the current reference against completed kick/zero trajectories and
compares its cubic Fourier transform directly with SALMON output. `screen`
requires identical dual-W metadata and matching exact optical tensor, and
records W errors/timings/unfinished response solves. `spectra` requires matched
16384-step kick/no-field trajectories, subtracts no-field current, applies
SALMON's cubic window without sigma0 correction, compares1–6eV Re/Im and the
first stable peak against a75% time window, and writes figures/derived tables.
The plotted Im axis starts at zero. The declared numerical comparison is
peak shift<=0.01eV, peak height difference<=1%, and Im L2 difference<=1%; it
is a same-window approximation test, not a physical convergence criterion.


### Adopted validation preset

The user accepted rank384, representative-q batch4 and internal q symmetry
for the validated Si Nr12³/Nk=Nq4³ static-W RT comparison on2026-11-22.
`samples/hsex_reduced_w_si/` contains reusable fresh-GS, occupied-only RT and
zero-field inputs. Use this preset and the paired finite-window spectrum as
the accuracy baseline for subsequent optimization. The sample disables the
extra reference-W computation during ordinary use; enable it for comparison.
Acceptance is specific to the measured observable/conditions, and does not
turn unfinished response solves into strict residual convergence. The full
source projector remains replicated; current small-grid performance is not
a speedup over the optimized hybrid path. These limitations remain part of
the adopted option's documentation.


### Finite-q density Gamma cell for frozen W

`hsex_frozen_gamma='q2'` selects an experimental density-only closure with
`hsex_w_solver='frozen'`, positive `hsex_frozen_rank`,
`yn_hsex_frozen_optical='n'`, and the existing occupied-only RT update controls.
The default `'shell'` reproduces the existing selection, including the separate
stationary `yn_hsex_frozen_optical='y'` path. The two corrected closures cannot
be selected together.

Each nonzero-q density matrix response is solved at eta, eta/2 and eta/4.
The existing `hsex_w_eta_tol` gates extrapolation of chi0 before dielectric
inversion and the change of W. Then `hsex_w_head_tol` independently gates
nested full-box-matrix fits of q^2 W, using the 26/52 nearest points with all
same-radius ties. Samples on a Nyquist boundary are excluded. The current
power-of-two frozen mesh requires Nk=Nq >=8 in each direction; Nk2/Nk4 is
rejected before Hamiltonian actions rather than silently using a finite shell.

The full spherical Gamma-cell average replaces the q=0 box matrix once,
including its regular term when the singular coefficient vanishes. The old
scalar head is not added. Complex odd wings, Hermiticity, positive dielectric
matrices, fitted positivity, and eta/q convergence remain checked. The q=0
response is not reconstructed from a residual current or by forcing a pumped
Drude term to zero. The source is the frozen occupied projector, without extra
propagated orbitals, instantaneous eigenstate projection or thermal replacement.

Only the selected small-q matrices are gathered on root; the full Nq box
kernel is not replicated for the fit. Final data retain their existing pair
ownership. The finite-q solver still uses its existing owner service, so this
option is not a claim of full large-k performance scaling. A failed check
aborts before a new W is published. Some coherent snapshots need not have an
admissible real static limit and can therefore be rejected.

The spherical cell and the existing discrete-R Fourier convention are retained.
This is not an exact continuous Fourier integral of the complete mini-BZ for
every R. Actual pumped-Si q/eta/grid convergence requires separate measurements;
the native and stationary-Si checks do not establish that physical convergence.

## Streamed nonstationary frozen-W k/q workers

`hsex_frozen_parallel='pairs'` selects local Krylov solves after streaming only
one source and one target k frame to each worker. `hsex_frozen_nproc_q` (default
1) divides the k-only MPI communicator into q groups; the remaining ranks per
group divide the source k sum. For example 8 MPI with q workers2 means2q×4k.
This is a work assignment within the existing communicator, not additional MPI
ranks. It supports empty original k owners, inactive workers and empty final
box-pair owners. Frame preparation and result transfer are collective; all
Krylov H actions are local. Each newly generated H-basis candidate is
orthogonalized once; previously rejected dependent candidates are not revisited.
The legacy expansion order stays available for comparisons. A failed worker
rejects the update collectively.

Frames contain the occupied RT orbitals, frozen local potential, the same A,
rebuilt target-k nonlocal pseudopotential phases, and the source/target HSEX
sigma and projector phases. No full all-k density matrix is replicated. RT
continues to propagate exactly the original occupied orbitals. The original
Hamiltonian and original-equation residual tolerance certify each solve.
Only converged response blocks enter the existing k sum, q/−q closure,
eta checks, dielectric inversion and Gamma-cell procedure. The final W remains
box-pair distributed. Legacy owner-service mode remains the default/reference.

`hsex_frozen_precondition='kinetic'` optionally generates expansion vectors by
applying the positive Fourier inverse `(1 Ha + |G+k+A|^2/2)^−1` to H times each
new basis vector. This changes the trial subspace only, not the response
operator, occupation, shift or convergence threshold. The default is `none`.
The FFT helper supports general Nr, including12; non-power-of-two line lengths
use the existing direct-transform fallback. A nonorthogonal cell and nonzero
k+A are retained. Numerical breakdown and nonfinite callbacks reject the
update. The preconditioner does not guarantee fewer basis vectors or faster
execution; the small512-point Si benchmark reaches full bases with or without
it, so its benefit must be measured at larger grids rather than assumed.

This worker path currently requires scalar CPU, uniform full k weights,
Nk=Nq power-of-two meshes and k-only MPI. Stationary optical closure is kept
on the legacy path; pair workers can use the existing shell closure or the
finite-q `hsex_frozen_gamma='q2'` closure. +U frame transfer is explicitly
rejected because its additional k-dependent operator has not been validated.
The per-q chi matrices and dielectric assembly/inversion still reside on root.
Thus the expensive k/q solves are distributed, but this is not complete
large-Nq memory distribution or a demonstration at Nk32³. Source-pair workspaces
and full trial bases can still dominate memory.

Native tests cover1/4MPI q/k assignments, empty owners, coherent snapshots,
owner-local H and preparation failures, invalid layouts, NaN/Inf preconditioner
rejection, explicit original-equation residuals and the complete q2 Gamma
pipeline. The latter uses a zero-response analytic fixture and is not a
physical pumped-Si convergence test. Small real Si uses8MPI×1OMP/BLAS and the
same nonzero-A1fs snapshot for W comparison. A small shell-closure performance
fixture is not evidence for an experimental E1 peak or optical-limit accuracy.

## Explicit Hermitian static closure for coherent snapshots

`hsex_frozen_closure='strict'` remains the default, including rejection of
non-Hermitian static response above `hsex_w_eta_tol`. The experimental
`'hermitian'` choice defines an additional approximation, rather than raising
that threshold: after the existing q/−q reality check and pairing, use
`chi_H=(chi+chi^dagger)/2` before Coulomb symmetrization and dielectric inversion.
This preserves `Re[x^dagger chi x]` for every trial field of the paired response.
The discarded part is not thereby demonstrated to be negligible or unphysical.
This is a quadratic static surrogate; it is not the full causal response,
a two-time screening theory, or a proven accurate strong-excitation W.

The raw relative defect `||chi-chi^dagger||_F/||chi||_F` is still printed for
all eta samples. `HSEX_FROZEN_STATIC_PROJECTION` records half that defect,
the discarded Frobenius fraction, not an estimated W error. q reality,
nonfinite norms, certified linear-solver residuals, dielectric positivity,
eta convergence/W changes and final W reciprocity still reject updates.
The mode uses eta,eta/2,eta/4, with the existing chi-level extrapolation and
acceptance estimate BEFORE inversion, including when the shell Gamma fixture
is selected. Failed candidates are never published. No automatic fallback
or thermal/eigenstate projection is introduced; four occupied RT orbitals
remain four. The option requires full frozen W and excludes stationary optical
closure. It works with both owner-service and streamed pair workers.

For comparisons, prepare a common equilibrium W using the same closure/eta
rules and feed it to BOTH fixed-W and updated-W runs. The extra eta extrapolation
can change equilibrium W relative to the older finite-eta fixture; attributing
that change to the pump would be incorrect. Gamma/head/grid convergence remains
a separate question; Nk=Nq2³ shell tests do not validate a metallic optical limit.

The real-quadratic preservation and Hermiticity are algebraic properties of
this projection, not claims established by a literature implementation.
TD-aGW/static-COHSEX background: [Many-body enhancement of HHG in monolayer
MoS2](https://www.nature.com/articles/s41467-024-50534-3). That work does not
validate the specific instantaneous coherent projection introduced here.

The initial pumped-Si validation of the explicit projection did not produce an
accepted strong static W: reducing eta worsened the eta-limit indicator, and
finer high-field cases rejected dielectric assembly/inversion. Do not adopt
this option as an established strong-excitation fix. A gapped coherent two-level
model, independently checked against its spectral sum, retains a real1/eta
pole after Hermitian projection. Static-limit existence requires more than
Hermiticity. Details: calc2026/si-static-closure-20261123/report.md.

## Finite causal window on pump-induced differential screening

The experimental `hsex_frozen_window='gaussian'` or `'exponential'` mode
keeps the static equilibrium screening intact and windows only its pump-induced
change. `none` is the default. `hsex_frozen_window_time` is the one-sided
integral area A in the selected input time units (default 0.5 fs). The exponential
is exp(-s/A); the Gaussian is exp(-s²/(2 sigma²)) with sigma=A sqrt(2/pi).
Thus A=0.5 fs means sigma=0.398942 fs for the Gaussian, or eta=hbar/A=1.3164 eV
for the exponential resolvent. Comparing equal areas does not make their
frequency responses identical.

At fresh zero-A RT initialization the driver computes three static-eta channels
and one window channel using shared H-action bases. It caches the paired,
Hermitian equilibrium static and window responses on root. The static W_eq must
match the supplied W0 to relative tolerance at most 1e-6. During propagation:

    chi_eff = chi_eq_static + (Herm(chi_window(t)) - Herm(chi_window(eq)))
    W_out = W_eff + (W0 - W_eq)

The subtraction precedes addition; no-pump updates preserve W0 to roundoff.
Initialization never publishes a replacement W. An incompatible W0 is rejected;
a boxed W alone cannot reconstruct its generating full-G response. The initial
state must pass the equilibrium stationarity check. This implementation requires
full frozen W, Hermitian closure, pair MPI workers, a positive response rank,
and occupied-only scalar CPU propagation. Restart caches, independently supplied
LDA reference W, and the stationary optical option are not supported in this mode.
Failed or partially initialized caches require explicit reset; layout/parameter
mismatch is rejected collectively before H actions.

The window is finite by definition, rather than an eta-to-zero extrapolation at
every update. Only the equilibrium reference uses static eta extrapolation.
`WINDOW_DONE` reports a normalized Duhamel response-error certificate, not a
Gaussian equation residual. `WINDOW_CHI_ERROR_BOUND` propagates absolute raw-chi
Frobenius bounds including the cached reference; it is not a relative delta-chi
or W error estimate. q reality, dielectric positivity, finite values and final
W reciprocity remain mandatory. Root stores the response cache; large-system
memory and scalability still require validation.

The empirical Drude relaxation time 0.5 fs motivates a baseline auxiliary window.
It does not validate that same decay for every interband channel. No collision
term or damping is applied to the propagated occupied RT orbitals, and this mode
is not a transport relaxation model. Width and window shape are physical model
choices to test. Small Nk=Nq2³ shell tests demonstrate implementation behavior,
not converged optical screening, metallic Gamma limits or experimental E1.

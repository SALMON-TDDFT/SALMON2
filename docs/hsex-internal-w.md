# Internal microscopic screening and W updates

Production screening is now callable in SALMON. Python tools remain independent
reference/analysis programs and are not launched by the Fortran production path.
The input variables below belong to `&functional`.

| Variable | Default | Meaning |
|---|---|---|
| `yn_hsex_prepare_w` | `'n'` | Generate `hsex-w.bin` after a converged GS calculation. |
| `yn_hsex_head_only` | `'n'` | For finite-temperature PZ SOC GS screening, solve only the inversion-symmetric small-q stencil and apply the usual metallic-head acceptance test. |
| `hsex_w_boxes` | `6` | Auxiliary boxes per axis for GS preparation; must divide each grid dimension. RT inherits the existing W box layout. |
| `yn_hsex_update_w` | `'n'` | Recompute W internally at `hsex_screen_step` after TD projection. |
| `hsex_w_qgrid` | `0,0,0` | Screening q mesh; zero triplet selects electronic k mesh. Each dimension must divide k mesh. |
| `hsex_w_cutoff` | `3 Ha` | Dielectric plane-wave energy cutoff, entered in the selected energy units. |
| `hsex_w_eta` | `0.001 Ha` | Positive numerical broadening, entered in selected energy units. RT computes eta and eta/2. |
| `hsex_w_head_tol` | `0.05` | Relative guarded metallic-head fit/nested-stencil threshold. |
| `hsex_w_eta_tol` | `0.02` | Relative eta convergence threshold, including charge and wing channels. |

GS preparation supports static KLI and TB-mBJ, closed-shell periodic full uniform
k meshes and k-only MPI. It retains the real-space Hamiltonian, finite-grid box
factors and local-field dielectric matrix. It does not require a large conduction
manifold. Normal GS partial diagonalization still uses its configured empty states.

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
spectral convergence. The default integration fixture now uses eta=0.001 Ha;
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

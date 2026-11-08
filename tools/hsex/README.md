# Real-grid RT-HSEX on a static KLI-SIC ground state

The two modes target weak transverse impulses in three-dimensional, closed-shell
insulators on a full, uniform k mesh. Both keep the initial LSDA + KLI XC
potential and update the Hartree potential from the propagated density. Select
`kli_mode='hsex_hartree'` for this control or `'hsex'` to include the induced
statically screened exchange. Propagation uses the ordinary real-space SALMON
midpoint predictor/corrector, with occupied states only. Neither mode is an
independent-particle approximation. They both include induced Hartree.

Use a converged static KLI restart with extra conduction bands for screening.
The Python preparation tool reads a single native complex-double SALMON
wavefunction checkpoint (single-file, full mesh, native endianness), eigenvalues
in Hartree and KLI metadata. For the 12³ grid with 4³ k points:

    OPENBLAS_NUM_THREADS=1 python3 tools/hsex/prepare_w.py GS RT/hsex-w.bin \
        --mesh 4 4 4 --boxes 3 --cutoff 2 --empty 16

The output must be new; put it in the working directory of the HSEX calculation.
Use the same `directory_read_data` for the two RT runs, four occupied states
for Si, and the usual `tddft_response`, `ae_shape1='impulse'`, `trans_longi='tr'`,
`propagator='middlepoint'`, `yn_predictor_corrector='y'`, integer occupations
and k-only MPI decomposition. A reference-potential, energy, cell, ion and
k-order check rejects incompatible screening files. Only the static screened
interaction and the occupied GS restart are needed for propagation.

The preparation uses the spin-degenerate Adler–Wiser static polarizability,
including both occupation directions at each q, and inverts the complete
symmetrized dielectric matrix within the spherical `--cutoff` (Ha for
|q+G|²/2). The kernel is W = sqrt(v) inv(epsilon) sqrt(v), including off-diagonal
local fields. Box densities are normalized to one and integrated on the SALMON
grid. The q=0 head is replaced by the equal-volume sphere average of 4π/q²
times the inverse macroscopic dielectric function from the shortest nonzero q
star. The angular average of the singular wings is set to zero. This head rule
is for cubic insulators. The finite-q value printed at preparation is not an
extrapolated optical dielectric constant.

For the exchange only, normalized, disjoint piecewise-constant functions on
`--boxes`³ primitive-cell blocks provide auxiliary real-space projectors. This
is a density-density approximation. It is distinct from a Wannier-orbital
approximation and from the microscopic full-grid exchange at a finite block
size. Ordinary orbitals, their density, and the Hartree solve remain on the
fine grid. Check sensitivity to block size, dielectric cutoff, empty-band count
and k mesh before assigning quantitative exciton parameters.

Specifically, define p_k as the block integral of each propagated cell-periodic
orbital multiplied by exp(i(k+A/c)·r). Define D_R = sum_k exp(-i(k+A/c)·R)
p_k p_k† / Nk minus its field-free GS value. The closed-shell exchange energy
per cell is -sum_Rlm W(l0,mR)|D_R(l,m)|². The Hermitian orbital derivative is
used at each Taylor action; the midpoint corrector averages incoming and
predicted D at the same field. Peierls phases couple both the finite blocks
and the lattice translations to A. The exchange part of the printed current
is the explicit field derivative of the same energy at fixed grid orbitals,
divided by cell volume in the SALMON current convention. The reference
subtraction gives zero exchange action on the initial state at zero field.
This is static RPA screening, without a time-dependent W or a qKI quasiparticle
Hamiltonian.

For example, check the nonlocal current with one matched step of 1e-6 au with
and without exchange, using a 1e-4 au kick:

    python3 testsuites/hsex/verify_current.py GS ON OFF

For spectral analysis place the two kicked runs in `hartree` and `sex` and
their zero-field controls in `hartree0` and `sex0`, then use:

    MPLBACKEND=Agg python3 tools/hsex/compare.py RUN_ROOT epsilon

The script uses a common propagation time, subtracts each zero-field current,
and applies the SALMON cubic window. Its output contains the real and imaginary
dielectric functions over 1–6 eV.

## Empty-state-free Sternheimer screening

In the static KLI GS (or a field-free initialization from its frozen restart),
set `yn_out_kli_hamiltonian='y'` in `&functional` in a fresh working directory.
This optional diagnostic exports the exact SALMON grid action at each k to
`kli-hamiltonian-*.bin` by applying hpsi to unit grid columns. For a frozen
restart use zero impulse, a single RT step and occupied states only. The files
belong in calculation scratch. They contain the nonlocal pseudopotential as
well as the kinetic and local terms; an existing file is never overwritten.

    OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/hsex/prepare_w.py \
      GS RT/hsex-w.bin --mesh 4 4 4 --boxes 3 --cutoff 2 \
      --hamiltonians GRID_OPERATOR_DIRECTORY --workers 8

This backend needs SciPy. It uses the occupied GS eigenvalues and wavefunctions
and solves the projected static Sternheimer equations for microscopic Fourier
perturbations. A Hermitian tridiagonal grid reduction is reused for every q,
occupied-energy shift and right hand side at each target k. A shifted
positive-definite solve and its explicit residual are checked before using
the response. Both occupation directions enter the irreducible susceptibility.
The dielectric inversion, local fields, cubic head prescription, box projection
and RT file format are shared with the band-sum backend. This is a **direct
reference solver**: storage for exported H is O(Nk Ng²), reduction costs
O(Nk Ng³). It avoids empty eigenpairs but does not have the large-grid memory
scaling of a matrix-free Sternheimer iteration. Normal GS partial diagonalization
can continue to use empty bands; extra states specifically for W are unnecessary.

## Time-dependent TB-mBJ reference

For the TB-mBJ comparison use `xc='tbmbj'` in GS and both RT inputs. In GS
use `kli_mode='off'` with `yn_out_kli_hamiltonian='y'` and request a restart.
This produces a TB-mBJ-specific screening reference and the grid operators;
prepare W with the same Sternheimer invocation above. The occupied GS
states determine W. SALMON can keep extra bands for GS diagonalization.
Use identical `cval` conventions in GS and RT. With the default negative
value SALMON evaluates c from the instantaneous density.

In RT, `kli_mode='hsex_hartree'` updates TB-mBJ and Hartree;
`kli_mode='hsex'` also updates the induced screened exchange at fixed GS W.
The predictor evaluates the native TB-mBJ potential on the predicted grid
orbitals (including density, kinetic density and current); the corrector
uses the average local potential, followed by the corrected endpoint update.
The HXT1 screening file is distinct from HXR1 for static KLI. A TB-mBJ XC
potential is not a functional derivative of the printed XC energy; assess
the dynamics with current, orthogonality and time-step/field tests.
The extra induced exchange is a chosen response approximation on top of the
TB-mBJ response, not the qKI + HSEX model of Cistaro et al.

For a time-reversal-invariant static reference, add `--time-reversal` to the
Sternheimer command. It verifies the supplied H(k), occupied projectors and
occupied energy operators against their conjugate partners. Only one member
of each matching k/−k pair is solved; the q and G indices of the partner
are reconstructed including zone-boundary folding. Omit the option for
symmetry-breaking references and meshes without explicit opposite k vectors.

During RT the full k mesh is retained, including under a driving field.
An all-to-all transpose redistributes projected density-matrix elements from
k owners to auxiliary-pair owners. Each owner executes the existing HSE
batched FFT plans, with the shifted mesh and vector-potential phases factored
out explicitly. The inverse transpose restores the exchange action to its
original k owner. The lattice part of the current derivative is contracted
with the density in translation space.

## Si pump–probe with screening updated after the pump

Use a closed-shell TB-mBJ GS with extra conduction bands, on the full mesh.
Start every branch from that GS and the same equilibrium hsex-w.bin. SALMON
propagates only occupied orbitals in RT. Set nstate to the number of GS bands
for the separate instantaneous eigensolve. Place the projection step at the
end of the pump, before the first nonzero probe field.

1. Run a pump with theory='tddft_pulse', kli_mode='hsex', projection_option='td'
   and hsex_screen_step set to the last pump step. Leave hsex_probe_w='none'.
   The projection solves the endpoint Hamiltonian with the same frozen endpoint
   exchange operator acting on all trial vectors. Files hsex-pump-*.bin contain
   the instantaneous eigenpairs, projected populations and grid Hamiltonian.
   They are scratch files local to this calculation. Use a fresh directory.
2. Recompute the static box integrals. For example, for an 8×8×8 mesh:

       python3 prepare_w.py GS W-pump.bin --pump PUMP --eta 0.001 --mesh 8 8 8 --boxes 6 --cutoff 3 --workers 8

   The occupations describe a dephased band ensemble. Pump RT and probe RT
   retain the coherently propagated orbitals; only W uses the diagonal
   populations. This version performs one static update, offline. It writes
   one dense grid Hamiltonian per k to scratch, so disk space scales as
   16 Nk (Nx Ny Nz)^2 bytes plus the projected orbitals. The solver includes
   interband and finite-q intraband occupation differences, including both
   directions of a transition. Eta in Hartree regularizes the static retarded
   resolvent. Each k is solved explicitly. Check projection completeness,
   k-point and eta convergence for quantitative carrier screening. Snapshots without a demonstrated insulating manifold require at least six q
   points per axis and use a guarded smooth-metal head (described below).
3. Replay the identical pump in four directories named fixed, fixed0, updated
   and updated0. All four continue to the same final time. Use the same delayed
   weak z impulse as field 2 in fixed and updated. Omit field 2 in fixed0 and
   updated0. For the fixed pair leave hsex_screen_step=0. For the updated pair
   specify the same positive step as in the snapshot and hsex_probe_w pointing
   to W-pump.bin. The update checks the instantaneous energy levels and
   populations before continuing. Use out_projection_step longer than nt to
   avoid unnecessary intermediate diagonalizations. This replay procedure
   preserves the GS reference density for induced exchange across the switch.
4. Run python3 compare_pump.py ROOT STEP OUTPUT to subtract each pump-only
   current and transform the difference using the same probe time window.
   OUTPUT is a filename stem for the dielectric and current tables and figure.

The operator after the switch is the existing induced-exchange approximation
with W replaced in its entirety, i.e. minus W times the difference between
the current and GS density matrices. This comparison includes the action on
the pump-induced density already present when the probe starts. It retains
TB-mBJ and Hartree updates and does not add a separate Coulomb-hole correction.

A concrete Si pump-only input is supplied as `si_pump.inp`. Place the Si
pseudopotential and equilibrium hsex-w.bin in the run directory. The example
expects the matching ground state in ../gs. Its 400-step pump spans 20 atomic
time units. To reproduce the short four-branch check set nt=2448 and use the
field-2 impulse or no probe as described above. The resulting probe window
(2.48 fs) tests the implementation, not the location of E1. Inspect the
omitted population relative to the excited population, not just the absolute
electron count, before interpreting a pump-induced shift.

For a distributed screening run, use a shared input directory and launch the
same command on all ranks, for example
`mpirun -np 8 python3 prepare_w.py GS W.bin --mesh 8 8 8 --q-mesh 4 4 4 --boxes 6 --cutoff 3 --hamiltonians GS --mpi --time-reversal`.
MPI requires mpi4py and one worker per rank. The optional q mesh must divide
the electronic mesh axis by axis. Every electronic k contributes. On a coarse
q grid the code subtracts the box projection of the bare Coulomb interaction
scaled by the shortest-star inverse dielectric head, periodically interpolates
the remainder trilinearly, and restores the analytic term on the electronic
difference grid. The q=0 cell is integrated at the final mesh spacing.
Converge the q mesh for each physical observable and excited carrier density.
The default uses the full difference grid without interpolation.

New equilibrium and pump exports store the lower Hermitian triangle in
Fortran column order (KSH2, TSH2, HSP2). Legacy full matrices are accepted.
Each MPI rank reconstructs only its assigned target matrices for solving.

## Fractional populations, metallic head and Drude response

Pump screening evaluates eta and eta/2 in one target-k reduction and uses the
latter. It rejects a per-q relative Frobenius W difference above
`--eta-tolerance` (default 0.02), nonpositive static dielectric matrices and
unresolved metallic heads. Eta is a numerical regulator, **not** a physical
Drude scattering rate. Taking q to zero at fixed eta suppresses intraband
screening; converge eta at the sampled nonzero q before refining q.

When the snapshot cannot establish a common, completely filled contiguous
occupied band manifold and a positive indirect gap, the full box-projected W
is fitted near zero to
a Cartesian quadratic on two inversion-paired stencils (at least 26 and 52
points). Their fit residuals and cell-average difference must be below
`--head-tolerance` (default 0.05). Matrix, uniform charge and both
charge-neutral wing errors are normalized separately, so large neutral
local-field modes cannot hide unresolved carrier screening. The same channel
checks apply to eta convergence; only numerical roundoff floors are shared.
The polynomial is integrated over the final
reciprocal mesh cell, including anisotropic quadratic terms. Coarse-q
interpolation uses the finite fitted kernel and never restores an insulating
1/q² tail. The insulating rule is used only when all k points share that filled
contiguous manifold and the projection includes an empty level demonstrating
a positive indirect gap. Integer electron/hole pockets therefore use the
finite-q path, not an insulating Coulomb tail. Integer snapshots without empty
projection levels cannot establish a gap and also use the finite-q path.
Diagnostics are written beside the new W as `W_FILENAME.screening.json`.
The finite-q path is a trigger for validation, not proof that a
smooth metallic limit exists: flat or inverted populations may fail.

This is an explicit smooth-limit approximation. Passing nested fits cannot
prove convergence or handle every directionally nonanalytic nonthermal
response. Independently refine electronic k and screening q grids and reduce
eta. On insufficient grids the tool stops without writing a usable W; do not
loosen tolerances merely to obtain a file. The metallic path may need finer
meshes than equilibrium Si and retains an additional W array for eta checks.

RT still propagates the same occupied orbitals. Total current already includes
the nonlocal HSEX phase derivative; no phenomenological Drude term is added.
The projected intraband diagnostic now includes the derivative of that same
instantaneous HSEX operator, holding its density kernel fixed. Projection may
use extra eigenstates but does not add propagated orbitals. Collisionless
carrier acceleration is represented by RT; physical collisions, a dynamic W
and the coherent off-diagonal pump density in the screening solver are outside
this static, dephased approximation. A zero-field/pump-only subtraction and
sufficient probe duration remain necessary to separate the low-frequency
response from interband peaks.

## CPU parallelism and reviewed current implementation

The projected lattice velocity differentiates the retained original R-space
exchange kernel before Hermitianization. Reconstructing that kernel from the
Hermitian k-space action loses finite-field Nyquist terms on even meshes.
`test_lattice_velocity.py --build BUILD --ranks 2` exercises the actual public
initialization/stage/current path against an independent finite difference.

HSEX density products and operator actions use ZGEMM. Independent local-k work
uses OpenMP with private scratch (reused over each thread's k range), capped
at the number of locally owned k points. MPI calls remain outside parallel
regions. BLAS thread limits are scoped and restored using the existing vendor
controls. FFT packing/unpacking and phase multiplication use independent OMP
tiles/rows above a small-work threshold; forward scatter walks contiguous k
FFT rows. The original kernel storage is reused between stages.

On Apple with GNU Fortran and an OpenMP OpenBLAS backend, concurrent BLAS calls
are conservatively disabled: a full Si regression reached end SALMON but
crashed in the vendor's shutdown routine. A minimal reproducer did not establish
the underlying cause; mixed OpenMP runtimes were present. The fallback keeps
BLAS calls on the parent while retaining MPI, FFT and packing OpenMP. Sequential
OpenBLAS builds also avoid concurrent calls because locking is not guaranteed.
Fujitsu and other supported backends retain local-k OpenMP; platform-specific
speed and compilation still require validation. Do not remove these guards
based only on a successful one-step run.

For timing, compare the corrected baseline with optimized OMP1 and OMP4 using
identical inputs, serial vendor BLAS, idle compute resources and repeated runs.
The optional `test_lattice_velocity.py --build BUILD --ranks 2 --benchmark 80
--threads 1,4 --report REPORT.json` isolates stage/update and four Hamiltonian
actions, and saves all final action elements for numerical comparison. It is
a synthetic throughput test, not a physical Si spectrum.

## Distributed W preparation and bounded reuse

The existing `--mpi` option now keeps q/-q partners on a q owner after the
k reduction. That owner assembles/inverts the dielectric matrix and projects
it onto boxes. Only small diagnostics and the two metallic-head stencils are
gathered; eta validation walks individual matrices rather than copying whole
rank-owned W arrays. Head/snapshot root-validation failures are propagated
before the next broadcast. The command-line exception handler still aborts
MPI on other rank-local failures.

A typed, batched all-to-all changes ownership from q to contiguous box-pair
intervals. Pair owners perform the same periodic interpolation and normalized
FFT; metallic cell replacement occurs before the FFT. Large pair arrays use
anonymous temporary mappings and FFT workspace is slab-bounded. Set TMPDIR to
suitable fast scratch on the computing node when available. Mapped pages and
filesystem caches still affect RSS: the workspace limit is not a hard process
RSS cap. The final file is streamed through root in bounded R blocks, checked
for realness/finiteness and R/-R transpose reciprocity, then published without
overwriting an existing output. HXR1/HXT1/HXP1 layout is unchanged. Output I/O
remains streamed root I/O rather than parallel MPI-IO.

For boxes 6³ and electronic mesh 32³, one full complex W is 22.78 GiB. Pair ownership
removes that entire-array requirement from the root: at 512 ranks each pair-owned
array is about45–46 MiB (plus local q matrices, solver work and caches). Aggregate
storage is not eliminated. The current Sternheimer solve still retains all-q
susceptibility accumulators on k workers and uses dense grid reductions;
packed H storage at 12³, Nk=32³ is roughly 729 GiB per snapshot. Matrix-free hpsi-based
Sternheimer remains a separate next step.

Box factors use exact products of finite discrete sums, not continuum sinc
approximations. Integer phases and pump metadata/orbitals use separate bounded
process-local LRUs (64 MiB payload each, at most 512 entries). File identity/change
stamps invalidate pump entries. Root population validation requests metadata
only, avoiding serial orbital reads. No dense Hamiltonians or open descriptors
are cached. Cache capacity can increase RSS for small cases, even though the
large W distribution reduces per-rank memory requirements.

Validation on 2026-11-05: old/new complete W files agree within 1.7e-16 on synthetic
full-pipeline equilibrium/pump tests, including independent q meshes and 3 MPI.
Same 3 MPI, three fresh-process median times on a 6³-grid, k=4³, q=4³ synthetic fixture:
equilibrium 1.048→0.972 s; pump 2.383→1.731 s. A separate 216-box, q=4³→k=8³ poststage
case gives 0.642→0.478 s and median root RSS 510→220 MiB with identical output bytes.
Aggregate RSS was not reduced in that poststage test. Forms-only Si 12³/3 Ha
microbenchmark is 23.44× faster; this is not the full-W speedup. These short local
measurements do not establish Si production or Fugaku performance.

### Numerical-library thread scaling for W preparation

The Python W preprocessor is MPI-parallel; setting `OMP_NUM_THREADS` alone
neither parallelizes its Python loops nor proves that NumPy/SciPy use OpenMP.
Inspect both packages' `show_config()` before choosing library thread controls.
The local Python wheels measured on 2026-11-05 use Apple Accelerate, independently
of SALMON's Fortran BLAS linkage. Accelerate chooses its own parallel work;
`VECLIB_MAXIMUM_THREADS` is a requested ceiling, not evidence of active threads.

Fresh-process, interleaved three-repeat tests compared nominal eight-core
configurations with MPI binding disabled. A synthetic 12³-grid operator,
k=q=2³, 6³ boxes and 3 Ha cutoff gave these median full-W preparation times
(excluding process startup and imports):

| MPI ranks × library thread ceiling | Equilibrium W | Pump W |
| --- | ---: | ---: |
| 8 × 1 | 2.833 s | 2.869 s |
| 4 × 2 | 2.998 s | 3.250 s |
| 2 × 4 | 4.532 s | 5.011 s |

Summed process CPU time / wall time was approximately 7.9, 4.0 and 2.0,
respectively, including MPI waiting. Additional library threads were not
materially utilized by this workload. Prefer MPI-first (8 × 1) for this tested
local environment; these are synthetic performance tests, not Si spectra or
Fugaku scaling results. The test explicitly set Accelerate's calling-thread
model to single-threaded for ×1 and automatic for ×2/×4 via BLASSetThreading,
and queried BLASGetThreading on every rank. Environment ceilings were set before
imports. On Apple this API control is needed to reproduce the strict ×1 test;
the generic OMP environment is not equivalent. Other providers require their
own controls. Re-measure with the actual Fugaku libraries and larger k/q meshes
before selecting a production MPI/OpenMP balance.

### Optional outer-q OpenMP backend

`python3 tools/hsex/build_q_openmp.py /absolute/new/build-directory` builds a
separate shared library, audits source and application objects, and writes its
SHA256 to `audit.json`. Existing build directories are refused. The Python
loader verifies that hash before loading the library. Keep the record beside
the library. This does not rebuild or replace SALMON. On macOS the helper uses
GNU Fortran/OpenMP and Accelerate (BLASSetThreading requires macOS 15 or later).
On other platforms set `FC`, `Q_OMP_FLAGS`, and `Q_BLAS_FLAGS` for an OpenMP
compiler and a **sequential, reentrant LP64** BLAS/LAPACK. External application
object arguments and linker response files are refused by this helper. Custom
library flags must identify vendor BLAS/LAPACK libraries only. Fujitsu cross-compilation was verified on 2026-11-05 with `frtpx` 4.12.2,
`Q_COMPILE_FLAGS="-O3 -Kpic -Cpp -Nalloc_assign"`,
`Q_OMP_FLAGS="-Kopenmp -Nfjomplib"`, and `Q_BLAS_FLAGS="-SSL2"`
(the sequential BLAS/LAPACK variant). Source and object prohibition audits passed.
Compute-node validation job 53643551 exited before Python initialization: the
Spack volume was not mounted. Diagnostic job 53644632 confirmed that
`/vol0004/apps` was unavailable. Fugaku jobs using these Spack packages must
include `#PJM -x PJM_LLIO_GFSCACHE=/vol0004`. The corrected 4 MPI × 12 OpenMP
validation was submitted as job 53645828 (queued on 2026-11-05);
compile success alone is not production qualification. On Fugaku, the user
requires **at least 12 OpenMP threads per MPI rank**: use `--q-threads 12`
(or higher), allocate/bind at least that many cores per rank, and verify the
reported native team size. The wrapper rejects a team smaller or larger than
requested, including an OpenMP runtime cap. Keep `OMP_DYNAMIC=FALSE` and ensure `OMP_THREAD_LIMIT`
is not below the requested team size. Mac 8×1 / 4×2 / 2×4 measurements are not
Fugaku launch recommendations. Begin Fugaku scaling with 12 q threads fixed
and vary MPI ranks; then compare higher thread counts with matched resources.

Add `--q-kernel /absolute/build-directory/libhsex_q.dylib --q-threads 4` to
`prepare_w.py` (Linux library suffix: `.so`). The backend is optional; omitting
it retains the Python/SciPy reference. It requires a Sternheimer backend and
`--workers 1`, and may be combined with `--mpi`. For example use two MPI ranks
with four q threads per rank. Disable rank binding to a single core or provide
an explicit rank placement with four cores per rank. Limit vendor BLAS threads
before Python starts to avoid oversubscribing the sequential Python stages.
Accelerate's native worker setting is additionally forced to single-threaded
and restored per calling worker. This is genuine outer-q OpenMP, distinct from
increasing BLAS's thread ceiling.

The target-k tridiagonal reduction and rotation are computed once and shared.
Python reads snapshots and packs normalized source orbitals/phases in batches;
one OpenMP region spans each batch's q loop. Workers allocate private projection,
rotation, solve, and Gram scratch once per batch and reuse it for assigned q.
MPI calls stay outside that region. The packed-input budget is 64 MiB, except
for one indivisible q item. This is **not** an RSS cap: source arrays, packed
copies, output, shared target matrices, caches and per-thread scratch coexist.
The batch has at most twice the requested thread count. Packing/I/O, target-k
reduction, and remaining dielectric/FFT stages are not made OpenMP-parallel by
this option, so full-W speedup depends on their share of work.

For native validation set `HSEX_Q_KERNEL` to the audited library and run
`testsuites/hsex/test_q_openmp.py`, `test_pump_screen.py`, and
`test_prepare_distributed.py`. The latter two also support `--mpi`; they check
fractional intra/interband response against complete band sums and whole W
files against the reference, respectively. The optional native unit test skips
explicitly if this environment variable is absent. Residual, positivity,
finite-value, eta and metallic-head checks remain enabled.

Measured on 2026-11-05: at 1728 grid points, the isolated synthetic q kernel
(32 q, 16 G, 4 bands) improves from 0.990 to 0.325 s for equilibrium and
1.093 to 0.363 s for paired-eta pump response at 1→4 threads (three-run medians).
This is a kernel speedup, not full W. In 48 full-W runs on the 6³/12³ synthetic
fixtures, the original 8-MPI reference remains fastest. All full-W outputs agree
within 4.21e-13. Native unit, full-W and fractional-response numerical tests also
pass at 12 threads on Mac, including rejection of a runtime-capped smaller team.
No Fugaku performance claim follows from these local tests. Tests accepting
`HSEX_Q_THREADS` use 2 by default; set it to 12 for the full-W/fractional checks.

The fixed-core comparison above must not be used to conclude that adding OMP
at fixed MPI count is ineffective. A subsequent fixed-1-MPI test of the same
native backend (three fresh-process medians, 1/4/12 q threads) gave full-W times:

| Synthetic fixture | Mode | OMP 1 | OMP 4 | OMP 12 | 1→12 speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| grid 12³, k=q=2³ | equilibrium | 7.893 s | 6.003 s | 5.661 s | 1.39× |
| grid 12³, k=q=2³ | pump | 8.923 s | 6.285 s | 5.795 s | 1.54× |
| grid 6³, k=q=4³ | equilibrium | 3.165 s | 1.459 s | 1.111 s | 2.85× |
| grid 6³, k=q=4³ | pump | 5.302 s | 2.145 s | 1.421 s | 3.73× |

All 36 files have zero numeric difference from their one-thread native reference.
This comparison adds cores rather than exchanging MPI ranks for threads. It was
run on an 18-core Mac; it does not establish Fugaku scaling. The 12³ fixture has
only eight q tasks, limiting usable q concurrency, and its target reduction is
still serial within each rank. For Fugaku, retain the intended MPI decomposition
and add at least 12 q threads per rank with matching resources.

### Original-H response refinement (ABI 2)

The native q library now exports `hsex_q_abi()` returning 2. Rebuild older native
libraries; a matching recorded SHA alone does not establish API compatibility.
The loader checks both. The pump response refines the existing lower-triangle
Hermitian operator and retains fractional populations; it does not replace H
with supplied approximate eigenstates. See `docs/hsex-internal-w.md` for the
convergence checks, memory/time cost, and independent high-precision tests.

The historical job 53645828 mentioned above subsequently failed while loading
the Python shared-library backend (`zptsv_` resolution). It is not a validated
execution path. The SALMON-internal implementation links the native solver
directly; its validation and recorded executable hashes are separate.

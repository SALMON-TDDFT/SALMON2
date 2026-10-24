# Ni primitive-cell Mac comparison

User requests: local Mac test, fcc primitive cell, grid width roughly Rps/3
with allowance for skew axes, electronic temperature300K, laser1.55eV/10fs.

Use the generated Ni SOC+NLCC V2 candidate unchanged for both functionals.
At a=3.52Angstrom the primitive columns are (0,1.76,1.76),
(1.76,0,1.76),(1.76,1.76,0). Rps_min=2.8bohr gives Rps/3=.493899Angstrom.
An8^3 grid has axis spacing .311127Angstrom and maximum singular grid step
.44Angstrom;10^3 gives .248902 and .352Angstrom. First use2^3 k points and
24 spinor states for16 valence electrons. This is a bounded test, not a
k-point-converged Ni prediction. Check a finer mesh after stable propagation.

Keep GS checkpoint data transient in /private/tmp; do not collect GS data.
Only RT comparison, inputs, numerical validation and execution provenance go
to canonical calc2026. Earlier580K trials are superseded, not final results.
Each functional must reach a self-consistent initial state before its RT run.
An LSDA checkpoint may seed target-MGGA SCF; direct LSDA-to-MGGA RT would instead
be a functional quench and is not the requested own-state response.

The first real-material test exposed that raw orbital-derived gKS coefficients
bypass ordinary rho mixing. Add SCF-only linear damping of the full local,
tau and current matrices with mixing%mixrate. Preserve instantaneous physical
energy/band/current diagnostics by applying damping only before solve_orbitals;
refresh physical coefficients before energy and final output. Require full
operator residual as well as the ordinary density residual. RT remains fully
instantaneous and unmixed except for its time midpoint predictor/corrector.
Verify tau-only fixed point damping and a charge-invisible spin-current residual.

Build optimized Mac CPU serial and MPI executables in new temporary build dirs;
retain mandatory source/configure/build/object prohibition audits and record
SHA256 before use. MPI decomposition is by k point within the Mac. Validate
numerical component tests before material calculations. Use identical pulse,
geometry, occupations, discretization and propagation settings; compare a
matched no-pulse baseline. Default test intensity1e11W/cm^2, linear x
polarization, dt=.0005fs, eighth-order Taylor midpoint predictor/corrector.
A six-order dt=.001fs attempt failed norm preservation near1fs and is excluded.
The no-pulse stability/half-step test precedes acceptance of the long run.
No physical demagnetization claim before assessing baseline drift and
SCF/grid/timestep limitations.

The real off-axis LSDA test also exposed an existing simple_dm error: after
diagonalizing the local spin density, it overwrote eigenchannel populations
with laboratory-frame diagonals, discarding transverse polarization. Remove
that overwrite; test a pure transverse spinor whose density eigenvalues are
1,0 although laboratory diagonals are1/2,1/2. For conventional CPU simple_dm
SCF additionally require the full spin-density-matrix fixed-point residual;
undo its exact linear alpha suppression after history initialization. This
prevents a charge-only convergence decision from hiding spin relaxation.
The full residual helpers reject NaN and Infinity with scalar IEEE checks
and collective rejection. Both final states use a bracketed initial native-output spin direction near (1,2,1)
and subsequently relax freely with no constraint; this finds stationary states,
not a global-minimum proof.

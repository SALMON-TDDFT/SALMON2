# Unexcited Si: matched impulse and finite-frequency Sternheimer diagnostic

The user approved comparison at identical discretized H, Hartree/TB-mBJ
linearization and current operator, followed by independent-particle chi0.
The known artificial branch moves down in energy on grid refinement.

Implement a real-linear, distributed Arnoldi representation of the occupied
tangent equation d(delta psi)/dt=K delta psi+b A. K includes the central
directional derivative of the native density/Hartree/TB-mBJ routines. Solve
(s-K)x=b at complex frequencies and verify full Arnoldi residuals. No empty
states enter the response or occupied-only reference propagation.

Run the same native Taylor propagator with predictor/corrector, current and
potential routines as RT, entirely inside a converged GS calculation. Keep
four occupied states after the eight-state GS solve. Disable GS output and
checkpoints, and return before the GS output section. This avoids creating
new restart/GS payloads, addressing the automatic-review rejection of the
initial separate-GS/restart job command. Existing production jobs are untouched.

The user subsequently authorized GS collection explicitly. The implemented
comparison therefore reuses normal GS restarts and an ordinary physical RT
run. The finite-frequency diagnostic is a separate invocation after RT
initialization and returns before physical propagation. This supersedes the
proposed embedded-GS propagation above.

Matched impulse window: approximately20 fs. Native finite-frequency poles are
not inferred from finite-time Fourier maxima. Nr16 requires a smaller Taylor
step; dt=.05 au failed conservation at step17, so dt=.0125 au is used there.
The user confirmed trans_longi='tr' for the known artifact.

First validate the numerical resolvent against a known oscillator and its
time-domain transform. Match actual-Si current derivatives to finite-A native
current. Check directional-difference step convergence and Krylov residuals.
Compare self-consistent and frozen-potential impulses at fixed Nk=2^3 and
Nr=8^3,12^3,16^3. Track pole position and weight, without subtracting a known
artifact or enforcing an insulating DC response. Density-coupled finite-q
chi0 must be distinguished from velocity/current contact artifacts.

The final optical solver uses a real rational Krylov basis, preconditioned by
static occupied-only Sternheimer inverses. Distributed real QR supplies the
weighted defect norm. Every frequency is rechecked with the original native
operator. Frequency solves and finite-time companions have separate outcomes.
The time reconstruction solves in the full right-eigenvector basis so that
degenerate modes do not depend on vendor-specific left/right pair scaling.

An impulse-direction current derivative is symmetry-null. Test that direction
and its imaginary orbital quadrature together, scaling their errors with the
nonzero derivative rather than dividing the null derivative by a tiny floor.
The relative threshold remains 1e-5; finite-value rejection covers both tests.
The native physical RT and finite-frequency diagnostics remain separate runs.

The user identified the intended low-frequency artifact with the residual
constant current discussed immediately after Fig.4 in the SALMON code paper
(Noda et al., CPC235,356,2019). The transverse comparison is therefore assessed
using both raw current and stationary density chi0, not solely finite-energy
local maxima. Self-consistent optical and independent density-head diagnostics
completed at Nr8/12/16, Nk2; an additional Nr12,Nk4 optical comparison passed.
Direct density q=0 body and finite-q channels completed at Nr8/12. An optional
Nr16 finite-q/body extension was stopped after the three optical-head tests
completed; its unfinished data are excluded. No production W was replaced.

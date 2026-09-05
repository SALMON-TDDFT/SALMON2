# Fixed-Density DG Interface Continuation Design

## Goal

Determine whether the poor Hybrid divided-SCF convergence is caused by applying the full SIPG interface operator from the first Schwarz epoch.  Reuse the converged ordinary-DC density and potential, keep them fixed, and introduce the complete SIPG interface operator continuously from zero to full strength.

This first implementation is a diagnostic non-self-consistent DG ground-state route.  It must not be presented as a converged Hybrid density-SCF result.

## Evidence and hypothesis

The simple- and Pulay-mixed Si64 runs show nearly the same behavior.  The Schwarz residual starts near `1.9e3` and remains near `2.0e3` through epoch 1000, although the common-occupation electron defect is normally at roundoff.  Density mixing reduces the density-change measure substantially but does not reduce the orbital residual.

The working hypothesis is therefore that the fixed operator seen by the Schwarz iteration, especially the full-strength DG interface term, dominates the instability.  This is not yet established: failure at zero interface strength would instead point to the Schwarz solver, overlap metric, or preconditioner.

## Operator continuation

Use

\[
H(\lambda)=H_{\mathrm{volume}}[\rho_{\mathrm{DC}}]+\lambda H_{\mathrm{SIPG}},
\qquad 0\leq\lambda\leq1.
\]

`H_SIPG` is the complete interface contribution, including its diagonal and off-diagonal jump, consistency, and penalty contributions.  Scaling the complete block by one real scalar preserves the assembled Hermitian pairing.  The volume kinetic term, local and nonlocal potential terms, and overlap matrix are not scaled.

The local Schwarz preconditioner must be assembled from the same `H(lambda)` used by the operator callback.  Using a full-strength interface diagonal in the preconditioner while applying a partial interface operator would make the diagnostic ambiguous.

## Schedule and controls

The first route uses a fixed schedule.  The configured Hybrid mixing rate is interpreted by this diagnostic route as the DG interface increment:

\[
\lambda_{k+1}=\min(1,\lambda_k+\alpha),
\]

where `alpha` is the configured rate and the default is `0.2`.  The schedule includes `lambda=0` and reaches `lambda=1` exactly even when `alpha` does not divide one.  A nonpositive or greater-than-one rate is rejected collectively.

At each continuation point, run at most `dg_hybrid_fragment_cg_steps` Schwarz CG steps, normally three, and carry the accepted coefficients into the next point.  The basis inventory, rank-to-fragment mapping, occupations, and fixed DC potential remain unchanged.  Rollback retains the last accepted coefficients and does not advance lambda.

Density Pulay/Broyden history is not used or updated in this mode.  A future adaptive-lambda route or short post-continuation density-SCF route is outside this first diagnostic implementation.

## Data flow

1. Restore or compute the ordinary-DC seed with the existing exact rank/fragment compatibility rules.
2. Build the projected WF+PW basis, immutable volume operator, SIPG interface rows, overlap rows, and Schwarz schedule once.
3. Extract the ordinary-DC core density and construct the corresponding total local potential once.
4. Initialize the Schwarz state at `lambda=0`.
5. For every continuation point, assemble the matching local preconditioner, perform the capped Schwarz update, assign the common 300 K occupations, and record diagnostics.
6. At `lambda=1`, assemble the full operator and perform exactly one terminal generalized LCFO diagonalization.
7. Publish the terminal state without a post-LCFO density update.

## Diagnostics and acceptance

Each continuation point reports lambda, accepted CG steps, residual, orthogonality defect, common-occupation electron defect, energy or Rayleigh trace, and the norm of the scaled interface action.  It also reports whether the step was accepted or rolled back.

The diagnostic interpretation is:

- If the residual is already large or fails to decrease at `lambda=0`, the DG interface term is not the primary cause.
- If the residual is controlled at `lambda=0` but rises with lambda, the first lambda interval showing the rise localizes the interface-coupling problem.
- Reaching `lambda=1` is not by itself proof of density self-consistency.  The result remains a fixed-density DG/LCFO state.

Collective disagreement in lambda, schedule position, basis generation, operator fingerprints, or rank-to-fragment mapping is fatal.  Nonfinite diagnostics or rejected Schwarz updates trigger collective rollback and stop the diagnostic with the last accepted lambda reported.

## Testing

Unit and MPI tests must establish:

- the fixed schedule includes zero and reaches one exactly;
- all ranks use the same lambda;
- the Hamiltonian action and local preconditioner scale the same complete SIPG block;
- `lambda=0` removes every SIPG interface contribution and `lambda=1` reproduces the existing full operator;
- the overlap and volume/nonlocal terms are invariant with lambda;
- DC density and potential callbacks are not updated during continuation;
- each point respects the configured CG-step cap and carries accepted coefficients forward;
- rollback does not advance lambda;
- exactly one generalized LCFO solve occurs after reaching one;
- no post-LCFO density update is introduced.

The Si64 diagnostic runs compare residual-versus-lambda and energy-versus-lambda using the same reused DC seed, MPI size, rank-fragment map, basis inventory, and CG cap.

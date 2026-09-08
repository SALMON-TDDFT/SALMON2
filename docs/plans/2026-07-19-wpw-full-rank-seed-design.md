# WPW Full-Rank Deterministic Seed Design

## Problem

The production Si64 calculation retains 160 states.  The current initial
vectors use terms of the form `sin(A_i+B_j)` and `cos(C_i+D_j)`.  The angle
addition identities restrict the W block to rank at most four and the P block
to rank at most four.  Consequently the global initial metric has effective
rank eight and `dg_metric_factor` returns `100 + (160 - 8) = 252` before the
first matrix-free eigensolver update.

## Design

Generate each complex coefficient from a deterministic stateless mapping of
the basis namespace (W or P), stable global basis ID, and retained-state
index.  The mapping includes an ID-by-state cross term, so it cannot factor
into a fixed small number of row and column functions.  Use bounded integer
modular arithmetic and map two independently salted residues to centered real
and imaginary components.  No random state, global gather, or communication
is introduced.

Place the generator with the matrix-free SCF implementation and expose one
initializer accepting owned W and P stable IDs.  `lcfo_flux` supplies its
existing ownership IDs and only fragment roots populate coefficients, exactly
as in the current ownership contract.

## Failure Handling

The initializer rejects shape mismatches, nonpositive IDs, and nonpositive
state counts.  The existing collective algebra preflight and rank-revealing
metric factor remain authoritative; this change does not weaken the metric
cutoff or convert a rank failure into success.

## Verification

First add a regression using a Si64-sized distributed layout and 160 retained
states.  It must demonstrate that the old formula has effective rank eight
and that the production initializer has all 160 metric eigenvalues above the
production relative cutoff.  Then run the matrix-free SCF MPI tests, bounded
operator tests, source contracts, full MPI/EigenExa build, and a fresh short
Si64 preflight.  A production checkpoint run is permitted only after the
preflight passes the first algebra step without rank collapse.

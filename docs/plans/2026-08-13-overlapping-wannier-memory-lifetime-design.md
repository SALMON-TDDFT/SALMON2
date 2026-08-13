# Overlapping-Wannier Memory Lifetime Design

## Goal

Reduce the Si64 overlapping-Wannier per-rank peak without changing the MPI
decomposition, physical basis, Wannier gauge, gradients, or checkpoint data.

## Evidence

For the Si64 fixture, each rank owns a `28^3` buffered grid, a `16^3` core,
and 384 retained states.  After the core values and gradients are extracted,
the old buffered values and gradients remain allocated until the end of the
character-sector pipeline even though they are not read.  They account for
about 514.5 MiB per rank, or about 4.0 GiB over eight ranks.  Smaller duplicate
lifetimes include the 24 MiB `global_seed_values`/`w90_anchors` copy and an
8 MiB occupied LCFO source retained after its last use.

## Design

### Stage 1: shorten lifetimes without changing arithmetic

Immediately after copying the initial buffer values and gradients into their
core arrays, deallocate the old buffer arrays.  The existing post-character
core-to-buffer materialization remains the sole producer of the final buffer.
Move, rather than copy, `global_seed_values` into `w90_anchors` after its last
seed-space use.  Release assembly-only Wannier inputs immediately after M/A
assembly, and release the occupied LCFO source after its last density check.

This stage changes allocation lifetime only.  Numerical operations and their
order remain unchanged.

### Stage 2: stream the Gamma transform

Replace full-size `new_values` and `new_gradients` temporaries in
`apply_dg_w90_gamma_transform` with one-point value and gradient work vectors.
For each spatial point, compute all transformed states from the unchanged old
point vector, then assign the completed point back.  Preserve the existing
canonical ordering/sign selection and all MPI receipts.

Stage 2 is kept separate because it changes the implementation of a numerical
transform, even though the intended result is identical.

## Safety and verification

A static lifetime checker fixes the ordering requirements and rejects renewed
full-buffer retention or duplicate seed/anchor allocation.  Existing focused
MPI fixtures compare values, gradients, fingerprints, and rank independence on
1/2/4/8 ranks.  The monitored Si64 MPI8 run supplies process RSS evidence.  No
build or additional MPI test is run concurrently with the active Si64 job.


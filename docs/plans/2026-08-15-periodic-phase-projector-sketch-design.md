# Periodic-Phase Projector Sketch Design

## Problem

`align_dg_w90_character_sectors_by_periodic_phase` currently fingerprints the
aligned sector by reconstructing every row of the full global projector.  For
each of `N` global rows it broadcasts an `m`-component sector row, allocates a
global `N`-component projector row, and performs an `MPI_Allreduce` of that
row.  This costs O(N^2) communication and dominated the Si64 post-Wannier run.

## Design

Replace the full-projector stream with a small fixed set of deterministic
complex probes.  For every probe `v`, compute

```
u = A^H v
w = A u = A A^H v
```

where `A` is the metric-weighted aligned sector frame.  `u` requires one
length-`m` `MPI_Allreduce`; `w` is evaluated locally and streamed in global-row
order into the fingerprint.  The probe phases depend only on the global row ID,
so the result is invariant under MPI decomposition and local row ordering.
Because `A A^H` is unchanged by `A -> A U`, the receipt remains invariant under
internal target-frame rotations.

Use the same bounded deterministic probe family and tolerance quantization as
the existing joint periodic-center canonicalizer.  Guard every value before
integer conversion and reject collectively on nonfinite or out-of-range sketch
values.  Remove the `projector_row(global_row_count)` allocation and update the
workspace receipt accordingly.

## Verification

- Add a source-route regression that rejects the old global projector-row
  allocation and global-length projector `MPI_Allreduce` in this routine.
- Keep and extend the W90 MPI fixture to require the same alignment fingerprint
  for target-frame rotations and MPI 1/2/4/8 decompositions.
- Run the focused W90 MPI test, route checks, `git diff --check`, and a serial
  Release build (`-j 1`).


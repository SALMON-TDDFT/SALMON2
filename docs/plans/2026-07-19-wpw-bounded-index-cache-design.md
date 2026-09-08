# WPW Bounded Operator Index Cache Design

## Problem

The full-basis Si64 preflight passes deterministic-seed rank validation but
spends nearly all sampled CPU time in `apply_blocks`.  Every H/S application
repeats `find_id` linear searches for both endpoints of every WW, WP, and PP
sparse entry.  These mappings are invariant for the lifetime of a bounded
operator epoch, so recomputing them inside every LOBPCG batch adds avoidable
work to the dominant production loop.

## Design

Keep all existing global stable-ID arrays as the authoritative sparse format
and checkpoint/provenance representation.  Add six private-to-the-operator
integer position arrays corresponding to WW row/column, WP W/P, and PP
row/column endpoints.  Resolve and validate them once while constructing the
candidate bounded operator.  Publish the candidate only if every endpoint is
valid collectively.

`apply_blocks` uses the cached positions directly and retains the existing
sparse-entry order, vector update order, MPI owner schedules, and floating
point summation order.  Rebuilding an operator creates new caches; releasing
it deallocates them through default derived-type assignment.

## Verification

Add a source contract that rejects `find_id` from the apply loop and requires
all six caches.  Keep the two-rank dense-oracle H/S fixture as the functional
equivalence test, including rebuild and stale-epoch failure paths.  Run the
matrix-free wide-spectrum fixture, full build, and a fresh full-basis Si64
preflight.  Compare the time from `Hamiltonian matrix: done` to the first
algebra outcome against the preserved pre-cache preflight.

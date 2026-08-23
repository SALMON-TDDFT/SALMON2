# Nonlocal Projector Range and Symmetry Diagnostic Design

## Purpose

Determine whether broad Wannier tails produce material overlaps with
pseudopotential projectors on neighboring or remote fragments, and whether the
total-system nonlocal operator respects affine generator operation 5.  The
diagnostic must not replace, symmetrize, or otherwise alter the production
Hamiltonian.

## Evidence and scope

The direct total-system `hpsi` projection is Hermitian to `6.6e-15`, but its
nonlocal covariance residual is `0.7855` for operation 5.  The former fragment
projector assembly is also Hermitian but uses a different partition and
normalization convention, so its raw matrix magnitude is not a valid reference.
The next diagnostic therefore works directly from `dc%ppg_tot`, the same
projector tables used by production `hpsi`.

## Architecture

Add a standalone MPI diagnostic routine in a focused DC module.  It consumes
row-owned Wannier values, their physical grid IDs and centers, the total-system
projector grid, atomic geometry, and one affine operation.  A cached physical-ID
redistribution places bounded Wannier tiles on `dc%mg_tot`.  For every global
nonlocal channel, ranks accumulate
`<beta_a,lm|w_i>` from their local projector support and combine it collectively.

The routine reports only reductions, not the full overlap table:

- projector contribution summed by center-fragment, adjacent-fragment, and
  remote-fragment distance class;
- maximum and total remote fraction over Wannier functions;
- operation-5 atom-map residual;
- norm and maximum defect between overlaps of symmetry-related Wannier/atom
  projector pairs;
- count of channels whose symmetry partner is missing or ambiguous.

Distance classification uses periodic minimum-image distance between the
Wannier center and projector atom.  Fragment classification is derived from the
existing `2x2x2` physical fragment geometry, not MPI rank arithmetic.  The
operation-5 atom permutation is obtained from the atomic positions, species,
integer rotation, and fractional translation with a one-to-one periodic match.

## Memory and communication

Process Wannier functions in tiles of at most 16.  Store at most
`tile_width * Nprojector` complex overlaps plus bounded redistribution buffers.
Use one collective reduction per tile, never one collective per orbital or
projector.  Do not allocate a full-system orbital array and do not retain the
overlap table after its receipts are accumulated.

## Safety and interpretation

The routine is diagnostic-only and runs once after the final Wannier basis and
centers are available, before Hybrid-SCF.  Failure to construct an exact atom or
projector partner is reported and stops the diagnostic; it does not silently
average the Hamiltonian.  A large remote fraction indicates that fragment-local
projector assembly is invalid for the current basis.  A large operation-5 pair
defect with complete partners instead identifies a total-system projector or
symmetry-action inconsistency.

The temporary raw comparison against `assemble_ow_nonlocal_rows` is removed
from every Hamiltonian build once this diagnostic is connected.

## Verification

Unit MPI tests use synthetic short-range projectors and Wannier functions with
known local, adjacent, and remote contributions.  Tests cover exact operation
pairing, missing partners, rank-independent results on MPI 1/2/4/8, bounded
workspace, and tile-width equivalence.  The Si64 acceptance run uses MPI 8 and
`OMP_NUM_THREADS=1` without a time cutoff and records the new receipts before
the existing strict covariance gate.

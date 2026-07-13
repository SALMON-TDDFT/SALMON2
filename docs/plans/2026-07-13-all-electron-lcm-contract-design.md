# All-electron Local Chern Marker contract design

## Scope

Fix both scalar/collinear and spin-orbit Local Chern Marker paths so that the
reported marker is the all-electron marker.  Reject unsupported non-idempotent
occupations and fail collectively when the dual-state overlap matrices cannot
be inverted reliably.

## Occupied projector

The occupied projector is built from unweighted occupied orbitals.  Occupation
weights must not be applied before Lowdin orthonormalization because that
transformation cancels every positive diagonal weight.  Instead, each sharp
spin channel contributes with its electron multiplicity:

- non-SOI, `nspin=1`: multiplicity 2;
- non-SOI, `nspin=2`: multiplicity 1 for each spin channel;
- SOI spinors: multiplicity 1.

The selected occupations must equal the corresponding multiplicity within a
small numerical tolerance.  Fractional, smeared, and restricted odd-electron
occupations are outside the idempotent-projector contract and terminate with a
clear diagnostic.  Empty orbitals remain excluded by the existing cutoff.

## Ill-conditioned dual overlaps

The S1 and S2 reciprocal condition estimates must exceed `1e-10`.  The orbital
root reports the label, reciprocal condition estimate, and norm.  It then
broadcasts failure to every rank in the orbital communicator so all ranks stop
consistently before inverse coefficients are distributed or marker output is
written.  No regularized inverse is substituted silently.

## Verification

Source-contract tests cover the multiplicities, removal of pre-Lowdin
occupation scaling, fail-closed fractional occupations, collective condition
failure, and preservation of the zero-local-occupied-rank fix.  The final gate
is the MPI/EigenExa build plus the focused LCM regression checks.


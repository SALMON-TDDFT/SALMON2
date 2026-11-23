# Built-in material adapter implementation plan

Approved M4 continuation, first LD/SBE subset. DFT/TDDFT and full M5/M6 stay open.
Add explicit builtin_material_registry() and material_from_catalog() in a separate
module, reusing config.MaterialRegistry and actual LorentzDrude/SBEInitialState.
LD parameters require explicit atomic_gaussian units and existing pole values.
SBE parameters require atomic units, local relative directory/sysname and existing
SALMON text-loader dimensions/lattice/electron options. Record three file hashes;
rehash around loading and compare with resolved catalog before model publication.
No discovery, native initialization, compile, network or MPI during resolution.

Return ConfiguredMaterial with stable material_id and definition_digest, original
model accessible as model. Its new_state forwards existing arguments unchanged,
returns identity-carrying state which delegates numerical methods/context cleanup.
LD E and SBE A/c semantics stay distinct. No point/run ID, custom backend factory,
initial_state recipe selection, serialized restart or output propagation claimed.

1. Fail tests before module creation: LD/SBE resolve/build, input rejection,
   modified SBE file rejection, separate states and direct-model numerical equality.
2. Implement adapters and copied identity records without touching native kernels.
3. Executable LD recipe example; document SBE text recipe/local staging; JP/EN/AI
   and public inventory/packaging updates.
4. Read-only review, local tests and installed-wheel sample. No Actions.

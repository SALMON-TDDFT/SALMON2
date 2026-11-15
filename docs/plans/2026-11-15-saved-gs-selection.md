# Saved GS selection by initial_state_id

Approved M4c continuation. Add tddft_from_catalog(catalog, initial_state_id,
base_directory=..., executable=..., max_steps=..., threads=..., timeout=...).
The selected entry must reference a DFT material and source object with exactly
format=salmonpy.gs/1, relative path and manifest_sha256. No automatic prepare,
network, arbitrary SALMON directory, GS/2 conversion or generic backend selection.

Validation: anchored manifest; existing GroundStateResult integrity/convergence;
exact configured scientific signature including ordered fractional positions;
scientific sections of actual gs.inp and actual pseudo bytes; same executable hash.
Use existing GS->TDDFT loader, then verify copied snapshot hashes against original
GS inventory and recheck pinned manifest. Close candidate on any post-load failure.
Propagate initial_state_id and source manifest digest alongside material identity
into resulting model/state. Do not change direct legacy calls or promise MPI.

Implementation: failing fixture tests, resolver in salmonpy/initial_states.py,
identity fields on existing wrappers, paired guide + AI JSON + example, local
review/full tests/sdist->wheel. Final sample generates a short GS, explicitly pins
it, reopens by ID with sockets disabled, checks independent TDDFT and existing Si
reference. Use existing audited runtime, one thread/worker and max two workers.
No Actions, no full runtime build. Preserve source records; no GS payload collection.

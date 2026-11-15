# DFT material configuration and GS-to-TDDFT identity plan

Approved M4 continuation. Register salmonpy.dft/dft alongside existing LD/SBE.
Parameters: explicit local structure path/format/frame, pseudopotential path/Z/
local channel list, XC, states, kmesh and grid counts OR spacing value/unit.
Use existing StructureInput, grid_from_spacing, Pseudopotential and dft_from_atoms.
Record normalized values, grid report, neutral electrons and hashes; revalidate
files at model creation. Optional ASE dependency only when this backend resolves.
Do not infer XC/potentials or download, initialize MPI, compile, or execute at resolve.

ConfiguredMaterial.run delegates DFT.run and returns ConfiguredGroundState with
material ID/digest; result.tddft_initial_state wraps the existing verified GS path
and returns identity-carrying TDDFT model. Model contexts close only owned models
with a close method. The wrapper does not persist identity to disk or implement
arbitrary saved-GS selection. Existing native/runtimes and A/c->j semantics remain.

1. Tests first: explicit count/spacing input equivalence, hash changes and invalid
   parameters, no runtime in resolution, ID forwarding across mocked existing GS
   boundary. Do not call mocks a native integration test.
2. Add DFT resolver and wrappers; schema/docs/AI/sample update.
3. Read-only review; source/wheel tests; installed offline preparation. Reuse an
   audited existing local runtime for short Si GS/TDDFT and independent workers,
   with one thread each and synthetic/smoke scope explicit. No new CI or compilation
   of full runtime; no broader XC convergence claims.

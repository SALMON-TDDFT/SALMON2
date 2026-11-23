# Structure and grid preparation plan

Approved continuation of material settings M2/M3 and offline preparation.
Add optional ASE-based StructureInput (detached Atoms, explicit file format/frame,
file SHA256 and metadata) and immutable GridSelection. Structure files are local
CIF or extxyz, parsed with ASE parallel=False and no filename @ indexing; no
implicit cell standardization, repetition or coordinate wrapping. CIF asymmetric
units necessarily expand per ASE's crystallographic reader; document this.
Grid helper accepts explicit angstrom/bohr width, uses ceiling plus minimum and
explicit integer alignment, reports actual width. Initially axis-aligned positive
orthogonal periodic cells only; no rounding based on MPI or convergence claims.
DFT and TDDFT ASE entry points accept either old grid counts or new spacing/unit,
retaining existing electron/XC/PP and runtime guards. Add preparation example
binding a staged structure file to material ID via config's explicit registry.

1. Failing tests for file/frame/hash handling, independent Atoms, grid boundary/
   alignment/units and rejection cases; old/new DFT input_text equality.
2. Implement salmonpy.structure and width support in salmonpy.ase.
3. Local example through preparation bundle, docs JP/EN and AI/public API updates.
4. Review; focused/full tests, source distribution→wheel, installed offline example.
No Actions, no native numerical changes, no automatic fetching or production run.

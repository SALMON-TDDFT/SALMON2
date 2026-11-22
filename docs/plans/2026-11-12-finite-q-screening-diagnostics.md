# Finite-q screening diagnostics

Approved scope: record every computed finite-q dielectric head and compare
screened and unscreened Si response at matching k mesh and time window.

The driver already reduces the finite-q inverse heads across MPI ranks. Print
those existing values on the communicator root after successful assembly. Record
the q index, fractional coordinates, minimum Cartesian |q+G| (bohr^-1), and the
local-field-inclusive inverse head. Retain separate optical-limit diagnostics.
No dielectric solve, interpolation, W normalization, or RT Hamiltonian changes.

1. Extend the existing pipeline test to require one record per finite q,
   root-only output, correct reciprocal-space radii, and finite inverse heads.
   Observe failure with the unchanged driver.
2. Add the root-only Fortran output and document its units and definition.
3. Run the source prohibition audit, compile out of source, audit application
   objects, and run the reference comparison with 1 and 4 MPI, 12 OMP.
4. Inspect existing Fugaku results and unchanged validated executable hashes.
   Reuse the converged Nr12/Nk16 GS on Fugaku, collect no GS payload, and obtain
   missing finite-q diagnostics and W-free kicked/zero-field RT as needed.
   Keep existing fixed jobs and executables unchanged.
5. Compare currents with their zero-field controls, transform using the same
   SALMON window and trajectory duration, and save spectra and provenance under
   the canonical calc2026 result directory. Do not compare peak heights across
   different windows or label TD-TBmBJ as IPA/RPA.

Defer box/cutoff/empirical-W scans until the finite-q and W-free comparisons
identify what requires further convergence. At Nr12, boxes8 is invalid.

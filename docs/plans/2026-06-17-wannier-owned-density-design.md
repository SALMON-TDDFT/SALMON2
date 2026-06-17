# Wannier-Owned DG Density Design

DG-RT should assign localized orbitals to fragments by their Wannier centers, not by cutting wavefunctions at fragment core boxes.  A fragment owns the Wannier functions whose centers lie in its core region.  The real-space values of those Wannier functions are still allowed to extend into neighboring fragment buffers, and density contributions are summed to the global grid owner through the existing sparse density exchange.

The first implementation keeps the current fragment-basis density path as a diagnostic fallback and adds a Wannier-owned density path.  DC export writes local Wannier metadata including centers.  RT reads the metadata, reconstructs Wannier real-space amplitudes from `local_wannier_coef` and buffered fragment basis functions, accumulates occupied density on all grid points covered by those functions, and uses the existing density communication maps to add contributions to the authoritative grid.

Initial states remain Flux-inclusive fragment-local block Hamiltonian eigenstates.  External impulse excitation is applied as a first-step length-gauge rectangular electric field, not as a direct coefficient kick.

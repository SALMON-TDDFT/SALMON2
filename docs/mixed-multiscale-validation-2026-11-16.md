# Mixed endpoint LD/SBE regression — 2026-11-16

Apple Silicon, OpenMPI, mpi4py, OMP1/BLAS1. Source-script and installed-wheel runs use 3 ranks (one field, LD point10, SBE point20), dt=.02,800 steps, transverse z polarization. LD and SBE parameters are synthetic. All A and j values match the independently executed serial field -> raw LD/SBE schedule exactly. Maximum absolute currents are 1.9125200721218774e-5 and 2.7166724704943676e-5 a.u.; both responses are nonzero. Communications are measured separately; a timing failure does not redefine calculation success.

Source archive -> wheel -> installation -> MPI module execution passed. Outputs (point/material mapping, A/j histories, per-rank communications) are in the persistent local evidence directory salmonpy-state-inputs-20261116/mixed-final-installed. Result SHA256: `454a63fb7915ad0ac5ee7e5a31eab7d6ccb3d3de1d0299c4de3248c004189d41`.

An initial invalid x-polarized trial produced zero response. The sample was corrected to transverse z and now requires nonzero current at every point, preventing vacuous parity. A separate follow-up is needed for Multiscale1D validation of longitudinal components in prefilled incident history when boundary samples are all zero; this sample does not fix that validator.

The optional managed TDDFT third point is prepared in the example but not executed in this validation. Internally parallel TDDFT, native MPI E/H with field-local LD, restart of the whole mixed session and physical convergence are not established by this result. The demonstration uses one rank per LD point only for the small endpoint test, not the intended production LD decomposition. No GitHub Actions or push.

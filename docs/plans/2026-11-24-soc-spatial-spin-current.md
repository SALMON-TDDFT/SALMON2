# CPU spatial SOC spin current

The model uses the existing Pauli-current definition. For each projector, sum the ordinary overlap over its spatial group. Keep the coordinate-weighted Pauli overlap local, contract with the summed ordinary overlap and the existing projector weight, then use the existing RKO sum exactly once. This avoids replicating the nonlocal contribution once per spatial rank. There is one new spatial-group collective per orbital/k point, only when spin observables are evaluated; no application-wide vote is added.

The charge-current calculation already updates wavefunction halos before RT spin gradients. The spin observer reuses those valid halos. Fixed-spinor tests supply explicit valid halos; x/y Ni spatial tests exercise the real charge-current/observer chain.

Serial and spatial CPU calls share the two contraction helpers. No mutable module cache or state resource is introduced. CPU capability is declared by the local C query `sp_tddft_soc_spatial_protocol()=1`; OpenACC returns0. The Python held-GS, initial-bundle and registry routes require that declaration for spatial SOC.

## Numerical finding

An independent Pauli-matrix oracle did not match the old GNU15.2.0/Mac O2 CPU spin-current calculation. A copy of the old formula compiled at O0 agreed with the oracle; O2 differed. Diagnostic executables/logs are retained. The compiler-internal cause is not claimed as established. Splitting projection and coordinate contraction avoids the observed sensitivity without changing the physical formula or globally disabling optimization.

Old-library negative tests are reproducible with `test_soc_spin_current_spatial` compiled using `-DSOC_SERIAL_ONLY`. New tests compare all twelve nonlocal components against the matrix oracle, and local microscopic currents and local+nonlocal band currents against spatially undistributed results. They cover MPI1/2/4, OMP1/2, an empty local projector support and independent spatial replicas.

Actual Ni software tests use4 SCF iterations, not a converged physical material. Primitive input:18 electrons,36 states,300 K,2x2x2 k points. Cubic fcc input:4 atoms,72 electrons,144 states with the same temperature/mesh and16-cubed grid. Electronic datasets are never exported. Charge-current and magnetization tables remain exact relative to the prior fixture; spin-current values are corrected, not forced to match the old result. Pauli0 agrees with native matter current Jm.

GPU, Fujitsu, changed-topology SOC disk restart, real-material convergence and fast-HSE SOC are outside this result. Normal/performance comparisons distinguish checked numerical recovery from explicitly protected performance rollback. No backward-data conversion or reader-side spin-current correction is introduced.

The corrected pseudo_pt_current_so file and test also pass an O3 compile with MPI1/2/4 and OMP1/2. Remaining solver objects in that targeted experiment retain O2; this is not a whole-runtime O3 validation. Spatial2×orbital2 recovery also passes on Ni, including Pauli0/Jm consistency.

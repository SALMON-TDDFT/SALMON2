# Native HSE spinor orbital stability regression

`check.py LOG` checks the three scalar stages emitted by
`SALMON_HSE_STATIONARITY=solve`: initial orbitals, production `solve_orbitals`,
and final Ritz rotation. The diagnostic uses a frozen physical Hamiltonian,
fixed fractional occupations and direct exchange (`hse_taylor4_full`). It
returns before SCF and exports no wavefunctions.

Use a local, existing Ni restart; do not collect GS data for this test. Copy
only its input and symmetry file to a fresh temporary directory, reset the
restart step, read `wfn`, use `nscf_init_redistribution=0`, and set `ncg=32`.
Run the audited executable with two k ranks, then check its local stdout log.
The existing initial state must be nonstationary (>1e-6 Ha); both output
stages must have occupied-density commutator residual <1e-6 Ha, Gram error
and projected-Hermiticity defect <1e-10. All metrics must be finite.

Recorded RED on the pre-fix solver: residual 2.89384e-3 -> 3.67172 ->
2.40239 Ha. GREEN on the production blocked solver: 2.89384e-3 ->
9.54035e-8 -> 9.54035e-8 Ha. The independent diagnostic `blocked` mode
provided the same result before changing the production solver.

This is a targeted material regression, requiring the transient local
restart. It does not establish SCF convergence, RT conservation, general
basis convergence, or portability to another compiler. It is intentionally
separate from the standalone synthetic `unit_exx_stationarity` metric tests.

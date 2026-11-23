# Owned k-point, pseudopotential and XC inputs

User direction: separate shared k-point/pseudopotential/XC input. Preserve native update order and fast-HSE results; never implement active-global swapping.

1. Add value-only configs and deep-copy them into the relevant owner.
2. Extract explicit k generation including symmetry; retain compatibility wrapper and route HSE mesh/symmetry consumers.
3. Extract the complete radial PP loader including format readers and parser scratch flags.
4. Retain functional identity and hybrid/rVV10 settings; pass pre-SCF evaluation phase explicitly to calc_xc. Route native exchange selection/scaling consistently.
5. Run independent lower-level probes, interleaved full-context comparisons, original fast-HSE regression, Python and sdist-to-wheel verification locally.

Scope boundary: PP grid construction and atom mapping, global SCF phase and algorithm controls, communicator/file-unit ownership and initialization of unlike full materials remain. Full native TDDFT keeps its common grid/k/XC guard. New helpers are internal Fortran interfaces, not public Python/C APIs, and retain native fatal-I/O/error behavior. Initialized Libxc/HSE handles must never be copied.

# Explicit real-space grid initialization

Continue the approved shared-input extraction. init_dft_system still derives cell,
k vectors and atoms from global input; pseudo/HSE also retain dependencies. Do not
relax the whole-context common-grid/XC guard yet.

Extract init_grid_whole into a compatibility wrapper that captures input and a
configuration-only implementation. The value configuration holds periodicity,
requested point counts, spacing and consistency-check lengths (rsize normally,
al_em for Maxwell). Preserve grid rounding, bounds, halo and index ordering.
Explicit entry validates finite/positive geometry, spacing/point consistency and
integer bounds before allocation, returning failure without changing the grid.
Legacy wrapper retains native fatal-error behavior for invalid input.

Test periodic and isolated, point/spacing modes, legacy vs explicit, two shapes,
global-input poisoning and invalid-call retention. Fresh audited full build,
whole TDDFT regression and native-driver reference, Python and wheel checks.
Document that this is lower-level geometry extraction, not heterogeneous TDDFT.
No Actions.

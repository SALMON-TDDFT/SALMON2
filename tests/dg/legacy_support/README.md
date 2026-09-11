# Test-only legacy references

These modules are excluded from the SALMON production build. They support
standalone numerical regression tests, not selectable GS or RT routes.

- `dg_hybrid_scf`: old global fixed-basis SCF controller.
- `dg_hybrid_divided_scf`: old divided SCF callback driver.
- `dg_hybrid_continuation_state` and `dg_hybrid_continuation_controller`:
  old continuation state/scheduling tests. The controller imports the live
  publication policy from `src/gs/dc/dg_hybrid_publication_policy.f90`.
- `dg_hybrid_real_space_residual`: standalone residual diagnostic reference.
- `rt_dg_overlapping_wannier`: old V3 coefficient-RT numerical reference.

Test runners compile these explicitly. Do not add them to production CMake
source lists. Current production uses fragment-local Hybrid GS, terminal
LCFO refinement and the v5/Exp RT path. Old selectors fail explicitly rather
than converting old checkpoints or selecting a different calculation.

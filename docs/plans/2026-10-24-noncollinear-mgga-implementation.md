# Corrected noncollinear MGGA implementation

Approved objective: implement the corrected 2023 noncollinear BR-CS functional
and compare spin trajectories with LSDA; independently generate and validate
Ni SOC+NLCC KY potentials, using larger radii to reduce internal oscillations.

## Repository

Implemented in the existing clean KY/SOC-enabled checkout
`.worktrees/codex/cartesian-soc-fix`, branch codex/fix-cartesian-soc-symmetry,
not the older SALMON2 checkout originally proposed. No historical source
snapshot, fixed executable or production job was modified.

## Completed implementation

1. `nc_mgga23_kernel.f90`: energy and 36 forward-AD partial derivatives,
   implicit BR root derivative, corrected correlation Eq. (6) checked against
   the primary 2024 APS erratum PDF. gamma=.8, tau without factor1/2.
2. `nc_mgga23.f90`: spinor density/tau/current accumulation, Cartesian
   derivative and MPI halo handling, full variational matrix gKS action.
   Keeps all off-diagonal spin components and uses explicit system payload.
3. XC selector `nc_mgga23`, input compatibility validation, formulation log,
   CMake registration and no LSDA fallback. No new namelist variables.
4. Hamiltonian complex path and isolated band-energy correction. Operator
   coefficients midpoint-averaged together. LSDA SOC matrix also averaged
   when predictor/corrector is enabled, for a comparable integration order.
5. Post-impulse operator refresh; initial/restart magnetization sample;
   native spin-file columns retained, additional local-norm/commutator file.
6. JSON/CSV/optional PNG comparison tool with unit/time validation and
   explicit first-sample reference time. Rotation and compensated initial
   magnetization handled without a spurious normalized loss.

## Validation

See canonical calc2026/nc-mgga23-validation-20261024 for build/test records
and synthetic comparison output. Tests exercise the kernel, orbital energy
variation, Hermiticity, spin rotation, full hpsi and Taylor integration,
restart-like coefficient reconstruction, actual writer and XC dispatcher.
MPI tests cover both orbital and x-domain decompositions; domain test uses
all four finite-difference radii and real halo exchanges against full grid.
No physical GS data or material demagnetization benchmark was collected.

## Explicit limitations and follow-up

- Nonzero NLCC is rejected. The available frozen core density alone does not
  define a consistent kinetic-density/current extension. Therefore the new
  Ni SOC+NLCC file is NOT a valid input for this MGGA selector yet. It is a
  separately validated LSDA candidate. Do not strip its core data to bypass
  this guard. A physically justified core extension needs its own validation.
- OEP/KLI, accelerators, DC, moving ions, forces, frozen-XC and non-midpoint
  propagators are not supported for this selector. Existing other selectors
  retain their routes.
- Extra torque columns describe the combined LOCAL multiplicative XC
  commutator only, not separate x/c torques or the complete gKS spin balance.
- Synthetic no-SOC tests establish numerical conservation, not the accuracy
  of real material SOC dynamics. Full checkpoint-file round-trip, physical
  material accuracy and Fujitsu compilation remain unverified.
- Local input documentation is supplied; an external SALMON-DOCS checkout
  was not modified. No merge, push or remote production submission performed.

The Ni radius exploration, rejected ghost-bearing candidates, independent
radial screening and remaining transferability limits are documented in
utility/crps_soc_nlcc/NI_VALIDATION.md and its canonical calc2026 records.

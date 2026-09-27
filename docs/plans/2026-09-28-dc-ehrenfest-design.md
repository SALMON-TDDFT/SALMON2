# DC Ehrenfest dynamics: corrected primary scope

Status: user correction on 2026-09-28 supersedes the finite-temperature BOMD
priority for laser-excited dynamics. This records the direction; it is not a
claim that a moving-basis force derivation or implementation is complete.

## Electronic state and initialization

Prepare the initial DC ground state by SCF; finite electronic temperature may
be used to prepare initial occupations. During laser excitation, propagate the
non-equilibrium electronic state with real-time TDDFT and ions with Ehrenfest
dynamics. Do not fit an instantaneous electronic temperature, reconverge a
ground state each step, or redistribute populations with a Fermi solve each
step. The new dc_thermal module remains useful for initial SCF and a separate
equilibrium BOMD route; E-TS is not the driven trajectory's energy functional.

For closed, unitary propagation the initial density-matrix eigenvalues are
preserved, while populations in an instantaneous orbital basis can change.
MLWF rotations must transform the entire density matrix consistently; treating
fractional unequal occupations as an unchanged diagonal matrix under arbitrary
rotations would change the state.

## Existing route and gaps

The conventional RT route already has yn_md integration in
src/rt/time_evolution_step.f90 and src/rt/md.f90. The PBEh LCFO route currently
supports fixed nuclei and impulse response, with explicit yn_md rejection in
src/io/inputoutput.f90. Reuse the existing ionic integrator only after a native
DC/LCFO force adapter and moving-state update have been validated. Removing an
input guard alone is not implementation.

1. Audit the RT density/state representation, LCFO overlaps, global atom map,
   force ownership, and the existing split ionic update. Specify an internally
   consistent electronic/nuclear discretization before changing propagation.
2. Derive instantaneous-state forces from the same discrete dynamics/energy.
   Separate explicit projector derivatives from moving LCFO basis, overlap and
   core-partition terms. A finite difference that reconverges SCF is a BO force
   reference, not a validation of a non-equilibrium Ehrenfest force. The missing
   BO adjoint correction from prior audits must not be transplanted unchanged.
3. Couple moving fragment images, projector updates, LCFO overlap/basis transport
   and total-ion forces. Treat basis connection terms and gauge transport
   consistently. Preserve MLWF/ACE acceleration, update exchange for the actual
   time-dependent density matrix and invalidate geometry/state-dependent caches.
   Evaluate rVV10 from the instantaneous assembled density in the existing
   adiabatic functional approximation. Begin with full MLWF support.
4. Connect the existing split ionic integrator with real-time propagation, then
   extend beyond the fixed-nuclei impulse-only excitation. Reject unsupported
   configurations until their force, state-transfer and field coupling work.

## Verification contract

- One-fragment/full-buffer limits against a supported conventional reference;
  frozen-state nuclear variations for explicit force terms, plus independent
  checks of moving basis/overlap contributions.
- Field-free conservation of electronic plus ionic energy, charge, orbital
  orthogonality and density-matrix spectrum; refine electronic and ionic steps.
- During a pulse compare total-energy changes with external-field work using
  the chosen gauge and boundary conditions. Do not demand energy conservation
  while the laser is doing work, or minimize E-TS during excitation.
- Translation, fragment-face crossing, buffer convergence, MPI decomposition
  and exchange/dispersion consistency. Static force diagnostics alone do not
  certify moving partitions.
- Initial target is mean-field Ehrenfest dynamics. Electronic thermalization,
  decoherence, surface hopping and a two-temperature model are separate models,
  not implicit consequences of this implementation.

Native DC-MD and PBEh LCFO moving-ion guards remain until these gates pass.

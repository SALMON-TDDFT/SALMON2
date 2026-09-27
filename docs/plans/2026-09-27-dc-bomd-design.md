# DC-BOMD design proposal

Status: approved by user continuation. Force certification precedes MD guard removal.
Branch: pbeh40-rvv10-water-md at db8841d2.

## Initial scope

Fixed-cell, orthorhombic, unpolarized PBEh40 / PBEh40+rVV10 Born–Oppenheimer
NVE dynamics. Each ionic step reconverges DC-SCF. Preserve MLWF gauge transport,
ACE and total-density rVV10 with the existing backend selection. Begin with
exx_mlwf_radius=0; finite source support needs its own force certification.
This does not reuse fixed-nuclei LCFO time propagation as an ionic integrator.
No restart, NPT/stress, or electronic nonadiabatic dynamics in this milestone.

## Audit of current code

- main_dft_md never constructs/passes s_dcdft; initialization_dft_md calls the
  conventional force routine. main_dft explicitly skips forces for yn_dc=y.
- calc_total_energy_dcdft accumulates core-restricted kinetic, nonlocal and
  exchange energies, then adds total-grid electrostatics/XC. Fragment SCF
  convergence alone must not be treated as proof that the assembled energy is
  stationary with respect to all fragment orbitals.
- init_fragment builds atom lists once and drops original global atom IDs.
  s_dcdft has grid maps but no fragment-atom-to-global-atom map.
- ne2mu_core updates weighted occupations; its final electron residual is not
  surfaced as a fatal failure. Finite electronic temperature requires an
  explicit energy/free-energy/force convention before admitting fractional
  occupations to an energy-conservation test.

## Force and propagation contract

1. Establish central-difference forces from converged static DC energies at
   two displacement sizes. This is a validation oracle only, not a scalable
   production force algorithm. Separate SCF, grid, buffer and displacement
   errors; compare one fragment, full-cell buffers and genuinely truncated
   buffers with conventional analytic forces where equivalent.
2. Implement a native DC force adapter with global atom IDs and periodic-image
   shifts. Compute ion-ion and local terms on the total grid; assemble nonlocal
   projector terms with the same partition/occupation conventions as energy.
   Derive/check fragment-response contributions instead of assuming that sums
   of ordinary fragment Hellmann–Feynman forces differentiate DC energy.
   Require finite-difference agreement before connecting this adapter to MD.
3. Select and document the stationary thermodynamic potential for occupation
   updates. Check electron count and chemical-potential convergence collectively.
   Use integer-occupation insulating fixtures for the first force comparison;
   fractional occupations require consistent entropy/response treatment.
4. Store ionic positions, velocities and masses once on system_tot. Velocity
   Verlet updates total ions; synchronize fragment images and rebuild affected
   projectors/atom lists at each step. Fixed real-space cores remain in place.
   Preserve orbital/MLWF data where valid and invalidate atom-dependent caches.
   Boundary crossings must neither lose atoms nor double-count forces.
5. A failed SCF, localization, electron-count or force check rejects the ionic
   step before the second Verlet update. Native yn_dc=y/theory=dft_md guards
   remain until force and moving-fragment tests pass.

## Validation gates

- Native projector and partition force against DC energy finite differences,
  with/without rVV10; conventional and one-fragment limits.
- Buffer convergence with actual truncated fragments, MPI decomposition parity,
  FFT backend parity, translated geometry and fragment-boundary crossings.
- Short NVE runs at dt and dt/2 from identical initial conditions, with tighter
  SCF comparison. Record drift, charge, force mismatch and per-step convergence.
- Water/solution scale tests follow the small-fixture certification; no claim
  of production scaling from a reference finite-difference driver.

## Static audit evidence

The existing H4 fixture converges for all eight displaced calculations (2 ranks,
SCF density threshold1e-10). For buffer4 (full-cell fragment support) and buffer2
(truncated fragments), central differences at0.002/0.001 bohr agree within
5e-6 eV/bohr. See ../results/pbeh40-rvv10/dc-md-force-audit.json. This establishes
usable energy derivatives for this fixture, not analytic force correctness,
MD support, or water convergence.

# Real-space PBEh Ehrenfest dynamics from DC initial states

Status: user correction on 2026-09-28: SALMON evolves wavefunctions on the real-space mesh. LCFO is used only to reconstruct the initial DC state on that mesh, not as a time-propagation basis. Supersedes the prior moving-LCFO-basis direction.

## Representation and dynamics

Prepare a DC ground state, reconstruct occupied wavefunctions onto the total real-space mesh using the existing conventional-from-DC reader, then release them from the LCFO subspace. Use SALMON's real-space Hamiltonian, Taylor4 predictor/corrector and existing ionic Verlet steps. MLWF/ACE accelerate exchange evaluated from evolving mesh orbitals. Do not reoccupy instantaneous states thermally or reconverge SCF during excitation. Initial SCF temperature is distinct from electronic dynamics.

## First native implementation

Use full MLWF support, periodic orthorhombic unpolarized occupied orbitals, k-only parallelism already supported by native MLWF exchange. Start with an impulse excitation and NVE ions, with pseudopotentials updated every step. Functional/run provenance remains checked by the DC reconstruction reader. The separately available projected LCFO optical response is unchanged and remains fixed-ion.

A nuclear Verlet drift produces endpoint positions. Build local/nonlocal pseudopotentials at midpoint positions for real-space electron propagation, and rebuild at endpoint positions for instantaneous-state force and energy. Include both initial and endpoint force evaluations. Exchange and rVV10 use the actual time-dependent mesh state/density. No moving-basis connection term belongs in this mesh representation.

## Gates and limits

Compare stationary nuclei versus moving ions, verify ions actually move, charge and fixed occupations, energy drift and trajectory refinement at dt,dt/2,dt/4. Use actual nonlocal projectors, not H-only local-projector fixtures alone. Preserve unsupported-path input checks. Finite-duration laser pulses require a separate external-field-work/endpoint-energy validation; do not enable them by inference from an impulse test. Long water trajectories, distributed real-space/orbital exchange scaling, truncated MLWF support and direct fragment dynamics remain later milestones.

This route uses DC for initial-state preparation and then conventional real-space Ehrenfest propagation. It does not claim that independently propagated truncated fragment forces have been certified.

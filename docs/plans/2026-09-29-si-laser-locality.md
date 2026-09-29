# Si laser locality timing plan

User requested immediate measurement without new per-orbital diagnostics.
Reuse existing solver/output; do not modify the physics code for this experiment.

- Si 4x1x1 conventional cells, 32 atoms, Gamma; mesh64x16x16, a=10.26 bohr.
- Same DC-LCFO preparation for PBE/HSE06/PBE0: four cores,8 x buffer points,
  64 fragment states,128 total LCFO states,300K DC chemical potential matching.
  Reconstruct64 occupied real-space orbitals for fixed-ion RT.
- MPI4 OMP4 BLAS1, runs sequential, saved executable and input hashes.
- Acos2, x polarization,3.1eV, I_wcm2_1=1e13 W/cm2,4fs width and1fs postpulse.
  dt=.08 au initially; probe .04 before production. No zero-field subtraction.
- Hybrid RT source ACE,99.9% adaptive radius, interval5,maxiter100,tolerance1e-7.
  Sparse-source experimental route stays off. Existing transported-gauge and ACE guards retained; invalid adaptive gauge can stop RT.
- Record runtime/RSS,existing radius maxima,pair counts,localization status,
  current/energy/norm. Do not infer orbital energy specificity or free-electron
  formation from these aggregate outputs. Small transverse box limits localization.
- Check GS convergence,short RT and time step first; then run finite pulses.
  Save all failures; do not overwrite prior results. No automatic restart on failure.

Runner and immutable results: work/si-laser-4x1x1 (workspace level).

GS trial log: initial/default Broyden alpha=.75 oscillated; preserved and stopped. Explicit alpha_mb=.2 reaches convergence. All final runs use alpha_mb=.2, mixrate=.2 and preconditioning. RT spatial decomposition is1x2x2 (required by the existing EXX pencil route), shared by all functionals.

Preflight: all three DC GS converged below1e-8. GS elapsed PBE10.47s,HSE06 24.86s,PBE0 30.50s. Initial16 steps(dt=.08) versus32 steps(dt=.04) currents differ by1.11e-9,4.28e-11,3.34e-10au respectively. This verifies only pulse onset, not full-pulse convergence. Production launched sequentially,2584 steps each. Record output-arrival timestamps for approximate half-fs timing bins; sampling/buffering limits bin precision. Final rank RSS uses child process resource usage.

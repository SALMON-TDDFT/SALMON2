# LCFO RT efficiency implementation plan

User-authorized optimization of WF updates, communication and matrix operations.
Preserve the existing physical approximation and impulse-first refresh schedule.

1. Split coefficient-space MLWF tracking from real-space source construction.
   Retained ACE steps track frames only; keep predictor rollback/cache semantics.
   Add a transport-only regression before implementing the new API.
2. Apply ACE using local LCFO rows: reduce F_local^H C_local, then form only
   -F_local times the reduced overlap. Avoid gathering all coefficients and
   redundantly multiplying the full global output matrix on every rank.
3. Fuse projection of local Hamiltonian output with addition of exchange in
   coefficient space; reconstruct once. Keep full-matrix fallback for invalid ACE.
   Test distributed action against dense ACE on arbitrary complex targets and
   unequal row partitions, plus projected-Hamiltonian equivalence.
4. Run native regression and Si128 before/after at radius9, ACE intervals1/4,
   dt0.16, 16steps, out_rt_energy_step=10. Isolate the energy-output change by
   running the old binary with exactly the same new input. Record wall/RT timing
   and E_inf with the existing common definition. One MPI job at a time.
5. Append results to the current note, preserve the paused long run and resume it.

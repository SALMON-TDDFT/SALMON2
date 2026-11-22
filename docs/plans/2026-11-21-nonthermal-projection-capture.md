# Nonthermal RT screening with adaptive projection capture

User selected actual RT occupations and requested implementation. Preserve the
same endpoint A, TB-mBJ/HSEX Hamiltonian, local RT potential and real-space
occupied-only propagated orbitals. No Fermi refilling or population rescaling.

Add capture metrics before W generation/publication. Measure each propagated
orbital's unrepresented norm and weighted electron loss, excited electrons and
holes. Start from existing nstate GS auxiliary orbitals. If capture is inadequate,
expand only projection workspaces, initialize additional orthogonal real-space
auxiliary vectors from missing RT components then deterministic grid vectors,
and reconverge the same frozen Hamiltonian. Bound expansion by
hsex_projection_max_states (0:auto up to min(2*nstate,grid dimension)); initial
nstate remains the GS restart size. hsex_projection_capture_tol defaults1e-8 in
orbital norm. Stop collectively at the cap, with previous W preserved.

Scope is scalar CPU, full-grid/k-only MPI, as existing internal pump screening.
Ordinary GS/TD projections without a screening snapshot retain their behavior.
Expansion beyond GS state count uses new auxiliary guesses, never extra propagated
states. Rebuild projection-only info, halo cache, CG scratch and eigenvalue arrays;
keep RT state/info/occupation arrays unchanged. GPU automatic growth is rejected
at the actual growth point until device mappings are supported; captured paths
retain their existing behavior. Restart reinitializes projection from GS and may
regrow at a later update; active W restart format remains unchanged.

Tests: independent capture/seed Fortran unit, finite/NaN/Infinity/overcapture
rejection; actual Si no-growth and forced growth from deliberately small initial
projection spaces, cap rejection before W save, same endpoint A/propagated count,
fixed-H population equivalence and restart currents. Audited native/full builds,
MPI/OMP checks including4MPI12OMP; no GS payload collection. New inputs in
functional namelist, defaults/read/broadcast/variables/ranges and manual updated.
Existing code/worktree and binaries are preserved; no automatic commit.

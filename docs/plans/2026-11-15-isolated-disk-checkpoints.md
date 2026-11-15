# Isolated E replay checkpoints

Continue approved F3c: public save_checkpoint/restore_checkpoint on isolated
states, accepted E and six named dipole observations only; GS supplied separately.
Reuse the strict bounded Zarr serializer through explicit private layout selection.
Keep periodic schema and layout unchanged. Isolated schema is
salmonpy.isolated-tddft-checkpoint/1; arrays E0, dipoles0, dipoles,
E_history, dipole_history. Input width 3, observable width 6. Bind initial files,
executable, dt, capacity, threads and sample position via existing settings hash.
No executable path is loaded from archive. Candidate replay remains atomic.

Test rebuilt-model continuation, corruption/budget/position mismatch rejection
before spawning, and periodic checkpoint regression. Validate actual local E
worker continuation, update paired guides, sample and AI JSON. No new physics,
MPI restart, runtime publication, or Actions.

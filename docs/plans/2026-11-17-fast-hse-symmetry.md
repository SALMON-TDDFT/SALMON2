# Fast-PR HSE symmetry is the primary path

Project lead direction: prioritize fast-PR HSE; symmetry is central. Defer the
standard k-point port. Preserve separate standard/reference physics, never make
conventional HSE a supported engine or change fast symmetry to match standard.

First implementation slice: optimized EXX refresh still checks the full grid
against global num_rgrid when symmetry is enabled. Replace this with a per-system
full-grid check using ngrid, primitive vectors and hgs; keep actual grids/maps and
existing exchange arithmetic unchanged. New test uses two independent symmetry
owners with different cell/omega/group/orbitals, unrelated poisoned globals,
interleaved refresh/action and close/reopen. Compare reduced representative
exchange actions/energy to the same fast kernel on independently enumerated full
k mesh. Test ACE and full-action profiles OMP1/2. No FFT/MPI owning handle clones.

Fast core supports symmetry but tddft_context still refuses symmetry and current
numeric checkpoints refuse nonempty symmetry maps. Do not remove those guards in
this slice; full TDDFT/FDTD needs field-little-group checks, complete map/history
checkpointing and MPI fault/restore tests. Document exact capability boundaries.

Verify old guard failure first; new external builds and HSE regression, source/
wheel/native hashes, paired guides/AI/checklist. Persistent local evidence, no
Actions; advisory budgets only, molecular work paused.

Implementation finding: native_owner_supported also excluded symmetry entirely.
Allow only explicit serial-k CPU/full-q symmetry owners after matrix-layout
admission; retain automatic GS binding and all TDDFT/checkpoint prohibitions.
Reject bad owned-grid/action layouts before state mutation. Two independent MPI
processes with COMM_SELF do not establish internal symmetry k-MPI2.

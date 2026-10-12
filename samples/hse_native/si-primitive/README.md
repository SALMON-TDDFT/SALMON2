# Nonorthogonal Si primitive-cell check

Two atoms, eight electrons and eight states, 12^3 grid, full 2^3 k mesh.
Use separate `gs` and `rt` working directories. Copy `gs.inp` and `rt.inp`
as their input files. Place `Si_rps.dat` from the existing Si sample in the
parent directory. Run GS before RT, with eight MPI ranks and one OpenMP
thread per rank. Set `SALMON_FACTOR_HISTORY=learned` for learned ACE;
omit it for ordinary ACE. The RT example has K=8, automatic warmup and
192 steps of 0.02 atomic time units.

These are small implementation checks. Converge the spatial grid, k mesh
and time step before physical production. DC and SR localization are off.

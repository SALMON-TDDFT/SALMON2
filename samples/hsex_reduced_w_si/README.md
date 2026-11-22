# Adopted finite-k Si response-compression preset

The user accepted this setting on2026-11-22 after the paired16384-step RT
spectra matched the same-GS hybrid reference. Use `gs.inp`, `rt.inp` and
`rt_zero.inp` as reusable input examples, copied to fresh scratch directories.
The solver remains explicitly selected through `hsex_w_solver='reduced'`.

| Setting | Value |
|---|---|
| Real-space grid |12³|
| Electronic / screening mesh |4³ /4³|
| XC / W |TB-mBJ / static occupied-state response|
| Boxes / cutoff |6³ /3Ha|
| Shared rank / representative-q batch |384 /4|
| Internal q symmetry |validated positive fractional-axis little-group permutations|
| GS / propagated states |8 /4|
| RT timestep / steps |0.05au /16384|
| Pulse |z kick1e-4au; separate zero-field control|

Create `gs`, `rt` and `rt-zero` under a scratch calculation directory. Copy
`gs.inp` to `gs/input`, and the RT inputs to the corresponding `input` files.
Copy `samples/exercise_04_bulkSi_gs/Si_rps.dat` into each run directory; its
SHA256 is8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c.
After GS/W completion, link `../gs/hsex-w.bin` as `hsex-w.bin` in each RT
directory. Both RT inputs read `../gs/data_for_restart/` and propagate only
four occupied orbitals in real space. Use the project's audited executable
and verify its recorded hash before reuse. Mac validation used8MPI×1OMP;
Fugaku runs must use at least12OMP and require separate capacity/performance
validation. Keep GS payloads on scratch; collect only derived results under
canonical calc2026.

The ordinary sample does not pay for an additional reference W. For an
accuracy check, enable `yn_hsex_reduced_reference='y'` in GS to also generate
`hsex-reference.bin` from the identical in-memory GS/Hamiltonian. Run a second
kick/no-field pair with that file linked as `hsex-w.bin`, then compare using
`tools/hsex/analyze_reduced_rt.py`.

Acceptance is by the paired physical observable: same cubic time window,
zero-field subtraction, first stable peak shift<=0.01eV, height difference<=1%,
and1–6eV Im-epsilon relative L2 difference<=1%. Response residual status2 is
not renamed convergence. The accepted Si measurement has W error0.0070%,
peak-height difference0.0069%, Im-curve difference0.0182%, and identical
2.765eV grid peak. This does not establish convergence for other systems,
k meshes, stronger excitation or metallic/SOC screening.

The current implementation is accurate here but is slower than the optimized
hybrid/dense reference on this small grid. Retain this preset and the paired
spectrum as the accuracy baseline when optimizing basis reuse, projections,
BLAS and parallel layout. This sample is not a speedup claim.

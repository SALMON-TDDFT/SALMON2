# Stationary Si: impulse versus finite-frequency Sternheimer

Prepare the same scalar TB-mBJ Si GS at Nr=12^3, Nk=2^3 with eight states
(four occupied, four empty) and the Si KY pseudopotential. Place its normal
restart under `../gs/data_for_restart/` and `Si_rps.dat` beside this input.
Use a tightly stationary GS; the diagnostic rejects a reference orbital rate
above 1e-6 Ha. The response and physical RT both retain only four occupied
states. This small k mesh is a solver check, not a converged Si spectrum.

Run this input on eight k-only MPI ranks with one OpenMP/BLAS thread. Its
`both` kernel writes independent and self-consistent current responses, density
chi0 probes, reduced operators and separately validated time companions,
then returns without physical propagation. For a shorter optical-only check,
select `hsex_sternheimer_kernel='selfconsistent'`.

In a separate directory run the same input with
`yn_hsex_sternheimer_optical='n'` for the ordinary impulse. Compare the current
with the companion on identical `nt,dt`, polarization, damping and SALMON cubic
window. Use both signs of a small impulse to check linearity. Do not compare
an infinite-time resolvent peak directly with a finite-window maximum.
Require `STERNHEIMER_OPTICAL_PASS`, and additionally `STERNHEIMER_WINDOW_PASS`
for a companion-time comparison, plus normal `end SALMON` completion.

Raw optical current includes its contact/DC term. The separate insulating
optical density-head conversion is restricted to an unexcited stationary
insulator; do not subtract a physical pumped/metallic Drude contribution.
The density probes do not publish a new full local-field W.

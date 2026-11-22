# Frozen-H single-q probe

This is an experimental response validation, not a new production W mode.
Use a previously validated Si GS restart with grid8^3, k2^3,8GS states,
and its unchanged HXT1 `hsex-w.bin` (box2). Keep the GS data in its existing
calculation directory. Copy the pseudo potential and initial W into a new
scratch directory and edit `directory_read_data` in `rt_probe.inp`.

The example propagates the4occupied orbitals for1fs under a weak4fs pulse.
At step800 it probes q=(1/2,0,0) without eigenstate projection. Rank128 is
an internal response-space cap. Inspect both the induced-density reference
error and the equation residual; a capped diagnostic is not accepted
convergence. Set the field amplitude to0 for an equilibrium control.

`yn_hsex_frozen_reference='y'` enables a validation-only512-point full-H
spectral reference. Set it to'n' to time the matrix-free diagnostic alone.
H-action ownership is MPI-distributed, but the coordinator is serial and
does not establish production scaling. No updated W or GS payload is produced.
Follow the repository source/object eigensolver audit before using a newly
built binary; reuse unchanged validated binaries by their recorded SHA256.

## Full-q update comparison

`rt_update.inp` is the experimental full-q version: 4 occupied RT orbitals,
2 fs propagation under a weak4 fs pulse, `hsex_w_solver='frozen'`, and endpoint
updates every1 fs. It requires Nq=Nk and power-of-two full k meshes.

Generate a matched initial W first: use the same RT input with zero field,
`nt=1`, `hsex_screen_step=1` and `hsex_screen_interval=0`. Keep the old initial
W and GS data in place for this preparation. The output
`hsex-frozen-updated-000000001.bin` has HXT1 layout and unchanged TB-mBJ
generating metadata. Use a copy of that output as `hsex-w.bin` in both new
comparison directories. Then run `rt_update.inp` for the updated trajectory.
For its fixed-W control set `yn_hsex_update_w='n'`, `hsex_w_solver='auto'`,
`hsex_frozen_rank=0`, `hsex_screen_step=0` and `hsex_screen_interval=0`;
retain4 orbitals, the same initial W and identical pulse/time settings.
Use zero-field counterparts when interpreting weak differences.

Every response must satisfy the conservative residual upper bound. In the
default strict closure, static Hermitian/reciprocal defects and assembled-W
eta dependence are logged and must pass `hsex_w_eta_tol`; no unconverged W is adopted. This mode retains
all finite-q local fields but uses the minimum-finite-q-shell insulating
head closure, not an analytic optical/metallic limit. It is an adiabatic
snapshot approximation, not general two-time RPA. A checkpoint is published
before the active W is replaced and endpoint exchange quantities refreshed.

## Coherent optical probe

`rt_coherent_probe.inp` samples a nonstationary occupied snapshot at 1 fs.
It applies the electric-impulse two-resolvent response, including the
nonlocal/contact derivatives, without an instantaneous-state projection.
Edit the existing GS restart path, provide its matching corrected-optical
initial W and pseudopotential, and run a separate zero-field control.

This diagnostic keeps the main RT and W unchanged. It reports IP optical
response at q=0 (not a requirement that the Bloch mesh contain Gamma).
It does not supply a general nonstationary static-W Gamma closure. Inspect
the two equation residuals, charge conservation and actual finite-A
derivative checks, and compare the low-frequency response to the zero-field
reference before interpreting Drude terms.

`rt_update_density_gamma.inp` selects the finite-q density Gamma closure,
`hsex_frozen_gamma='q2'`, for Nr12^3 and Nk=Nq8^3. Prepare a matching GS and
corrected initial W (boxes6, cutoff3 Ha); this input is not compatible with
the older Nk2 sample restart. It is an experimental update example, not a
claim that its eta or q mesh is converged. A failed static/eta/q check aborts
before replacing W. Inspect `HSEX_FROZEN_DENSITY_LIMIT_PASS` at every accepted
update. For weak-pulse comparisons, run the same input with updates off and a
separate zero-field control. Physical RT retains four occupied states.

The density route extrapolates chi0 at each finite q, then fits the complete
q^2 W box matrix and replaces the spherical Gamma cell once. It does not
remove a physical pumped Drude term. Nk2/Nk4 cannot supply the nested small-q
stencil and are rejected before any attempted screening update. The old
`shell` setting and the separate stationary optical option remain available.

### Streamed k/q workers

`rt_update_pairs_small.inp` is the Mac performance fixture: Nr8³, Nk=Nq2³,
8MPI×1OMP/BLAS,2q groups×4k workers,4occupied RT orbitals and one update at1fs.
It needs a matching existing GS and `hsex-w.bin`; it does not generate a new GS.
The finite-shell Gamma convention is kept solely for equivalence/performance,
not for E1 or optical-limit accuracy. Compare `hsex_frozen_parallel='legacy'`
(with q workers1) against `pairs`; optionally select `kinetic` preconditioning.
Each H-basis candidate is orthogonalized once in pair mode, while the original
equation still certifies convergence. Full physical q2 Gamma uses the separate
larger-grid sample and needs additional memory/cost validation.

### Experimental Hermitian static projection

`hsex_frozen_closure='hermitian'` explicitly projects the paired chi onto its
Hermitian part before inversion, with3eta samples and retained eta/positivity
guards. The raw discarded fraction is logged, and can be large for strong
coherent excitation. This is a defined additional approximation, not a claim
that the full nonstationary response is static or reciprocal. Keep `strict`
as the reference. Prepare a common initial W with these same rules before
comparing fixed and updated trajectories; do not attribute its equilibrium
finite-eta correction to the pump. The small-shell sample remains a diagnostic.

Initial strong-Si tests of this projection still fail eta convergence or
dielectric assembly/inversion. It is not adopted as a strong-excitation fix.
The coherent two-level regression demonstrates a real1/eta pole that survives
Hermitian projection; specify the auxiliary screening time/coarse-graining
before interpreting a static limit of a coherent state.

`rt_update_window_small.inp` tests the finite-window differential-screening mode.
Supply the same static-extrapolated Hermitian equilibrium W0 used by the fixed-W
control; the driver rejects a reference mismatch greater than 1e-6. It requires a
fresh zero-A equilibrium GS, pair workers, and excludes independent LDA W and RT
restart. The example uses Gaussian area A=0.5 fs (20.6706866682456 au), hence
sigma=0.398942 fs. Select `exponential` for decay time A with the same area.
Only the pump-induced chi change is windowed; equilibrium W and the coherent
occupied RT propagation are retained. This is an auxiliary screening model,
not a Drude collision term. Compare window shape and area before production use.

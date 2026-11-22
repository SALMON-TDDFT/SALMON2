# Explicit Hermitian static projection implementation plan

Goal: define an opt-in static approximation for nonstationary coherent response, without raising the existing strict threshold.

Keep strict default and raw-response solver/certified residual. New hsex_frozen_closure='hermitian' takes the Hermitian part after q/−q reality validation, before dielectric inversion. This is a defined conservative quadratic surrogate, not a demonstrated exact nonequilibrium W. Record discarded norm (half ||chi-chi†||/||chi||); do not claim a general dissipative/Drude or two-time theory. It preserves the real quadratic form of trial complex fields. No density/eigenstate projection or extra RT orbitals.

Use eta,eta/2,eta/4 and existing chi-level extrapolation even with the small shell fixture. Retain eta/W change, dielectric positivity, q reality, final W reciprocity and finite guards. Failure never publishes W. Strict mode retains original path. Exclude stationary optical option with the new mode; q2 closure remains optional and unchanged.

1. Failing native helper projection/quadratic/IEEE tests, then helper in hsex_frozen_limit.
2. Wire optional full-W closure, namelist/default/broadcast/dump/checks and adapter. Add full MPI stationary/coherent/partition tests and rejected eta/positivity tests.
3. Source/application-object prohibition audits, native GNU/OpenBLAS and Accelerate tests; HSE-on/off build checks.
4. Si strict regression and E1e-3/E1e-2 pilot at1fs, keeping eta tolerances. If passed extend fixed-versus-updated trajectories and update-interval checks. If rejected diagnose and report actual gate; do not raise it.
5. Review, manuals, sample and canonical results; no GS payload, no existing remote job/binary modifications.

Background: TD-aGW/static-COHSEX is used in https://www.nature.com/articles/s41467-024-50534-3 . That does NOT establish the explicit coherent Hermitian projection proposed here; treat it as an experimental additional closure.

## Completed outcome

The explicit Hermitian option is experimental and remains opt-in. Native algebraic quadratic/idempotence/IEEE, projected negative-epsilon, coherent scalar pole,1/4MPI projected and strict regressions, q2 Gamma and OMP checks pass. OpenBLAS/Accelerate pole results agree. HSE-on/off builds and source/application-object audits pass. Review added explicit nonfinite q-reality rejection. No new default or threshold increase.

Si strict weak W is bitwise identical to the previous validated binary. A common projected/extrapolated equilibrium W was prepared (relative change1.67246e-4 due eta rules). Strong E1e-3 pilots at coarse eta=.01,.005,.0025,.00125 still reject eta: indicators3.717%,7.555%,14.905%,25.358%. E1e-2 rejects eta at .01/.005 and dielectric assembly/inversion(stage33) at .0025/.00125; minimum eigenvalues are not logged and are not claimed. No rejected W is published. Strong post-update trajectories were therefore not launched.

A gap1Ha two-level coherent model with a real diagonal density vertex matches the independent spectral reference and has chi_H≈-1.13/eta, establishing that Hermitian projection alone need not create a finite static limit. The next physical choice is time-window/coarse-graining of the auxiliary screening density or causal finite-frequency W; neither is implemented by this experiment. Full causal strong-screening accuracy is not established.

Canonical results: calc2026/si-static-closure-20261123. No GS/restart payload, remote-job or fixed-binary changes; no commit requested.

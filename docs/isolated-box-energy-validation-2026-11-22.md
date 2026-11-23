# Molecular / LD / field box-energy validation — 2026-11-22

The opt-in observer reads a live pre-step backup and the accepted field without changing either. Fortran computes the six boundary fluxes, two definitions of field energy, LD work and external-current work. The query reduces 20 scalar values once within the owning field group. The test carries its result in the existing accepted-E gather. No extra application-wide collective was added.

Scope: existing C2H2 PZ ground state, 64 field steps, outer MPI4, one rank per field context, isotropic LD, dt=0.003 au, inclusive zero-based inner box (13,13,13)..(27,27,27). No new GS was calculated. Results and permitted worker copies reside in the existing calc2026/salmonpy-c2h2-feedback-20261122 campaign. The subsequent real MPI decomposition checks are recorded below; the manufactured kernel test separately covers partition additivity.

At the last fresh electronic energy (electronic step60 / field step61):

| Quantity | Hartree |
|---|---:|
| Modified discrete field energy change | 5.107273611394702 |
| LD work (storage plus dissipation) | 0.007353325279850679 |
| External-current work | 8.380216853189417e-6 |
| Outward six-face flux integral | -5.114635316891417 |
| Accumulated modified field residual | -9.137667543175966e-15 |
| Maximum single-step modified field residual | 5.1986193128072955e-14 |
| Molecular-plus-field residual, including removed own-vacuum work | 5.859442012558545e-8 |
| Grid source work minus mapped point work | 0 |
| Conventional time-centered field energy change | 5.128070452959051 |
| Joint conventional-energy residual | 0.020796900158768814 |

The two field energies are intentionally distinct because of Yee time staggering. Closure of the modified discrete energy is not continuum physical-accuracy certification. The molecular discrepancy is about 0.70% of its work; zero-field drift is retained, not subtracted. The conventional residual is also retained.

Source continuous and installed-wheel rejection/retry (step30) have identical saved physical arrays and all six-face records. The arrays also exactly match the previous midpoint run without energy observation. Source and installed-wheel suites each ran688 tests with one optional skip, no failures. NLCC XC variation and force tests passed with ranks1/2 (and their OpenMP variants). Source and application-object prohibition audits passed. Final Git-tracked-only prepared source matches the built source byte for byte.

Evidence: canonical campaign `box-energy-summary.json`, `box-energy-artifacts/`, `box-final-source-ld-continuous/`, `box-final-wheel-ld-retry/`. Library SHA256: `4abccd676d621c75eb78d1bb0849dae54165ab1d5f68db0f0cd0743a46917d8e`. Maintained solver inventory SHA256: `6765473bad4c8c948e3a059e6dd7af643a4e0f781033e7a4fc37f9fd2bbecbaf`.

Unsupported diagnostic profiles warn and return unavailable; they do not prevent the simulation. Tensor/subpixel media and larger or multi-node field decompositions remain open. No physical CI campaign was launched.

## Real FDTD domain decomposition

The same audited field library and existing C2H2 GS were reused. Each context now has two field ranks and one bookkeeping rank, so the two contexts use outer MPI6. The electronic worker remains one process/thread. No new solver build or GS calculation was needed.

| Field grid | All physical arrays versus one field rank | Cumulative discrete field residual (Hartree) |
|---|---|---:|
| 2 × 1 × 1 | exactly equal | -5.731616144620343e-15 |
| 1 × 2 × 1 | exactly equal | -2.7536964551399806e-15 |
| 1 × 1 × 2 | exactly equal | 4.846934613298165e-16 |

Every saved step's six-face records, both field energies, LD/source work and residuals agree with the one-rank result at rtol=2e-10, atol=2e-13 Hartree. Reduction-order rounding is expected for these integrated quantities. The mapped point lies next to or on the decomposed cell boundary, and the observation box crosses that boundary. The LD region is retained.

The installed-wheel run with grid (2,1,1) rejects interval30 and retries. Its physical arrays **and every six-face record** are exactly equal to the source continuous run with the same grid. Quantum histories also equal the one-field-rank baseline. The electron rollback still replays accepted history; this is not a direct quantum restart claim.

Reproduce the authorized fixture with `tools/test_molecular_feedback_mpi.py feedback BUILD_JSON NEW_OUTPUT`, the existing `--isolated-executable`, `--isolated-validation`, `--isolated-bundle` and `--isolated-input` arguments, plus:

```sh
mpiexec -n 6 python tools/test_molecular_feedback_mpi.py feedback BUILD_JSON NEW_OUTPUT \
  --isolated-executable WORKER --isolated-validation RECEIPT \
  --isolated-bundle GS_BUNDLE --isolated-input input.inp \
  --environment ld --source-amplitude .1 --energy-box 13 13 13 27 27 27 \
  --field-grid 2 1 1 --inject-reject-step 30
```

Use outer rank count `2*(NX*NY*NZ+1)` for `--field-grid NX NY NZ`. Omit injection for the continuous run. Output directories must be new. Evidence is in canonical campaign `box-energy-mpi-summary.json` and the `box-grid-*` cases. These are local development checks, not a multi-node benchmark or an added physical CI job.

# Molecular / LD / field box-energy validation — 2026-11-22

The opt-in observer reads a live pre-step backup and the accepted field without changing either. Fortran computes the six boundary fluxes, two definitions of field energy, LD work and external-current work. The query reduces 20 scalar values once within the owning field group. The test carries its result in the existing accepted-E gather. No extra application-wide collective was added.

Scope: existing C2H2 PZ ground state, 64 field steps, outer MPI4, one rank per field context, isotropic LD, dt=0.003 au, inclusive zero-based inner box (13,13,13)..(27,27,27). No new GS was calculated. Results and permitted worker copies reside in the existing calc2026/salmonpy-c2h2-feedback-20261122 campaign. This does not demonstrate real multi-rank field decomposition; the manufactured kernel test covers partition additivity.

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

Unsupported diagnostic profiles warn and return unavailable; they do not prevent the simulation. Tensor/subpixel media and real multi-rank field-face verification remain open. No physical CI campaign was launched.

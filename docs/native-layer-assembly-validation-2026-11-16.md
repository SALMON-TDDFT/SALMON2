# Native layer assembly validation — 2026-11-16

Base: b5b2053c; following changes were tested in the state-inputs checkout. Local Apple Silicon/macOS, Python3.14, host OpenMPI, one OpenMP/BLAS thread, at most four MPI ranks. GitHub Actions was not triggered. Persistent evidence directory: `~/.codex/outputs/salmonpy-state-inputs-20261116` (developer evidence, not a required user path).

## Input regression

Multiscale1D previously checked incident polarization at boundary sample times but not the negative-time history prefilled in vacuum. A short x-polarized pulse with vacuum(300,50),spacing5,dt.02,10steps reproduced acceptance despite nonzero forbidden prefill. It now raises InputError before state use. The corresponding z-polarized pulse remains nonzero and accepted. Valid-update arithmetic is unchanged.

## Configuration and packaged execution

Five layer-planning tests cover exact manual material-grid parity, ceiling thickness, repeated-material point identities, immutable reports, budgets without subsampling, invalid IDs/units/fields, execution rank contracts and LD neighborhood rejection. Planning does not prepare material states.

Built sdist and then wheel using `python -m build --no-isolation`, two build workers. Installed the wheel and executed examples with `python -I -m`:

- Four ranks: field1 + SBE1 + TDDFT2, selected from the installed mixed_layers.json. Native LD holds two nine-cell materials. Si starter GS is PZ,16³,2³k;128 steps,dt.003. Generated layout preserves material order and maps to transport IDs1/2 rather than manual10/20.
- Compared the generated layout to the prior manual native fixture after explicitly checking point material order: A/c,E,j maximum differences all **zero**. Native material update code was unchanged.
- Independent serial replay: SBE exact; TDDFT maximum j error **7.657177843178875e-19**, within the existing tolerance.
- Two ranks: field1 + SBE1, same LD regions, no TDDFT layer. Passed with deliberately nonexistent TDDFT runtime/bundle paths. Independent SBE replay exact. This guards against unnecessarily preparing TDDFT or an empty TDDFT MPI group.
- Installed-wheel Python suite:442 tests,1 skipped, no failures. Source Python suite against the installed kernel:442 tests,1 skipped, no failures. Manual synchronization and whitespace checks pass.
- Additional review found no remaining important issue after fixing conditional TDDFT preparation.

Evidence: `native-layers-final`, `native-layers-sbe-only`, `layers-final-comparison.json`, `layers-final-{installed,source}-tests.log`, `layers-release-build.log`.

Artifact SHA256:
- sdist: `e04c7136a34ebffe5e109ecdcab975099eba00fcaca9043a611f144b8b812858`
- macOS arm64 cp314 wheel: `e15e20706dc353d3be91e5fa990c778ea9b8519848974236817dbbe6e662e2f7`

## Limits

Finite rectangular geometry, one independent electronic point per electronic cell, explicit Bohr lengths and current native admission constraints. Not an infinite-film homogenization. The sample has four predefined model IDs; the planner accepts explicit user bindings. Rank budgeting does not estimate/enforce memory (aggregate memory unmeasured). No implicit MPI launch, GS preparation, asset download, representative-point reduction, unequal timesteps or restart. Repeated material bindings produce separate states through the existing session contract; this test is not unrestricted native reentrancy. Physical convergence and production scaling remain outside this connection regression; practical TDDFT still defaults to 1D recommendations.

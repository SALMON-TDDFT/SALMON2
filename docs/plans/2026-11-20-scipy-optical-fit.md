# SciPy optical fitting implementation plan

Approved by the user after the proposed SciPy/native division.

Goal: reuse SciPy bounded least_squares without adding a runtime dependency.
Architecture: optional optimizer="scipy"; native residual/seed/statistics; no silent fallback. Preserve native default for compatibility, use SciPy by default in explicit experimental-data fitting tool. Keep original catalog and all trials until comparison completes.

1. Add synthetic optional-SciPy contract tests; verify failure before implementation.
2. Add audited scalar-IEEE native evaluation/seed ABI; integrate bounded TRF multistart and optimizer/version/termination provenance.
3. Refit the retained CC0 tables explicitly; curate without hiding failures or stiffness.
4. Sync Japanese/English guides and AI JSON; run source suite and sdist-to-wheel installed checks. No physical-data fitting in CI.

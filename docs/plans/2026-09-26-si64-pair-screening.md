# Si64 fragment Wannier pair-screening validation

Goal: Validate locality screening on the user's historical Si64, 2x2x2 DC/DG fixture against full exchange.
Architecture: Start from frozen localized density factors Q (including sqrt(f/2)), measure symmetric pair-density overlap, and bound omitted exchange energy using max(V_G)=pi/omega^2. This is a diagnostic, not a production SCF operator. Test Hermiticity and negative definiteness of the resulting construction metric before considering ACE.
Tech stack: Python/NumPy diagnostic; existing Fortran exact Wannier backend remains the reference.

## Fixture and scope
- Use the original Si64 DG geometry, pseudopotential, grid, buffer and retained-state count once located. Do not substitute an invented input and call it the historical benchmark.
- All three axes are fragmented, so current DC restrictions require Gamma. Each core contains 8 Si atoms; buffer size determines actual fragment cost.
- Existing OneDrive/SALMON2 checkout is dirty and remains untouched. Continue the previous dc-hse-mlwf-ace branch in work/SALMON2-dc-hse, cloned locally from the clean committed implementation because the project root is not a Git repository.
- No speedup or SCF certification is inferred from frozen tests.

## Work
1. Add failing tests for zero-budget exactness, symmetric pair selection, discarded-energy bound, and complex action/metric diagnostics.
2. Implement samples/dc_hse/pair_screening.py for Gamma localized factor fixtures. Always keep self pairs; rank off-diagonal pairs by a rigorous L2 upper bound. Count ordered FFT actions consistently.
3. Evaluate full and screened Q actions on the same input; report energy error, theoretical bound, action error, metric Hermiticity, and metric minimum eigenvalue. An arbitrary target-dependent screened action is not a linear Fock operator; do not install it into LCFO or SCF.
4. Write an explicit NPZ input contract and tests. Retain references to raw input and units.
5. Await the historical Si64 input path, then reproduce its DC baseline and export converged localized factors. Sweep error budgets, buffers and localization convergence separately. Production pair pruning needs additional derivation/testing if the metric tests fail.

## Progress
- Historical fixture located through the referenced chat in SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement/tests/dg/data/si64_overlapping_wannier_rt.
- Preserved cell 20.52 bohr, total grid 32^3, buffer 6, retained states 400. Added coordinate and pseudopotential hashes.
- Frozen Gamma diagnostic and independent dense real-space normalization/complex-gauge tests added. No production approximation is enabled.
- Independent reviewer confirmed conjugation, spin factor and Parseval energy bound.
- Full-HSE convergence, real Si64 factor export, physical error sweeps and production pair-screened SCF remain outstanding.
- All nine new/existing compiled/Python unit tests pass (including five screening tests). The historical Si64 input was accepted by the committed MPI HSE executable; all eight first exchange refreshes completed, then the 120-second pilot limit stopped it before SCF convergence.

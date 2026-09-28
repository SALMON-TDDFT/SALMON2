# Reproducible distribution measurements

Goal: measure native exchange and LCFO dense-solver time and per-rank peak RSS, separately from scientific water/solution validation.

Use existing production kernels and peak_rss.c. Each measurement runs in a fresh MPI process set; no replicated oracle arrays are retained. Native synthetic orbitals are generated analytically on local grid rows and orbital columns. Verify ACE source reproduction and decomposition-independent exchange energy. Dense Hermitian input is generated one small column block at a time and discarded; verify the solver's distributed eigensystem diagnostics. Report baseline/max-rank/sum-of-rank-peak RSS and phase times, with repeated runs and complete metadata. Sum-of-rank peaks is not simultaneous node RSS.

Start with locally manageable deterministic workloads and 1/4/8 ranks. Strong-scaling comparisons keep workload and algorithm fixed. Do not infer giant-water performance from synthetic states or an all-eigenpair synthetic matrix. Do not add memory caps or change production algorithms based on a tiny-case speedup. Record measured regressions as well as gains. Existing exchange and dense comparison oracles remain the correctness reference.

Completed: 32³ grids with 16/64 synthetic wave-packet orbitals and a 1024-dimensional complex matrix; three fresh-process repetitions per layout (81 launches). Results and limitations are in `docs/hybrid-distribution-measurements-ja.md`. Localization reaches the three-iteration limit; these are not converged-water timings. Production code is unchanged.

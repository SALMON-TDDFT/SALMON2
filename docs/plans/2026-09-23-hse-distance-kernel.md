# Common distance kernel exchange implementation plan

**Goal:** Evaluate a shared electron-separation cutoff without truncating individual Wannier orbitals, and measure accuracy, ACE compatibility and runtime on Si 4³.

**Architecture:** Stream primitive-cell density-matrix row blocks. Transform the occupied density matrix over the uniform k mesh to cell translations, multiply by the original sampled periodic screened kernel (optionally zeroed outside a minimum-image radius), and transform back to obtain a common per-k exchange action. Keep the scalar kernel real and inversion symmetric. All sources remain intact. No MLWF gauge optimization is needed for this alternative. Retain the full-radius case as an exact independent control.

**Tech stack:** Python, NumPy, existing HSE samples and ACE, unittest.

## Approved design and scope

The user approved the common relative-distance cutoff after discussion of locality and ACE amortization. Fixed ions, Si8, primitive grid12³, k4³. No k-convergence study. The existing full-support response job is left unchanged. This design uses kernel locality and lattice translation symmetry rather than imposing orbital-specific supports, which failed ACE in the preceding experiment.

## Tasks

1. Add `test_distance_exchange.py`: missing API must fail first. Verify full-radius Bloch parity on a shifted/permuted k grid, finite-cutoff Hermiticity on arbitrary targets, fixed-source linearity, occupied-unitary covariance, energy finite differences and invalid geometry/cutoffs.
2. Add `distance_exchange.py`: blocked density-matrix algorithm, optional radius in bohr, exact discrete periodic kernel normalization, no dense full-supercell operator. Skip primitive column blocks that cannot intersect the radius. Return separate construction/application timings and workspace sizes. The radial mask is a hard cutoff for an initial fixed-radius experiment, not a redefinition of analytic HSE.
3. Benchmark full radius and several cutoffs on the persistent HSE ground state against the optimized full MLWF action. Verify occupied metric Hermiticity and positivity through unmodified ACE, and ACE interpolation. Record repeated timings, energy/action error and rebuild/apply break-even.
4. If compatible, expose an optional backend on the independent HSE functional and verify its energy gradient. Run a short zero-field/impulse PT-CN pilot with fixed cutoff and separately report initial-state stationarity error (GS was solved without cutoff). Do not switch the ongoing reference run.
5. Review, full sample tests, save raw results and report the achieved scope. Do not equate changed-kernel results with unmodified HSE or claim k scaling from one mesh.

## Scaling assumptions

At fixed primitive grid G and occupied count b, streamed density-matrix formation is O(Nk G² b), transforms O(Nk log Nk G²), and action O(Nk G² b_target), with row-block working memory O(Nk block G). A finite radius can skip impossible spatial blocks; for a radius exceeding the primitive-cell minimum-image diagonal it cannot skip primitive columns, although translation entries are zero. This must be measured honestly rather than attributing all speedup to truncation. ACE application remains O(Nk G b²).

## Ledger

- Plan approved through user's “やってみましょう”; implementation proceeds inline with test-driven-development.
- Existing TDCDFT feature branch is used. Changes are isolated to sample files and optional backend selection, not the running process.
- Tasks1–2 complete: missing-API tests observed red, then direct Bloch/finite-cutoff Wannier, gradient, Hermiticity, covariance and onsite tests green. Density-matrix blocks use the same periodic sampled kernel.
- Task3 complete: full exchange7.692s versus12.212s MLWF, relative agreement7.55e-15. Rc16 error3.25e-5 and0.0160meV/atom; ACE passes. Rc6 kernel is indefinite despite occupied ACE passing.
- Task4 complete: optional backend wired after failing integration call, then real zero/impulse Rc16 pilots passed. Impulse25.688s versus37.823s old pilot. Added a red/green positivity gate at HSE integration; diagnostic apply remains allowed for negative-kernel benchmarks. No MLWF update required for this backend.
- Task5 complete: independent phase/normalization review passed; added mesh3 regression and qualified default-backend docstrings. 67 tests pass. Results and remaining limitations recorded in docs/results/si-hse-distance-kernel/README.md.
- Ruling: report speedup as primarily algorithmic, not sparsity from cutoff. Rc10/16 do not skip primitive columns and dense k FFTs still process zeroed translation entries. No unsupported larger-k speedup claim.
- Ruling: do not interpret the one-step cutoff trajectory as a spectrum; original full-HSE GS is not stationary for Rc16 (zero-field residual9.63e-6Ha). A consistent cutoff GS is needed before a long response run.

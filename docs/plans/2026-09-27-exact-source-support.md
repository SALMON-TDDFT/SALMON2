# Exact source-support exchange loops

Goal: reduce zero-region memory traffic in the dominant exchange construction without changing the operator, FFT grid, cut radius, Taylor4, or MPI layout.

Design: construct exact nonzero source row indices once per translated source. When fewer than half of grid rows are nonzero, form pair products and accumulate exchange only on those rows. Clear the full FFT work array before each pair. Dense sources retain the original contiguous loops. Count actual product/accumulation grid points for regression and measured work reduction; no magnitude threshold. The internal dense switch is for direct reference testing.

Alternatives: sparse propagation changes approximation/gauge semantics and needs a separate R_prop convergence study; increasing ACE/U intervals changes error; this exact kernel optimization is a bounded step under the approved efficiency work.

1. Extend direct-convolution exact-pair test for compact versus dense execution, work counters and finite values, including tiny components, empty/disjoint support, FFT batches1..8 and OMP1/2/4.
2. Implement compact source loops with disjoint target-column ownership; keep dense fallback for broad support. Build full HSE binary.
3. Native direct/default MPI regression; preserve ACE counts, current/density/energy.
4. Sequential before/after benchmarks with identical C32/C64/C128 GS and inputs. Retain only if useful; single samples remain provisional.
5. Review, update root development note and bundled data, commit/push accepted result. Sparse propagation and dense U costs remain explicitly unresolved.

Verification ledger: initial regression failed on missing compact-support API. Direct-convolution tests now pass for batch1..8/OMP1,2,4, with product and accumulation counters, tiny finite components and disjoint support. A second rectangular two-k-point case checks translated compact sources against independent direct convolution and the dense path. Six existing Wannier/Bloch-reference tests and direct/default native MPI regressions pass. Read-only review found no blocker; its accumulation-counter and translated-source test requests were implemented. Identical IEEE behavior for nonfinite input is not claimed. Six sequential before/after runs completed. C32 RT21.388→20.506s, C64 40.171→32.216s, C128 52.892→52.032s. Each single sample; variable background load. C32-normalized weak efficiency improves for C64 but not C128. Product points94,371,840→51,339,264 and accumulation points60,817,408→39,763,968 per fragment; FFT count7424 unchanged. No new accuracy approximation; current/density/energy comparisons pass. Keep the exact work reduction but do not claim general weak-scaling improvement.

# Symmetry-aware exchange worker layout

Approved by the user on 2026-10-28: enumerate MPI factor pairs, minimize the maximum grid-pair work per rank, then maximize active ranks; retain explicit overrides.

The r dimension is the actual `nfirst` orbit count computed from the selected symmetry operations. The r-prime dimension contains all `ng` grid points. Neither is the number of k-point owners. Work ownership of k points remains unchanged.

Add a small integer-only layout module. For every factor pair of the exchange communicator size, reproduce the worker's bounded row batch and column chunk, including its 64 MiB density scratch cap. Score the exact cumulative blocked-cyclic row and column counts, not a fractional average. Choose minimum maximum pair count, maximum active ranks, then fewer column groups as a deterministic tie break. Report the chosen geometry, representative count and predicted work range; this is not a runtime speed optimum claim.

Input value zero means automatic for each worker setting; positive values preserve explicit selection. Negative values and non-divisible manual group counts remain collective errors. Defaults become zero. The selector runs only after `nfirst` is known, before allocating scratch and splitting communicators. It preserves the contraction, kernel, symmetry operations, finite guards, output ownership and checkpoint format.

Validate deterministic choices against an independent exhaustive assignment oracle, memory-cap/overflow/prime/partial/empty cases, and 2/4/8-rank scalar and spinor actions versus explicit layouts and direct exchange. Preserve NaN/Inf and rank-mismatch rejection. Audit source before compiling and every linked application object before execution.

Current remote jobs keep their validated executable. Nr16 Nk8 job 52847533 receives only an explicit equivalent configuration while QUE: row16/column80/tile80. Existing Nk16 job 52911342 uses row27/column1024/tile1024. No production binary replacement or scheduler operation is part of implementing automatic selection.

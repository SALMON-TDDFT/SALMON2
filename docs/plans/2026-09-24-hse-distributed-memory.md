# Distributed native HSE exchange memory

User requests distributing everything practical, not just a memory audit.
Implement inline in the existing TDCDFT branch. Preserve the physical kernel,
ACE equations and integration methods; avoid changing the active reference job.

Design: retain source, target, action and phase arrays only for owned k points.
For each round of spatial row tiles, build density tiles on k owners using
BLAS. All-to-all transpose to row owners, perform existing k FFT/kernel/inverse
FFT, transpose back to k owners, then apply to local target orbitals. Equal
padded messages support uneven contiguous k layouts. Every rank executes the
same collective schedule, including idle row owners in the final round.
Only tile buffers contain all k; reduce tile width with process count to bound
these buffers. The small sampled kernel and k ordering remain replicated
because every row owner needs them. Avoid global source/action reductions.
Free Taylor initial/midpoint factors immediately after final propagation.

Alternatives: streaming broadcasts use less tile memory but excessive repeated
collectives; retaining global orbital gathers is simpler but violates the goal.
The two all-to-all transposes retain efficient BLAS/FFT without orbital gathers.

1. Add independent Python/full-kernel parity MPI probe for distributed exchange,
   arbitrary targets, shuffled shifted k mesh, uneven layout, partial/idle rows.
   Verify failure before implementation.
2. Add distributed kernel API with exchange callback, local phase initialization;
   use existing serial reference as oracle, not as the production fallback.
3. Connect native refresh and full-Taylor audit using local arrays; remove all-k
   source/action allocations and release Taylor stage temporaries.
4. Build, numerical regression, MPI1/4/8 physical parity and restart; compare
   old/new peak RSS and RT timing sequentially from the same initial state.
5. Document achieved memory reduction, communication/time tradeoff and remaining
   replicated buffers. Commit verified work; publish alongside existing notes.

Completed: independent kernel MPI1/3/8 parity and rank-local invalid-data tests;
79 standard numerical tests, six Si/TDCDFT CTests; native MPI1/4/8 and bitwise
restart. Fresh read-only review found no blocker after collective validity
and finite checks were added. Sources/targets/actions/phases remain k-local.
Initial b=2 doubled Taylor time; final max(1,min(16,64/P)) amortizes collectives.
Observed MPI8 summed RSS reductions: Taylor 36.8%, PT-CN 34.9%, full Taylor
57.6%; step times increase 39%,22%,31%. Final old/new orbital relative errors
below 1.78e-14. Full audit no longer constructs unused initial/midpoint ACE.
Evidence: docs/results/si-hse-distributed-memory/. Small sampled kernel and
whole-k row tiles remain replicated/row-owned work; no ideal 1/P memory claim.

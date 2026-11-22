# Static finite-q RHS batching

Measured Si Nr12³ Nk=Nq4³ boxes6³ cutoff3Ha,8MPI1OMP: after form-factor reuse,
W169.04s, finite-q solves110.72s, assembly14.53s. Batch static occupied-band
charge vertices through the existing occupied projection and target rotation,
then reuse unchanged per-band tridiagonal solves, residuals, Gram and spin
weights. Pack only nonzero source occupations. Preserve q OpenMP distribution,
thread-local BLAS policy and pump/refinement arithmetic. Bound additional wide
buffers to64MiB/thread, except the unavoidable one-band minimum. Do not change
physical parameters, solver defaults, MPI layout or existing executables.

Alternatives: matrix-free finite-q was slower at this size; Gram Hermitian BLAS
has lower potential because target rotation scales N²G. RHS batching is the
bounded first measurement. Test batched versus one-band calls, mixed q counts,
occupations, projection and eta/pump compatibility. Record optional local
performance comparison against repeated one-band calls without CI timing limits.
Reject or narrow optimization if actual Si W fails tolerance2e-12+2e-11|W| or
end-to-end/phase timing shows no useful improvement. Audit every build and object;
use new build directories and executable hashes. No GS payload collection.

## Validation details

Native Accelerate: nb7 forces full and partial four-band panels, with q-dependent
G counts, skipped occupations, a fractional weight, OMP1/2 and an all-zero source
set. Independent diagonal-resolvent and repeated single-band comparisons pass.
Generated capacity fixtures retain active counts80/73/40 and pad mg to1024/1400/
2048, driving the memory-bound panel sizes to3/2/1; all pass. Fixtures are offline
test variants only, compiled with -fmax-stack-var-size=65536 to avoid placing
large program-local reference arrays on the main stack; production sources and
flags are unchanged. All linked fixture and kernel objects were audited.

Scalar uneven5MPI, spinor, optical/local-field/zero-limit and inverse IEEE
regressions pass. Unchanged original-H refinement passed an independent90-digit
oracle at eta5e-4/1e-6/1e-9/1e-13 with4MPI×12OMP (maximum tolerance ratio7.85e-4).
The unguarded standalone sequential-OpenBLAS OMP2 probe failed before the change;
macOS Accelerate passes. The full SALMON build has the existing
HAVE_HSEX_APPLE_GNU and HAVE_EXX_OPENBLAS_THREADS guards enabled. Do not claim this
preexisting standalone-backend limitation was fixed by RHS batching.

First actual Si8MPI measurement: finite-q solve110.7156→104.3307s;
whole W169.0387→163.7743s. W-header unchanged and full W within2.22e-16 of the
dense reference; optical epsilon and Gamma quadrature unchanged. A second
reference/new sequential pair is run because this small timing difference needs
checking against runtime variation.

Second Si pair: reference161.7788s (solve105.9449s), batch161.1747s
(solve103.3170s). Two-run means: whole165.4088→162.4745s (1.0181×),
finite-q108.3303→103.8239s (1.0434×). Both pairs improve, but broad scaling and
Fujitsu performance remain unmeasured. Retained as a modest bounded optimization;
do not quote the standalone1.2× as full Si acceleration. Actual full SALMON
4MPI×2OMP also completed, with team2 and W numerical equivalence verified.

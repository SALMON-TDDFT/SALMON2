# Frozen response k-sum and exact residual contraction

**Goal:** Continue the approved frozen-response validation by comparing the full
k-summed induced density, then reduce exact residual-evaluation work.

**Architecture:** Fold each k-pair density into the canonical q-periodic gauge,
sum with Ng/(Nk*volume), and compare to the same full-H reference. Independently
check its Fourier contraction against the existing physical trace. Keep both
per-pair and summed errors. No production W replacement is introduced.

Residual factors [Q,HQ,Q,U] and [Q,Q,HQ,V] can be combined exactly into
[Q,HQ,U] and [Q,HQ,V], with coefficient blocks [zC,C,0;-C,0,0;0,0,-S].
The Frobenius norm needs only QR's R, so omit construction of unused Q.
Retain the previous roundoff safety bound. No basis-admission tolerance or
acceptance threshold is relaxed.

1. Add failing MPI test for optional k-summed density errors and trace closure.
2. Implement fold/sum and independent-reference comparison; test serial and
   MPI including empty ranks and failure exits.
3. Compare reduced exact residuals to explicitly formed small-grid residuals
   at capped and converged ranks, then implement the contraction reduction.
4. Audit/build focused tests and HSE SALMON; retain old validated binaries.
5. Repeat matched equilibrium and laser-snapshot cases at rank128 and a strict
   convergence cap, recording density error, residual, H actions, time and RSS.
6. Save derived results under canonical calc2026 without GS data. Do not claim
   scalable production W or full spectral agreement from these diagnostics.

## Validation outcome

Completed the native small-model and MPI1/4 checks, including independent
explicit residuals at converged and capped ranks. The nonzero-wrap MPI case
closes the canonical Fourier density contraction against the physical trace.
Full Si equilibrium and weak-laser snapshots used Nr8^3, Nk=Nq2^3,
q=(1/2,0,0), 8 MPI x 1 OMP. Ranks64/96/128 map accuracy versus cost;
equilibrium cap512 achieves the unchanged 1e-6 equation-residual criterion.
Response values, basis sizes and H-action counts agree with the previous
implementation. Three interleaved rank128 runs establish the timing comparison.
No production W is published; response-space growth remains unresolved.

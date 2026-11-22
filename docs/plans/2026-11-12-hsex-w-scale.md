# Screened-exchange strength diagnosis

Use the requested lambda=0.5, 0.7, 1 comparisons at fixed Nr12/Nk16,
boxes6/cutoff3 and the same existing microscopic W and GS. Keep TB-mBJ updated.
The existing all-q diagnostics remain independent of lambda.

Add dimensionless `hsex_w_scale` to `&functional`, default 1, range [0,1].
Apply it once to -W(P-P0) when forming the exchange kernel. All actions,
current derivatives and induced energy then share the scale; neither W files
nor reference densities are pre-scaled. Record scale in HXC3 RT checkpoints;
read HXC2 as scale 1 and reject changed scale on restart collectively.

1. Extend the actual HSEX operator/current test with lambda cases and checkpoint
   roundtrip/mismatch/legacy and nonfinite rejection. Observe missing interface.
2. Add namelist/default/broadcast/log/validation and the single kernel multiplier;
   update checkpoint metadata and both manuals.
3. Audit source before builds and all application objects before new execution;
   verify GNU serial/MPI, actual action/current/energy identities and restarts.
4. Inspect the existing off/HSEX and finite-q jobs. Reuse full-duration matching
   data when available. Build a separate Fugaku executable; preserve existing
   fixed jobs and binaries. Validate 4 MPI x 12 OMP before lambda RT comparisons.
5. Save only RT and extracted screening results/provenance in canonical calc2026.
   No GS collection. Do not infer excessive binding from lambda nonlinearity
   alone, or insufficient screening from a small inverse dielectric head.

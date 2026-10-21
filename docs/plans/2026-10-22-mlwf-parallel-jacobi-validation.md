# Experimental distributed parallel MLWF Jacobi

User approved the independent-pair round design on 2026-10-22. Existing fixed jobs and binaries are preserved. New branch mlwf-parallel-jacobi-20261022 is experimental pending native Si512 validation.

Core: block-cyclic matrices remain distributed; batch 2x2 blocks and angles per round; process-row/column mate exchanges apply commuting sparse rotations. Identical rotation criterion, functional and residual threshold; PZHEEVD seed and LCFO driver replace prohibited PZHEEV. MPI errors fail-stop before outputs are consumed. Empty tile padding initialized. Final gauge unitarity checked at1e-10.

Local numerical validation: MPI1/4/6/8 N1,3,8,17; final MPI8 also N65 block32, reversed ranks, forced skips/NaN rejection. Si64 elongated8x1x1 sameGS/128orbitals/8MPI: both10sweeps, finalgradient serial2.00465e-8/parallel1.30027e-8 with1e-7 threshold. Spread610.2916039695147/610.291603969441. Full HSE06 action on2 fixed probes old/new relativeL2 2.99681e-15; R10/norm.9999 masked old/new7.93163e-7. Approx-vs-full relativeerrors .001363827221/.001363822806. Four own-gauge workerroute comparisons are zero and all8ranksexit0. Those comparisons alone are not a cross-gauge accuracy test.

MLWF-only wall time in the instrumented pair2.976701s/.260823s (~11.41x); single local pair, not a native-scaling claim. Diagnostic initialization and EXX action export are excluded from this timer. Native measured RSS+Huge and speed pending.

SourceSHA bf1ed7cae2bcc21da244545263777f8c01597d793de0fb7eb91a33ccb316baf1. Native binarySHA c19c76055d06d74b59d3d3a392d9e9e27f60582dc5139ab3c568eaf00d3cdc78. Native52614422 submitted once64nodes256MPI12OMP6h, new algorithm round test then sameSi5121024MLWF/R10 comparison; full preflight514GS sizeSHA and all actual linked applicationobject checksPASS. Native precision/capacity acceptance pending. Original reference52613213 unchanged. No automatic retry or threshold relaxation.

See workspace work/mlwf-parallel-jacobi-20261022 for reproducible compile/link commands and Calc2026/mlwf-parallel-jacobi-20261022 for nonGS raw outputs/actions/provenance. No GS download. Cloud sync unverified. No push or production promotion.

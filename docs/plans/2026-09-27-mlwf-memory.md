# Initial MLWF memory reduction

Goal: remove replicated six-direction link matrices and bound construction scratch without changing the initial localization problem. This is a prerequisite step, not certification of the 8000-atom production case.

Design: form N x min(64,N) overlap tiles by ZGEMM with conjugate transpose, reduce only to spatial root, and fill the same six-direction root array. Store three independent directions, reconstruct each negative direction as the adjoint on root. Non-root ranks allocate no global link array. Preserve snapshot wire format, MV algorithm, seed, tolerances and runtime settings. Release link storage immediately after localization. Avoid changing the root seeding algorithm in this bounded step.

Alternatives: lowering occupation count or localization accuracy changes the physics; merely reducing to root still leaves N² local construction scratch. Tiling bounds local scratch to O((grid+N)*64), while preserving exact matrix products to roundoff.

Validation: MPI2/4 independent dense reference with complex data, asymmetric 3D positions, empty local rows, tile tails and varied widths; assert root-only storage and explicit allocation counts. RED absent module/API, then GREEN. Native direct-WF regressions and same-GS C128 before/after including binary initial links, current, density and energy. Read-only review. Publish measured/theoretical allocation reduction and explicitly retain the production hold: root QR/gather, MV workspace and RT replicated orbital matrices remain.

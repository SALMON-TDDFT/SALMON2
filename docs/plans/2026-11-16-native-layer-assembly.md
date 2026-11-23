# Native mixed layer assembly

User-approved scope: fix incident validation first; build mixed calculations from material IDs, order, thickness and parallel settings. Preserve existing native updates. Completed implementation is `salmonpy.layers.assemble_native_layers`.

- A layout references explicitly resolved material objects, without duplicating their physics parameters. Existing catalog/initial-state loading stays separate. User response implementations enter through MaterialBinding.
- Versioned layout uses Bohr, propagation axis, origin cell, ordered uniquely named layers, material IDs/null vacuum, thickness and transverse cell counts.
- Axial thickness rounds upward to whole grid cells; report input and actual thickness and exclusive bounds. No representative-point subsampling or implicit gaps.
- One independent electronic state/group per occupied electronic cell; repeated material IDs reuse immutable initialization, not evolving state. LD remains on the field group and retains native E/H update order.
- Validate ranks, backend execution, source support, absorbing boundaries and LD/electronic neighborhoods before preparing native states. Rank budget is not a memory cap.
- Return calculation, bindings, points, MPI plan and detached JSON report. Do not launch MPI, download assets or prepare GS during planning.
- Installed sample accepts `--layers` and `--rank-budget`. Its built-in model dictionary is a fixture; public assembly accepts arbitrary explicit bindings.

Acceptance: longitudinal prefill regression; geometry parity, ceiling and ID/rank rejection tests; sdist-to-wheel; installed mixed four-rank equality to manual layout; SBE/LD-only two-rank execution without TDDFT files; source/installed Python suite. Physical convergence, whole-session restart and large-scale memory estimates remain separate.

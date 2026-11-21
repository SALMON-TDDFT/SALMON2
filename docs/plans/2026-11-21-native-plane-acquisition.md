# Native FDTD plane acquisition (2026-11-21)

Prior integration committed/pushed ac0354f7. New profile is fresh batch FDTD.
PlaneMonitor axis/nearest-grid position generates an explicit sidecar asset;
management-node preparation remains data staging only. main_maxwell opens/closes
per-fe diagnostic resources. eh_calc appends only after existing saved E/H are
prepared. Owning ranks alone write fixed-endian-detectable stream tiles for the
plane and preceding normal plane; no new collectives, no all-point sampling.
Embedded contexts do not open these resources, preserving transactional IO/state
contracts and numerical checkpoint bytes. They warn if given batch-only planes.

Python verifies clocks/coverage/geometry, assembles tiles and calls the native
periodic Yee2-edge E/4-face H spatial interpolation. H is already at E time;
initial fresh fields are recorded at0. Final sample is the last observation
stride boundary. FDTDResult reference/budget methods use collected plane areas;
material work/storage remain explicit, never assumed from insufficient data.
Old runtime or observation IO/resource failure warns and retains completed
calculation/point outputs. No auto-build on import or preparation.

All3 axes and two simulated tile owners cover colocation/clock/incomplete-data
contracts. Actual native grating: monitoring off/on and vacuum controls reproduce
128 old point-output files bytewise. Full EM checkpoint bytes also match.
Serial/transverse MPI2 plane results/budgets match (~0.0107% residual). Another
normal MPI split places the plane/neighbor on distinct ranks: nonzero E/H are
bit-identical. This is local validation, not optical convergence/multi-node proof.

Native source/app objects were audited; unchanged prior executable reused by
SHA only. Live inventories and all source-profile references stay synchronized.
Git-only export/preparation/build is used, preserving history/fixed binaries.
Results: canonical calc2026/salmonpy-auto-native-planes-20261121. No GS collected
and cloud synchronization not asserted.

Final source and installed sdist→wheel suites: 678 tests, one optional skip,
no failures. Installed wheel rereads serial/MPI native streams and reproduces
both budgets. Forced native observation IO failure preserves original point
output and marks diagnostics unavailable. Prior ac0354f7 automatic LD wheels
run 44867826995 succeeded. Validation report/logs live in the canonical result
root above; follow-up implementation remains uncommitted at this checkpoint.

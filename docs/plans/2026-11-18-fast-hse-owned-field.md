# Owned fast-PR symmetry field admission

Continue the project lead's fast-PR/symmetry priority. Implement a status-returning
uniform A/c little-group check on each initialized value-only symmetry map,
plus an exx_native wrapper tied to the explicitly bound system owner. No global
sym_sub fallback. Validate map readiness/layout/finite Cartesian rotations; use
normalized A direction for scale-independent tolerance1e-10, zero permitted only
after validating the map. Return0 compatible,1 invalid/unready,2 incompatible.
State/model/map arrays are read-only and no handles/caches are cloned or updated.

Wrapper rejects stale model geometry/operation identity; disabled symmetry only
permits finite fields with a consistent bound non-reduced owner. Wire the local
check before TDDFT step numerical backup/propagation, agree existing admission
status collectively. Preserve current symmetry prohibition in rt_open and in
checkpoint support until complete symmetry persistence/MPI tests exist.

Fail-first native probe: independent C2z/mirror-x/rotated-cell groups, global
poisoning,zero/tiny/huge/nonfinite fields,tolerance boundaries,incomplete/malformed
maps,stale/closed/wrong owners and state preservation. Existing exchange core16
cases and non-reduced mixed/protected step/direct restore remain unchanged.
New complete override exx_symmetry file,strict manifest chain,new external build,
source/wheel tests,paired guides/AI/checklist,persistent evidence,no Actions.
No general initialization recovery/public C ABI/physical-convergence claim.

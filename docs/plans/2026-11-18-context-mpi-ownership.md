# Per-context initialization MPI ownership

Next bounded slice of the approved state/resource ownership backlog. Tag child
communicators created during synchronous rt_open with a unique attachment-local
owner ID. Context close releases request caches/FFT/XC owners first, then its
children, preserving sibling state communicators and the host parent. Reuse
registry slots with distinct owner tags; old owner close never frees new children.
Default unscoped/standalone paths remain attachment-owned. No general arbitrary
input/threaded reentrancy or runtime-step resource guarantee is implied.

Density symmetry construction currently keeps the last context's r communicator.
Use stable COMM_GROUP_SELF for its existing k-only size-one r subgroup, or rebind
explicitly to the live owner's r communicator on context entry; select the
smallest change that preserves original arithmetic and generalized grid paths.
Test actual MPI ownership counts and repeated allocation/close, mixed sibling
propagation/save/restart, symmetry and full-k, old baseline trajectories. Update
paired guides/AI and newest checklist. Explicit new external build/source dirs,
manifest updates, source/installed/package checks; no Actions. Larger different-
input isolation remains open rather than claiming complete reentrancy.

# Standard embedded TDDFT initialization MPI ownership

Port the proven fast-PR child-owner tracker into the bundled standard
communication module. Keep the standard one-native-lifetime/process constraint.
Capture children from standard embedded TDDFT initialization through its initial
current calculation; close request caches and XC, discard pending restore payloads,
then release that state's children before deallocation. Host/attachment resources
remain owned by their caller. No changes to numerical update order or checkpoint
format, no claim of general reentrancy, no process/MPI failure recovery.

First show old runtime lacks the resource capability. Raw host-MPI probe records
children before/during/after close while attachment stays live; exercise pending
restore candidate close too, verify host collective and old repeat-init refusal.
Compare old/new traces with PZ/TB-mBJ/Libxc-PBE under orbital/r-space MPI splits.
Verify direct restart/protected ordinary mixed FDTD/LD/SBE/TDDFT and full-k HSE
mixture. Fresh explicit external source/build, manifests/shared-kernel byte parity,
source/wheel/install tests. Paired guides/AI/checklist. No Actions.

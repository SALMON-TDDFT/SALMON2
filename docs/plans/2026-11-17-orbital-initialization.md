# Orbital layout initialization inputs

Continue approved initial-state ownership extraction with an internal value-only
s_orbital_initial_config and s_orbital_initial_layout in initialization_sub.
Capture calculation mode,periodicity,spin/SOC,state/electron counts,temperature and
k-mesh/file sentinel before downstream initialization. Explicit resolve selects
spin count,orbital count and real/complex storage without global access; output
is assigned only on successful resolution. Orthogonality comes from resolved
stencil. HSE uses its owned occupation input/effective kmesh and captures the same
legacy hybrid predicate; standard retains its existing nonhybrid predicate.
Keep legacy field/file side effects and occupation evaluation order unchanged.
No full GS reentrancy, electron-capacity certification or stable ABI claim.

Native table tests cover GS Gamma/non-Gamma/nonorthogonal/isolated/SOC/hybrid and
RT negative/finite-temperature/even/odd/spin-resolved counts; capture followed by
poisoned globals,alternating outputs,invalid-result preservation and file sentinel.
Expected missing-symbol compile first. Standard/HSE/RT-profile source preparation,
new external explicit builds,source/wheel tests,GS MPI2 and PBE/HSE mixed old/new
fault/retry/direct-restore regression. Update JP/EN/AI/checklist,retain evidence,
no Actions,OMP/BLAS1,warning-only memory budget.

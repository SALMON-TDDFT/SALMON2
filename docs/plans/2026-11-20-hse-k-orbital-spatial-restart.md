# Fast-PR HSE multi-k orbital/spatial bridge and restart

Extend hse_kspace compute3 to orbital groups, preserving full/reduced-k and
legacy routes. Gather accepted rows/columns over the r/o group for refresh,
run the existing fast k-distributed exchange, then scatter local rows/columns.
Exchange scratch is replicated across r/o ranks, not a scaling claim. Apply
orbital-distributed dense ACE per local k slice using existing orbital overlap
collectives; no new exchange approximation. Integer occupied spin-pair CPU
profile, same owned orthogonal cell; reduced-k retains existing1/2 k groups.

Mixed repartition3 already maps five-dimensional ownership and HSE row/column/k
caches. Extend admission without replaying SCF/density/exchange. Verify full and
reduced-k rank changes2↔8 and uneven orbital blocks, exact cut/zero bootstrap,
continuation1e-9/1e-10, late donor and protected different-input retry, prior
profiles, source/object audits, wheel/sample, JP/EN/AI/list. Do not dispatch CI.
Packed/tiled SCF storage remains separate; built-in SCF still generates dense ACE.

Accepted: full2→4/2→8/8→4,orbital2→4,uneven6→4,reduced2→8/8→4 and
pre-PBE2→8. Physical cut interiors exact, bootstrap0, continuation within
unchanged1e-9/1e-10. Protected full/reduced different-input owner failure/continue/
retry exact; previous full/reduced spatial/Gamma/ordinary results exact. Source/
installed637 tests each636pass/1skip; installed reduced8→4 and public sample8→4
pass. Source/object audits and source/wheel/docs/header/native-byte parity
verified. Campaign source/wheel GS/restart export explicitly authorized. Evidence:
calc2026/salmonpy-hse-korspace-record-20261120. No CI dispatched.

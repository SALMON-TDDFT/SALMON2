# Approved Si middle five-intensity design

Human approved option1: one middle job, 10 concurrent RT cases, 5472 nodes each
(54720 nodes total), 4 MPI/node, OMP12, 24h. Each of five intensities
1e10, 10**10.5, 1e11, 10**11.5, 1e12 W/cm2 has pump-only and pump-kick.
Photon energy1.55eV, Nr16/Nk16, existing 8cycle z pulse,1fs wait,20fs response,
dt.03au/K16/nteach8/history3; projection off; validated balanced46caf executable.
Shared GS computed once on one case partition, verified before 10 RT launches.
Explicit disjoint vcoordfiles; all-rank layout/memory/finite/SVD/phase/checkpoint
checks; no repeated source/object audits for unchanged binary. Preserve evidence.
Prepare/verify new job completely before cancelling old QUE/unstarted52968000.
Cancel once only if still unstarted; submit new middle once with receipt/locks.
If old already starts, do not cancel/modify automatically. Other jobs unchanged.
No auto retry/resubmit; final incomplete checkpoint EC86 differs from complete0.
No GS/checkpoint/wfn payload downloaded. NonGS records canonicalcalc2026, same96.

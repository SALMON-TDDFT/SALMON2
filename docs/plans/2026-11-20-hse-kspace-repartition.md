# Fast-PR HSE multi-k / spatial SCF bridge and restart

Keep the existing fast full-k exchange algorithm. Store wavefunctions and dense
ACE/cache rows locally. Explicit hse_kspace profile gathers full-grid scratch
within each k group's spatial communicator for exchange refresh, runs the fast
k-distributed kernel and returns row-local accepted caches/factors. Exchange
FFTs/work are replicated over spatial ranks, not a distributed-FFT speedup.
Initial profile: CPU full uniform mesh, no symmetry/orbital split, scalar full
occupied spin pairs and orthogonal cell. Old profiles keep their route.

Mixed repartition3 maps rows and k slices from authenticated rank files, retains
SCF/mixing and HSE state, rebinds communicator resources and publishes collectively.
No exchange/density/SCF reconstruction of the saved state. Verify changed k and
spatial rank counts in both directions, exact cut interiors, bootstrap0 and five
iterations vs source-initialized continuous8 within1e-9/1e-10. Verify late donor
failure preservation, independent/protected owners and previous full-k results;
source/object audits, Python suite, source/wheel/public sample, docs/JSON/list.
No routine CI. Results in calc2026; obtain campaign-specific GS export permission
if automatic review requires it.

Accepted with final source-v5:2↔4,2→1,4→8,changed-axis4→4 and pre-PBE2→4;
physical wave/density cut interiors exact, bootstrap0, continuation within
unchanged1e-9/1e-10. Protected different-input mixed owners and failure/continue/
retry exact; previous k-HSE/Gamma-mixed/ordinary results exact. Source/installed
637 tests each636pass/1skip, wheel4→2 and fixed8k public sample4→2 pass.
Source/object audits and source/wheel/docs/header/native-byte parity verified.
Earlier failed builds and outputs retained. Main initial failure was the input
checker's automatic Gamma Wannier selection before native refresh, now avoided
only for attached full-k spatial SCF. Non-MPI guard retains CLI compatibility.
Evidence: calc2026/salmonpy-hse-kspace-record-20261120. No CI dispatched.

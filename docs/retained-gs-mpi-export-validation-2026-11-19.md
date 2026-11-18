# Retained GS collective portable export — 2026-11-19

Regular-mesh DFT-owned retained GS now exports through one subgroup-collective
API. Native lower bounds and physical density bounds determine tiles; no rank
numbering is assumed. Replicated occupations/eigenvalues and owner/input/options
must match. Exact global k/band/grid coverage rejects overlap/gaps before transfer.
The root receives one tile per collective, assembles the wavefunction and writes
existing GS/2; workers return None. Root's global-wave allocation and publication
failures propagate collectively; native owners are never mutated by this helper.
Incomplete new archives are removed by the writer; existing paths are preserved.

No blanket memory cutoff. Root still needs full global wavefunctions, local
copies, one received tile, density and writer buffers. Rank-ordered communication
and pairwise tile checks need larger-scale calibration. Export alone does not
implement native GS-handle/full mixing-history checkpoints or changed-topology
restoration. Arbitrary MPI/process failure is outside this contract.

Actual MPI4 synthetic tests cover five k/orbital/spatial/mixed decompositions,
including uneven band ranges. All saved file hashes equal the serial reference
when synthetic provenance is fixed. Eight injected failure cases cover bounds,
replicated arrays, overlap, incomplete coverage, native metadata, selected options,
root allocation and root write; every rank receives the failure, source arrays
remain unchanged and a clean retry matches serial. Existing output is preserved.
Fourteen cases pass both from source and from the installed final wheel.

Source623 and installed623 unit methods pass with one skip each. Both real and
complex synthetic Gamma archives remain byte-identical to b66b97fd. The sample,
JP/EN/AI guides, package manifests and bundled MPI runner are synchronized.
Wheel built from sdist offline; bundled-kernel source and all application-object
prohibition audits pass, as do native wheel hashes. Only the existing LD/SBE
wheel kernels were rebuilt; no electronic SALMON runtime was rebuilt or used.
No new SCF, no actual native GS export, and no real GS arrays were collected.
These tests establish MPI/Python assembly, not native GS/SCF/TDDFT physical parity.
Measurement platform: Apple Silicon/GNU/OpenMPI/Python3.14; other hosts unverified.

Selected logs/hash manifests/wheel only (no GS arrays):
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-retained-gs-mpi-export-20261119
Local writes verified; cloud synchronization unverified. No CI/push this slice.
Preceding b66b97fd remains pushed; current multi-k/MPI work is uncommitted.

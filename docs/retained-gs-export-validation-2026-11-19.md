# Retained GS portable export — 2026-11-19

Preceding mixed-HSE assembly was committed as `2340787d` with `[skip ci]`.
This subsequent Python-only slice adds `NativeGSState.export_archive` for
serial Gamma DFT-owned retained states, using existing GS/2 save/open machinery.
No native source, numerical update order or ABI was changed.

Synthetic tests use a noncubic8×10×12 physical grid, asymmetric halos, two bands,
real/complex orbital buffers, occupations/eigenvalues and an unconverged flag.
Wavefunction values survive portable save/open exactly, with correct axes/halo
removal; reconstructed density agrees within1e-15 relative tolerance (complex
absolute value versus real/imaginary squared operations differ in rounding).
Source buffer mutation and owner close preserve archived values. Readback is
immutable; existing output and invalid layout rejection preserve source storage.
Unconverged SCF warns rather than blocks. Material ID, native identity/runtime,
convergence and local channels are validated in optional GS/2 provenance.

Source suite:622 methods, one skip, all remaining pass. Focused native-GS/archive
suite:17 methods pass. The optional sample's `--help` passed; its SCF was not run.
JP/EN guides, AI search metadata/embedded text and checklist are synchronized;
manual sync and git whitespace checks pass. No new wheel validation is claimed.

No actual retained GS was exported, no new SCF was performed, and no real GS
restart/orbital/density data was collected. Real runtime export/TDDFT parity is
not established by these synthetic tests. Portable array readback is not native
GS handle restoration or a full SCF checkpoint: mixing history/native resources/
forces/mixed density are absent. Multi-k/MPI retained export, imported-owner
physical identity and changed-topology restore remain future work.

Selected evidence only:
`/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-retained-gs-export-20261119`.
Local files verified; cloud synchronization unverified. No native build, push or
GitHub Actions was performed. Existing HSE regression evidence remains separate.

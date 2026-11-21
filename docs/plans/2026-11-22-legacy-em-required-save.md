# Legacy EM required-save recovery

Implemented checked per-file staging, mandatory sticky status, and one save-phase MPI agreement. No backward-compatibility fallback is maintained for unpublished Fortran helpers. Preserve failed staging evidence; failed multi-file sets are not valid archives.

The old implementation overwrote the existing payload before a sidecar directory error. The negative fixture reproduced this; the new implementation retains that payload and returns status. Text and real/complex arrays of dimensions1–5 match the old normal bytes (40 files at MPI1 and66 at MPI2 per distribution). Source/installed ownership/remapping fixtures and679 Python tests pass with one optional skip. Real LD+FDTD normal outputs and four native source/wheel MPI2 scalar/tensor controls match prior bytes. Batch failures return nonzero without hanging. No GS is collected.

Evidence: calc2026/salmonpy-legacy-em-save-recovery-20261122/summary.json. Remaining R3 scope includes selected startup/helper and restart-reader I/O and auxiliary allocation handling. This does not certify arbitrary hardware/filesystem failures, Fujitsu, multi-node execution, fsync durability, or atomic multi-file archives.

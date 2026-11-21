# Required batch EM restart-directory preparation

Prior work was committed/pushed as a7073552. The next R3 item exposed an undefined optional status when filesystem.create_directory received an existing directory. An existing test fixture with sentinel status37 reproduced the bug before this fix. The routine now explicitly returns zero for existing directories. Batch FDTD passes creation status and replaces its creation barrier with one owning-group failure agreement; root-only STOP is avoided. No compatibility fallback, extra ordinary-step communication or new numerical calculation is added.

GNU macOS MPI1/2: fresh/existing/symlink-to-directory LD+FDTD outputs and raw restart files match previous successful controls. Conflicting files and symlinks to files remain unchanged and batch ranks return nonzero without timeout. Source/wheel ownership and writer tests,679 Python tests(each one optional skip),and four native scalar/tensor MPI2 direct-save/restore regressions pass. This is not a writability preflight or a guarantee against OS/MPI failures. No GS data is collected.

Evidence: calc2026/salmonpy-restart-directory-recovery-20261122/summary.json. Follow-up remains uncommitted. Next R3 work remains selected restart-reader/helper I/O and auxiliary allocation failures; R4 native ABI stability follows before group3.

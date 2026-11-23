# Fast-PR HSE k-only host-MPI direct continuation — 2026-11-16

## Result and scope

Implemented experimental `NativeFastHSEBackend`: host-owned MPI communicator, one material point/group, one native lifetime/process on the main thread. Only the pinned fast-PR HSE engine is supported. Validated k-only groups of one/two ranks; orbital and real-space decomposition remain rejected. No MPI init/finalize, compilation, download or native loading occurs on import. Native calls temporarily change cwd.

Local environment: Apple Silicon, macOS26, Python3.14, GNU15/OpenMPI, OpenBLAS, FFTW, Libxc5. Each MPI rank used one OpenMP/BLAS thread. At most four active compute ranks. This is not a4GiB memory-ceiling measurement, cross-platform validation or scaling benchmark.

## Coupled numerical acceptance

Installed the wheel built from the source distribution. Used native E/H FDTD40³, two field-local LD materials, one SBE point, one HSE point with internal k ranks2: field1+SBE1+HSE2=4 outer ranks. Microscopic HSE input: real orthonormal synthetic orbitals,8³ grid,4×1×1 k mesh,two bands,four electrons. This is not a converged material GS. No initial-current subtraction, no internal driving pulse; transverse external A/c only.

| Fast propagator | Total / split / resumed calls | Full A/c,E,j histories | Maximum current difference from independent serial fast-PR response |
|---|---|---|---|
| Taylor4, learned ACE |256 /127 /129|bit-identical|4.063047641389428e-16 au|
| PT-CN, learned factors |256 /127 /129|bit-identical|1.7049079162334557e-16 au|

The point-ID and serial-reference arrays also matched byte-for-byte between continuous and resumed outputs. dt=.001 au. Fresh-job LOAD restores accepted step/current/A and dynamic numerical payloads directly; the communication counters on every rank show129 new material/dispatch/collection calls. Initialization reconstructs static exchange/FFT resources but performs no elapsed propagation. Independent serial replay runs afterward solely for validation. Update formulas and coupling order are unchanged.

Saved state includes field/PML,local LD P/J/Q,old fields,SBE density matrix,HSE orbital buffers/density/potentials/RT and ACE/learned/PT-CN histories,integrated A and lagged current/clock. Group size/rank,k ranges and array layouts are validated. Rank-local opaque MPI/FFT/XC handles are freshly built, never serialized. Archives were approximately54.77MB, including JSON histories. Same binary,input,bundle,topology and identical shared files are mandatory.

## Communication measurements

Field-rank accumulated wall seconds for256 continuous steps:

| Propagator | A distribution | Material phase including waiting | j collection |
|---|---:|---:|---:|
| Taylor4 |0.000722|1.915895|0.009227|
| PT-CN |0.000815|2.962447|0.010538|

These phase timings include synchronization/waiting and are not isolated network latency or strong-scaling measurements. Logical point traffic is48 bytes of A and48 bytes of j per step; MPI metadata and internal HSE collectives are not included in that byte estimate. Resource measurement errors are separate from numerical failures.

## Additional checks and corrected defect

- Native HSE1 versus HSE2 coupled256-step runs also passed: maximum absolute A/c,E,j differences Taylor `(1.71e-17,1.33e-16,4.02e-16)` and PT-CN `(1.97e-17,1.53e-16,4.62e-16)` au, below the explicit comparison tolerances rtol2e-10/atol2e-12. MPI reductions need not be bit-identical across decompositions.
- Native two-rank HSE direct save127→256: restored current/A and continuation current/energy bit-identical. `tools/test_fast_hse_mpi_checkpoint.py` exercises fresh jobs and per-rank shards.
- Truncating only rank1's shard and swapping rank0/1 shards return collective failure without hanging or advancing Python step/A/current. Tests bypass manifest authentication deliberately to exercise native validation. Coupled archives additionally verify all hashes before native LOAD.
- Existing serial fast-PR direct checkpoint regression:16 continuation cases and24 rejection cases passed, including byte equality of internal traces. Footer/trailing-byte validation now precedes any numerical commit, so these early rejections retain accepted step0.
- Python source and installed-wheel suites:470 tests,1 skip, both successful. Source-distribution→wheel build succeeded. Shared LD/SBE numerical kernels were unchanged.
- Found an allocation-dependent false restore rejection: the exchange kernel's distance-index table only initializes each axis's active rectangle, while checkpoint comparison included unused padding. Capture now zeroes padding and copies initialized elements only; comparison uses only active rectangles. The actual solver table and numerical operations are unchanged.

All shards are validated before HSE numerical commit. A late second-read/allocation/I/O failure can still invalidate contexts after partial commit; abort the job. This does not provide cross-group rollback, changed-topology restore, unrestricted reentrancy, orbital/rgrid MPI, GPU or generic backend persistence.

## Reproduction and evidence

Use the internal-k-MPI section of `salmonpy/COUPLED.md` / `COUPLED.en.md`. Prepare fixtures explicitly with `tools/prepare_fast_hse_coupled_fixture.py PSEUDO NEW_DIR --kpoints 4 [--propagator ptcn]`; build with `tools/build_fast_hse_mpi.py`. Run continuous,save127,resume jobs using the same256-step configuration and `--hse-ranks 2`. The serial fast-PR worker is the independent reference, never a conventional HSE engine.

Local evidence root: `/Users/otobetoshihito/.codex/outputs/salmonpy-state-inputs-20261116`.

- `hse-inner-mpi-build-v2/build.json`: native build/source/adapter receipt. Library SHA256 `d6bddfc02f4f351491ca647e065685d160627017fb0b61addb6be10f56e23947`.
- `hse-inner-final-comparison.json`: byte comparisons,archive inventories/hashes and communication records.
- `hse-inner-final-{taylor,ptcn}-{continuous,save,resume}` and corresponding logs.
- `hse-inner-negative-{save,resume,truncate,swap}.log`.
- `hse-inner-serial-regression-v2/result.json`.
- `hse-inner-source-tests.log`, `hse-inner-installed-tests.log`, source/wheel build logs.

No remote Actions were triggered for this validation.

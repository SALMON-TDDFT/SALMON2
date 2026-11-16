# Fast-PR HSE in MPI / FDTD with direct persistence

User goal: connect fast-PR HSE to MPI/FDTD, then restore full internal coupled state without elapsed propagation replay, matching continuous execution. Only fast-PR HSE is an execution engine.

Design implemented: preserve the existing serial fast-PR Fortran context in a dedicated serial process. An explicit worker builder links audited fast objects and SALMONpy adapters into a new output; import never builds/downloads/initializes MPI. Each single-rank material group owns one worker, exchanges endpoint A/c and charge j, and can directly SAVE/LOAD numerical state. FDTD retains its native MPI communicator and local LD update order. This is outer MPI coupling, not HSE internal MPI.

Worker protocol saves stdin separately from the legacy Namelist reader with POSIX /dev/stdin, reads inputfile, emits prefixed flushed replies, and supports STEP/SAVE/LOAD/CLOSE. Provider pins immutable bundle identity; backend pins worker/origin hashes, stages independent state copies and keeps host cwd unchanged. Reject bad profiles, paths, hashes, concurrency, missing input and timeout. Return actual initial current and no background subtraction.

The coordinator admits explicit direct-state callbacks alongside existing native FDTD/TDDFT and SBE paths, authenticates shared same-topology archive, restores native worker step/current and Python A/adapter clocks. One worker per rank is enforced. No elapsed time evolution is replayed; fresh static initialization and independent reference validation remain separate.

Acceptance completed: same3k time-reversal real-orbital synthetic fixture; learned ACE Taylor and learned PT-CN; native field1/2 + SBE1 + serialHSE1 + two local LD materials. Four configurations each continuous256 versus save127/new-process129: full A/E/j bit-identical, worker command logs128..256, independent response replay exact. Field1/2 trajectories also identical. Source/installed459 tests,1 skip each; packaging and guides synchronized. Full evidence: docs/fast-hse-mpi-coupled-checkpoint-validation-2026-11-16.md.

Remaining: HSE inner-MPI/native-host profile, stable generic checkpoint ABI, cross-topology restore, broad HPC/platform/performance/physics validation and memory estimates/caps. No conventional HSE or implicit fallback is added.

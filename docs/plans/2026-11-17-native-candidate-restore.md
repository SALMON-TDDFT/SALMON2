# Native candidate restore

Extend the existing same-binary numerical checkpoint ABI with explicit modes:
0 save, 1 legacy immediate load, 2 stage detached candidate, 3 publish prepared
candidate, 4 discard candidate. An exported checkpoint_protocol query identifies
this revision without attempting a load against an old library.

Candidates own numerical allocatables only. Communicators, Libxc handles, FFT
plans, output units and caches remain owned by the initialized runtime. HSE stages
its validated numerical exchange/history payload once; publication does not
reopen the checkpoint. Staging an advanced original is supported; the legacy
fresh-load restriction remains for mode1. There is one pending candidate per
native lifetime. Hosts pause all evolution during this transaction.

All ranks stage, then agree success over the host communicator. Any handled
failure discards every attempted candidate without touching numerical originals.
Publication runs only after the global vote and uses prepared memory, with no
archive I/O or user callbacks. Allocation failure/process death/failed MPI is
outside the handled-error contract. Required immutable kernels must already
exist; staging must not initialize exchange into the original context.

Verify real field+standard-TDDFT+two-rank-fast-HSE groups: save at a split, advance
originals, fail a later shard read, compare before/after native state bytes and
continue originals. Retry and compare future fields/currents against continuous
execution. Retain existing checkpoint formats and regression checks.

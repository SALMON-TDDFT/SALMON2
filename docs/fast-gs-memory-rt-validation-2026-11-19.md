# Retained GS to fast-PR HSE memory initialization — 2026-11-19

Scope: local Apple Silicon, GNU15.2/OpenMPI/Libxc/OpenBLAS, scalar periodic fast-PR
HSE06, k-only MPI1/2, OpenMP/BLAS1. Synthetic8³/2-band fixtures test implementation
parity, not physical convergence. No fresh GS/SCF calculation or GS dataset
collection. No GitHub Actions runs were launched.

## Implementation

`NativeGSState.tddft_model` dispatches on the runtime profile inside the existing
GS/TDDFT adapter. The fast bridge copies a value-owned seed into the model registry;
initial and subsequent independent RT owners use it instead of restart_rt. Static
RT/HSE initialization proceeds normally. The model's seed outlives the source GS.
No GS restart files are copied/written/read for this handoff and no accepted steps
are replayed. This uses memory for retained copies; it is not zero-copy.

Subsequent GS imports park/restore active RT input/symmetry workspaces. Both source
profiles use the standard timer's existing initialize-once guard. The first mixed
GS/RT test exposed duplicate timer allocation in the previous fast path; failed
build/runtime logs are retained as evidence. A later no-MPI compilation exposed
an MPI-only public declaration, now guarded by USE_MPI.

New builds query the compiled memory protocol and declare gs_abi.memory_rt=1;
old runtimes are not rebuilt implicitly. Existing standard GS/RT calls remain
unchanged. Fast session.run, force/eigenvalue computation and archive export
remain pending, as does ordinary-XC RT migration into the fast source profile.

## Controls

The existing `tools/test_memory_gs_mpi.py --fast-hse` splits2*size ranks into a
memory-seed group and an existing-file group. It compares16 steps of electric
charge-current and total energy with rtol3e-13/atol3e-15. Each group drives two
independent states with opposite A. Full-k drives all three components with a
rotating xy part; reduced C2z uses symmetry-preserving z. At step8 a direct checkpoint is loaded into a
fresh owner; its subsequent currents/energy must match continuous execution
exactly. Source GS closes and its staged source bundle is removed before further
state creation. Interleaved GS import, unused bootstrap cleanup, altered first
state dt/capacity and live children after model close are checked. RT work contains
no wfn.bin and its configured GS read directory does not exist.

Full/reduced learned PT-CN and learned Taylor are covered. The C2z-reduced fixture
uses z drive; full-k has no universal z-only constraint. The zero internal impulse
is explicitly assigned a harmless z direction to satisfy the existing input parser;
it never projects or changes externally supplied A.

`tools/test_hse_registry_mpi.py` compares the new normal file-initialized registry
with the preceding fixed runtime, internalMPI2/overallMPI6. Currents, total energy,
exchange history, fresh state and direct-resume results remain exact. These controls
use independent old binaries and do not infer normal parity from a self-comparison.

Existing retained_gs sample import mode optionally accepts --current-steps 16,
showing source-GS close followed by two independent current calculations. It is
an execution demonstration; small dt/pulse/grid choices do not certify optics.

Records and distribution checks are retained under the local canonical results
store at calc2026/salmonpy-fast-gs-memory-rt-20261119. Only logs, receipts, checksums,
results summaries and distribution/binary artifacts are collected there; GS,
wavefunction and direct-state checkpoint payloads are excluded. Cloud sync and
other OS/compilers are unverified.


## Final local checks

- Source and installed suites:628 methods each, one external-runtime skip each.
- Four final memory/file controls: max current difference0 and max energy
  difference0; direct resume exact, SCF0/replay0. MPI2/4 overall as above.
- Final old-runtime registry comparison and GS-import independence: MPI6 pass.
- Same prepared source builds fast-HSE MPI and HSE-disabled serial. Mandatory
  source audits and all application-object audits pass, including host adapters.
- sdist→wheel succeeds. Installed source preparation reproduces all258 src hashes.
  Installed retained_gs sample and MPI4 memory/file/restart controls pass.
- Installed FDTD/LD/SBE/two-input-HSE mixed128 is exact against the preceding
  mixed baseline. Direct save8→resume120 matches continuous128 without replay.
- JP/EN guide snapshots and AI JSON synchronized. No test-method/module/sample
  count growth. Compiler/platform and synthetic-physics limits above still apply.

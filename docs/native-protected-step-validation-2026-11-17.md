# Native protected mixed-step validation — 2026-11-17

Scope: the explicit whole-step transaction protocol, separate from individual
native numerical recovery and disk checkpoint protocols. Implementation and usage
are in [NATIVE_STEPPING.en.md](../salmonpy/NATIVE_STEPPING.en.md) and its Japanese
pair. Molecular work remains paused.

## Local campaign

Apple Silicon macOS, host capacity 64GB / 18 cores; at most seven MPI ranks,
OMP/BLAS one thread per rank, sequential runs, 48GiB advisory budget. Budgets only
warn. This campaign does not establish Linux/Intel Mac or physical convergence.
Production libraries were explicitly built into new external directories from
manifest-verified sources; shared LD/SBE kernels remain unchanged.

Ten 32-step cases passed with the final production builds:

- Native E/H field2 + local two-material LD + SBE1 + electronic2: PZ,
  automatic/fixed TB-mBJ, Libxc-PZ, Libxc-PBE, fast-PR HSE Taylor and PT-CN.
- Field2 + local LD + SBE1 + Libxc-PBE2 + fast-PR HSE2 (seven ranks),
  Taylor and PT-CN. HSE and PBE use separate initial-state bundles and identities.
- Multiscale1D field1 + two independent Libxc-PBE2 groups (five ranks).

At attempts 1, 13 and 32, the harness rejects after actual field, SBE and last-rank
electronic advances and exercises real electronic overflow; native field overflow
is also exercised at attempt 13. It compares accepted numerical checkpoint bytes,
Python array identities/values/clocks/flags and every ordinary `.data` stream before
and after rollback. No propagation kernels are mocked. All valid explicit retries
and final A/E/j arrays match uninterrupted execution exactly. All ordinary `.data`
streams also match byte-for-byte. The 1D grid/history comparisons are exact.
Five cases additionally match unchanged pre-transaction standard/fast-PR libraries
(PBE, HSE Taylor/PT-CN and both seven-rank mixtures), including ordinary output bytes.

Nested native transactions, inconsistent MPI modes and checkpoint staging during
an active step backup reject collectively. Three-rank negatives verify old-library
rejection, an actual HSE output subdirectory rejected before advance, and a failed
native rank remaining failed across its group. The latter deliberately does not
claim original-state preservation. A seven-rank protected checkpoint-restore
failure/continuation/retry regression also passes with protected stepping enabled;
checkpoint publication uses freshly constructed targets for replaced owners.

## Verification and reproduction

- Source and independent installed-wheel suites: 545 tests, one skip each.
  Installed seven-rank PBE/HSE and two-point 1D fault/retry campaigns pass; arrays
  and all ordinary output streams match uninterrupted execution exactly.
- Source distribution to wheel and native SHA256 validation pass.
- `tools/test_native_protected_step_mpi.py`: ordinary mixed sample arguments,
  optional source-only `--kernel`; `--multiscale` selects the 1D sample.
- `tools/test_native_step_admission_mpi.py`: old-protocol, output-rejection and
  unsafe-lifetime negatives. Use launcher timeouts and matching audited libraries.
- `tools/test_mixed_protected_restore_mpi.py`: accepts `--protect-steps`.
- Paired human guides and AI JSON are synchronized. No remote Actions were used.

Persistent local evidence is retained under the Codex outputs directory as
`salmonpy-native-protected-step-20261117`: production build receipts/libraries,
source/adapters, execution recipes, test logs, comparisons and checkpoint/file
hash inventories. These local libraries are not a published portable runtime.

## Limits

The implementation restores numerical memory directly; it does not replay elapsed
accepted propagation. Retrying recomputes the rejected current step. FFT plans,
communicators and other resource owners remain live rather than being cloned.

Output capture first flushes accepted buffered records. Rollback closes audited
RT units, truncates private regular single-link root `.data` files, and reopens the
same units. Tentative appends may be externally visible: no concurrent readers or
writers, moving work directories or durable crash-atomic publication is promised.
Unsupported output profiles reject rather than silently disabling user requests.
I/O failure, process/MPI failure, allocation exhaustion, initialization/uncontrolled
STOP, unreturned native calls, general reentrancy and stable public C ABI remain
outside this contract. CLI errors still terminate when the host leaves them
unhandled; no automatic retry or timestep adjustment is introduced.

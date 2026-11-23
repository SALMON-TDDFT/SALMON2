# RT context and output ownership validation — 2026-11-18

This slice moves fast-PR protected-step buffers and pending direct restore payloads
from process globals to the context that owns the numerical state. It does not
lift the public Python/C one-native-lifetime limit. All tests are local; no
GitHub Actions run is required. Synthetic HSE fixtures certify state isolation,
not physical convergence. Ordinary results are compared with the preceding build.

## Reproduction

Explicitly prepare/build a new fast-PR runtime using the current adapters, then:

```sh
python tools/test_fast_hse_context_owners_mpi.py BUILD FULL_K_FIXTURE NEW_OUTPUT
```

The tool verifies the library receipt/current adapter, explicitly compiles native
probes, and runs two independent material groups at internal MPI1/2 and OpenMP1/2.
The same immutable model within each group has A(dt=.001,capacity32) and
B(dt=.002,capacity34); different groups use different HSE omega. A uses z; B
uses a rotating xy field. Independent controls set their legacy dt/capacity too,
so agreement does not merely reproduce the same unintended global configuration.

Three modes (Taylor, learned Taylor, learned PT-CN), four execution layouts and
eight sequences cover isolated A/B, ordinary interleaving, simultaneous protection,
staged restore/commit/cancel, close A while B's restore is pending, opening B
while A's output is closed, failed reconnect followed by close or retry. Each B
continues after A closes. All current/energy traces and final rank-local numerical
probes agree byte-for-byte with isolated execution. Ordinary and late-open cases
also preserve complete RT/energy output bytes. The staged-restore sequence does
not truncate text logs by itself: host journaling performs output rollback.

## Failures reproduced before correction

- Two live contexts can interleave normally, but the second simultaneous capture
  fails status14 because the first owns the global step reservation.
- Closing A's outputs lets B reuse A's former positive unit numbers. Reopening A
  on those numbers silently reconnects B's streams. The small output probe fails
  STOP5 before correction.
- Discarding transaction metadata after a reconnect failure lets close(A) close
  B's reused unit. Numeric B traces can still match while output is lost: one
  reproduced B RT stream had472bytes instead of12979bytes. This is why actual
  output bytes are now compared, alongside numerical arrays.

Reopen now selects an available positive unit and updates the owning output
object. Failure retains numerical/output metadata collectively. Close uses the
transaction's closed flags, not stale unit numbers. The byte-identical standard
helper and embedded TDDFT call sites are kept synchronized; ordinary standard
close also closes its actual dynamically allocated energy unit.

## Remaining boundary

Only non-symmetric same-model internal contexts sharing one communicator are
exercised here. Global initialization/symmetry/input, timer and MPI attachment
ownership and the synchronous serializer array bound remain separate tasks.
No stable multi-handle C API, thread safety, repartitioned checkpoint or general
heterogeneous same-process execution is promised. Checkpoint encoding, numerical
kernels, and time-update formulas/order are unchanged.

## Recorded results

Apple Silicon/macOS, explicit GNU/OpenMPI/Libxc/OpenBLAS builds:

- 12 configuration cases,96 MPI invocations and192 owner comparisons passed.
  A further three-case tool smoke includes the direct unit-reuse helper probe.
- Standard Libxc-PBE internal-MPI2 output-close probe: old ordinary close fails
  STOP7 (energy stream left connected); final ordinary/protected close passes,
  including preservation of the host output opened during rollback.
- Four existing native protected mixed-step regressions pass against preceding
  continuous references: PBE, Taylor HSE, learned PT-CN HSE and combined PBE/HSE,
  with5/7 outer ranks. Final normal trajectories retain previous results.
- Reduced-k learned Taylor/PT-CN deep numeric failure atstep13: original state
  and output are retained, then retry matches an independent continuous group.
- Source and installed-wheel suites: 590 tests each, one skipped, no failures.
  Source-distribution-to-wheel and installed native hashes pass. JP/EN/AI guide
  synchronization and source inventories pass.

HSE library SHA256:
`9291c7caa65172a0c3100c76d4b5bebfd74acada74ff10bcdb5edba4afae66bb`.
Standard library SHA256:
`ed6c3dee8f8229bff711ba903e251e693cfa08b1f0b89d9f03cf0067a46f807d`.
Evidence bundle: local `salmonpy-rt-context-owners-20261118` output directory;
receipts retain compiler/link commands, source hashes and original build paths.

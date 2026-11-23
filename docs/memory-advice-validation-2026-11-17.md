# Advisory memory budgets — 2026-11-17

## User decision and implementation

The user replaced budget-excess stopping with warnings. No memory admission stop,
OS limit, scheduler limit or automatic reduction of physics was implemented.
MemoryEstimate/MemoryReport/memory_advice provide a common byte-count interface
for built-in and custom backend estimates. None/unknown are distinct from zero;
the communicator's rank estimates are summed, not treated as per-rank budgets.

The ordinary mixed sample's optional --memory-budget-gib records pre-initialization
heuristics, initialized/segment-end process peak RSS, and protected-restore advice.
Restore includes extra4*native payload sizes+32*JSON sizes+snapshot allowances,
above observed rank process peak RSS, plus128MiB/rank reserve. The multipliers are
heuristics, not a proved heap bound. Native field includes all LD/PML state;
standard TDDFT and fast-PR HSE payloads and independent SBE snapshots are counted.
Static initialization workspace gaps and excluded child processes remain unknown.

Process probing is lazy standard-library Linux/macOS peak RSS. No MPI init,
compile, native state query/allocation probe, download or psutil dependency on
import. Probe failures return unavailable data; warning-filter exceptions and
report-output failures are separate diagnostics and do not invalidate computation.

## Local verification

All jobs run sequentially with4 MPI ranks, OMP_NUM_THREADS=1 and
OPENBLAS_NUM_THREADS=1. Existing audited native v4 libraries were reused unchanged.
No remote Actions or new native runtime build was used.

- Source and installed-wheel suites each pass509 tests, one optional skip.
  Tests cover coexistence arithmetic, unknowns, warning exceptions, payload-size
  estimation without reading arrays, report write failure and child uncertainty.
- Real MPI4 API test: over-budget warning on root only, aggregate4GiB comparison,
  different rank budgets reported unknown, unavailable probes, warning-filter
  error capture; host communication remains usable.
- Real native field+two LD materials+SBE+orbital-MPI2 Libxc-PBE TDDFT with an
  aggregate4GiB budget: original byte preservation on late restore failure,
  old-state continuation and exact retry65–67 all pass. Installed restore
  estimates are about1.373–1.509GiB including candidate/read copies and reserve;
  candidate allowances are about0.306–0.312GiB across the job. Initialized advice
  is0.911GiB. Segment-end process high-water RSS sums to about0.709GiB; adding
  the0.5GiB reserve yields the1.209GiB segment-end estimate.
- The pre-initialization known-component estimate is3.284GiB, but initialization
  workspaces remain unknown. It is explicitly not classified as assured compliance.
- Identical-sized native work with only the advisory threshold lowered to1MiB
  and a rank3 RSS probe forced to fail still completes, preserves originals and
  retries exactly. The failed probe stays JSON null/unknown; excess only warns.
  This does not claim the actual job fits1MiB. Its three continuation frames
  equal the4GiB case exactly (maximum error0).
- Installed native fast-PR PT-CN HSE k2+SBE+two field-local LD materials also passes
  with4GiB/4 cores: restore estimates1.202–1.321GiB, strict128→131 retry parity,
  and end-of-segment advice1.109GiB including reserve.
- A fresh ordinary installed CLI job resumes the4GiB checkpoint64→68 with memory
  advice enabled; all fields/clocks/currents at65–67 equal continuous exactly
  (maximum error0). The scientific archive does not bind the advisory budget.
- Source-distribution-to-wheel build/install, public API inventory and JP/EN/AI
  synchronization pass. Independent review found excluded child memory missing
  from restore reports; that was fixed, unit-tested and rechecked.

Persistent logs, JSON advice, scientific references and final wheel:
`/Users/otobetoshihito/.codex/outputs/salmonpy-memory-advice-20261117`.

## Limits

RSS is a per-process high-water mark, not current use or simultaneous whole-job
peak. Summing rank peaks may double-count shared pages and different-time peaks;
MPI launchers/child workers are excluded from the probe. The sample reports a
reserve and child uncertainty, not a measured full-machine limit. Only the local
Apple Silicon environment was executed here. Windows/other unsupported probes
return unavailable; callers may provide their own estimates.

Future native workspaces, larger systems and arbitrary external backends need
their own calibration/declared estimates. This does not certify capacity,
scientific convergence, Windows native MPI support or a universal4GiB upper bound.

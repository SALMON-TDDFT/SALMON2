# Fast-PR symmetry exchange checkpoint — 2026-11-18

Internal same-binary V2 extends the admitted scalar periodic CPU/full-q/internal-k1
exchange core. The public symmetric TDDFT open prohibition remains; this is not a
TDDFT/FDTD propagation or coupled-resume claim. Conventional HSE is not supported.

The checkpoint contains value-only symmetry map/cache plus existing numerical
source/ACE/history and kernel identity. An initialized matching fresh owner supplies
trusted static array layouts. Reads compare static map values, stage optional bounded
cache and numerical values, then install. No FFT plans, C pointers or MPI ownership
are copied/replaced. Fresh operator construction is required; accepted propagation
steps are not replayed. Raw streams rely on caller-authenticated archives/config.

Local Apple Silicon evidence (no GitHub Actions):

- Fail-first actual old-library run rejects symmetric checkpoint write (STOP101).
- 32 cases: full/ACE × star/orbits × phase-cache setting0/1 × OMP1/2 × MPI1/2
  processes, each independently owning COMM_SELF cores. Distinct per-process
  cell/screening/call schedules; C2z/mirror-x owners; explicit absent/present caches.
- Validate-only and staged reads preserve live bytes. Installation preserves FFT
  owner identity. Restored action/energy and subsequent refresh/action archive
  bytes match continuous core execution. Closing/reopening the owner also matches.
- Truncated streams, different operation identity and phase partition reject
  without mutation. Additional invalid-cache-range/NaN write guards pass OMP2.
- Old/new non-reduced ACE and full-action V1 checkpoint bytes exact; full-k action
  reference bytes unchanged. No exchange algorithm changed.
- Existing MPI2 GS outputs, five-rank PBE/HSE Taylor/PT-CN, seven-rank mixed32-step
  protected failure/retry remain exact. Mixed protected direct restore save16,
  fail18 preserves all original field/LD/electronic state bytes; continuation and
  three future retry steps exact, zero accepted-step replay.
- Source/installed545 Python tests (one skip each), sdist-to-wheel/native hashes,
  installed mixed/1D fault-retry and packaged native probe retained in evidence.

Evidence: `~/.codex/outputs/salmonpy-symmetry-checkpoint-20261118/`.
These small implementation regressions do not establish physical convergence.
Memory budgets remain advisory. Next: symmetric TDDFT source/state integration,
then internal-k MPI partitioning and symmetric mixed protected failure/resume.
Learned/PT-CN symmetry stays rejected.

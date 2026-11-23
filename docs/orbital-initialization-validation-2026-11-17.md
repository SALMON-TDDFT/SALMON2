# Initial orbital layout ownership — 2026-11-17

Continuation of approved GS/TDDFT initialization ownership extraction. Internal
value-only `s_orbital_initial_config` captures mode,periodicity,spin/SOC,
state/electron counts,temperature and k-mesh/file selection. The same explicit
`resolve_orbital_layout` in standard and fast-PR chooses nspin,no and real/complex
storage without shared input reads. Existing derived-type layouts/public APIs
are unchanged; no whole solver/resource handles are copied.

## Integration and limits

Capture occurs before downstream setup. Resolved stencil orthogonality is passed
explicitly. Fast-PR uses owned occupation settings/effective system kmesh and the
same initially captured legacy hybrid predicate. Standard file k input captures
the legacy -1 mesh sentinel; initial occupation uses captured local electron/SOC
values. Existing orbital/occupation evaluation order and compatibility side
effects remain. No elapsed propagation or new physical approximation is added.

For valid inputs, GS Gamma/non-Gamma/nonorthogonal/isolated/SOC/hybrid rules and
RT finite/negative-temperature branches retain their results. Positive odd
counts use integer ceiling without nelec+1 overflow; spin totals use wider
integer arithmetic. Bad mode,nonfinite temperature or unusable selected layout
reject before assigning the result. This is not full electron/occupation-capacity
certification, SCF-input ownership or general recoverable native initialization.
Remaining k/symmetry/communicator and other input dependencies still apply.

## Native verification

`tests/native/test_orbital_initialization.f90` fails against old modules (missing
symbols), then passes standard/fast-PR/RT-profile builds at OMP1/2 with probe bounds
checking. Manually specified cases cover GS Gamma/non-Gamma/file-sentinel/
nonorthogonal/isolated/SOC/hybrid, RT even/odd/SOC/spin-resolved counts and both
sides of the -1e-12 threshold. Capture followed by poisoned global inputs and
alternating two different outputs leaves each independent. Invalid resolution
preserves the original result; capture never opens the selected k file.
Cell/grid/ion probes also pass all three profiles at OMP1/2. Read-only review
found no blocker; source bodies of the explicit selector/capture match.

Strict fresh preparation verifies standard/HSE/RT-step source inventories and
before/after manifest links. New external libraries recompile initialization_sub
and reuse unchanged matching-profile objects. RT uses its matching earlier full
RT module tree; no cross-profile module reuse. Existing binaries remain fixed,
all commands/source hashes/core reuse receipts are retained. Shared LD/SBE kernel
bytes remain unchanged. No remote Actions.

## Runtime regression

Local Apple Silicon64GB/18-core host, sequential MPI runs,OMP/BLAS1, at most seven
ranks,48GiB advisory budget. Budget excess only warns.

- MPI2 standard Libxc-PBE GS warm-start three-iteration regression: newly generated
  energy/info,eigen,k and force files exactly equal the unchanged old executable.
  Comparison files are removed before running; this is not physical convergence.
- Five-rank field2/local LD + SBE1 + electronic2 with PBE TDDFT, fast-PR Taylor
  and PT-CN; seven-rank independent PBE2/HSE2 mixture:32-step first/warm/final
  late/deep rejection, original state/output preservation and valid retry pass.
  Final arrays and all ordinary streams match unchanged normal baselines exactly.
- Seven-rank direct save16/late restore failure18: original state bytes retained,
  original continues,retry and three future frames exact without accepted-history
  replay. Existing identity/admission negatives pass.
- Source/fresh installed suites545 tests,one skip each;source-distribution→wheel/
  native SHA256 and guide synchronization pass. Installed seven-rank mixed and
  five-rank1D two-TDDFT-group fault/retry arrays and all ordinary streams exactly
  match unchanged normal baselines.

The varied-layout tests isolate the selector, not full physical calculations with
all combinations of spin/SOC/temperature or multiple GS contexts. Tested full
runtime profiles retain the scopes from preceding validation. RT-profile probes
are not a new general RT physical/restart validation campaign. Molecular work is
paused; general initialization STOP/I/O/allocation/process/MPI failures remain
outside recovery.

## Evidence

Persistent local bundle:
`/Users/otobetoshihito/.codex/outputs/salmonpy-orbital-initialization-20261117/`.
Includes source snapshots,new libraries/receipts,compilation/probe commands,old/new
comparisons,output inventories,direct-restore evidence,distributions and test
logs. Local paths record provenance,not user installation requirements or a
published portable runtime release. JP/EN API guides and AI JSON state the same
limited internal ownership boundary.

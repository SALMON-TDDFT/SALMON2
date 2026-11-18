# Active state-isolation checklist

Continuation of handoff §7.3–4 and the approved GS/TDDFT ownership work. This is not a replacement for the wider project backlog. Historical plans record their completion-time scope; later evidence below extends it. Direct mixed-session persistence, non-PZ inner-MPI validation and protected ordinary mixed restoration are completed within the recorded scopes. Wider isolation, resources and distribution continue as separate objectives.

## Latest actionable queue — 2026-11-18 (after functional test cleanup)

This section takes precedence over all earlier dated status/check counts below.
Historical unchecked bullets are scope records, not additional independent tasks.
Verified GS handoff/HSE registry/test cleanup committed as cd0841fe [skip ci].
Ordinary HSE assembly2340787d and Gamma portable exportb66b97fd committed and
pushed to salmonpy/codex/salmonpy-state-inputs-recovery with [skip ci].
Regular-mesh/MPI export committed asd67ce847 and pushed with [skip ci].
Referenced-k-file import fix and declared-input metadata are uncommitted;
existing-GS/native checks below. No new push or Actions this slice.

### Completed within measured scopes

- [x] Explicit different-input fast-PR HSE model registry, complete symmetry/density
  workspace switching, Python leases and direct restart; see latest validation below.

- [x] Standard TDDFT: same-process different-input model registration/activation,
  Python registry and ordinary material-ID/multi-point integration. Different
  cell/grid/k/XC/decomposition and independent/direct-resume controls verified.
- [x] Standard retained GS: independent native/Python owners. Existing-GS load
  without SCF; owned in-memory GS→TDDFT seed, independent RT state creation after
  GS close, interleaved GS loading/RT stepping. This slice does not test new SCF
  interleaved with RT; no new GS/reference data is authorized.
- [x] GS→TDDFT existing-file current parity and direct no-replay resume in5 cases:
  PZ MPI1/orbital2/spatial2, Libxc-PBE orbital2, TB-mBJ spatial2. Analytic
  real/complex empty-band redistribution also passes. Imported convergence is
  unknown and uncomputed forces are rejected.
- [x] Same-model fast-PR HSE full/reduced-k independent states, symmetry Taylor/
  learned/PT-CN, internal MPI1/2, mixed FDTD and protected direct resume/retry.
- [x] Functional test consolidation:620 unit methods (one skipped),83 files;
  LD/SBE and existing-GS MPI runner7 cases pass. Removed redundant old entries.
  Selected49→43 unit methods retain the identical581 traced runtime lines;
  line tracing is not branch coverage or a physics certificate.
- [x] Relevant JP/EN/AI guides, examples, test removal mapping and canonical
  diagnostic records updated. Mixed-HSE source/wheel regressions are620; retained portable export source
  regression is622 (one skip), with no new wheel run for that Python-only slice.

### Remaining work (do not double-count historical notes)

1. [ ] Commit/review and push the currently verified changes when requested;
   retain CI quota-saving policy. HEAD is d67ce847; preceding work pushed; referenced-k-file import fix uncommitted.
2. [x] Connect the different-input fast-PR HSE registry to ordinary material-ID/
   ResponseSession/rank-plan assembly. Grouped4/separate6 and protected aggregate
   save8/resume24 are exact. Specialized SOC/+U/DC remain separate; molecular work paused.
3. [ ] Complete retained GS/SCF checkpoint and portable archive/version migration:
   regular-mesh DFT-owned GS/2 collective export is implemented (synthetic validation only);
   physical GS identity, broader PP format/content coverage, inactive legacy
   checkpoint canonicalization, changed subgroup/decomposition restoration.
   Portable GS/2→new retained GS initialization now works through explicit conversion
   and compatible MPI splits. This is not full SCF/cache or RT checkpoint repartitioning.
   Existing same-topology RT direct restore is already complete in measured scopes.
4. [ ] Remaining shared initialization/I/O/timer/unscoped MPI resources and broader
   recoverable initialization/allocation/I/O contracts. Current numerical rollback
   and protected restore do not imply general STOP/OOM/MPI/process recovery or
   threaded/concurrent GS/RT. New-SCF/RT interleaving remains unverified.
5. [ ] Larger electronic grids/more points and communication/resource calibration;
   Linux/Intel Mac/other runtime builds, packaging and matching distribution
   evidence for the new native features. Budgets remain warning-only.

Smaller/future: band-iteration warning policy; optional interpolation/subcycling
and E-enabled runtime catalog publication. Converged physical benchmarks remain
explicitly separate from implementation parity; new SCF/GS generation needs
applicable authorization. No conventional HSE requirement or routine CI dispatch.

Evidence: [memory bridge](../memory-gs-tddft-validation-2026-11-18.md),
[test runner](../component-tests.md), [functional removal map](../functional-test-cleanup-2026-11-18.md).

## Retained standard GS handles and Python — 2026-11-18

- [x] Stable heap owners for GS numerical arrays, XC and child communicators;
  default temporary GS and optional retained model run share the original SCF.
- [x] Generation-tagged native create/info/shape/read/close/resources and Python
  NativeGSBackend/session/state. Immutable material ID/calculation and fixed
  runtime/decomposition. Requested arrays are read-only copies, valid after close.
- [x] Ten real MPI/OMP cases, PZ/Libxc-PZ/PBE, real/complex orbitals, different
  cells/grids/atoms/electrons/k; interleaved A/B and independent arrays exactly
  match. Close/reuse/stale/mixed tokens, sibling preservation and resource release pass.
- [x] Native SCF flag snapshots preserve exact convergence semantics; unfinished
  SCF warns and retains. Existing Fortran GS9 cases match old status/eigen/info/k.
- [x] Darwin NumPy/Accelerate complex-BLAS ABI interposition reproduced and fixed
  with two-level library binding and actual shared libmpi COMMON sentinels.
  Host Fortran probes use the same link helper. Linux link path unchanged.
- [x] Rank-local work visibility/content and chdir failures agree before native
  collectives. Existing TDDFT currents and direct resume stay exact in MPI6 controls.
- [x] Paired retained-GS guide, public C header, Python sample and AI JSON updated.
- [x] Source/final-wheel635 tests each (one skip). Installed MPI6 real/complex
  visibility-rejection test, public MPI2 sample, different per-state SCF cap/
  tolerance test and true-convergence-flag test pass. Array values agree with
  independent controls at the same MPI partition; no physical-convergence claim.
- [x] Selected hash-only/diagnostic evidence, approved runtime/wheel and explicitly
  unapproved failed-library evidence saved to canonical
  calc2026/salmonpy-retained-gs-20261118; no GS restart/orbital/density collection.
  Local OneDrive writes verified; cloud sync unverified. No Actions dispatched.
- [x] In-memory standard GS→TDDFT handoff, existing-GS import and independent RT states;
  see 2026-11-18 memory bridge evidence below.
- [ ] Retained GS checkpoint; HSE/SOC/+U/DC GS, concurrent GS/TDDFT and broader
  GS initialization recovery remain separate.

## Standard GS model-owned input boundary — 2026-11-18

- [x] Internal s_gs_model.load/run/clear retains owned parsed values and calls
  existing SCF with explicit selected inputs; rejected unsupported replacement
  preserves the model. Default main_dft arithmetic/entry remains unchanged.
- [x] Caller input/calc_mode/unit views, actual symmetry switches, optimization
  and XC scratch switch/restore. Model-only GS retires communication requests
  before local buffers and releases scoped child groups after each run.
- [x] Nine actual GS cases: PZ/PBE, Libxc-PZ/TB-auto, PBE/TB-fixed × MPI1/orbital2/
  spatial2,8³/12³ cells/atom/electron/k controls. A→B→A and independent old status,
  energy/residual/iterations/eigen/info/k output bytes exactly match. Six repeats
  retain zero live child groups and stable warmed management slots.
- [x] Two existing TDDFT registry/control MPI6 cases: interleaved currents and
  direct no-replay resume unchanged. Source628 tests (one skip); paired guide/AI
  metadata and reusable development probe updated. Final runtime/installed probe
  repeats all9 GS cases; installed628 tests (one skip), installed TDDFT MPI6 and
  legacy default GS spatial2 match old controls. Evidence: canonical
  calc2026/salmonpy-standard-gs-inputs-20261118; cloud sync unverified.
  No Actions; no GS restart written.
- [x] Retained standard GS/Python handles are implemented in the newer section
  above. Full specialized GS contexts, HSE/SOC/+U/DC and concurrent GS/TDDFT remain separate.

## Ordinary heterogeneous layer groups — 2026-11-18

- [x] model_groups maps explicit group IDs to material-ID lists; normal layer
  planning preserves ID/point order and per-material input/internal decomposition.
- [x] NativeTDDFTBackend(registry_model=True) and session-owned group resources
  register one standard model per binding; state/factory/resource close order.
- [x] Local LD + SBE + two PBE/two TB-mBJ points: separate10 versus combined4
  MPI ranks,128-step A/E/j exact. Direct save63/resume65 exact, no replay.
- [x] Late last-material restore failure preserves all original states/outputs;
  retry succeeds. Paired manuals/AI JSON/sample group file updated.
- [x] Source/final wheel628 tests each (one skip); installed128-step continuous
  and final-wheel direct resume65 steps match source continuous exactly.
  Evidence: canonical calc2026/salmonpy-layer-registry-20261118 (no GS collected;
  cloud synchronization unverified). No Actions used.
- [x] Two unsupported native HSE owner cases rejected during planning; custom
  backends keep their existing callback contracts and explicit grouping.

## Heterogeneous standard TDDFT registry — 2026-11-18

- [x] Explicit native generation-tagged model registration/spawn/close; input and
  initialized views, optimization arrays and XC workspace owned per model.
- [x] Python open_registry/open_model, deterministic collective model selection
  and cleanup; closing one owner leaves live siblings and other models intact.
- [x] Eight heterogeneous MPI/OMP cases: different cell/atoms/grid/k/XC/dt/splits,
  independent old-runtime controls and direct resume agree exactly. Fifteen legacy
  default-model MPI comparisons also pass. No GS computed or collected.
- [x] Source and installed wheel:621 tests each (one skip); installed MPI sample
  and six-rank independent comparison pass. Source and wheel capacity tests let
  model B reach56 steps beyond model A capacity48 without cross-contamination.
- [x] Paired registry guides, public C header, AI JSON and sample synchronized.
  Source/object audits pass; no Actions used. Selected evidence (no GS) stored in
  canonical calc2026/salmonpy-heterogeneous-models-20261118; cloud sync unverified.
- [x] Explicit model_groups in ordinary material-ID plans builds a shared standard
  registry via ResponseSession; per-ID inputs/decomposition and per-point state retained.
  Built-in fast-PR HSE heterogeneous grouping remains unsupported.
- [ ] Full GS contexts and remaining specialized SOC/+U/DC/process resources;
  concurrent threads, HSE registry and portable changed-topology restart are separate.

## Standard parser workspace/default ownership update — 2026-11-18

- [x] Model-owned s_input_parser: scratch units, namelist status/coordinate flags,
  unit-dependent fresh defaults and final selected scalar-TDDFT values. Native
  model construction uses its owned parser. Legacy no-workspace paths/MD named
  coordinate temporary files are preserved; synchronous value view remains.
- [x] Actual parser A→B→A and reused workspace, omitted XC/k/pulse/mixing,
  changed units/atoms/Cartesian-external-reduced-absent sources, two interleaved
  scratch sources and handle release at MPI1/2. Fresh legacy diagnostic logs
  match exactly; existing named temporary-file sentinels remain untouched.
- [x] Five XCs × MPI1/orbital2/spatial2:15 new/old trajectories, both-point RT
  numerical columns and energies match exactly, including direct restart.
- [x] Normal two-PBE/two-TB-mBJ mixed layout: independent/shared/protected direct
  resume and late failed-restore/retry match; source/installed617 tests each,
  one skip. Source-distribution→wheel and installed parser tooling verified locally.
- [x] Heterogeneous registration/activation and Python handles completed in the
  registry section above. Concurrent threads and specialized resources remain separate.

## Standard native model-input ownership update — 2026-11-18

- [x] Model-owned initialization record (227 scalar/array settings including atom
  coordinates/species, units, PP/XC/cell/grid/k and effective dt/parallel inputs),
  deep allocation ownership and explicit `open_initial_state(owner,...)` boundary.
- [x] Actual Fortran two-input capture/clone/resize/release tests. New/old traces
  match exactly at PZ, Libxc-PZ/PBE, automatic/fixed TB-mBJ × MPI1/orbital2/spatial2;
  direct resume and independent sibling lifetimes retain parity (15 cases).
- [x] Scalar IEEE restored-A checks retain all6 original elements; two-rank NaN/Inf
  corruption at each element rejects before publication (12 cases), valid resume exact.
- [x] Standard registry/different-model activation and heterogeneous Python handles
  completed in the newer registry section. Remaining specialized/common ownership
  is tracked only in the latest queue above.

## Different TDDFT input audit — 2026-11-18

- [x] Different inputs can use separate native MPI groups. Normal sample now accepts
  an offline per-material-ID catalog for bundle/input/decomposition/shared-point
  width; serial controls, asset identity and resource advice follow each material.
- [x] Two distinct Si XC inputs (PBE/TB-mBJ), two points each, internal MPI2:
  shared6 vs separate10 ranks match bitwise; protected63→128 direct resume
  matches continuous with zero replay, including late-rank failure/retry.
- [x] Source and installed614 tests each (one skip); installed catalog restart
  ignores unused fallback paths/flags, and old no-catalog direct archive resumes.
  Installed unequal per-material MPI (PBE2/TB-mBJ1) agrees within2e-12 absolute
  tolerance (largest current difference8.61e-19); no new GS/native build required.
- [x] Standard same-process heterogeneous model registry/activation, Python wiring
  and differing-input current/restart validation are completed in the newer registry
  section. Fast-PR HSE heterogeneous input sharing and broader common ownership
  remain in the latest queue; do not count the standard registry twice.
- [ ] Generic allocation/I/O/MPI failure recovery, canonical archive migration,
  larger-scale calibration and additional platform distribution remain separate.

## Standard Python shared-model update — 2026-11-18

- [x] `NativeTDDFTBackend.open_model` / `NativeTDDFTModel`, one immutable model/work
  copy, independent point states, common dt/capacity/decomposition.
- [x] Opt-in normal layer/session grouping; independent/shared and direct resume
  are bitwise equal for LD/SBE/two-TDDFT128-step mixed run (inner orbital MPI2).
  Handled late-rank restore failure retains all owners; retry matches continuous.
  The combined LD/SBE/TDDFT2/HSE2 layout also matches independent execution
  and direct protected restart (10 independent vs6 shared ranks,128 steps).
- [x] Five standard XC profiles × MPI1/orbital2/spatial2:15 public-model
  independent/direct-continuation cases pass; final output-collision retry and
  mixed-epoch shard rejection verified in three MPI configurations.
- [x] Per-state charge-current sign, checkpoint current metadata admission and
  owner-local output journaling; default legacy create_state preserved.
- [x] Paired guides, AI catalog and installed module/example/tool inventory updated.
- [x] Standard different-input sharing completed via explicit registry; common model
  parameters remain per-model, not per-state. Specialized/common work is tracked
  in the latest queue.

## Earlier shared-model audit — 2026-11-18 (superseded by latest queue)

The dated sections below retain the scope at their completion time. Older statements such as "occupation remains shared" or "multiple field ranks remain open" are superseded by the later evidence, not separate unfinished tasks.

### Actionable status — Python shared-model audit

This is an earlier queue retained as evidence. Use the latest actionable queue
at the beginning; do not count historical unchecked items again.

- [x] Delivery throughbd27f0ed is committed and pushed to the independent
  SALMONpy branch. Ownership9f07fb54 and C handles8fe6cd15 are included.
- [x] Low-gradient Libxc HSE correction, full/reduced-k parity and internal MPI1/2.
- [x] Cell/atom/grid/k-star/representative weights, PP references/KY records and
  legacy GS structural admission. This is not complete physical GS authentication.
- [x] General-cell native initialization and actual Cartesian x/y/z/circular A.
- [x] Symmetry Taylor strict/learned and PT-CN off/learned with direct history resume.
- [x] Per-context step/restore buffers and generation-checked same-model C handles.
- [x] Explicit Python same-model shared model: independent dt/capacity/state/output,
  one work copy/attachment, idle recreation and direct save/resume. Independent
  and preceding singleton trajectories match in six MPI/OpenMP cases.
- [x] Fast-PR HSE rt_open MPI children are context-owned and retired at close;
  bounded registry reuse, independent sibling/control, full/reduced-k MPI1/2
  verified. Unscoped/common ownership remains separate.
- [x] Standard embedded TDDFT initialization/initial-current children retire at
  close, including pending payload cleanup. PZ/Libxc-PZ/PBE/TB-mBJ auto/fixed,
  orbital2/rgrid2 match old traces; one lifetime/process stays in effect.
- [x] Standard same-model native C handles: per-context protection/pending/XC/A/
  serializer ownership, generation tokens, original-GS recreation and direct
  restart. Fixed common dt/capacity/decomposition. Python shared factory and
  ordinary standard multipoint session integration are now implemented opt-in.
- [x] Latest source/installed suites:617 tests each, one skip each; sdist→wheel,
  hashes, installed MPI/example, paired guides and AI catalog verified locally.

Remaining work, grouped rather than duplicating historical bullets:

1. [x] Connect the explicit shared-model factory to normal ResponseSession /
   rank-plan assembly: multiple local macropoints, LD/SBE/TDDFT/HSE mixtures,
   FDTD coupling and aggregate protected save/resume. Opt-in grouping is implemented; default ordinary plans remain unchanged.
   Taylor/PT-CN, PBE mixture, distributed field and protected failed-restore retry pass.
2. [x] Multiple reduced-k/symmetry HSE contexts in one immutable model/process.
   Taylor strict/learned and PT-CN off/learned, MPI1/2, independent controls,
   full/reduced-k and direct restart verified. Shared mixed restore/retry covered
   for LD/SBE/HSE128-step compatible local fields; arbitrary PBE fields may break
   retained symmetry and require full-k, as recorded in the validation.
   Different symmetry groups within one process remain in item3, not enabled here.
3. [ ] Standard native C handles now connect to an opt-in Python shared model and normal
   multipoint sessions (see SHARED_TDDFT). Model input records/explicit initialization
   ownership are complete; remaining: broader GS/TDDFT different-input activation,
   unscoped/common MPI resources, common initialization/I/O/timers and specialized
   SOC/+U/DC resources. Fast-PR and standard-TDDFT initialization-child retirement are complete.
   Same-model HSE token identity is complete within its documented boundary.
4. [ ] Extend defined recovery beyond the verified numerical/staged-restore
   subset: initialization, allocation/I/O and unsupported native error paths.
   Process/MPI failures do not acquire a general retry guarantee from these tests.
5. [ ] Canonical/interoperable archives and version migration, physical GS
   identity, non-KY content coverage, arbitrary-backend identity and changed-MPI
   topology restore. Existing same-topology direct restart and GS/2 are complete
   in measured scopes. Canonicalize inactive legacy checkpoint workspace too.
6. [ ] More points/larger electronic grids and communication scaling, calibrated
   resource advice, additional native platforms and matching runtime distribution.

Band-iteration non-convergence policy remains a separate smaller audit item.
Interpolation/subcycling is optional later scheduling work, not a prerequisite
for fixed-dt mixed calculations. Converged physical benchmarks remain material /
user-owned. Molecular work is paused; budget excess warns; conventional HSE and
routine Actions are not pending requirements.

Evidence: [normal shared mixed-session validation](../validation/2026-11-18-shared-hse-mixed-session.md) [shared symmetry validation](../validation/2026-11-18-shared-symmetry-hse.md), and [context MPI retirement](../validation/2026-11-18-context-mpi-ownership.md), and [standard TDDFT MPI retirement](../validation/2026-11-18-standard-tddft-mpi-ownership.md), and [standard native handles](../validation/2026-11-18-standard-tddft-handles.md).

### Fast-PR HSE priority

Fast PR remains the only HSE engine, with symmetry central. Reduced-k symmetry,
actual-A admission, general cells and learned Taylor/PT-CN persistence are verified
on the established internal-MPI1/2 routes. Explicit C/Python same-model sharing is
now available for full/reduced-k states of one immutable symmetry model. The factory
is integrated into normal mixed point/rank assembly. Next address remaining
different-input ownership and remaining common/context resources.
Same-model full/reduced-k multi-context symmetry is complete within the measured
scope; different-input or different-group sharing remains open. Standard-solver symmetry arithmetic remains
a distinct port; no physical-convergence or cross-platform certificate is implied.

### Completed within measured scope

- [x] GS mixing, normalization and occupation inputs; HSE electron-count/temperature/state-capacity routing; TDDFT electron-count diagnostics. Broader initialization and shared controls remain below.
- [x] Mixed LD/SBE/managed Si TDDFT, with exact serial/MPI field/current comparison.
- [x] Native E/H field-local two-material LD + SBE + internally parallel TDDFT (1/2 orbital ranks), serial response replay and pure-LD regression.
- [x] Two-way x/y/z FDTD decomposition with SBE1/TDDFT1: exact A/E/j parity to a single field rank. Installed-wheel samples and JP/EN/AI instructions updated.
- [x] Existing configured material adapters, serial/MPI point identity and managed TDDFT point replay are available in the documented M4/M5/M6 subsets. Do not describe them as entirely missing.
- [x] Same-topology direct field/PML/LD/SBE/TDDFT restart: PZ, TB-mBJ automatic/fixed coefficient and Libxc-PBE include native TDDFT orbital MPI2;63 saved/65 resumed steps match continuous128 bit-for-bit. Libxc-PZ is validated with TDDFT1. Generic archives and other HSE partitions remain below.

- [x] Fast-PR HSE serial Fortran context direct checkpoint: ACE/exact/learned/PT-CN histories and RT buffers, fresh-process6/7→12,32 cases overOMP1/2; listed internal traces match. Native-host/k-only inner-MPI HSE and outer serial-worker coupling are covered below.

- [x] Fast-PR HSE outer-MPI/FDTD coupling via one serial worker per material rank: field1/2+SBE1+HSE1 with local LD, Taylor/learned PT-CN direct127→256; all four trajectories match continuous exactly. The separate native k-only inner-MPI entry is covered below.

- [x] Fast-PR HSE native host-MPI k-only groups1/2, field1+SBE1+HSE2 with two field-local LD materials. Taylor/learned PT-CN direct127→256: all A/E/j histories match continuous bit-for-bit, only129 new step calls. HSE-only continuation current/energy also match exactly; rank1 truncation and swapped shards are collectively rejected. Ordinary serial HSE regression16 cases remains exact. Evidence: `../fast-hse-kmpi-validation-2026-11-16.md`. Orbital/rgrid HSE, changed topology and general native-step recovery remain open; handled restore failures are covered by the later protected-restoration evidence.

- [x] Common checkpoint admission/execution and backend-owned isolated-candidate restore. Ordinary native mixed direct resume now protects field/PML/all local LD, SBE and native TDDFT/fast-PR HSE against handled staging failure; retry and fresh CLI continuation verified. Source/wheel502 tests, one skip each. Latest evidence: `../mixed-protected-restore-validation-2026-11-17.md` (fe286bc1).
- [x] GS/2 Zarr-v3 array archive, strict Gamma export/import and opt-in serial multi-k bridge exist within documented periodic/unpolarized/fixed-ion scopes. Do not list GS archive v2 as wholly missing. MPI/spin/DC/general portability remain separate extensions; see salmonpy/GS_ARCHIVE.en.md.

### Delivery status and remaining mixed-material work

1. [x] Delivery: native/distributed-field samples and evidence pushed in b5b2053c, layer setup in01d43cc1, direct persistence in1851c533, non-PZ persistence inea22fe15. No routine Actions.
2. [x] Reject longitudinal incident components in Multiscale1D prefilled histories, including zero boundary-sample cases. Regression and valid transverse preservation verified.
3. [x] Configuration-to-native-mixed setup: explicit material IDs/models or bindings, ordered geometry/thickness and rank budget via assemble_native_layers. Initial-state preparation remains explicit, not an automatic arbitrary backend factory.
4. [x] Same-topology supplied mixed fixture direct checkpoint/restart: native field/PML/local LD/SBE/TDDFT numerical state and coupling history/clocks/point IDs captured together; no elapsed propagation replay. PZ and non-PZ measured scopes below. Generic user-backend archives and changed topology remain open. Handled mixed restore failures are now protected as recorded below.
5. [ ] General native failure contract: scoped semilocal TDDFT, fast-PR HSE and native E/H numerical recovery plus explicit aggregate numerical/ordinary-output step rollback are implemented below; initialization/unsupported native errors and process-terminating paths remain open. Handled mixed restore staging failures already preserve originals and permit continuation/retry; do not list that completed subset as missing. The CLI still exits on an unhandled error.
6. [ ] Remaining GS/TDDFT initialization, ionic/I/O/timer/DC/SOC/+U/localized-exchange ownership and stable native handle interface. Required where the chosen backend/capability needs it; full reentrancy is not a prerequisite for already tested MPI/process-isolated mixtures.
7. [ ] Remaining scaling/resource work: larger numbers of points/LD cells, calibrated backend memory models and communication scaling. Aggregate CPU memory advice including restore candidates is now implemented below. Budget excess only warns by user decision; hard enforcement is not a pending requirement. Four-core runs do not certify a universal4GiB cap.
8. [ ] Additional environment/topology validation: field2+SBE1+TDDFT2 (5 ranks) and96³ field12+SBE1+TDDFT2 (15 ranks) are now validated below; new native mixed paths still have local Apple Silicon evidence, not fresh Linux/Intel Mac validation or a published matching runtime catalog. Windows native support/delivery policy is distinct from existing lightweight wheel support.
9. [ ] General coupled archives: interoperable archive layout/version migrations, full arbitrary-backend identity propagation, changed-topology restore and complete historical native output streams. Common admission/execution and protected restore already exist; the current native sample remains same binary/input/topology. Unsupported XC/resource profiles and HSE orbital/rgrid decomposition require their own validation.
10. [ ] Optional later scheduling: interpolation/subcycling and heterogeneous dt. Current mixed paths use explicit fixed clocks and matching dt.

Molecular feedback/energy-schedule work remains paused. Physical convergence for each user's material is user-owned; do not turn it into a mandatory developer sweep. This audit does not close the broader historical backlog.

## Completed within documented scopes

- [x] Independent native HSE resources and TDDFT interleaving regression.
- [x] Per-state timestep, capacity, OpenMP selection, XC workspaces and stencil tables.
- [x] Value configurations for cell/grid, k-point, radial pseudopotential and XC inputs; scalar PP-grid ownership.
- [x] Owned XC phase and pre-SCF counters; serial GS and real DC/MPI regression.
- [x] SCF exchange flags/localization checks routed to the selected owner.
- [x] Ordinary eligible GS creates/releases its HSE owner through result/checkpoint output.
- [x] Two complete same-input GS calls in one process match original separate-process output, including generated restart payloads. See `../gs-lifecycle-validation-2026-11-15.md`.

- [x] Retained GS storage and internal run/close entry; selected retained-state arrays checked while another same-input GS runs/closes. Explicit input and broader independence remain open.

- [x] Value-owned GS SCF iteration cap, convergence threshold and pre-SCF controls, with semantic validation and poisoned-global regression. Other shared inputs remain open.

## Still open in the existing isolation objective

- [ ] Complete GS context with explicit input/state, rather than the global-input `main_dft` entry; different-input independence remains unverified.
- [ ] Remove remaining guarded process settings and shared initialization dependencies from full TDDFT; common grid/k/XC constraints remain.
- [ ] Own MPI communicators, DC material/communicator state and remaining fixed I/O/timer resources.
- [ ] Isolate specialized SOC/+U, localized exchange and remaining unsupported resources before extending the independent path.
- [ ] Replace process-terminating native errors with a defined recoverable contract and validate failure atomicity.
- [ ] Expose the independent solver through the planned stable handle ABI, with save/restore/lifetime guarantees.

Existing regression cases are not physical-convergence or cross-platform evidence. Checkpoints matching on output do not prove TDDFT import/restart compatibility. Do not mark a whole item complete based on a narrower fixture. Keep day-to-day validation local to conserve Actions quota.

The approved direction for per-material settings and user backends is specified in `2026-11-15-material-config-backend-design.md`. It guides existing input extraction; the configuration-only milestone below does not complete solver input isolation.

## Material configuration milestones (approved; alongside input isolation)

- [x] M1: configuration-only `salmonpy.materials/1`, explicit custom resolver registry, immutable per-material definitions, initial-state ID reference checks; local tests and packaged example.
- [x] M2 (ASE Atoms and local CIF/extxyz; other formats remain unsupported): structure adapter for ASE Atoms and explicitly referenced files, including supercells; preserve cell/atoms, source/format/frame/hash; no implicit supercell generation.
- [x] M3 (positive axis-aligned periodic cells; reject other width geometries): per-material maximum grid spacing with units → ceiling counts plus backend constraints, report actual spacing; retain explicit counts, reject dual specification. Validate nonorthogonal/isolated geometry separately.
- [ ] M4: built-in LD/SBE/DFT/TDDFT configuration adapters with model-specific schemas, units, provenance and initial-data compatibility; preserve existing A/c→j and LD E contracts.
- [ ] M5: ID-carrying bindings and independent per-point states; connect existing serial and MPI adapters, capabilities and execution plans.
- [ ] M6: preserve material/point/initial/run ID mappings in output and restart; validate hashes and backend compatibility.

M1 is not full simulation schema support. M2 uses MateriApps as a discovery portal, with actual structure data coming from the chosen provider or user file. These milestones do not relax native shared-input guards or promise physical convergence.

M1 evidence: `../material-config-validation-2026-11-15.md`. Structure/grid adapters were subsequently added within M2/M3 scope below.

- [x] M1b: explicit local input-file staging and portable hash/provenance manifest; offline verification without backend execution. M2/M3 structure and grid semantics use explicitly registered adapters, not automatic dependency discovery.

M1b evidence: `../offline-input-preparation-validation-2026-11-15.md`.

M2/M3 evidence: `../structure-grid-validation-2026-11-15.md`. Nonorthogonal/isolated width resolution and automatic scientific backend preparation remain unsupported.

- [x] M4a: explicitly registered LD/SBE parameter adapters and factory, SBE source hashes/revalidation, existing E versus A/c contracts. DFT/TDDFT adapters remain M4.
- [x] M5a: material ID and definition digest carried by configured LD/SBE models and state proxies. Point/run IDs, MPI and automatic output/restart propagation remain M5/M6.

M4a/M5a evidence: `../builtin-material-validation-2026-11-15.md`.

- [x] M4b: DFT configuration resolver with structure/PP hashes, XC/grid/states/kmesh and neutral electron checks; UPF element check. PP/XC provenance remains user-verified.
- [x] M5b: configured DFT result → verified GS-to-TDDFT bridge → independent worker states retain material ID/digest. Configured warm-start results unwrap into existing compatibility checks.
- [x] M4c (GS/1 scope): explicit initial_state_id selection with pinned local source, exact configured-DFT compatibility and executable hash checks; ID/digest propagate to independent TDDFT states. No implicit prepare or general restart conversion.

M4c evidence: `../saved-gs-selection-validation-2026-11-15.md`.

M4b/M5b evidence: `../dft-material-config-validation-2026-11-15.md`. M5/M6 point/MPI/disk identity work remains open.

- [x] M5c (serial scope): PointRun binds ordered point IDs, owns one state lifetime per run/point, carries material/initial identity, and accepts explicit custom Python models with record/new_state. No hidden-sharing certification or MPI allocation.
- [x] M6a (postprocessing scope): outputs/1 archives retain point-run/1 metadata and explicit array point axes; readback validates identity consistency and axis lengths. Caller supplies ordering; this is not restart restoration.

M5c/M6a evidence: `../point-run-validation-2026-11-15.md`. MPI and restart aspects of M5/M6 remain open.

- [x] M5d (explicit MPI identity): ResponseSession/from_multiscale validate user point/run/material/initial records against an explicit integer transport map, agree metadata before callbacks, expose immutable group-local IDs. Numerical and native capability contracts unchanged.
- [x] M6b (MPI postprocessing): collective output saves identity and rank/group ownership mapping; field-leader I/O errors are collectively propagated, reader validates mappings. No restart restoration or implicit provider certification.

M5d/M6b evidence: `../mpi-point-identity-validation-2026-11-15.md`. Remaining M5/M6 work includes automatic model/provider integration and restart identity/compatibility.

- [x] M5e (endpoint-A subset): bind_recorded_material/response_identity derive binding and MPI point metadata from configured SBE or recorded user models; immutable records, pre-prepare/create drift checks and session agreement. One rank per material group. Dedicated LD/TDDFT providers and collective native capability integration remain open.

M5e evidence: `../recorded-material-binding-validation-2026-11-15.md`. Empty-message callback exceptions now produce nonempty collective diagnostics; bounded create/step fault cases pass.

- [x] M5f (dedicated serial-native providers): configured LD endpoint and managed TDDFT factories retain immutable IDs, explicit LD conversion and verified TDDFT runtime; TDDFT initial snapshot is borrowed without extra prepare-time copy. Installed three-rank LD/PZ-Si currents match direct paths.

M5f evidence: `../configured-native-binding-validation-2026-11-15.md`. Collective native provider capabilities, restart identity and full native isolation remain open.

- [x] M6c (managed TDDFT point disk resume): PointState wraps existing replay checkpoint with exact run/point/material/initial identity and inner-manifest hash; mismatch rejected before replay. Installed Si rebuilt-model restore/continuation matches uninterrupted currents; cross-point restore leaves target unchanged.

M6c evidence: `../point-checkpoint-validation-2026-11-15.md`. MPI whole-session restart, other backend disk formats, branch semantics and in-memory point snapshot identity remain open.

## Common field input (user-approved A/E/H extension)

- [x] F1: immutable local A/c/E/H packet with None distinct from zero; per-component time/position/unit/sampling, TimeSpec-based requirements and declared gauge.
- [x] F2 (local adapter): explicit legacy endpoint A/E compatibility and custom multiport step_fields; declared backend clocks checked; native LD/SBE direct parity. No implicit conversions or magnetic-physics claim.
- [x] F3a (experimental standard runtime): dedicated E1/STEP_E protocol, PZ isolated fixed-ion midpoint, dipole observables; FieldPacket driver and reference/interleaved/history-replay checks. Validation-only tools, not an installed public loader.
- [x] F3b (local E1 scope): binary-bound validation loader, public isolated native initial-state snapshot, independent E states and in-memory atomic replay restore.
- [x] F3c1 (local replay): isolated E/dipole disk checkpoints with rebuilt-model
  compatibility checks and candidate-worker atomic restore. See
  `../isolated-disk-checkpoint-validation-2026-11-15.md`.
- [ ] F3c2: published E-enabled runtime catalog; self-field physics remains caller-selected.
- [x] F4a: explicit MPI A/E/H packet transport and per-point required-port,
  position/time/unit agreement; current-sign/output declaration checks; legacy A
  regression. See `../mpi-field-packet-validation-2026-11-15.md`.
- [x] F4b1 (point samples): immutable SampledField and sample_packet with explicit
  convex spatial weights and exact/bracketed-linear temporal sampling. Manufactured
  affine fields, separate port clocks and off-stencil rejection verified.
- [x] F4b2a: local detached A/E samples from native E/H and MultiscaleField,
  preserving separate clocks and mapped positions; native 1/2-field-rank checks.
  Evidence: `../field-extraction-validation-2026-11-15.md`.
- [x] F4b2b-H: explicit collective raw-H accessor at (n+1/2)dt, owned-face
  averaging, old ABI rejection and rank-local preflight failure propagation.
  Evidence: `../native-h-sampling-validation-2026-11-15.md`.
- [x] F4b2c (exact mapped samples): PacketCoupledStepper with native/1D port
  and current clock admission, legacy-A parity and synthetic E/H feedback.
  Evidence: `../packet-coupling-validation-2026-11-15.md`.
- [ ] F4b2d: general interpolated/subcycled schedules and explicit adapters for
  existing electronic/dipole providers. No implicit normalization or magnetic model.

F1/F2 evidence: `../field-input-validation-2026-11-15.md`. These do not complete native E/H physics or MPI transport.

F3a evidence: `../isolated-worker-validation-2026-11-15.md`. RT step-config overlay rejects E1; HSE/embedded routes are not E1-enabled.

## One-way molecular radiation

- [x] R1: conservative explicit dipole-to-macro-cell polarization/current map, immutable atomic-unit contract and molecule/point IDs.
- [x] R2: precomputed dipole trajectory to NativeEHField, explicit initial hold and midpoint current times; collective preflight before updates, no native numerical change.
- [x] R3: native impulse sign/deposition, one/two field-rank parity and recorded C2H2 one-way smoke test.
- [x] R4a: finite-window analytic point-dipole refinement and source-work comparison; finest field error 4.3%, not production convergence.
- [ ] R4b (development): explicit self/local-field selection and two-way molecular coupling contract, with sign/time/conservation regressions. Production boundary/flux/physical convergence is user-owned, not a pending developer sweep.

Evidence: `../dipole-radiation-validation-2026-11-15.md`.

R4a evidence: `../radiation-refinement-validation-2026-11-15.md`. Self-field design: `2026-11-15-molecular-self-field.md`.

F3b evidence: `../isolated-public-api-validation-2026-11-15.md`. Physics convergence
is the user's task; no further convergence sweeps are part of this API work.

## Backlog audit after the sample review (2026-11-15)

The unchecked M4–M6 headings are umbrella milestones, not statements that their
subtasks are unimplemented. Preserve the scoped completion records above:

| Milestone | Implemented scope | Remaining scope |
|---|---|---|
| M4 | LD/SBE adapters, configured periodic DFT, GS/1 selection and managed TDDFT bridge | Broader initial-data/model compatibility, including isolated configuration integration; do not infer support from the local E loader |
| M5 | Serial per-point states, recorded/custom bindings, MPI identity and dedicated serial-native providers | Collective native provider capability integration; full native state/input isolation |
| M6 | Serial/MPI postprocessing, managed TDDFT point resume, ordinary native mixed-session restart identity | Other backend disk formats and branch semantics; process-local point identity completed below |

- [x] S1: audit all 39 Python sample entry points; run 13 lightweight examples,
  two preparation examples, serial two-layer FDTD and three three-rank MPI cases.
  Other input-dependent workflows received CLI checks only, not fresh scientific runs.
- [x] S2: repair six MPI help entry points, update paired sample indexes and add
  a regression that blocks mpi4py. Full source suite: 409 tests, one skip.
- [x] S3 (Apple Silicon/macOS/Python 3.14): build from sdist, install a fresh
  wheel, run 409 tests (two explicit skips), all 39 sample execution/help checks
  and archive-to-dipole handoff; public native isolated replay/independence/recovery
  also passes. No Linux/Windows/Intel Mac validation or publication in this step.
- [x] Delivery: read-only independent review and local commit of the isolated
  E/radiation/sample changes; no required correctness fixes identified. Review
  inspected code/tests; execution evidence remains the separate validation records.
  No push or Actions run is implied by local completion.

Sample evidence: `../examples-audit-2026-11-15.md`.

Recommended continuation within the approved work: F4b explicit sampling/native
coupling integration; see `2026-11-15-native-field-sampling.md` for the audited
clock/accessor gaps. F3c1 and F4a additions have completed local review and
validation and are recorded in a separate local commit (no push/Actions).
Publishing an E-enabled runtime is F3c2 and requires its own
packaging/platform evidence. The original native GS/TDDFT isolation objective
above remains open; the process-isolated Python APIs do not complete it.

Physical convergence for each user's material/geometry is explicitly user-owned.
Developer work supplies numerical regression, interface checks and usable samples;
it does not continue the radiation refinement study as a prerequisite for API work.
This is an active-work audit, not a declaration that the wider historical backlog
(non-collinear/Ehrenfest, broader XC/platform coverage, etc.) is complete.


- [x] F4b2d-A: explicit existing A/j binding adapter, endpoint/Yee mappings,
  independent SBE states and two-point managed TDDFT legacy parity in 1D/native
  E/H. Evidence: `../electronic-packet-adapter-validation-2026-11-15.md`.
  General interpolation/subcycling and dipole adapters remain open.


- [x] R4b-selection: pure explicit molecular drive selection, own-source vacuum
  reference provenance and exact sample checks, immutable records and examples.
  Evidence: `../molecular-field-selection-validation-2026-11-15.md`.
- [ ] R4b-schedule: causal native/reference oscillator loop and source-work
  accounting before attaching isolated TDDFT. R4b as a whole remains open.


- [x] R4b-causal-fixture: separate native/reference oscillator loop, illumination
  baseline, source-work diagnostics and temporal-error isolation. Evidence:
  `../molecular-feedback-schedule-validation-2026-11-16.md`.
- [ ] R4b-energy: held-drive work discrepancy is first order and remains large;
  design and verify centered coupling before isolated TDDFT attachment. Native
  energy/flux conservation and material-scattering validation remain unclaimed.

## User priority update — 2026-11-16

Molecular work is paused at the user's explicit request. Preserve R4b evidence
and open limitations, but do not continue molecular feedback, energy-schedule
refinement or isolated molecular extensions until requested again.
Return focus to complete GS/TDDFT explicit input ownership, collective native
backend contracts, recoverable errors/handle ABI and whole-session restart.
Existing GS/2 serial multi-k conversion and process-isolated replay are already
implemented in their documented scope; do not list them as wholly unimplemented.

- [x] GS metric selection: convergence_metric joins owned GS SCF settings;
  all five residual definitions and restart outputs match the original fast-HSE
  reference with poisoned shared controls. Evidence:
  `../gs-convergence-config-validation-2026-11-16.md`.
  Mixing/normalization inputs and general material initialization remain open.

### 2026-11-16: GS mixing settings

- [x] Capture method/rates/Pulay memory/auto adjustment/no-mixing stage per retained GS state.
- [x] Preserve capacity21 and formulas/order; reject malformed explicit values before activation.
- [x] Compare 13 two-state cases (different mixing coefficients/history length) with original solver outputs, restart bytes and residuals.
- [x] Retained settings/adaptive rates/history independence; separate halving/reset check.
- [x] Japanese/English runtime and AI guides synchronized.
- Remaining: normalization electron counts, initialization/material and parallel/timer inputs; concurrent/general different-material GS not implemented. DFT-MD global auto switch remains.
- Evidence: `docs/gs-mixing-config-validation-2026-11-16.md`. Molecular work stays paused.

### 2026-11-16: GS normalization inputs

- [x] Capture normalization nelec/nelec_spin in the owned GS config.
- [x] Preserve spin fallback and DC state-count precedence; reject nonpositive non-DC rho_dne divisor before activation.
- [x] Original equality across five metrics/14 cases; second-state divisor doubling across4 cases changes only its residual.
- [x] Capture normalization in retained trace and synchronize Japanese/English/AI instructions.
- Still shared: occupation/gap/material initialization inputs and parallel/timer controls. No general concurrent GS claim.
- Evidence: `docs/gs-normalization-config-validation-2026-11-16.md`.

## Recovery 2026-11-16
Occupation snapshot and tests reconstructed after temporary checkout loss; fresh persistent evidence in docs/gs-occupation-recovery-validation-2026-11-16.md. HSE shared-input work follows this checkpoint.

### HSE occupation routing (2026-11-16)

HSE nelec/temperature/nstate now resolve from each system's occupation snapshot. Native scalar full-support two-owner checks poison global values and compare actions/energies. Other shared controls and the outer TDDFT configuration guard remain. See docs/hse-occupation-inputs-validation-2026-11-16.md for measured scope; thermal/Wannier simultaneous owners are not established.

### TDDFT electron-count diagnostic ownership (2026-11-16)

- [x] subdip electron counts use system occupation input, retaining spin precedence, strict bounds and legacy fallback.
- [x] Four original-output comparisons, eight bounds, full TDDFT eight-case regression and packaged native probe pass.
- [x] Japanese/English/AI guides synchronized.
- Shared ionic/I/O/parallel controls and native STOP remain; outer TDDFT guards are not removed. Evidence: docs/dip-occupation-inputs-validation-2026-11-16.md.

## Mixed-material priority (2026-11-16)

The user explicitly targets LD/SBE/TDDFT in one multiscale run. See 2026-11-16-mixed-multiscale.md. LD+SBE endpoint parity, managed three-way TDDFT, native field-local LD and distributed-field mixed fixtures are now verified within the scopes below. Preserve the molecular pause and don't equate full reentrancy with mixed-backend integration.

- [x] Mixed managed-worker LD/SBE/TDDFT: real Si starter GS,4 ranks,100 steps, nonzero current at all three points and exact serial/MPI A/j parity. Evidence: docs/mixed-three-material-validation-2026-11-16.md. Native internally parallel TDDFT and field-local LD were subsequently verified below.

- [x] Native mixed fixture: field-local two-material LD + SBE + TDDFT internal orbital MPI (1/2 ranks), original serial replay, LD-only batch parity and LD vacuum contrast. Installed wheel runs4 ranks. Scope/evidence: docs/native-mixed-material-validation-2026-11-16.md. Multiple field ranks were subsequently verified below; scaling and mixed-session restart remain open.

- [x] Mixed FDTD field decomposition: x/y/z2-way field MPI with SBE1/TDDFT1,128 steps, exact A/E/j parity to field1. Sample exposes --field-grid and --compare. Simultaneous2field+2TDDFT ranks,scaling and full-session restart remain unverified. Evidence: docs/distributed-field-mixed-validation-2026-11-16.md.

### Mixed layer input entry (2026-11-16)

- [x] Reject longitudinal incident content in Multiscale1D's prefilled history, not only boundary samples; preserve valid transverse history.
- [x] Public native layer assembly from explicit material IDs/models, order, thickness, transverse counts and rank settings. Preserve actual thickness, point/material identities and backend provenance in the report.
- [x] Installed4-rank layout matches manual A/E/j exactly; installed2-rank SBE/LD needs no TDDFT files. Source and wheel442 tests pass (one skipped each).
- [x] Japanese/English coupling guides, AI task catalog/API inventory and packaged sample synchronized.
- Remaining: full mixed-session restart, memory estimation/cap, production scaling, simultaneous field2+TDDFT2 configuration within a larger resource budget, general native reentrancy. Layout planning does not automatically prepare materials or launch MPI; arbitrary user models enter as explicit bindings. Molecular work remains paused.
- Evidence: docs/native-layer-assembly-validation-2026-11-16.md; design:2026-11-16-native-layer-assembly.md.

- [x] Repeated-material layer validation: independently evolved Si×2 and SBE×2 in separate groups, four ranks each; distinct point IDs/drives/currents, exact independent serial replay for all three electronic points. Reproduction recipes added to both coupling guides and AI JSON. Evidence: docs/repeated-mixed-material-validation-2026-11-16.md. Full coupled-session restart remains the next acceptance milestone; native same-process reentrancy remains separate.

- [x] Sample-level replay restart: save at64, reconstruct field/local-LD/electronic histories in a new four-rank job, continue to128; full A/E/j exactly matches continuous execution. Inputs/library/GS/topology and saved histories are checked. Evidence: docs/native-replay-restart-validation-2026-11-16.md.
- [x] Direct native field/PML/LD/SBE/PZ-electronic internal state save/restore in the supplied mixed fixture;64 and63 step restart matches continuous128, with64/65 propagation calls. Generic user-backend archives/other XC and changed-topology restart remain open. See docs/native-direct-checkpoint-validation-2026-11-16.md.

### Non-PZ direct checkpoint (2026-11-16)

- [x] CheckpointV2 XC identity/derivative flags/input cval validation; restore scalar cval and rebuild Libxc handles from original configuration.
- [x] TB-mBJ automatic/fixed cval,Libxc-PZ,Libxc-PBE and PZ regression: per-XC converged Si SCF, continuous128 vs direct63+65; A/E/j bit-identical in all five cases. Three ranks(field1/SBE1/TDDFT1), no elapsed-time replay.
- [x] Source/wheel446 tests (one skipped each), packaging, source manifest update, English/Japanese/AI guides.
- Subsequent campaign below validates TB-mBJ/Libxc-PBE orbital MPI2. Fast-PR HSE native direct checkpoint, other XC, changed topology and generic backend persistence remain open. Molecular work stays paused. Evidence: docs/xc-direct-checkpoint-validation-2026-11-16.md.

### TB-mBJ / Libxc-PBE inner-MPI direct checkpoint (2026-11-16)

- [x] ea22fe15 delivered to independent SALMONpy branch.
- [x] Four outer ranks(field1/SBE1/TDDFT2), verified actual nproc_ob=2: TB-mBJ automatic/fixed cval and Libxc-PBE continuous128 versus direct save63/new-job restore65, all full A/E/j histories bit-identical.
- [x] Both electronic rank files, completion inventory and all hashes checked; clocks/lagged current and independent serial response pass.
- [x] Current audit reconciled: old delivery, incident validation, layer setup and sample direct restart are completed; broad GS ownership/native failure/resources/generic archive/environment work remains open.
- Evidence: docs/xc-inner-mpi-direct-checkpoint-validation-2026-11-16.md. No extra native implementation was needed. Field2+SBE1+TDDFT2 is still five ranks and outside budget. No unrestricted reentrancy or physical convergence claim.

### HSE engine constraint reaffirmed

Only fast-PR HSE is an eligible execution engine. Do not add conventional HSE or comparison baselines as runtime options. A future direct checkpoint must include the fast engine's ACE/factor-history/native resources and its dedicated propagator; adding hse06 to the semilocal middlepoint whitelist is insufficient.

### Checkpoint verification I/O (2026-11-16)

- [x] Partition payload SHA256 work across ranks, agree canonical manifest and all errors before restoration, retain shared-directory contract and inventory checks.
- [x] Installed Libxc-PBE innerMPI2 direct63→65 remains bit-identical to continuous128. Checksum-streamed bytes82,191,025 vs328,764,100 at4 ranks; no memory or wall-time guarantee.
- [x] Source/wheel450 tests (one skip), package and guide synchronization pass. Evidence: docs/checkpoint-verification-io-validation-2026-11-16.md.
- Fast-PR HSE direct persistence remains open; no conventional HSE engine is introduced.

### Fast-PR HSE direct internal persistence (2026-11-16)

- [x] Internal Fortran rt_checkpoint/rt_accepted_step, native numeric payload and HSE ACE/factor/PT-CN/source/action/history serialization; retain fresh FFT resources, never serialize opaque handles.
- [x] Thirty-two fresh-process cases:2/3k,ACE/learnedACE/full-action Taylor and learned PT-CN,split6/7,OMP1/2; currents,total energy and selected internal trace bit-identical, only6/5 remaining calls.
- [x] Truncated archive, changed settings and trailing-byte errors rejected; late-load failure invalidates stepping. Prior eight two-context interleaving/lifetime/configuration regressions pass.
- [x] Only fast-PR HSE objects used. No conventional HSE engine or comparison baseline engine added.
- Remaining: host-MPI/coupled HSE, public Python/C persistence API, broader symmetry/localized/distributed/GPU resources and physical convergence. Existing HSE-open notes above are historical; serial internal checkpoint is complete within this scope. Evidence: docs/fast-hse-direct-checkpoint-validation-2026-11-16.md.

### Fast-PR HSE outer MPI / FDTD direct restart (2026-11-16)

- [x] Explicitly audited fast-only worker builder and flushed STEP/SAVE/LOAD protocol; actual initial charge current query, no baseline subtraction, no host cwd changes or import-time compile/MPI init.
- [x] Public FastHSEMaterialProvider / FastHSEResponseBackend, one rank/group owning a serial child; immutable bundle/worker identities, independent work copies, timeout and cleanup.
- [x] Same-epoch native field/PML/LD, SBE, HSE worker internal state plus Python counters/current/A; generic worker checkpoint callback checks step/current/A and refuses unsupported many-workers-per-rank collisions.
- [x] Installed-wheel Taylor/learned PT-CN × field1/2, externalMPI3/4 ranks: save127/new-job remaining129 through256. A/E/j bit-identical to continuous; field1/2 also bit-identical. Real worker command logs contain exactly128..256 after LOAD127.
- [x] Source and installed459 tests (one skipped each), sdist-to-wheel, public API and English/Japanese/AI instructions synchronized.
- Remaining: HSE inner MPI/native-host integration, other topologies/resources, generic public checkpoint archives, recoverable native errors, memory caps and production/physical convergence. Outer-MPI/FDTD HSE connection is now complete within the serial-worker profile; older HSE-coupling-open notes are historical. Molecular work stays paused.
- Evidence: docs/fast-hse-mpi-coupled-checkpoint-validation-2026-11-16.md.

### Backlog check — 2026-11-17

- [x] Fast-PR HSE host-MPI/k-only direct continuation delivered in ace8a557. Historical HSE-inner-MPI-open notes above are superseded within this scope, not separate unfinished tasks.
- [x] Native layer planning admission: optional backend `supported_group_sizes` positive unique integer tuple, fast-HSE `(1,2)`. Reject unsupported/malformed declarations before MPI or material preparation. Undeclared custom backends retain existing execution/decomposition checks. This is not certification of their numerical or restart capabilities.
- [x] Direct ResponseSession admission and cross-group declaration agreement, shared with layer planning; A-packet wrapper forwards declarations.
- [x] Supplied mixed direct checkpoint: collectively compare persistence routes/modes, exact MPI ownership and field/electronic decomposition before LOAD; same-size different layouts rejected. Legacy manifest/1 requires original reader; arbitrary snapshot-only/replay backends are not admitted.
- [ ] Next: generic archive I/O/capability migrations and archive portability; retain custom backend freedom without assuming undeclared capabilities.
- [ ] Remaining primary work: native recoverable-error/handle contracts, remaining GS/TDDFT global ownership, generic checkpoint archives and changed topology, resource estimates/caps and scaling, platform/runtime distribution evidence. See current audit items5–10.
- [ ] Orbital/rgrid HSE is distinct from completed k-only HSE; no claim of support. Molecular work remains paused.

Validation for parallel admission: failing regression reproduced acceptance of an unsupported3-rank collective group before the fix; valid2-rank planning and malformed declarations now pass their checks without invoking provider.prepare. Source and installed-wheel suites each pass471 tests (one skip); sdist→wheel succeeds. Persistent artifacts: `/Users/otobetoshihito/.codex/outputs/salmonpy-parallel-admission-20261117`. No native kernel, update order or numerical formula changed; no new remote Actions run.

Direct-session admission validation: nine real three-rank collective rejections before callbacks, canonical tuple-order agreement and legacy undeclared trajectories verified in source/installed environments. Installed native FDTD+LD+SBE+HSE-k2 256-step outputs remain byte-identical. Source/wheel473 tests pass (one skip each). Evidence: `../response-admission-validation-2026-11-17.md`. This closes rank-admission agreement only; restart capability compatibility and broader native failure/ownership work remain open.

### Restart admission — 2026-11-17

- [x] Native-direct/2 sample manifests bind persistence formats/modes, exact group/point ownership and inner field/k/orbital/rgrid layouts; collective comparison precedes every LOAD. Fresh-candidate and within-group agreement checks included.
- [x] Require actual native payload filenames in inventory and reject collisions; unsupported save leaves no output directory. Old archives lacking compatibility metadata require their original reader.
- [x] Eight three-rank MPI preflight failures have zero LOAD calls. Final installed HSEk2 Taylor/PT-CN,Libxc-PBE orbital2 and field2+serialHSE direct resumes remain byte-identical to continuous with only129/65 new calls.
- [x] Source/wheel476 tests (one skip each),sdist→wheel and JP/EN/AI synchronization. Evidence: `../restart-preflight-validation-2026-11-17.md`.
- Remaining: generic public checkpoint contracts and portable/changing MPI layouts; full native error atomicity and state/input ownership; resource/cross-platform work. The sample's preflight does not complete these broader items.

### Common checkpoint admission API — 2026-11-17

- [x] Public immutable CheckpointCapability and serial/host-MPI check_checkpoint_compatibility; ordinary group mappings and explicit custom model identity, no SALMON context requirement.
- [x] Shared declaration,mode,saved-record,exact partition and optional required-file admission; earlier caller errors can be agreed with local_error. No file I/O, MPI initialization or backend callbacks in the API.
- [x] Existing native sample bridges lifecycle/ABI checks to the public validator, retaining direct-compatibility/1 record shape and legacy native-direct/2 archives.
- [x] Packaged user-model JSON persistence/independent continuation example and paired/AI guide entry. Generic archive writers,changed-topology redistribution and native rollback remain separate work.

Common-API evidence: `../common-checkpoint-api-validation-2026-11-17.md`. Source/wheel485 tests (one skip), public MPI3 ten rejection cases, existing coordinator eight rejection cases and a complete user-model17+15 restart pass. Four native continuation routes read unchanged earlier native-direct/2 archives and remain byte-identical to continuous execution, with129/65 new calls only. Public declaration/admission is complete within this scope; generic I/O,redistribution,resource controls and native rollback remain open.

### Common checkpoint execution — 2026-11-17

- [x] CheckpointTarget binds backend-owned zero-argument save/restore callbacks; execute_checkpoint admits before callbacks and returns opaque results locally without MPI serialization.
- [x] Synchronized group-local component slots support unequal point counts; callback failure stops later slots collectively. Completed changes are not rolled back; backend collective progress remains its responsibility.
- [x] Custom JSON example and native field/electronic/snapshot routes use the shared entry; completion-marker publication errors are also agreed collectively in the native coordinator.
- Remaining: portable generic archive layouts,redistribution,aggregate resource limits and native failure recovery/ownership. The common dispatcher prescribes no payload encoding.

Execution evidence: `../checkpoint-execution-validation-2026-11-17.md`. Source/wheel489 tests (one skip),MPI unequal-group/opaque-result/failure-slot tests and collective completion-publication failure pass. New common-entry saves and direct resumes across HSEk2 Taylor/PT-CN,PBEorbital2 and field2/serialHSE remain byte-identical to continuous with129/65 new calls. Pre-existing native-direct/2 archives remain readable. Generic callback execution is complete within this contract; archive interoperability/migration,redistribution and native rollback remain separate.

### Protected candidate restore — 2026-11-17

- [x] Common `restore_checkpoint` and `CheckpointState`: independent candidates, collective readiness vote, single local bundle publication; no destructive fallback.
- [x] Late-failure cleanup preserves originals; cleanup failures reported separately. Successful replacement returns previous objects for explicit host retirement.
- [x] Real four-rank MPI/custom backend failure injection and successful retry, actual LD/SBE future-result equality, source/wheel495 tests (one skip each), sdist→wheel, JP/EN/AI synchronization.
- [x] Experimental native numerical candidate ABI added below; old direct-only backends remain rejected. Current ordinary native coordinator is still explicitly non-atomic.

Evidence: `../atomic-restore-validation-2026-11-17.md`. Backend isolation and healthy collective progress are required; process crashes, hidden shared mutable state and topology redistribution are outside this guarantee.

### Native numerical restore candidates — 2026-11-17

- [x] Standard native FDTD/TDDFT and fast-PR HSE stage detached numerical payloads; publish after collective success without checkpoint re-read, or discard while preserving originals.
- [x] Group mode/readiness agreement and HSE context ownership; older binaries rejected by candidate protocol query. One pending candidate per native lifetime.
- [x] Real four-rank feedback tests with field2, TDDFT2 and HSE2 layouts: late HSE shard failure preserves every original shard byte-for-byte; originals continue, retry16-step A/E/j agrees exactly with continuous. Taylor and installed-wheel PT-CN agree with previous native binaries; legacy HSE immediate save/resume also passes.
- [x] Paired/AI documentation distinguishes numerical candidates and explicit native publication from the generic holder and the unchanged ordinary non-atomic coordinator.
- [x] Ordinary mixed-native coordinator now defaults to prepared native/snapshot candidates and Python/session metadata publication after global agreement; details below.
- [ ] Stable public native ABI, allocation-failure recovery, unrestricted multi-context native runtime and changed-topology restart remain open.

Evidence: `../native-candidate-restore-validation-2026-11-17.md`. Cache validity flags are deliberately reconstructed; do not claim full post-success electronic archive byte identity from current/field equality.

### Ordinary mixed-native protected restoration — 2026-11-17

- [x] Default direct resume stages field/PML/all field-local LD, SBE handles, native TDDFT/fast-PR HSE and Python arrays/clocks/adapters/ownership before publication. Advanced paused sessions can restore and retry.
- [x] Handled staging failure preserves original numerical bytes and live references/flags; originals remain usable. Snapshot retirement occurs after success with separate cleanup diagnostics.
- [x] Restore/protected mode agreement and history parsing/validation precede LOAD. Older ABIs/direct-only backends reject without fallback; explicit unprotected flag retains the legacy route.
- [x] Four-rank source/wheel tests cover field2, TDDFT2 and HSE2 (Taylor/PT-CN) with two field-local LD materials plus SBE; late rank failure preserves all originals, three subsequent frames after retry match continuous exactly. A fresh ordinary installed CLI job also resumes64→68 and matches continuous65–67 exactly.
- [x] Source/wheel502 tests (one skip each),sdist→wheel,eight existing MPI preflight rejections and JP/EN/AI synchronization. No remote Actions.
- [ ] Stable public native ABI,changed MPI topology,arbitrary backend cloning and native allocation/process-failure recovery remain separate work; this does not finish unrestricted multi-context ownership.

Evidence: `../mixed-protected-restore-validation-2026-11-17.md`.

### Checklist reconciliation after fe286bc1 — 2026-11-17

- [x] Corrected the current audit's stale cross-group restore/open statements; ordinary protected restore is complete for the measured supported paths.
- [x] Distinguished completed common checkpoint API and GS/2 subsets from broader archive/identity/MPI portability work. M4–M6 remain umbrella milestones with completed subtasks, not wholly missing features.
- [x] Kept molecular work paused and preserved historical completion records. Ordinary native output-file truncation, process/OOM recovery and general step failure atomicity were not marked complete.
- Suggested next practical priority: aggregate memory estimation and budget admission, including candidate coexistence during protected restore. Stable native context/ABI ownership and supported-platform runtime evidence remain parallel broader objectives. This reconciliation adds no feature or new numerical test claim.

### Advisory memory budgets — 2026-11-17

- [x] Public MemoryEstimate/MemoryReport/memory_advice and lazy peak-RSS probe; custom backends supply bytes/basis/unknowns. Host-MPI sums across the supplied communicator; missing values are not zero.
- [x] Ordinary mixed CLI --memory-budget-gib: pre-initialization heuristic, initialized/segment-end RSS, protected restore original/candidate/read/copy/JSON/snapshot allowances plus reserve. Budgets are runtime advice, not physical input.
- [x] User changed the action from stopping to warning only. Excess never rejects computation or restore; measurement, warning-filter and report-output failures stay diagnostic. OS/scheduler hard caps are not installed.
- [x] Source/wheel509 tests (one skip each),sdist→wheel,MPI4 public arithmetic/unknowns/warning-filter cases and actual4GiB/4-core mixed restore. Tiny advisory budget plus a failed rank probe still completes and retries exactly. JP/EN/AI synchronized; no remote Actions.
- [ ] Broader backend-specific calibration, children/unsupported-platform measurement, aggregate production scaling and capacity predictions remain open. No universal heap/RSS upper bound or physical-convergence guarantee is claimed.

Evidence: `../memory-advice-validation-2026-11-17.md`. The earlier recommendation for budget admission is superseded by the user's explicit warning-only decision.

### 2026-11-17 ordinary mixed-session identity

- [x] M6d: ordinary native mixed sample carries run, logical point, material, definition and initial-state identity plus MPI ownership and LD cell mapping into direct/replay checkpoints and normal output archives. Explicit configured LD/electronic records are retained.
- [x] Reject identity mismatch before native load; agree identity across ranks before checkpoint writes. Old unidentified archives require explicit legacy identity mode; protected restore remains default.
- [x] Four-rank TDDFT and fast-PR HSE k-only two-rank tests preserve all original state after injected failure and resume exactly without replay; fresh installed-wheel TDDFT continuation also matches exactly.

Evidence: `../mixed-identity-validation-2026-11-17.md`. M6 remains an umbrella milestone: arbitrary backend formats, branch semantics and general in-memory identity remain open. Memory budgets remain advisory warnings.

### 2026-11-17 point snapshot identity

- [x] M6e: PointState process-local snapshots wrap backend tokens with exact run/point/material/initial identity, schema and version. Reject mismatches and closed targets before callbacks. No arbitrary token conversion, implicit cloning or branch semantics.
- [x] Real LD/SBE same-point future currents match exactly after restore; cross-point rejection leaves the other state unchanged. Matching IDs do not bypass backend dt/model checks. Custom opaque payloads remain backend-owned.
- [x] Updated point_ld sample and JP/EN/AI guides. Source/wheel517 tests (one skip), sdist→wheel and installed example pass. No solver or MPI code changed; no remote Actions.

Evidence: `../point-snapshot-identity-validation-2026-11-17.md`. General branch APIs, arbitrary disk payloads, changed topology and native ownership/recoverable errors remain open. Molecular work stays paused.

### 2026-11-17 recoverable standard TDDFT admission

- [x] Standard embedded TDDFT protocol1 rejects statuses10..14 collectively on the attached communicator before initialization latch or numerical update. Invalid/nonfinite/stale/rank-inconsistent inputs, lifecycle/capacity and pending-candidate checks; scalar IEEE checks only.
- [x] Python initialization retry only for declared unanimous pre-setup10/11; step InputError only for declared unanimous pre-update rejection. Old/unknown binaries retain conservative RuntimeFailure; successful native initialization remains one lifetime.
- [x] Two independent2-rank groups, raw ABI13 negative cases with numerical checkpoint byte preservation, Python init retry and valid8-step old/new current parity. Four-rank ordinary mixed64-step old/new A/E/j exact, protected-failure/retry regression exact. Source/wheel519 tests (one skip),sdist→wheel and JP/EN/AI sync; no Actions.
- [ ] General native failure contract remains OPEN: deep initialization/step STOP, nonfinite results after mutation, allocations, native FDTD/fast-HSE equivalents, whole-session step rollback and process/MPI failures are not solved by admission.

Evidence: `../native-tddft-admission-validation-2026-11-17.md`. The current audit item5 remains partial; do not mark general native failure recovery complete.

### 2026-11-17 expanded Mac resource campaign

- [x] User supplied host capacity64GB/18 cores. Preserve4GiB/4-core beginner profile. Campaign uses48GiB advisory budget, at most15 ranks, OMP/BLAS1 and sequential jobs.
- [x] Sample --field-shape: even counts>=40, default40³ and centered physical point positions preserved. Field dimensions now explicit in comparison/result metadata; loop timing excludes setup/output/serial replay. Electronic GS grid/k sampling unchanged.
- [x] Previously unverified5-rank field2+SBE1+TDDFT2 protected restore passes; native fast-PR HSE k2 also passes with field2 and SBE1.96³(13.824×40³ cells) field1/2/4/8/12, total4/5/7/11/15 ranks, two repetitions each: all128-step A/E/j arrays exactly equal to baseline.
- [x]96³/15-rank protected late-failure rollback preserves all original native bytes and Python states; original continues, retry and3 future frames exactly match uninterrupted execution without replay. Normal summed rank peaks<=1.83GiB; failure-injection path about5.1GiB, not simultaneous RSS or a universal bound.
- [x] Source/wheel524 tests (one skip),sdist→wheel, sample/controller and JP/EN/AI guides. No native numerical source changes, no remote Actions.
- [ ] Larger electronic grids/k meshes, many material points, scaling beyond this local fixture and cross-platform/runtime delivery remain open. Native post-update error recovery is unchanged.

Evidence: `../expanded-mixed-scaling-validation-2026-11-17.md`.

### 2026-11-17 standard TDDFT deep numerical recovery

- [x] Native per-step in-memory numerical backup and state-owned stage status. Collective finite checks after propagation/density/Hartree/XC/endpoint; predictor-corrector errors return collectively with balanced timers. Protected invalid-propagator error returns status instead of STOP. Production code owns no new mutable recovery globals.
- [x] Rollback copies original numerical payload/XC state and reconstructs PP phases without replay. Preserve original overlap flags, unlike external disk LOAD normalization. Failed attempts write no ordinary RT records. Unsupported protected profiles reject before mutation.
- [x] Real overflow and separately compiled one-rank fault stages20..24; PZ, automatic/fixed TB-mBJ,Libxc-PZ/PBE internal orbital MPI2, plus PBE k/rgrid MPI2. Every rejected numerical checkpoint is byte-identical; retry16-step currents and accepted energy/current records match the corresponding old binary/decomposition exactly. Odd/even buffers and two independent MPI groups also pass.
- [x] Installed production/fault paths and 128-step normal5-rank FDTD2+LD+SBE1+TDDFT2 A/E/j/serial replay exact to old binary. Existing5-rank protected checkpoint/failure/retry still passes.
- [x] Source/installed525 tests (one skip), sdist→wheel and JP/EN/AI sync. No remote Actions. Test-only injection symbol absent in production.
- [ ] This closes the measured TDDFT-group numerical failure subset, not whole mixed-step rollback, allocation/MPI/process failures, third-party STOP, deep initialization recovery, HSE/FDTD equivalents, stable native handles or general reentrancy. Per-step copying/checking has resource overhead; budgets remain warning-only.

Evidence: `../tddft-deep-recovery-validation-2026-11-17.md`.

### 2026-11-17 fast-PR HSE deep numerical recovery

- [x] Separate HSE recovery protocol1; per-step electronic/XC and ACE/history/cache/counter numerical backup, with no copied FFT/MPI owners. Handled collective finite/operator/PT-CN/history errors restore accepted state without replay. Failed attempts append no ordinary RT records.
- [x] Internal k-only MPI2 Taylor/PT-CN, first/warm exact/warm prediction/final-step faults, genuine overflow and PT-CN nonlinear rejection, independent two-group MPI4. Payload bytes preserved; twelve32-step current/energy trajectories and accepted output match the old fixed binary exactly.
- [x] Five-rank FDTD2+LD+SBE1+HSE2 normal128-step A/E/j exact to old; protected late-failure restore/continuation/retry regression passes without replay. Raw exhausted-capacity ABI rejection tested after a reviewer-found defect was fixed.
- [x] Source/installed526 tests (one skip),sdist→wheel/hash/layout checks, JP/EN and AI JSON sync. No remote Actions.
- [ ] General initialization/third-party STOP, allocation/process/MPI failures, native FDTD deep failure and whole mixed-step rollback remain open. k-only HSE recovery does not establish arbitrary native reentrancy or physical convergence.

Evidence: `../fast-hse-deep-recovery-validation-2026-11-17.md`. Earlier entries describe their historical scope; the supported HSE-group subset is now closed, not the general failure contract.

### 2026-11-17 native E/H FDTD numerical recovery

- [x] Field recovery protocol1 and checked endpoint ABI. Per-step numerical field/PML/all-local-LD/current/clock backup; collective E/LD/H/sample/integrated-A/ordinary-observation gates precede scientific output. Handled errors restore without propagation replay and preserve Python A/E/step/failed.
- [x] MPI1/2 and independent2x2 groups, first/warm/final one-rank faults, actual finite-input overflow, raw mixed-entry/capacity rejection and three-step observation stride. Rejected bytes/output unchanged, valid retry exact to old accepted32-step A/E/H and observation records.
- [x] Five-rank Libxc-PBE TDDFT and fast-PR HSE mixed128-step A/E/j exact to old; both protected restore/failure/continuation/retry regressions still pass. Source/installed528 tests (one skip), sdist→wheel/hash/layout and JP/EN/AI sync. No Actions.
- [ ] General initialization/TTM/LR/ASE/ART/frequency-observation paths, allocation/process/MPI/uncontrolled STOP/output-device failures and whole mixed-step rollback remain open. Field-group recovery does not complete general native failure atomicity or multi-context ownership.

Evidence: `../fdtd-deep-recovery-validation-2026-11-17.md`. Current audit item5 now includes scoped standard TDDFT, fast-PR HSE and native field numerical recovery; broader work remains open.

### 2026-11-17 common explicit protected-step entry

- [x] Public StepTarget/guarded_step/StepRecoveryError with opaque local tokens, group key/release-policy agreement, unequal target counts, partial-capture release(None) and reverse rollback/cleanup. Unsafe rollback never claims preservation; cleanup failure after success reports possible acceptance.
- [x] CoupledStepper.step_protected requires all local field/material owners and additionally restores coupling A/J/clocks/adapters/measurement metadata. Incomplete rollback/cleanup invalidates the response and field. Ordinary step unchanged; no automatic native target inference.
- [x] Six-rank Multiscale1D + real LD/SBE + collective custom two-point group: source/installed first/warm/final custom failures, LD and post-field mutations preserve original signatures; retry80 frames exactly match ordinary execution without replay. Missing owner/release policy/partial capture and rollback/release failure negatives pass.
- [x] Source/installed534 tests (one skip),sdist→wheel/hash/layout, installed example/manual search and JP/EN/AI sync. No solver changes or remote Actions.
- [ ] Native aggregate numerical/output transactions remain open. Individual E/H/TDDFT/fast-PR HSE recovery and direct checkpoint protocols do not establish whole-step atomicity. User callbacks must explicitly cover all effects; MPI/process/native STOP/deadlock/OOM remain outside this contract.

Evidence: `../protected-step-validation-2026-11-17.md`. This closes a common-API/light/custom subset, not the broader native failure item5.

### 2026-11-17 native aggregate protected steps

- [x] Separate explicit native whole-step protocol and native_step_targets connect native E/H/PML/all local LD, SBE, standard TDDFT, fast-PR HSE and Multiscale1D to the common protected-step API. Numerical/XC/ACE/history backups exclude FFT/MPI owners and release all allocated snapshot arrays.
- [x] Ordinary private RT/point-output rollback closes/truncates/reopens audited writers; checkpoint staging and step backups exclude each other. Unknown/failed lifetimes fail closed. Other output profiles reject before advancing; no hidden fallback or compilation.
- [x] Ten final 32-step MPI cases cover PZ/TB-mBJ/Libxc, fast-PR Taylor/PT-CN, separate PBE+HSE seven-rank mixtures and two-point 1D. First/warm/final late/deep rejection preserves accepted numerical bytes and output; explicit retries match uninterrupted arrays and every ordinary stream exactly. Five unchanged-library comparisons also match exactly.
- [x] MPI ownership/admission/unsafe-lifetime negatives, seven-rank protected checkpoint-restore regression, source/installed545 tests (one skip each), source-distribution-to-wheel/hash validation, installed seven-rank/1D fault-retry campaigns and paired JP/EN/AI guides. No Actions.
- [ ] Broader audit item5 remains open for initialization/unsupported profiles, I/O/allocation/process/MPI failures and unreturned STOP calls. This is not durable or observer-invisible file publication; private output ownership is required. General native reentrancy and physical convergence are separate work.

Evidence: `../native-protected-step-validation-2026-11-17.md`. Earlier open aggregate-step notes are superseded only within this measured scope. Memory budgets remain warning-only.

### 2026-11-17 enlarged protected mixed-step measurements

- [x] Existing installed wheel/production libraries:96³ Maxwell + two local LD materials + SBE1 + Libxc-PBE2 + separate fast-PR HSE PT-CN2. Field2/4,total7/9 ranks,64 accepted steps,two repetitions each with/without whole-step protection. All eight arrays and every ordinary stream exactly equal across repetitions/decompositions/protection modes.
- [x] Local protected/plain mean loop ratios1.198(field2)/1.082(field4); protected maximum summed rank peak RSS2.692/4.417GiB. Peaks are not simultaneous totals or heap bounds. Unknown resource allowances remain unknown;48GiB budget advisory only,OMP/BLAS1,no Actions.
- [x]96³ seven-rank step13 late-electronic and actual field-overflow rejection preserve full numerical payload hashes,Python references/clocks and ordinary output bytes. Explicit retry64-step arrays/streams equal uninterrupted execution; zero elapsed-history replay. Diagnostic summed rank peaks3.943GiB.
- [x] Paired JP/EN resource notes and AI JSON synchronized;11 manual tests pass. No solver/API/packaging changes or rebuild required. Recipes,commands,measurements,arrays/streams and diagnostic harness retained persistently.
- [ ] Larger electronic grids/k meshes, many material points/LD cells, transferable backend calibration, cross-platform evidence and output/initialization failure contracts remain open. These measurements do not prescribe production rank choices or certify physical convergence. Molecular work remains paused.

Evidence: `../protected-mixed-scaling-validation-2026-11-17.md`. Audit item7 is reduced only by this measured subset; no universal memory model is claimed.

### 2026-11-17 repeated material definitions and independent points

- [x] New packaged64³ multipoint layout:four field ranks,SBE4 independent points,Libxc-PBE2 points with two orbital ranks each,fast-PR HSE PT-CN2 points with two k ranks each,total16 ranks. Drude/Lorentz1152 cells each remain field-local.
- [x]64-step uninterrupted/protected arrays and every ordinary output stream exact; independent serial response checks pass at all8 points. Same-material points have distinct drives/currents. Direct31→64 restart advances only33 native steps,all full-history arrays exact; logical/material/LD identities preserved.
- [x] Warm13 failure in just one HSE point group and actual field overflow restore all groups' original payload hashes/Python states/outputs; explicit retry exact without elapsed-history replay. Source/installed545 tests(one skip),sdist→wheel/native hashes and paired JP/EN/AI synchronization pass. Packaged-layout installed restart verified. No solver/API changes or Actions.
- [ ] Many more points, larger electronic grids/k meshes, transferable capacity models and arbitrary in-process native multi-context ownership remain open. Existing one-point-per-group architecture and same-topology restart remain; molecular work stays paused. Budget excess only warns.

Evidence: `../multipoint-mixed-validation-2026-11-17.md`. Audit item7 now includes this8-point/2304-LD-cell subset, not general production scaling.

### 2026-11-17 ion initialization ownership

- [x] Standard and fast-PR full-source value-owned ion input (atom/element counts,Cartesian/reduced positions,species). Explicit initializer validates allocation/shape/lower bounds/species/finiteness before changing state and owns ion arrays. No global reads/writes or resource-owner copies in explicit helper; legacy wrapper retains downstream Rion normalization.
- [x] Two different captured configurations,poisoned globals,alternating initialized systems and malformed inputs:standard/HSE native probes OMP1/2 pass. Lower-bound review finding fixed with rejection regressions. No arbitrary full-GS different-input independence claim.
- [x] MPI2 GS fresh three-iteration energy/eigen/k/force outputs exact to old;32-step five-rank PBE/HSE Taylor/PT-CN and seven-rank combined protected step/fault/retry results exact. Seven-rank direct restore/failure/continuation/retry and identity negatives pass without replay.
- [x] Manifests/provenance refreshed,new external builds,source/installed545 tests(one skip each),source-distribution→wheel/native hashes,installed seven-rank mixed/five-rank1D protected retries and paired JP/EN/AI guides. No Actions or shared-kernel changes.
- [ ] Remaining cell/k/symmetry/communicator/ionic dynamics/I/O/timer/DC/SOC/+U ownership,stable handles and general initialization/allocation-failure recovery remain open. Legacy wrapper still normalizes shared Rion for downstream code. This only closes the captured ion-input/owned-array slice of audit item6.

Evidence: `../ion-initialization-validation-2026-11-17.md`.

### 2026-11-17 standard cell/grid initialization ownership

- [x] Port existing fast-PR cell/grid value types and explicit resolvers to standard solver/optional RT-step structures. Both profiles have capture-only cell_config_from_input. Existing derived-type layouts,valid cell precedence/arithmetic,grid epsilon/parity and compatibility al normalization preserved.
- [x] Standard/HSE cell/grid/ion probes OMP1/2,global-poisoning and frozen periodic/isolated/Maxwell grid references pass;invalid input preserves original geometry/indices. Standard/RT-step/HSE source preparation and manifest chains verified;shared kernels unchanged. Optional RT-step full fresh4-worker build and matching-profile cell/grid/ion OMP1/2 probes pass;no full RT-profile dynamics claim.
- [x] MPI2 fresh GS energy/eigen/k/force exact;five-rank PBE/HSE Taylor/PT-CN and seven-rank mixed protected failure/retry exact. Seven-rank direct restore failure/original continuation/retry pass without replay. Source/installed545 tests(one skip each),sdist→wheel/native hashes and JP/EN/AI guides synchronized. No Actions.
- [ ] Compatibility wrappers still consume shared inputs;whole grid coordinates,parallel ownership,k/symmetry/FFT resources and full heterogeneous GS/TDDFT initialization remain open. Explicit resolver ownership does not establish arbitrary native reentrancy or stable handles. Molecular work paused;budget warning-only.

Evidence: `../cell-grid-initialization-validation-2026-11-17.md`.

### 2026-11-17 initial orbital layout ownership

- [x] Standard/fast-PR capture mode,periodicity,spin/SOC,state/electron counts,temperature and kmesh/file sentinel into explicit value input. Shared selector resolves nspin/no/real storage without globals. HSE uses owned occupation/effective kmesh/legacy hybrid predicate;standard initial occupation consumes captured electron/SOC values. Existing type layouts/public API unchanged.
- [x] Native GS Gamma/non-Gamma/file/nonorthogonal/isolated/SOC/hybrid and RT even/odd/SOC/spin-resolved/temperature-threshold cases,global poisoning,alternating results and rejection preservation pass standard/HSE/RT OMP1/2. Cell/grid/ion probes also pass. Strict manifests/preparation,new external matching-profile builds and shared-kernel equality verified.
- [x] MPI2 fresh GS energy/eigen/k/force exact;five-rank PBE/HSE Taylor/PT-CN and seven-rank mixed failure/retry exact. Direct restore failure/original continuation/retry exact without replay. Source/installed545 tests(one skip each),sdist→wheel/native hashes,installed seven-rank mixed/five-rank1D fault-retry and JP/EN/AI updated. No Actions.
- [ ] Full occupation-capacity/SCF-input ownership,remaining k/symmetry/communicator dependencies and arbitrary heterogeneous native GS/TDDFT initialization remain open. This closes only the initial orbital-selection slice;no new native initialization-failure/handle guarantee. Molecular work paused,budgets advisory only.

Evidence: `../orbital-initialization-validation-2026-11-17.md`.

### 2026-11-17 fast-PR HSE symmetry primary path

- [x] Project lead: prioritize fast-PR HSE and symmetry. Standard k-point port deferred before tracked solver changes; distinct standard/fast-PR arithmetic preserved. No conventional HSE engine.
- [x] Explicit symmetry EXX owners bounded to CPU scalar periodic full-q/internal-k1/factor-off/Taylor-ACE;test-only full-action comparison. Independent symmetry maps/phases/FFT resources;owned full-grid count/per-axis admission before refresh/action. Automatic GS binding unchanged;learned/PT-CN symmetry/internal-k2 still reject.
- [x]16-case C2z/mirror-x4-rep vs independently enumerated full8-k fast kernel comparison: full/ACE,star/orbits,OMP1/2,MPI1/2 processes each COMM_SELF. Different per-process settings/call counts; finite values,interleaving,source update,grid/layout rejection,close/reopen. Max action difference~2.80e-16,energy~1.12e-16;old full-k reference bytes exact.
- [x] Existing NON-reduced GS/mixed/protected retry/direct restore remain exact;source/installed545 tests(one skip each),sdist→wheel/native/object audits,packaged MPI2-process/OMP2 probe,installed mixed/1D fault-retry and paired JP/EN/AI scope guides. No Actions,budgets warning-only.
- [ ] Next priority: owned actual-A little-group checks,complete symmetry-map/history/cache direct/protected persistence,then reduced internal-k MPI and mixed fault/resume. Current TDDFT/FDTD/checkpoint symmetry prohibitions remain. Core proof does not establish symmetric TDDFT or physical convergence. Molecular work paused.

Evidence: `../fast-hse-symmetry-validation-2026-11-17.md`.

### 2026-11-18 owned fast-PR uniform-field admission

- [x] Read-only state-owned Cartesian little-group check for uniform A/c; normalized nonzero direction/tolerance1e-10, finite/readiness/geometry/operation identity checks and status0/1/2. No shared symmetry fallback. TDDFT local pre-backup/propagation hook with existing collective admission.
- [x] Native independent C2z/mirror-x/rotated-cell, zero/tiny/huge/nonfinite fields, tolerance boundary, malformed/stale/closed owners and preservation; existing16 exchange cases and non-reduced GS/mixed/protected retry/direct restore remain exact.
- [ ] Complete symmetry-map/history/cache persistence; reduced internal-k MPI and symmetric mixed failure/retry/resume. Symmetric TDDFT/FDTD/checkpoints remain prohibited. No conventional HSE engine; molecular work paused; budgets advisory only.

Evidence: `../owned-hse-field-validation-2026-11-18.md`.

### 2026-11-18 fast-PR symmetry exchange checkpoint

- [x] Internal exchange-core V2 saves value-only map/cache and existing numerical history. Fresh-owner static layout/identity match, bounded cache load, staged/validate-only reads and value-array moves; FFT/C/MPI ownership untouched. Non-reduced V1 bytes unchanged.
- [x]32 full/ACE/star/orbits/cache0-1/OMP1-2/MPI1-2 COMM_SELF-process cases; independent C2z/mirror-x, reopen, exact action/energy/future checkpoint bytes, truncation/group/partition rejection preservation; cache-range/NaN guards. Existing non-reduced mixed/protected retry/direct restore exact.
- [ ] Symmetric TDDFT state integration and coupled direct/protected resume; reduced internal-k MPI and symmetric mixed failure/retry. Existing symmetric TDDFT open prohibition stays. No claim of learned/PT-CN symmetry or physical convergence; molecular work paused; budgets advisory.

Evidence: `../symmetry-checkpoint-validation-2026-11-18.md`.

### 2026-11-18 internal fast-PR symmetry TDDFT

- [x] Explicit internal Fortran rt_open allow_symmetry opt-in for scalar CPU/full-q/internal-k1/factor-off Taylor, one initialization/process. Ordinary worker/C/Python defaults stay denied. Remaining process-common symmetry flags/operations/lattice/inverses guarded before step/current/checkpoint/staged commit.
- [x] Full8/reduced4-k C2z/mirror-x synthetic TDDFT; ACE/full-action, MPI1OMP1/MPI2OMP2 independent COMM_SELF processes with different screening/12-13 step schedules. Odd/even direct/staged resume and final numerical checkpoints exact, zero accepted-step replay; corrected-A retry exact. Drift/truncation/trailing/default/learned/PT-CN/internal-k2/second-lifetime rejection verified. Symmetry archive overlap flags canonicalized only in copies; live flags/non-reduced V1 unchanged.
- [ ] Public native-C/Python opt-in and field-admission retry contract, protected symmetric mixed/FDTD testing; reduced internal-k MPI partitions/failure/restart. No claim of arbitrary heterogeneous same-process initialization, learned/PT-CN symmetry or physical convergence. Molecular work paused; budgets advisory.

Evidence: `../symmetry-tddft-validation-2026-11-18.md`.


### 2026-11-18 public symmetry and protected mixed connection

- [x] Backward-compatible C default init plus explicit symmetry init/query1. NativeFastHSEBackend allow_symmetry=True is lazy, identity-bound and group1-only; old/missing protocols and inappropriate inputs reject before lifetime. Declared unanimous finite34 yields state-preserving InputError; corrected A retry verified in independent2 groups and fresh direct restart.
- [x] Ordinary sample --hse-symmetry/--hse-reference-bundle, no arbitrary GS conversion or field projection. Synthetic phase.23 full/reduced inputs: field1/2 + two local LD + SBE1 + PBE-TDDFT2 + HSE1, five/six ranks; before/after late/deep/field-admission failure bytes and ordinary streams preserved. Protected restore save16/fail18 preserves all originals, retry3 future frames exact; fresh CLI16→32 continuation equals continuous, no accepted replay.
- [ ] NEW numerical issue: Libxc5.2.3 HSE near-zero-gradient derivatives are sensitive; phase0 failed full/reduced current comparison~5.7e-10, not fixed by dt refinement. Pure Libxc diagnostic reproduces finite but large derivative variation. No smoothing/clamping/tolerance change/library patch applied. Phase.23 is an explicit integration fixture, not a production remedy. Resolve with a mathematically consistent verified treatment before claiming robust physical accuracy.
- [ ] Reduced internal-k MPI2, learned/PT-CN symmetry, broader initialization/physical convergence remain open. Default remains denied and molecular work paused. Resource budgets advisory only; no Actions.

Evidence: `../public-hse-symmetry-validation-2026-11-18.md`.


## 2026-11-18 — Libxc HSE small-gradient cancellation

- [x] Stable same-formula WPBEh(0) energy/derivatives in squared reduced gradient;
  screened quadratic limit below roundoff, untouched PBE correlation.
- [x] Native regression fails on the old library;240-case density/screening sweep,
  independent80-digit unscreened/screened energy differences, polarized/core
  derivatives and low-density screening boundary checks pass.
- [x] Original phase0 full8/reduced4-k weak-field regression, OMP1/2; FDTD
  coupling, failure/rollback, protected restore and fresh-process direct resume.
- [x] Ordinary PBE byte equality, fast-PR Taylor/PT-CN/mixed roundoff agreement.
- [ ] Reduced internal-k MPI, broader initialization and converged physical benchmarks.

Historical unresolved sensitivity notes above describe the old binaries. Both
serial full-k reference and MPI runtime need explicit rebuilding from new sources.
Raw Libxc diagnostic sensitivity remains expected; no library binary is patched.


## 2026-11-18 — Reduced-k symmetry internal MPI2

- [x] Explicit one/two k-rank profile; separate capability rejects old libraries.
- [x] Rank-local star/phase cache checkpoint bounds and authenticated same-partition restore.
- [x] MPI consensus preserves recoverable status34; rejected field leaves state unchanged.
- [x] Phase0 full8/reduced4-k comparison; two independent2-rank save7→32 resume traces/final bytes.
- [x] Ordinary mixed sample supports --hse-ranks2; same-decomposition direct resume.
- [ ] Broader partitions/initialization and converged physical validation; checkpoint repartitioning is not supported.


## 2026-11-18 — Symmetry initialization admission

- [x] Read-only canonical sym.dat parsing and finite/retained-group checks before native lifetime/work copy.
- [x] Disabled-axis filtering and fast-PR signed-permutation/identity/duplicate/closure rules; Fortran D/repeat syntax.
- [x] One malformed rank per2-rank group: collective rejection, no lifetime, repaired provider starts normally; three fault kinds.
- [x] Normal current/energy traces retained; own save/restart still direct without temporal replay.
- [ ] Cell isometry, atomic/species/grid/k-star/restart and general Fortran initialization failures remain native hardening work.

### Next actionable backlog (historical unchecked notes above are superseded)

1. Broader initialization admission/recovery for selected native assets; no general STOP recovery yet.
2. Learned/PT-CN symmetry has not been admitted by the factor-off Taylor work.
3. Broader ownership/handles (different-input GS/TDDFT, MPI/I/O/timers and specialized resources).
4. Additional environments, calibrated resource/scaling models and interoperable archives.

Physical convergence remains material/user-owned; molecular work stays paused.
Larger MPI partitions and checkpoint repartitioning are not implied by MPI2 results.


## 2026-11-18 — Cell/atom/grid initialization admission

- [x] Effective cell and grid follow native priorities/units/rounding and integer bounds; finite volume/inverse; current native profile remains axis-aligned/au-only.
- [x] Explicit Cartesian/reduced inline/file atoms, counts/species/finite values; case-preserving bundle-local paths, no new dependency.
- [x] Symmetry cell metric, bounded affine full-grid checks and unique same-species atomic permutations.
- [x] MPI2 groups reject one bad rank before lifetime/work copy, repair/retry, unchanged normal current/energy and own direct resume.
- [ ] General-cell native initialization (axis-permuted fixtures hit ACE init STOP); k-star, pseudopotential-content/restart-identity and general native recovery remain open.

This supersedes the cell/atom/grid subset of older unchecked initialization notes.
GS/TDDFT ownership, learned/PT-CN symmetry, environment/scaling and portable archive work remain separate.

## 2026-11-18 — Warning-first restriction audit

- [x] Layer rank budget is advisory; all points and required ranks are preserved.
  Actual communicator-size/decomposition validation remains an execution requirement.
- [x] Finite GS occupation range/electron-sum discrepancies warn without mutation.
- [x] Removed default full-orbital finiteness rescan from state creation; size checks remain.
- [x] Memory budgets were already advisory. JP/EN coupled guides and AI JSON synchronized.
- [x] Six asset admission faults on two independent MPI2 groups: collective rejection,
  repair/retry and bit-identical current/energy; fresh-process own direct resume.
  Mixed LD/SBE/field/reduced-HSE and ordinary PBE/Taylor/PT-CN regressions pass.
- [x] Native KY trailing columns, continued records and implicit signed exponents
  accepted; non-KY content coverage and legacy GS physical identity remain open.
- [x] SCF threshold failure warns and preserves usable finite results, with explicit
  converged=False and residual/tolerance/iterations; warm start/TDDFT remain usable.
- [ ] Band-iteration convergence policy remains separate from SCF handling.

Local source suite: 581 tests, one skip. Source-distribution-to-wheel and native
library hash check passed. No GitHub Actions used. No numerical kernel changed.

SCF warning-first validation: actual Si one-iteration SCF, reopen/force access and
2-step TDDFT current completed with finite outputs. No native kernel changes.

## 2026-11-18 — k/PP/GS structural admission

- [x] Generated uniform mesh, native odd Gamma behavior, affine retained stars,
  representative count/weights and bounded native ownership checks.
- [x] PP species references/readability/suffixes and plain scalar KY records.
- [x] Legacy global GS metadata/occupations/byte lengths; optional explicit deep scan.
- [x] Finite occupation range/sum differences remain warning-only and unmodified.
- [ ] Physical GS identity requires a richer archive; legacy records lack it.
- [ ] Non-KY content admission and general native initialization recovery remain open.

General cells and learned/PT-CN symmetry remain subsequent native work, unchanged
by these structural checks. Resource budget advice does not enforce a rank/memory cap.

## 2026-11-18 — General-cell and actual-field fast-PR initialization

- [x] Root cause: scalar orbit gauges/k maps/wrap phases assumed aligned axes;
  rotated reduced-k ACE failed Hermiticity while full-k/reference action worked.
- [x] Complete replacement uses full cell/Cartesian maps only for non-diagonal
  cells; positive diagonal arithmetic stays unchanged. Original sources retained.
- [x] General-cell geometry uses native column norms and actual cell metric.
- [x] Per-owner A/B symmetry comparison no longer falsely rejects2.9e-18 roundoff.
- [x] Host-attached zero-impulse initialization has no fixed-z assumption; actual
  A[3] checked per-step, with capabilities identifying corrected native cores.
- [x] Rotated,left-handed,xy-sheared,z-permuted full/reduced MPI1/2 comparisons;
  independent material groups and direct7→32 resume. x/y/circular-xy also covered.
- [x] Original diagonal trajectory bytes retained. Nonbinary diagonal cells
  match old optimized kernels with the admission-only fix.
- [ ] Symmetry learned/PT-CN, arbitrary space groups, broader process ownership,
  larger partitions/environments and portable archives remain next tasks.

Historical axis-aligned-only and fixed-z restrictions above describe older
binaries; use the current guide and corrected explicitly built runtime.

## 2026-11-18 — Symmetry histories and PT-CN

- [x] Admit existing Taylor strict/learned and PT-CN off/learned with explicit
  native history capability; old factor-off Taylor requirements unchanged.
- [x] Rank-local representative history payloads retained; no new restart format.
- [x] PT-CN electron target slices active weights to match reduced rocc shape.
- [x] Read-only local diagnostics prove actual history readiness/prediction use.
- [x] Full8/reduced4-k, MPI1/2, independent material groups, splits3/7/13 direct
  resume and byte-identical future current/energy/history diagnostics.
- [x] Learned Taylor/PT-CN deep numeric overflow after training preserves state
  bytes and output files; retry matches an independent continuous baseline.
- [x] FDTD/LD/SBE/learned-HSE mixed examples match full-k fast-PR reference.
- [x] Nonuniform27→15-k learned PT-CN with uneven MPI ownership, own resume;
  circular-xy learned PT-CN full/reduced comparison.
- [x] Native resource/checkpoint probes: MPI2 processes with independent cores,
  OMP1/2, full/ACE, stars/orbits, optional cache presence across fresh owner.

Optional phase-cache allocation is not MPI identity; its trusted phase span still
must match before allocation. Learned CN interval1 is an existing method contract,
not an extra accuracy restriction. Historical learned/PT-CN-unavailable notes
above are superseded. Remaining next work: broader state/input/I/O/MPI/timer
ownership and stable native handles, larger environments/scaling and archives.

## 2026-11-18 — Context-owned RT transactions and restore candidates

- [x] Move protected numerical/XC/HSE buffers, output metadata and pending direct
  restore candidates into each fast-PR `tddft_context`; no whole-context copies.
- [x] Internal same-model non-symmetric A/B contexts: distinct dt/capacity and
  z/circular-xy waveforms; isolated vs interleaved numerical probes and traces.
- [x] Simultaneous transactions, independent rollback, simultaneous staged
  candidates, commit/cancel, close A with B's candidate/step still live.
- [x] Reopen acquires available output units and updates its owning descriptor;
  rolled-back close skips units already reused by another owner.
- [x] Reconnect/flush failure retains ownership metadata collectively until
  success; missing-file close/retry leaves B's outputs intact.
- [x] Taylor, learned Taylor/PT-CN; internal MPI1/2, two material groups,
  OpenMP1/2. Compare current/energy plus final orbital/density/potential/exchange
  probes; compare ordinary/late-open RT logs as well.
- [ ] Continue remaining initialization/symmetry/input/I/O/timer/MPI ownership
  and generation-checked C handles. Public one-native-lifetime stays in force;
  this does not claim heterogeneous or threaded same-process TDDFT support.

The serializer allocation limit is still synchronous process-local scratch.
Numeric restart format and time-update formulas/order are unchanged. Host output
journaling remains responsible for file-byte rollback. See
`docs/rt-context-ownership-validation-2026-11-18.md` for scope and verification.
- [x] Standard embedded close releases its actual energy output unit; MPI2
  normal/protected-close probes preserve host units. Source/installed590-test
  suites (one skip each), sdist→wheel/hash and mixed/deep regressions pass.

## 2026-11-18 — Generation-checked same-model C contexts

- [x] Version1 host-MPI fast-PR C header and opaque int64 slot/generation tokens;
  independent original-GS states with per-state dt/capacity/A/failure status.
- [x] Dynamic registry growth without whole-context copies; closed tokens remain
  invalid after slot reuse; root close leaves child numerical/restore state live.
- [x] Collective invalid/mixed token rejection and mismatched initial dt/capacity
  rejection before model lifetime consumption; corrected initialization retries.
- [x] Independent simultaneous transactions/candidates, recoverable numerical
  failure, staged restore, fresh direct9→32/33 continuation without replay.
- [x] Learned Taylor/PT-CN: two material groups, internal MPI1/2, OpenMP1/2;
  independent current/energy/history and complete numeric HSE payload equal to
  preceding singleton controls. External C sample compiles and runs on MPI2.
- [ ] Python shared-model leases/direct use of the new handle ABI, heterogeneous
  same-process inputs/symmetry, common initialization/I/O/timers and per-context
  MPI subgroup retirement. Existing Python one-lifetime behavior remains.
- [ ] Canonicalize inactive legacy checkpoint workspace before portable archives.
  Examined periodic fixed-ion files contain unused tpsi overlap flag, Qp_e tail,
  projection E_old, isolated Poisson counts and +U energy. Do not mistake whole
  archive byte differences for propagated-state differences. No format change
  in this slice; raw diagnostic archives remain private local evidence.

See `docs/fast-hse-native-handles-validation-2026-11-18.md` and the JP/EN
`NATIVE_HANDLES` guide. This is bounded handle identity/dispatch, not completed
arbitrary native reentrancy. Child MPI communicators currently retire at detach.

## 2026-11-18 — Python shared full-k HSE model

- [x] Push preceding ownership/C-handle commits to independent SALMONpy branch.
- [x] backend.open_model + NativeFastHSEModel.new_state: one verified GS work
  copy/attachment, independent state dt/capacity/A/current/clock/output labels.
- [x] State-bound ctypes shims preserve underlying native signatures; leases
  retain the real parent context until closed model + last closed state.
- [x] Native model protocol1 allocates after zero live states; no hidden numeric
  template. Unattempted empty model releases its unused Python reservation.
- [x] Owner-scoped output journals and public native_step_targets preserve
  sibling outputs during rollback; same-owner transaction exclusion retained.
- [x] Learned Taylor/PT-CN, MPI1/2 and OMP1/2: independent/mixed/preceding
  singleton trajectories identical. Save8→A32/B34 resumes without replay;
  corrupted fresh restore can be retired without changing original states.
- [x] Shared model close stops creation while a survivor continues; last close
  detaches. Idle-model fresh recreation is tested twice in isolated controls.
- [ ] Automatic multi-point ResponseSession/rank-plan integration for this
  explicit factory, heterogeneous/symmetry shared inputs and per-context MPI
  communicator retirement remain separate; ordinary default MPI plans unchanged.

See SHARED_HSE JP/EN and examples/shared_fast_hse.py. Numerical source kernels
and checkpoint format unchanged; model allocator reuses the existing constructor.

## 2026-11-18 — retained GS to TDDFT bridge

- [x] Owned orbital/occupation/input seed at the existing RT restart-read boundary.
- [x] Python `NativeGSState.tddft_model` and existing-GS `session.load` entry.
- [x] Existing-file vs memory current parity; sibling/source-close independence.
- [x] Direct checkpoint continuation without replay; MPI orbital/spatial checks.
- [x] Synthetic real/complex empty-band redistribution; source regression (637 tests).
- [x] Final installed-wheel637 tests (one skip), installed MPI2 public sample and
  canonical hash/diagnostic evidence: calc2026/salmonpy-memory-gs-tddft-20261118.
  Local OneDrive write verified; cloud synchronization unverified; no Actions dispatched.
- [ ] General subgroup repartitioning, HSE and concurrent/threaded GS/RT ownership.

No new SCF/reference GS generation is authorized in this validation slice.

## 2026-11-18 — component regression consolidation

- [x] Shared LD/SBE point-MPI test body, preserving model-specific checks; redundant
  individual entry scripts and their model-selection compatibility branch removed.
- [x] Explicit handles/memory-GS matrix runner; independent MPI lifetimes, per-case logs,
  timing/status index, dry-run and failure reporting; no automatic builds/SCF/Actions.
- [x] Four LD/SBE pre/post JSON records exactly equal, including currents/cleanup counts.
- [x] Existing-GS MPI1/orbital2/spatial2 checks pass through the common runner.
- [ ] Other native HSE/TDDFT scenarios remain separate where reference data or lifecycle differs.

See docs/component-tests.md. Unit/error-case coverage remains unchanged.

## 2026-11-18 — unit-test pruning after consolidation

- [x] Repeated module-import contracts folded into public API import checks,
  preserving explicit optional-MPI-unavailable and ASE/f90nml lazy-import checks.
- [x] Five existence-only smoke tests removed; inventory or actual operation tests remain.
- [x] GS link tests moved into native-GS tests; standard/HSE shared-model classes
  grouped without removing numerical/ownership/failure assertions.
- [x] Final626 tests pass (one skip), down from637; files87, down from89.
  Evidence: calc2026/salmonpy-test-consolidation-20261118/unit-pruning.
  No new native builds, GS runs or GitHub Actions.

## 2026-11-18 — functional test consolidation

- [x] LD numerical/ownership/input tests grouped in test_ld.py; SBE numerical/import/
  ownership checks in test_sbe.py, with independent failure diagnostics retained.
- [x] Six superseded test methods removed; unique immutable-settings/parser-type/
  accepted-state checks transferred, repeated worker-start injection removed.
- [x] Selected49→43 tests execute the same581 production lines (not branch coverage).
- [x] Full620 tests pass, one skip; unit files83. MPI/physics/other-backend tests retained.
- [x] Removal mapping and selected logs/JSON recorded in canonical
  calc2026/salmonpy-test-consolidation-20261118/functional-cleanup; no GS/CI/build.

See docs/functional-test-cleanup-2026-11-18.md.

## 2026-11-18 — different-input fast-PR HSE registry

- [x] Owned input/unit values, complete symmetry and density-map workspace;
  no copying of numerical/XC/FFT/MPI owners; fresh coordinate parse allocations.
- [x] Negative generation model tokens and positive state handles; register/spawn/
  close and per-state model activation; old single/same-model behavior retained.
- [x] Explicit Python open_registry/open_model, collective selection/cwd admission,
  deterministic close and retained live children; idle child preserves reservation.
- [x] Six MPI/OMP/general-cell/density-map cases: old-independent current/energy/
  history, fresh state and direct7→24 resume exactly equal, no time replay.
- [x] Density-map sharing bug reproduced before fix and resolved on original input.
  Scalar IEEE NaN/+Inf/-Inf checks cover all6 old/new A entries, one-rank fault.
- [x] Old same-model protected mixed probe matches8 stored traces; source620 and
  installed620 tests pass (one skip each). Installed MPI2 sample and strict old-
  independent comparison pass; JP/EN/AI/C header/example and wheel synchronized.
- [x] Source/all-object prohibition audits; canonical evidence (no GS collection)
  calc2026/salmonpy-hse-model-registry-20261118. Local cloud sync unverified; no CI.
- [x] Ordinary heterogeneous HSE material-ID/rank-plan assembly completed below;
  changed topology, general failures, platforms/physical validation stay separate.

See docs/hse-model-registry-validation-2026-11-18.md.

## 2026-11-18 — ordinary heterogeneous HSE material groups

- [x] Commit preceding verified work as cd0841fe [skip ci]; no push/Actions.
- [x] registry_model=True/open_group session resources and model_groups layer
  assembly; matching library/root/symmetry and rank checks, existing defaults kept.
- [x] Offline material-ID catalog/layout templates and mixed_hse_models example:
  native FDTD/LD + SBE + two different HSE models, independent per-point states.
- [x] Fixed symmetry-file unit1001 collision and unique registered initial labels;
  protected output capture retained. Current/energy/history still equal old control.
- [x] Grouped4 vs separate6,32-step A/E/j and material energy columns exact;
  aggregate protected save8 +24 resumed steps exact, zero replay.
- [x] Source620/installed620 tests (one skip each), installed example exact;
  guides/AI/templates/source inventories and canonical selected evidence updated.
- [x] Source/application-object audits; no GS collected/new SCF, no Actions.

Evidence: docs/hse-layer-assembly-validation-2026-11-18.md and canonical
calc2026/salmonpy-hse-layer-assembly-20261118. Local cloud sync unverified.


## 2026-11-19 — retained GS portable array export

- [x] Commit mixed HSE assembly as2340787d [skip ci]; no push/Actions.
- [x] NativeGSState.export_archive connects serial Gamma DFT-owned arrays to GS/2;
  physical density bounds remove asymmetric orbital halos; density reconstructed.
- [x] Material ID/runtime/GS identity/convergence/local-channel provenance validated;
  warning-only unconverged exports, new destinations only, source state unchanged.
- [x] Synthetic noncubic complex arrays roundtrip exactly; reopened arrays independent
  after source changes/close, malformed layout rejected. Source622 tests, one skip.
- [x] JP/EN/AI guides and explicit opt-in serial sample updated; sample SCF not run.
- [ ] Real retained-GS export/TDDFT parity, retained multi-k/MPI export, imported-owner
  physical metadata and native handle restore remain unverified/not implemented.
- [ ] Complete SCF checkpoint (including mixing) and changed-decomposition restore.

No new SCF, no real GS data collected, no native rebuild, no CI. Portable array
readback is not native GS restoration. See docs/retained-gs-export-validation-2026-11-19.md.


## 2026-11-19 — retained serial regular-k export

- [x] Push preceding Gamma exportb66b97fd and mixed HSE2340787d to salmonpy remote;
  remote branch SHA verified. [skip ci], no manual Actions dispatch.
- [x] Immutable DFT regular full mesh export supports2/3/6/8 k points: exact
  wavefunction/occupation/eigenvalue payloads, native x-fastest ordering, uniform
  weights and weighted-electron density; independent readback after source close.
- [x] Compare both synthetic real/complex Gamma archives against b66b97fd:
  all saved-file SHA256 values identical, not only decoded numerical values.
- [x] Source623 tests, one skip; guides/AI metadata and opt-in sample synchronized.
- [ ] Real native export/TDDFT parity, imported-owner physical identity, retained MPI
  export/native restoration/full SCF mixing checkpoints remain future work.

No new SCF, real GS collection, native rebuild or wheel rebuild. Evidence:
docs/retained-gs-multik-validation-2026-11-19.md. Local storage only; sync unverified.


## 2026-11-19 — retained GS collective array export

- [x] Native bounds drive k/orbital/spatial tile extraction with halo removal;
  exact nonoverlapping global coverage and replicated metadata checks before transfer.
- [x] Root-only publication/return, one rank tile per collective, collective allocation/
  write error propagation; preserved GS arrays and existing outputs, retry succeeds.
- [x] Actual MPI4 synthetic runner: five partition cases and collective bounds/
  replicas/overlap/gap/metadata/options/allocation/write/existing-output checks.
- [x] Source623 tests (one skip), Gamma archive file parity withb66b97fd; JP/EN/AI
  metadata/examples and packaging updated. No new SCF/real GS collection/electronic build.
- [x] Final sdist/wheel source+object audits/hash checks; installed623 tests (one
  skip) and installed actual-MPI4 all14 cases pass. Audited wheel LD/SBE rebuilt.
- [ ] Actual native GS export/TDDFT parity, imported-owner physical metadata,
  native GS handle restore/full SCF mixing checkpoint/changed-topology restoration.

See docs/retained-gs-mpi-export-validation-2026-11-19.md for final wheel/results.


## 2026-11-19 — existing-GS referenced k-table import

- [x] Commit/push regular-k/MPI exportd67ce847 to salmonpy remote [skip ci].
- [x] Copy/hash/retain native file_kw assets with relative directories; preserve
  quoted ! in filenames. Missing required files fail before native initialization.
- [x] Shared existing-GS runner --external-k: PZ2-k MPI orbital2/spatial2,
  Libxc-PBE8-k orbital2, TB-mBJ8-k spatial2 and unchanged regular control.
  Currents, independent states and direct no-replay restart are bit-exact.
- [x] Missing k-table while RT is active leaves GS resources/current unchanged;
  subsequent valid loading and RT continue. No new SCF/reference GS generated.
- [ ] Imported canonical physical model/portable export, full native GS/SCF mixing
  restoration and changed-topology restore remain; no claim of physical convergence.

See docs/gs-external-k-import-validation-2026-11-19.md for final records.

- Final import slice: source624/installed624 tests (one skip each), installed native
  external-k orbital2/file/current/resume/missing-asset comparison pass. Offline
  sdist/wheel audits/hashes pass; original electronic library reused by recorded SHA.
  Evidence: canonical calc2026/salmonpy-gs-k-table-import-20261119 (no GS arrays).


## 2026-11-19 — local declared GS input descriptions

- [x] GSInputMetadata/describe_gs_input and retained input_metadata snapshots;
  optional f90nml metadata extra, no dependency/native/MPI work on ordinary load.
- [x] Explicit cell/coordinate units/species/grid/k/XC/PP descriptions, unknown
  defaults/unsupported fields remain unknown. Source ID is not physical identity.
- [x] Captured input/assets and convergence remain available after GS close,
  detached immutable/JSON results; imported provenance/convergence remain unknown.
- [x] Synthetic units/indexed species/nonorthogonal/external coordinate/unknown/
  dependency/mutation tests; existing PZ/PBE/TB-mBJ MPI current/resume controls pass.
- [ ] Full physical GS identity/native effective model binding/imported portable
  export, SCF-mixing checkpoint and changed-topology restoration remain separate.

See docs/gs-input-metadata-validation-2026-11-19.md for final wheel/results.

- Final metadata slice: source626/installed626 tests (one skip each), installed
  native MPI4 and public metadata example MPI2 pass. Import loads no parser/MPI;
  offline wheel source/object/hash audits pass. No new SCF/electronic rebuild/CI.
  Evidence: canonical calc2026/salmonpy-gs-input-metadata-20261119 (no GS arrays).


## 2026-11-19 — portable arrays to retained native GS

- [x] load_archive: matching collective selections/environment, root-only explicit
  conversion into new caller-owned bundle, receipt/hash verification and GS load.
- [x] Immutable archive_origin retains input archive ID/convergence claim separately
  from current imported convergence(None); default ID inherited, explicit ID honored.
- [x] Known marker/integer-record mismatch rejects; legacy unknown compiler/source
  evidence warns. No new restrictions on ordinary GS run/load or implicit compiler.
- [x] Audited/probed converter uses unchanged extracted I/O from existing validated
  serial runtime. No solver rebuild/new SCF/real GS collection.
- [x] Synthetic native orbitals/occupations exact; MPI2 orbital/spatial/k contexts
  and PBE orbital compare currents to independent file route; direct8+8 resume,
  independent states, conversion failure ownership and retry exact.
- [ ] Complete native SCF mixing/cache/iteration restore and arbitrary RT checkpoint
  repartitioning remain separate. Archived energy/eigenvalues are not native records.

See docs/portable-retained-gs-validation-2026-11-19.md for final wheel/evidence.

- Gamma complex-array import explicitly stages an equivalent k table before GS
  initialization, avoiding native DFT's real allocation. Actual Gamma orbital/
  spatial readback/current/independent/direct-resume controls pass. Full-SCF
  restoration remains separate; no imaginary component is discarded.

- Final portable-load slice: source627/installed627 tests (one skip each), six
  actual native synthetic cases (PZ/PBE, Gamma and two-k) plus installed PBE/Gamma
  and repeated shared-session public sample pass. Converter/wheel audits/hashes
  pass. No new electronic build/SCF/real GS collection/CI. Evidence:
  calc2026/salmonpy-portable-retained-gs-20261119. Complete SCF caches remain open.


## 2026-11-19 — GS/2-imported retained re-export

- [x] Capture small source metadata/eigen/geometry/k arrays, no wavefunction cache;
  GS/2-imported live owners re-export after source removal. Legacy models not guessed.
- [x] Native orbitals/occupations and derived density separated from carried source
  model/sampling/eigenvalue/energy; new validated reexport provenance, con=None.
- [x] Complex reads preserve signed-zero bits without arithmetic; NaN/Inf rejection
  remains collective. Computed-DFT export arithmetic/old provenance preserved.
- [x] Source synthetic MPI Gamma/PZ/PBE/k/orbital/spatial/nonuniform save/load/current
  comparisons and active-RT partial-write failure/source-owner preservation/retry pass.
- [ ] Full native SCF mixing/cache/iteration snapshot and arbitrary RT checkpoint
  migration remain separate; reimport creates a new model with its own identity.

See docs/gs-reexport-validation-2026-11-19.md for final wheel/evidence.

- Final re-export slice: source628/installed628 tests (one skip each), seven source
  actual-MPI synthetic cases and installed PBE/Gamma/nonuniform cases pass. Public
  sample repeats with reexport_to and legacy current hash unchanged. Computed
  archive b66b97fd file parity, source removal, signed-zero and nonfinite checks pass.
  Evidence: calc2026/salmonpy-gs-reexport-20261119; no real GS arrays/SCF/CI.

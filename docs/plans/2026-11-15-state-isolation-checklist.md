# Active state-isolation checklist

Continuation of handoff §7.3–4 and the approved GS/TDDFT ownership work. This is not a replacement for the wider project backlog. Historical plans record their completion-time scope; later evidence below extends it. Direct mixed-session persistence, non-PZ inner-MPI validation and protected ordinary mixed restoration are completed within the recorded scopes. Wider isolation, resources and distribution continue as separate objectives.

## Current audit — 2026-11-17 (takes precedence over historical status notes)

The dated sections below retain the scope at their completion time. Older statements such as "occupation remains shared" or "multiple field ranks remain open" are superseded by the later evidence, not separate unfinished tasks.

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
5. [ ] General native failure contract: standard semilocal TDDFT deep numerical recovery is now implemented below; initialization/unsupported native errors and process-terminating paths remain open. Handled mixed restore staging failures already preserve originals and permit continuation/retry; do not list that completed subset as missing. The CLI still exits on an unhandled error.
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

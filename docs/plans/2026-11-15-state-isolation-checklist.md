# Active state-isolation checklist

Continuation of handoff §7.3–4 and the approved GS/TDDFT ownership work. This is not a replacement for the wider project backlog. Historical plans record their completion-time scope; later evidence below extends it. New XC coverage is deferred until the existing isolation work is addressed.

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
- [ ] F4b: explicit field-position/time interpolation adapters and native coupling
  driver integration; do not infer molecular dipole-to-current conversion.

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
| M6 | Serial/MPI postprocessing identities and managed TDDFT per-point disk resume | MPI whole-session restart, other backend disk formats, branch and in-memory point identity semantics |

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

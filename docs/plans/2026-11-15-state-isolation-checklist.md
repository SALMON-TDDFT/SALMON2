# Active state-isolation checklist

Continuation of handoff §7.3–4 and the approved GS/TDDFT ownership work. This is not a replacement for the wider project backlog. Historical plans record their completion-time scope; later evidence below extends it. Direct mixed-session persistence, non-PZ inner-MPI validation and protected ordinary mixed restoration are completed within the recorded scopes. Wider isolation, resources and distribution continue as separate objectives.

## 現在の状態（2026-11-23 作業・検証照合）

この節と末尾の最新検証記録を現在の状態として扱います。過去の日付の
未完了チェックや「未コミット」は、その時点の記録であり、独立した残件として
数えません。未完了の物理条件を、今回の局所的な検証だけで完了扱いにしません。

| 状態 | 項目 | 確認した範囲 |
|---|---|---|
| 完了 | 実用上選定した異常復帰・操作契約 | 主なLD作業領域、入出力、電子状態・FDTD・GSの呼び出し契約。任意のOS/MPI障害や全ての小さな確保までの保証は対象外 |
| 完了 | 現行CPUのNLCC認証・単位照合 | PZ/Libxc-PZ/PBE局所試験、MPI2/4、固定Si/PBE比較。今回のソース・バイナリ向け記録も更新。換算定数は変更なし |
| 完了 | 初期推定の入力準備・保持GSからの受け渡し | 準備と実行を分離し、保持GSの波動関数・占有数を独立コピー。前回実装は8f3a7663までpush済み |
| 完了 | アーカイブ・自作readerからの初期推定 | 共通配列をFortranへ直接渡す。PZ/高速HSEのMPI、独立性、拒否後の再試行、既存経路との一致を確認。本項の実装・文書を今回のcommitへ収録 |
| 完了 | 最新パッケージ・分子接続検証 | 二つのFDTD＋分子の直接再開、非対角背景・テンソルLD・界面平滑化・Ohmic散逸の収支観測、各部品の性能モードを7d27edbaまでcommit/push済み。分子E駆動も通常／性能・直接復元で一致。今回の外側MPI対応はwheel707件（706成功・任意スキップ1）。CI結果は未確認 |

### 直近の実装残件（2026-11-23 外側MPI検証後）

- [x] LD／SBE、直列の場、通常TDDFT、高速PR版HSEの単一／共有／登録状態、周期系・分子の別プロセス型の性能モードを、各検証記録の範囲で実装。
- [x] 外側A→J MPI結合の性能モード。自作物質モデルで3/5ランク、実E/H＋LD＋2点TDDFTで3/6ランクを確認。後者は96ステップ、通常の連続実行と明示的保護復元・再試行が一致。wheel707件中706成功・1スキップ。
- [x] **高速PR版HSE＋実E/Hの外側性能モード検証**。Taylor/PT-CN×全/縮約k点×3/6ランク、ソース24・wheel24条件で完了。末尾の検証範囲を参照。
- [ ] **A/E/Hの場の受け渡し・分子の外側性能スケジュール**。分子ワーカー単体の対応と区別する。
- [ ] **別プロセス型高速HSEの経路・必要性の確認**。共有／登録型を根拠に対応済みとせず、不要な経路は増やさない。
- [ ] **保護付き混在計算の規模別計測**。多数点・長時間のメモリ・通信・観測I/O、実際の性能差を測る。

性能モードは明示的に選択し、既定は検査付き。明示的な保護更新は退避・復元を維持する。性能モードの通常更新に同じ異常復帰を約束しない。過去の「次：分子」等は後続の完了記録で更新されている。

### 残る作業と対象拡大

**分子の双方向FDTD結合は、C2H2・PZの開発試験まで接続済みです。** 真空・LD環境で独立E応答との一致と取消・再試行を確認しました。電子状態の取消はワーカー内保存点から直接復元します。短時間の電子数保存は確認済みです。量子系の仕事量は観測・比較と選択電場の中点時間配置まで確認済みです。光源・吸収層を除いた内側の箱について、六面流束を含む場・等方LD・分子の離散エネルギー収支と実MPI2分割を確認済みです。単一分子用の一般利用向け入口はMolecularFeedbackとして実装済みです。正の背景誘電率・透磁率と対角LDの収支観測まで完了しました。背景誘電率が対角の場合は非対角LD極の収支観測も実装・ソース検証済みです。非対角の背景誘電率の収支観測も接続・検証済みです。LD界面平滑化の収支観測も実装・検証済みです。Ohmic伝導の散逸観測も実装・検証済みです。広い物理条件での収束確認は残件です。以下はそれと区別する運用・対象拡大の残件です。

1. **施設MPI・多ノードでの実証**：サイトのランチャー、実際のジョブ単位検証記録、ランク配置・多ノードの検証。Mac上の少数ランク試験を代用しない。
2. **大規模化の計測**：マクロ点・ランク・観測間隔に対するメモリ、通信、ファイル数・I/O負荷の測定。集約・バッファ化は測定上の必要に応じて実装する。
3. **他環境・配布の検証拡大**：重い計算ライブラリと初期推定経路について、Linux・Intel Mac等および施設環境での構築・配布・実行証拠を追加する。既存LD/SBE wheelの検証とは区別する。
4. **Agのフィット改善**：局所的な形状誤差と最適化終了を改善し、現状のneeds_refinementを根拠なく解除しない。
5. **追加対象の設計・検証**：追加保存形式、GPU、現在の独立状態経路で未検証の物理条件。必要性と対応範囲を決めて進め、全機能への包括的対応を主張しない。

高速HSEの電子数診断と実材料の格子収束調査は利用者指定で後回し、
分子は2026-11-22の利用者指示で再開（双方向FDTD結合を先行）。TTMは削除済みです。格子/k点の暗黙の補間は追加しません。
メモリ予算は警告とし、SCF未収束時は結果・残差を保持して確認を促す方針を維持します。
重いCIを追加で手動起動せず、ソース・目録・配布物の必須整合性検査は省略しません。

## Current actionable queue — 2026-11-22 (resource inventory and operational follow-ups)

This section takes precedence over historical dated checkboxes below. Historical
unchecked scope statements are not independent tasks; retain them as evidence.
Molecular work resumed2026-11-22: actual C2H2 PZ E-worker/FDTD wiring and retry are tested in vacuum/LD. Short-run measured electron-number conservation is verified. Quantum-work observation/comparison and selected-field midpoint timing are verified. Six-face inner-box energy balance is verified for scalar LD, including real two-partition FDTD tests along x/y/z. Ohmic-conductivity energy observation is now verified below. Physical convergence remains open. LD subpixel energy observation is now verified below. Non-diagonal background observation is now integrated and verified below. Off-diagonal LD pole strengths with diagonal background are source-verified; wheel validation of that increment is also verified below. Positive diagonal electric/magnetic backgrounds and diagonal LD energy observations are now verified below. MolecularFeedback, live direct cancellation and process-exit molecular direct disk restart are implemented and verified in the latest entries below. Coordinated direct disk persistence of both fields, the molecule and explicit accepted step indices is now verified for MPI6/same-layout C2H2 PZ; see the latest entry. Memory budgets warn; SCF nonconvergence retains
results/residual diagnostics and pauses the downstream workflow for user review
(the explicit continuation policy is documented below). GS saving/reuse is the normal workflow. Collected GS validation datasets require the applicable campaign-specific user authorization; existing grants do not authorize arbitrary new campaigns. Store authorized results in calc2026. Do not manually dispatch optional heavy CI. Critical source/inventory/packaging commits must not skip required integrity CI.

### Completed and verified

- [x] Practical tensor-update workspaces: call-owned drive/solver/off-diagonal
  buffers allocate before physical mutation and reuse the existing input gate.
  Known24 returns RuntimeFailure with the live field still usable; no E/H,LD,
  clock,output or checkpoint mutation/replay. No extra collective or physical
  checkpoint payload. GNU MPI1/2 source failure/retry controls pass at early,
  middle and later steps. Production/testing results and direct midpoint data
  match the previous library. Standalone32-step tensor/subpixel MPI2 results
  match the previous executable. Source/wheel679-test suites pass(one optional
  skip each); two independent2-rank groups also pass workspace failure/retry.
  No blanket backup/OS/MPI failure guarantee.
  Record: calc2026/salmonpy-tensor-step-workspace-recovery-20261122.

- [x] Off-diagonal LD grid coefficient storage:5D weights,4D diagonal values and
  integer material indices allocate before the existing material-vector gate;
  values initialize at their original position. Partial STAT failure returns1120
  after field-group agreement and owned cleanup, allowing constructor retry.
  No startup/step collective added, arithmetic and disk format unchanged.
  Source GNU MPI1/2 production/testing controls and repeated partial failures
  pass; scalar/diagonal/tensor fields and direct midpoint data match5f255105-era
  controls byte-for-byte. Source/wheel679-test suites pass(one optional skip
  each). Subpixel normal/direct-restore source/wheel controls also match the
  previous library on GNU MPI1/2; its allocation-failure coverage is unchanged.
  Subpixel allocation and tensor propagation workspaces remain separate scopes.
  Record: calc2026/salmonpy-ld-tensor-grid-allocation-recovery-20261122.

- [x] LD coefficient preparation: material list, six common1D/2D and three
  directional2D/3D arrays now return1120 on partial allocation failure. The
  existing directional status broadcast becomes a max reduction; scalar LD
  carries failure to the existing material-vector gate. No new startup/step
  collective. GNU MPI1/2 source controls and repeated failure/retry cases pass;
  production/test fields, histories and direct midpoint checkpoints are bitwise
  identical to the previous e3666160-era controls. No additional input restriction
  or old-format fallback. Memory budgets remain warning-only.
  Record: calc2026/salmonpy-ld-coefficient-allocation-recovery-20261122.

- [x] Removed old build/preparation positional arguments, Python chdir runtime
  rescue, old mixed-checkpoint schemas and identity rescue, inferred GS restart
  flags, SCF rebuild/1 and missing codec/revision rescue, implicit HSE settings,
  old EXX input aliases and missing RT reference-density rescue. Current custom
  backends, numerical reference kernels and explicit execution modes remain.
  Evidence: calc2026/salmonpy-remove-backward-compatibility-20261122.

- [x] Directional LD input-array allocation failures: diagonal/full tensor
  partial allocation returns1120 before array broadcast; group cleanup and
  constructor retry are verified with GNU MPI1/2. Production/testing scalar,
  diagonal and rotated-tensor histories and direct midpoint checkpoints match
  bitwise; scalar/tensor prior controls also match. Existing startup broadcast
  is replaced, adding no collective. Source/wheel evidence is recorded separately
  in calc2026/salmonpy-ld-material-allocation-recovery-20261122.
  No backward-compatibility path is added. Broader workspace failures remain open.

- [x] Covered ordinary-XC/fast-PR HSE SCF memory and same-binary/model/topology disk continuation, residual CSV/plots, protected numerical staging and common disk binding. Ordinary-XC pure k/orbital/spatial and mixed-layout changed-rank restart are measured. Dense fast-PR HSE full/reduced-k × orbital × spatial mixed restart is measured below. Same-source rebuilt binaries and exact declared metadata-only source transitions are verified; arbitrary changed-numerics/cross-platform migration and packed/tiled/sparse formats remain outside the measured scope.

- [x] One maintained fast-PR solver tree and shared default preparation/build path;
  ordinary XC does not require HSE dependencies. Historical source evidence retained.
- [x] Ordinary TDDFT measured electron number: read-only Fortran density integral,
  process protocol and collective native/shared/heterogeneous-handle interfaces.
  State independence, direct/replay restore and prior-current comparisons verified.
- [x] dt campaigns: norm target1e-5; finite-value failures end only the failed trial.
  Recovery-capable MPI returns the pre-existing gross-norm failure instead of STOP.
- [x] Primary grid criterion: current relative RMS error<=5% at norm-qualified dt,
  common physical times and finer-grid reference; confirmed SCF-converged pairs
  stop before0.1A. Scalar energy comparison remains supplementary.
- [x] Mac/internal-MPI2 Si16³/~2fs: SCF converged,dt0.2 fails at step8;dt0.1
  completes with maximum relative electron loss4.22e-10. This is not an optical
  or general-grid convergence certificate. Existing12³→16³ current error18.83%
  does not pass5%;12³ SCF was unconverged, so that comparison is provisional.
- [x] Latest completed communication verification: source/wheel668 tests, one
  optional skip; status-bearing A/J transport, failure-only reporter Bcast,
  MPI4/6 protected retry and direct resume, paired manuals/AI JSON verified.
- [x] Optical duration/window/k diagnostics and reference-frequency LD-only
  subpixel modules completed below; historical pending optical checkboxes are superseded.
- [x] Operational FDTD input-only staging and site launcher argv verified below.
- [x] Fresh DFT and explicit native RT input-only staging verified below.
- [x] Recent FDTD cold-init recovery is verified within recorded scopes:
  early grid/boundary/time1101–1104, explicit MPI/grid/mapping1105/1106/1107, automatic factorization1117,
  material-ID CUBE1115/1116, LD/subpixel file1110/1113 and coefficient/special-mode
  returns1111/1112/1114. Source/installed GNU MPI1/2 controls match prior fields,
  point/material/shape data and direct checkpoints. Automatic-factorization
  failure→explicit retry is additionally measured on local MPI7.
  Tiny reader staging-allocation1116 has not been fault-injected; this is not
  blanket allocation/OOM/process-failure recovery. Bloch/tensor evolution remains
  unsupported; safe rejection/retry does not implement that mode.
- [x] Owned file units for ordinary point/plane/shape, common EM restart/copy,
  LR/Fourier-plane and ASE/ART output. Caller units600/601/602/700/701 are checked
  in the relevant separate fixtures. Six spectral cases × GNU MPI1/2 compare
  occupied/unoccupied/prior-runtime outputs exactly; periodic LR output now closes.
- [x] TTM cancelled/removed by the project lead: implementation, automatic
  activation, build references, examples, tests and utility removed. Obsolete
  input is inert. Current sdist/wheel exclude TTM and historical source archives;
  frozen original records/archives remain local evidence, not supported engines.
- [x] Native core3D allocation STAT failures return1118 after field-group
  phase agreement; partially allocated material/core E/H arrays and MPI caches/
  child communicators are released. Source/installed GNU MPI1/2 tests inject
  partial material and field failures, retry twice each and directly resume.
  Normal production/test/wheel fields, observations, shape and checkpoint bytes
  match the previous control. External-CUBE MPI2 controls also match. Production
  has no injection hook. Two startup field-group agreements, no extra step
  collectives; memory-budget excess still warns. LD-pole/PML/other auxiliary
  allocation, physical OOM/process/MPI failures remain outside this bounded scope.
- [x] Initial LD-state allocation1119:24 owned pole/current/polarization/
  coefficient/history arrays return STAT failure. Source/installed GNU MPI1/2
  scalar/off-diagonal-tensor tests fail after partial4D and3D allocation,
  release resources and retry twice each; combined core/LD fixtures also
  retain1118 behavior.32-step normal/direct-resume results match prior controls
  byte-for-byte. Shared native/solver kernels remain identical; production has
  no injection hook. One LD startup field-group agreement,zero new step
  collectives. Material/coefficient scratch,tensor auxiliary,subpixel-specific
  reallocation,PML/observation allocation and OS/MPI failures remain excluded.
- [x] Material/boundary vector allocation1120:three material and six E/H
  boundary coefficient vectors return STAT failure before coefficient use or
  communication. Partial-state cleanup and protocol-gated Python retry verified
  in source/installed GNU MPI1/2 scalar/off-diagonal-tensor fixtures. The material
  gate now shares the existing LD-state agreement (retaining1119);only one new
  startup field-group agreement with LD,zero new step collectives. Combined
  core1118/LD1119/aux1120 regressions and direct continuation match prior bytes.
  PML bulk3D coefficients were already covered by1118;small automatic PML
  workspaces are not a separate bulk-allocation implementation gap.
- [x] Optional point-output I/O:Fortran header/append open/write/close status
  checks warn and continue physical evolution. Only owned successfully opened
  units are closed;warning is bounded per state/output-owner rank. No new MPI
  communication. Real directory/read-only-file failures and destination repair
  measured on source/installed GNU MPI1/2;later rows resume,missing rows remain
  explicitly incomplete. Fields,shape/info and physical checkpoints match normal
  controls;normal scalar/tensor output bytes remain unchanged. Caller file unit
  lifetime is verified through initialization,steps,direct restore and finalize.
  Other outputs and filesystem-wide failures remain separate scoped work.
- [x] Optional plane-integral output:Fortran header/append open/write/close
  status checks warn and continue without skipping any existing plane/MPI sums.
  Owned warning state,no extra collective,caller file units retained. Source/
  installed GNU MPI1/2 real directory/read-only faults and repairs preserve
  fields/checkpoints;unaffected yz/xz planes and point files remain complete.
  Normal scalar/tensor integral bytes match prior controls. Missing/partial rows
  are reported as incomplete,not replayed. Other output families remain scoped work.
- [x] Optional final observation-info output:obs0_info.data open/write/close
  errors warn and proceed through normal state/MPI cleanup. Completed physics
  is kept separate from missing animation metadata;no data deletion or extra
  collective/recomputation. Source/installed GNU MPI1/2 real directory/read-only
  faults preserve time series,integrals and checkpoints;caller files survive.
  Normal scalar/tensor metadata and other output bytes match prior controls.
  Required restart errors and other output families retain separate policies.
- [x] Optional plane-snapshot I/O:open/write/close failures warn and retain
  every existing plane/component aggregation and workspace release. State-owned
  warning flag,no extra MPI collective. Source/installed GNU MPI1/2 directory/
  read-only faults continue physical evolution;repair plus explicit continuation
  after direct restore writes the affected frame. Normal scalar/tensor output
  bytes and physical checkpoints match prior controls;caller files survive.
  Snapshot staging allocation and frequency-resolved maps remain separate scopes.
- [x] Batch frequency-resolved plane-map I/O:header/row/close failures warn
  and continue remaining windows,planes,energies and cleanup. Fourier sums and
  MPI order unchanged,no new collective. Source/installed GNU MPI1/2 compare
  two energies×three planes×two windows exactly to prior controls. Eight real
  directory/read-only fault cases preserve other11 maps/time series and caller
  units600/601/602;map units are checked closed. Native coupled unsupported
  return15 remains unchanged;native spectral transactions/restart are not claimed.
- [x] Batch LR/ASE/ART optional spectra:Fortran open/write/close errors warn
  and retain Fourier processing,per-window integrals/MPI sums and cleanup.
  Separate window units/status/names;no extra collective or global selector.
  Source/installed GNU MPI1/2 isolated/periodic/non-unit-background LR,ASE,ART:
  20 normal-configuration comparisons match prior bytes;40 directory/read-only
  faults preserve all other outputs and caller units600/601/602. Spectral units
  are checked closed. Native unsupported return15 and required-save policy unchanged.
- [x] Existing native required-checkpoint save contract verified without numerical
  changes: source/installed scalar/tensor GNU MPI1/2, rank-local existing-file,
  directory and missing-parent failures return group errors and preserve live
  field/LD bytes, clocks and observations. Existing target retained; failed partial
  sets quarantined as evidence, never published as archives. Repair/retry and
  isolated direct restore require zero evolution. 72 state comparisons and72
  normal-control file comparisons match exactly. No multi-file atomic-save claim;
  disk-full/close-only faults and legacy batch restart writers remain separate.
  Record: calc2026/salmonpy-fdtd-required-save-verification-20261122/summary.json.
- [x] Legacy batch required EM writers: exclusive sibling staging and checked
  open/write/close/rename; explicit required sticky io_status, no unpublished-API
  compatibility fallback. One owning-group failure agreement per metadata/state
  save phase, none on ordinary steps; batch fails nonzero rather than claiming
  saved data. Original overwrite-before-sidecar-failure reproduced, then fixed.
  Source/installed GNU MPI1/2 helper failures repair/retry;106 normal payload files
  per distribution match old bytes. Caller units600/601/700/701 and remap reads
  remain intact. Real LD+FDTD batch metadata/payload failures terminate without
  hangs,normal outputs match old controls. Four native scalar/tensor source/wheel
  MPI2 regressions preserve prior field/checkpoint bytes and direct restoration.
  Per-file publication only: failed sets can mix generations; preserve a separate
  completed directory. Legacy readers/directory setup, disk-full/close-only faults,
  fsync and atomic multi-file publication are not claimed.
  Record: calc2026/salmonpy-legacy-em-save-recovery-20261122/summary.json.
- [x] Required batch EM restart-directory creation:existing directories now
  return defined status0 (old undefined output reproduced). Root creation uses
  explicit status and replaces the existing barrier with one owning-group
  agreement;all ranks fail nonzero on creation errors. No new collective count
  or step communication. GNU MPI1/2 fresh/existing/directory-symlink LD+FDTD
  controls match prior data and restart bytes;file and file-symlink collisions
  preserve originals and terminate without hangs. Source/installed ownership,
  normal writer-byte controls,679 tests(one skip),and four native MPI2 direct
  save/restore regressions pass. Not a directory-writability preflight;later
  writes remain checked by writers. Committed/pushed as06bf7cea.
  Record: calc2026/salmonpy-restart-directory-recovery-20261122/summary.json.
- [x] Required legacy EM restart readers:mandatory collective io_status with
  checked read-only open/read/close. Shared status precedes payload;local errors
  agree within solver group;remap combines3 metadata broadcasts into1 status
  header and gates workspace/payload transfer. No ordinary-step communication
  changes or compatibility fallback. Source/wheel94 cases each on GNU MPI1/2:
  missing/truncated/invalid-text/directory errors retain caller MPI/unit600 and
  explicit repaired retry reads poisoned buffers. Read-only archives accepted;
  existing1D–5D real/complex ownership/remap tests pass. Actual LD+FDTD64→128
  resume and1↔2 repartition match continuous raw states;normal prior-runtime
  controls are unchanged. Shared input/output repartition failure originally
  overwrote layout metadata;metadata publication now waits for every old payload
  to read successfully and failed input sets retain identical contents. Later
  multi-file saves are not atomic. Copied observation text differs only in old
  trailing blanks;physical values and raw payloads are exact. Both679-test suites
  pass(one skip),plus four native scalar/tensor MPI2 save/restore regressions.
  Failed legacy buffers/partial models are invalid,not atomically restored live
  states. Scope is detected I/O errors,not arbitrary corruption/OS/MPI recovery.
  Record: calc2026/salmonpy-restart-read-recovery-20261122/summary.json.
- [x] Current user/AI documentation readability:104 Markdown documents reviewed,
  explanatory internal labels removed while callable names/paths remain intact.
  Further15-guide pass corrected obsolete LD state-reconstruction and unconnected
  Python-model claims against implementation,translated the Japanese onboarding,
  and added task/state-creation navigation with developer details folded separately.
  Manual11/public-API5 tests,local links and executable-example checks pass.
  Records:calc2026/salmonpy-documentation-readability-20261122/{summary,followup-summary}.json.
- [x] Reader I/O and documentation integration: published as23cb1342.
- [x] Directional LD allocation recovery and backward-compatibility removal
  were published as e3666160.
- [x] LD coefficient scratch recovery was published as5f255105.
- [x] Off-diagonal LD grid-coefficient allocation recovery published as46589f24.
- [x] Tensor-update workspace recovery was integrated in the recorded earlier commits. Retain required lightweight source/packaging checks; no optional heavy CI dispatch.
- [x] Latest source/installed suites:679 executed,one optional skip,no failures.
  Paired current manuals and AI JSON synchronized. These latest local records
  do not supersede older MPI4/6, multi-platform or physical-convergence evidence.
- [x] Integration checkpoint: verified grid/layout/input/output/TTM and required
  save changes consolidated for the requested2026-11-22 commit/push. Required
  source-integrity CI remains enabled; no optional heavy workflow is dispatched.
  Published commit:a7073552 to salmonpy/codex/salmonpy-state-inputs-recovery.
- [x] Initial-guess preparation/live/archive/custom-reader integration is completed below. Selected group5 R3/R4 and CPU group3 scopes are completed. Broader site/layout/multinode evidence remains in the current queue; groups1/2 stay deferred.

### Next work after this reconciliation

1. **Group5 R3 — selected practical scope completed:** supported input/read/save
   returns and principal LD initialization/tensor-update allocation paths are
   verified above. Do not keep extending this task to every small allocation or
   standalone STOP. Remaining subpixel initialization/reconfiguration, observation
   staging, compiler temporaries and state-backup allocation are outside the
   selected guarantee; revisit only for an identified practical need. Arbitrary
   OS OOM, trapped FP, broken MPI and process termination are not recoverable
   contracts. No restriction on valid/custom inputs or budget-based hard stop.
2. **Group5 R4 — selected current-operation scope completed:** common host,
   electronic-state, FDTD and GS declarations/bindings are organized and locally
   verified. GS source/wheel SCF continuation/protected restore is completed.
   Optional diagnostics and wider profile/platform/site evidence remain separate;
   these tests are not a blanket stable ABI or arbitrary-failure guarantee.
3. **Group3 — selected CPU scope completed:** maintained-runtime
   NLCC certification covers PZ/Libxc-PZ/PBE local variation/core-force tests,
   MPI2/4 ×OMP1/2, and the fixed Si/PBE application reference. Unit audit preserves
   every conversion constant and numerical result; unit_contract exposes existing
   conventions. Source/wheel679 tests(one skip each) and evidence-reader7 pass.
   Other materials/GPU physical accuracy is not certified by these checks.

NLCC/unit and retained-GS work are committed through8f3a7663; archive/custom
initial guesses are published in7009381c. Those initial-guess implementations
are complete within their recorded scopes. Do not reopen their historical
unchecked rows as separate tasks. New GS campaigns require applicable user
permission; prior grants remain campaign-specific.

The active molecular stage resumed2026-11-22. MolecularFeedback, the PZ isolated
worker, live direct rollback and molecular-only direct disk restart are now verified
within the latest recorded scope. Coordinated direct restart of both fields, the molecule and accepted step indices
is now verified for the recorded MPI6 scope. Tensor/subpixel/background energy-observation coverage is completed below.
Next: outer-MPI performance integration and physical convergence. Other
operational work remains site-MPI/multinode/receipt filesystem evidence, scale
and observation I/O, Ag fit refinement, optional storage/GPU and wider physical
profiles. Preparation-only DFT/TDDFT/FDTD and customizable launcher construction
already exist. Groups1/2 remain deferred; molecular work is active. TTM is deleted.
Historical unchecked entries below are scope-at-the-time notes, not additional
independent tasks or a mandate to generalize every possible configuration.

### Current-format policy and historical evidence

Older completed source-pair migrations and old-schema records below remain
historical evidence. They do not authorize restoring removed compatibility paths.
Current SCF rebuild matching requires schema2, codec identity and state revision;
old missing-revision transition declarations are removed. Current mixed direct
and replay writers use direct/4 and replay/2. Distinct current numerical methods,
custom backend contracts and single/shared/registered execution remain supported.
The remaining ABI task is to define and verify the current caller/runtime contract,
not to maintain older releases. No automatic old-format rescue is planned.

### User-directed execution order (2026-11-19)

Active order:4→5→3. Groups1and2 are deferred at the user's request; do not
prioritize fast-HSE norm wiring or a current-grid campaign now. Molecular work
remains paused. Group6 follows separately with the existing CI budget constraint.
Group4 SCF history, numerical checkpoint/live continuation and ordinary mixed
repartition are complete within the measured scopes. Fast-PR HSE dense Gamma orbital/spatial mixed ownership and restart are also
measured below, together with full/reduced-k orbital/spatial mixed restart and limited rebuilt/source-pair compatibility. Packed/tiled formats and optimized distributed exchange FFT remain separate extensions. The next active group is5, then3; do not manufacture new storage formats merely to extend group4. Keep authorized validation results in calc2026.

### Remaining groups (not a count of individual implementation patches)

Group5 now has concrete R0–R4 scope in [resource inventory](2026-11-20-group5-resource-inventory.md). R0 scalar IEEE compliance is implemented and GNU-tested (Fujitsu untested). R0–R2 are implemented and verified within documented profiles; remaining group5 milestones are the still-unhandled R3 wider semantic/numerical-init failures (including allocation failures outside the verified coordinate-return path) and R4 native ABI stability. R3 preflight, parsed admission, detected namelist read errors and early eh_init grid/boundary/time/CFL returns are verified; scheduled-job startup declaration errors now record and exit collectively. R4 declared contracts and explicit isolated protocol probes are implemented. C/ISO_C_BINDING Fortran host state lifecycles are verified for standard PZ and fast-PR HSE full/reduced-k on Mac MPI1/2, including independent-state controls and rejected-call state preservation. Stable ABI negotiation, broader profiles/platforms and the remaining host operations are still open; lifecycle success is not blanket ABI certification. Management-node checking remains an optional independent user action. Do not treat every process-scoped selector as a new isolation task.

1. [ ] DEFERRED by user: extend measured-norm diagnostics to fast PR HSE. Isolated E-drive norm
   wiring belongs to paused molecular work and is not the next active priority.
2. [ ] DEFERRED by user: establish an actual current-grid campaign with two SCF-converged grids
   passing5%, then separate k-point/response checks. Add unambiguous PP cutoff
   extraction; explicit species/channel radius is required until that exists.
3. [ ] Certify the new maintained runtime's NLCC-force capability with appropriate
   physical evidence; audit legacy SALMON eV/Hartree versus Python CODATA without
   silently changing existing results.
4. [ ] Dense ordinary-XC and fast-PR HSE full/reduced-k orbital/spatial direct
   SCF repartition is measured below. Remaining: packed/tiled/sparse state
   migration, arbitrary changed-numerics source-version transitions and cross-platform/portable-format migration. Same-source rebuilt-runtime and exact declared metadata-only source-pair
   SCF restart are measured below with explicit native-format/codec/state-revision proof. Current
   built-in SCF generates dense ACE; do not add an artificial format switch only
   to fill a checkpoint checklist. Same-topology RT remains complete in scope.
5. [ ] Retained GS timers are now owner-scoped and independently queryable;
   transient DOS/PDOS/band/AVS file units and model-derived output conversions
   are isolated below. Retained GS now uses a native path scope without host chdir.
   Standard/shared/registered TDDFT, fast-PR HSE and FDTD now require directory
   protocol3 through one common helper, including GS→RT setup, without host chdir.
   HSE initialization settings are explicitly caller-owned; the implicit fallback
   has been removed.
   Remaining initialization/
   other I/O/MPI ownership (native unit30 fixed; propagation MPI child ownership
   is now state-scoped below). Standard TDDFT/fast-PR HSE singleton and
   handle timers are owner-scoped. Remaining stable public native
   ABI and wider initialization/allocation/process-failure recovery. Existing
   rollback does not guarantee arbitrary STOP/OOM/MPI-failure or threaded reentrancy.
6. [ ] Larger grids/material counts, resource/communication calibration, matching
   runtime distribution and Linux/Intel Mac/other-platform evidence. Respect the
   Actions budget; no routine CI dispatch.

Latest restart evidence: calc2026/salmonpy-scf-declared-record-20261120,
calc2026/salmonpy-scf-rebuild-record-20261120,
calc2026/salmonpy-hse-korspace-record-20261120,
calc2026/salmonpy-hse-symspace-record-20261120,
calc2026/salmonpy-hse-kspace-record-20261120,
calc2026/salmonpy-hse-mixed-record-20261120,
calc2026/salmonpy-scf-mixed-record-20261120 and
calc2026/salmonpy-hse-spatial-record-20261120. Convergence evidence:
calc2026/salmonpy-tddft-measured-norm-20261119 and
calc2026/salmonpy-current-grid-criterion-20261119. Local writes are verified;
cloud synchronization is not confirmed. Detailed historical records follow.

## Historical queue — 2026-11-18 (after functional test cleanup)

This was the active queue at that checkpoint; the current queue above supersedes it.
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
  diagnostic records updated. Current source/installed-wheel regression is628
  (one skip each); latest standard SCF-progress slice preserves this method count.
  Historical620/622 counts below describe their dated validation slices.

### Next priority — source consolidation and SALMON convergence (2026-11-19)

User direction: fast-PR HSE is the maintained solver basis. Audit effective
source/feature differences and migrate standard-only interfaces before removing
duplicate maintenance trees. Reduce archive-helper fragmentation and curate three
beginner examples while preserving public API. Design only so far; no deletion.

User convergence recipe: PP radius/3 → about0.1 angstrom with few fixed k points,
then a few-fs impulse dt scan from about0.2au downward. Stop a NaN/Inf trial
immediately; report nonbreaking and value-converged dt separately. k convergence
is a later independent check; SCF/resource warnings do not halt the campaign.
See 2026-11-19-source-consolidation-and-convergence-design.md for acceptance and
unresolved radius/tolerance choices. Full SCF continuation is deferred behind
this priority; previously verified pending changes are retained.

### Remaining work (do not double-count historical notes)

1. [x] Verified external-k import, declared input metadata, portable retained-GS
   load and source-labelled re-export committed/pushed as f6005fb9 [skip ci].
   Subsequent work remains uncommitted until requested; no routine CI dispatch.
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


## 2026-11-19 — collective retained-GS array admission

- [x] Invalid quantity types and differing supported quantities agree across the
  subgroup before native access; fixed requests can retry without owner changes.
- [x] Existing unit contract extended (no new test method); real synthetic MPI4
  PZ orbital2 and Libxc-PBE spatial2 preserve GS bytes/resources and active RT,
  match independent currents and direct resume without replay or new SCF.
- [x] Paired guide/AI JSON synchronized; historical completed slices remain below
  their dates and must not be counted as new outstanding work.
- [ ] Full SCF iteration continuation remains: retained mixing storage exists,
  but current model.run refuses active storage and always enters initialization.
  Restart must rebuild resources and restore explicit iteration/mixing controls.
  See 2026-11-19-scf-checkpoint-followup.md; this fix does not implement it.

- Final check:628 unit methods (one skip), existing SHA-verified wheel LD/SBE
  kernels reused without rebuild. Synthetic MPI4 orbital-PZ/spatial-Libxc-PBE
  admission/retry/current/direct-resume controls pass; orbital current hash equals
  f6005fb9's measured baseline. Docs/AI sync/diff checks pass. Evidence:
  calc2026/salmonpy-gs-array-admission-20261119 (logs/hashes only, no GS arrays).
  No Actions or new SCF. Local storage verified; cloud synchronization unverified.


## 2026-11-19 — owner-scoped SCF convergence progress (standard)

- [x] Move previous density/potential arrays, normalization, no-mixing counter and
  convergence flag into s_mixing; capture final iterations/residual/rion_update.
  Existing fresh-entry reset, arithmetic and native convergence policy retained.
- [x] Complete dependent-object rebuild and source/object prohibition audits;
  synthetic ownership/copy/reset/bounds and extracted five-criterion convergence
  bit parity in MPI1/2 x OMP1/2; iteration body unchanged beyond ownership/capture.
- [ ] Continue active SCF without reinitialization; serialize complete numerical
  context and rebuild resources in independent restore candidates. This slice
  does not implement a complete SCF checkpoint or public resume entry.

- Final progress slice: source628/installed628 (one skip each); installed native
  MPI4 import/current/direct-resume and installed ownership/extracted convergence
  runner pass. PZ orbital2/Libxc-PBE spatial2 current hashes equal prior runtime.
  Runtime/standalone probe object audits and offline sdist→wheel/source inventories
  pass. Current optional RT-step profile synchronized to the new base/type layout.
  Evidence: calc2026/salmonpy-scf-progress-20261119, no GS arrays/new SCF/CI.
  New shared-runtime SHA6950875d; GNU/macOS only, not Fujitsu/platform validation.


## 2026-11-19 — consolidation first batch

- [x] Effective-source inventory validates ordered complete-file replacements and
  preserves per-file origins/hashes without mutation. Standard670/fast-HSE725
  fully match actual prepared trees:573 same,77 different,20/75 unique paths.
- [x] Consolidate three GS archive implementation modules into _gs_archive.py,
  retaining public imports and separate GS→TDDFT lifecycle. Private module count
  reduced by2; function ASTs unchanged except private naming/self-import.
- [x] Beginner examples: LD→synthetic SBE→Si execution check; JP/EN/AI synchronized.
- [x] One maintained fast-PR-capable solver tree, compatible no-HSE/HSE profiles (completed d86d5692):
  migrate standard-only owners/bridges and reconcile77 effective differences.
  Current trees/references/fixed binaries preserved; no blanket copy/delete.
- [x] SALMON-specific grid/dt planning, reporting and synthetic campaign tests (first callback/native-diagnostics stage);
  PP radius definition/explicit tolerance contract per approved design.

- Final first-batch verification: source628/installed628 test methods (one skip
  each), source and installed real-MPI14 archive cases pass. Native PZ orbital2
  current hash/direct resume unchanged; public LD/SBE examples pass, Si CLI help
  checked without a new SCF run. Function AST parity and source/installed/prepared
  inventory equality pass. Offline sdist→wheel, native checksum and JP/EN/AI sync
  verified. Current Python module count68/private17 (was70/19).
  Evidence: calc2026/salmonpy-consolidation-first-batch-20261119, no GS collection
  or Actions. Native single-source migration/convergence execution still pending.


## 2026-11-19 — fast-source SBE sharing and HSE-off build path

- [x] Fast explicit routing references existing standard SBE types/kernel/wrappers
  and CMake list; no new file copies. Native/solver byte equality checked by a
  common preparation helper; history/reference-only trees unchanged.
- [x] HSE-only Taylor/time-evolution recovery/history references guarded; ordinary
  numerical checks preserved. HSE-preprocessed statements equal the previous
  version except NVTX source-line labels; no-HSE preprocessor has no EXX calls.
- [x] Same final source builds HSE ON and OFF. OFF omits HSE FFTW configuration
  and linkage. Source/object prohibition audits retained.
- [x] Full maintained-source unification (completed d86d5692; following counts historical):18 standard-only paths/74 differing
  effective paths remain; default GS/TDDFT/EH preparation is not switched.

- Final shared-SBE slice: same-source HSE ON/OFF builds audited; no-HSE executable
  has no FFTW linkage.16 SBE configurations in each of HSE/OMP1,HSE/OMP2,no-HSE
  and HSE/MPI2 compare to the unchanged validated old executable and Python.
  OMP1 files bit-exact; OMP2 within3e-13 relative/3e-15 absolute tolerance.
  HSE Taylor/PT-CN learned full/reduced-k internal-MPI2 independent current/
  energy/history and direct no-replay resume are exact to old controls.
  Source628/installed628 methods (one skip each), installed preparation source
  hashes/native MPI/SBE example, offline sdist→wheel and JP/EN/AI sync pass.
  Evidence: calc2026/salmonpy-hse-shared-sbe-20261119 (no GS arrays/new SCF/CI).
  Full source unification/convergence driver and other platforms remain open.


## 2026-11-19 — fast-source FDTD/LD/EH sharing

- [x] Six existing standard Maxwell paths now serve fast explicit routing too:
  CMake,LD,EH kernel,EH calculation,main and embedded E/H. No extra source copy.
- [x] Fast build publishes field_abi only after probing the compiled protocol;
  NativeEHField accepts that record while retaining legacy embedded records.
  Separate executable-build.json binds batch/MPI FDTD executable SHA.
- [x] HSE ON/OFF builds;76 observation files across single/multiple LD/dielectric/
  magnetic/ID/order/axis/impulse cases match the current standard reference bytes.
  MPI4 multi-LD HSE/on-off,magnetic rank corners and native-EH MPI2 failure/retry
  controls pass;A/E/H equal the standard library.
- [x] Same library serves separate field/material groups in the ordinary mixed
  sample:local LD,SBE,two full-k HSE model inputs,128 continuous steps versus
  old standard field,and direct8→128 resume match with zero replay.
- [ ] Standard GS/input/TDDFT source migration and remaining71 differing paths;
  15 standard-only paths include five alternative HSE files NOT to migrate,
  plus bridge/output/audit functions already present elsewhere. Counts are not
  a list of15 missing physics features.


## 2026-11-19 — shared IO prerequisites for GS/ordinary TDDFT

- [x] Fast preparation shares the existing standard IO list, model-input value
  storage and RT stream management. Remove identical adapter step_output source;
  builders reuse solver object or compile the same source for compatible older
  builds. Include each linked application object in the audit.
- [x] Native model-input capture/activate/deep-copy/unallocated-array controls pass.
  HSE's additional complete settings snapshot remains separate and unchanged.
- [x] Same prepared source builds fast HSE with host MPI, standalone serial, and
  MPI without HSE. Fix unconditional MPI-only symmetry/child-owner references in
  serial paths without changing MPI numerical expressions or field validation.
- [x] Two different full/reduced HSE model inputs under internal MPI2 match old
  independent current/energy/history and direct resume; ordinary four-rank
  LD/SBE/two-model mixed128 and direct8→128 controls match with zero replay.
  Serial Taylor/PT-CN worker current/energy and direct restart match old controls.
- [ ] Migrate parser workspace/GS owning handles/ordinary TDDFT into the fast
  tree; sharing an input value record does not implement those lifecycles.
  Current inventory587 same/70 differing/13 standard-only/75 fast-only; five
  alternative standard HSE modules remain excluded from porting.

No new SCF/GS collection or Actions; only synthetic regression fixtures used.

Final IO slice: source628/installed628 (one skip each), sdist→wheel, installed
source-hash reproduction and standalone worker linking checked. Validation:
calc2026/salmonpy-hse-shared-io-20261119. No new Python module/test method.


## 2026-11-19 — fast model-owned parser workspace

- [x] Port standard workspace ownership to fast inputoutput: scratch units,
  coordinate-source flags, namelist diagnostics, common defaults/value records.
  Registered HSE models use an isolated parser; legacy batch keeps named files.
- [x] Preserve complete fast defaults, unit initializer, KY spin resolution and
  input validation; keep deprecated local MLWF aliases explicitly initialized.
  Broadcast isolated parser metadata on all ranks; keep file handles root-local.
- [x] Actual MPI2 A→B→A, reused workspace, split/interleaved scratch inputs,
  omitted defaults, mixed units/coordinates and previous temporary-file controls
  pass. All variables.log entries equal the unchanged old fast library; standard
  parser regression still passes. HSE-specific inverse-length/GPU-batch inputs
  survive reading and reset to native defaults when omitted.
- [x] MPI2 heterogeneous HSE registry (six ranks total) current/energy/history,
  direct resume and fresh independent state match old controls. LD/SBE/two-model
  four-rank mixed128 and direct8→128 match with zero replay.
- [x] Same source HSE-on MPI/serial and HSE-off MPI builds pass source/object
  prohibition audits. Serial Taylor/PT-CN24 current/energy/direct resume match;
  no-HSE SBE16 full-application configurations match byte-for-byte.
- [ ] Retained GS numerical owners and ordinary embedded TDDFT lifecycle in the
  fast tree. Common parser values do not replace the complete HSE input snapshot.
  Global activation is still serialized; do not infer threaded native reentrancy.

Source628 methods pass (one skip). No new Python module/test method or new GS/SCF.
Validation records:calc2026/salmonpy-hse-owned-parser-20261119. No Actions used.

Final parser slice: installed628 methods (one skip), offline sdist→wheel,
installed source-hash reproduction, installed MPI2 parser and MPI6 registry
controls pass. JP/EN/AI snapshots synchronized. Local result storage verified;
OneDrive cloud synchronization and other platforms/compilers not confirmed.


## 2026-11-19 — retained fast GS import and common seed

- [x] Fast GS context import-only entry reads existing compatible wavefunction
  restart without SCF/warmup/GS checkpoint output. Skip irrelevant SCF checks,
  resolve unused initial-no-mix sentinel locally, keep buffer/mixing validation.
- [x] Context-owned input snapshot, MPI child ownership, cache release and detached
  common gs_seed export. No HSE plan/pointer exported; common seed is not a full
  model/checkpoint. Share existing standard gs_seed source, no extra kernel copy.
- [x] Internal MPI2/serial HSE four XC pairs, including two HSE states with
  different screening/coordinates; no-HSE MPI2 two non-HSE pairs. Finite orbitals,
  input/array independence, fresh import equality, close isolation, no child leaks.
  Existing restart files remain unchanged. All tests use synthetic fixtures.
- [x] Common real/complex seed redistribution with empty bands passes in MPI2.
  Normal SCF iteration/output body byte-exact to preceding source; new SCF runs0.
  HSE registry independent current/energy/history and four-rank mixed128/direct
  save8→resume120 match old controls with zero replay. Source628 methods pass,
  one skip; existing test methods/files extended rather than multiplied.
- [ ] Python fast GS handles/NativeGSBackend connection; ordinary embedded TDDFT
  initialization from the retained seed and complete model/resource lifecycle.
  Native Fortran import alone does not implement those public workflows.

Current inventory588 identical/70 different/12 standard-only/75 fast-only paths
(fast733). Five alternative standard HSE modules are not migration targets.
Evidence:calc2026/salmonpy-hse-gs-import-20261119. No new GS collection/SCF/Actions.

Final GS-import distribution slice: installed628 methods pass (one skip),
sdist→wheel/source inventory reproduction and installed MPI2 four-pair GS import
controls pass. JP/EN/AI guide snapshots synchronized. Probe legacy compilation
is retained for earlier builds; real/complex seed redistribution also passes.
Normal fast GS omitted no-mixing sentinel resolution remains to audit with the
public GS lifecycle migration; this import slice resolves that unused sentinel
locally and does not run a new normal SCF. Local evidence only; cloud sync and
other platforms/compilers unverified.


## 2026-11-19 — Python fast GS import handles

- [x] Existing GS Python API accepts fast host build records explicitly declaring
  gs_abi/1 import mode. Compiled protocol probe required; no implicit build/load
  during import or backend construction. Standard GS behavior retained.
- [x] Native stable pointer slots and monotonic tokens: import/info/array/close/
  resources. Density/occupations/orbitals are copied read-only; metadata not
  established by import is None. Computed-only requests reject before native use.
- [x] Registered GS reads switch/restore symmetry scratch, preserve sym.dat and
  normalize unused RT-learning/pre-SCF controls without changing XC/screening.
  Source GS and active retained owners remain unchanged.
- [x] MPI6/internal2 full/full and full/reduced GS arrays match independent groups;
  closing/fresh import, stale/mixed tokens and zero child-resource leaks pass.
  Different HSE/Libxc-PBE and TB-mBJ/Libxc-PZ inputs pass the same controls.
  Existing retained_gs sample import mode runs in MPI2 with preserved IDs/nulls.
- [ ] Fast GS SCF creation, energy/force/eigenvalue computation and portable export;
  ordinary-XC GS→TDDFT initialization remains to connect; HSE06 memory initialization
  is completed in the following slice. Current fast GS
  public path imports existing compatible scalar periodic wavefunctions only.

No new Python module/test method/sample file, no new SCF or Actions. Evidence:
calc2026/salmonpy-python-fast-gs-20261119. Remaining standard/fast source migration
and approved convergence workflow are not completed by this handle connection.

Final Python fast-GS slice: source628/installed628 methods (one skip each),
installed MPI6 full/reduced ownership and MPI2 sample, source-hash reproduction,
sdist→wheel and JP/EN/AI synchronization pass. Same-source no-HSE MPI and serial
HSE builds also pass. New SCF0; no GS datasets or Actions. No extra Python module,
sample file or test method. The subsequent HSE06 GS→TDDFT connection is recorded
below; ordinary-XC RT and computed GS lifecycle migration remain pending.


## 2026-11-19 — retained GS to fast-PR HSE memory initialization

- [x] Optional value seed through RT initialization; HSE models retain independent
  orbitals/occupations for later states, without GS restart files, SCF or RT replay.
  Static potential/XC/exchange initialization retains the existing update order.
- [x] Existing NativeGSState.tddft_model connects HSE06/atomic-unit scalar periodic
  imports to the fast host registry; k-only internal MPI1/2. Initial bootstrap
  cleanup, alternate dt/capacity and source-close lifetimes are explicit.
- [x] Sequential subsequent GS import parks/restores active RT settings and
  symmetry workspaces. Common timer initialization now shares standard's existing
  once-only allocation guard; old normal initialization behavior is preserved.
- [x] Synthetic full/reduced learned Taylor/PT-CN memory/file current and energy
  comparison; source GS close/removal, independent sibling/fresh states and exact
  direct resume after8 with no replay. Normal HSE registry current/energy/history
  remains exact against the preceding fixed runtime in internalMPI2 controls.
- [x] Existing retained_gs sample gains optional --current-steps; JP/EN and AI task
  catalog document scope. No new Python module, test method or sample file.
- [x] Ordinary PZ/TB-mBJ/Libxc RT migrated in the HSE-disabled fast profile below.
- [ ] Fast computed GS SCF/energy/forces/eigenvalues/export and final source-tree merge.
  Existing standard-runtime features remain supported. This HSE bridge alone does
  not complete the single-tree migration or approved convergence workflow.

Evidence: ../fast-gs-memory-rt-validation-2026-11-19.md. No fresh GS/SCF collection
or GitHub Actions. Same-source comparison inventory is now589 identical,
69 different,12 standard-only,75 fast-only (standard670/effectivefast733).


Final memory-RT slice: source/installed628 (one skip), installed MPI4 controls,
MPI2 existing sample, MPI6 old-runtime registry/GS-import comparisons pass.
Same-source HSE-MPI/no-HSE-serial builds/audits and sdist→wheel/source reproduction
pass. Installed mixed128 matches the fixed prior baseline exactly; save8→resume120
also matches continuous128 with replay0. These are local synthetic controls.


## 2026-11-19 — ordinary XC on the fast source profile

- [x] Share embedded TDDFT, endpoint worker and MPI bridge source with the fast
  prepared tree; no copied second TDDFT adapter. Compile-time source-profile branch
  consumes fast value-owned GS seeds and omits obsolete global scratch switching.
- [x] USE_HSE=OFF host build advertises source_profile=fast-pr, ordinary_rt_abi=1
  and import-only gs_abi/1. Existing NativeGSBackend.load/tddft_model and ordinary
  NativeTDDFTBackend/registry APIs are connected. No SCF or GS-file output.
- [x] PZ/TB-mBJ/Libxc-PZ/PBE, k/orbital/spatial innerMPI2: current parity against
  fixed preceding standard runtime, exact sibling/fresh/direct-resume checks.
  PZ/TB-mBJ current difference0; Libxc at most5.421010862427522e-20au in these controls.
- [x] Different-XC GS import preserves active RT; explicit external k assets,
  missing-asset admission and optional input metadata pass. File-registry tests
  cover different XC/dt/capacity/MPI partitions and56 accepted steps, exact against
  independent owners, with direct restart and stale/mixed-handle rejection.
- [x] File initialization captures API dt/capacity/MPI values after worker_configure;
  the first mixed-input file test reproduced old parser-snapshot order failure.
- [x] Same-source fast-HSE MPI build and fixed prior HSE registry comparison retain
  exact current/energy/history/restart. Standard MPI build/memory connection also
  remains exact against its prior fixed binary. JP/EN/AI and existing sample updated.
- [ ] Computed GS (SCF/energy/eigenvalues/forces/export) migration and final tree
  consolidation. One-session ordinary/HSE registry mixing remains separate; normal
  mixtures can still use separate material MPI groups. Physical convergence remains
  to validate; no new GS dataset collection or Actions is part of this slice.

Current source comparison:592 identical,69 different,9 standard-only,75 fast-only;
standard670, effectivefast736. These are source-provenance counts, not physics claims.
See ../ordinary-fast-rt-validation-2026-11-19.md for controls and measured limits.


Final ordinary-RT slice: source/installed628 (627 pass,1 skip); sdist→wheel and
installed-source261-file parity; installed GS/PBE current and spatialMPI2 resume
controls pass. FDTD/LD/SBE/PZ128 vs independent serial control passes; direct
63→65 restarts match all128 A/E/j frames exactly. Serial communication cleanup
and MPI-only GS compilation are corrected without enabling serial host-MPI APIs.
Different cell/XC model registry ownership remains exact. No new SCF or Actions.

Final admission check: HSE06 on the HSE-disabled frontend rejects collectively
before native initialization and preserves GS/RT resources/current. Ordinary RT
XC checks precede the worker call; valid GS arrays remain usable. Coverage is
folded into existing unittest/MPI tests; test-method count remains628.

### Computed ordinary GS on the fast source (2026-11-19)

- [x] HSE-disabled run capability, native convergence metadata and five array kinds.
- [x] Orbital/spatial MPI2 retained A+B versus independent A/B: exact arrays; close/reuse isolation.
- [x] Nonconvergence warns; large-threshold control reports native convergence.
- [x] Export assembly checked in RAM; no GS archive/restart dataset collected.
- [x] Previous fast executable PZ/PBE scalar/text regression matches. Old standard unconverged PBE differs by 3.4e-7 Ha; do not relax its comparison.
- [x] Fresh HSE-enabled build regression after the GS mixing metadata addition.
- [x] Computed GS to ordinary RT lifecycle regression.
- [x] Installed-wheel validation: sdist→wheel; source/installed628 (627 pass,1 skip); computed GS→current sample and spatial MPI6 direct resume exact.
- [x] HSE-enabled computed GS (new declared run capability, full-k/k-only native profile).
- [x] Final source-tree consolidation (2026-11-19 unified-source stage).

### Computed GS→ordinary RT continuation

- [x] Six-rank orbital/spatial controls: computed PZ/Libxc-PZ and PBE seed ownership, interleaved SCF, source GS close, independent current equality.
- [x] Direct 8→16 restart exact, no replay; transient RT checkpoint files removed.
- [x] Existing retained_gs sample now supports --current-steps in computed mode, with three-component A drive; MPI2 execution passes.
- [x] Fresh HSE-enabled build and all-object audit; full learned Taylor/reduced learned PT-CN MPI6 currents, energy and history match the fixed previous build exactly; direct restart exact.

### Computed fast-PR HSE GS (2026-11-19)

- [x] HSE06 DFT descriptor and declared native run capability; conventional-HSE fallback and unconnected batch HSE reports are not used.
- [x] Explicit exx_pre_scf_threshold=0 in new HSE descriptors: run HSE itself, not unfinished PBE warmup. Nonconvergence warns.
- [x] Darwin fast-HSE builder shares ordinary mpi_link_command; avoids NumPy Accelerate ZDOTC complex-return ABI collision. Failed release/debug binaries and logs retained.
- [x] Fresh build, source/all-object audits, MPI6 full-k different-input retained owners, direct8→16 RT resume, interleaved SCF and output assembly in RAM. All independent arrays/currents exact.
- [x] HSE GS info/eigen text output matches prior fixed fast standalone executable byte-for-byte. No GS restart/array dataset collected.
- [x] HSE computed sample --xc hse06 --current-steps16 passes MPI2. Ordinary Libxc-PZ/PBE spatialMPI6 remains exact.
- [x] Installed wheel/sampleMPI2/HSE MPI6 direct restart, source/installed628 (627 pass,1 skip), JP/EN/AI/tool byte parity; old-HSE RT with matched linking is exact.
- [x] Final current source-tree consolidation and batch-report migration (2026-11-19).

### Fast GS batch reports and optional restart output (2026-11-19)

- [x] Fast-source fixed-ion non-DC GS status/force report migration, preserving numerical operation order.
- [x] Public PZ/Libxc-PBE batch runs without restart arrays; fixed-fast energy difference0 and force print-precision comparison; SCF nonconvergence warns.
- [x] Per-run save_restart=False, actual restart_saved metadata, backward saved=True default, warning-only unexpected output; clear missing-array admission to restart consumers.
- [x] Existing lifecycle --batch controls and dft_si --no-restart/--max-iterations sample; source/installed8-atom4-iteration execution, no GS datasets collected.
- [x] Fresh fast serial/no-HSE and MPI/HSE builds, source/object audits, HSE CLI report comparison, MPI6 retained GS→RT/direct restart and full/reduced learned-history regression.
- [x] Source/installed628 tests (627pass,1skip), sdist→wheel, JP/EN/AI synchronization; no new unittest methods.
- [x] Batch HSE fast-PR provenance connection and final maintained source-tree consolidation.
- [ ] Audit legacy SALMON eV/Hartree versus modern Python CODATA conventions without silently changing existing calculations.

### NLCC parity gate before solver-source retirement (2026-11-19)

- [x] Reproduce missing fast Libxc NLCC energy weighting (~0.122 derivative mismatch); port only the missing weight, retaining fast state/spin/gradient implementation.
- [x] Port sampled-core force using owned species/PP boundary/Vxc, V1 cutoff and V2 full support/skew/periodic images; forward Vxc at all five callers.
- [x] Scalar IEEE used-entry validation, V2 unused-padding acceptance; finite/NaN/Inf/overflow controls.
- [x] Prescribed-density variation/force OMP1/2; MPI2/4 spatial + orbital-replica rejection. V2 compensation reduces1.82e-13 to2.50e-17 without loosening MPI1e-13 limit.
- [x] No-NLCC PZ/PBE/HSE parity, independent GS/RT/direct restart, full/reduced learned/PT-CN history regression; fresh HSE-off/on source/object-audited builds.
- [x] Builder adapter sources select prepared-tree files, not a hard-coded standard checkout. Preparation-level shared-file references remain.
- [x] Source/installed628 (627pass,1skip), sdist→wheel, installed HSE MPI6, installed preparation736-file byte parity and JP/EN/AI synchronization.
- [x] Default fast-source preparation/provenance, canonical source placement and retirement of duplicate current trees; synthetic archive converter gate passes.
- [ ] New-source runtime NLCC capability certification remains a separate gate.

### Unified maintained source and default build (2026-11-19)

- [x] One current736-file solver tree; historical sources/overrides preserved as checked archives, excluded from distributions.
- [x] Default and HSE preparation share verified copies; default public build disables HSE.
- [x] Fast-specific31 initialization inputs captured/activated per model; native ownership probe passes.
- [x] Fresh default serial, ordinary MPI and fast HSE MPI builds; source/all-object audits retained.
- [x] Six ordinary-XC MPI cases exact against prior fast build; HSE full/reduced learned Taylor/PT-CN exact, direct resume no replay.
- [x] HSE computed GS independent ownership and16-step direct resume exact; public GS batch scalar comparison passes.
- [ ] Broader converged-material validation, new-source NLCC runtime certification and other-platform execution remain separate gates.

- [x] Installed default-source preparation,736-file byte parity; wheel excludes duplicate historical archives.
- [x] Source/installed628 tests (627pass,1skip), prescribed NLCC OMP/MPI2/4, archive converter synthetic I/O and host-MPI bridge regressions.
- [x] Restore existing process-worker parser/endpoint wiring in the unified source; keep fast HSE numerical paths.

Final closure: source/installed628 (627pass,1skip), installed736-file and JP/EN/AI byte parity; default serial and ordinary/HSE MPI audits/builds; ordinary six-case currents and HSE learned/history/direct restart exact. Installed mixed128-step serial control and direct63→128 A/E/j resume match exactly without replay. No GS dataset collection or CI. Evidence: calc2026/salmonpy-unified-source-20261119.

### Norm-based convergence campaigns (2026-11-19)

- [x] Unified-source commit d86d5692 pushed to independent SALMONpy branch; original history unchanged. OneDrive Git mmap issue avoided through verified local object transfer; no CI.
- [x] Grid radius/3→0.1 angstrom plans, ceiling/deduplicated grids, explicit radius definition/species/PP hashes; finer limits allowed.
- [x] Independent dt trials from0.2au, fixed impulse and common physical duration; norm tolerance default1e-5 after user clarification.
- [x] Fortran scalar-IEEE sample checks, norm loss/gain/error and optional trace comparisons. Python handles planning/dispatch/records.
- [x] Tolerance exceedance reports without stopping; numerical breakdown stops only that trial; runtime/cleanup failure separated. SCF nonconvergence warns/provisional.
- [x] Wire ordinary built-in TDDFT electron-number measurement and memory GS/RT callbacks; execute native MPI2 Si campaigns without collecting GS datasets. Fast-HSE/isolated measurement remain separate.
- [ ] Automatic PP cutoff-radius selection remains deferred: require an explicit physical radius rather than guess channel semantics.

Campaign verification: source and installed631 tests (630pass,1skip); final provenance update covered by3 affected tests in source/final wheel. Installed planning sample and Fortran below/above1e-5 loss/gain/NaN/Inf controls pass; old LD/SBE regressions remain passing. JP/EN/AI byte parity and sdist→wheel verified. No new GS/native TDDFT norm measurement campaign or CI. Evidence: calc2026/salmonpy-convergence-campaigns-20261119.

### Measured TDDFT norm and early grid stopping (2026-11-19)

- [x] Commit/push5ba6406d before continuing; skip-ci retained.
- [x] Add read-only Fortran occupied-density integral to ordinary native lifetime and per-state-handle C ABI; collective internal-MPI calls. Process worker optional norm payload preserves old protocol.
- [x] Default run_dt_scan reads state.electron_number; no configured-nelec substitution. Old runtime normal use remains valid, missing norm reports explicitly.
- [x] Confirmed adjacent SCF-converged grid pair ends the scan before0.1A; full-scan override and provisional evidence retained.
- [x] Real MPI2 Si8/grid12 andgrid16/2fs trials. Grid16 GS converged,dt0.2 fails step8; finer trials continue;dt0.1 max relative electron loss4.22e-10, meeting1e-5. Grid12 SCF unconverged stays provisional.
- [x] Recover the existing gross-electron-count STOP as status21 in recovery-capable native stepping; thresholds unchanged, accepted state restored, failing trace excluded. Raw batch behavior unchanged.
- [x] Six XC/decomposition controls: measured norm/query isolation/direct restart and independent current equality to fixed prior build; process current/replay norm exact.
- [ ] Fast PR HSE/isolated E-drive measured-norm API, PP cutoff parser and broader physical/grid/k/response validation remain.

Final measured-norm verification: source/installed632 tests (631pass,1skip), sdist→wheel and installed native-MPI example pass; JP/EN/AI and implementation byte parity. Si16³ GS converged, ~2fs:dt0.2 returns numerical error at step8 while remaining trials continue;dt0.1 max electron-loss fraction4.22e-10. No GS datasets or CI. Canonical evidence: calc2026/salmonpy-tddft-measured-norm-20261119.

### Current-based grid criterion (user clarification)

- [x] Primary grid API uses each grid's norm-qualified dt trace; default relative current RMS tolerance5%, finer-grid reference and common physical times. Fortran computes errors.
- [x] Confirmed SCF-converged passing pair stops before0.1A; scalar-energy API remains supplementary, unresolved dt excluded.
- [x] Synthetic10% rejection / passing pair / inclusive5% / zero-crossing / zero-reference controls. Existing Si12³→16³ trace comparison18.83%, not passing and provisional due coarse GS.
- [ ] A physically converged current-grid campaign with two SCF-converged grids meeting5% is not established by the existing pair.

Current-grid criterion verification: source/installed633 tests (632pass,1skip), sdist→wheel and JP/EN/AI/implementation byte parity. Error statistics use Fortran only. Existing Si12³/16³ waveforms differ by18.83%, provisional and not accepted at5%. No new GS data or CI. Evidence: calc2026/salmonpy-current-grid-criterion-20261119.

### Group4 first step: shared-source SCF progress ownership

- [x] Reproduce missing owned-progress type in current fast-source modules with the existing native ownership probe before implementation.
- [x] Restore previous density/local potential, electron normalization, convergence/no-mix counts and iteration/residual/ion-update capture into s_mixing%scf_progress. Preserve existing mixing%scf_converged consumers.
- [x] Entire iteration body identical after reversing ownership substitutions; scalar arithmetic/order unchanged. Synthetic MPI1/2×OMP1/2 convergence and copy-isolation probes pass.
- [x] Computed PZ/PBE spatialMPI2 GS/independent arrays/16-step RT and direct resume equal the unchanged preceding runtime; no GS dataset collection.
- [ ] At this ownership-only checkpoint full active-SCF continuation was pending; the in-memory stage is implemented below. Portable checkpoint/version and changed decomposition remain open.

Shared-source SCF progress validation: normal and fast-HSE builds source/object-audited; native progress/convergence MPI1/2×OMP1/2 and installed probes pass. PZ/PBE and fast-HSE computed GS→16-step independent/direct-resume controls exact; HSE learned/PT-CN history exact. Source/installed633 tests (632pass,1skip), sdist→wheel and JP/EN/AI byte parity. Full active-SCF continuation remains next. No GS datasets or CI. Evidence: calc2026/salmonpy-scf-progress-shared-source-20261119.

### Group4: active SCF memory snapshot and continuation

- [x] First commit/pushaba74f9b; no CI.
- [x] Pause accepted SCF boundaries; preserve convergence/no-mix progress and mixing history; continue without init or replay.
- [x] Native numerical payload: orbital/CG/occupations/mixing, density/spin/potentials, Poisson/reciprocal/XC and effective post-init inputs; static handles remain live.
- [x] Distinguish declared GS→RT input from effective SCF continuation settings. Fix interleaved-model force mismatch from derived Ewald cutoff input restoration.
- [x] Six PZ/Libxc-PZ/PBE MPI1/2 orbital/spatial/k controls, no-mix0/3, same-owner snapshot3/5 and mixed-token/budget rejection. Numeric arrays/info match continuous new and fixed previous runtime exactly.
- [x] Linear8 steps: initialization1/executed8. Restore executes0 iterations; intentional future branches are counted.
- [ ] Remaining numerical coverage: symmetry/SOC/+U/DC/ionic-motion continuation; portable disk/process/checkpoint version and changed-decomposition restore remain group4 work. Fast-PR HSE full-k memory continuation is completed below. No completed disk checkpoint claimed.

Active SCF validation: source/installed634 tests (633pass,1skip). Six native ordinary-XC controls exact; PBE kMPI2 continued GS→16-step TDDFT currents equal both continuous new and fixed preceding runtime. Fresh-path arithmetic/copy probes pass. Fast-HSE computed GS→RT and learned/PT-CN history regression exact; this does not enable paused HSE SCF. Source and wheel example restore with initialization1/executed8 unchanged by restore. GS data not collected, CI not invoked. Evidence: calc2026/salmonpy-scf-live-resume-20261119.

### Group4: fast-PR HSE live SCF continuation

- [x] Commit/push previous ordinary-XC work as78c5bcd1 with skip-ci.
- [x] Keep explicit HSE numerical ACE/cache/source/action/occupation/history/support flags independent in SCF snapshots. Never copy opaque FFT/backend owners; rebind the restored system to its live HSE owner.
- [x] Preserve mixing age and pre-SCF first iteration/count/phase across segments. No GS initialization, no-mix replay or exchange evaluation in restore.
- [x] Actual MPI1/2 full-k and Gamma controls, interleaved different materials, snapshot3/5, continuous new/fixed preceding arrays/forces/info exact. Sixteen TDDFT currents from continued GS exact too.
- [x] PBE→HSE switch across cuts exact against fixed continuous runtime with a control-flow tolerance. Tight-threshold target HSE warns/retains (no false convergence) against the new continuous runtime; preceding binary rejects this case and is not claimed to pass it.
- [ ] Same-process memory scope only. Symmetry/SOC/+U/DC/ionic motion, portable disk/process restart, topology change and common direct integration completed below; candidate integration follows below.

No GS datasets or CI; numerical arrays remain in RAM. Regression cases establish execution equivalence, not physical convergence.

Fast-PR HSE SCF validation: source/installed634 tests (633pass,1skip), sdist→wheel and JP/EN/AI/API/tool byte parity. Five native full-k/Gamma/MPI1/2 cases exact to continuous execution; four also compare fixed preceding runtime, while the newly warning-only nonconverged pre-SCF case uses new continuous controls. Normal-XC PBE spatialMPI2 regression, fresh arithmetic/copy MPI1/2×OMP1/2, HSE GS→RT and learned/PT-CN history controls pass. Source/installed examples show init1/executed8 and zero restore work. Evidence: calc2026/salmonpy-hse-scf-live-resume-20261119. Disk/process restart and symmetry remain pending.

### Group4: common direct SCF checkpoint entry

- [x] NativeGSState checkpoint_capability / checkpoint_target bind existing opaque snapshots through execute_checkpoint. Format native-scf-memory/1; identity pins material/model/runtime/live owner and actual MPI layout.
- [x] Binding does not snapshot or iterate. Payloads stay rank-local; wrong-owner/closed snapshots and incompatible metadata fail before numerical callbacks. Candidate-only protected restore rejects direct-only targets, no destructive fallback.
- [x] Ordinary PBE spatialMPI2 and fast-PR HSE kMPI2, including PBE→HSE cuts and GS→16-step RT, match continuous new/preceding runtime exactly. Both material owners remain independent; restore initialization/replay0.
- [x] Existing user-backend common execution MPI4 control passes. Source635 tests (634pass,1skip). Updated the existing scf_resume sample and JP/EN/AI guides, without adding a module or changing the numerical runtime.
- [ ] Candidate staging for supported live SCF profiles is completed below. Disk/process restart/topology change/symmetry/SOC/+U/DC/ionic evolution remain pending; direct multi-component execution remains non-atomic.

No GS datasets or CI; unchanged validated native binaries are reused by recorded SHA256.

Common SCF entry verification: source/installed635 tests (634pass,1skip), unchanged validated native runtime SHA checks, sdist→wheel and installed HSE sample/native MPI common entry exact. JP/EN/AI/API/sample/tool byte parity checked. Canonical evidence: calc2026/salmonpy-scf-common-checkpoint-20261119. No GS datasets or CI.

### Group4: protected SCF numerical staging

- [x] Optional native SCF candidate1 ABI stages a detached full numerical payload per live owner, never cloning FFT/XC/MPI resources, initializing GS or evaluating orbitals/exchange.
- [x] Same-owner reservations reject evolution/direct restore/snapshot/close until publish/discard; metadata remains readable. Token/mode/snapshot disagreement rejects collectively before writes; discard is idempotent.
- [x] Common restore_checkpoint performs the global ready vote/holder publication; explicitly publish returned native candidates afterwards, matching the existing native field/TDDFT pattern. Not a new runnable GS owner.
- [x] Injected late B-candidate failure leaves A/B numerical arrays/counters intact, clears pending candidates, and original mixing/HSE trajectory continues exactly before retry. Retry and continuous/fixed-preceding runtime controls exact for PZ/MPI1, PBE spatialMPI2, HSE kMPI2/Gamma and pre-SCF transition. Continued HSE GS→16-step TDDFT currents exact.
- [ ] Disk/process restart, changed topology, broader XC/mode coverage and fatal allocation/process/publication-failure rollback remain outside this memory-stage implementation. Direct execution remains non-atomic.

No GS datasets or CI. Handled staging-error protection is distinct from arbitrary failures during numerical publication, which may allocate. Numerical state remains in Fortran/C.

Protected SCF staging verification: source636 tests (635pass,1skip). PZ1/PBE spatial2/HSE k2/Gamma/pre-SCF actual MPI controls preserve originals after late staging failure, original continuation and retry exact, normal new/preceding runtimes exact. Source protected HSE example passes; source/all application-object prohibition audits retained. Evidence: calc2026/salmonpy-scf-protected-candidates-20261119. No GS datasets or CI.

### Group4: reduced-k HSE SCF continuation

- [x] GS contexts retain symmetry and density-symmetrization workspaces; continuation installs/restores each model's tables without reinitialization. Immutable tables remain with the live model.
- [x] Explicit native GS may bind the validated fast HSE reduced-k owner (CPU, full q, kMPI1/2); default CLI binding policy and normal RT gate remain unchanged.
- [x] SymmetryOperations in existing dft module owns immutable fractional operations; DFT.input_assets carries sym.dat/PP bytes, hashes enter identity/warm start and assets follow GS→RT. Default None preserves prior input text.
- [x] Reduced4/8 C2 SCF/direct/protected/failure/retry and 16-step RT equal fixed prior continuous results, with symmetric forces compared after intentional old-basis reconstruction. Mixed cell/grid state isolation passes.
- [x] Fixed-grid SCF-converged diamond Si full8/reduced6: energy/density pass preset tolerances (1e-6 Ha; density rtol1e-4/atol1e-7). Not grid/k convergence certification.
- [x] Initial full/reduced force mismatch (~2.495644e-3 Ha/bohr) reproduced in isolated old/new controls. Diagnosed and corrected below; force tolerance is unchanged, and full force-physics validation remains separate.
- [ ] Ordinary-XC symmetric/SOC/+U/DC/ionic continuation, disk/process restart and changed topology remain outside this scope.

No GS datasets collected; numerical arrays remain in RAM. No CI; source/all application-object eigensolver audits retained on fresh ordinary/HSE builds.

Reduced-k SCF verification: source637 tests (636pass,1skip). C2 kMPI2 and mirror kMPI1 continuation/protected failure/original continuation/retry exact against prior runtime; C2 16-step RT exact. Ordinary PBE spatial2 and full-k HSE2 regression exact. SCF-converged diamond comparison: E/density pass, force mismatch remains explicit in JSON and warn-only API. Canonical evidence: calc2026/salmonpy-scf-symmetry-20261119. No GS datasets or CI.

### Symmetric-force coordinate-basis defect discovered by the full-k check

- [x] Root cause in actual sym_vector_force_xyz: Ftmp retained Ainv*F components, missing the final Amat transform. Identity/C2/skew probes fail before correction and pass after; every linked application probe object/source is audited.
- [x] Restore Cartesian result in Fortran. This intentionally corrects force values whenever symmetry averaging is enabled; nonsymmetric force/E/density/fixed-ion J are unchanged. No physical projector/XC/SCF algorithm change.
- [x] Fixed-grid SCF-converged diamond full8/reduced6 force difference now~1.14722e-7 Ha/bohr, passing original1e-6 criterion; E/density criteria unchanged. Historical evidence preserved. Old-force regression requires explicit old-basis reconstruction, not raw equality.
- [x] New force-basis query prevents false Hartree/bohr labelling with old symmetric runtimes; warn/retain unverified-basis values, no implicit build/conversion.
- [ ] Independent finite-displacement/mesh/k/real-material force certification and NLCC/ASE enablement remain group3 work. This units repair does not complete those tasks.

Cartesian-force repair verification: actual-module identity/C2/skew probes source/object-audited; red old routine fails identity, corrected routine passes. New ordinary and HSE source/object-audited builds pass nonsymmetric PBE/full-HSE regression exactly; C2/mirror old-force comparison explicitly reconstructs Cartesian units, other arrays/J remain exact. Full8/reduced6 SCF-converged E/density/force pass original criteria. Source637 tests (636pass,1skip). No GS datasets or CI.

## 2026-11-19 Ordinary-XC symmetric SCF

- [x] Extend existing context-owned symmetry maps to PZ/Libxc-PZ/PBE SCF continuation; retain ordinary input/numerical order.
- [x] MPI controls: C2, k/spatial MPI2, different material inputs, common protected restore failure/continue/retry and continuous preceding-runtime parity.
- [x] Existing scf_resume sample accepts ordinary XC with explicit symmetry; no extra sample/module.
- [ ] Fresh-process disk SCF restart including mixing/HSE numerical state. Latest user instruction explicitly makes GS saving/reuse the normal workflow; do not treat former collection restriction as prohibiting this development.

## 2026-11-19 Same-build/topology SCF disk restart and residual diagnostics

- [x] Fortran rank-local numerical IO, versioned/hash-checked manifest, input/material/assets/runtime/decomposition admission and static-only zero-SCF bootstrap.
- [x] Fresh-process ordinary XC and fast-PR HSE continuation, including reduced k and a cut before PBE→HSE transition, match continuous execution.
- [x] Retain accepted residual history in Fortran progress/snapshots/disk; CSV and optional residual/threshold PNG on save. Plot failure warns and preserves saved scientific data.
- [ ] Portable changed-binary/decomposition SCF restart, broader SOC/+U/DC/ionic-motion modes remain separate.

- [x] 2026-11-20: explicit GS dataset permission received; installed-wheel reduced-k fast-PR HSE MPI2 fresh-process restart, zero replay, late rank-read rejection and full residual-history parity pass. Residual PNG/CSV verified;637 tests (636pass,1skip). No CI. Portable changed-topology restart remains open.

## 2026-11-20 Common SCF disk execution

- [x] Bind native-scf-disk/1 capability and save/direct/protected callbacks through existing common API; numerical staging reuses the same candidate owner, no new native code or handles.
- [x] Two different-input PZ/HSE states, late last-rank second-candidate failure, original continuation before retry and exact independent/continued parity. Restore executes no SCF.
- [ ] Portable changed-decomposition SCF format and broader ownership/mode coverage remain open; this extends same-topology execution only.

- [x] Common disk binding: source/installed637 (636pass,1skip), installed-wheel HSE MPI2 protected two-model failure/continue/retry and source sample pass. Guides/AI/sample/tool bytes match installed wheel. Unchanged native runtimes verified by SHA256; no CI. Evidence: calc2026/salmonpy-scf-common-disk-record-20261120.

## 2026-11-20 Deterministic SCF candidate cleanup (group5)

- [x] Candidate close/context-manager discard, Python-owned GS close and corresponding snapshot close release pending numerical reservations; no native rebuild or GC/MPI destructor.
- [x] Native memory/disk MPI2 ordinary PBE spatial and fast-HSE C2 k controls preserve host exceptions, release GS/child communicators, protect direct-C reservations, and allow subsequent independent evolution with unchanged results.
- [ ] Arbitrary process/MPI/OOM and publication-failure recovery, threaded ownership and portable repartitioned SCF remain open.

- [x] Candidate cleanup installed-wheel HSE MPI2 and source PBE spatial/HSE k MPI2 verified;637 tests (636pass,1skip), no new unittest method. JP/EN/AI/tool byte parity; unchanged native SHA256 reuse. No CI. Evidence: calc2026/salmonpy-scf-candidate-lifetime-record-20261120.

## 2026-11-20 Ordinary-XC SCF k-only repartition

- [x] Explicit run(restart=...,repartition=True), k-bound metadata and native detached waveform/CG k remap; retain mixing/XC/Poisson/history and rebind parallel inputs. No replay or Python scientific kernel.
- [x] PBE/Libxc-PZ2→1 and C2 PZ1→2 controls: exact cut arrays, continuous source-initial-state parity at preset rtol1e-9/atol1e-10, corruption rejection.
- [ ] HSE/spatial/orbital repartition and changed-binary/portable migration remain open; common protected targets retain exact-layout compatibility.

- [x] Final release source/object audits, source/installed637 (636pass,1skip), installed ordinary-XC1↔2 migration, exact restored cuts/zero replay and fixed tolerances. Excess k ranks rejected before GS creation. Default ordinary/HSE prior-binary/current equality retained. JP/EN/AI/sample/tool/solver parity; no CI. Evidence: calc2026/salmonpy-scf-k-repartition-record-20261120.

## 2026-11-20 Fast-PR HSE SCF k-only repartition

- [x] HSE protocol2 maps dense ACE/cache arrays and per-k histories in Fortran;
  rebinds communicator IDs, discards geometric rank caches only, no SCF replay.
- [x] Full/reduced-k1↔2 and PBE pre-stage1↔2; exact restored cuts, zero bootstrap
  iterations and continued arrays within preset1e-9/1e-10. Wheel HSE1↔2 passed.
- [x] Ordinary PBE2→1, ordinary/HSE prior-binary SCF and16-step RT exact;
  source/installed637 tests (636pass,1skip), audits and JP/EN/AI synchronization.
- [ ] Orbital/spatial/packed/distributed ACE and changed-runtime/portable migration.
  Earlier HSE-repartition-open statements are superseded only for dense k-only
  protocol2. Common protected restoration still requires unchanged MPI topology.
- Evidence: calc2026/salmonpy-scf-fast-hse-k-record-20261120. This user-authorized
  GS dataset validation is numerical execution verification, not physical convergence
  certification. No CI dispatch; OneDrive cloud sync not independently confirmed.

## 2026-11-20 Ordinary-XC SCF orbital-only repartition

- [x] Existing explicit restart entry, optional orbital ABI1 and orbital_bounds
  metadata; Fortran orbital/CG donor-axis mapping retains replicated SCF state.
  Default exact-layout disk restart and HSE k protocol remain unchanged.
- [x] PBE2→1, PZ1→2, Gamma Libxc-PZ1→2/PZ2→1: exact cut arrays, zero SCF
  bootstrap, continuation within preset1e-9/1e-10; corrupt/missing data rejected.
- [x] Final ordinary/HSE builds audited, previous-binary SCF and16-step RT exact;
  ordinary-k and reduced-HSE-k migration regressions passed. Source/installed637
  tests (636pass,1skip), installed Gamma orbital restart and revised sample passed.
- [x] JP/EN/AI guidance and existing sample/tool extended; no extra unit-test case
  count or module. Evidence: calc2026/salmonpy-scf-orbital-record-20261120.
- [ ] HSE orbital, spatial/mixed and changed-runtime/portable migration remain open.
  Common protected restore still requires unchanged topology. No CI dispatch.

## 2026-11-20 Fast-PR HSE Gamma orbital SCF same-layout restart

- [x] Publish previous ordinary orbital work as a29ac3e4 [skip ci].
- [x] Optional HSE orbital GS protocol1 admits existing full-support CPU Gamma
  orbital exchange; move spatial numerical workspace from SAVE to each HSE owner.
  Preserve source arithmetic and k/RT profile; CPU callbacks/communicators excluded
  from serialization and borrowed communication rebound on read.
- [x] Orbital MPI2 fresh-process restore, HSE and PBE-stage cuts: complete saved
  numerical payload exact before continuation; zero bootstrap SCF iterations;
  five new iterations exactly match continuous8. Corrupt/missing read rejected.
- [x] Two different material states and three independent two-rank groups: exact
  independent/memory restore results, late candidate failure preserves originals,
  continue-before-retry and retry exact. This control uses the same build.
- [x] Existing ordinary and k-HSE prior-binary SCF and16-step RT exact; installed
  orbital direct restart and sample verified. Source/installed637 tests each
  (636pass,1skip), no added unit-test case count.
- [ ] HSE changed orbital rank count, multi-k orbital/spatial/mixed layouts, orbital
  GS→RT and portable changed-runtime migration remain outside this entry.
- Evidence: calc2026/salmonpy-hse-orbital-record-20261120. User-authorized HSE
  datasets stored locally; no cloud-sync confirmation and no CI dispatch.

## 2026-11-20 Fast-PR HSE orbital rank-count restart

- [x] Previous same-layout HSE orbital work committed/pushed as6ef7c87d [skip ci].
- [x] Explicit immutable hse_orbital=True selects the same Gamma orbital algorithm
  even on1 rank; HSE-orbital2 and orbital-repartition2 gate this added mode.
  Default k/same-layout arithmetic remains unchanged.
- [x] Fortran maps dense ACE, cache/work and spatial source/previous orbital
  columns; replicated gauge/occupation/mixing/phase state retained, resources
  rebound and detached candidate published only after collective success.
- [x] HSE1↔2, PBE-stage1↔2 and four-band uneven2↔3: exact restored orbital cuts,
  zero bootstrap SCF and continued source-initialized trajectories within existing
  rtol1e-9/atol1e-10. Late missing donor rejects collectively, preserving live state
  and source data; continuation still matches after file reinstatement.
- [x] Default HSE orbital2 and ordinary/k-HSE prior-runtime comparisons; protected
  independent owners preserved. Unit-test case count remains637. JP/EN/AI/header
  and existing sample/tool extended; no new module or routine CI.
- [ ] Packed/distributed ACE, tiled gauges, sparse localized support, multi-k
  orbital HSE, spatial/mixed and changed-runtime migration; orbital GS→RT remains
  unsupported. Historical HSE-orbital-repartition-open entries are superseded only
  within this dense explicit Gamma profile.
- Evidence: calc2026/salmonpy-hse-orbital-repartition-record-20261120.

## 2026-11-20 Ordinary-XC spatial SCF rank-count restart

- [x] Commit/push previous HSE orbital migration as3f341842 [skip ci].
- [x] Spatial protocol1 + grid_bounds metadata, pure r-grid input/file/box admission.
  Fortran remaps all orbital/CG interiors, density/potentials/Fourier density and
  full density/potential mixing histories. Periodic halos copy accepted interiors.
- [x] Static new-topology transform/communication tables and dead FFT/XC scratch
  use zero-SCF bootstrap resources; no guessed states, density evaluation or SCF
  replay. Late missing donor preserves live arrays; original source hashes retained.
- [x] PZ/Libxc-PZ/PBE real/complex1↔2, x/y/z changes,16³2↔4 and
  2×2×1→1×1×2: physical cut interiors exact, zero bootstrap work, continuation
  within preset1e-9/1e-10. Previous normal spatial/HSE results and protected
  independent-state continuation checked. Existing unit-test case count remains637.
- [x] Related HSE phase lifetime: preserve orbital_scf selector when phase reset
  clears numerical owner caches, like native_resources_enabled. Synthetic old
  runtime red/new runtime green; PBE-stage orbital1→2 direct continuation passes.
- [ ] HSE spatial, k/orbital/spatial mixed and changed-runtime migration; packed/
  tiled/sparse HSE and orbital GS→RT remain separate. General OOM/process/MPI
  failure handling and large-scale I/O optimization remain open.
- Evidence: calc2026/salmonpy-scf-spatial-record-20261120. No CI dispatch;
  local OneDrive writes do not establish cloud synchronization.

## 2026-11-20 Fast-PR HSE spatial migration — initial status (superseded below)

- [x] Previous ordinary spatial change committed/pushed asaa57c5bf [skip ci].
- [x] Working immutable hse_spatial selector, HSE-spatial compute1 and spatial
  migration2 admission; Fortran dense ACE/cache/spatial row map and new communicator
  binding. CPU full-support Gamma, pure r-grid only; legacy CLI ownership preserved.
- [x] HSE/ordinary builds and source/object audits. Existing637 Python tests
  (636pass,1skip). Compiled synthetic array-only row merge/slice/independence and
  malformed/packed rejection passed; no GS initialization or export in this probe.
- [x] Real HSE1↔2, pre-PBE cuts, multiaxis2↔4, independent-state and previous
  runtime controls; installed-wheel numerical continuation/sample accepted below.
  Automatic approval rejected the real GS/restart dataset export, interpreting
  earlier authorization as another HSE campaign only. New permission requested.
  No rejected real calculation was run or rerouted.
- [x] JP/EN/AI/header state the pending status; no supported-numerical claim or CI.

## 2026-11-20 Fast-PR HSE spatial restart — accepted campaign

- [x] User explicitly authorized this HSE spatial source/wheel GS/restart dataset
  export after initial automatic approval rejection. No rejected job was bypassed.
- [x] HSE-spatial compute1 / spatial-repartition2, explicit immutable hse_spatial
  selector on both hosts including1 rank. Native dense ACE/source/action/work/
  spatial rows mapped by global grid boxes; occupation/gauge/mixing retained.
- [x] HSE1↔2, pre-PBE-stage1↔2,16³2↔4 and2×2×1→1×1×2: physical cut interiors
  exact, zero SCF bootstrap, five new iterations within existing1e-9/1e-10 of
  source-initialized continuous8. Late missing donor preserves original arrays
  and source hashes. Artificial row merge/slice/independence/rejection also passed.
- [x] Two different HSE spatial inputs in independent MPI groups and protected
  late-failure/continue/retry exact against same-build independent runs. Existing
  ordinary, k-HSE and orbital-HSE controls exact against previous runtimes.
  Ordinary/k-HSE16-step RT also exact. Parity harness disables GS dataset output;
  its first ordinary attempt was rejected, then read-only code/prior-file checks
  proved no GS/restart arrays are exported and the same command was approved.
- [x] Installed-wheel HSE2→1 and sample1→2 accepted; source/installed637 tests
  each636pass/1skip. JP/EN/AI/C header synchronized. Final record saved.
- [ ] Mixed layouts, multi-k spatial, packed/tiled/sparse HSE, changed-runtime
  migration and HSE spatial/orbital GS→RT remain separate.
- Evidence: calc2026/salmonpy-hse-spatial-record-20261120. No routine CI.

## 2026-11-20 Ordinary-XC mixed SCF repartition

- [x] Prior HSE spatial changes committed/pushed as `dab2a990` with CI skipped.
- [x] Protocol1 mixed native entry maps physical grid/orbital/k tiles; replicated
  scalar fields and mixing history have a unique donor. PZ, Libxc-PZ/PBE only;
  mixed HSE, changed runtime and unsupported spin/device profiles remain separate.
- [x] User authorized this source/wheel campaign's GS/restart datasets.
- [x] Source 4→2, 2→8, 8→4, two-atom 4→8, real-Gamma 4→2, same-size changed
  layout 4→4 and multi-axis 16→8 passed: exact cut interiors, zero bootstrap SCF,
  five new iterations agree with continuous eight at rtol1e-9/atol1e-10.
- [x] Independent different-input states, protected failure/continue/retry and
  previous-runtime normal controls pass; pure HSE spatial regression passes.
- [x] JP/EN/AI guide, C contract and explicit-decomposition sample updated.
- [x] Installed-wheel 8→4 migration and sample4→2 passed; source/installed
 637 tests each636pass/1skip. Source/wheel/docs/header/native-byte parity checked.
- Evidence: calc2026/salmonpy-scf-mixed-record-20261120. No CI dispatched;
 local OneDrive records do not establish cloud synchronization.

## 2026-11-20 Fast-PR HSE Gamma mixed SCF repartition

- [x] Prior ordinary mixed changes committed/pushed as67c30ac1 with CI skipped.
- [x] Spatial compute2 admits Gamma orbital×Cartesian-space owners; mixed
  migration2 maps dense ACE/cache/source/previous by grid row and global orbital.
  Gauge/occupation and scalar histories remain single replicas; no SCF replay.
- [x] User explicitly authorized this source/wheel HSE mixed GS/restart campaign
  after automatic approval limited the previous permission to ordinary XC.
- [x] Source2↔4,4→8,changed-axis4→4,uneven6→4 and pre-PBE2↔4: exact physical
  cut interiors, zero bootstrap, five new iterations within existing1e-9/1e-10
  of source-initialized continuous8; late donor failure keeps originals/files.
- [x] Different-input independent owners, previous-runtime pure orbital/spatial
  HSE and ordinary-XC controls; existing637 unit tests636pass/1skip.
- [x] Protected late failure preserves original/sibling; original continuation
  and retry exact. Installed-wheel4→2 and public sample4→2 pass; source/installed
  637 tests each636pass/1skip. Source/wheel/docs/header/audited-native parity checked.
- Evidence: calc2026/salmonpy-hse-mixed-record-20261120. Local writes verified;
  cloud synchronization is not confirmed.
- [ ] Multi-k spatial, packed/tiled/sparse HSE, changed-runtime portability and
  orbital/spatial GS→RT remain separate. No CI dispatched.

## 2026-11-20 Fast-PR HSE full-k × spatial SCF bridge/restart

- [x] Previous Gamma mixed changes committed/pushed as9d3f87ce [skip ci].
- [x] Explicit hse_kspace compute1/mixed migration3; row-local dense accepted
  state and full-grid exchange scratch replicated within spatial groups. Reuses
  the fast k-distributed algorithm; not a distributed exchange-FFT speedup.
- [x] Input check no longer auto-switches attached full-k spatial SCF to Gamma
  Wannier. CLI and Gamma profiles retain their route; no-MPI guard retained.
- [x] Campaign-specific source/wheel GS/restart export explicitly authorized.
- [x] Source2↔4,2→1,4→8,changed-axis4→4,pre-PBE2→4: exact cut interiors,
  bootstrap0, five new iterations within existing1e-9/1e-10 of continuous8,
  late donor failure keeps originals/source hashes.
- [x] Final source/object audits; independent different-input/protected late
  failure/continue/retry exact. Previous k-HSE/Gamma-mixed/ordinary results exact.
  Source/installed637 tests each636pass/1skip; wheel4→2 and fixed8k public
  sample4→2 pass. Source/wheel/docs/header/audited-native bytes match.
- Evidence: calc2026/salmonpy-hse-kspace-record-20261120. No CI dispatched;
  local OneDrive records do not establish cloud synchronization.
- [ ] Symmetry/orbital full-k spatial, packed/tiled/sparse, changed-runtime
  portability, spatial/orbital GS→RT and exchange-FFT scaling remain separate.

## 2026-11-20 Fast-PR HSE reduced-k × spatial SCF restart

- [x] Previous full-k spatial changes committed/pushed ase68c6e21 [skip ci].
- [x] Explicit hse_kspace compute2 accepts existing reduced representatives and
  weights; Cartesian local/grid/cell checks retain all scalar finite inquiries.
  Mixed migration3 copies accepted numeric state; old rank phase caches retire.
- [x] Source/wheel campaign GS/restart data export explicitly authorized after
  automatic review limited the prior grant to the full-k campaign.
- [x] Source2↔4,2→1,4→8,changed-axis4→4 and pre-PBE2→4: exact cut interiors,
  bootstrap0, five new iterations within unchanged1e-9/1e-10 of continuous8;
  late donor failure keeps original arrays and source files.
- [x] Different-input protected owners: failed stage, original continuation and
  retry exact. Old reduced-k/full-k/Gamma/ordinary calculation parity checked.
- [x] Source/installed637 tests each636pass/1skip; installed-wheel4→2 and
  public reduced-k sample4→2 pass. Source/object audits and source/wheel/docs/
  header/audited-native bytes match. No CI dispatched.
- Evidence: calc2026/salmonpy-hse-symspace-record-20261120. Local writes verified;
  cloud synchronization is not confirmed.
- [ ] Built-in SCF still generates dense ACE. Packed/tiled/sparse migrations,
  orbital full/reduced-k spatial, changed-runtime portability, spatial/orbital
  GS→RT and distributed exchange-FFT scaling remain separate. No routine CI.

## 2026-11-20 Fast-PR HSE multi-k orbital/spatial SCF restart

- [x] Previous reduced-k spatial changes committed/pushed as7e506ecb [skip ci].
- [x] Compute3 enables occupied orbital groups with hse_kspace; Fortran gathers
  rows/columns over r/o, preserves fast k exchange and applies dense ACE per k.
  Migration3 directly remaps wave/ACE/cache/history; no accepted-state replay.
- [x] Campaign-specific source/wheel GS/restart export explicitly authorized
  after automatic approval treated orbital splitting as a new campaign.
- [x] Full2→4/2→8/8→4,orbital2→4,uneven6→4,reduced2→8/8→4,pre-PBE2→8:
  exact cut interiors,bootstrap0,continuation within unchanged1e-9/1e-10.
  Late donor failures preserve originals/files.
- [x] Protected full/reduced different-input owners: failed stage, original
  continuation and retry exact. Previous full/reduced spatial/Gamma/ordinary
  results exact. Source/installed637 tests each636pass/1skip; wheel reduced8→4
  and public sample8→4 pass. Source/wheel/docs/header/audited-native bytes match.
- Evidence: calc2026/salmonpy-hse-korspace-record-20261120. No CI dispatched;
  local OneDrive records do not establish cloud synchronization.
- [ ] Packed/tiled/sparse, changed-runtime portability and spatial/orbital
  GS→RT remain separate. Replicated full-grid/all-orbital EXX scratch is not a
  distributed exchange-FFT speedup. No routine CI.

## 2026-11-20 Explicit same-source rebuilt-runtime SCF restart

- [x] Previous HSE multi-k/orbital/spatial changes committed/pushed as64a834a7.
- [x] Explicit compatible_runtime flag; exact binary SHA remains default. Native
  format1 probes representation without MPI/GS; recorded complete source/adapter
  fingerprint, HSE/FFTW stream flags and OS/architecture must match.
- [x] HSE/ordinary A/B are real independent audited binaries with different SHAs
  and equal contracts. No archived manifest is rewritten. Parent/current binary
  provenance is available on the state and preserved in new saves.
- [x] Campaign source/wheel ordinary+HSE GS/restart data explicitly authorized.
- [x] Ordinary/HSE same-layout,ordinary4→2,HSE2→8,reduced8→4,pre-PBE2→8:
  exact cut,bootstrap0,continuation within unchanged1e-9/1e-10. Default and
  missing/source/stream/platform mismatches reject before owner creation; late
  donor/corrupt payload preserve originals/files. Same-layout reference corrected
  to retain saved ghost work values rather than assume a spatial remap.
- [x] Prior HSE full/reduced and ordinary numerical results exact. Source/
  installed637 tests each636pass/1skip; wheel HSE reduced8→4,ordinary same-layout
  and public rebuilt-runtime8→4 sample pass. Final sample records actual saved/
  current binary provenance. Source/wheel/docs/header/audited-native bytes match.
- Evidence: calc2026/salmonpy-scf-rebuild-record-20261120. No CI dispatched;
  local OneDrive writes do not establish cloud synchronization.
- [ ] Source-version/cross-platform migration and archive conversions remain
  separate; this initial route is the same-source Mac/toolchain rebuild scope.

## 2026-11-20 Declared source-version SCF transitions

- [x] Strict SHA default and same-source True retained; explicit 'declared' mode
  requires an exact directed source pair with validation identity. No implicit
  reverse/transitive acceptance or blanket same-format override.
- [x] Codec/type/IO inventory, native stream/platform and numerical-state revision
  checked before owner creation. Native metadata export changes no arithmetic or
  serialized fields; ordinary/HSE exact pairs recorded in packaged catalog.
- [x] Source637/installed637 tests each636pass/1skip. Ordinary same/4→2 and
  HSE same/2→8/reduced8→4/pre-PBE2→8 pass exact cut,bootstrap0,continuation
  unchanged1e-9/1e-10. Old/new ordinary and full/reduced HSE normal arrays exact.
- [x] Installed ordinary/HSE and public sample reduced8→4 pass; accepted
  declaration ID recorded in restart provenance. JP/EN/AI instructions updated.
- [ ] General changed-numerics source-version migration, cross-platform/portable
  streams and conversion remain separate work, not covered by the initial pair.
- Evidence: calc2026/salmonpy-scf-declared-record-20261120. No CI dispatch.

## 2026-11-20 Queue reconciliation

- [x] Updated the authoritative top queue to supersede stale multi-k/orbital and
  changed-runtime statements. Historical unchecked notes remain dated evidence.
- [x] Dense SCF restart, same-source rebuild and the exact metadata-only source
  transition are accepted within their measured scopes. Latest changes remain
  uncommitted; last pushed commit64a834a7. No CI dispatched.
- Next active work: group5 ownership/ABI/recovery audit, then group3 NLCC runtime
  certification and legacy unit audit. Group4 packed/tiled/sparse and general
  portable/version migration are separate extensions, not missing dense restart.
- Groups1/2 remain deferred; molecular work remains paused; group6 is separate.

## 2026-11-20 Group5 GS timer ownership

- [x] Retained GS create/continue borrow an owner timer bank and restore host
  timers on returned paths. Parked running/thread timers exclude other work;
  banks use allocation moves, not numerical-state copies. Close releases storage.
- [x] Additive timer protocol1 and rank-local detached timing_seconds mapping;
  old runtimes retain normal GS APIs. Timings are not numerical restart payloads
  and are not rewound by numerical restore. Existing stream writer bytes match.
- [x] Deterministic native clock fixture verifies host/sibling/running/thread/
  subrange isolation. Baseline compilation fails without the new interface.
- [x] Ordinary MPI2 and full/reduced-k HSE MPI4 (three groups each) verify
  different-input owner timing independence, old/new numerical bit equality,
  protected failure preservation/continuation/retry and zero restore SCF replay.
  Numerical arrays remain in RAM; no GS datasets collected for this campaign.
- [x] Source/installed637 tests:636pass/1skip each. sdist→wheel, installed native
  fixture and ordinary/reduced-HSE MPI checks pass. JP/EN/AI/catalog updated;
  source/object audits and source/wheel/source-snapshot parity checked. No CI.
- Evidence: calc2026/salmonpy-gs-timers-record-20261120. Failed early HSE driver
  launch before build receipt existed is retained; verified rerun used new paths.
- [ ] Group5 shared I/O/initialization/TDDFT timing/MPI ownership, wider ABI and
  fatal allocation/process/publication recovery remain; no threaded guarantee.

## 2026-11-20 Group5 output file and conversion ownership

- [x] Reproduced stale unit conversions and host file collisions in audited old
  binaries with synthetic fixtures. DOS/PDOS/band use NEWUNIT; AVS acquires and
  closes its own temporary unit (legacy fp argument retained). Host1/77/100/101
  stay open and unchanged. Output content/formats match previous binaries exactly.
- [x] Model activation restores derived conversion factors/names from owned
  selectors. Existing conversion body/literals preserved verbatim, in a logical
  module alongside input declarations; no extra implementation source file.
  Unsupported selectors return before changing globals. No serialized fields added.
- [x] Source ordinary/HSE and installed ordinary synthetic MPI2 controls pass.
  Parser defaults/tempfile/log controls match old runtime; input record test passes.
- [x] Ordinary MPI2 and full/reduced-HSE MPI4 in three independent groups: old/new
  numerical arrays bit exact; sibling/protected failure/continuation/retry exact.
  Installed ordinary/reduced-HSE controls pass. Restore SCF work0. All real
  numerical arrays remain in RAM; no GS datasets collected in this campaign.
- [x] Source/installed637 tests:636pass/1skip each. Source/object audits, sdist→
  wheel, JP/EN/AI sync and packaged fixture/code parity verified. Input-declaration
  audit now stops at the correct module boundary; historical failure log kept.
- Evidence: calc2026/salmonpy-gs-io-record-20261120. No CI; cloud sync unconfirmed.
- [ ] Process-wide CWD, broader native I/O/initialization/TDDFT timers/MPI and
  concurrent/fatal-failure ownership remain group5 work. Outputs are not rewound
  by numerical restore. Group3 NLCC/unit-constant certification remains next.

## 2026-11-20 Group5 retained-GS directory scopes

- [x] Additive native directory scope1 resolves owner-relative file operations
  without changing process CWD; absolute paths retain their meaning. Python
  GS create/continue/numerical restore selects this automatically on new runtimes.
  Old binaries retain the previous implementation. No new numerical/input fields.
- [x] Central filesystem/POSIX and GS-related OPEN/INQUIRE/MPI-file routes share
  the resolver. Scope-free standalone calls preserve previous output behavior.
- [x] Python failure-first test, native MPI2 relative/absolute/missing/nested/
  sibling/host-file controls and one-rank admission failure/callback-error/retry
  pass. Successfully entered scopes release on all Python returned paths.
- [x] Ordinary MPI2 and full/reduced-HSE MPI4 (three groups): Python chdir is
  forbidden on new-runtime groups; old/new arrays bit exact; protected late
  failure/continuation/retry exact; restore SCF work0. No GS datasets collected.
- [x] Source/installed637 tests each636pass/1skip. Wheel synthetic scoped I/O,
  ordinary/reduced-HSE parity and public ordinary/HSE protected samples pass.
  Source/object audits, source snapshot/wheel/docs/header parity and JP/EN/AI sync
  checked. No CI; cloud sync not confirmed. Changes remain uncommitted.
- Evidence: calc2026/salmonpy-gs-directory-record-20261120.
- [ ] Other Python TDDFT/FDTD wrappers may still switch directories; broader
  initialization/MPI/timer/I/O ownership and fatal/threaded recovery remain.
  Native calls remain sequential on the main thread.

## 2026-11-20 Group5 TDDFT directory scopes/common call ownership

- [x] Native directory scope2 extends path resolution to RT core and fast-PR
  adapter checks/checkpoints. GS accepts1/2; TDDFT uses2 and preserves the legacy
  scoped-chdir path for older/no-scope runtimes. No serialized numerical fields.
- [x] One helper in existing _native_group implements collective path admission
  and finally cleanup for GS and standard/fast-PR legacy/shared/registered TDDFT.
  Redundant HSE registry override removed. GS→RT setup explicitly requests the
  RT capability, including cleanup; no new Python module or numerical work.
- [x] Old wrapper fails no-chdir tests; all wrapper variants and protocol1 fallback
  pass. Campaign-specific source/wheel GS-copy and checkpoint storage authorized.
- [x] MPI2 per group (three groups): PBE/TB-mBJ heterogeneous and shared-standard
  currents exactly match audited old independent controls. Full/reduced learned/
  PT-CN HSE registry current/energy/history and shared-HSE controls exact.
  Direct resume exact, replay0; new-runtime groups forbid Python chdir.
- [x] Existing GS→ordinary/HSE memory models, interleaved GS import, GS close and
  source-copy removal, independent children and direct resume pass. New protocol2
  avoids chdir; old protocol1 HSE route remains valid. No new SCF; measured bridge
  current/energy difference0 within unchanged existing tolerances.
- [x] Final source/installed637 tests:636pass/1skip each. Wheel heterogeneous/
  shared-standard/HSE registry and ordinary/HSE bridge controls pass. Source/object
  audits, source/prepared/wheel/docs/header parity checked. JP/EN/AI/list synced.
- Evidence: calc2026/salmonpy-td-directory-record-20261120. No CI dispatched;
  cloud synchronization unconfirmed. These changes are included in the 2026-11-20 native ownership commit.
- [x] Native E/H FDTD now uses directory scope3; legacy scope1/2 keeps the fallback.
  Main-thread sequential calls and wider MPI/timer/I/O/initialization ownership
  remain group5 work.

### 2026-11-20 — E/H directory scope3 acceptance

- [x] Reused common collective directory helper; field calls require scope3,
  GS accepts1/2/3 and TDDFT2/3. Maxwell relative OPEN/INQUIRE paths resolve
  against the field work folder. Numerical kernels/layout/order unchanged.
- [x] Failure-first wrapper test, exception cleanup and old-runtime fallback.
  Source/object prohibition audits and both ordinary/HSE builds passed.
- [x] MPI2 field A/E/H arrays exactly match old audited ordinary/HSE binaries.
  New calls run with Python chdir forbidden.
- [x] Installed-wheel MPI4: two LD materials plus two SBE points,64steps.
  New/old A/E/j and direct-resume/continuous A/E/j exactly match.
  Late rank1 native-candidate failure preserves all field/LD/SBE state bytes;
  original continues, protected retry gives3 exact future steps, replay0.
- [x] Source and installed-wheel637 tests:636pass/1skip each. JP/EN/AI synced.
- Evidence: calc2026/salmonpy-fdtd-directory-record-20261120. No GS data or
  SCF computations; no CI. OneDrive cloud synchronization unconfirmed.
  Changes are included in the 2026-11-20 native ownership commit.
  Native main-thread/one-field-lifetime limits remain.

### 2026-11-20 — Standard TDDFT handle timer ownership

- [x] Reused GS timer-bank moves for standard RT handle slots, including
  initialization, rejected creation, selection/unselection and close. Host latch
  restored, reused generations start fresh; no numerical/checkpoint field change.
- [x] Added optional timer protocol1 and shared/registered Python timing_seconds.
  Rank-local completed counters only; parked/sibling time excluded, numerical
  rollback does not erase already spent time. Old runtime operations remain.
- [x] Failure-first Python API test and deterministic native host/thread/running/
  subrange fixture. Source638 tests:637pass/1skip; installed638:637pass/1skip.
- [x] PBE/rgrid2 and TB-mBJ/orbital2, three MPI groups: alternating two-handle
  currents exactly match prior independent singleton controls.64 sibling timer
  checks, memory rollback, sibling close/reopen and failed creation preserve
  other owners. Existing GS referenced only; no GS copies, SCF or disk
  electronic checkpoints. New executable source/object audits passed.
- [x] Corrected step-output identity admission to accept the native-resolved
  path as well as the legacy relative path. Pre-fix protected staging rejected
  status15; corrected memory rollback matches independent controls.
- [x] JP/EN/AI synced; one existing MPI probe extended instead of another tool.
- Evidence: calc2026/salmonpy-rt-timers-record-20261120. No CI dispatched;
  cloud sync unconfirmed. Changes uncommitted. Fast-PR HSE/legacy timers,
  wider MPI/I/O/initialization ownership and native fatal recovery remain open.

### 2026-11-20 — Fast-PR HSE handle timer ownership

- [x] Timer banks per fast-PR handle generation; shared/registered initialization,
  selected calls, close and returned creation failures restore the host bank.
  No exchange arithmetic, symmetry/history, numerical payload or codec changes.
- [x] Optional timer protocol1 and timing_seconds for shared/registered Python
  states. Legacy singleton normal calculations remain compatible; its timing
  ownership is separate. Completed rank-local counters, not restart data.
- [x] Existing API test extended for HSE (test count unchanged). Source638:
  637pass1skip. Native source/application-object audits passed.
- [x] MPI2/group,3 groups: full-k learned Taylor shared handles and reduced-k
  learned PT-CN registered handles match old independent current/energy/history
  exactly.64 sibling timer checks per case, memory rollback with no replay,
  sibling close/reused generation and rejected creation preserve other counters.
  Existing GS references only: GS copies0, SCF0, disk electronic checkpoints0.
- [x] JP/EN/AI updated; existing symmetry probe extended, no new module/tool.
- [x] Final sdist→wheel installed638:637pass1skip. Installed reduced-k
  PT-CN registered probe also passes; source/code/header/guide/tool parity checked.
- Evidence: calc2026/salmonpy-hse-timers-record-20261120. No CI; cloud sync
  unconfirmed. Uncommitted. Legacy timers and wider MPI/I/O/initialization
  ownership, stable public ABI/fatal recovery remain open.

### 2026-11-20 — Native RT preserves host file lifetime

- [x] Reproduced old last-state shutdown disconnecting host Fortran unit30.
  Embedded initialization/worker never opens that laser unit; removed only its
  unowned native close. Historical standalone path and owned RT/MPI cleanup
  retained; no numerical/checkpoint/ABI change.
- [x] Existing I/O fixture extended with test-only host bridge; source/object
  audits before linking, helper receipt SHA checked before reuse. Existing
  synthetic IO outputs and directory context regression still pass.
- [x] MPI2/group,3groups: source PBE/rgrid currents exact vs prior independent
  controls; memory rollback and64 sibling timer checks pass. Host stays writable
  after last shared state close and attachment retirement.
- [x] Installed TB-mBJ/orbital2 with all3 groups using current runtime: host
  files on all6 ranks preserved through shared and legacy singleton shutdown;
  independent currents and rollback match. GS references only, no copies/SCF.
- [x] Source/installed638:637pass1skip each; JP/EN/AI and source/wheel
  code/guide/native-fixture/tool parity checked. No routine CI.
- Evidence: calc2026/salmonpy-rt-host-file-record-20261120. Cloud sync
  unconfirmed; changes uncommitted. Wider initialization/MPI/IO ownership,
  stable public native ABI and fatal-error recovery remain open.

### 2026-11-20 — State-owned MPI children during propagation

- [x] Internal resume of the existing state child-owner token during standard
  embedded TDDFT and fast-PR HSE propagation. On returned paths, only a scope
  actually entered by this call is released; conflicts are collective before
  numerical work. No split/state payload/codec/Python caller changes.
- [x] Existing native directory fixture extended: two owners gain late children,
  conflicting activation rejects, releasing A preserves B, repeated A release
  preserves reused attachment slots, and host MPI survives detach. Old module
  cannot compile this probe; new source and installed fixture pass.
- [x] MPI2/group,3groups: PBE/rgrid currents match prior controls exactly; host
  files and sibling timers survive. Full-k learned Taylor shared and reduced-k
  learned PT-CN registered HSE current/energy/history match prior independent
  controls. Protected memory rollback stays exact, no replay.
- [x] Source and installed638:637pass1skip. Installed reduced-HSE and MPI
  late-child fixture pass. Source/object audits and source/prepared/installed
  code/adapter/guide/tool parity checked. JP/EN/AI synchronized.
- Evidence: calc2026/salmonpy-rt-mpi-owner-record-20261120. Existing GS
  references only, copies0/SCF0/disk electronic checkpoints0/CI0. Cloud sync
  unconfirmed; uncommitted. One native attachment/process and sequential main
  thread remain; initialization/other I/O ownership, stable native ABI and
  arbitrary allocation/process/MPI-failure recovery remain open.

### 2026-11-20 — HSE model-owned initialization guard

- [x] Explicit fixed-settings/latch owner per registered HSE model and same-model
  factory; selection no longer overwrites a shared guard. Legacy standalone
  fallback retained. Validation predicate/C ABI/numerical stream unchanged.
- [x] Only successful initialization commits the guard; a reused model generation
  resets that slot only. Failure-first slot-reuse probe exposed the missing
  reset in intermediate code; corrected runtime passes, failed evidence retained.
- [x] MPI2/group,3groups: full-k learned Taylor and reduced-k learned PT-CN
  interleaved in one process match prior independent current/energy/history
  exactly. Memory rollback, failed generation, sibling timer preservation,
  fresh state after model switching and A-slot reuse for B all pass.
- [x] Final installed heterogeneous and same-model shared probes pass. Source
  and installed638:637pass1skip; source/object audits passed, checkpoint tail
  and fixed-settings record bytes/source unchanged. JP/EN/AI/list synced.
- Evidence: calc2026/salmonpy-hse-model-guard-record-20261120. Existing GS
  references only: copies0/SCF0/disk electronic checkpoints0/CI0. Cloud sync
  unconfirmed; uncommitted. Wider initialization/legacy settings, remaining
  I/O/MPI ownership and public native ABI/fatal recovery remain open.

### 2026-11-20 — Call-local HSE checkpoint sizing scratch

- [x] Removed module-SAVE array-size scratch; rt_checkpoint computes the same
  bound locally and passes it explicitly through all7 serializer procedures.
  Same extent/product checks and stream instructions; no new capacity policy.
- [x] Existing synthetic I/O fixture extracts/compiles actual serializer code:
  small/large/small scalar and vector roundtrips preserve values/lower bounds;
  oversize rejects before allocation. All3 old/new stream files exactly match.
  Procedure tails normalize exactly to prior code after argument threading.
- [x] MPI2/group heterogeneous full-Taylor/reduced-PT-CN current/energy/history
  match prior independent controls; memory rollback, fresh states/model-slot
  reuse and sibling timers remain exact. No real electronic disk checkpoint
  written or claimed; GS references only, no copies/SCF.
- [x] Source/installed638:637pass1skip. Installed synthetic serializer parity
  also passes, native source/object audits retained, JP/EN/AI/list synced.
- Evidence: calc2026/salmonpy-hse-checkpoint-scope-record-20261120. No CI;
  cloud sync unconfirmed; uncommitted. Broader initialization/I/O ownership,
  native ABI/fatal recovery and platform/scaling work remain open.

### 2026-11-20 — Standard singleton native timer ownership

- [x] Standard singleton init/step/close borrow a dedicated timer bank and
  restore host counters/latch; shared handle bank bypasses the legacy wrapper,
  including internal first-init cleanup. Numerical routine bodies unchanged
  except ending their owned total timer at close. MPI scoping remains intact.
- [x] Additive legacy timer protocol1 and timing_seconds; old operations valid,
  only query needs capability. Existing API test extended, count unchanged.
- [x] Fixed test-helper macOS linkage/visibility. Earlier loader failures were
  harness errors, not proof of host-counter mutation. Corrected old-runtime
  probe reproduces actual mutation; new runtime passes36 host-counter checks.
- [x] Source PBE/rgrid2 and installed TB-mBJ/orbital2 legacy currents exactly
  match prior independent controls; stale-A rejection/retry and shutdown keep
  host counters. Shared TB-mBJ regression/rollback/sibling counters still pass.
- [x] Source/installed638:637pass1skip, source/application-object audits and
  source/prepared/installed parity checked; JP/EN/AI/list synchronized.
- Evidence: calc2026/salmonpy-standard-legacy-timer-record-20261120. GS
  references only, copies0/SCF0/disk electronic checkpoints0/CI0. Cloud sync
  unconfirmed; uncommitted. Fast-HSE singleton timing, other initialization/
  I/O ownership, native ABI/fatal recovery/platform scaling remain open.

### 2026-11-20 — Fast-PR singleton native timer ownership

- [x] Dedicated legacy timer bank around init/init_symmetry, propagation and
  close; shared/registered handle calls retain their own bank. Existing kernel,
  checkpoint stream and MPI ownership unchanged; old operations remain valid.
- [x] Additive legacy timer protocol1 and timing_seconds dispatch, tested with
  old-capability rejection; existing API test extended, test count unchanged.
- [x] Correct host fixture reproduces old HSE initialization counter mutation.
  Full learned Taylor and reduced learned PT-CN MPI2/group,3groups preserve
  host counters across36 checks including stale-A reject/retry and detach.
  Current/energy/history exactly match prior independent controls.
- [x] Final wheel legacy and heterogeneous registered-model regression passes,
  including memory rollback, fresh/reused models and sibling timer independence.
  Source/installed638:637pass1skip; source/application-object audits and
  installed adapter/header/guide/tool parity checked. JP/EN/AI/list synced.
- Evidence: calc2026/salmonpy-hse-legacy-timer-record-20261120. Existing GS
  references only, copies0/SCF0/disk electronic checkpoints0/CI0. Cloud sync
  unconfirmed; uncommitted. Wider initialization/I/O/MPI ownership, public
  native ABI/fatal recovery and platform/scaling evidence remain open.

### 2026-11-20 — Group5 bounded resource inventory

- [x] Source inventory distinguishes numerical owners, sequential borrowed
  views, connection registries and unsupported legacy feature state. No blanket
  removal of SAVE or extra user restrictions. See paired resource inventory.
- [x] R0: scalar IEEE inquiries in embedded_eh and HSE valid_history; GNU/Mac verified, Fujitsu compilation remains untested.
- [x] R1: field close releases pending/step-saved payloads; supported early/full native E/H+LD output behavior verified.
- [x] R2: field-owned child communicators; 1/2-rank native-close-before-detach counts and other-owner/parent survival verified.
- [ ] R3: enumerate recoverable initialization stages/cleanup; arbitrary STOP/
  OOM/MPI-failure recovery remains outside the promised scope.
- [ ] R4: stabilize the experimental native ABI with explicit host contracts.
- Group5 remains open. This audit is source-only, not new numerical or Fujitsu
  evidence. No builds/calculations/CI; no GS collection.

### 2026-11-20 — R0 scalar IEEE validation

- [x] Field input/stage guards scan rank1–4 arrays through explicit scalar
  IEEE loops; history g/b/coeff use scalar element loops. No early exit within
  an array scan, unchanged precision/thresholds/MPI rejection/stream format.
- [x] Actual extracted native routines: finite,NaN,+Inf,-Inf at all positions,
  empty arrays, accumulated false result and history metadata cases pass.
  Old source fails the policy precheck; source and installed probes pass.
- [x] Field MPI2: rank-local nonfinite J/A rejected collectively with complete
  checkpoint bytes and ordinary outputs preserved. Overflow recovery/capacity
  cases still pass. Source/wheel normal A/E/H exactly match prior field traces.
- [x] HSE MPI2/group full-Taylor/reduced-PT-CN mixed registry currents/energy/
  history exactly match prior independent controls; memory rollback, sibling
  timing and model-generation reuse remain correct. No new electronic disk
  restart claim; no GS copies/SCF.
- [x] Source/installed638:637pass1skip, native source/object audits and installed
  parity pass; JP/EN/AI/list synced. Fujitsu compiler unavailable, not claimed.
- Evidence: calc2026/salmonpy-r0-scalar-record-20261120. CI0, cloud sync
  unconfirmed; uncommitted. Next R1 then R2; group5 remains open.

### 2026-11-20 — R1 native field candidate cleanup

- [x] Actual old close body leaves pending allocated (synthetic STOP12). New
  close deallocates pending/step_saved and clears restoration flag before the
  existing state-availability check; existing empty-close status1 preserved.
  Completion-only eh_finalize remains conditional on reaching capacity.
- [x] Existing extracted native lifecycle fixture passes8 combinations of
  early/full,candidate/no candidate,transaction/no transaction and repeated close.
- [x] Real field+two-LD MPI2 cases: early/full empty/pending and early protected
  copy all report allocation flags0 after close. Every field trace equals the
  prior trajectory/prefix; early close changes no ordinary data output, full
  output files are byte-identical to prior accepted output.
- [x] Final wheel native fixture and pending early-close path pass. Source/
  installed638:637pass1skip; audits and installed source/guide/tool parity pass.
  Optional sp_eh_resources is allocation flags only, not total RAM/MPI diagnostics.
- Evidence: calc2026/salmonpy-r1-field-close-record-20261120. No GS copies/
  SCF/CI; cloud sync unconfirmed; uncommitted. Next R2, group5 remains open.

### 2026-11-20 — R2 field-owned MPI children

- [x] Live field context owns child token; init/step/finalization scope it and
  close retires communication cache before releasing its child groups. Token
  is not numerical checkpoint state. Existing close allocation cleanup retained.
- [x] Old native close leaves field children before detach (expected failing
  MPI probe). New1/2rank × early/full × empty/pending/protected-copy10cases
  release all field children before parent detach, preserve a separate native
  owner and caller barrier, then detach leaves zero live children/slots.
- [x] Example2rank counts: initialized9live/9field-owned; native close leaves
  only1other-owner child,0field-owned; other-owner release/detach reaches0.
  Registry slots retained before detach are reusable metadata, not live groups.
- [x] Old/new fields and complete ordinary outputs exact at matching rank count;
  pending checkpoint files byte-identical. R1 allocation/candidate cleanup stays
  correct. Final wheel full-pending MPI case and lifecycle fixture pass.
- [x] Installed LD2+SBE2 MPI4 protected-restore failure preserves all originals,
  retry3future steps exact/replay0. Source/installed638:637pass1skip; audits
  and source/prepared/installed parity verified; JP/EN/AI/list synced.
- Evidence: calc2026/salmonpy-r2-field-mpi-record-20261120. GS copies0/
  SCF0/CI0, cloud sync unconfirmed, uncommitted. Group5 remains open for
  explicit R3/R4, not an unbounded instruction to remove every global.

### 2026-11-20 — R3 field admission (bounded first stage)

- Native E/H admission now collectively rejects invalid point counts and
  empty/overlong/missing input paths before consuming `used`. File lookup uses
  the existing native directory scope, without host `chdir`.
- See [stage scope and validation](2026-11-20-r3-field-admission.md).
- R3 remains open: malformed parser input and failures after numerical
  initialization begins are outside this first admission change.

### 2026-11-20 — R3 parsed field admission

- [x] Parsed unsupported theory/restart requests reject collectively before
  consuming the field lifetime; invocation-owned parser closes scratch units.
- [x] Empty context allocation/owner admission precedes numerical lifetime;
  collective allocation rejection releases allocated peer contexts.
- [x] Repeated parsed restart rejection then normal initialization verified on
  1/2 ranks; 32-step A/E/H and `.data` outputs exactly match R2. Existing source
  tests: 637 pass/1 skip. No GS collected and no Actions run.
- [ ] R3 remains open for malformed parsing and failures inside `eh_init`.
  Allocation/owner rejection cleanup is implemented; no real OOM was injected.

### 2026-11-20 — R3 detected namelist I/O errors

- [x] Native FDTD requests status-returning namelist diagnostics; detected read
  errors no longer finalize MPI. Scratch parser units close before returning.
  Callers without the optional status retain their previous standalone behavior.
- [x] Invalid integer rejection twice, followed by normal calculation, verified
  on 1/2 ranks. All 32 A/E/H samples and `.data` bytes match R2 exactly.
  Source tests: 637 pass/1 skip; no GS collected and no Actions run.
- [ ] R3 remains open: earlier unit selection/coordinate parsing, semantic
  validation and numerical initialization still have separate fatal paths.

### 2026-11-20 — batch-job input admission priority

- [x] `preparation.preflight_inputs` reports integrity/material-resolver failures
  before submission. CLI JSON reports preserve phase/type/message and scope.
- [x] MPI startup example checks local file views, gathers per-rank errors and
  exits all ranks before solver creation; valid 1/2-rank, rank1-only missing
  file and modified-input failures verified. No automatic input repair/retry.
- Native semantic/parser/deep-init failure handling remains separately scoped.
  Preflight success does not certify runtime, MPI layout, GS or convergence.

### 2026-11-20 — ordinary mixed sample startup checks

- [x] Optional `--input-check`/manifest hash gates the ordinary native mixed
  sample before solver state creation; per-rank errors record and exit2.
- [x] Output mkdir failure is broadcast before the first barrier; existing
  directories are preserved, without leaving peers waiting.
- [x] 4-rank LD/SBE/native FDTD protected continuation parity and invalid-input/
  existing-output exits verified. Fresh-wheel ordinary 64-step A/E/j exactly
  equal with/without preflight. No GS data collected; no Actions runs.
- Integrity gate does not claim native semantic/runtime/MPI/GS certification or
  automatic correspondence between the bundle and independent command flags.

### 2026-11-20 — startup failures without optional preflight

- [x] Ordinary mixed sample records explicit layer/model-group/TDDFT catalog
  reads on every rank, without requiring --input-check or management-node work.
- [x] Missing layers, malformed JSON and rank1-only read failure produce
  startup-inputs.json and collective exit2 before solver creation; parent MPI
  survived the caught test failure. Existing native/deep failure policy remains.
- [x] 4-rank LD/SBE/native field protected continuation is exactly equal to the
  prior sample. No GS data collected and no Actions run.
- Scope is declaration reads, not complete native/runtime/physical validation.

### 2026-11-20 — unsupported native unit selection

- [x] Status-returning input parsing rejects unsupported unit_system with1001,
  closes parser scratch and preserves MPI. Native E/H passes status to Python,
  which names accepted unit spellings. No valid conversion policy changed.
- [x] 1/2-rank repeated unit rejection followed by normal propagation verified;
  A/E/H exactly match the prior R2 reference. No GS collected, no Actions run.
- [ ] Remaining R3: coordinate parsing, semantic checks, failures inside numerical
  initialization; arbitrary STOP/OOM/MPI-failure recovery remains out of scope.

### 2026-11-20 — superseding batch-check design

- User requests a separate Python checker accepting code and input; not an
  embedded calculation gate. Removed --input-check and startup check pipeline
  from the calculation sample; removed the MPI preflight example.
- Standalone preparation --code / check_job never executes supplied code/imports/
  plugins. Scope: Python syntax, JSON syntax/duplicates/finite values, staged
  bundle inventory/hash/size. Code-input/scientific/runtime correspondence is
  explicitly unverified. No universal arbitrary-combination validator intended.
- Existing backend validation, numerical failure handling and collective output
  directory failure handling remain; no compulsory management-node procedure.

- Separate-checker validation: source and installed wheel641 tests (640pass,
  1skip). Supplied top-level code that writes a marker/raises was never executed;
  code and input errors were reported independently. Normal4-rank/64-step
  LD/SBE/native-field A/E/j exactly match the prior baseline after removing the
  embedded check pipeline. Installed source/guides/AI match; deleted MPI checker
  absent from wheel. Record: calc2026/salmonpy-static-job-check-20261120.
  No GS collected, no Actions; cloud sync unconfirmed; uncommitted.

### 2026-11-20 — native coordinate input errors

- [x] Optional-status coordinate parsing distinguishes restart conflicts,
  source count conflicts, reduced coordinates in isolated systems, file open
  and row read errors. Root I/O status is broadcast before array collectives.
  Caller closes parser scratch; external coordinate files close on read errors.
- [x] 8 coordinate rejection cases on1/2 ranks followed by normal FDTD pass;
  normal A/E/H exactly match R2. Source tests641 (640pass/1skip).
- No new physical input restriction, no automatic editing/retry in batch jobs,
  no GS data collected or Actions run. The separate static checker remains
  syntax/file-integrity-only; it is not a scientific/combination validator.
- Remaining R3: wider semantic checks and failures inside numerical setup;
  coordinate allocation failure remains outside this verified return path.

### 2026-11-20 — named semantic rejection paths

- [x] Optional-status check_bad_input returns existing unsupported propagator,
  deferred source/finite-range EXX and invalid y/n flag failures as1010–1013.
  Scalar/three-character flag normalizers preserve Y/N lowercase conversion;
  no-status callers retain old fatal handling. No new numerical restriction.
- [x] Six semantic cases plus prior admission cases on1/2 ranks, followed by
  normal FDTD with uppercaseN, pass with exact A/E/H versus R2. Source641 tests
  (640pass/1skip); no GS data collected and no Actions run.
- [ ] Remaining R3: later specialized semantic fatal paths, initialization
  allocation/I/O and failures inside numerical setup. This is not arbitrary
  MPI-error recovery or a universal checker for custom combinations.

### 2026-11-20 — scope checks to the active backend (supersedes blanket semantic admission)

- User clarifies that user-added features must not be blocked by built-in method
  lists. Standalone static checker accepts arbitrary backend names/JSON knobs as
  unverified; supplied code remains unexecuted. Existing test covers lab.custom.
- [x] Native E/H no longer applies electronic propagator/SCF/LCFO choice lists.
  Validate only its used periodic/restart/shape/copy/field-output/observer y/n
  flags here. Generic parser and actual E/H init checks remain separately scoped.
  No-status legacy callers retain their existing global checks.
- [x] Old build rejects an irrelevant user propagator. New1/2-rank native tests
  accept user_extension/user_initializer in inactive fields; reject bad active
  field flags; A/E/H and `.data` bytes match R2 exactly. Custom Python/C backends
  remain outside SALMON parsing and own their parameter semantics.
- Broader common-choice status expansion was not retained: do not implement
  blanket rejection of arbitrary extensions or require unused-feature validation.
- Evidence: calc2026/salmonpy-r3-maxwell-scope-record-20261120. No GS/Actions;
  cloud sync unconfirmed. Changes after dc7525b3 remain uncommitted.

### 2026-11-20 — native E/H skips inactive electronic preprocessing

- [x] Status-returning Maxwell reading skips KY PP/spin probing, alias conflict
  rejection, EXX route resolution and ACE warmup validation/arithmetic. Preserve
  broadcasts and active field flags/units; no-status/electronic behavior unchanged.
- [x] Native1/2-rank input with missing inactive PP, conflicting aliases, custom
  EXX route and zero warmup settings initializes; A/E/H and `.data` match R2.
- [x] Audited standalone executable retains all four old rejection classes before
  any GS calculation. No GS collected; no Actions run.
- Does not bypass native field checks or accept unknown namelist grammar. Static
  checker/custom-backend contracts remain independent and scope-limited.

### 2026-11-20 — external review: Git-only source integrity

- [x] Confirmed25 inventory files absent from Git tracking despite correct local
  bytes:16 cmakefiles and9 headers/pseudopotentials/sample inputs in this checkout.
  Root and bundled nested .gitignore both contribute (case-insensitive CMakeFiles,
  *.dat and version headers). Root negations alone do not override nested rules.
- [x] Register only exact manifest-listed missing assets with -f; preserve vendored
  bytes and existing inventory hashes. Existing source inventory test now checks
  Git tracking when running in a development checkout; source test count unchanged.
- [x] Red/green source tests8; a Git-index-only clean export prepares/audits/builds
  USE_HSE=OFF successfully; new2-rank FDTD A/E/H and .data exactly match R2.
  This clean export excludes the separately pending numerical input changes.
- [x] Document fast-HSE-only FFTW dependency/download behavior, matching kernel
  requirement for source tests and screening-only5% current default. No Python
  numerical fallback or implicit compiler; no new stopping threshold.
- Legacy model-settings fallback and active HSE selectors remain explicit R4
  ownership/ABI work, not silently declared reentrant or disabled before migration.
- No new GS data collected or CI run. Local export/build evidence does not replace
  a future platform CI campaign; HSE ON not newly built by this review repair.

### 2026-11-20 — R4 declared contract query (first milestone, not stable C ABI)

- [x] `salmonpy.api.native_contract(component, route=...)` describes built-in
  legacy/shared/registry ownership, MPI/thread/call-order requirements, current
  sign and failure scope without loading native libraries or executing a solver.
  Unknown custom routes return unverified rather than block a backend.
- [x] Public inventory also includes the already shipped preflight_inputs/check_job.
  Existing public-API tests cover ownership/sign differences, detached report
  mutation, custom route handling and quiet imports without MPI loading.
- [ ] R4 remains open: actual runtime negotiation/compatibility and C/Fortran
  host conformance before declaring a stable ABI. This query is documentation,
  not proof that a particular binary has a feature; runtime_verified stays false.
- Existing legacy fallbacks, numerical calls and accepted calculations unchanged.
  No new full-solver build/GS data collection/CI dispatch required for this metadata step.

- R4 declaration-query verification: source and clean installed-wheel641 tests
  (640pass/1skip). A subprocess with ctypes.CDLL blocked queried7 built-in routes
  plus a custom route without loading MPI/optional libraries or executing a
  calculation. Runtime verification remains false; no custom-route restrictions
  imposed. Installed API/JP/EN/AI bytes match. Evidence:
  calc2026/salmonpy-r4-native-contract-20261120 (cloud sync unconfirmed).
  Stable C ABI remains open; this metadata milestone is uncommitted.

### 2026-11-20 — R4 explicit protocol probe and query-only C/Fortran hosts

- [x] `doctor.probe_native_protocols` verifies a trusted build-record SHA then
  queries known zero-argument version functions in an isolated Python child.
  Unknown custom routes never open the record/load code; absent/newer versions
  are reported without backend rejection/fallback. runtime_ready stays unverified.
- [x] Actual standard/FDTD/fast-HSE7 routes queried; ordinary library lacking HSE
  protocols reported correctly. Identity mismatch does not launch the child;
  simulated child failure reports loading as unknown without killing the host.
- [x] C/ISO_C_BINDING Fortran hosts query4 versions per ordinary/HSE library;
  values agree and MPI stays uninitialized. Mandatory source/all-linked-object
  audits precede helper links. No solver states or GS data are created.
- [x] Host compile found missing HSE checkpoint/recovery query declarations in
  the handle header; add declarations only, without changing library functions.
- [ ] R4 remains open for full model/state/step lifetime conformance and stable
  ABI negotiation. Protocol-version matching alone is not physics/profile,
  parameter-signature or checkpoint/topology certification.
- The explicit probe loads native code and may run loader constructors: trusted,
  audited libraries only. Static check_job remains strictly separate/nonexecuting.

- R4 probe evidence: source and installed641 tests (640pass/1skip); standard/HSE
  C/Fortran query-only hosts and installed HSE host runner pass. All return1 and
  MPI remains uninitialized. Installed CLI/doctor/API/header/C fixture match source.
  Canonical record: calc2026/salmonpy-r4-runtime-probe-20261120; cloud sync
  unconfirmed. No solver state, GS data or Actions run. Still uncommitted.

### 2026-11-20 — R4 standard-PZ functional C/Fortran host lifecycle

- [x] C/Fortran MPI1/2 paired-state lifecycle matches independent controls; bad-old-A rejection preserves accepted state, stale generation remains rejected after slot reuse, fresh state reproduces its baseline, caller MPI survives detach.
- [x] Reuse common Darwin MPI link policy in host runners. Same Fortran object relinked with shared MPI COMMON sentinels fixes the host-only failure; no numerical solver change. Record link commands and SHA256.
- [x] Common directory-scope C declarations installed through native_context.h. Existing GS read in place and byte-checked unchanged; no GS copy/electronic checkpoint/Actions run.
- [ ] Full R4 remains open: fast-HSE functional hosts, additional profiles/platforms and stable ABI negotiation.
- Evidence: calc2026/salmonpy-r4-state-host-confirmed-20261120 (cloud sync unconfirmed).

- Verification: source and sdist-built installed wheel each ran641 tests,640pass/1skip. Installed C/Fortran host lifecycle and independent controls pass at MPI1; source matrix passes MPI1/2. Installed fixture/header discovery and common MPI link helper verified. 日英APIとAI JSON同期検査通過。未コミット。

### 2026-11-20 — R4 fast-PR HSE functional host lifecycle

- [x] Reuse the C/ISO_C_BINDING Fortran lifecycle fixtures for fast-PR HSE via a compile-time profile; no duplicated host programs or numerical solver changes. HSE signatures/status codes remain distinct from ordinary TDDFT.
- [x] Full-k HSE MPI1/2: paired currents match independent controls and C matches Fortran; sibling states remain unchanged, bad-old-A preserves accepted state, stale handles remain invalid after slot reuse, fresh state matches the original baseline, caller MPI survives detach.
- Existing GS read in place, byte-checked unchanged; no GS collection or electronic disk checkpoints. Record: calc2026/salmonpy-r4-hse-state-host-20261120 (cloud sync unconfirmed). No Actions dispatch.
- [ ] Stable ABI remains open: symmetry-reduced HSE functional hosts, more XC/platform profiles and negotiation coverage. This is state-lifecycle interoperability, not physical convergence certification.

- Verification: source and sdist-built installed-wheel suites each641 tests (640pass/1skip). Installed HSE C/Fortran MPI1 lifecycle/independent controls pass; source HSE matrix passes MPI1/2. Common-fixture ordinary PZ MPI1/2 also passes; all8 traces byte-identical to the preceding ordinary-PZ run. JP/EN/AI guide sync and diff checks pass. Changes remain uncommitted.

### 2026-11-20 — R4 symmetry-reduced fast-PR HSE host lifecycle

- [x] Common C/ISO_C_BINDING Fortran host fixtures select symmetry-aware initialization without duplicating programs. Runner --state-symmetry references existing sym.dat and drives the fixture's retained z direction; full-k xyz driving stays unchanged.
- [x] Symmetry-aware PT-CN/learned exchange,8 full k points reduced to4 representatives: MPI1/2 C/Fortran paired currents match independent controls; bad-old-A preserves state; stale generation rejected; fresh state reproduces baseline; caller MPI survives detach. Runtime logs confirm the orbits algorithm and representative k points.
- [x] Existing GS, PP and symmetry data read in place and unchanged; no GS collection or electronic checkpoints. Evidence: calc2026/salmonpy-r4-hse-reduced-host-20261120 (cloud sync unconfirmed). No Actions dispatch.
- [ ] Stable ABI still requires remaining negotiation/profile/platform coverage. This fixture is not physical convergence certification or a general all-direction symmetry guarantee.

- Verification: source and sdist-built installed-wheel suites each641 tests (640pass/1skip). Installed reduced-HSE C/Fortran MPI1/2 all8 lifecycle/control cases pass; logs confirm4 representatives. PZ/full-k HSE MPI1 controls,4 traces each, remain byte-identical to preceding runs. Wheel guides/AI JSON/host fixtures/runner match source bytes; sync/diff checks pass. Changes remain uncommitted.

### 2026-11-20 — commit/push queue reconciliation

- Completed in this commit: native E/H checks only active field settings and skips inactive electronic preprocessing; current-grid5% screening scope/FFTW dependency documentation; declared native contracts, explicit trusted-runtime protocol probes, common directory header, C/Fortran query and functional lifecycle hosts, shared Darwin MPI link policy.
- R4 measured lifecycle: standard PZ and fast-PR HSE full/reduced-k, Mac MPI1/2, independent controls and rejected-call state preservation. JP/EN/AI guides and packaged host fixtures synchronized. Historical 'uncommitted' notes above describe their recording time and are superseded by this commit.
- Next active work remains group5: R3 later initialization/allocation failure returns/cleanup, then R4 remaining operation/negotiation and platform coverage. Group3 NLCC runtime certification and legacy unit audit follows. Group1/2 deferred; group4 portable/packed migration and group6 scaling/distribution are separate extensions.
- Routine Actions remain disabled by commit skip marker; no new GS collection or production-job changes.

### 2026-11-20 — R3 early native E/H initialization returns

- [x] Existing grid/boundary/time/CFL failures return1101–1104 with optional status; standalone termination remains unchanged. Clean up communication caches before owned children, destroy failed context and permit corrected native retry. No new restrictions or numerical changes.
- [x] Old build reproduces grid-error host termination. MPI1/2 repeats all4 failures twice; contexts/live children absent after returns, registry slots reused; successful32-step A/E/H and .data outputs exactly match the previous build.
- [ ] R3 still open for failures inside grid/MPI setup and later specialized initialization, I/O and allocation stages; arbitrary STOP/OOM/MPI recovery is outside the guarantee.
- Record: calc2026/salmonpy-r3-geometry-record-v2-20261120 (cloud sync unconfirmed). No GS data collected; no Actions dispatched.

- Distribution verification: source and sdist-built installed-wheel suites each641 tests (640pass/1skip). Installed FDTD MPI1/2 repeats the11 admission attempts and normal32-step propagation; all A/E/H samples and .data bytes match the unchanged reference. Standalone executable preserves all4 original diagnostics. Installed solver/host test/JP/EN/AI bytes match source. Native audits and guide-sync/diff checks pass. Changes remain uncommitted.

### 2026-11-20 — external review follow-up

- [x] Separate explicit native_host from optional returned status; default/explicit false preserve ordinary validation. GS-free MPI1/2 regression covers all3 routes.
- [x] Add existing-GS-only PZ host fixture preparation, small GS input starter and exact theory/files/decomposition documentation. No GS generation/copy/download; host MPI1/2 controls pass.
- [x] Expose verified from_runtime delegation, move the simple FDTD example to load_runtime, provide FDTD reading order and actionable existing-output messages. Existing load_runtime callers retain their entry.
- [x] One-time legacy-settings fallback warning; explicit owners stay quiet, compatibility stays allowed. MPI1/2 tests use no GS or numerical state. Keep new instance state in contexts/slots.
- [x] Expose paired per-axis periodic boundaries and native neighbor communication, retaining centered coordinates. Dielectric/LD slab Fresnel references run in C; MPI1/2 and transverse-probe comparisons pass (max relative complex error0.0554%/0.0375%). Existing isolated32-step PML/LD output matches previous bytes.
- [x] Fix three whitespace typos and accept build-record directories with clear missing-build.json diagnostics.
- [x] Configure required inventory + Git-only build + GS-free ABI/parser CI for critical paths; no skip marker for those commits. Keep Linux sdist automatic, heavy OS/physics campaigns PR/manual. No Actions dispatched here; CI outcome remains unmeasured until a push/run.
- [ ] Common grid/MPI setup minimum-domain and later initialization failures remain R3 work; thin MPI trial used too few cells and hit the existing grid initializer STOP. Retain the failure record; do not claim arbitrary init recovery.

- Final review evidence: Git-index-only export inventory8/8, prepare/build and C/Fortran query contracts pass. Source/installed suites641 tests (640pass/1skip), final alias/output-attribute targeted checks pass; wheel critical sources/JP/EN/AI/examples match bytes. New PZ fixture MPI1/2 all8 cases, GS-free parser/warning MPI1/2, periodic dielectric/LD Fresnel MPI1/2, old isolated FDTD and fast-HSE current parity pass. CI YAML checked locally; no remote Actions run or GS collection. Record: calc2026/salmonpy-review-followup-record-20261120 (cloud sync unconfirmed). Uncommitted; next critical commit must not use skip-ci.

### 2026-11-20 — optical-feature request and CI scope

- [x] Record the new feature request separately from the preceding bug/operations review. Physical reproduction/convergence is not required automated CI. NLCC/TDDFT/full-runtime campaigns are manual-only; retain automatic inventory, tracked-only build, packaging and GS-free host/API checks. No Actions dispatched.
- [x] Transverse-periodic FDTD/Fresnel validation from the earlier review is already implemented and recorded; do not schedule its physical reproduction in CI.
- [ ] New optical queue: public TB-mBJ interface -> conductivity/dielectric observable contracts -> linearity/sum-rule diagnostics -> epsilon_infinity in native LD -> constrained fitting/provenance -> time/window sensitivity, post-hoc k contribution and measured cost advice. See 2026-11-20-optical-observables-plan.md for dependencies and exclusions.
- [x] Subpixel smoothing has completed its LD-only reference-frequency implementation and local synthetic validation (2026-11-21; details below). It was a separate design/validation phase after the periodic reference. Do not promise second-order convergence or a coefficient-only implementation before checking tensor/dispersive requirements. FEM/DG replacement, automatic nstate/k planners and copyrighted optical tables are excluded.

- [x] 2026-11-20: public serial TB-mBJ DFT/cval, ASE structure bridge and material JSON; coefficient preserved through GS/TDDFT/archive and restart checks. Synthetic native XC probes agree bitwise with Libxc ON/OFF. No new GS collected; real-material precision and native parallel TB-mBJ SCF are not claimed. Optical observable work remains next.

### 2026-11-20 — HPC readiness review (baseline 05827a7e)

The external Linux review reports host-MPI evidence up to five ranks. Larger-scale correctness, site MPI compatibility and scaling remain unverified; an inherited SALMON kernel is not by itself proof for the embedded Python workflow. These are queued tasks, not newly validated capabilities. Finish the current serial TB-mBJ work first; prioritize job-scoped validation receipts before GPU work.

- [x] **Job-scoped runtime validation receipt core:** explicit job/token-bound RuntimeValidation, one full preparation and authenticated metadata reuse, executable SHA256 retained, invalid proof/changed metadata warning with complete-verification fallback. Local POSIX MPI4/data-only58MiB evidence, byte/stat/walk/time measurements; see2026-11-21-job-runtime-receipts.md.
- [ ] **Receipt site validation:** shared parallel-filesystem freshness/inode identity, multi-node/site MPI and higher-worker-count evidence remain unverified. No local-cache timing extrapolation to FEFS.
- [ ] **HPC launch modes:** document and expose (a) preparation/analysis-only Python with the site's job script launching SALMON directly, and (b) a collectively called mpi4py/embedded path on caller-owned communicators. Keep nested launcher protection; removing the MPI-rank guard alone is not a solution. Preserve caller MPI after detach. Support site-specific argument construction rather than assuming every launcher accepts `-n`.
- [ ] **Multi-node and multi-rank evidence:** audit tool-only rank limits separately from product limitations; parameterize appropriate tests without promising unsupported layouts. Test k-point/orbital/spatial combinations, independent state ownership and checkpoint repartition on site MPI and multiple nodes. Record compiler/MPI/node layout, correctness and scaling separately. Large physical campaigns remain manual/site jobs, not required GitHub CI.
- [ ] **Observation I/O at scale:** measure output volume and file count versus ranks, points and stride; design optional aggregation/buffering or sharded output while preserving observation metadata and existing small-run interfaces. Do not assume byte volume alone measures filesystem cost.
- [ ] **GPU prerequisites (after receipt work):** declare validation per feature/backend/precision rather than a single GPU-supported flag; assign device allocations and cache ownership to each slot/model, define stream/concurrency behavior, and set precision-specific comparison tolerances before implementation. Audit GPU eigensolver call paths under the existing prohibited-driver policy. Preserve explicit rejection of known unimplemented numerical paths without blocking unrelated custom backends.
- [ ] Reconcile this review with current code before implementation: its 05827a7e timing and rank observations are baseline evidence, not measurements of later changes. No multi-node or GPU validation is claimed on this Mac.

- [x] 2026-11-20: single-polarization optical_conductivity/dielectric_function, explicit total-E/charge-current/Gaussian-AU contract, masked weak-drive/DC bins; native Fortran ratios and finite-band/linearity diagnostics on supplied data. Synthetic phase, sign, amplitude scaling and finite-band integral tests; no material accuracy certification or extra solver runs.
- [x] Native LD epsilon_infinity completed; fitting/provenance work is recorded below. HPC remains after the optical queue.

- [x] 2026-11-20: LD epsilon_infinity (default1), native continuous permittivity/D/background-energy and FDTD coefficients; JSON/provider/config/checkpoint identity propagation. Preserve pole P/J and default update order. Non-unit-background legacy impulse spectra warn/omit without stopping evolution. Source synthetic and native MPI1/2 Fresnel/default parity/direct checkpoint checks passed; real-material certification is not claimed.
- [x] Constrained LD fitting/provenance completed below; remaining optical work is duration/window/k post-analysis, then subpixel design. HPC remains deferred.

- [x] 2026-11-20: native passive LD fitting, direct epsilon/nk/masked-optical inputs, explicit pole counts and deterministic starts, optional initial model/background, residual/optimizer separation, data/file SHA256 and literature reference, calibrated band and warning-only excitation checks. Synthetic tests cover fit accuracy, passive failures, source identity and input independence. Local optimization and finite-band calibration do not certify a global optimum or real materials.
- [x] Superseded/completed on2026-11-21: duration/window/k diagnostics and LD-only reference-frequency subpixel modules; see their dedicated plans.

- [x] 2026-11-20 user explicitly authorized acquiring/retaining legally redistributable experimental optical data and fitting a reusable list. This supersedes the earlier exclusion of bundled Johnson-Christy numerical data for the selected CC0 database subset only. Preserve original measurements, license, pinned commit, SHA256, DOI and conditions; no journal articles/figures. Five materials,22 retained candidates; choose defaults using residual and time-resolution screening, retain warnings and optimizer nonconvergence. All experimental refits stay explicit/local, not required CI.

- [x] CC0 subset data library:5 sources,22 stored fit trials,4 screened defaults. Experimental comparison exposed Ag local weighted errors up to38.6% despite4.34% global RMS; mark needs_refinement, retain warning-only explicit reuse. Al stiff low-pole fits stay available but are not defaults. Preserve raw/license bytes across platforms with attributes.
- [ ] Improve Ag local-shape fitting/optimizer termination before presenting it as a default; do not relax/hide screening to obtain a pass. All new experimental-data matching/refits remain explicit local tools, not required CI.

## SciPy optical fitting (2026-11-20)
- Implemented optional optimizer="scipy" and fit extra; native objective/seed/statistics, SciPy bounded TRF and derivative estimates. No import-time install or fallback.
- Local CC0 study:44 SciPy trials plus22 retained native trials; selected SciPy fits locally terminated, Ag remains needs_refinement. Smallest screened defaults Au6/Cu6/Al4/Si6. Raw data preserved.
- Source and sdist-to-wheel installed validation:649 tests passed (one optional skip) in each; final artifact/record checks follow. No experimental matching in required CI.

## Optical oxides and tensor direction (2026-11-20)
- Added pinned CC0 scalar n,k films SiO2/Al2O3/TiO2/Ta2O5/HfO2 and22 SciPy trials (88 retained total); preserved prior five-material fit records. No n-only zero-absorption substitution. TiO2 residual passes but optimization-budget warnings remain.
- Common tensor representation in ALL modes supersedes3D-only scope. Implemented immutable measured/model(...,3,3), native isotropic embedding, compressed LDTensorModelConfig, FDTD/standalone/MPI-provider acceptance. Unknown directional response uses one independent scalar on the diagonal.
- Source tests651 pass, one optional skip. Audited unchanged FDTD binary:scalar/tensor inputs and observations bitwise equal for xyz polarizations, MPI1/2 (six paired cases). Wheel validation follows.
- REMAINING:independent diagonal axes/off-diagonal native propagation and fitting/rotations across modes; do not claim general anisotropy from isotropic tensors.

## Reciprocal tensor material response (2026-11-20)
- Implemented native PSD/SPD matrix validation, tensor epsilon, D/background energy and persistent coupled P/J updates. Isotropic kernel path retained. TensorLorentzDrude supports separate diagonal axes and reciprocal off-diagonal couplings; no automatic rotation or tensor fitting.
- Direct process-local snapshots restore without replay; overflow rejection preserves accepted P/J/E and permits retry. Independent states do not interfere.
- Local source suite653 tests passed, one optional skip. Actual MPI3 with two different tensor points:40 endpoint steps bitwise equal to independent native controls; caller MPI survives.
- REMAINING: general anisotropic native FDTD field/coefficient update with Yee-grid interpolation and checkpoint/schema wiring across field modes. General tensors are rejected by existing native FDTD rather than projected away. Source-to-sdist-to-wheel installation verified:653 tests passed, one optional skip; wheel MPI3/two-point controls and direct snapshot continuation match. Previous isotropic kernel responses match bitwise in three500-step controls.

## Native diagonal E/H and waveplates (2026-11-21)
- Principal-axis diagonal tensors connected to native E/LD coefficients, preserving P/J/Q ordering and existing history arrays. Explicit ld-diagonal input asset, hash binding, matching-runtime check and both batch/host-MPI staging. No new process-global material arrays.
- First prototype:old/new scalar, diagonal isotropic limit and MPI1/2 outputs bitwise equal. Dispersive distinct-axis controls match separate scalar native FDTD to below3e-18 absolute field difference.
- Synthetic quarter/half-wave plates:normalized |S3|0.99669 and S2=-0.999697 respectively; Stokes-norm discrepancy vs Fresnel slab0.00305/0.00563 at spacing8. This is not real-material or complete convergence validation. Refinement verified:error0.000760 at spacing4 (4.01-fold reduction). Git-only739-file export/build and source/installed-wheel655-test suites passed (one optional skip each). Native MPI2 direct checkpoint A/E/H and RT outputs match continuous execution bitwise, invalid/missing assets return1110 collectively. Final record:salmonpy-diagonal-waveplate-record-20261121 in the canonical results store.
- Native carrier/analytic Stokes diagnostic and explicit --run sample added. Off-diagonal Yee field coupling, automatic rotations and all-mode propagation coverage remain.

## Three-drive tensor extraction/fitting (2026-11-21)
- [x] Native matrix inversion of three independent J/E spectra; DC/weak/rank-deficient bins masked. Native PSD-factor tensor fitting with common poles, fixed initial SPD background, optional SciPy, input/model provenance, source and installed-wheel657-test suites passed (one optional skip). Synthetic negative off-diagonal roundtrip weighted RMS below1e-14. Updated paired LD guides, AI task and ld_fit --tensor example.
- [x] Reciprocal off-diagonal native Yee coupling, dispersive cross-drive, model/geometry-bound checkpoints and rotated waveplate/MPI validation completed below. See2026-11-21-offdiagonal-fdtd-tensor-fit.md; real-material certification is separate.

## Reciprocal off-diagonal native E/H completed (2026-11-21)
- [x] Full symmetric SPD background and PSD shared-pole LD matrices on material-owned Yee edges; paired means and transpose means, implicit E PCG, pole J before original P/Q order; native numerical work. Old scalar and diagonal results remain bitwise equal.
- [x] Finite iteration exhaustion warns/continues. MPI numerical rejection rolls back and permits retry. Checkpoint compares matrices/rates/edge ownership/weights/grid/dt/options before commit; late-file and model failures preserve all state. Material on only one MPI rank does not deadlock. Continuous and direct-restarted A/E/H and raw observations are bitwise equal; no replay.
- [x] Synthetic rotated30deg quarter/half waveplates agree with Fresnel Stokes within0.00305/0.00563; refined quarter0.000760. Two-pole rotation errors below3e-15; all-pair3D axis covariance below2e-21. Installed synthetic three-drive epsilon->tensor-LD->FDTD passes. No new real-material TDDFT/GS claim.
- [x] Source and wheel658 tests, one optional skip each; Git-only739-file audited build and wheel source parity; paired guides, AI JSON and existing waveplate/recovery examples updated. No remote CI started. Evidence and audited local runtime retained in canonical calc2026/salmonpy-offdiag-final-record-20261121.


## LD-only interface preparation and optical A4 (2026-11-21)
- [x] LD-only center-frequency subpixel material coefficient initialization;
  native module/shared copies, compatible runtime/assets, response-point exclusion,
  passive tensor lift and unchanged ADE order. Default disabled preserves legacy.
  Fixed-time-step curved refinement, homogeneous/slab/MPI2 checkpoint evidence
  recorded in2026-11-21-reference-subpixel-smoothing.md and canonical results.
- [x] A4 observation-only correction in Fortran; independent raw/corrected
  conductivity and dielectric spectra, masks, fit band/coefficients/assumptions.
  Native solver/current feedback unaffected. Paired manuals and AI JSON updated.
- [x] Duration/window/post-hoc k analysis implemented in the Fortran observation kernel (2026-11-21). HPC
  validation receipts and launcher/site-scale work retain their earlier queue.


## Stored optical trajectories (2026-11-21)
- [x] Common-frequency prefix/window comparison and already weighted k-current
  reduction with optional total-charge-current consistency check; no solver,
  no new GS, no k scan or implicit extra weighting. Public analysis APIs,
  paired guide/AI JSON and NPZ sample connected. Native sums and Fourier work.
- [ ] Remaining operational queue: receipt site validation, HPC launcher
  modes, multi-node evidence and observation-I/O/large-data memory advice.
  Larger material-accuracy studies and Ag fit refinement remain separate.


## Runtime validation I/O (2026-11-21)
- [x] Optional per-job receipt with caller-owned secret/ID, immutable bundle
  fingerprints and mandatory executable hash; no global cache, source/object
  startup audit, solver launch, MPI initialization or automatic build.
- [x] Passed explicit context through process TDDFT provider/backend; unchanged
  material identities/current algorithm, existing NLCC capability checks preserved.
  Default and unsupported metadata platforms retain full verification. Custom
  backends and native build/library loaders are not forced into this mechanism.
- [x] Local MPI4/three reuses per rank, unchanged executable checks, no tree walk,
  safe fallback/corruption rejection and I/O measurements. Paired runtime guides,
  AI JSON and preview-first validation-only MPI example updated.
- [ ] Next core HPC task: site-launch argument construction/preparation-only
  versus caller-owned MPI entrypoints. Receipt multi-node/FS proof stays separate.


## Minimize outer communication (2026-11-21)
- [x] Project policy recorded: no unnecessary world diagnostics/health barriers;
  keep physical sums/halo exchange and consistent failure/retry guarantees.
- [x] Receipt-only example setup bcast + root-only final gather; no replicated
  allgather, health allreduce or barrier. Native E/H Python admission/result
  agreements coalesced6->2 within the field group, native internals unchanged.
- [x] Owner-only A/field-packet routing implemented (2026-11-21);
  preserve physical timestep and mixed snapshot/retry semantics; reduce duplicate
  outer agreements only when consistent failure ordering can be proven.


## Owner-only response inputs (2026-11-21)
- [x] Numeric A Scatterv/per-owner packet scatter, only local material subgroup
  broadcasts, no data to empty field group or singleton rebroadcast. Cached
  point/count/permutation mapping; root order/interleaved/different counts covered.
- [x] J owner gather/order preserved. Global checkpoint array shape retained,
  material A history updates owned rows only; old full-array checkpoints accepted.
- [x] MPI6 custom routing, MPI4 LD and native field+LD+two SBE points, MPI3 packet
  admission/failure, MPI6 protected mixed step and old/new direct restore checked.
  No new GS or real-material TDDFT/HSE/multinode claim.
- [x] Native field spatial J ownership implemented; see the spatial-current routing plan.
- [x] Adjacent outer field-status/response-admission votes coalesced; retain pre-field, pre-H and post-material decisions. See the coalesced outer admission plan.


## Spatial current routing (2026-11-21)
- [x] Native owner map from actual Yee deposit edges; boundary duplication and
  zero-point field ranks; packed J step ABI with legacy full-input fallback.
- [x] CoupledStepper/PacketCoupledStepper route owned J within the field group;
  reporting leader retains globally ordered output/checkpoint J.
- [x] Local MPI4/6 owner and boundary tests, preceding-binary/legacy-path bitwise comparisons,
  owned/full-row direct resume and protected failure/retry; source/wheel668 tests
  (one optional skip). See 2026-11-21-spatial-current-routing.md and canonical record.


## Coalesced outer coupling admission (2026-11-21)
- [x] A coupling outer agreements4→3; packets without H7→4 and with H7→5.
  These are hierarchical agreements, not individual MPI calls or speedup ratios.
- [x] Explicit upstream-error argument; standalone signatures and custom public
  response overrides retain their previous path. No persistent error flag.
- [x] H-sampling barrier decision, recipient validation, native reductions and
  protected-step transaction votes retained. Fault-injection tool follows owned-J
  field entry and skips native-electronic probes if no native electron owner exists.
- [x] Source/wheel668 tests (one optional skip), MPI4 native late/deep failure
  preservation/retry and direct resume, MPI6 custom protected mixed states,
  packet pre/post-H failures and rank-local custom response dispatch verified.
  Evidence in the dedicated plan and canonical outer-admission record.


## Status-bearing A/J transport (2026-11-21)
- [x] Managed coupling exchanges input status with owner A/packets and gathers
  material status with J; field returns carry status plus owned current.
- [x] Reporter-only explanatory Bcast occurs only after every rank knows the
  error origin; normal steps never use the owned failure-only communicator.
  Caller input-communicator lifetime and MPI finalization remain independent.
- [x] Separate hierarchical votes A3→1, packet no-H4→2, H5→3; leader collective
  calls5→4,6→5,7→6. Keep post-J acknowledgement, pre-H and transaction votes.
- [x] Final source/wheel668 tests (one optional skip), MPI4 protected numerical
  failure/retry and direct resume, MPI6 custom states and MPI3 freed-parent/root
  order validated. Normal reporter Bcast0; failed A/J Bcast1 from selected rank.
  Seventeen A/E/J comparisons bitwise equal; evidence in the status transport plan.
- [x] Reported lightweight CI failure repaired in test-only commit8d5d90f5;
  run44852912505 passed. No manual platform/physics matrix dispatched.


## HPC FDTD operational path (2026-11-21)
- [x] prepare_fdtd writes only native input/assets plus prepared-not-run manifest;
  no executable, subprocess, MPI initialization, compilation or submission.
- [x] FDTDMPI.command/rank_arguments construct site argv without shell evaluation;
  preserve default -n and existing nested-launch rejection.
- [x] Explicit two-rank local launch with -np matches the old -n path: identical
  input bytes and four native FDTD/LD .data files. No GS collected.
- [x] Source/wheel669 tests (one optional skip), preparation-only installed CLI,
  two-rank default/custom native launch and tensor sidecar parity verified.
  Dedicated plan and canonical record retain the bounded evidence.
- [ ] Remaining HPC scope: generalized initial-guess staging and real site
  launcher/runtime/multi-node evidence. Fresh DFT/RT staging is complete below.
  General init/ABI/NLCC work remains in its own queue.


## DFT/TDDFT input-only preparation (2026-11-21)
- [x] prepare_dft stages fresh native input plus PP/symmetry bytes, preserving
  default GS saving and existing run/warm-start/HSE execution contracts.
- [x] prepare_tddft streams only supplied native input/explicit local files;
  no GS discovery/snapshot, namelist changes, worker, MPI or job submission.
- [x] Short native DFT input/status/forces parity with restart output disabled;
  TDDFT current/energy parity reading existing GS in place. No new GS records;
  existing GS hashes unchanged. This is parity evidence, not physical convergence.
- [x] Source/wheel671 tests (one optional skip), management-only CLI, preceding
  wheel DFT parity and native RT current/energy parity verified. Existing GS
  unchanged and no new GS records collected. Dedicated plan/record retained.
- [ ] Actual site/multi-node launch evidence and generalized warm-start staging
  remain separate. Group5 R3/R4 and group3 are not closed by this work.


### 2026-11-21 SCF未収束時の確認待ち

- [x] バッチDFTとnative GSの反復上限未収束をデータ保持＋確認待ちへ変更。
- [x] 残差CSV・任意PNGと警告を案内。意図的な部分SCFは通常継続可能。
- [x] 日英マニュアル・AI向けJSONの自動継続案内を更新。
- [x] 実計算で従来出力との不変性を確認。Si8の8反復・6種類の出力が従来wheelとバイト一致。GS保存なし。記録: calc2026/salmonpy-scf-review-20261121/record.json。


### 2026-11-21 FDTD energy observations

- [x] Native per-macro-point charge j dot total E work, separated E/H R/T and unwindowed fluences.
- [x] FDTD budget with volumes/area/field storage/other outflow; warnings without stopping runnable solver.
- [x] Uniform-port scope, reversible-storage accounting and no extra MPI communication documented in JP/EN/AI catalog.
- [x] Nonuniform transverse-periodic plane area integration added via area_weights/surface_flux; fresh batch native plane collection/Fortran colocation connected; embedded/restarted streams remain separate.


### 2026-11-21 Real coupled energy-budget validation

- [x] Commit/push 58edd163; automatic lightweight LD workflow 44864514394 succeeded.
- [x] Actual LD+E/H FDTD and vacuum control. Residual 0.0108%; half-spacing/dt 0.00270%.
- [x] Reconstructed LD current verified against native saved final current (rtol1e-10/atol1e-18).
- [x] Actual 2-point/8-rank TDDFT+Weyl FDTD, source GS reused unchanged. Delayed incidence residual 0.000711%.
- [x] Prefilled-wave initial-energy accounting, finite coverage caveat and software/physical scope documented.
- [x] Reusable source-only observation analysis tool and C reference; no physical runs added to CI.
- Further precision search stopped at user's accepted accuracy.


### 2026-11-21 Periodic transverse-plane energy

- [x] Native local cross product before area integration in time and frequency.
- [x] Existing R/T/budget and NPZ sample accept nonuniform planes; uniform path preserved.
- [x] Analytic phased diffraction (mean E=0 but finite flux), patch additivity and finite checks.
- [x] Real native lossless periodic grating: area-budget residual0.0107%; field means miss10.1% of incident energy.
- [x] JP/EN/AI guidance distinguishes full plane integration from automatic collection and order separation.


### 2026-11-21 Automatic batch native plane acquisition

- [x] Commit/push ac0354f7 for prior plane integration.
- [x] Owner-local native plane+neighbor streams; no new collective, no numerical update changes.
- [x] PlaneMonitor/PlaneData and Fortran3-axis spatial colocation; H clock preserved.
- [x] FDTDResult reference R/T and energy budget; missing diagnostics preserve completed calculation.
- [x] Real grating128 point files and EM checkpoints byte-identical to prior executable; serial/transverse MPI2 budgets agree.
- [x] Normal split with planes on different ranks: nonzero colocated E/H bit-identical.
- [x] Live source manifests/profiles and JP/EN/AI/sample synchronized.
- [x] Source and installed sdist→wheel suites:678 tests, one optional skip; installed wheel reproduces serial/MPI budgets. Native diagnostic I/O failure preserves completed calculation.
- [x] Explicit batch restart stream reading/joining: serial/MPI2 checkpoint1200→2400 E/H bit-identical, point values equal, final EM checkpoints byte-identical, energy budget unchanged.
- [x] Installed-wheel suite678 tests, one optional skip; restart-segment budgets reproduce source results. Existing mock startup timeout isolated from injected step timeout without changing production timeouts.
- [x] Embedded explicit read-only plane ABI and committed-only sampling: one field-group gather, immutable rank0 snapshots, no automatic global communication/native writes.
- [x] Pending step/restore candidates cannot sample or publish. Optional no-overwrite archives warn on failure; existing observations and native state stay usable. Mixed sample opts in with --planes-x/--plane-stride and saves one final archive.
- [x] LD serial/MPI2 monitor off/on, protected rollback/retry and direct restore preserve traces and point-output bytes; prior batch E/H and flux agree (max field difference below8e-19). Active checkpoint payloads agree. Source/installed679 tests, one optional skip; installed native MPI2 matches source.
- [x] Unused fe%coo padding now zero-initialized before unchanged coordinate arithmetic. Old NaN/Inf poison test fails, repaired GNU probe passes periodic/odd/even cases. Serial/MPI2 currents, E/H and point bytes match prior runtime; full checkpoint bytes now match after protected/direct restore on every rank. Historical validation records retain the earlier padding difference. Source/wheel679 tests, one skip; installed SDK probe and native MPI2 full-byte comparisons pass. Existing CI gains only a small synthetic check, with no manual dispatch. Batch restart launching remains explicit native/site input.


### 2026-11-21 R3 local MPI grid initialization return

- [x] Commit9eb6f2d6 records plane restart/committed embedded observations and deterministic coordinate padding; no push requested/performed here.
- [x] Optional status in init_grid_parallel recognizes the same existing minimum across its already replicated domain table; no extra collective/new restriction. Native E/H propagates1106 before cache creation; standalone status-absent behavior unchanged.
- [x] Explicit sp_eh_init_recovery_protocol1 + unanimous1106 lets Python release the failed native lifetime reservation after detach. Corrected constructor/new work directory retries on the same caller communicator.
- [x] Old validated binary reproduces process STOP (exit0 with no validation record). New GNU MPI1/2 tests repeat low-level and public rejections twice, including last-rank-only short domain, release all children, retain parent MPI, then complete32 steps/direct resume. Old normal A/E/H arrays and point files match exactly.
- [ ] Other MPI mapping/setup, allocation, specialized initialization and R4 ABI stability remain open; no arbitrary STOP/OOM/process/MPI recovery claim. Groups1/2 and molecular work remain deferred.
- [x] Source/installed-wheel679 tests, one optional skip; installed SDK MPI2 repeats public/raw initialization retries and direct resume, matching source arrays/output. Normal midpoint checkpoint bytes also match old controls on every rank. JP/EN/AI synced.
- Records:calc2026/salmonpy-fdtd-grid-init-recovery-20261121; no GS collection or manual Actions dispatch.


### 2026-11-21 R3 explicit MPI layout/mapping returns

- [x] Old validated runtime reproduces host termination for explicit unsupported Maxwell MPI factors, with exit0 and no completed validation record.
- [x] Optional status in init_communicator_dft returns from actual unimplemented mapping/count branches; preserve status-absent legacy behavior and Fujitsu/custom branch selection. No outer mapping whitelist added.
- [x] Native E/H forwards1105/1107 through the existing cold cleanup/protocol, alongside1106. Python releases failed attachment/reservation only for unanimous known statuses/protocol1. No new collective or valid-input restriction.
- [x] Existing E/H regression tool consolidates27 admission/retry checks per GNU MPI1/2 layout, child/caller lifetime, normal32-step/direct-resume comparisons; old A/E/H, point bytes and checkpoints exactly match. Existing alternate orbital mapping remains usable on MPI2.
- [ ] Automatic layout factorization, Fujitsu topology, later allocation/IO/specialized initialization and R4 stable ABI are still open. Groups1/2, molecular work and group3 physical certification are not implicitly resumed.
- [x] Source/installed-wheel679 tests, one optional skip. Installed SDK MPI2 repeats27 checks and direct resume; fields, point bytes and checkpoints match source/old controls. Paired manuals/AI JSON synchronized.
- Records:calc2026/salmonpy-fdtd-mpi-layout-recovery-20261121; no GS collection or manual Actions dispatch.

### 2026-11-21 R3 material-ID shape input recovery

- [x] Native owned-unit streaming CUBE reader returns1115 for supported IO/data/grid failures,1116 for reported tiny staging allocation failure. Preserve conversion and standalone path; scalar finite checks include partial records.
- [x] External-input initialization only: one field-group agreement before halo exchange; partial-rank read failure returns collectively. Known cold failures release cache/children/reservation for retry.
- [x] GNU MPI1/2 retries/direct resume match prior generated-shape field, point and checkpoint bytes. Host unit600 remains owned through rejected input attempts.
- [x] Ordinary point/plane observation and generated-shape output unit ownership addressed below. General OOM, device/MPI/process failure and specialized initialization remain open.
- Records:calc2026/salmonpy-fdtd-shape-input-recovery-20261121. No GS collection or manual CI dispatch.
- [x] Source and isolated sdist→wheel suites679 tests, one optional skip. Installed SDK MPI2 repeats the full shape-input rejection/retry test and exact source field/point/checkpoint bytes. JP/EN/AI synchronized.

### 2026-11-21 ordinary E/H observation file-unit ownership

- [x] Reproduce host-owned unit600 closed by the prior validated runtime after normal32-step/direct-resume execution; preserve failed-run evidence.
- [x] Local newunit ownership for point/plane headers, point append, final observation information, plane snapshots/integrals and material-ID CUBE writer. Preserve helper signatures, filenames, formats, loop order and communications. Unimplemented shape-output formats no longer close a caller-owned unit.
- [x] Audited maintained source and all application objects; Git-only inventory export prepared and built in a new external directory. Preserve old runtime and records.
- [x] Legacy restart/helper units isolated below; linear-response/spectral/ASE/ART and TTM file units remain to audit/isolate. No claim of general output IO-error/OOM recovery; Fujitsu execution is not available here.
- Records:calc2026/salmonpy-fdtd-output-unit-20261121. No GS collection, push or manual CI dispatch.
- [x] GNU MPI1/2 preserve caller unit600 throughout normal initialization,32 steps,direct restore/finalize. Prior fields/checkpoint/shape bytes and point/plane outputs match exactly (77 data files with planes enabled). Source/isolated-wheel679 tests,one skip; installed SDK MPI2 matches source. JP/EN/AI synchronized.

### 2026-11-21 common Maxwell restart/helper file-unit ownership

- [x] Prior audited runtime passes unoccupied synthetic control but closes a host-owned unit in the occupied variant (STOP81). Preserve failed run.
- [x] Use owned newunit for text/binary single/all restart data and metadata and remapping reads; copy helper owns its two units. Keep signatures/paths/bytes/array order/communication; no new collective.
- [x] Existing native I/O fixture consolidates real/complex1–5D, text, same-layout/cyclic target-domain remapping, inactive payload and copy checks. GNU MPI1/2 preserve units600/601/700/701; control/occupied and prior-runtime bytes identical.
- [x] New Git-only export prepared/built with source and application-object prohibition audits. Normal LD+FDTD/direct-resume MPI2 is unchanged; no GS collection.
- [x] Linear-response/spectral/ASE/ART file units isolated below; TTM remains open. No general IO-error recovery, full legacy batch restart matrix or Fujitsu/multi-node execution claim.
- Records:calc2026/salmonpy-fdtd-restart-io-20261121. No commit/push or manual CI dispatch requested/performed.
- [x] Source/final installed-wheel679 tests,one optional skip. Installed SDK synthetic MPI1/2 helper tests and normal native MPI2/direct resume match source and old controls exactly. Fixed fixture-only writer/read phase synchronization; no solver collective added. Initial failed fixture/launcher evidence retained. JP/EN/AI synchronized.

### 2026-11-21 spectral output file-unit ownership

- [x] Prior validated runtime completes small LD+FDTD without host files but loses occupied unit600 at spectral finalization (STOP81); evidence retained.
- [x] Local owned newunit for LR trajectories/spectra and Fourier planes; two local handles for simultaneous ASE/ART window variants. Explicit close for periodic LR response file; no new global selector or communication.
- [x] Updated maintained inventory/profiles and audited Git-only source export plus all application objects before linking new library.
- [x] Source and installed SDK:6 cases × GNU MPI1/2, occupied600/601/602, named spectral files closed, expected outputs finite/nonempty; all control/occupied/prior-runtime and source/wheel output bytes match. Normal MPI2/direct-resume fields, point/plane/shape/checkpoint bytes also match prior runtime.
- [x] TTM ownership work cancelled by removal below. General IO-error recovery remains part of R3; no full material/convergence, multi-node or Fujitsu certification. Embedded spectral-mode restrictions unchanged.
- Records:calc2026/salmonpy-fdtd-spectrum-io-20261121. No GS collection, commit/push or manual CI dispatch.
- [x] Source/installed679 tests,one optional skip; JP/EN/AI JSON synchronized. Test-input construction failures retained and corrected without solver restrictions or new communication. TTM was queued at that checkpoint, then cancelled by removal below.

### 2026-11-21 TTM removal by project lead

- [x] Supersedes earlier pending TTM ownership work: remove unpublished TTM from maintained solver, build paths, examples, test suite and utility. Remove automatic activation and native dependency. Preserve frozen historical evidence.
- [x] Git-only clean build and source/object audits; old TTM input inert;6 spectral cases × MPI1/2 preserve caller files and previous outputs. Final installed SDK normal MPI2/direct-resume fields,point/plane/shape/checkpoint bytes match. Source/wheel679 tests,one skip. sdist/wheel exclude TTM implementation and archived source payloads; original local archives preserved. No GS collection or manual CI dispatch.
- [x] Corrected sdist force-include overriding archive exclusion; restrict runtime includes to maintained sources/manifests/docs. Current JP/EN/AI manuals and constraints synchronized. Earlier TTM-development tasks are cancelled by removal, not deferred.

### 2026-11-21 checklist reconciliation after TTM removal

- [x] Reviewed current implementation, recorded results and top queue. Fold recent R3 and file-unit work into completed scope; TTM removed, not a pending implementation.
- [x] Preserve active order5→3 after measured dense group4 completion; groups1/2 and molecules remain paused. Stable ABI, general supported initialization/IO cleanup and site evidence are still open.
- [x] Distinguish completed input-only staging/launcher construction from remaining site validation and initial-guess staging. No new calculation, GS collection, commit/push or CI in this reconciliation.

### 2026-11-21 automatic MPI factorization return

- [x] Reproduce prior native host termination on automatic7-rank assignment; no numerical/partition restriction added.
- [x] Optional helper status plus native1117 uses existing cold cleanup and status agreement; recognized unanimous/protocol1 lets Python retry. Legacy status-absent behavior and C signature retained; no new collective.
- [x] Helper MPI1/2 automatic choices match prior controls;1/2/7 status/continuation checks pass. Native MPI7 raw/public retries, child cleanup, caller survival, explicit32-step/direct-resume fields/points/checkpoint bytes match old controls.
- [x] Source/installed679 tests,one optional skip. Installed SDK helper MPI1/2/7 and native MPI7 repeat13 rejection/retry checks; fields,point/shape/checkpoint bytes match source/old controls. JP/EN/AI synchronized; TTM remains absent from sdist/wheel. Remaining R3: allocation/IO/specialized initialization; Fujitsu topology/MPI/process failure are outside this proof. R4 remains open.
- Records:calc2026/salmonpy-fdtd-auto-layout-recovery-20261121. No GS collection or manual CI dispatch.

### 2026-11-21 material-input cold rejection lifetime release

- [x] Reproduce native-cleaned1110 but Python lifetime remains consumed. Reuse the unchanged source/object-audited runtime by SHA256; no Fortran build/extra startup audit.
- [x] Recognize unanimous1110/1113 only with recovery protocol1; identify material errors, release failed Python lifetime/attachment. No numerical/input/communication changes or broader unknown-error retry.
- [x] Existing tool adds3 material modes with missing/header/truncated/nonfinite failures repeated twice; GNU MPI1/2 complete public rejection/cleanup/retry,32-step/direct-resume tests. Pre-change Python and installed SDK controls verify fields,point/material assets and direct checkpoint bytes exactly. Source/installed679 tests,one optional skip. Shared native library unchanged; JP/EN/AI synchronized.
- [ ] At that milestone coefficient returns were open;1111/1112/1114 are completed below. Remaining R3 is early-return Python classification plus unchecked supported allocations/IO and other unhandled returns; no arbitrary OS/MPI/process-failure recovery claim. R4 remains open; groups1/2 and molecules stay deferred.
- Records:calc2026/salmonpy-fdtd-material-input-recovery-20261121. No GS collection or manual CI dispatch.
- [x] Completed6 mode/rank combinations;8 rejected-input retries each in source and installed SDK. Runtime SHA256 verified, no solver rebuild or communication added. TTM remains absent from sdist/wheel.

### 2026-11-21 handled coefficient/mode cold-rejection retry

- [x] Reproduce1111/1112/1114: native cleans state/children but Python retains failed lifetime. Existing audited binary reused with recorded SHA256; no native source/binary changes.
- [x] Recognize only unanimous1111/1112/1114 with init recovery protocol1; improve diagnostics. Preserve formulas, finite checks and unsupported-mode boundaries; no new collective or silent input edits.
- [x] Existing fixture combines reader and coefficient/mode rejection; GNU MPI1/2 per3 material modes complete repeated rejection, parent/candidate cleanup, valid32-step/direct restore. Pre-change/installed SDK fields,point/material assets and checkpoint bytes match exactly. Source/installed679 tests,one optional skip; JP/EN/AI synchronized.
- [ ] Remaining R3: unhandled initialization/allocations/IO; arbitrary trapped FP/MPI/process/device failures not claimed. Stable ABI R4 remains open. Groups1/2 and molecules remain deferred.
- Records:calc2026/salmonpy-fdtd-initial-coefficient-recovery-20261121. No GS collection or manual CI dispatch.
- [x] All6 mode/rank combinations passed in source and installed SDK with10 repeated reader/coefficient rejections each. Parent MPI and resources preserved; valid retry/direct resume match old controls byte for byte. Solver SHA256 unchanged; no solver rebuild/communication added. TTM remains excluded from sdist/wheel.

### 2026-11-21 checklist reconciliation after handled coefficient returns

- [x] Re-read current native/Python error paths and latest passed summary records; mark known1111/1112/1114 recovery complete, not broad specialized-mode support.
- [x] Identify next bounded gap: native early1101–1104 returns versus Python recognized cold-return set. Keep larger allocation and output-IO recovery separate; historical STOP listings do not each become a new task.
- [x] Active sequence remains group5 R3→R4, then group3 NLCC-force/unit audit; groups1/2 and molecular work deferred, TTM deleted. Operational/site/large-scale and optional storage extensions remain separate.
- [x] Latest source/wheel evidence679 tests,one optional skip; six material/rank combinations with repeated reader/coefficient rejection and exact old controls. Recent implementations remain uncommitted. No new calculation, build, GS collection, commit/push or CI in this checklist review.

### 2026-11-21 early native grid/boundary/time cold-return retry

- [x] Reproduce clean1101 with consumed Python lifetime before fix; existing native already returns1101–1104 and cleans cold state/children.
- [x] Recognize only unanimous1101–1104 with init recovery protocol1, add actionable reasons. Preserve input/CFL rules, native binary and existing communication; no automatic repair/resubmission.
- [x] Existing grid-init/geometry fixture combines raw and public paths: GNU MPI1/2 each35 rejected-input/retry checks, parent/child cleanup, valid32-step/direct restore. Field/point/shape/checkpoint bytes match prior valid controls.
- [x] Source/installed679 tests,one optional skip; installed SDK MPI1/2 each35 rejection/retry checks and direct restore. Old/source/installed field,point/shape/checkpoint bytes exactly match. JP/EN/AI synchronized. Remaining R3: unhandled supported startup, bulk allocation/IO; arbitrary OOM/FP/MPI/process/device failures not claimed. R4 remains open; groups1/2 and molecules deferred.
- Records:calc2026/salmonpy-fdtd-early-init-recovery-20261121. No GS collection or manual CI dispatch.
- [x] Runtime SHA256 unchanged; no solver rebuild, input/CFL rule change or collective added. TTM remains excluded from sdist/wheel; recent changes still uncommitted.


## 2026-11-21 — bounded core FDTD allocation recovery

Implemented `allocate(...,stat=...)` for material-ID/scratch and split-E/H,
coefficient and optional saved-E/H 3D arrays on the native status path. All ranks
agree before entering the next communication/setup phase; normalize nonzero STAT
values to a common failure flag (including negative implementation codes). Cold
failure1118 reuses existing owned cleanup and the protocol-gated Python lifetime
release. Standalone status-absent behavior is preserved. No budget-based stop.

Validation record: `salmonpy-fdtd-core-allocation-recovery-20261121/summary.json`
in the canonical calc2026 results store. GNU/macOS MPI1/2 source and installed
wheel: two partial-allocation sites, two attempts each; parent MPI survives,
owned resource counts return to zero.32-step controls and direct checkpoint
restore (zero replay) compare byte-for-byte with previous runtime outputs.
External material CUBE MPI2 uses the same injected failure/retry/control checks.
Production and testing builds were source/object-audited; production exports no
fault hook. Source and installed suites each679 run,one optional skip,no failures.
Paired manuals and embedded AI guide text synchronized. No new GS,CI dispatch,
commit/push or confirmed cloud synchronization. This is deterministic STAT-path
injection, not a real OS-OOM campaign, Fujitsu or multi-node certification.


## 2026-11-21 — initial LD state allocation recovery

The initial LD allocation path now checks STAT for18 pole4D arrays and6
current/old-E3D arrays, retains original allocation/zero/arithmetic order and
returns a normalized failure flag. An owning-field-group agreement returns1119
before subsequent coefficients use incomplete state. The native owner releases
partially allocated components/caches/communicators; Python requires unanimous
known1119 plus the declared init-recovery protocol before retry. No new global
selectors or step collectives; memory-budget excess still warns.

Validation: `salmonpy-fdtd-ld-allocation-recovery-20261121/summary.json` in canonical
calc2026. GNU/macOS single-node MPI1/2, scalar and off-diagonal-tensor LD; both
source and installed wheels inject failures in the second4D pole/history array
and second3D FDTD-current array, twice each. Existing1118 cases are retained.
After cleanup caller MPI is alive and owner/resource counts return zero. All12
new production/testing/installed controls match the preceding runtime's fields,
material/observation data and direct checkpoint bytes;32-step direct continuation
requires zero replay. Both Python suites679 run,one optional skip,no failures.
Shared kernels are byte-identical; production/test source/object audits pass
and production contains no injection hook. The wheel enables preprocessing of
this shared kernel explicitly; no import-time compilation. Paired manuals and
AI guide text updated. No GS collection or optional CI dispatch.

Only the first LD-state allocation is covered. Material/coefficient scratch,
tensor auxiliary and subpixel-specific reallocation/PML/observation allocations,
OS kill/physical OOM/broken MPI remain outside this scope. State partial cleanup
is not a universal initialization-failure or Fujitsu/multi-node certificate.


## 2026-11-21 — material/boundary vector allocation and PML inventory

PML uses the existing split-field3D coefficient arrays already covered by1118.
Its helper has only small automatic workspaces for the existing eight layers;
no new bulk PML allocation was needed. Do not turn every automatic local array
into a recovery framework or claim OS/memory/MPI failure recovery.

Added STAT/partial-state recovery for three material coefficient vectors and
six E/H boundary vectors, returned as1120. The material gate coalesces with the
previous LD admission:allocation/initialization/arithmetic order is preserved,
1119 is retained,and only the boundary agreement is newly added when LD is
used. Without LD,both material and boundary gates are needed. No new step or
application-wide collective and no memory-budget stop.

Final record: `salmonpy-fdtd-aux-allocation-recovery-v2-20261121/summary.json` in
canonical calc2026. GNU/macOS single-node MPI1/2 source/installed scalar and
non-diagonal-tensor LD:partial material/boundary failures,owned cleanup and retry,
plus retained core/LD cases. All12 new production/testing/installed controls
match preceding runtime fields,shape/observations/sidecars and checkpoint bytes;
32-step direct continuation has zero replay. Both Python suites679 run,one
optional skip,no failures. Production/test source and object audits pass;
production exports no fault hook. Paired manuals/AI guides synchronized.

The prior uncoalesced validation and failed fixture-preparation attempt remain
preserved as records,not the final implementation. No new GS or optional CI
execution. R3 remains bounded:next inspect supported output-IO failures and
selected auxiliary buffers,then R4 and group3. No commit/push was requested.


## 2026-11-21 — optional point-observation I/O failure handling

`obsN_at_point_rt.data` headers and appends now check open/write/close status in
one private Fortran helper. Optional diagnostic errors emit a bounded warning
and continue physical evolution,not an automatic solver abort. Successful
open owns the unit;close retains any earlier write error. A subsequent append
retries after the path/header is repaired. Missing/partially written rows are
not reconstructed. Warning state is owned by the field context and omitted
from physical checkpoints. No new MPI communication or import-time behavior.

Record:`salmonpy-fdtd-point-output-recovery-20261121/summary.json` in canonical
calc2026. The prior runtime stops on a real directory collision (retained red
record). New source/installed GNU macOS MPI1/2 covers header/append directory
collisions and read-only writes (actual IOSTAT5007).18 source/installed/caller-unit
fault cases continue32 steps,release owned resources,retain caller MPI,repair
future output and directly restore with zero replay. Normal scalar/tensor point
output and all physical checkpoints match previous bytes. Failed-output cases
match normal fields/shape/info/checkpoints;post-repair point rows match the exact
expected suffix with deliberately missing rows. Caller file lifetime is verified
through initialization,steps,direct restore and finalization on MPI2.

Both suites679 run,one optional skip,no failures. Native source/object audits,
Git-only source preparation/build,and sdist→wheel installation passed. Paired
manuals and AI guide text synchronized. No new GS or optional CI dispatch.
Close errors are checked,but a close-only failure was not injected. This is not
protection for plane-integral/info/restart output,arbitrary filesystem failures
or automatic reconstruction of diagnostics. Continue bounded R3 work,then R4
and group3;do not overgeneralize optional-output warnings to required restarts.


## 2026-11-21 — optional plane-integral output I/O

Plane-integral headers/appends now check IOSTAT in a private Fortran writer;
optional output failure warns rather than aborting evolution. The field context
owns the warning flag;physical checkpoint layout/numbers are unchanged. Existing
integration and MPI summations remain in their original order,including all
subsequent planes after root write failure. No new MPI communication. Successful
opens own their units;failed writes still close them and preserve the first error.
The public Fortran routine keeps existing source calls valid with an optional
trailing warning flag;no new process-global selector.

Record:`salmonpy-fdtd-integral-output-recovery-20261121/summary.json` in canonical
calc2026. Old runtime stops on a real directory collision (retained red record).
New source/installed GNU macOS MPI1/2 cover18 header/append directory/read-only/
caller-unit fault cases. Read-only writes return actual IOSTAT5007. Filesystem
repair resumes later rows;deliberately missing rows are not fabricated. Normal
scalar/tensor integral/point outputs and physical checkpoints match old bytes.
Failed xy output leaves yz/xz files and physical fields/checkpoints unchanged;
32-step zero-replay direct continuation and caller-file/MPI lifetimes pass.

Both suites679 run,one optional skip,no failures. Source/object prohibition audits,
tracked-only source preparation/build and sdist→wheel installation pass. Paired
manuals/AI guides synchronized. No GS collection or optional CI dispatch.
Close errors are handled but close-only failure was not injected. Plane snapshots,
frequency-space maps,info/restart outputs and arbitrary filesystem-wide failures
remain outside this bounded scope. Continue bounded R3,then R4 and group3.


## 2026-11-21 — optional final observation information

The obs0_info.data animation helper file now checks open/write/close IOSTAT.
A failed optional metadata write emits SALMONPY_OBSERVATION_INFO_WARNING and
continues normal field/state/MPI finalization. Preserve the first error,close
only a successfully opened owned unit,and do not delete earlier results or
silently reinterpret failed required restart saves. No new collective,global
state,checkpoint layout or numerical calculation. No automatic retry after
finalization;check missing/incomplete metadata before visualization/analysis.

Record:`salmonpy-fdtd-info-output-recovery-20261121/summary.json` in canonical
calc2026. The corrected red fixture reproduces a prior-runtime Fortran directory
collision at finalization;the initial fixture-listing failure is also preserved.
New source/installed GNU macOS MPI1/2 cover10 real directory/read-only/caller-file
fault cases (write IOSTAT5007),completed32-step trajectories and zero-replay direct
continuation. Normal scalar/tensor info/output bytes match prior controls.
Faulted metadata leaves time series,integrals,shape and physical checkpoints
unchanged;state/MPI counts return zero and caller files/parent MPI survive.

Both suites679 run,one optional skip,no failures. Final native source/object
audits,Git-only preparation/build and sdist→wheel install pass. Added guards
fit standard132-column free-form lines without relaxing compiler limits;initial
build evidence remains preserved. Paired manuals/AI guides synchronized.
No GS collection or optional CI dispatch. Close-only failure was not injected.
Remaining R3 work:snapshots/frequency maps,other optional metadata,required
restart I/O and selected auxiliary buffers;then R4 and group3 in the active order.


## 2026-11-22 — optional plane-snapshot I/O

Time-resolved E/H plane snapshots now check open/write/close IOSTAT and warn
instead of aborting valid physical evolution. Keep original plane aggregation,
unit conversion,component order and MPI sums;failed root writes never skip a
later collective. Workspaces are released as before. Warning state is field-owned
and excluded from physical checkpoint layout;no new global flag or collective.
The optional final warning argument preserves existing Fortran source calls.

Record:`salmonpy-fdtd-snapshot-output-recovery-20261122/summary.json` in canonical
calc2026. The old runtime stops at a real snapshot-directory collision. New
source/installed GNU macOS MPI1/2:10 directory/read-only/caller-file fault cases
continue,repair the destination and perform zero-replay direct restoration.
The validation then explicitly advances17/18 again;that continuation writes
frame18 successfully and matches the previously obtained physical trajectory.
All90 sampled frame files and other outputs match normal controls. Normal
scalar/non-diagonal-tensor outputs/checkpoints match prior runtime bytes.
Caller-file/state/parent-MPI lifetimes pass. Observation stride6 limits test
file count;this is not a new default or output-size restriction.

Both suites679 run,one optional skip,no failures. Native source/object audits,
Git-only preparation/build and sdist→wheel install pass. Paired manuals and AI
text synchronized. No GS collection or optional CI dispatch. Read-only faults
were measured on a non-root Mac;close-only failure was not injected.
Missing/partial frames are not automatically reconstructed. Staging allocation,
frequency maps,required restart I/O and arbitrary filesystem failures remain
separate scopes. Continue bounded R3 work,then R4 and group3 as already ordered.


## 2026-11-22 — batch frequency-resolved plane-map I/O

Added IOSTAT checks to frequency-map header/row/close writes. Root failure warns
once per finalization and continues every remaining window/plane/energy and
existing Fourier aggregation. Close only successfully opened owned units and
retain the first error. No new global selector,collective or physical numerical
operation. Missing/partial maps are not reconstructed or retried after finalize.

Scope is the SALMON batch FDTD path. The native E/H API already rejects these
accumulators with15;those guards are intentionally unchanged. This work does not
certify native protected stepping or direct spectral-state save/restore.

Record:`salmonpy-fdtd-frequency-output-recovery-20261122/summary.json` in canonical
calc2026. Prior runtime stops on a real map-file directory collision. New source/
installed GNU macOS MPI1/2:two energies,three planes,two windows;all12 normal maps
and other outputs match prior runtime bytes. Eight directory/read-only fault
cases preserve other11 maps,point time series and info output. Actual read-only
write IOSTAT5007;caller units600/601/602 retain before/after content and all tested
map units are closed before returning to the host. Host/MPI programs finish.
Existing consolidated spectrum probe is reused;focused mode suppresses redundant
time snapshots and limits the trajectory to64 steps. No new Python test module.

Both suites679 run,one optional skip,no failures. Final runtime and new host
probe source/object audits,Git-only preparation/build and sdist→wheel install
pass. Standard132-column free-form source retained. Early environment-import
failure,old-runtime red run and an unused pre-edit build remain preserved as
evidence;only the final source-verified runtime is used for positive results.
Paired manuals/AI text synchronized. No GS collection or optional CI dispatch.
Close-only failure was not injected;non-root Mac permission semantics measured.
LR/ASE/ART output,required restarts,staging allocation and arbitrary filesystem
failure remain separate scopes. Continue bounded R3,then R4 and group3.


## 2026-11-22 — batch LR/ASE/ART optional spectral I/O

Added shared private open/close helpers within eh_finalize,checking IOSTAT and
closing only owned opened units while preserving the first error. LR trajectory
and response writes are independently guarded;failed raw I/O never skips Fourier
processing/subsequent tables. ASE/ART windows keep separate units,status and names
through their energy loops,so one failed file does not contaminate the other.
Physical formulas,normal formatting and every existing integration/MPI sum remain
unchanged. A local per-finalization warning flag bounds messages;no new global
selector,collective or numerical recomputation. No public/native C ABI expansion.

Record:`salmonpy-fdtd-optional-spectra-recovery-20261122/summary.json` in canonical
calc2026. The old runtime stops on a real ART directory collision. New source/
installed GNU macOS MPI1/2 cover isolated LR,periodic LR,non-unit-background raw
LR,ASE and ART:20 normal-configuration comparisons match prior byte payloads;
40 directory/read-only fault cases retain all other output bytes,caller files
and spectral-unit cleanup. Actual read-only write IOSTAT5007. Programs complete
MPI normally. The non-unit-background legacy normalization warning/omission is
preserved. Existing consolidated host probe reused;no new Python test module.

Both suites679 run,one optional skip,no failures. Runtime/host source and object
audits,Git-only build and sdist→wheel install pass. Standard132-column source
retained. Paired manuals/AI text synchronized. No GS or optional CI dispatch.
Native E/H diagnostic rejection15 is unchanged;native spectral protected stepping/
restart is not certified by this batch test. Missing/partial spectra are not
fabricated/retried;close-only failure was not injected. Non-root Mac permissions
measured. Required restarts,initialization/helper I/O and selected auxiliary
allocation remain bounded R3 work,then R4 and group3 in the agreed order.


## 2026-11-22 — current common C/Fortran host descriptor

- [x] Previous verified tensor step workspace changes committed/pushed as `f4841bcd`.
- [x] Added local `sp_host_abi` version/type-width/communicator/ownership/thread contract query. No MPI initialization, solver initialization or per-step collective. Missing/mismatched descriptors reject before attachment; no old-version fallback. Single/shared/registry are current public route names.
- [x] Ordinary and fast-PR HSE rebuilt with mandatory source/object audits. C/Fortran probes verify short-buffer rejection without writes, unchanged tails, equal descriptors and MPI remaining uninitialized. These are contract probes, not HSE physical trajectories.
- [x] Local MPI6 split into groups1/2/3: rank-local missing/version/width/status errors reject collectively before attachment; subsequent normal attachment and custom A/J response succeed.
- [x] Source suite679 tests,one optional skip. LD/FDTD scalar,diagonal and rotated-tensor MPI1/2 normal/direct-resume outputs match prior controls byte-for-byte. Results: canonical calc2026 `salmonpy-host-abi-contract-20261122/comparisons.json`. No GS collection.
- [x] Japanese/English API/runtime guides and AI JSON synchronized; manual synchronization and whitespace checks pass.
- [ ] R4 remaining: operation-specific signature/lifetime guarantees and wider profile/platform evidence. The common descriptor is not blanket stable solver ABI certification. New changes remain uncommitted; no optional CI dispatched.

Continue bounded R4 operation contracts,then group3 NLCC runtime certification/unit audit. Groups1/2 and molecular work remain deferred. Multinode/Fujitsu/GPU and full HSE physical validation are not established by these local tests.


## 2026-11-22 — common host operation contracts

- [x] Common descriptor changes committed/pushed as `26dca23c` to independent SALMONpy.
- [x] One operation declaration supplies exact ctypes signatures and detached `native_contract(...)["host_operations"]` documentation for attach/detach/info/barrier/sum/library. Resolve all six symbols before executing any operation; missing functions leave binding declarations untouched. Custom backends remain unconstrained.
- [x] Public C header now declares all six common operations and states communicator/buffer ownership, collective/local scope, main-thread ordering, lifetime and returned-error limitations. No Fortran arithmetic, solver initialization, timestep reduction or checkpoint payload changed.
- [x] Local source/installed MPI6 groups1/2/3 reject one-rank missing operation symbols before attachment, then attach/reattach and compute manufactured independent A/J responses successfully.
- [x] Audited C host against ordinary and fast-PR HSE runtime: MPI1/2/3 range rejection, duplicate attach rejection, rank/size, sum, barrier, truncated MPI text, detach/reattach and surviving host MPI pass. Query-only C/Fortran probes pass without MPI initialization. No GS data collected.
- [x] Japanese/English API and AI JSON synchronized. Source679 tests,one optional skip,pass. Wheel-based normal/direct-resume LD/FDTD scalar/diagonal/tensor MPI1/2 match prior bytes.
- [ ] Solver-specific operation signatures, all state/transaction error guarantees and broader platforms are not certified by the common operation suite; continue bounded R4 then group3 NLCC certification/unit audit. No arbitrary MPI/process failure recovery claim.

Final installed validation details are recorded in canonical calc2026
`salmonpy-host-operations-20261122/summary.json`. These operation-contract changes
are uncommitted; groups1/2 and molecular work remain deferred. No heavy CI dispatch.


## 2026-11-22 — TDDFT / fast-PR HSE state-operation contracts

- [x] Common host operation changes committed/pushed as `a2411deb` to independent SALMONpy.
- [x] Current single/shared/registry required state signatures now have one declaration per component. `native_contract(...)["state_operations"]` exposes ordered parameter names, exact type tokens, in/out direction, ownership, collective scope, endpoint/current conventions and checkpoint/transaction modes. Python factories bind those same declarations; no import-time loading or new timestep communication. All required symbols resolve before binding; custom backends are unconstrained.
- [x] Preserve component differences: ordinary single C returns matter J, handles return charge j; HSE also returns total energy. TDDFT inherits model dt/capacity/decomposition while HSE spawn receives dt/capacity. Transaction mode2 is standard release vs HSE reopen/flush/release with retention on failure. Raw checkpoint load is not whole-session atomicity.
- [x] Existing consolidated signature tests cover all missing required entries without executing state operations or leaving partial type declarations. No extra Python test module/count. Japanese/English guides, AI JSON and C header ownership/error notes synchronized.
- [x] Existing C/Fortran state host reused with optional `--abi-state --state-checkpoint`. Ordinary PZ and fast-HSE MPI1/2 C/Fortran currents agree; independent points, rejected old A, stale generations and caller MPI lifetime remain correct. C loads direct checkpoints into new states and compares subsequent currents exactly. GS is referenced by symlink and byte-checked unchanged; no GS calculation or collection.
- [x] Source vs prior installed Python, ordinary/HSE × shared/registry × MPI1/2:8 profiles,16 executions, two independent drives each; currents match bitwise and direct resume after step6 matches continuous steps7–12, no replay. Result directory: canonical calc2026 `salmonpy-state-operation-contracts-20261122`.
- [ ] Broader XC/symmetry/platform evidence, every optional diagnostic/GS import and arbitrary MPI/process/invalid-pointer failure recovery are separate; do not label the complete solver ABI stable from this scoped verification. Continue bounded R4 then group3 NLCC certification/unit audit. Groups1/2 and molecules remain deferred.

Final source/distribution checks and installed matrix are recorded in that
directory's summary.json. This state-contract work remains uncommitted; no
optional heavy CI dispatched.


## 2026-11-22 — FDTD operation contract and next queue

- [x] TDDFT/fast-HSE state contracts committed/pushed as `ab987567`.
- [x] Added current `field_operations` report with ordered names/types/directions, C/Fortran layout, current sign/time normalization, local vs collective scope and checkpoint/transaction modes. New `native_fdtd.h` ships with sdist/wheel. Python creation/update/routing/H/plane reads and protected-step/restore entry bindings use the same declaration. Optional features stay checked at use; custom backends are unchanged.
- [x] Explicit distinctions: one-based int32 grid triples vs Python zero-based; C (points,3) float64 equals Fortran (3,points); J parameter is physical charge j, no cell-volume multiplier. current_mask/plane_read are local; H sampling and stepping/save/restore are collective. Field transaction mode3 is no-op, not electronic stream reopening.
- [x] Source and installed679 tests each pass,one optional skip. Both source and wheel run scalar/diagonal/off-diagonal LD at MPI1/2 with normal/direct-midpoint resume; all compared fields/output/checkpoint payloads match prior controls byte-for-byte. No solver arithmetic/update order/state format changed; no extra MPI collectives.
- [x] New header compiled with -Wall/-Werror; typed references to all current operations link against the hash-verified audited runtime. Protocol query host runs without field/MPI initialization. This link/query test does not certify every optional observation profile.
- [x] Paired API manuals/AI JSON synchronized; final wheel bytes match current source. First matrix's tensor setup failures came from a missing test-only LD kernel link before setup, not the field update; logs preserved. Prepared link and fresh v2 source/wheel matrix pass. No GS collection, no heavy CI dispatch.

Record:canonical calc2026 `salmonpy-field-operation-contracts-v2-20261122/summary.json`.
Changes here are uncommitted. Next bounded R4 task is GS operation declarations/
bindings,then group3 maintained NLCC certification and unit audit. Broader
platform/site/optional-profile proof remains an operational queue; do not extend
this into every small allocation or add restrictions on usable/custom inputs.
Groups1/2 and molecular work remain deferred; TTM remains deleted.


## 2026-11-22 — GS operation contracts

- [x] Current GS create/import/info/array/close, active SCF, numerical snapshots, candidate, disk bootstrap/read and explicit repartition declarations organized in `native_contract("gs",route="single")["gs_operations"]`. Python core/SCF/disk/candidate/repartition/import/format bindings share those declarations. Single names the interface, not a one-state restriction. Preserve optional capability queries and custom backends.
- [x] GS info/read are collective because gs_lookup agrees tokens; resource counts and pure protocol/format queries are local. Named argument types/directions, Fortran array extents, snapshot ownership and candidate/disk modes documented in Japanese/English/AI guides. No Fortran arithmetic or format changes.
- [x] Source679 tests pass,one optional skip. C header compiles/links with audited host object and unchanged verified runtime; protocol/format query runs without MPI/GS initialization.
- [x] Source PZ/fast-PR HSE internal MPI1/2 comparison, mixing-history snapshot restore, sibling preservation and protected candidate failure/retry pass; MPI2 scalar result reports equal prior installed frontend. No wavefunction/density/restart files collected. This is an8-iteration regression, not physical convergence.
- [x] Existing regression tool now has explicit `--allow-unconverged-fixture`; default production residual-review pause remains. Initial nonconvergence pauses retained as evidence.
- [ ] Installed-wheel SCF numerical matrix awaits explicit permission: auto-review rejected it under the current GS-collection restriction. Actual test trees contain standard GS scalar outputs (eigenvalues/forces/residual/logs), even though wavefunction/density/restart writing is disabled. Do not describe these tests as producing no GS-related output. Preserve existing evidence; do not bypass rejection.

Installed unit/query/package evidence is recorded separately from the blocked
SCF matrix. GS contract changes and the preceding FDTD contract remain
uncommitted. Once the authorized installed matrix finishes, proceed to bounded
NLCC certification/unit audit; broader platform/site proof is separate.


### GS wheel validation completed after explicit output permission

2026-11-22:user authorized eigenvalue/force/residual/log saving for this wheel
validation. The previously blocked installed matrix now passes PZ and fast-PR
HSE,internal MPI2 ×3 groups. Source/installed reports agree exactly;numerical
arrays,mixing history and snapshot restore equal continuous execution. Late
restore failure preserves originals and retry matches;restored iterations0.
File inventory contains prepared input/PP and authorized standard outputs/reports.
Wavefunction/density/GS restart files absent;created restart directories are empty.
This closes the preceding pending installed matrix. Preserve initial rejection
as historical evidence;no broader future GS-saving permission is inferred.
Summary/inventory:canonical calc2026 `salmonpy-gs-operation-contracts-20261122`.
This short deliberately unconverged regression is not physical-convergence
certification. FDTD/GS contracts remain uncommitted. Next:bounded NLCC/unit audit.


## 2026-11-22 — maintained NLCC certification and unit audit

- [x] Previous FDTD/GS contracts committed/pushed as `a2f0a9a4`.
- [x] Current audited fast-PR source/executable:prescribed XC variation PZ/Libxc-PZ/PBE with/without core at OMP1/2,core-force6 geometries,and MPI2/4 ×OMP1/2 both core-table formats pass. Initial unsupported MPI1 selection preserved;same hash-verified helper rerun on prescribed2/4 ranks. No prohibited driver executed.
- [x] Fixed Si8/PBE whole-application comparison authorized explicitly by user for this NLCC campaign:energy error0,maximum force error9.772998e-16 Hartree/bohr. Only input/energy-force-residual reports/log retained;wavefunction/density/restart datasets absent. This is implementation regression,not physical convergence.
- [x] Fixed evidence reader bugs:core parser had silently ignored geometries4–6;MPI parser expected one success line instead of both current table formats. Tests reject a bad sixth geometry,missing formats,nonfinite errors and noncollective rejection.
- [x] Updated certification to the actually validated maintained source tree. Require tree identity in build/receipt/native evidence;old/missing identity rejects. MPI-linked certified build additionally requires its own complete MPI-force proof. Current CMake source binding used;no assumed old serial-source directory. Completed binary/reference-bound application record may be reused explicitly,without SCF rerun. No automatic certification from a feature name.
- [x] `unit_contract()` documents unchanged solver input27.21138505 eV/Hartree,Python energy output27.211386245988,geometry/ASE0.52917721067 angstrom/bohr. Relative energy convention difference4.3951750e-8;separate legacy+U27.2116 noted outside declared profile. No conversion constants,numerical formulas,reference data or input interpretation changed.
- [x] Japanese/English API and AI JSON synchronized. Source679 tests pass(one optional skip);validation-tool7 tests pass. Final installed unit/capability checks recorded separately below/in summary.

Record:canonical calc2026 `salmonpy-nlcc-current-audit-20261122`. Validated runtime
and receipt copied byte-identically under validated-runtime for reuse;original
build/test hashes retained. NLCC/unit changes uncommitted. Do not reopen every
small allocator/optional operation as a blocker. Broader site-MPI/multinode,other
platforms,scale/observation I/O,initial-guess staging,Ag fit and optional storage/
GPU remain separate queue. Groups1/2 and molecules remain deferred;TTM deleted.


NLCC/unit final installed verification:source/wheel679 each pass(one optional
skip),evidence-reader7 tests pass. Final installed unit_contract AI search,
byte-identical saved runtime+receipt loading and real-NLCC-UPF ASE forces
admission pass without another SCF run. Source/installed API/runtime/capability/
reference/manual files match. Current fixed Si SCF converged at72 iterations,
residual7.60746e-12<threshold1e-11;this is not physical discretization convergence.
Summary.json records exact SHA256 and limited CPU/XC scope. No new optional CI
dispatched. Changes remain uncommitted.


## 2026-11-22 — explicit warm-start input preparation

- [x] NLCC/unit work committed/pushed as `5ee9b6c8`.
- [x] `prepare_dft(...,initial_guess=GroundStateResult)` stages only info.bin,occupation.bin,wfn.bin. Shared DFT.run compatibility check and namelist renderer preserve existing run input bytes. Source remains unchanged;no automatic discovery,SCF/runtime load,compile or submission. Required executable SHA256 and source manifest digest are recorded for the separate compute-node job.
- [x] Stage verifies copied seed hashes before publication; incompatible models or changed copied bytes leave no final prepared directory. Existing unconverged-SCF review and wavefunction-only semantics remain;this is not mixing-history restart. Fresh preparation is unchanged.
- [x] Consolidated existing DFT/preparation tests use synthetic files,compare old run/prepared input bytes,and inject a changed staged wavefunction. No real GS copied,collected or computed. Sample `prepare_electrons` has an explicit --initial-guess flag;paired preparation guides,DFT guide,AI guidance and searchable task/API index synchronized.
- [ ] Remaining initial-guess generalization:retained GS/fast-HSE live-state and other archive/custom-backend seed handoff. Do not infer changed-grid/kmesh interpolation or HSE support from this batch-file staging step.

Source/wheel suite and package details are recorded in canonical calc2026
`salmonpy-guess-preparation-20261122/summary.json`. Current changes uncommitted.
Site/other-platform/scale validation,Ag fit,optional formats/GPU remain separate;
groups1/2 and molecular work deferred. No optional heavy CI dispatch.


## 2026-11-22 — retained GS / fast-PR HSE initial guess

- [x] Preparation-only work committed/pushed as `4d6c2b59`.
- [x] Common GS C contract and Python `initial_guess` copy only wavefunctions/occupations from a live owner. Target density/potentials, mixing and HSE resources start independently; temporary seed storage is released after initialization.
- [x] Same-session/model/layout checks permit atomic displacement; source and siblings remain unchanged on rejection, target continuation, and source close. This is not exact SCF restart or changed-grid interpolation.
- [x] Git-only export: all728 solver inventory entries tracked and matching; source/object eigensolver audits and fast-HSE build passed.
- [x] Final source/wheel679 tests each passed (one optional skip each). Final PZ/fast-HSE k2 MPI6 and Gamma PZ MPI3/fast-HSE spatial MPI6 passed; wheel PZ/HSE MPI6 passed. Normal cold-start arrays match prior runtime exactly; source/siblings and receiver-after-source-close checks pass. C typed declaration/link/query passes. Evidence: canonical calc2026/salmonpy-gs-direct-guess-20261122/summary.json. No GS datasets collected; finite developer fixture is not physical convergence. Current direct-handoff changes remain uncommitted.
- [ ] Archive/custom-backend initial-guess handoff remains separate. New-binary NLCC recertification, site/multinode/other-platform/scale work and optional formats/GPU remain separate. Groups1/2 and molecules remain deferred.


## 2026-11-22 — NLCC evidence for retained-GS initial guesses

- [x] Current source/binary NLCC consistency: XC variation maximum3.061e-10, core-force variation5.951e-8;2/4-rank×1/2-thread replication maximum5.552e-17. Source/object eigensolver audits passed before linking each native test.
- [x] Fixed Si PBE whole-application comparison: energy difference0;maximum force difference9.773e-16 Hartree/bohr. Only input/scalar reports retained; no wavefunction/density/restart data. Native records and binary-bound receipt saved under canonical calc2026/salmonpy-gs-guess-nlcc-20261122; validated-runtime holds identical executable bytes for explicit reuse.
- [x] Advance the NLCC source identity only after these checks; existing source-inventory test now detects a stale capability identity. Identity mismatch reproduced before update; focused source/receipt tests pass afterward. Numerical formulas/constants/reference values unchanged. Daily CI checks identity/integrity only; no material calculation added.
- [x] Final source/wheel679 tests each pass(one optional skip each). Installed wheel loads the byte-identical saved executable and current receipt; real NLCC UPF ASE force admission passes without another SCF. Evidence:canonical calc2026/salmonpy-gs-guess-nlcc-20261122/summary.json. Fixed Si SCF converged at72 iterations,residual7.60746e-12<threshold1e-11. No optional CI dispatch. Current changes remain uncommitted.
- [ ] Archive/custom-backend initial guesses and broader site/platform/scale/GPU evidence remain separate; groups1/2 and molecules deferred.


## 2026-11-22 — common archive/custom-reader initial guesses

- [x] Previous retained-GS/NLCC work committed/pushed as8f3a7663;approved direct-array design recorded in67917b56(local).
- [x] GSInitialGuess.from_arrays/from_archive/from_provider shares one synchronous rank-local C/Fortran reader contract. Explicit DFT source model and GS/2 conventions; immutable metadata/array copies, no inheritance or intermediate GS record files. Native bounds avoid guessed MPI mapping and mandatory full-state gather/broadcast. Existing setup agreement includes mode/metadata hashes, with no additional per-step communication.
- [x] Fortran owns numerical copy/finite/kpoint/weight checks and fresh SCF density/potential/mixing/HSE initialization. No complex truncation for real Gamma. Existing retained-state route and cold results remain identical. Reader exceptions, NaN/Inf, sampling mismatch and MPI mode disagreement preserve live owners; retry matches controls. Allocated communicator-slot capacity is not confused with live-resource ownership.
- [x] Source and installed-wheel680 tests each pass(one optional skip each);typed C declaration/link/query, tracked-only728-file inventory build and forbidden-eigensolver audits pass. PZ/HSE k2 MPI6, Gamma PZ MPI3 and Gamma HSE spatial MPI6 pass. Synthetic GS/2 route matches direct arrays; archive deletion/provider-buffer mutation after initialization do not change continuation. No real wavefunction/density/restart datasets saved. Wheel scalar outputs explicitly authorized; source/whole-app NLCC outputs deleted after temporary comparison.
- [x] New-binary NLCC local/MPI/fixed Si reference checks pass;energy difference0,max force difference9.773e-16 Hartree/bohr. Capability source identity and binary-bound receipt refreshed only after validation. Paired guides and AI task/API catalog synchronized; new module packaged in wheel. Evidence:canonical calc2026/salmonpy-array-guess-20261122/summary.json. Implementation changes remain uncommitted;no optional heavy CI dispatched.
- [ ] Remaining broader work: site/multinode/other-platform/scale evidence, optional formats/GPU and Ag fitting; no implicit grid/k interpolation or arbitrary solver-profile certification. Groups1/2 and molecules remain deferred.


## 2026-11-22 — publication and queue reconciliation

- [x] Retained-GS/NLCC implementation8f3a7663 already pushed to the independent SALMONpy branch.
- [x] Approved design67917b56 and verified common archive/custom-reader implementation are included in this requested publication. Current queue at the top supersedes historical unchecked entries about initial guesses and pending selected R3/R4/CPU-NLCC integration.
- [x] Final pre-commit checks: source inventory and receipt binding, initial-guess/public API/manual tests, guide synchronization and whitespace checks. Previously completed source/wheel680-test and MPI numerical records remain the numerical evidence; no redundant physical campaign or optional heavy CI dispatch.


## 2026-11-22 — molecular work resumed

The project lead explicitly resumed molecular work and selected two-way FDTD
coupling first. This supersedes the historical molecular pause, not other deferred
groups. Approved scope: single fixed molecule, length-gauge E drive; verify a
Fortran linear oscillator with centered time/work and own-vacuum subtraction
before attaching isolated TDDFT. Plan:2026-11-22-molecular-centered-feedback.md.
No real GS data collection is authorized or needed for the oscillator stage.


## 2026-11-22 — molecular centered-feedback foundation completed

- [x] Public dipole→P/j mapping and own-source subtraction use the bundled Fortran kernel; existing API/atomic units, conservative weights, read-only ownership and regression values preserved. No import-time build/MPI or old numerical fallback. Direct IEEE NaN/Inf/overflow rejection preserves inputs; mapping tests cover conservation, superposition and static offsets.
- [x] Fortran development midpoint oscillator, accepted-state work and normalized summaries pass independent analytic second-order and reversed-work-sign controls. This is a test substitute, not molecular TDDFT.
- [x] Two independent fields use trial/rollback/corrected-current updates after one explicit hold. Vacuum own-source subtraction and synthetic LD scattering pass. dt=.003/.0015/.00075 at fixed duration satisfy normalized work tests; dt=.003 temporal mismatch falls from held3.598e-11 to2.67e-27 au. No full Maxwell/flux conservation or overall coupled-order claim.
- [x] Ordinary FDTD with accepted j reproduces accepted E without trial/rollback. Public one-way mapping passes1/2 field ranks (2/3 total); centered fixture uses4 total ranks with one rank per field. Source/wheel E,p,j,statistics and one-way histories match exactly.
- [x] Source/wheel681 tests pass(one optional skip each), standalone oscillator2 tests and tracked-only lightweight build pass. Source/object prohibition guards remain; unchanged field library reused by recorded SHA256. Paired guides, AI instructions/JSON and sample description synchronized. Artifacts/code patch saved in canonical calc2026/salmonpy-molecular-centered-20261122. No real GS data or optional heavy CI dispatch. Implementation uncommitted; planbd45e553 local.
- [ ] Next: connect isolated length-gauge E→charge-dipole TDDFT to this trial/accept schedule; verify quantum-state preservation/retry and endpoint-E/clock alignment against the current worker/runtime. Real-GS collection still needs applicable permission.
- [ ] Wider physical/media/multi-molecule/platform/GPU/multinode validation remains separate; no automatic production refinement campaign. Molecular work is active; groups1/2 remain deferred.

Evidence:canonical calc2026/salmonpy-molecular-centered-20261122/{summary.json,validation.md}.

## 2026-11-22 — isolated E-worker / two-field connection

- [x] Private single-molecule transaction passes predicted endpoint E to the existing
  isolated worker, maps charge-dipole differences in Fortran and validates actual
  selected endpoint E before committing. No change to the solver's update order.
- [x] Initial field hold is explicit; electronic time=field time-dt, initially zero
  selected E required. No second Python midpoint average. Existing isolated and
  periodic public interfaces are unchanged.
- [x] MPI4 native-field + explicitly synthetic E-worker cases pass in vacuum and LD.
  Independent replay of the supplied E history reproduces dipoles exactly. Injected
  endpoint rejection/retry gives identical full saved E/p/j/drive/time histories to
  continuous runs. These are transport tests, not molecular TDDFT physics evidence.
- [x] Fix discovered in rejection test: field transaction rollback is single-use;
  retire predictor token and capture a fresh accepted-attempt token. Electronic
  cancellation uses existing accepted-history replay, not direct/no-recompute restore.
- [x] Source Python suite682 tests passes(one optional skip); paired molecular
  guides and AI catalog synchronized. No new real GS collection or optional CI.
- [ ] Real isolated TDDFT/FDTD validation needs an existing GS bundle plus validated
  isolated E-worker binary/receipt. None identified in this phase; requested existing
  paths. Quantum norm/work, full-field energy and a public reusable two-way driver
  remain open. Groups1/2 remain deferred.


## 2026-11-22 — authorized real C2H2 / FDTD connection

- [x] User authorized this C2H2 campaign's GS/wavefunction/density/restart and worker
  copies. GS generated once: PZ, fixed ions,24×24×32 at0.35Å,123 iterations,
  residual1.0780440e-10. Saved under canonical calc2026/salmonpy-c2h2-feedback-20261122.
- [x] Reused audited binary7f8e0c7f… after SHA256 check. Two batch RT references
  match exactly. Isolated E1 worker/reference maximum dipole error3.8553e-12 e·bohr;
  independent workers and accepted-history replay match exactly.
- [x] Actual C2H2 connected to two-field MPI4 vacuum/LD fixture;64 steps at.003au.
  Independent replay of the supplied E drive matches dipoles exactly. Injected
  rejection at interval30 rolls fields and electrons back, retries and matches all
  continuous E/p/j/drive/time arrays exactly. Electronic rollback replays history.
- [x] Fixed a real-MPI startup defect: independent single-process workers inherited
  parent MPI connection/rank environment and hung. Filter these variables in the
  child only; preserve parent MPI/environment and scheduler allocation IDs. Unit
  test reproduces before/fixes after. macOS/OpenMPI evidence only.
- [x] Source suite683 tests passes(one optional skip); paired guides/AI JSON updated.
  No optional heavy CI. Source and package changes are not committed yet.
- [ ] Remaining molecular scope: measured norm/quantum work/full electromagnetic
  energy, physical convergence, public reusable two-way workflow, broader platforms.
  The previously missing real-GS/E-worker inputs are now available for this campaign.
  Groups1/2 remain deferred. No approval implied for unrelated GS campaigns.

- [x] Installed wheel683 tests pass(one optional skip); real vacuum/LD retry histories exactly equal source continuous runs. Evidence:docs/isolated-fdtd-validation-2026-11-22.md and canonical campaign records.

## 2026-11-22 — isolated measured electron number

- [x] Previous molecular feedback committed/pushed as b7657d4a. No manual heavy CI.
- [x] E2 emits Fortran's already integrated occupied electron count with six dipole
  values; no extra global reduction. Public scalar and frame observation, seven-value
  replay snapshots, revision2 evidence/checkpoint schema. No old-format fallback.
- [x] Authorized C2H2 GS reused: actual E2/reference dipole agreement unchanged;
  memory/disk replay and independent states preserve measured electron number.
  MPI4 LD/FDTD continuous vs reject30/restore/retry: all arrays bit-identical,
  including electron count. Earlier E/p/j/time arrays unchanged. Relative electron
  loss9.9476e-15<1e-5 diagnostic target. Do not stop solely for exceeding target.
- [x] Git-only source build and source/object guards pass. New-source NLCC local XC
  and core-force derivatives MPI1/2×OMP1/2 pass before source capability pin update;
  no new GS or whole-application NLCC force certificate. Source/wheel684 tests each
  pass(one optional skip). Paired guides/AI JSON updated. Follow-up uncommitted.
- [ ] Molecular remaining: quantum work/full electromagnetic energy, physical
  convergence and general reusable two-way workflow; broader sites/platforms.
  Electronic restoration is still accepted-history replay. Groups1/2 deferred.

Evidence: docs/isolated-electron-number-validation-2026-11-22.md and canonical
C2H2 campaign's electron-number-summary.json (local OneDrive, cloud sync unconfirmed).

## 2026-11-22 — molecular energy/work diagnostic

- [x] Electron-number work committed/pushed as416b1a83. No manual heavy CI.
- [x] E3 exposes internal energy and freshness without changing existing energy
  cadence, physical arithmetic or adding a reduction. None at stale endpoints.
  Revision3 evidence/replay checkpoints; no old-format fallback.
- [x] Fortran held-endpoint and midpoint dipole work, sparse energy differences,
  independent zero-field drift and finite/overflow rejection. Analytic sign/time
  controls pass. Actual C2H2 reference energy difference<=2.132e-14Hartree.
- [x] Authorized existing GS reused in MPI4 LD/FDTD. Old E/p/j/Ne arrays unchanged;
  source continuous/wheel injected reject30/retry histories exactly match, including
  work/energy. Source/wheel686 tests(one optional skip each); memory/disk replay
  and freshness/None are checked with actual C2H2. Work-observation follow-up uncommitted.
- [ ] Energy agreement is NOT closed: default weak signal is below zero-field drift.
  Explicit amplitude.1 control has0.897% held-work residual, but28.49% midpoint-work
  residual. Next resolve molecular/field time placement, then full Maxwell/LD/flux
  energy balance. No automatic background subtraction or residual-based stopping.
- [ ] Physical convergence/general public two-way driver and broader site/platform
  evidence remain separate. Electronic restore still replays accepted history.
  Groups1/2 deferred. No new GS campaign authorized by this diagnostic.

Evidence:docs/isolated-work-validation-2026-11-22.md and canonical C2H2 campaign's
energy-work-summary.json. Local OneDrive storage; cloud synchronization unconfirmed.

## 2026-11-22 — selected-field midpoint work alignment

- [x] Fortran selected left/right E average; midpoint is the development coupling
  default, explicit endpoint is the numerical control. Standalone E3/Fortran solver
  unchanged. No fallback or additional MPI communication.
- [x] Accept updates left E only after actual-field validation. Cancel preserves
  left E and quantum state; retry histories exactly match continuous execution.
  Original endpoint-control physical/energy/work arrays remain bit-identical.
- [x] Actual authorized C2H2 GS reused, MPI4 vacuum/LD,64 steps,.003au,amplitude.1:
  applied and selected-field midpoint work numerically equal. Raw internal-energy
  residual vs selected midpoint work drops28.49%→0.695% (LD),0.688% in vacuum.
  Zero-field drift not subtracted. Corrected coupled trajectories intentionally
  change; no claim that midpoint and endpoint results are identical.
- [x] Total-field/own-vacuum/selected-field source work separately retained. Source
  continuous vs installed-wheel retry final histories match; independent quantum
  drive replay matches. Source/wheel687 tests(one optional skip); scalar finite,
  NaN/Inf and finite-large-average checks pass. Paired guides/AI JSON synchronized.
- [ ] Remaining: full spatial field/LD/PML/boundary energy balance, accounting for
  removed own-vacuum work/radiation reaction, and physical convergence/public workflow.
  Point-work alignment alone is not a total electromagnetic conservation certificate.
  Electronic restore still replays accepted history. Groups1/2 stay deferred.

Evidence:docs/isolated-midpoint-validation-2026-11-22.md and canonical C2H2 campaign
midpoint-work-summary.json. Follow-up still uncommitted; no manual heavy CI.

## 2026-11-22 六面エネルギー収支の接続

- [x] Fortranの読み取り専用観測：六面流束、場の離散／時刻中心化エネルギー、LD・外部電流の仕事。人工解で符号・領域分割の加法性と非有限値拒否を確認。
- [x] C2H2＋等方LDの実接続：観測前との全保存配列が完全一致。ステップ30取消・再試行を通した六面収支を取得。自己場仕事と量子エネルギーを含む集計もFortran。
- [x] この追加観測を含むwheel導入・実結合試験。ソース・wheelとも688件実行、各1件任意スキップ。連続／ステップ30取消再試行の全保存配列・六面観測が完全一致。
- [x] 実FDTDの複数場ランクでの面境界検証。x/y/z各2分割（外側MPI6）で場1ランクと分子の全保存配列が完全一致、各時刻の六面収支は丸め誤差内で一致。同じ分割のwheel取消・再試行は全保存配列・六面記録とも完全一致。人工解での領域分割加法性とは別に実測。
- [ ] テンソル・界面平滑化・任意の背景誘電率への観測拡張。現状は警告して観測不可を返し、計算は継続可能。

今回の実測は量子仕事の約0.70%の残差を残す短時間開発検証。通常の時刻中心化エネルギーの残差も保存し、離散保存量と区別する。一般利用の分子結合ドライバ・広い物理条件の検証は引き続き残件。

## 2026-11-22 リスト照合（5d6892c7後）

49707ce3・5d6892c7は独立SALMONpyへpush済み。古い節の「未コミット」「全エネルギー未着手」は、その時点の記録であり現時点の残件ではない。

優先候補は、開発試験にある分子双方向結合を一般利用できる入口へ整理すること。現在の電子状態の取消は履歴再生なので、再計算なしの直接復元は別項目として残す。収支観測のテンソル・界面平滑化・任意の背景誘電率への拡張と、広い物理条件の収束確認も未完了。施設・多ノード、大規模化計測、他環境配布、Agフィット改善は上記キューのまま。GPU等は将来対象、利用者指定の高速HSE電子数診断・実材料格子収束は後回しを維持する。

## 2026-11-22 単一分子の公開結合入口

- [x] 既存dipole_sourceモジュールへMolecularFeedback/MolecularFrameを追加。呼出元が用意する二つの場・一つのE駆動分子状態・MPI・DipoleMapを接続。startの初期保持とstepの中点予測／取消／電流注入／受理を公開。
- [x] 開発用MPIドライバの重複した結合ループを公開API呼出へ置換。数値処理は既存Fortranのまま。正常時の全体通信は従来の4回／区間、復元結果の全体合意は失敗時だけ。退避記録の代表ランク単独失敗は更新前の群内合意へ伝える。
- [x] 既存C2H2／等方LD・MPI6：公開API連続、wheel取消再試行、代表ランク退避記録失敗後の再試行が、以前の全保存配列・六面記録と完全一致。完全サンプルの観測値・エネルギー鮮度も一致。
- [x] 最終wheel・真空例・全体テストを照合。ソース・wheel各689件実行、各1件任意スキップ。真空の従来全配列、LDの従来全配列・六面記録、完全サンプルの観測・鮮度が完全一致。API・日英ガイド・AI JSON・サンプルのwheel収録バイトも一致。
- [x] 本計算と自己真空参照の格子・境界・採取／電流配置を初期化時の既存合意で照合。不一致は更新前に拒否。

一分子・一点・自己真空場差し引きの範囲。分子状態・場の生成終了とファイル保存は呼出元が担当し、SCF・ビルド・ダウンロード・再試行を暗黙に実行しない。電子状態の取消は履歴再生であり直接復元ではない。残件は直接復元、テンソル等の収支観測拡張、物理条件の収束、施設・多ノード・他環境、大規模化、Agフィット改善。

公開結合入口の記録：docs/molecular-feedback-api-validation-2026-11-22.md、既存C2H2キャンペーンのpublic-feedback-release-summary.json。今回の変更は未コミット。追加の物理CIは起動していない。

## 2026-11-22 孤立系ワーカー内の直接退避（第一段階）

- [x] PZワーカーのmain_tddft所有の一状態退避を実装。波動関数両バッファ・密度・ポテンシャル履歴・Poisson・RT履歴・エネルギー・XC作業量・外場ベクトルと更新番号を直接コピー。新しいモジュール大域セレクタは追加しない。
- [x] 初期0／奇数7／偶数8／容量20で直接復元。プロセスIDは同一、復元による時間発展追加0回。再試行後の9観測量と主要二出力が連続実行と完全一致。別ワーカーの状態は不変。二重退避／退避なし復元は状態を変更せずステータスを返す。
- [x] Git追跡ファイルだけから準備・ビルドし、ソースと全アプリケーションオブジェクトの禁止ソルバー監査を通過。NLCC XC・力の局所MPI1/2検証後にソース識別子を更新。既存C2H2 GSを再利用し、新規SCFなし。
- [x] 通常E3応答とバッチ基準の比較誤差は従来と同じ。Python公開API／FDTDの取消はこの段階では履歴再生のまま。
- [ ] 公開Python snapshot／FDTD取消へ直接退避を接続する。
- [ ] プロセス終了後にも使える直接ディスク保存／復元。

本節はワーカー内メモリ退避の開発検証であり、公開APIの直接保存完了とは数えない。公開API互換動作・source/wheelの最終確認と記録は完了。ソース・wheel各689件実行、各1件任意スキップ。新ワーカーを使った既存のLD＋FDTD取消再試行も、以前の全保存配列・六面記録と完全一致。

ワーカー内直接退避の確定記録：docs/isolated-savepoint-validation-2026-11-22.md と既存C2H2キャンペーンのnative-savepoint-final-summary.json。公開Python snapshot／FDTD取消への直接接続と耐久ディスク保存は未完了。今回の変更は未コミット。追加の物理CIは起動していない。

## 2026-11-22 公開保存点と分子FDTDの直接取消

- [x] IsolatedSavepointとsavepoint／restore_savepoint／release_savepointを公開。数値状態はFortran所有の保存点から直接コピーし、Pythonは出力位置と呼出契約を管理する。履歴再生のsnapshot／ディスク保存は別の明示APIであり、直接復元の代替として暗黙に呼ばない。
- [x] MolecularFeedbackの取消を直接復元へ接続。正常時の外側通信は区間あたり4回を維持。releaseの失敗も既存の受理通知に載せる。
- [x] C2H2＋LD、外側MPI6・場2分割。source連続／取消再試行と最終wheel取消再試行の全保存配列・六面記録が従来の連続実行と完全一致。同一PID、受理63更新、取消試行1更新、復元追加0更新。
- [x] 一時的な出力巻戻し失敗後の明示的完了、古い／別状態の保存点拒否を確認。ワーカー終了時は再生せず利用不能・original_preserved=Falseを通知し、呼出元MPIは存続。プロセス消失からの状態保持は保証しない。
- [x] ソース／最終wheel各690件、任意スキップ各1件。日英ガイド・AI JSON・サンプルの収録バイト一致。検証記録はdocs/isolated-python-savepoint-validation-2026-11-22.md。
- [ ] 次段階：プロセス終了後も使える直接ディスク保存／復元。現在の履歴再生アーカイブとは区別する。

既存の許可済みGSを再利用し、新規SCFなし。変更は未コミット。追加の物理CIは起動していない。

## 2026-11-22 分子の直接ディスク保存／再開

3af42587は独立SALMONpyへpush済み。以下はその後の変更で、未コミット。

- [x] save_direct_checkpoint／restore_direct_checkpointを公開。Fortranが波動関数両バッファ・一時軌道・密度／ポテンシャル・RT履歴・Poisson・エネルギー・XC作業状態・外場・更新番号を直接ストリーム保存。出力ファイルとPython観測／履歴メタデータは別に保持。履歴をSTEP_Eで再計算しない。
- [x] 形状／dt／セル／完了印、ファイルSHA256、入力・実行ファイルの一致を照合。新しい候補プロセスで読み込み、出力再接続を完了してから切替。破損・候補失敗で元のPID・時刻・観測値を保持する。旧形式への救済分岐なし。
- [x] 許可済みC2H2 GSを再利用。0／7／8／20ステップで元プロセスを終了後に再開し、復元直後の時間発展0回、続行した9観測量と出力ファイルが連続実行に完全一致。PZ・固定イオン・ワーカー1ランクの範囲。
- [x] 配列保存処理を通常TDDFTと共通化。数値処理文は既存と同一で、所有するサイズ上限の参照だけ関数化。Git追跡ファイルだけから新規ビルドし、禁止ソルバーのソース／全オブジェクト監査を通過。NLCC XC・力MPI1/2、GS不要のC／Fortran問い合わせ契約を確認。
- [x] ソース・最終wheelとも691件実行、各1件任意スキップ。wheelで実プロセス再開と破損候補拒否を実測。サンプルの途中保存・別実行での再開が、連続実行の全保存配列と完全一致。日英ガイド・AI JSON・サンプルの収録バイト一致。
- [ ] 次候補：分子の保存を本計算の場・自己真空参照場・各計算の時間ステップ番号と処理段階と協調させた一括ディスク再開。今回の分子単体保存と区別する。

検証文書：docs/isolated-direct-disk-validation-2026-11-22.md。
既存キャンペーンのdirect-disk-verified-summary.jsonとdirect-disk-verified-artifactsに記録。
新規SCF・追加の物理CIは実行していない。異なる実行ファイル／容量／コンパイラをまたぐ可搬性、
孤立系の内部MPI、広い物理条件の収束は今回の保証範囲に含めない。

## 2026-11-22 最新リスト照合（直接ディスク再開後）

現在の実装・実測記録と照合済み。旧節の未完了項目は重複して数えない。

- 完了：分子の公開双方向結合、直接取消、プロセス終了後の分子単体の直接ディスク再開。
- 次候補：本計算の場・自己真空参照場・分子・各計算の時間ステップ番号と処理段階を揃えた一括保存／再開。
- 物理／観測：テンソル・界面平滑化・任意背景誘電率の六面収支観測、分子結合の物理条件の収束確認。
- 運用／対象拡大：施設MPI・多ノード、規模別メモリ／通信／I/O計測、Linux・Intel Mac等の重いライブラリ配布・実行証拠、Agフィット改善。
- 別枠：GPU・追加保存形式などは必要性に応じた設計／検証対象。現在の必須完了条件に自動追加しない。
- 利用者指定の後回し：高速HSEの実測電子数診断と実材料の格子収束キャンペーン。

3af42587はpush済み。今回の直接ディスク再開と文書・リスト更新は未コミット。
検証集計はdirect-disk-verified-summary.json、ソース／wheel各691件・任意スキップ各1件。

## 時間ステップ番号による処理の記述

「結合時刻」「時間位置」という略記は使わず、各計算の時間ステップ番号と
電場・電流を渡す処理を明示する。現在のMolecularFeedbackでは次の対応となる。

| 処理段階 | 本計算・自己真空参照のFDTD | 分子TDDFT |
|---|---|---|
| 受理済み状態 | ステップn | ステップn−1 |
| 電場の予測 | 一時的にn→n+1へ進め、電場を取得 | n−1を保持 |
| 分子の更新 | 予測した電場を渡す | n−1→nへ進め、分極電流を返す |
| 電流を反映した場の更新 | ステップnへ戻し、分極電流を使ってn→n+1へ進める | nを保持 |
| 受理完了 | ステップn+1 | ステップn |

初期保持後はFDTDがステップ1、分子がステップ0。
一括保存／再開の次段階では、受理完了後の両FDTDと分子のステップ番号、
次の更新で参照する保存済み電場を揃えて扱う。未受理の予測途中と区別する。
これは現行処理の説明であり、一括ディスク再開の実装完了を意味しない。

## 2026-11-22 二つのFDTD・分子の一括直接保存／再開

54c43c2b（分子単体の直接ディスク再開）は独立SALMONpyへpush済み。
以下はその後の実装で、未コミット。

- [x] MolecularFeedback.save_checkpoint／restore_checkpoint。本計算FDTD n・自己真空参照FDTD n・分子 n−1を保存し、両場と分子の候補を検査後に切替。次の分子更新で使う電場、A/E採取値、出力ファイル、任意の呼出元記録も保持。
- [x] Fortranの読み取り専用sp_eh_step_indicesで、現在と候補の実ステップ番号を確認。保存・復元時だけ利用し、通常更新の通信は変更しない。
- [x] 許可済みC2H2 GS・PZ固定イオン、MPI6（各FDTD群2ランク）。30／30／29で保存して元ジョブを終了し、別ジョブで30／30／29へ直接復元。分子の追加時間発展0回、残り34更新だけで最終64／64／63へ到達。
- [x] 連続実行・変更前の場ライブラリ・最終wheel再開の全64フレームと六面収支が完全一致。最終FDTDのE/H・LD・PML・電流の全ランク保存バイナリと点出力ファイルも完全一致。
- [x] 分子の波動関数・密度等の有効な保存配列と観測値は一致。独立ジョブ間では未使用Qp_e、一時軌道の更新フラグ、Poisson作業配列の未使用末尾にバイト差があることを明示。保存ファイル全体の一致とは区別して記録。
- [x] FDTD5／分子4まで進んだ復元先へ参照場準備失敗・分子準備失敗・保存番号とFortran候補番号の不一致を注入。全ランクで元の数値状態の全バイトと番号を保持し、再試行後の結果は連続実行と完全一致。復元中は場・分子の更新呼出を禁止した検証でも成功。
- [x] Git追跡ファイルだけから準備・ビルド、ソースと全アプリケーションオブジェクトの禁止ソルバー監査。NLCC XC・力MPI1/2とGS不要のC／Fortran問い合わせ確認。ソース・最終wheel各692件、各1件任意スキップ。日英ガイド・AI JSON・サンプルの収録バイトも照合。
- [ ] 次の主な残件：テンソル・界面平滑化・任意背景誘電率への六面収支観測拡張、分子結合の物理条件の収束確認。施設・多ノード・大規模化・他環境配布・Agフィット改善は既存の別枠。

記録：docs/molecular-bundle-checkpoint-validation-2026-11-22.md、既存キャンペーンの
molecular-bundle-summary.jsonとmolecular-bundle-artifacts。追加SCF・物理CIの手動起動なし。
同一入力／実行ファイル／ライブラリ／容量／MPI分割の範囲であり、任意のOS・MPI障害や
切替時の致命的確保失敗への原子性、レイアウト変更、他の分子XCまでの証明ではない。

## 2026-11-22 正の背景ε・μと対角LDの六面収支

f4cc1b03（一括直接保存／再開）はcommit済み、利用者の今回の指示にpushは含まれないため未push。
以下の観測拡張は未コミット。

- [x] Fortranの読み取り専用観測で、実際の電場・電流・磁場の更新係数から正の対角エネルギー重みを取得。物質IDで係数を推定しない。LDの半陰的寄与を相殺して背景εを取得し、µでHの保存量を重み付けする。ソルバー、LD更新順、保存状態と通信回数は変更しない。
- [x] 人工解：非等間隔格子、方向別ε・µ、LD／外部電流の仕事、全方向の静的エネルギー、領域分割加法性、0／NaN／Inf係数の拒否を確認。失敗時の出力は保持。
- [x] 実FDTD：誘電体ε=2.5・µ=1.4、等方LD背景ε=4、対角LD背景ε=(2,3,4)、x成分だけ1の対角LD、従来の単位背景。全64ステップの観測が利用可能。最大離散残差は約2.89e-14～4.13e-14 Hartree。
- [x] 観測有無のE/A履歴と全ランク最終E/H・LD・PML・電流保存バイナリが完全一致。従来ライブラリとの対角LD状態と標準背景の観測も完全一致。既存C2H2・MPI6の全64フレーム／収支も不変。
- [x] 場MPI1/2はE/A完全一致・収支は丸め誤差内。ステップ32保存から別ジョブで直接再開した全履歴・収支・最終状態は連続実行と完全一致。ソース・最終wheelでも同じ結果。
- [x] Git追跡ファイルのみの準備・ビルド、共有カーネルのバイト一致、禁止ソルバー監査。NLCC局所MPI1/2も確認。ソース・wheel各693件、任意スキップ各1件。日英ガイド・AI JSON・収録バイトを同期。
- [ ] 次：非対角テンソルと界面平滑化の空間補間を含むエネルギー／流束の導出と実装。Ohmic伝導の散逸項も別途必要。これらが未対応でも観測だけ警告し、計算は継続する。

記録：docs/diagonal-field-energy-validation-2026-11-22.md。
calc2026/salmonpy-diagonal-field-energy-20261122/summary.jsonとartifacts-verified。
場だけの新規検査ではGS／SCFなし。分子の回帰は許可済みC2H2 GSを再利用。
時刻中心化した通常のエネルギー残差は離散保存量の残差と別であり、物理精度の保証と混同しない。

- Non-diagonal LD pole energy observation with diagonal background: implemented and verified with field MPI1/2, observation invariance and exact direct restart (see tensor-pole-field-energy-validation-2026-11-22.md). Non-diagonal background and subpixel energy remain pending.

## 2026-11-22 最新リスト照合（非対角LD極の観測後）

この節を現在の残件として扱い、過去の未完了チェックを重複計上しない。

- [x] 背景が対角の非対角LD極：全電流を使った仕事量観測、場MPI1/2、観測有無・変更前ライブラリとの状態一致、ステップ32からの直接再開一致。ソース693件（1件スキップ）。
- [x] 今回の非対角LD極の変更を含むwheelをsdistから構築し、別環境へ導入。ソース版と全64ステップの履歴・収支・全ランク最終状態バイトが一致。ソース保存点から別MPIジョブで直接再開した結果も一致。wheel693件（1件スキップ）。
- [ ] 非対角の背景誘電率：領域境界をまたぐ構成則のエネルギー交換を導出・実装・検証する。
- [ ] LD界面平滑化の六面収支観測：分散と空間補間を含む保存量を評価する。
- [ ] Ohmic伝導の散逸量を収支へ加える。
- [ ] 分子双方向結合の広い物理条件での収束確認。現在のC2H2/PZ短時間試験を一般の精度保証にしない。
- [ ] 運用・対象拡大：施設MPI／多ノード、大規模時のメモリ・通信・I/O計測、Linux／Intel Mac等の重いライブラリ配布・実行、Agフィット改善。

高速HSE電子数診断と実材料格子収束キャンペーンは利用者指定で後回し。
GPU・追加保存形式は必要性に応じた別枠。後方互換・TTMは追加しない。
f4cc1b03はローカルcommit済み・未push。その後の対角背景／非対角LD極観測と文書更新は未コミット。

### 非対角背景の独立会計カーネル（2026-11-22）

- [x] Fortranで旧／新の非対角電場エネルギーと、Poynting流束と区別した境界交換量を評価。対称な非局所演算子の人工解で、中点仕事量との恒等式、領域分割の加法性、全領域での交換相殺、NaN／Inf拒否時の出力保持を確認。694件（1件スキップ）。
- [ ] 実FDTDの背景補間と観測出力への接続、およびMPI／直接再開での比較。独立カーネルの検証だけで非対角背景の対応完了とはしない。

導出：docs/plans/2026-11-22-cross-background-energy.md。検証ログ：calc2026/salmonpy-cross-background-energy-20261122。新規GS・SCF・物理CIなし。

## 2026-11-22 最新リスト照合（非対角背景の実FDTD接続後）

- [x] 実際の背景補間演算子をFortran観測へ接続し、非対角エネルギーと `constitutive_outward_energy` を評価。Poynting六面流束とは別に保持し、残差・集計には両方を含める。現行観測18成分／集計11成分に統一し、旧契約への分岐なし。
- [x] 材料内部を切る領域の64ステップ：最大離散残差6.09e-13 Hartree、構成則境界交換は最大4.31e-7 Hartree。場MPI1/2、x/y/z分割が丸め誤差内で一致。
- [x] 観測有無・変更前ライブラリとの最終保存状態が完全一致。ステップ32の保存から別ジョブで直接再開し、全履歴・全ランク最終状態が連続実行と完全一致。従来の等方／対角の既存観測値・状態は不変、新項目は0。
- [x] 現行目録からのソース準備、ソース／全アプリケーションオブジェクトの禁止ソルバー監査、sdistからのwheel構築・別環境導入、最終wheelでの実行／再開。ソース・wheel各694件、任意スキップ各1件。日英ガイド・AI JSON・収録ソースのバイト一致。
- [ ] 次：LD界面平滑化の収支観測、Ohmic散逸の観測、分子結合の広い条件での物理収束。施設MPI／多ノード・規模計測・他環境配布・Agフィット改善は別枠の残件。

非対角背景の観測時のみ、担当MPI群内で準備合意と補間の境界通信が増える。通常更新・全体通信・保存状態は変更しない。新規GS／SCF・物理CIの手動起動なし。
記録：docs/cross-background-energy-validation-2026-11-22.md、calc2026/salmonpy-cross-background-energy-20261122/fdtd-summary.json。今回の変更は未コミット。

## 2026-11-22 最新リスト照合（LD界面平滑化の観測後）

- [x] LD界面平滑化の実係数を六面収支観測へ接続。背景の非対角境界交換をPoynting流束と分離し、平滑化された全LD電流の仕事量を評価。既存の初期化済み判定を再利用し、追加の全体通信・保存形式変更なし。
- [x] 格子からずらした曲面LD、544平滑化サンプル、参照／形状フォールバック0。最大離散残差7.41e-13 Hartree。格子整列界面も確認。物理的な広帯域精度・格子収束とは区別する。
- [x] 場MPI1/2とx/y/z分割が一致。観測有無・変更前ランタイムの最終保存状態が完全一致。ステップ32から別MPIジョブで直接再開した履歴・収支・最終状態も完全一致。既存テンソル／等方の対照は不変。
- [x] ソース／オブジェクトの禁止ソルバー監査、sdistからのwheel構築・別環境導入、最終wheelからの直接再開、ソース・wheel各694件（任意スキップ各1件）。日英マニュアル・AI JSONとwheel収録バイトを照合。
- [ ] 次：Ohmic伝導の散逸量を収支へ追加。分子結合の広い条件での物理収束、施設MPI／多ノード、規模計測、他環境配布、Agフィット改善は既存の別枠。高速HSE電子数と実材料格子収束キャンペーンは利用者指定で後回し。

記録：docs/subpixel-field-energy-validation-2026-11-22.md、calc2026/salmonpy-subpixel-field-energy-20261122/summary.json。新規GS・SCF・物理CIなし。今回の変更は未コミット。

## 2026-11-22 最新リスト照合（Ohmic散逸の観測後）

- [x] Fortranで実更新係数から導電率を取得し、LD半陰的寄与を差し引いて `ohmic_work` を評価。LD分極への仕事・Poynting流束・構成則の境界交換と分け、残差／集計に含める。現行観測19成分／集計12成分。後方互換分岐なし。
- [x] 人工解で正の散逸、ゼロ導電率、LD二重計上の防止、領域分割加法性、NaN／Inf・負係数の拒否時の出力保持を確認。
- [x] 実FDTDの導体・導体＋LD・単位背景導体で64ステップを確認。最大離散残差1.87e-14〜4.47e-14 Hartree。混在計算の場MPI1/2、x/y/z分割、32ステップからの別ジョブ直接再開が一致。観測有無・変更前ランタイムの状態は完全一致。
- [x] ゼロ導電率の等方・非対角背景・界面平滑化の既存観測値／状態は不変。新項目は0。通常更新・保存形式・通信回数は変更なし。
- [x] 禁止ソルバーのソース／オブジェクト監査、sdist→wheelの別環境導入、ソース・wheel各695件（任意スキップ各1件）、日英ガイド・AI JSON・wheel収録ソースの一致を確認。
- [ ] 残る主要作業：分子結合など広い物理条件の収束確認、施設MPI／多ノード、規模別メモリ・通信・I/O計測、他環境の重いライブラリ配布・実行、Agフィット改善。高速HSE電子数・実材料格子収束キャンペーンは利用者指定で後回し。追加保存形式・GPUは必要性に応じた別枠。

記録：docs/ohmic-field-energy-validation-2026-11-22.md、calc2026/salmonpy-ohmic-field-energy-20261122/summary.json。新規GS／SCF・物理CIの手動起動なし。今回の変更は未コミット。

## 2026-11-23 最新リスト照合（観測の資源負荷測定）

- [x] 開発用検査に格子数・観測間隔・段階別計時・安全なRSS／出力ファイル計測を追加。入力不正はMPI初期化前に拒否。計測値の全体集約は終了後1回で、計算核・通常更新は変更なし。
- [x] 単一Macの40³／72³、場1／4ランク、観測off／毎ステップ／8ステップ間隔を各2回、計24件。最大5 MPIランク、1スレッド。72³・場1ランクの観測時間中央値2.733→0.361秒、点出力10,520→1,448 bytes。一方、個別ピークRSS合計は約1007MiBでほぼ不変。実行時間のばらつきが大きく、速度向上率や最適ランク数は断定しない。
- [x] 同じ格子／MPI配置の物理状態は一致。保存した観測最大値2個の差を区別し、それ以外の保存バイトは一致。資源計測失敗を注入しても計算完了・状態一致。間引いた値は該当する毎ステップ観測と完全一致。
- [x] 日英MEMORYガイド・AI用JSONを更新。ソース・wheel各696件（任意スキップ各1件）、最終wheelの実行・sdist内の測定コードと文書の一致を確認。
- [ ] 残件：実際の保護付き混在計算での多数点・長時間の資源測定、施設MPI／多ノード・共有ファイルシステム、他環境の重いライブラリ配布・実行、広い物理条件の収束確認、Agフィット改善。高速HSE電子数・実材料格子収束キャンペーンは利用者指定で後回し。GPU・追加保存形式は別枠。

記録：docs/field-observation-profile-2026-11-23.md、calc2026/salmonpy-field-observation-profile-20261123/summary.json。今回の変更は未コミット。新規GS／SCF・遠隔ジョブ・CI手動起動なし。

## 2026-11-23 パフォーマンスモード（FDTD＋内蔵LD）

- [x] NativeEHField(execution_mode="performance") を追加。既定はchecked。毎ステップの全配列走査・A/J一致検査・検査用MPI合意・自動復元用コピーを省く。値の正しさ・呼出順は呼出元の責任とし、自動復元／再試行を保証しない。初期化・SHA・禁止ソルバー監査、最小限のローカルABI検査、物理通信、ソルバー処理、指定出力は維持。
- [x] 非ゼロ電流、通常／空間所有点配送、MPI1/2/4、導体＋LD、界面平滑化、直接保存／別ジョブ再開、明示的観測の一致を確認。通常／性能の履歴・最終保存状態は完全一致。
- [x] 72³・4場ランク・64更新を各2回。更新時間中央値6.710→5.379秒、個別ピークRSS合計1567.8→1050.8MiB。局所的な短い試験であり一般化しない。
- [x] 全730目録資産の追跡・一致、Git-only準備／ビルド、ソース／オブジェクト監査、wheel導入・697件（任意スキップ1件）、C/Fortran問い合わせ・NLCC局所MPI1/2、日英COUPLED・AI JSONの同期を確認。
- [ ] 次：独立LD／SBE、その後TDDFT／高速HSE・混在結合外側へ、同じ明示的な実行モードを展開。現時点で場以外の検査が外れるとは扱わない。

記録：docs/performance-field-mode-2026-11-23.md、calc2026/salmonpy-performance-field-20261123/summary.json。今回の変更は未コミット。

### 2026-11-23: standalone LD/SBE performance mode

- Implemented per-state and batch `execution_mode="performance"`; checked remains default.
- LD scalar/tensor arithmetic preserved; per-step backup/scans omitted.
- SBE directly evolves owned density without trial allocation; A drive unchanged.
- Independent states and explicit snapshot/restore tested against checked mode.
- Initial/library/output validation retained; automatic failed-step rollback is not promised.
- Validation: source and installed wheel each 699 tests, 698 passed / 1 optional skip; synthetic comparisons exact. Installed scalar/tensor LD and SBE examples matched. No remote CI or GS collection. See `docs/ld-sbe-performance-validation-2026-11-23.md`.

## 2026-11-23 直列の場＋SBE性能モード

- [x] 4種類の場の状態と `maxwell.couple` に明示的な実行モードを接続。初期化・契約検査は維持し、更新中の有限値走査・時刻比較を省略。数式・演算順は変更なし。
- [x] 非ゼロ・多成分の駆動による4種類の場＋SBE、2点の1D配置で通常モードとのA/J完全一致。場の明示的snapshot/restore後も一致。自作部品は性能モードへの明示的対応を要求し、通常モードは変更不要。
- [ ] 次：TDDFT／高速HSE内部とMPI結合外側への性能モード展開。今回の直列対応でこれらの検査まで省かれるとは扱わない。施設MPI・多ノード等の既存残件も維持。
- [x] 最終ソース／別環境wheel各700件（699成功、任意スキップ1件）。同梱例の横波・縦波・3点1Dを各4000ステップ、time/A/J全配列が完全一致。日英FIELDS・AI JSON・例のwheel収録バイト一致。記録：docs/serial-coupling-performance-validation-2026-11-23.md。未コミット。

## 2026-11-23 通常TDDFT・単一状態の性能モード

- [x] NativeTDDFTBackend(execution_mode="performance") を単一状態の呼出元MPI経路へ接続。Fortran時間発展は通常経路と共通。更新時の全状態退避・診断合意・異常復帰走査を省き、初期化・明示的保存／復元・物理通信を維持。
- [x] 既存GSを読取専用のリンクで再利用。PZ／TB-mBJ／Libxc-PZ／Libxc-PBE × 内部1／2ランクで通常・性能・変更前ライブラリの電流／実測電子数が完全一致。メモリ内復元後も一致。新規GS・GSコピー・ディスク電子状態保存なし。
- [x] ソース目録・全Git追跡資産のみからの準備／新規ビルド、禁止ソルバー監査。Python契約テストで検査通信の省略・有限値検査の省略・失敗後の再利用拒否を確認。
- [ ] 次：共有／登録モデルの通常TDDFT、高速HSE、分子・別プロセス型、MPI結合外側へ性能モードを展開。今回の単一状態対応をこれらの対応済みとは扱わない。
- [x] 最終ソース／wheel各701件（700成功、任意スキップ1件）。各版8条件のMPI比較と明示的メモリ内復元が一致。wheel/sourceの配列・文書バイト一致。記録：docs/tddft-performance-validation-2026-11-23.md。未コミット。

## 2026-11-23 共有・登録TDDFTの性能モード

- [x] ハンドル更新に性能モードを追加。ローカルの世代・容量・有効期間の検査を維持し、診断合意と自動退避を省く。通常経路の合意・異常復帰は維持。
- [x] モデルの既定値と状態ごとの上書きを分離。通常セッション／登録管理器の再利用でも各物質設定の選択を伝達する。
- [x] 4種類のXCの共有状態、PZ↔TB-mBJ／Libxc-PZ↔PBEの異種登録を内部1／2ランクで独立対照と比較。交互更新・兄弟状態の保持・メモリ内復元後の結果が完全一致。既存GSは読取リンクだけ、新規GS・コピー・電子状態のディスク保存なし。
- [ ] 次：高速HSE、別プロセス型・分子、MPI結合外側への性能モード展開。施設MPI・多ノードなどの既存残件は維持。
- [x] 最終ソース／wheel各702件（701成功、任意スキップ1件）。各版12条件のMPI試験で独立実行・明示的メモリ内復元後の電流一致、兄弟状態の保持、閉じたハンドルの拒否を確認。wheel/sourceの配列・日英文書・AI JSONバイト一致。記録：docs/shared-tddft-performance-validation-2026-11-23.md。未コミット。

## 2026-11-23 高速PR版HSEの単一状態性能モード

- [x] NativeFastHSEBackend(execution_mode="performance") を単一状態の呼出元MPI経路へ接続。電子／交換状態の自動退避と復帰用診断合意・走査を省き、初期化・明示的保存／復元・必要な物理通信を維持。
- [x] 学習用データ収集の位置を通常／性能の両方でエネルギー出力前に維持。状態所有の一時フラグで診断設定と処理順を分離。縮約対称性とAの数学的な条件は維持。
- [x] Taylor／PT-CN × 全／縮約k点 × 内部1／2ランクの8条件で、通常／性能／変更前の電流・エネルギー・学習履歴が完全一致。明示的メモリ内復元後も一致。縮約で不適合Aを拒否して正常続行。
- [x] 既存GSはリンクから読み取り、新規GS・GSコピー・電子状態のディスク保存なし。Git-only準備／新規ビルド、ソースと全アプリケーションオブジェクトの禁止ソルバー監査。
- [ ] 次：高速HSEの共有／登録状態、別プロセス型・分子、外側MPI結合へ性能モードを展開。既存の施設MPI・多ノード等の残件も維持。
- [x] ソース／wheel各703件（702成功、任意スキップ1件）。各版8条件で電流・エネルギー・学習履歴とメモリ内復元後の一致、縮約対称性違反の更新前拒否を確認。wheel/sourceの配列・文書・サンプル・アダプターバイト一致。記録：docs/fast-hse-performance-validation-2026-11-23.md。未コミット。

## 2026-11-23 高速HSEの共有・登録状態の性能モード

- [x] 世代付きハンドルの性能更新を追加。ローカルな世代・容量・対称性条件と必要なモデル切替を維持し、診断用ハンドル合意を省く。通常更新の検査・異常復帰は維持。
- [x] モデルの既定モードと状態ごとの選択を分離。共有セッション・登録管理器の再利用でも各物質設定の選択を伝達。日英ガイド・AI JSON・共有／混在サンプルを更新。
- [x] 最終の共有8条件・異種登録4条件のMPI比較、wheel導入版の照合を完了。既存GSは読取リンクのみで再利用し、新規GS・コピー・電子状態のディスク保存は行わない。
- [ ] 次：別プロセス型・分子、外側MPI結合への性能モード展開。施設MPI・多ノードなどの既存残件は維持。
- [x] 最終ソース／wheel各704件（703成功、任意スキップ1件）。各版12条件で独立対照との電流／エネルギー／学習履歴、メモリ内復元、兄弟状態保持、スロット再利用後の旧世代拒否を確認。wheel/sourceの配列・日英文書・AI JSON・サンプル・アダプターは一致。記録：docs/shared-fast-hse-performance-validation-2026-11-23.md。未コミット。

## 2026-11-23 別プロセスの周期系TDDFT性能モード

- [x] TDDFTInitialState.new_state／crystal.respond／TDDFTResponseBackendに性能モードを接続。専用初期化コマンドと応答で実行ファイルの対応を確認し、黙って通常モードへ下げない。
- [x] 更新時の有限値検査と失敗後の自動再起動・履歴再計算を省略。初期化・通信形式／ステップ番号・容量・SHA照合・明示的出力を維持。通常経路と分子の経路は維持。
- [x] PZ／TB-mBJ／Libxc-PZ／Libxc-PBEで通常・性能・変更前の電流／電子数が完全一致。モード間snapshotを明示的な履歴再計算で復元して受取側モードを維持。ワーカー強制終了後の自動再起動なし、旧実行ファイルの拒否も確認。
- [x] 最終ソース／wheelの照合と検証記録を完了。実GSの収集を避けるため検証時だけコピー関数をGSリンクへ置換し、更新・通信・復元は製品コードを使用。GS内容は前後で同一。製品コピー処理は変更しない。
- [ ] 次：分子・高速HSEの別プロセス型、外側MPI結合への性能モード展開。施設MPI・多ノード等の既存残件も維持。
- [x] ソース／wheel各705件（704成功、任意スキップ1件）。各版4種類のXCで電流・電子数一致、明示的なモード間履歴再計算の復元、ワーカー停止時の自動再起動なし、旧ワーカー拒否を確認。wheel/sourceの配列・日英文書・AI JSON・サンプルは一致。記録：docs/process-tddft-performance-validation-2026-11-23.md。未コミット。

### 2026-11-23 分子E駆動ワーカーの性能モード

- [x] 専用の初期化・応答プロトコルと状態ごとのモード指定。
- [x] C2H2の20ステップで通常・性能モードの9観測項目が完全一致。
- [x] メモリ内復元とプロセス終了後の直接ディスク復元で継続結果が完全一致。直接復元後の時間発展回数は0。
- [x] GSは読み取りリンクで再利用し、内容不変。新しいGS計算・コピーなし。
- [x] この変更を含むwheelでの実計算確認。705件成功・1件スキップ。記録: `docs/molecular-performance-wheel-validation-2026-11-23.md`。

検証ログ: `/private/tmp/molecular-perf-native.log`。物理的収束や速度向上を示す試験ではない。

## 2026-11-23 外側A→J結合の性能モード

- [x] ResponseSessionの固定モードと全ランク一致の初期化契約。部品ごとのモードとは独立。
- [x] A/J有限値走査、有限値検査専用の群内合意、CoupledStepperの重複した開始合意を性能モードで省略。形状・所有者・必要な通信とコールバック例外の合意は維持。
- [x] 3/5ランク・2物質群・24ステップのMultiscale1D＋自作物質モデルで通常／性能のA/Jが完全一致。明示的な保護更新で途中例外後の状態保持・再試行を確認。ソース／wheelで各707件中706成功・1スキップ。
- [x] 日英ガイド・AI JSON・MPIサンプルを更新。記録：calc2026/salmonpy-outer-performance-20261123。GSなし、CI手動起動なし。
- [ ] 次：実際のE/H FDTD＋電子状態で外側モードを検証し、A/E/Hの場の受け渡しと分子の外側スケジュールへ適用範囲を拡張。今回の合成試験をその検証の代用にしない。

通常の性能更新には自動復帰を保証しない。明示的なstep_protectedは性能モードでも退避・復元を行う。速度向上の測定は別。本節と次節の変更を今回のコミットへ収録。

## 2026-11-23 実E/H＋LD＋2点TDDFTの外側性能モード

- [x] 既存ライブラリのSHA照合。GSは読み取りリンクで再利用し新規計算・コピーなし。
- [x] 3/6ランク、FDTD1/2分割、TDDFT各群1/2ランク、96ステップの通常連続実行・通常／性能の保護復元後がA/E/JおよびLD領域の非ゼロ電場波形で完全一致。ソース6条件・wheel6条件。
- [x] 65番目の更新中に電子状態を変化させた後の例外を全群へ返し、64ステップ終了時へ直接退避から復元して再試行。
- [x] wheel707件中706成功・1スキップ、日英・AI JSON・配布ツール一致。記録: docs/outer-eh-tddft-performance-validation-2026-11-23.md。
- [ ] 次：高速PR版HSEとの実E/H結合、場のA/E/H受け渡し・分子の外側性能スケジュール。今回のPZ合成軌道試験をこれらや物理収束の証拠にしない。

本節の検証・文書を今回のコミットへ収録。CIの追加手動実行なし。施設MPI、多ノード、規模測定などの既存残件は維持する。

## 2026-11-23 高速PR版HSE＋E/H＋LDの外側性能モード

- [x] 既存の高速PR版ライブラリをSHA照合し、既存GSを読取リンクで再利用。GS計算・コピー・電子状態ディスク保存なし。
- [x] 学習付きTaylor／PT-CN×全／縮約k点×3／6ランクのソース24・wheel24条件。2物質群、FDTD1／2分割、HSE各群1／2ランク、96ステップ。
- [x] 通常連続実行と、65番目の更新中の例外から64ステップ終了時へ直接復元・再試行した通常／性能モードで、A/E/J・LD領域波形・HSE電流／エネルギー出力が完全一致。出力97行で失敗更新は残らない。
- [x] 全／縮約k点とMPI分割間も記録した許容誤差内で一致。wheel707件中706成功・任意スキップ1。日英・AI JSON・検証ツールを更新。記録: docs/outer-fast-hse-performance-validation-2026-11-23.md。
- [ ] 次：A/E/Hの場の受け渡し・分子の外側性能スケジュール。HSE別プロセス型の必要性確認、規模計測・施設MPI／多ノード・他環境・物理収束等は既存の別枠。

未commit/push。重いCIの追加手動起動なし。今回のメモリ内復元試験をディスク／プロセス終了後の再開試験の代用にはしない。

# 項目5：共有資源の棚卸しと完了条件（2026-11-20）

この文書は現行CPU native経路のソース調査です。新しい実計算やコンパイラ検証の記録ではありません。
過去の未完チェックをすべて別タスクとして数えず、下記の具体的な資源単位で進めます。
本体の `SAVE` を機械的に全廃することは完了条件にしません。

## 現行の対象範囲

標準TDDFTのPZ/TB-mBJ/Libxc-PZ/PBE、固定イオン・周期系・A駆動、
高速PR HSEの検証済み全／縮約k点、GS保持状態、native E/H＋LDの結合です。
nativeはメインスレッドで逐次呼出し、ライブラリのMPI接続は1プロセスに1つ。
複数物質状態は同じ親コミュニケータの下で独立に持ちます。
FDTDの複数インスタンス／プロセスは今回の所有権整理に自動追加しません。
GPU・分子関連の再開や、未対応SOC/+U/DC等の制限解除も今回の作業に含めません。

## 資源単位の分類

| 資源・ファイル | 所有者／切替箇所 | 判定・残作業 |
|---|---|---|
| GSの数値・混合履歴、`gs/main_dft.f90`、GS storage/model | GSハンドル。継続／保存／破棄 | 測定済み範囲は完了。特殊プロファイルを勝手に広げない |
| 標準RT `rt_context`、`rt/embedded_tddft.f90` | ハンドルのcontext、pending/step_savedもcontext内 | 数値所有は対応済み。一般の初期化失敗復帰は別行 |
| 高速HSE `tddft_context`、`native_adapters/tddft_context.f90` | 軌道・XC・交換・学習履歴はcontext | 対応済み。checkpointサイズ作業値も呼出しローカル |
| HSE初期化の固定設定記録、`fast_hse_mpi.f90` | 登録モデルのinit_settings／同一モデル生成器 | 対応済み。世代再利用時リセット、失敗時不変 |
| `state/model/active_model`、handles/models、世代番号 | 標準/HSE C入口のレジストリと一時選択子 | 接続単位の管理資源として維持。数値配列の共有と混同しない |
| `salmon_global` と parser/model_inputs | 入力の実体はモデル。select_modelが既存kernel用の値を設置 | 逐次呼出し用互換ビュー。全フィールドの完全な引数化は将来の再入可能化 |
| `tddft_worker` embedded_dt/steps/decomposition/initial_A | embedded初期化用の互換ビュー。標準select_model/worker_configure | 既存プロファイルの逐次初期化では維持。汎用・並行初期化には明示引数化が必要 |
| `pending_gs`、GS→RT境界 | seedを登録呼出し中だけ借用し、戻りでnullify | 通常の戻りは解放済み。未返却STOP/OOMの回復とは区別 |
| `filesystem.io_root`、`_native_group._directory_call` | 呼出し範囲で借用。例外時はfinallyで解除 | 対応済み。旧runtimeの互換chdirは残す |
| timerの実行中バンク、標準/HSEのlegacy/handle banks | 生存状態のバンクと一時ホスト退避 | 検証済みプロファイルで対応済み。同時スレッド実行とはしない |
| `parallelization.nproc_*`、gap_time_system | 接続の親groupとその時計 | 同じ接続内で共有してよい。2つの親groupの同時接続は別設計 |
| communication.application_comm/embedded_comm、child registry | MPI接続。RTは初期化／更新時に状態ownerを指定 | TDDFT/HSEとFDTDは対応済み。FDTDはR2の1/2ランク測定で確認 |
| 対称性テーブル、symmetry/symmetry.f90・symmetry_rho.f90 | HSEモデルのworkspaceをselect_modelで交換 | 検証済みプロファイルは対応済み。GPU・他プロファイルへ一般化しない |
| XC/stencil作業配列、salmon_xc/code_optimization | system%xc_work／stencil、標準モデルのworkspace | 対応済み範囲を維持 |
| exx_native.rt_owned／互換aliases | bind_systemされた所有者を優先し、旧呼出しのみfallback | native绑定済み経路と未対応fallbackを区別。全legacy機能の書換えはしない |
| exx_native.wannier | 一部の局在交換／legacy経路 | 現行独立RTのwhitelistから外れる機能。今回有効化しない |
| DOS/PDOS/band/AVSの一時unit、step_output | 空きunit／状態のofl | 既知の固定unit問題は修正済み。任意I/O失敗の原子性は別行 |
| native FDTD state、step_saved、pending | `maxwell/embedded_eh.f90`の単一field寿命 | pendingのclose時解放はR1で修正・検証済み |
| native FDTD `eh_mpi_grid_sr`、fs%srg_ng | 初期化でinit_communicator_dftを呼ぶ | R2で状態owner指定・cache終了後の回収を実装。親・他所有者保持を測定済み |

## 次に実装する作業（順序付き）

### R0：指定済みスカラーIEEE検査規約への適合

`maxwell/embedded_eh.f90`の入力／段階検査と、
`native_adapters/hse_checkpoint.f90::valid_history`に、
`all(ieee_is_finite(array))`が残っています。現在のソースに対する確認事項です。
GNUで過去の試験が通ったことを、Fujitsu対応や規約適合の証明にしません。
次の関連ビルド前に明示ループ／scalar inquiry helperへ直します。

受入条件：全要素を走査し、精度・閾値・MPI全体の拒否を維持する。
有限値、NaN、Infと通常結果の回帰を確認する。Fujitsu実機／コンパイラが
利用できない場合は未検証と明記し、コンパイル済みとは報告しない。
新しい入力制限や収束停止条件は加えません。

### R1：FDTDの終了時に復元候補を確実に解放

`sp_eh_checkpoint(mode=2)`はmoduleのpendingをallocateします。
mode3/4は解放しますが、`sp_eh_close`はstep_savedとstateのみを解放し、pendingを残します。
Python側が通常discardしていても、native closeの契約で解放すべき資源です。

受入条件：候補あり／なし、早期終了／完走で数値状態の一時所有を解放する。
未完了の計算へ最終スペクトル等を出すためにeh_finalizeを強制実行しません。
途中終了時の出力unitについては実際のOPEN/CLOSE対応を確認し、未確認の漏れを断定しません。
連続実行・保護付き復元・LD結合の結果を維持します。

### R2：FDTDの子MPI資源をfield寿命に帰属

`sp_eh_init -> eh_init -> eh_mpi_grid_sr -> init_communicator_dft`は
TDDFTと異なりcomm_begin_childrenでfield ownerを設定していません。
`sp_eh_close`もcomm_release_childrenを呼びません。現在は接続のdetachが最後の回収地点です。

受入条件：fieldにownerを持たせ、生成・必要な更新範囲をそのownerへ帰属する。
通信cacheの終了後に子groupを回収し、親groupと他の所有者を維持する。
1／2ランクでinit/step/close前後のlive数を測り、早期終了・完走・候補保持を含める。
MPI分割やE/H/LDの更新順は変更せず、旧版との結果比較を行う。

### R3：初期化失敗の対象を段階ごとに定義

`initialization_rt`、`eh_init`と下位I/OにはSTOP／未捕捉のallocation/I/O失敗が残る。
最初に対応する失敗段階と後始末を一覧化し、戻せるものだけstatus契約を定める。
SCF未収束・メモリ予算超過は引き続き警告とし、計算可能なケースの停止条件を増やさない。
任意のMPI障害・OOM・外部ライブラリSTOPからの回復を完了条件にしない。

### R4：native ABIの安定化は独立した完了項目

現行ヘッダーは版付きでもexperimental/previewです。符号・時刻・所有権・
失敗時保証・寿命・能力照会と非対応範囲を統一し、C/Fortranホスト例で確認してから
安定性を宣言します。項目5全体を閉じる前に必要であり、この棚卸しだけで完了にしません。

## 項目5の完了の扱い

まずR0〜R2を閉じ、現在対応する計算経路の資源寿命を固めます。
R3の対応段階とR4の安定ABIが未完の間、項目5全体のチェックは開いたままにします。
全グローバル変数の排除、同時native実行、未対応物理機能の解放、
任意のプロセス障害回復は別の拡張であり、この作業へ際限なく追加しません。
その後の優先順5→3→6、グループ1・2の保留、分子の保留は維持します。

棚卸し時点ではソースの読取りと文書更新のみで、計算・ビルド・CIは実行していません。

追記：R0はその後実装・GNU/MacおよびMPI検証済みです。配列IEEE問い合わせを除去し、
有限値・NaN・正負Infと旧結果の一致を確認しました。Fujitsuコンパイルは未実施です。
R1の終了時候補解放とR2のMPI所有権も修正・検証済みです。次の対象はR3とR4です。詳細は2026-11-20-r0-scalar-ieee.mdを参照。

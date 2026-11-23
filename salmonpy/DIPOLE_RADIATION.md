# 分子双極子とFDTDの結合

## 単一分子の双方向結合

`MolecularFeedback` に、作成済みの二つのFDTD群と一つの孤立系TDDFT状態を渡します。全ランクで同じ順序で呼び、分子状態は `total_root` のランクだけに置きます。場の担当ランクでは `field`、それ以外では `None` を渡します。`reference_root` は、この分子だけが駆動する真空FDTD群の代表ランクです。この参照には入射波・他の源・物質を入れません。格子・境界条件・採取位置・電流の配置方法は本計算と揃え、開始前に照合します。MPI、分子状態、FDTD群の作成・終了は呼び出し元が担当します。

```python
from salmonpy.dipole_source import DipoleMap, MolecularFeedback

mapping = DipoleMap(('molecule',), (10,), [[1.]], [8.])
coupling = MolecularFeedback(
    world, field, state=state, mapping=mapping,
    total_root=0, reference_root=group_size,
    policy='environment_scattered_only', observer=observer,
)
first = coupling.start()  # 電流ゼロで場を1ステップ進める
for _ in range(steps - 1):
    frame = coupling.step()
    if world.rank == 0:
        records.append(frame.record())
```

`world` は呼び出し元のMPIコミュニケータ、`group_size` は各FDTD群と補助ランクを合わせた大きさです。`state` は代表ランク以外では `None`。`observer` は各場の `energy_observer(...)`、観測不要なら全ランクで `None` です。`start()` の戻り値も必要に応じて保存します。`MolecularFrame` は代表ランクだけに返り、他のランクでは `None` です。電場、電流、双極子、電子数、鮮度を反映した内部エネルギー、任意の六面収支を保持し、`record()` は独立した保存用辞書を返します。自動的なファイル出力は行いません。

結合は一分子・一点、実験室座標の双極子、中点E駆動です。初期の選択電場はゼロとし、非ゼロなら光源の距離・遅延を見直します。場の時刻は電子状態よりdt先です。自己真空場の差し引きは真空中の放射反作用も除きます。分子状態と場を結合の外から途中更新しないでください。

受理端点の不一致や、代表ランクだけで実行する `step(validate=callback)` の例外は、その区間を取り消します。全ランクで `StepRecoveryError` を捕まえ、`original_preserved=True` のときだけ同じ `step()` を再試行できます。自動再試行はありません。コールバックは観測の確認だけに使い、状態を変更したりMPI通信を行ったりしません。電子状態は直接退避で復元し、履歴再生しません。ワーカー喪失など復元を確定できない場合はoriginal_preserved=Falseとして返し、再試行を認めません。

完全な実行例は [molecular_feedback.py](examples/molecular_feedback.py) です。明示した既存GSと検証済みのE駆動計算プログラムを使い、SCF・ビルド・ダウンロードを暗黙に開始しません。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 \
  python -m salmonpy.examples.molecular_feedback FIELD_BUILD E_WORKER RECEIPT EXISTING_GS NEW_OUTPUT \
  --field-grid 2 1 1 --ld
```

この例の格子・時間長はC2H2の接続検証用です。任意分子の収束済み条件ではありません。

## 本計算・自己真空参照・分子の一括保存／再開

全ランクで`MolecularFeedback.save_checkpoint(directory, history=...)`を呼びます。
時間更新の受理後だけ保存し、本計算FDTDと自己真空参照FDTDはステップn、
分子TDDFTはステップn−1であることを確認します。予測更新の途中は保存しません。
FDTDの実ステップ番号はFortranから読みます。

```python
# 初回のジョブ。recordsはtotal_rootにだけ置く。
first = coupling.start()
records = [first.record()] if world.rank == total_root else None
while coupling.step_index < 30:
    frame = coupling.step()
    if world.rank == total_root:
        records.append(frame.record())
coupling.save_checkpoint('new-bundle', history=records)
# 各状態を終了し、MPIジョブを終了してよい。
```

```python
# 別のジョブで同じ入力・配置からFDTDと分子、couplingを新規作成する。
# start()は呼ばず、保存した番号へ直接復元する。
saved = coupling.restore_checkpoint('new-bundle')
records = saved['history'] if world.rank == total_root else None
while coupling.step_index < steps:
    frame = coupling.step()
    if world.rank == total_root:
        records.append(frame.record())
```

保存形式は`salmonpy.molecular-feedback-direct/1`。二つのFDTDのE/H・LD・PML・
電流状態を各場ランクがFortranで直接保存し、分子は直接ディスク保存APIを使います。
採取したA/E、次の分子更新に使う保存済み電場、出力ファイルも保持します。
`history`は任意の呼出元のJSON記録です。指定する場合は受理済みFDTDステップ数と
同じ長さのリストをtotal_rootだけから渡します。復元結果はtotal_rootだけに返り、
`field_step`、`molecular_step`、各担当ランクの`verified_steps`と`history`を含みます。

同じ実行ファイル／ライブラリ・GS入力・dt・容量・物質・境界・格子・採取位置・
MPI分割で再開します。SHA256と各設定を照合後、両FDTDと分子を独立した候補へ読み、
Fortranの候補ステップ番号と分子の更新回数0を確認してから切り替えます。
候補の準備失敗は全ランクへ伝え、元の状態を保持します。
任意のOS／MPI障害や切替時の致命的な確保失敗まで回復を保証するものではありません。
保存は新規ディレクトリのみで、既存保存を上書きしません。共有ファイルシステムが必要です。
通常の時間更新には通信を追加せず、保存／復元時だけ照合します。
復元中は元の状態と候補が同時に存在するため、追加のメモリと出力退避領域が必要です。
この機能には`sp_eh_step_indices`を持つ場のライブラリと、revision5の分子ワーカー検証記録が必要です。
読込み時にビルド・SCF・過去ステップの再計算は行いません。

完全サンプル`examples/molecular_feedback.py`に
`--stop-after 30 --checkpoint NEW_CHECKPOINT`を指定して保存・終了し、
別のMPIジョブで`--resume NEW_CHECKPOINT`を指定して再開できます。
両方とも同じ`--steps`（容量）・入力・MPI数を使い、出力先はそれぞれ新規にします。
`--ld --field-grid 2 1 1`ではMPI6です。分子はPZ・固定イオン・内部1ランクの実測範囲です。


## 一方向の放射源

分子の物理的な電荷双極子pを実験室座標、e·bohr単位で渡します。
DipoleMapは双極子総量を保存して明示したマクロセルへ分配します。

P_i = sum_m weights[i,m] p_m / cell_volumes[i]

j_i = sum_m weights[i,m] (p_m(new)-p_m(old)) / (dt cell_volumes[i])

体積はFDTDセルのbohr³単位の値で、TDDFTの真空を含む計算箱の体積ではありません。
重みは非負で、分子ごとの和を1とします。自動規格化しません。
分子IDは重複のない文字列、場の点IDは非負整数で、指定した順序を保持します。
複数分子を同じセルへ重ねられます。座標回転と実際の分子数の反映は呼出側で明示してください。
均質な分子集団ならP=n*<p>ですが、数密度・配向平均を計算箱から推測しません。
時刻は原子単位、jはe/(bohr²·原子時間)です。SIへの暗黙変換はありません。

```python
from salmonpy.dipole_source import DipoleMap
mapping = DipoleMap(('molecule-1',), (10,20), [[.25],[.75]], [8.,8.])
j = mapping.current([[0.,0.,0.]], [[0.,0.,.001]], dt=.003)
```

入力はコピーされ、返す配列は読取専用です。一定の双極子は電流を作りません。
その静電場を初期化する機能ではありません。電気双極子近似で、分子の空間電荷分布、
磁気・高次多極子放射、局所場・自己場の帰還は扱いません。

## FDTDのE/H更新への接続

`drive_native_dipoles(field, mapping, dipoles, dt=...)`は初期状態のNativeEHFieldを借用します。
場グループの全ランクから同じ引数で呼び、他のグループからは呼びません。
双極子は(samples,molecules,3)、時刻は0からn*dtです。
読取専用のE/j履歴、それぞれの時刻、分子ID・点IDを返します。
Eはt=0を含み、jは区間平均を中点時刻で表します。
マクロセルからYee格子への分配とFDTD更新は既存計算処理処理のままです。
Python側でもう一度Yee格子へ分配しないでください。

計算処理の初期化条件により最初の区間電流はゼロが必要です。
すぐ変化する既存波形を使う場合は、呼出側で初期双極子を1回繰り返します。
これによるdtの遅延をAPI内に隠しません。

```python
import numpy as np
from salmonpy.dipole_source import drive_native_dipoles
# p: n*field.dtの電荷双極子、形状(samples,molecules,3)。
# field: 点ID・セル体積が一致する初期化済みNativeEHField。
delayed = np.concatenate((p[:1], p), axis=0)
result = drive_native_dipoles(field, mapping, delayed, dt=field.dt)
```

場はlen(delayed)-1ステップ以上の容量を持ち、点の順序・体積が一致する必要があります。
時間刻みも一致させ、補間は明示的に行います。全波形を最初の更新前に検査し、ランク間で照合します。
backend異常時に計算全体を巻き戻す仕様ではありません。場のcloseは呼出側が担当します。
NativeEHFieldの明示的build record、呼び出し元プログラムが管理するMPI context、境界から離れた真空近傍、
平面波入力という制約は引き継ぎます。分子の放射だけを見る場合はパルス振幅をゼロにします。
入射場を残す場合は出力に両方の寄与が含まれます。境界・CFL条件はFDTD側で扱います。

一方向接続なので放射したEをTDDFTへ戻しません。放射反作用、自己場と局所場の区別は未対応です。
実験的孤立系E-workerでは、応答の最後の3成分である全電荷双極子を使います。
周期系TDDFTの物質電流の符号規則を転用しないでください。FIELD_INPUT.mdも参照してください。

## サンプルと検証

```sh
python -m salmonpy.examples.dipole_radiation
```

wheel内の軽量例は電流生成と総量保存だけを確認します。
ソース内のtools/test_dipole_radiation_mpi.pyは明示したembedded library recordと新しい出力先で
Pythonから直接呼ぶFDTDを実行します（合計2または3ランク）。tools/test_isolated_worker.pyが出力した
NPZを追加引数にすると分子の実測双極子も使えます。
検証記録はdocs/dipole-radiation-validation-2026-11-15.mdです。
短いC₂H₂波形は接続試験で、収束した分子発光スペクトルではありません。

## 細分化検証の現状

位置を固定した解析的な点双極子について、遅延電場の全成分と電流源の仕事を比較しました。
dt=.0015原子時間でdxを2→1→.75 bohrへ細かくすると、電場の相対誤差は46.5%→8.0%→4.3%、
源の仕事と連続体の放射エネルギーとの差は15.8%→2.5%→1.1%へ減りました。
閉曲面を通るエネルギー流の実測や、本番精度の収束保証ではありません。
時間刻みだけを半分にしても連続体解との一致は単調に改善しませんでした。
源セルのEは格子依存が強く、そのまま分子へ帰還する方法は未検証です。
再現ツールはソース内のtools/test_dipole_refinement_mpi.pyとtools/summarize_dipole_refinement.py、
条件と制限はdocs/radiation-refinement-validation-2026-11-15.mdにあります。
自己場差引きや真空の放射反作用モデルを自動的に有効化することはありません。

公開APIのisolated_responseサンプルが作った出力も、そのまま読み込めます。

```sh
python -m salmonpy.examples.dipole_radiation --archive ISOLATED_OUTPUT --cell-volume 8000
```

8000は使用するFDTDセルの実際の体積（bohr³）に置き換えてください。
単位・等間隔の時刻・分極電流生成と保存則を確認する例で、FDTD自体は起動しません。
電磁場の伝播計算は上記のMPI検証スクリプトに分けています。

## 分子へ渡す電場の明示選択

`select_molecular_field`は計算済みの電場`FieldValue`を選択する処理です。
双方向ソルバーは起動しません。方針を必ず指定します。

- `provided_drive`：利用者が選択済みの電場を保持し、provenanceの`meaning`に
  物理的な意味を記録します。参照場の併用は拒否します。
- `environment_scattered_only`：環境中の全電場から、対象分子自身だけによる
  真空参照場を引きます。外場・他分子からの場は残します。真空の放射反作用も
  除くため、完全な放射ダイナミクスや自然放出の計算とは呼びません。

両標本はE:au・bohr・point samplingとし、時刻・位置を厳密に一致させます。
参照のprovenanceには`source_ids=[molecule_id]`、`medium='vacuum'`、
`incident=False`が必要です。両provenanceの`footprint_id`、`sample_map_id`、
`own_source_history_id`も一致させます。これらは配置・観測・自己源履歴が同じという
利用者の宣言であり、別々のMaxwell計算の条件を検証済みという証明ではありません。
不足や不一致をゼロで補ったり、全電場へ暗黙に戻したりはしません。

戻り値`MolecularDrive`の`.field`と`.metadata`は不変です。`.record()`は入力・参照・
選択後のE、位置・時刻・方針・分子ID・出典を持つ独立したJSON互換データを返します。
補間・単位換算はしません。足し合わせた人工電場で動作を確認する例は次です。

```sh
python -m salmonpy.examples.molecular_field_selection
```

この選択処理自体は場を更新しません。振動子による中点結合の検証は下記の開発試験に
分けています。孤立系Eワーカーへの開発用接続は末尾を参照してください。

### 開発用の中点帰還検証（線形振動子）

2026-11-22に分子の作業を再開しました。ソースの
`tools/test_molecular_feedback_mpi.py`で、単一の線形振動子と二つの独立した場を
使います。全電場から、その振動子だけが作る真空参照場を差し引きます。
振動子更新、電流、仕事量と残差の数値処理はFortranです。これはTDDFTの代用となる
開発試験であり、実分子のモデルや公開の双方向TDDFTドライバではありません。

最初の1区間は電流ゼロで保持します。その後は各区間について、場を保存して終点の
電場を予測し、両端の電場を使って振動子を中点更新します。場を元へ戻し、補正した
電流で更新し直します。実際の終点電場を使って更新と仕事量を再検査し、予測と
整合した場合に受け付けます。試行更新の観測ファイルも巻き戻します。現在の試験は
各場1ランクと制御用ランクを使う合計4ランクで、場の大規模並列の実証ではありません。

```sh
# BUILDは既存の監査済みFDTDライブラリのbuild.json。自動ビルドしません。
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 2 python tools/test_molecular_feedback_mpi.py external BUILD NEW_BASE --dt .003 --steps 64
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python tools/test_molecular_feedback_mpi.py feedback BUILD NEW_RUN --baseline NEW_BASE --dt .003 --steps 64 --schedule centered
```

合成LD環境は両コマンドへ `--environment ld` を付けます。真空参照場にはLDと入射光を
置きません。基準の時刻刻み・ステップ数・環境・ライブラリが異なれば、更新前に拒否します。
`--schedule held` は旧アルゴリズムとの差を測る明示的な比較用です。救済分岐ではありません。

旧保持法の仕事量の時間ずれはdt=.003 auで約3.60e-11 auでした。中点帰還では同じ条件で
約2.67e-27 auになり、合成LD環境でも散乱帰還を残したまま中点仕事量の一致を確認しました。
絶対残差に加え、振動子エネルギーと仕事で規格化した残差を報告します。
これは観測点のEと双極子変化による仕事量検査です。全空間の場エネルギー・境界流束を
含む保存則、結合全体の二次収束、任意の物質や多分子への保証ではありません。
真空放射反作用は方針上除いており、自然放出も含みません。孤立系TDDFTの
試行状態・巻き戻しと端点Eの時刻を明示した接続です。

受け付けた電流列を通常のFDTD更新へ再入力する `replay` 比較でも、真空・LDとも
電場の履歴が完全一致しました。試行・巻き戻しが既存の場の更新を変えない検査です。

双極子から分極・電流への写像と自己場の差し引きも、同梱のFortranカーネルで計算します。wheel導入後はSALMON実行ファイルなしで利用でき、ソースから使う場合は軽量カーネルを明示的にビルドします。インポート時にコンパイルや通信はしません。


### 孤立系Eワーカーへの接続（開発用）

同じ試験に `--isolated-executable`、`--isolated-validation`、
`--isolated-bundle`、`--isolated-input` を指定すると、線形振動子に代えて
`IsolatedTDDFTState` を使います。必要なのは検証済みの孤立系Eワーカーと既存GS入力
一式です。周期系Aワーカーでは代用できません。SCF・取得・ビルドは自動実行しません。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python tools/test_molecular_feedback_mpi.py feedback BUILD NEW_RUN --isolated-executable E_WORKER --isolated-validation VALIDATION_JSON --isolated-bundle EXISTING_GS --isolated-input rt.inp --dt .003 --steps 64
```

場は最初の1区間を電流ゼロで保持します。電子状態の時刻ゼロをこの直後に対応させ、
`electronic_time` と場の時刻の差dtを明示して保存します。最初の自己場差引き後のEは
ゼロが必要です。非ゼロなら光源との距離・遅延を増やしてください。各区間の左端と予測右端の選択EをFortranで平均し、
その区間で一定のEとしてワーカーへ渡します。既存の時間発展順は維持します。
双極子の差から電荷電流をFortranで作ります。

場を巻き戻して電流を投入した後、実際の自己場差引き後の終点Eと予測値を照合します。
不一致なら場と電子状態を戻して拒否します。場の復元トークンは一度限りなので、予測の
取消後に確定試行を保護する新しいトークンを取ります。電子状態の取消は同じワーカーへの直接復元で、
時間発展を再計算しません。開発試験で `--inject-reject-step 30` を
付けると1回拒否して再試行します。任意の異常はこの試験ドライバではMPI終了となり、
一般の障害復旧APIではありません。

合成ワーカーに加えて、C2H2・PZ・固定イオンの実TDDFTを接続しました。電子格子は
24×24×32、幅0.35Å。許可を受けてGSを一度計算し、保存して再利用しています。
基準TDDFTとEワーカーの双極子差は最大約3.86e-12 e·bohrでした。真空・合成LDの
4ランク試験（各場1ランク＋制御ランク）で独立したE駆動応答との一致を検査します。
これは短い64ステップ・dt=.003 auの接続試験であり、光学応答の収束検証ではありません。
この接続試験に加え、後述の電子数・仕事量・内側の箱の六面収支を検証しています。

MPI版の1プロセスワーカーが親のMPI接続情報を引き継ぐと起動待ちになるため、子へ渡す
環境だけからOMPI/PMI/PMIX/PRTEの情報とランク識別子を除きます。親のMPIやジョブ割当ID
には触れません。この起動経路は今回のMac/OpenMPIで検証したもので、全サイトでの保証
ではありません。一般利用向けの現在の入口は冒頭の `MolecularFeedback` です。

E3では実測電子数も保存・復元照合し、Fortranで減少率を診断します。終点方式のC2H2＋LDの64ステップでは最大減少率約9.95e-15（目安1e-5）で、電子数を含む取消再試行履歴が連続実行と完全一致しました。短時間の診断であり、光学応答の収束保証ではありません。


E3の内部エネルギーは既存の観測間隔を保ち、更新済みの時刻だけ仕事量と比較します。
TDDFTが受け取る終点Eの区間一定仕事と、中点Eの仕事をFortranで別々に積算します。
過去の終点駆動方式の比較記録では、励起振幅.1で終点仕事との差が約0.897%、中点仕事との差が約28.49%でした。
既定の弱い励起は別測定した無電場ドリフトが支配し、精度を判断できません。
背景を自動的に差し引かず、残差で有限な計算を止めません。結合全体のエネルギー保存は
未確認です。時間配置を揃えた中点方式の結果は次節、内側の箱の六面エネルギー収支は後述します。


### 分子と場の中点時間配置

公開の孤立系結合は中点駆動です。自己場を差し引いた
左端・予測右端EをFortranで平均します。受理した右端Eだけを次区間へ継承し、取消時は
電子状態と左端Eを保持します。確定後の場から作る中点Eも照合します。自動的な終点方式
への切替はありません。終点駆動方式は過去の比較記録として残します。
単独の `IsolatedTDDFTState.step(E)` の意味やFortran本体は変えていません。

`drive_time`は場の時計で、中点方式では区間の中点、終点方式では右端です。
電子状態の`time`は更新後の区間終端です。`work_applied`と
`work_selected_field_midpoint`を区別して保存します。全電場・真空自己参照の仕事も
別記録します。自己場差引き後の仕事が一致しても、全電磁場との保存則とは限りません。

C2H2・PZ、励起振幅.1の中点比較では、適用Eによる仕事と選択電場の中点仕事が
数値的に一致し、内部エネルギーとの差は約0.695%でした。約28.49%だった終点方式の
時間配置差を解消しています。無電場ドリフトは差し引いていません。途中拒否後も
連続実行と全履歴が一致します。内側の箱の場・LD・六面流束による収支は次節で扱います。

## 箱の六面を使ったエネルギー収支（開発用）

`NativeEHField.energy_observer(lower, upper)` は、光源・吸収層・外端を除く箱を、両端を含む0始まりの格子番号で指定します。場の担当ランク全員が同じ順序で呼びます。更新前の退避状態を保持して更新した後、解放する前に `observer.sample_step()` を呼びます。場は進めません。場の担当群内の集計1回だけを行い、代表ランクに辞書、他のランクに `None` を返します。既存の受け渡しパケットに載せれば全体通信は増えません。

六面の外向き流束、場の保存量の変化、LDの仕事、外部電流の仕事をFortranで計算します。LDの仕事は分極に蓄えられるエネルギーと散逸の合計です。PML内の発熱を推定する方法ではありません。散乱光の側面からの流出も六面流束に含みます。

現在の検証範囲は真空、正の非分散背景誘電率・透磁率、等方LD・対角LDです。非負のOhmic伝導率も対象とします。LD媒質の透磁率は従来どおり1です。対角LDではx/y/zそれぞれの背景誘電率を使います。対応していない配置など、観測が利用できない場合は警告して `available=False` を返し、計算自体は止めません。

重みは物質IDから再構成せず、実際のYee更新係数からFortranで取得します。電場の係数をc1、電流係数をcJとすると、c4=−4πdt/cJ、背景誘電率はc4(1+c1)/2です。LDの半陰的な寄与が相殺されます。磁場は実際のH更新係数から透磁率を取得します。これにより、媒質境界でもソルバーが使った成分ごとの係数と一致します。標準の単位背景・等方LDでは従来の評価式を維持します。

ステップnの離散保存量は、各Yee成分についてE(n)²に背景誘電率を掛け、H(n−1/2)H(n+1/2)に透磁率を掛けた和です。時刻中心化した通常の磁場エネルギーも別に返します。これらの診断式だけを拡張しており、E/H・LDの更新順、チェックポイント状態、通信回数は変更しません。

開発検査`tools/test_diagonal_field_energy_mpi.py`はGSを使わず、場1／2ランク、正の背景誘電率・透磁率、等方／対角LD、32ステップ保存からの直接再開を確認します。実際の計算を実行するツールであり、管理ノード用の入力検査ではありません。非対角の背景誘電率には、空間補間を含む構成則の境界交換を追加します。LD界面平滑化の観測も、平滑化後の背景・分散応答を使います。

Yee格子のE/Hには半ステップの時刻差があります。`field_modified_*` はこの離散更新が保存する交差積のエネルギー、`field_physical_*` はHを時刻中心化した通常のエネルギーです。両方の残差を残し、離散収支の一致を連続体のエネルギー精度と取り違えないようにします。

自己場を引いた分子駆動では、場の変化＋六面流束＋LDの仕事＋分子内部エネルギー変化に、差し引いた真空自己場の仕事を加えます。ゼロ外場のドリフトは差し引きません。C2H2の短時間開発試験では離散場の累積残差は約9.1e-15 Hartree、分子を含む残差は約5.9e-8 Hartreeでした。分子の仕事に対して約0.70%であり、大きな入射エネルギーで規格化して誤差を小さく見せません。物理条件全般の精度保証ではありません。

実際のFDTD領域をx・y・zそれぞれ2分割する試験も通過しました（外側MPI6、場の各群2ランク＋管理1ランク）。分子の保存配列は場1ランク版と完全一致し、六面収支は集計順の丸め誤差内で一致しました。`tools/test_molecular_feedback_mpi.py` の `--field-grid NX NY NZ` で分割を指定し、全体には `2*(NX*NY*NZ+1)` ランクを割り当てます。同じ2分割での取消・再試行は、全保存配列・六面記録とも連続実行と完全一致です。電子状態の取消は直接退避です。多ノード・大きな分割数の実証とは区別します。

背景誘電率が対角の場合は、非対角のLD極強度にも対応します。仕事量には各成分の全LD電流を使い、成分間の応答も含めます。背景誘電率が非対角の場合も、実際の背景補間から領域境界をまたぐ交換量を評価します。

非対角背景の観測では `constitutive_outward_energy` を返します。これは領域外向きを正とする構成則の境界交換量で、`outward_energy`（六面のPoynting流束）と区別します。旧／新の場エネルギーには非対角寄与を加え、`modified_residual` と `physical_residual` にはこの境界交換も含めます。局所的に割り当てた非対角エネルギー単独が正とは限りません。

現行のFortran/C観測配列は19成分で、末尾2成分が構成則の境界交換量とOhmic散逸量です。集計配列は12成分です。旧契約の観測用ライブラリは使用せず、現在のソースから明示的に再ビルドします。非対角背景を観測する場合のみ、担当MPI群内で補間の境界通信と準備失敗の合意が追加されます。通常の時間発展、全体通信、保存データは変更しません。これは離散収支の検査であり、物理的収束の保証ではありません。

LD界面平滑化を指定した場合も、実際の更新係数と空間依存の背景補間を使って同じ19成分の収支を返します。非対角の背景が領域境界をまたぐ寄与は `constitutive_outward_energy`、分散応答の仕事量は全LD電流による `ld_work` に含めます。元のバルク材料係数に戻して評価しません。

この診断は、中心周波数で近似した平滑化モデル自身の離散収支を確認するものです。広帯域での材料応答や曲面の格子収束を保証しません。幾何形状や特異な参照条件で平滑化を適用できなかったセルは、計算が使った元係数で評価されます。`salmonpy-subpixel-summary.json` の適用数・元係数へ戻った数・参照誤差も確認してください。SBE／TDDFTは平滑化しません。

Ohmic伝導の散逸量は `ohmic_work` として返します。各辺の旧／新電場の平均から、1ステップあたりの `dt × 体積 × Σ σ E_mid²` をFortranで計算します。値は非負で、場エネルギーが失った分を表します。`ld_work`（LDの分極への仕事で、不可逆な熱だけではない）、外部電流への仕事、Poynting流束、構成則の境界交換と分けて記録し、残差と集計には全て含めます。

導電率は材料IDで推定せず、実際の電場／電流更新係数からLDの半陰的寄与を差し引いて取得します。係数の差し引きで生じる丸め誤差だけをゼロとして扱います。導電率ゼロの計算ではこの追加評価を省き、従来の観測値を保持して `ohmic_work=0` を返します。通常更新・保存形式・通信回数は変更しません。現行の場ライブラリと軽量カーネルは合わせて明示的に再ビルドしてください。

## 分子TDDFTの性能モード

`initial.new_state(dt=dt, execution_mode="performance")` で、E駆動の分子ワーカーも性能モードを選べます。既定は `"checked"` です。性能モードは毎ステップの有限値検査と自動回復を省きます。入力の形・容量、応答の種類、エネルギーの更新有無、明示的な保存・復元は引き続き検査します。`step_fields` の場の受け渡し契約も検査します。異常時の自動再試行は保証しません。

保存した状態に実行モードは含めず、復元先が選んだモードを保持します。ディスクから復元する候補ワーカーも同じモードで起動します。通常モードとの数値的一致を確認してから利用してください。この指定だけでは外側のFDTD結合全体が性能モードにはなりません。

### 分子の外側結合の性能モード

`MolecularFeedback(..., execution_mode="performance")` で、報告する場・電流配列の追加の有限値走査を省きます。場と分子のモードは独立で、必要なら `NativeEHField(..., execution_mode="performance")` と `initial.new_state(..., execution_mode="performance")` も指定します。サンプル `molecular_feedback` の `--execution-mode performance` は3か所をまとめて指定します。既定は検査付きです。

分子結合は、試行→場の復元→電流を含む受理更新を行うため、両FDTDと分子の退避を性能モードでも維持します。自己真空参照の差し引き、予測した電場と受理した電場の一致条件、時刻・MPI順序、明示的な検証コールバックも維持します。コールバックが更新を拒否した場合、`StepRecoveryError.original_preserved` と後処理エラーを確認してから再試行します。プロセス終了・通信障害・NaNからの一般的な自動復帰を約束するモードではありません。直接ディスク復元の候補ワーカーも復元先の実行モードを保持します。

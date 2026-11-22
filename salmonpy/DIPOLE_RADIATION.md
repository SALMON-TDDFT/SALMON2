# 分子双極子からFDTDへの一方向放射源

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
ゼロが必要です。非ゼロなら光源との距離・遅延を増やしてください。各区間で予測した
終点Eをワーカーへ渡し、既存の時間発展順を維持します。Pythonで再び中点平均は
しません。双極子の差から電荷電流をFortranで作ります。

場を巻き戻して電流を投入した後、実際の自己場差引き後の終点Eと予測値を照合します。
不一致なら場と電子状態を戻して拒否します。場の復元トークンは一度限りなので、予測の
取消後に確定試行を保護する新しいトークンを取ります。電子状態の復元は**履歴再生**で、
直接復元や再計算不要の保証ではありません。開発試験で `--inject-reject-step 30` を
付けると1回拒否して再試行します。任意の異常はこの試験ドライバではMPI終了となり、
一般の障害復旧APIではありません。

合成ワーカーに加えて、C2H2・PZ・固定イオンの実TDDFTを接続しました。電子格子は
24×24×32、幅0.35Å。許可を受けてGSを一度計算し、保存して再利用しています。
基準TDDFTとEワーカーの双極子差は最大約3.86e-12 e·bohrでした。真空・合成LDの
4ランク試験（各場1ランク＋制御ランク）で独立したE駆動応答との一致を検査します。
これは短い64ステップ・dt=.003 auの接続試験であり、光学応答の収束検証ではありません。
電子数・量子系の仕事量・全電磁エネルギー保存は未検証です。

MPI版の1プロセスワーカーが親のMPI接続情報を引き継ぐと起動待ちになるため、子へ渡す
環境だけからOMPI/PMI/PMIX/PRTEの情報とランク識別子を除きます。親のMPIやジョブ割当ID
には触れません。この起動経路は今回のMac/OpenMPIで検証したもので、全サイトでの保証
ではありません。一般利用向けの双方向ドライバは引き続き未公開です。

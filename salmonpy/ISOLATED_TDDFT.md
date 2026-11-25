# 孤立系TDDFT：E入力と双極子

PZ・固定イオン・非偏極・スカラー擬ポテンシャル・中点伝播を、明示的な別プロセスで使います。結晶のTDDFT/SBEは従来どおりA駆動です。返す量は電荷双極子であり、電流密度ではありません。格子・時間刻み・計算箱の収束は利用者が確認します。

## 計算プログラムと初期状態

```python
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState

runtime = load_isolated_runtime('salmon', 'validation/report.json')
with IsolatedTDDFTInitialState.from_salmon(
        'existing-bundle', 'rt.inp', runtime=runtime,
        max_steps=100, threads=1, timeout=60) as initial:
    with initial.new_state(dt=.02) as state:
        saved = state.savepoint()
        frame = state.step([0., 0., 1.e-5])
        print(frame.time, frame.dipole, frame.electron_number)
        state.restore_savepoint(saved)  # 同じワーカーへ直接コピー、時間発展なし
        state.release_savepoint(saved)
```

E・時刻は原子単位、双極子はe·bohrです。`dipole` はイオンと電子を合わせた全電荷双極子、`electronic_change` は電子双極子変化です。配列は読み取り専用です。初期E・時刻はゼロ。`new_state(position=...)` は場を評価する位置（既定(0,0,0) bohr）を指定し、分子座標の移動・回転は行いません。

読込みは64 KiB以下のローカル検証記録と、実行ファイルのSHA256を照合します。現行の記録形式は `salmonpy.isolated-worker-validation/5`。数値応答はE3の9成分のままですが、同一プロセス・時間発展追加0回の直接復元と、古い退避の拒否の実証を必須にしています。旧revision4以前の記録は受け付けません。読込み時にダウンロード・ビルド・検証計算・MPI初期化・ワーカー起動は行いません。記録はローカル検証の証拠で、電子署名や禁止ソルバー監査の代用ではありません。

既存の参照入力を使う `tools/test_isolated_worker.py EXECUTABLE EXISTING_REFERENCE NEW_VALIDATION` で、通常応答と直接復元を検証した記録を作成できます。これはRT計算を実行する開発用ツールです。新しいSCFは行いません。対応環境への配布が完了したという意味ではなく、手元の検証済み計算プログラムを明示して使います。

初期bundleには孤立系RT入力、`data_for_restart/info.bin`・`occupation.bin`・`wfn.bin`、相対パスで参照する擬ポテンシャル等を用意します。ローダーは不変の初期コピーを保持し、各状態が独立した作業コピーと子プロセスを持ちます。初期ローダーを閉じても既存状態は使えます。周期系GSの変換や孤立系GSの自動作成はしません。

## 同じ計算プロセス内での直接退避

`savepoint()` はFortran側の波動関数両バッファ、密度・ポテンシャル履歴、Poisson作業状態、双極子・外場履歴、エネルギー、更新番号を保持し、`IsolatedSavepoint` を返します。一つの状態につき一つだけ保持できます。`restore_savepoint(token)` は出力を元の長さへ戻して同じプロセスで復元し、時間発展を呼びません。トークンは復元後も残り、繰り返し復元できます。不要になったら `release_savepoint(token)` で明示的に解放します。解放だけなら現在の計算結果は変えません。

別状態のトークン・解放済みのトークン・再起動前のトークンは、復元前に拒否します。トークンはディスクへ保存する形式ではありません。退避を保持したまま履歴再生の `snapshot()/restore()` や `save_checkpoint()` は使えません。計算状態と退避の寿命は呼び出し元が管理し、状態を閉じれば退避も終了します。

`worker_statistics()` は明示的な問い合わせで、PID、実際の時間発展回数 `updates`、受理番号 `step`、直接復元回数 `direct_restores` を返します。通常の更新に全体通信を追加しません。復元だけでは `updates` は増えません。

直接退避中の失敗で、別ワーカーの起動や履歴再生は行いません。出力の巻き戻しに一時失敗した場合は `usable=False` として更新を止め、同じトークンで復元を完了できます。ワーカー喪失や復元結果を確定できない障害は使用不可として返します。`usable` を確認せずに更新を再開しません。任意のOS/MPI障害からの復帰を保証するものではありません。

退避の対象は現在検証している主要二出力です。未対応の空間出力を要求した場合は退避作成だけを拒否し、通常のE駆動計算は引き続き使えます。退避には電子状態一式を保持する分のメモリが必要です。

## 共通の場入力と観測

```python
from salmonpy.field_input import FieldPacket, FieldValue
frame = state.step_fields(FieldPacket(E=FieldValue(
    [0., 0., 1.e-5], time=state.time+state.dt,
    position=(0.,0.,0.), unit='E:au')))
```

この例は上記の生きている状態内で実行します。Eの時刻・位置・単位を照合し、未提供とゼロを区別します。A/Hは消費せず、AからEを暗黙に再構成しません。単独の `step(E)` は既存の定数区間Eを受ける積分順を維持します。

`frame.electron_number` はFortranが積分した実測電子数です。減少率1e-5は診断の目安であり、それだけで計算を停止しません。`frame.internal_energy` は既存の計算頻度で新しく評価された端点だけ値を返し、それ以外は `None` です。計算頻度を勝手に増やしたりゼロ外場ドリフトを自動で差し引いたりしません。

## プロセス終了後の直接ディスク再開

```python
with model.new_state(dt=dt) as state:
    for E in first_fields:
        state.step(E)
    state.save_direct_checkpoint('new-electronic-checkpoint')
# 元のワーカーは終了。同じGSと設定から新しい状態を作る。
with model.new_state(dt=dt) as resumed:
    resumed.restore_direct_checkpoint('new-electronic-checkpoint')
    assert resumed.worker_statistics()['updates'] == 0
    for E in remaining_fields:
        frame = resumed.step(E)
```

`salmonpy.isolated-direct/1`はFortranの数値配列・履歴・Poisson・XC作業状態、
受理ステップ、観測値、出力ファイルを保存します。各ファイルのSHA256と、
モデル・実行ファイルの一致を読み込み前に確認します。元のGSは別途指定し、
初期化時に読み込みます。SCFや経過済み時間ステップの再計算はしません。
別ワーカーで読み込みと出力再接続を完了してから切り替えるため、
候補の復元失敗時も元のワーカーを保持します。復元した状態の異常時に、
履歴再生へ暗黙に切り替えません。復元中は電子状態が追加で必要です。

実行ファイル・GS／入力・dt・最大ステップ数・位置・スレッド数は同一にします。
異なるコンパイラや配列配置をまたぐ可搬形式ではありません。
保存点を解放してからディスク操作を行います。max_bytesは任意の読み込み量制限で、
メモリ使用量の制限ではありません。既存ディレクトリは上書きしません。
対象は既存の確認済み出力モードです。ワーカー検証記録は実際のプロセス終了後の
再開を含むrevision5を要求し、古い記録は受け付けません。

応答サンプルで`--stop-after 7 --direct-checkpoint NEW_CHECKPOINT`を指定すると、
7ステップで終了して保存します。同じE NPZと設定で別途
`--resume-direct NEW_CHECKPOINT`を指定すると、残りだけを進めます。
これは分子単体の保存です。FDTD結合全体の再開には二つの場と結合時刻も
協調して保存する必要があります。

## 履歴の保存と再実行（直接復元とは別）

`snapshot()/restore()` は受理済みEと観測量の履歴を保存・再実行する入口です。同じ初期モデル、dt・位置・容量・スレッド数・バイナリが必要です。候補プロセスの観測量を照合してから交換し、復元失敗時は元の状態を保持します。費用は履歴長に依存し、検証中は新旧二ワーカー分の資源を使います。直接退避を使わない通常の `step()` の異常回復もこの履歴再実行方式です。

```python
state.save_checkpoint('new-history-checkpoint')
resumed.restore_checkpoint('new-history-checkpoint', max_bytes=16*1024*1024)
```

`resumed` は同じ初期データ・設定で呼び出し元が作成した状態です。ディスク形式は `salmonpy.isolated-tddft-checkpoint/3`。E0、初期／現在の9観測量、E履歴と観測履歴をハッシュ付きで保存するZarr v3互換形式です。GSは別途用意します。`max_bytes` は追加のファイル読み込み予算であり全RSS上限ではありません。これは波動関数の直接チェックポイントではありません。直接の数値保存には上記のsave_direct_checkpointを使います。

## 出力とFDTD

`python -m salmonpy.examples.isolated_response EXECUTABLE VALIDATION_JSON BUNDLE RT_INPUT E_NPZ NEW_OUTPUT` は、E[steps,3]とスカラーdtを原子単位で持つNPZから双極子アーカイブを作ります。`--checkpoint` は明示した履歴再実行用の保存です。応答・単位・ソースの出典を明示し、結果はソースツリー外へ保存します。

双方向FDTD結合は [分子とFDTD](DIPOLE_RADIATION.md) の `MolecularFeedback` と、完全例 `examples/molecular_feedback.py` を使います。取消に `savepoint/restore_savepoint/release_savepoint` を使い、電子履歴を再計算しません。真空自己場の差し引きは真空の放射反作用も除きます。

既存C2H2 GS・PZ、64場ステップの6ランク検証では、正常時は電子更新63回、1区間の取消・再試行時は64回で、各実行内のPIDは変わりません。復元による追加更新は0回、受理した全保存配列・六面記録は以前の結果と完全一致しました。これは短時間の接続検証で、光学物性の収束保証ではありません。

## 入力で指定する双極子源の近接場

SALMON入力の `num_dipole_source`、`cood_dipole_source`、`vec_dipole_source`、
`rad_dipole_source` を使う入口は、`from_salmon(..., native_near_field=True)` で明示します。
既定は無効です。この入口には、使用する実行ファイルについて近接場の比較結果を記録した
検証記録が必要です。従来の一様電場だけの検証記録では有効になりません。
C₂H₂の短時間検証で本体との双極子比較、独立状態、同じプロセス内の復元と、プロセス終了後の直接ディスク再開を確認しました。物理量の収束検証とは区別します。

対応範囲はPZ・固定イオン・スカラー擬ポテンシャル・1ランクです。
双極子源は最大2個、半径は正の有限値、波形は `Acos2` または `Ecos2` です。
座標・半径はSALMON入力の単位系に従い、Pythonの `step(E)` は引き続き原子単位です。
入力の双極子源の空間分布と時間波形を保持し、`step(E)` は**追加の一様電場**を渡します。
源の数が0なら一様電場だけです。源の入力や有効化設定が異なる状態間では復元できません。

これは入力で定めた双極子源の機能です。一般の空間依存E/Hを受け取る入口とは区別します。

## 空間依存E/Hを使う近接場モード

`from_salmon(..., near_field=True)` は、分子周囲のE/Hを経路積分して使います。
現在の対象は **PZ・固定イオン・非偏極・スカラー擬ポテンシャル・CPUの1ランク**です。
磁気結合は軌道に対するものです。スピン・SOCの磁気結合、分子からFDTDへの放射の戻しは含みません。
分子が置かれる領域は真空とし、Gaussian原子単位で `B=H` を使います。

### FDTD格子から直接渡す（展開しない方法）

`GridFieldSeries` は、保存したE/Hの格子データを直接補間し、Iwasaらの電場の経路積分へ渡します。
多極子の次数や減衰率を選ぶ必要はありません。現在は等間隔の3次元直交格子で、セル内は3重線形補間です。
経路が格子セルを跨ぐ場所で積分を分割します。磁場による軌道結合も同じ補間場から計算します。

```python
from salmonpy.near_field import GridFieldSeries

# E_samples、H_left、H_rightは、同じセル中心で保存したSampledField。
# H_leftはn-1/2、H_rightはn+1/2。相対時刻0の元のFDTD時刻も保持する。
grid = GridFieldSeries.from_yee_samples(
    E_samples, H_left, H_right,
    counts=(16, 17, 17), origin=(7., -15., -15.), spacing=(2., 2., 2.),
    molecule_position=(21.2, .7, .9), center=(0., 0., 0.), time_origin=.303)
grid.write('NEW_INPUT_BUNDLE/salmonpy-near-eh.dat')  # 保存GSを再利用する入力準備のみ

with IsolatedTDDFTInitialState.from_salmon(
        'NEW_INPUT_BUNDLE', 'rt.inp', runtime=runtime,
        max_steps=20, near_field=True) as initial:
    with initial.new_state(dt=.003) as state:
        for _ in range(7):
            state.step([0., 0., 0.])  # 追加の一様E。格子の場は保存時系列から取得。
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

位置・間隔はbohr、時刻は原子単位です。`counts`は各軸のサンプル点数で、FDTD全体の格子数とは区別します。
`origin`は最初のサンプル位置です。`molecule_position`は積分中心のFDTD座標、`center`は同じ点のTDDFT座標で、
Fortranが `FDTD位置 = molecule_position + (TDDFT位置-center)` として座標を対応させます。
どちらの位置も格子点と一致する必要はありません。座標軸は共通とし、回転は自動では行いません。原子だけでなく、TDDFT計算箱・差分の外側格子点と
非局所ポテンシャルの経路まで覆うパッチを用意します。格子外の推測や暗黙の周期折り返しはしません。
分子領域は真空としてH=Bを使います。

HをEの時刻へ揃えて保存済みなら `from_samples(E_samples, H_samples, ...)` を使います。
構造化配列を渡す場合は `GridFieldSeries(origin, spacing, times, E, H, ...)` で、E/Hは
`[時刻,z,y,x,成分]` です。成分ごとの格子位置が異なる場合は、`origin`を `[E/H,成分,座標軸]` の
配列で指定します。各成分は同じ点数・間隔を持つ共通のパッチに切り出します。
半ステップの構造化H配列には `from_yee_arrays(..., H_left_times=..., H_right_times=...)` を使えます。
いずれも時刻合わせ・空間補間・積分はFortranで行い、Hや端点を省略して推測することはありません。

入力サイズは時刻数×格子点数に比例するため、必要な局所パッチと区間を切り出します。
格子の時系列は相対時刻0から一定dtで、入力フレーム0〜max_steps+1を用意します。入力のdtとTDDFTのdtを揃えます。
ステップnではフレームnとn+1を使い、`step(E_extra)`で進めます。`step_near_field`は展開係数用なので使いません。
元のFDTD時刻は`time_origin`として入力ファイルにも記録します。この情報で分子の過去の時間発展を
再現したことにはなりません。初期状態と照射開始の条件を明示してください。
保存／再開では格子データ・座標対応と電子状態を保持し、同じ入力モデル・dt・計算プログラムで
`restore_direct_checkpoint`後に残りのステップを進めます。

準備用サンプルは `python -m salmonpy.examples.prepare_grid_field --help` です。準備だけで計算は始まりません。
直接格子入力の検証記録を持つ実行ファイルが必要です。磁気結合は元のIwasaらの電気相互作用とは
別に検証しており、スピン磁気結合や分子からFDTDへの放射の戻しは対象外です。
補間してもFDTDで解像できない場は復元できないので、格子と積分・分子応答の収束を確認します。

直接補間の解析場試験では、成分ごとにずれた格子で線形場・経路積分の絶対誤差は最大約2.1×10⁻¹⁷でした。
滑らかな場を格子幅0.5→0.25→0.125 bohrで調べると、電場の相対RMS誤差は6.60→1.46→0.328%、
結合ポテンシャルは3.42→0.646→0.0901%へ減少しました。実際のLD/FDTD保存場をC2H2へ渡す試験では、
分子位置をFDTD格子からずらし、ステップ7の保存から20までの全9観測量と出力が連続実行と完全一致しました。
これは接続と再開の検証であり、実ナノ構造・分子応答の物理的収束を保証するものではありません。

### フィット時に中心を決める

```python
from salmonpy.near_field import NearFieldSeries

# E_samples/H_samples: 同じ位置・時刻のSampledField。単位はE:au/H:au、位置はbohr。
near = NearFieldSeries.from_samples(
    E_samples, H_samples, center=(0., 0., 0.), radius=8., lmax=1)
print(near.center, near.relative_residual)  # 残差の列はE、H
near.write('GS_BUNDLE/salmonpy-near-eh.dat')  # 準備だけ。計算は始めない。
```

中心は分子中心などを明示指定することを勧めます。省略時はサンプル点の重心、半径の省略時は
そこから最も遠い点までの距離を、最初のフィットでFortranが決めます。
決まった中心・半径を全フレームと保存／再開で共有し、実行途中に選び直しません。
半径は基底の規格化尺度です。FDTDの分解能を細かくする操作ではありません。

E/Hの各成分を、半径方向の項を含む実Ylmで展開します。`lmax` は総次数dの上限です。
係数は `d=0..lmax`、各dで `l=d,d-2,...`、各lで `m=-l..l` の順に並び、
`(r/radius)^d Y_lm` を掛けます。d=l+2kなので半径因子r^(2k)も含みます。
項数は `P=(lmax+1)(lmax+2)(lmax+3)/6`、配列は `[時刻,P,xyz]` です。
電気ポテンシャルだけの入口とは係数の単位が異なります。
`lmax=0` は空間一様場、`lmax=1` は勾配を含む場です。1次には電気四重極までが含まれ、
一様Hは軌道磁気双極子に結合します。より高い次数も使えますが、有限次数で任意の場を厳密に表せるわけではありません。
フィットの相対残差が既定の5%を超えると警告し、計算可能な結果は返します。

FDTDのHが半ステップずれる場合は `NearFieldSeries.from_yee_samples(E_samples, H_left, H_right, ...)`
を使います。H_leftはn−1/2、H_rightはn+1/2、Eはnです。
位置・時刻の対応を確認し、Hの平均をFortranで計算します。初期の前後のHも明示的に渡してください。
欠けた端点の推測や外挿はしません。

準備用の例は `python -m salmonpy.examples.fit_near_field --help`。
実行例 `examples/isolated_response.py` には `--near-field` があります。
フィット・保存・importはTDDFTを起動しません。


分子位置と展開中心は実座標で指定し、FDTD格子点との一致は要求しません。
周囲の既知のE/Hサンプルから、その中心での展開を作ります。格子点に分子を移動したり、
最近傍の場だけを採用したりしません。現在は保存した場をTDDFT入力として使う一方向の処理です。
分子応答をFDTDへ戻す補間・電流配置・エネルギー整合性は別の設計課題です。
中心から分子の電子密度が広がる領域をサンプルが十分に覆うようにし、金属界面をまたぐ
単一の滑らかな展開や、サンプル領域外への外挿を精度保証と扱わないでください。

TipやBow Tieに必要な次数は未確定です。次数だけでなく、先端半径・ギャップ・分子位置・
フィット領域を変え、フィットに使わない点で絶対・相対残差を確認します。
現在の基底は指定した総次数までの完全多項式を表せます。ただしFDTDで解像できない情報は
復元せず、Maxwell方程式の整合性も自動では保証しません。独立した位置の残差とFDTD分解能を確認します。


従来の診断では角度成分だけの表現と完全多項式を比較しました。TDDFT入力も現在は半径項を含みます。
格子点での残差が小さくても、格子間の場や分子応答の精度は別途確認が必要です。
総次数8では各成分165係数になり、少なくとも165点の独立なサンプルが必要です。
点数・独立性不足は収束として扱いません。入力ファイルV3・直接保存V7で表現を識別し、
以前のファイルは作り直します。

半径項付きの接続検証では、C2H2の既存GSを使い、非ゼロの2次半径項を含むE/Hで20更新しました。
更新7で直接保存し、プロセスを終了してから再開した13更新の全9観測量は連続実行と完全一致し、
復元時の追加時間発展は0回でした。独立した解析場では総次数2・4・6・8の格子間点のE/H誤差は
最大約1.3×10⁻¹⁵ auです。これは表現・接続の検証であり、Tip/Bow Tieの物理的収束保証ではありません。


### 指数減衰を含む場

既知の方向に急減衰する場には、同じ半径項付き基底へ指数因子を掛けられます。

```python
near = NearFieldSeries.from_samples(
    E_samples, H_samples, center=(0., 0., 0.), radius=8., lmax=4,
    decay_vector=(0., 0., .25))  # 単位bohr^-1。正のz方向に減衰。
near.write('GS_BUNDLE/salmonpy-near-eh.dat')  # 入力準備だけ。
```

各E/H成分の基底は `exp(-decay_vector·(r-center)) (|r-center|/radius)^d Y_lm` です。
半径は多項式部分の規格化尺度で、減衰ベクトルの単位は変えません。
既定の `(0,0,0)` は通常の多項式展開です。E/Hで共通の減衰ベクトルを使い、中心・半径・次数と同様に、
全フレームと保存／再開で固定します。`NearFieldFrame.from_samples` にも同じベクトルを指定します。
減衰率の自動選択は行いません。別の点で選んだ率を使う場合も、最終検証点は分けて残差を確認します。
選ばれた率が物理的な減衰長を一意に示すとは限りません。単一因子は一般のTip/Bow Tieの全近接場を
表すものではなく、逆方向には増大するため、分子を覆う検証済み領域で評価してください。

場の近似と分子への結合積分はFortranで分離しています。指数因子を含む場合も、電場は経路上のEを、
磁場は経路パラメータを掛けたHを積分して結合ポテンシャルを作ります。指数積分では小さい引数の
桁落ちを避け、磁気の非局所ポテンシャル位相は指数依存に対応した数値積分を使います。
別の場表現を後から渡すための、経路上のサンプルから積分する共通処理もあります。
FDTDからの直接補間や双方向接続まで実装したという意味ではありません。


自作の場表現には `salmonpy.near_field.couple_path_samples(displacement, parameters, weights, E, H)`
からも共通積分を呼べます。中心からの変位はbohr、経路パラメータは0〜1、重みはその積分に対する重み、
E/Hは経路上の同じ位置にある `[点数,3]` 配列です。戻り値は電子のスカラーポテンシャルエネルギーとA/cです。
経路上の場と積分点・重みは呼び出し元が用意し、積分点数の収束も確認します。この関数だけではTDDFTは進みません。

### 係数を更新ごとに渡す

```python
with IsolatedTDDFTInitialState.from_salmon(
        'GS_BUNDLE', 'rt.inp', runtime=runtime, max_steps=80,
        near_field=True) as initial:
    with initial.new_state(dt=.003) as state:
        saved = state.savepoint()
        frame = state.step_near_field(near.frame(1), near.frame(2))
        state.restore_savepoint(saved)
        state.release_savepoint(saved)
        state.step_near_field(near.frame(1), near.frame(2))
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

更新nでは、SALMONの既存順序に合わせ、場のステップnとn+1を渡します。
各フレームには中心・半径・次数・時刻・ステップ番号・E/H係数が入ります。
`NearFieldFrame.from_samples(..., step_index=n, center=..., radius=...)` で1回分の空間E/Hから
係数を作ることも、`NearFieldFrame` に自作コードの係数を直接渡すこともできます。
隣接更新で共有する場のフレームを修正する場合は、その場を使う前のsavepointへ戻して更新し直します。

将来の場がまだない場合は、最初の `NearFieldFrame` の `write_live_input(path,max_steps=...,dt=...)`
で初期値を保持する入力を準備し、各更新で `step_near_field` を呼びます。
この入力は将来の場の予測ではありません。
通常の `step(E_extra)` は入力済みのE/H系列を使い、引数を追加の一様Eとして扱います。

保存は `savepoint` または `save_direct_checkpoint`、再開は対応する直接復元を使います。
ライブ係数はEだけの再実行履歴では保存できないため、近接場モードのreplay snapshotは提供しません。
係数系列・中心・半径・磁気用配列・非局所位相まで直接保存します。復元時に時間発展を再計算しません。

### 結合と確認範囲

Fortran側では、展開中心からの変位rについて
`V_e = r·integral_0^1 E(s r) ds`、`A/c = -r×integral_0^1 s H(s r) ds / c`
を使います。`V_e` は電子のポテンシャルエネルギーです。
対称化した運動エネルギーの磁気項とA²項を保持し、非局所擬ポテンシャルには
原子中心からの線積分位相を使います。これは[動径方向のゲージによるポテンシャル](https://arxiv.org/abs/physics/0204034)
に基づく、指定された古典場との結合です。

C₂H₂で電気四重極・磁気の効果、入力済み系列とライブ係数の一致、初期状態への復元、
プロセス終了後の直接再開を確認しています。これは短時間のソフトウェア検証です。
実材料の応答、展開次数、FDTD/電子の格子、dt、場の近似残差は目的に応じて確認します。E/Hを独立にフィットしただけでMaxwell方程式との整合性が保証されるわけではありません。
場のないGSから開始し、先頭から非ゼロのE/Hを与えた場合は、その初期状態からの時間発展です。
過去の照射履歴や場ありのGSを自動で作ることはありません。

## 保存した電場の準静的ポテンシャル近似

`salmonpy.near_field.ElectricPotentialSeries` は、同じ位置・時刻で取得した電場を、
実球面調和関数の正則な固体調和関数で近似します。数値処理はFortranです。
**一般の近接場のE/H入力ではなく、スカラーポテンシャルで表せる電気成分の近似です。**
Hの結合、磁気多極子、分子からFDTDへの放射の戻しはこの入口に含みません。
FDTDの格子で失われた空間情報も復元しません。

```python
from salmonpy.near_field import ElectricPotentialSeries

# electric_samples: SampledField('E', times, positions, values, 'E:au')
# positionsはbohr。Eの3成分の位置・時刻を揃えてから渡す。
potential = ElectricPotentialSeries.from_samples(
    electric_samples, center=(0., 0., 0.), radius=8., lmax=2)
print(potential.relative_residual)
potential.write('GS_BUNDLE/salmonpy-ylm.dat')  # 準備だけ。計算は開始しない
```

同じ入射波・位置・時刻の真空参照があれば、`reference=vacuum_samples` で差の電場を近似できます。
差の計算もFortranで行います。相対残差が既定の5%を超えた場合は警告し、値を返します。
電場がほぼゼロの時の相対値だけで判断せず、`electric_field(...)` と元の電場も比較します。
次数を上げても磁場や回転する電場成分の欠落を補えるとは限りません。

実関数は複素Ylmの[NISTの規格化](https://dlmf.nist.gov/14.30#E1)から、正のmは
`sqrt(2)*(-1)^m*Re(Ylm)`、負のmは対応する正のmの虚部を同じ係数で使います。
`l=1..lmax`、各lの中で`m=-l..l`の順です。`lmax=1..8`に対応します。
電子のポテンシャルは `V_e=sum(c_lm*(r/radius)^l*Y_lm)`、電場は `E=grad(V_e)`。
係数はHartreeで、定数ポテンシャルl=0は省きます。領域内に外部電荷がない近似を前提とします。

実行は明示的に別の操作で行います。対応する実行ファイルのYlm検証記録が必要です。

```python
with IsolatedTDDFTInitialState.from_salmon(
        'GS_BUNDLE', 'rt.inp', runtime=runtime, max_steps=80,
        electric_potential=True) as initial:
    with initial.new_state(dt=.003) as state:
        for _ in range(40):
            state.step([0., 0., 0.])  # このEは追加の一様電場
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

係数表の時刻は0から一定dtで、`0..max_steps+1`の行が必要です。
更新番号nではSALMONの既存処理に従いnとn+1の係数を参照します。
dt・入力・実行ファイルが同じ状態へ、通常と同じ`restore_direct_checkpoint`で復元できます。
FDTDの途中の時間窓を使う場合、元の時間ステップ・時刻と分子側の時刻原点の対応を保存してください。
場のないGSからその窓で開始する計算は、前の照射履歴を持つ分子の継続とは異なります。

3次元LD/FDTDの27点から分子へ渡し、80更新のうち40更新でプロセスを終了して直接再開した試験で、
連続実行と完全一致しました。一様場のYlm表現は本体と約4e-12 e*bohr以内で一致しています。
ただし実FDTD試験では電場の近似残差が大きく、近接場の物理精度を検証した結果ではありません。
空間依存E/Hを使う場合は、上記の近接場モードを使って電気・磁気を一緒に渡します。

# ASE Calculator：全エネルギーと本体の力

`salmonpy.ase_calculator.SalmonCalculator` は独立した任意モジュールです。
このモジュールのimportにはASEが必要ですが、SALMONpyの通常のimportやLDなどには不要です。
対応する物理量は `energy`（eV）と `forces`（eV/Å）。応力・free_energyは公開しません。
対応範囲は3次元周期・中性・固定イオンSCF・非スピン分極、PZ／Libxc-PZ／Libxc-PBEとスカラーKY／NC-UPF擬ポテンシャルです。
Libxcには対応したランタイムが必要です。擬ポテンシャルの生成汎関数・元素との整合を確認してください。TB-mBJの力は対象外です。

## NLCC付きUPFの力

**Libxc-PBE＋NLCC付きUPFの力は未検証です。構造最適化にはPZ＋KYを使ってください。**
擬ポテンシャルと汎関数を対応させ、格子・k点・SCFの収束も別途確認してください。

CalculatorはUPFの有効または不明な`core_correction`指定、もしくは`PP_NLCC`要素を検出すると、
汎関数にかかわらず力を無効にします。`atoms.get_forces()`と直接の力要求は
`PropertyNotImplementedError`となり、その要求によるSCFを開始しません。
エネルギー計算は引き続き利用できますが、力は読み出し・キャッシュしません。
不明なフラグやゼロ値だけの`PP_NLCC`要素も保守的に対象とします。

エネルギーを取得できることは、NLCCのエネルギーが検証済みであることを意味しません。
既存Libxcランタイムでは、NLCCを加えた密度と価電子のみのエネルギー重み・GGA勾配の混在が
確認されています。互換パッチにはXC整合性と明示的NLCC力の修正を追加していますが、
既存の実行ファイルがこの修正を含むと推測してはいけません。

このPython保護自体は本体の力を補正しません。検証済みランタイムを識別できる仕組みが
整うまでは、CalculatorのNLCC力の拒否を維持します。
低水準の`GroundStateResult.forces_hartree_per_bohr`から読み出せる値も、
NLCCの場合に検証済みの力とは解釈しないでください。

```python
from ase.build import bulk
from ase.optimize import BFGS
from salmonpy.ase_calculator import SalmonCalculator
from salmonpy.dft import Pseudopotential, SCFSettings
from salmonpy.fdtd import SalmonExecutable

atoms = bulk('Si', 'diamond', a=5.43, cubic=True)
atoms.calc = SalmonCalculator(
    executable=SalmonExecutable(executable_path, validated_sha256),
    directory='ase-calculations',
    pseudopotentials=[Pseudopotential.from_file('Si_rps.dat', atomic_number=14, local_channel=2)],
    electrons=32, states=32, grid=(32, 32, 32), kmesh=(2, 2, 2),
    scf=SCFSettings(max_iterations=500, tolerance=1e-12), threads=2,
)
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
# 格子・k点・SCF精度を目的の力の精度に対して確認してから使用。
# BFGS(atoms, logfile='opt.log', maxstep=0.01).run(fmax=0.05, steps=100)
gs = atoms.calc.last_ground_state
```

`electrons` は省略可能です。各構造の擬ポテンシャル価電子数の合計から決定し、明示値が異なればネイティブ計算前に拒否します。`states`・`grid`・`kmesh` は明示指定です。

上の格子はSiの例であり、任意の系や目標精度への推奨値ではありません。
実行ファイルにはGS収束報告と `gs_forces.dat` の `SALMONPY_FORCES_V1` 報告を含めます。
監査と実行環境の準備はDFT.md、追加依存はASE.mdを参照してください。

## 状態・保存・失敗

同じ構造へのenergy/forces問い合わせはASEのキャッシュを使います。位置・セル・元素・
周期境界・初期電荷/磁気設定が変われば、新しいGS計算を行います。初期電荷・磁気設定に
非ゼロ値があれば現在の非スピン分極入口では拒否します。既定では前のGSを流用しません。`warm_start=True` の場合は下記の条件で波動関数を再利用します。

各計算は指定ディレクトリ内の新規 `gs-...` に保存し、以前の結果を上書きしません。
このため長い最適化では保存量が増えます。結果の保持・削除は利用者が管理します。
同じCalculatorインスタンスを複数スレッドから同時使用することは想定しません。

モデル設定は構築時に固定し、`set(...)` による変更は拒否します。格子などを変更する場合は
新しいCalculatorを作って `atoms.calc` に設定します。未収束・失敗時は結果と
`last_ground_state` を消去し、古い力を返しません。失敗ログは計算ディレクトリに残ります。

`GroundStateResult.forces_hartree_per_bohr` で本体の力を原子単位で読み出せます。
このプロパティには新しい力報告ファイルが必要です。過去の報告なしGSを読み直せても、
そこからこのプロパティを取得できるとは限りません。

## egg-box効果と力の検証

実空間格子では、原子と格子の相対位置によって全エネルギーが揺れることがあります。
そのエネルギーを有限差分した力にも揺れが含まれます。本体の軌道微分による力を
有限差分値に置き換える実装にはしていません。

検証では次の三つを分けます。

1. 本体の力を正しい原子順・符号・単位で渡しているか。
2. 全原子の剛体並進に対するエネルギーの揺れが、格子間隔とともに小さくなるか。
3. 本体の力とエネルギー微分が、SCF精度・格子・変位幅に対してどう収束するか。

今回の変位Siでは16³の半格子並進で約0.27 eV、24³で約0.044 eVのエネルギー差があり、
1格子並進では元のエネルギーに戻りました。粗い格子の数値微分だけで本体の力の良否を
判定しないでください。実際の数値と条件は検証記録を参照してください。

ソース中の `tools/check_force_grid.py` で同じSi診断を再実行できます。
`tools/test_ase_calculator.py` は24³格子・固定セル・一部原子固定のBFGSを1ステップ実行し、
エネルギー低下と接続動作を確認するテストです。構造最適化の収束や物性精度の保証ではありません。

参考： https://pubs.rsc.org/en/content/articlehtml/2015/cp/c5cp00351b

## warm start（明示的に有効化）

`SalmonCalculator(..., warm_start=True)` を指定すると、直前の成功したGSと計算条件が一致する場合、波動関数を新しい原子位置のSCF初期値にします。セル・格子・k点・元素順・擬ポテンシャル・電子数・軌道数・XCが異なる場合は新規SCFです。元ファイルの破損はエラーとして扱い、黙って新規SCFへ切り替えません。失敗後は古いエネルギー・力・初期値候補を返しません。

混合履歴・旧原子位置は引き継がず、反復数をリセットします。状態と出力先は各計算で独立です。`warm_start` は既定Falseで、従来の更新経路を維持します。SCFの同一許容誤差で新規計算とのエネルギー・力を比較してから、最適化に使ってください。高速化は系によって変わります。

結晶TDDFTの基本は **transverse**（`trans_longi='tr'`）。Longitudinalは誘起場を扱う別の明示的な選択で、Calculatorの既定ではありません。ASEの固定セル最適化／Born–Oppenheimer MDと、電子・イオンを同時伝播するEhrenfestは別です。初版の目標には孤立系・non-collinear・Ehrenfestを残し、今回の対応済み範囲と区別します。

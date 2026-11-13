# SALMONpy：Pythonから使うSALMONの計算・解析部品

> Standalone layout: full SALMON sources are no longer bundled. Source-checkout
> runtime setup: `runtime/README.md`; use `tools/prepare_salmon_runtime.py` with
> the pinned external upstream checkout. Historical validation commands below
> describe the old source layout and are not current build instructions.

[English](README.en.md) · Topic: `getting-started` · Paired revision: `2026-11-13.1`
LD・SBEの物質応答、構造・光パルスの作成、解析を単独で利用できます。
対応する外部SALMONを用意すれば、DFT・TDDFT・FDTD・バンド計算にも接続できます。
現在はpreview版で、PyPI公開済みではありません。

## 最初の使い方

配布された、自分のOS・CPU・Pythonに対応する**修復済みwheel**を導入します。
wheel利用ではFortranコンパイラは不要です。import時のコンパイルも行いません。

```sh
python -m pip install "PATH_TO_COMPATIBLE_WHEEL.whl"
python -m salmonpy.doctor
python -m salmonpy.examples.ld_response
```

`PATH_TO_COMPATIBLE_WHEEL.whl` は実際のファイルに置き換えます。`python` は導入先環境の
Pythonです。Python 3.10以上が必要ですが、wheelのPython・OSタグにも適合する必要があります。
診断はライブラリの読み込みと依存パッケージの所在を確認します。インストールや計算は行いません。
追加依存の実動作、外部SALMONの機能対応、計算規模の上限を保証するものではありません。

## 目的から選ぶ

| やりたいこと | 入口 | 必要な準備 |
|---|---|---|
| 材料の誘電率・Eに対するP/J | `salmonpy.LorentzDrude` | wheel |
| 結晶のA→J（SBE） | `salmonpy.sbe` | wheel、実材料にはバンド・遷移行列データ |
| 光パルス・構造・解析・結果読み込み | `fields` / `geometry` / `analysis` / `outputs` | NumPy、描画にはMatplotlib |
| DFTの基底状態・力 | `salmonpy.dft` | 対応外部SALMON・擬ポテンシャル |
| 結晶のA→J（TDDFT） | `salmonpy.tddft` | 対応外部SALMON・整合する基底状態 |
| 光伝播（FDTD） | `salmonpy.fdtd` | 対応外部SALMON |
| バンド計算・描画 | `salmonpy.bands` | 対応外部SALMON・基底状態・Matplotlib |
| ASEの構造・経路・energy/forces | `salmonpy.ase` / `salmonpy.ase_calculator` | 上記に加え任意のASE依存 |

[14本のPythonサンプルと実行方法](examples/README.md)から始めてください。
各機能の対応範囲は [WORKFLOWS.md](WORKFLOWS.md)、[DFT.md](DFT.md)、
[ASE_CALCULATOR.md](ASE_CALCULATOR.md) に記載しています。
wheel内のLD/SBEと、別途準備が必要なSALMON本体は区別してください。

## AIに依頼する場合

まず [AI_USAGE.md](AI_USAGE.md) と `manual.json` を読み、既に回答された条件を再利用して、
不足する「目的・使ってよい総メモリ・CPUコア数」を確認します。

```sh
python -m salmonpy.manual 初めて
python -m salmonpy.manual レーザーによる電子の動き
python -m salmonpy.doctor --json
```

検索で一致しない場合も `task_catalog` と `getting_started` に案内があります。
`sample_files` が実行例、`beginner_intake` が資源上限の確認手順です。
環境診断で検出したCPU数を、そのまま利用許可と解釈しません。
メモリ自動見積もり・ハード上限の強制機能は未実装です。

外部SALMONを確認する場合は、監査済みの記録値を指定します。

```sh
python -m salmonpy.doctor --executable PATH_TO_SALMON --sha256 VALIDATED_SHA256
```

ハッシュ一致はファイルの同一性だけを示します。対応機能・バイナリ互換性や実行成功は
検証しません。準備・OS別検証状況は [PACKAGING.md](PACKAGING.md) を参照してください。

## LDのPython API

```python
import numpy as np
from salmonpy import LorentzDrude

material = LorentzDrude(
    plasma_frequency=0.5,
    strength=[0.8, 0.2],
    resonance=[0.0, 0.3],
    damping=[0.01, 0.02],
)
E = np.zeros((101, 3))
E[:, 0] = 1e-3
result = material.respond(E, dt=0.1)
# result.time: (101,), result.polarization/current: (101, 3)

state = material.new_state(dt=0.1, initial_field=E[0])
P, J = state.step(E[1])  # fields at t=0.1; returns physical P and J
```

All values use SALMON atomic/Gaussian units: angular frequency, damping and plasma
frequency in inverse atomic time, time in atomic time, E/P/J in the corresponding
atomic units. Do not pass eV energies or femtoseconds without conversion. Pole
strengths are dimensionless. Each pole obeys

```
dP/dt = J
dJ/dt = f * wp**2 * E / (4*pi) - gamma*J - w0**2*P
```

`respond` starts with P=J=0 at the first supplied E sample. Sample zero is the
initial condition, not an update. `step` takes E at the next time, using the stored
previous E. Time step is fixed for each state. Every call to `respond` starts a
fresh state; separate `new_state` objects do not share histories. Observed P/J are
summed over poles. They are not the auxiliary FDTD quantity Q.

`permittivity(omega)` evaluates the continuous model with exp(-i omega t):
`1 + sum(f*wp**2 / (w0**2 - omega**2 - i*gamma*omega))`.
The instantaneous background is fixed at 1; anisotropy, active media, fitting,
and unit conversion are not included. An exact real-axis pole raises ValueError.
Finite-step time propagation approximates this continuous response; reduce dt to
check temporal convergence. A single state must not be shared between threads.

The step adapter reconstructs a small native context each call; it is intended
for prototyping, not high-throughput grid coupling. Batch calculations reuse one
native context. FDTD has an external-executable batch wrapper. The built-in SBE adapter is `salmonpy.sbe.SBEInitialState`; see the loader and
independent-trajectory example in WORKFLOWS.md. A process-isolated PZ/TB-mBJ/Libxc-PZ/Libxc-PBE TDDFT adapter is available as `salmonpy.tddft.TDDFTInitialState`
with a separately built worker-enabled SALMON executable.
Crystalline TDDFT/SBE interfaces remain A-driven only.

See [AI_USAGE.md](AI_USAGE.md) for extension rules and the legacy FDTD boundary.

Input ownership, immutable resolved settings and failure behavior are specified in
[INPUT_CONTRACT.md](INPUT_CONTRACT.md). Inspect `material.config` and `state.config`
for all effective settings. Python 3.10 or newer is required.

## Execution environment

Commands use `python` for the interpreter in the active environment (Python 3.10+).
Use that environment's interpreter name or full path if different; do not assume
`python3` exists on Windows. Run repository-only commands from the checkout root;
`build_native`, `tools/` and tests are source-checkout commands. The examples are included in wheels and run with `python -m salmonpy.examples.NAME`.
Replace example executable paths and hashes with the selected platform's audited
SALMON executable (including `.exe` on Windows where applicable). `executable_path`
and `validated_sha256` in examples are caller-supplied values, not detected defaults.
Use new, writable output directories. Relative paths resolve from the process's
working directory, not the installed package directory.

Pure Python analysis needs NumPy; plotting additionally needs Matplotlib. On a
headless machine select its `Agg` backend before plotting (e.g.
`import matplotlib; matplotlib.use('Agg')`). No GUI is needed to save a figure.
Installed-wheel compatibility does not establish compatibility of an external
full SALMON executable. See PACKAGING.md for the version-specific validation matrix.


実行可能な用途別Pythonサンプルと必要な入力は [examples/README.md](examples/README.md) を参照してください。AI検索の各taskの `sample_files` からも辿れます。

## ソースから開発する場合

ソース構築のコンパイラ・BLAS要件は [PACKAGING.md](PACKAGING.md) を参照してください。
リポジトリのルートで明示的に構築します。再構築後はPythonを再起動します。

```sh
python -m salmonpy.build_native
python -m unittest discover -s salmonpy/tests -v
```

## 独立配布と日英ドキュメント

将来は独立リポジトリとして配布し、このブランチをSALMON本体へマージしない方針です。
`python -m salmonpy.manual --lang en "I want to run TDDFT"` と
`python -m salmonpy.doctor --lang en` で英語の開始案内・環境診断を利用できます。
AIカタログの詳細技術説明はまだ言語が混在しており、全リファレンスの翻訳完了ではありません。

公開版SALMON VRへの接続は [SALMON_VR.md](SALMON_VR.md) を参照してください。

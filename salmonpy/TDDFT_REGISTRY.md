# 同一プロセス内の異入力TDDFTモデル

[English](TDDFT_REGISTRY.en.md)

`NativeTDDFTBackend.open_registry(context)` で明示的に1つの呼び出し元プログラムが管理するMPI接続を作り、
`registry.open_model(initial, dt, decomposition=...)` で異なる入力のモデルを登録します。
`model.new_state(point_id)` は従来と同じA/c→電荷電流の状態を返します。
モデルごとに入力、XC、セル・格子・k点・原子・擬ポテンシャル参照、dt・容量・内部MPI分割を保持します。
同じコミュニケータを使うため、分割の積はその群のランク数に合わせます。
モデル内の入力は不変で、各状態は元のGSから独立に始まります。
物質ID・点IDは呼び出し元プログラムが付け、計算処理の不透明なモデル／状態トークンとは区別します。

```python
from salmonpy.embedded_tddft import NativeTDDFTBackend

# 検証済みの初期状態と呼び出し元プログラムのMPI群を用意する。
backend = NativeTDDFTBackend(build_record, work_root,
                            decomposition=(1, 2, 1, 1, 1))
registry = backend.open_registry(context)
models, states = [], []
try:
    a = registry.open_model(initial_a, 0.003)
    models.append(a)
    b = registry.open_model(initial_b, 0.004, decomposition=(1, 1, 1, 1, 2))
    models.append(b)
    sa = a.new_state("film-a/0")
    states.append(sa)
    sb = b.new_state("film-b/0")
    states.append(sb)
    ja = sa.step([1e-4, 0.0, 0.0])
    jb = sb.step([0.0, 1e-4, 0.0])
finally:
    for state in reversed(states):
        state.close()
    for model in reversed(models):
        model.close()
    registry.close()
```

保存した初期状態に適合する分割を選んでください。dtの違いは独立した時刻の例であり、
FDTDのサブサイクルを追加するものではありません。通常の層配置は既定ではIDを別群に置きます。明示的な model_groups 指定で同じ応答群への登録・配置も可能です（下記）。
明示的なregistryを使う呼び出し元プログラムは状態を直接呼べます。

`sp_tddft_registry_protocol()==1` に対応する明示的に構築・検証したCPU計算プログラムを使います。
周期・全k点（縮約なし）・非スピン偏極・固定イオンのPZ／TB-mBJ自動・固定c／Libxc-PZ/PBE、
transverse・middlepointが対象です。Pythonの初期Aは0で、x/y/z・回転するAも扱えます。
Aは直交座標のA/c、dtは原子時間、jは原子単位の電荷電流密度で、再び符号反転しません。
HSEは別のHSEエンジンであり、このregistryには登録しません。
import時のMPI初期化・ビルド・ダウンロードはありません。
入力・作業ファイルを変更せず、MPI群で共有する書込み可能なwork_rootを用意します。
作業コピーは状態ごとではなくモデルごとです。登録時に最初の状態を初期化して保持し、
最初のnew_stateで返します。使わずにモデルを閉じた場合も解放します。

群内の全ランクがメインスレッドから同じ順に呼びます。
registryは明示的な別経路で、legacy初期化や暗黙のmodel_newと混ぜません。
既定の単一モデル入口は維持します。モデル／registryのcloseは新規作成を止めますが、
既存状態は続行できます。全状態を閉じてから呼び出し元プログラムが管理するMPIを終了してください。
一方のモデルを閉じても他方は残ります。待機状態の取り出し前に群内のモデル選択を照合し、
終了順もトークンで統一します。扱える終了失敗では残りの所有者も終了を試し、再試行できます。

モデルは入力・既定値、初期化後の設定ビュー、初期A、最適化表・フラグ、XC作業配列を所有します。
配列はmove_allocで退避・復帰し、波動関数やLibxcハンドル・MPI資源をコピーしません。
状態は波動関数・密度・XC演算子・時刻・A・電流・出力と保護／復元用の別状態を所有します。
擬ポテンシャル形式のフラグは元素の初期化時にリセットします。
タイマー・MPI接続・同期的な互換ビューはプロセス共通の基盤として残ります。
並行スレッドや無制限な別ソルバー初期化ではなく、モデルの順次切替です。
一般の入力・割当・I/O・プロセス・MPI異常からの復帰は保証しません。

状態の `checkpoint(path)` と、新しい状態への `checkpoint(path, restore=True)` も使え、
過去のステップを再計算しません。ファイル名はランクごとに分け、同じバイナリ・入力・モデル条件・
dt・容量・MPI分割とファイルを呼出側で認証してください。単独ファイルは可搬アーカイブではありません。
既存の所有者単位の保護付きstep／復元は同じ正の状態ハンドルを使います。
通常計画からのregistry構築は明示的な model_groups 指定で利用できます。MPI分割を変える可搬再開は別作業です。

`python -m salmonpy.examples.tddft_models BUILD INITIAL_A INITIAL_B NEW_OUTPUT` を適合するMPI数で実行できます。
既定は軌道分割です。開発者用 `tools/test_heterogeneous_tddft_mpi.py NEW_BUILD OLD_BUILD INITIAL_A INITIAL_B OUTPUT --size 2 --b-split rgrid`
は6ランクを3群に分け、最初の群内の2モデルを別群の旧計算プログラム独立実行と比較します。
交互伝播・直接再開・兄弟終了・空モデル再作成・不正／古い／不一致ID・登録失敗・終了順を確認します。
既存Siの16³／8³格子、異なるセル・原子数・k点と5種類のXCを使い、新しいGSは作りません。
これは実装の一致検証であり、物性の収束保証ではありません。
基礎状態とパーサーは [共有TDDFT](SHARED_TDDFT.md)、C入口は [専用ヘッダー](native_adapters/standard_tddft_models.h) を参照してください。

## 通常の物質配置・混在計算から使う

`assemble_native_layers(..., model_groups={"film": ["silicon-pbe", "silicon-tbmbj"]})`
で、指定IDの全マクロ点を一つのMPI物質群へ配置できます。各IDのbindingには
`NativeTDDFTBackend(..., registry_model=True)` を指定します。同群内の標準TDDFTは
同じライブラリと作業ルートを使います。各物質の入力・内部分割は別々で、分割の積と
`ranks_per_point` は群のランク数に一致させます。各マクロ点の状態は独立です。

`ResponseSession` は群に一つのレジストリを作り、物質IDごとにモデルを登録します。
全状態・モデルの後にレジストリを終了します。LDは引き続きFDTD群内で更新します。
結合のdtは共通です。既定の配置は変わりません。同じ物質IDへ `points_per_group` と
`model_groups` を重複指定しないでください。群外のIDには従来の分割を使えます。

通常サンプルでは、既存の `--layers` と `--tddft-materials` に次を追加します。

```text
--tddft-model-groups mixed_tddft_model_groups.json
```

同梱のJSONは `{"silicon-models": ["silicon-pbe", "silicon-tbmbj"]}` です。
物質カタログの各IDへ既存GSバンドルを指定し、配置レポートの `required_ranks` に合わせて
起動してください。計算ノードでは外部取得せず、入力・擬ポテンシャル・計算プログラムを事前搬入します。
物質ID／点IDと実際のMPI配置は出力・一括チェックポイントにも記録されます。
変更したMPI分割での再開には対応しません。

自作の計算コードは従来の `prepare`／`open_group`／`create_state` を利用できます。
必要なら `SolverContext.group_resources` に群所有のclose可能な資源を登録できます。
共有資源は状態・物質factoryの後に、登録順で終了されます。群内の全ランクで登録順と
呼出順を揃えてください。組込みHSEの異入力モデルは別群へ配置します。

## 呼び出し元プログラムの作業場所

directory scope2/3対応runtimeでは、登録・状態生成・step・直接checkpoint・終了の
呼出しが呼び出し元プログラムの作業ディレクトリを変更しません。状態／モデルごとのworkフォルダ
を計算処理のパス解決に渡し、相対ファイルをその基準で扱います。例外やランクの
準備失敗時も範囲を解除します。旧runtime（scope1を含む）は従来の切替を維持します。
メインスレッドで順番に呼ぶ制約は継続し、数値保存の項目や更新順は変えません。

標準登録モデルの状態も、対応するtimer protocol1があれば `state.timing_seconds` を共有状態と同じ仕様で照会できます。[共有TDDFT](SHARED_TDDFT.md)を参照してください。

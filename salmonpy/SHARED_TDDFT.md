# 標準TDDFTのPython共有モデル

[English](SHARED_TDDFT.en.md) · Topic: `shared_standard_tddft` · Revision: `2026-11-18.1`

`NativeTDDFTBackend.open_model(initial, dt, context)` は `NativeTDDFTModel`
を返します。不変のGS・入力・擬ポテンシャルをモデル単位で読み、作業コピーとnative接続も
1つだけ持ちます。`model.new_state(point_id)` ごとに元のGSから独立した電子状態を作ります。
波動関数、密度・XC、A・電流・時刻、出力、保護・復元用状態は各状態が所有します。
point_idは利用者の識別子であり、nativeハンドルや物質定義とは別です。

標準handle/model protocol1対応の明示的に構築・検証したランタイムを使います。
import時のビルド・ダウンロード・MPI初期化・入力取得はありません。
変更のない検証済み実行物は記録済みSHA256で確認します。ホストがMPIを管理し、
物質群内の全ランクがメインスレッドから同じ順に呼びます。入力・作業ファイルは変更せず、
work_rootは初期状態のフォルダ外に用意してください。

基本コードは [English](SHARED_TDDFT.en.md) の実行例と共通です。
`provider.prepare(context)` で初期状態を用意し、`backend.open_model(initial, dt, context)`、
`a=model.new_state("film/0")`、`b=model.new_state("film/1")` と作成します。
各 `state.step(A)` でその点の電流jを得ます。
`a.checkpoint(path)` で保存し、新しい状態の `checkpoint(path, restore=True)` で直接復元します。
保存先親フォルダは事前作成し、MPIランクごとに異なる未使用ファイル名を指定します。
復元には同じバイナリ・入力・dt・容量・MPI分割と認証済みファイルが必要です。
単独nativeファイルの認証はホスト側の責任です。通常サンプルの一括保存なら由来検査も行います。
再開時に過去のステップを再計算せず、過去のテキスト出力の連結は行いません。

状態とモデルは明示的にcloseし、その後ホストのMPIを終了します。
モデルclose後は新しい状態を作れませんが、既存状態は継続できます。
1つの状態のcloseで兄弟状態は終了しません。モデルが開いていれば全状態close後の再作成も可能です。
1プロセスで1つのnative電子モデルの生存期間を扱う契約であり、異なる入力のモデルや並行スレッドを
同じプロセス内に自由に配置する機能ではありません。

対応範囲は周期系・全k点（対称性縮約なし）・非スピン偏極・固定イオン、
PZ、TB-mBJ（c自動／固定）、Libxc-PZ/PBE、transverseのmiddlepoint伝播です。
**dt、ステップ容量、初期A=0、並列分割、不変の物質入力はモデルで共通です。**
異なる物質は別MPI群に配置します。x/y/z方向・回転するAにも対応します。
Aは直交座標のA/c、dtは原子時間、返り値は原子単位の電荷電流密度jです。
handle経由の値は変換済みなので、再び符号反転しないでください。
テキスト出力列の符号はSALMONの既存定義に従います。
状態ごとのRT・電流とエネルギー出力に対応し、共通名を使う任意解析出力は従来の別状態経路を使います。
legacy初期化とhandle初期化は同一プロセスで混ぜられません。
初期化・メモリ確保・I/O・MPI異常一般について新たな復帰保証を設けるものではありません。

## 通常の複数点配置

`NativeTDDFTBackend(..., shared_model=True)` と
`assemble_native_layers(..., points_per_group={"si-tddft": 2})` で明示的に有効にします。
`ranks_per_point["si-tddft"]` は共有する物質群のランク数で、分割の積と一致させます。
物質ID、点ID、並び、厚み、LD領域を維持したまま、同じモデルの電子点を群にまとめます。
既定の1点1群・`create_state` は変わりません。ResponseSessionが `open_group` を呼び、
モデルと返された状態を管理します。

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD SERIAL_RUNTIME INITIAL_BUNDLE NEW_OUTPUT \
  --layers mixed_shared_tddft_layers.json --tddft-ranks 2 \
  --tddft-points-per-group 2 --steps 128 --dt .003 --protect-steps
```

JSONはインストール先 `salmonpy/examples/` からコピーするか絶対パスで指定します。
場1＋SBE1＋共有TDDFT2ランクと2種類のローカルLDを使います。
共有オプションを省く独立実行は6ランクです。
`mixed_shared_tddft_hse_both_layers.json` はHSEも2点置く例です。
適切なHSE引数と `--hse-points-per-group 2` を併用すると、TDDFTとHSEそれぞれで共有できます。
異なる物質間でモデルを共用する意味ではありません。
資源予算は警告で、実行可能な計算を止めません。実際のランク数と計画の不一致はエラーです。

`--direct-save-at 63` と、新しい実行での `--direct-resume CHECKPOINT` により、
同じトポロジーの場・LD・SBE・各TDDFT状態を一括で直接保存／再開します。
保護付き復元は全所有者の候補を準備してから公開します。準備途中の扱える失敗では元の全状態を保持します。
状態の電流メタデータもnative復元候補を作る前に検査します。
詳しくは [結合](COUPLED.md)・[分散実行](DISTRIBUTED.md) を参照してください。

## 検証範囲

`tools/test_python_shared_tddft_mpi.py BUILD INITIAL_BUNDLE NEW_OUTPUT --size 2`
は6ランク・3群で、異なる3成分Aを与えた共有2状態を独立2実行と比較します。
直接保存からの再開、NaN/Inf/複素Aの拒否、兄弟状態・空モデルからの再作成も検査します。
`--rgrid` で軌道分割から空間分割へ切り替えます。

PZ・Libxc-PZ/PBE・TB-mBJ自動／固定について、MPI1・軌道2・空間2で局所検証します。
通常のLD/SBE/TDDFT2点混在例では、独立／共有実行と直接再開、群内の最後のランクでの
復元途中失敗・再試行を比較します。小さい既存入力による所有権・数値一致の検証であり、
物性の収束保証ではありません。新しいGSデータの収集やFortranカーネル変更は行いません。
低水準の規約は [native C契約](STANDARD_HANDLES.md) を参照してください。

## 異なる物質入力：現在の対応と残作業

**別々のMPI物質群なら、異なる入力のTDDFTを同時に使えます。** 各群の中では同じ物質入力の
複数点を共有できます。通常サンプルの `--layers` と `--tddft-materials FILE` を組み合わせ、
物質IDごとにローカルの初期状態、入力ファイル、並列分割、共有点数を指定します。

```json
{
  "schema": "salmonpy.native-tddft-materials/1",
  "materials": {
    "silicon-pbe": {
      "bundle": "./pbe-initial",
      "input_file": "rt.inp",
      "decomposition": [1, 2, 1, 1, 1],
      "points_per_group": 2
    },
    "silicon-tbmbj": {
      "bundle": "./tbmbj-initial",
      "decomposition": [1, 2, 1, 1, 1],
      "points_per_group": 2
    }
  }
}
```

相対パスはカタログの場所から解決します。絶対ローカルパスも使えます。
管理・ログインノードなどでGS・擬ポテンシャル・入力一式を準備し、計算前に転送してください。
ローダーによる通信・ダウンロード・ビルドはありません。
既定は `input_file="rt.inp"`、分割 `[1,1,1,1,1]`、共有点数1です。
点数1は従来の別状態経路、2以上は上記の共有プロファイルを使います。
分割はk・軌道・x・y・zで、積がその物質群のランク数です。
入力が同じでも異なるIDは別モデル・別群に置きます。結合dtは全群共通で、サブサイクルは導入しません。

層配置の `material_id` がカタログ内のIDを選びます。
同梱の `mixed_distinct_tddft_layers.json` はSBE1点、PBE2点、TB-mBJ2点、真空と2種類のLDを配置します。
内部MPI2・共有点数2なら場1＋SBE1＋PBE2＋TB-mBJ2で6ランク、共有点数1なら10ランクです。
`distinct_tddft_materials.json` を入力の近くへコピーして参照を編集します。

```sh
mpiexec -n 6 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD SERIAL_RUNTIME DEFAULT_INITIAL NEW_OUTPUT \
  --layers mixed_distinct_tddft_layers.json --tddft-materials materials.json \
  --steps 128 --dt .003 --protect-steps
```

位置引数 `DEFAULT_INITIAL` は既存の `si-tddft`／HSE用の互換既定値です。
カタログの独自IDだけを使う場合は読みません。カタログの `si-tddft` は既定値を上書きします。
未使用の既定値・パスは再開条件に含めません。
入力ハッシュ、初期状態の識別情報、MPI配置、資源助言は各物質の実際の入力を参照します。
`sbe-point`、`hse-point`、`drude`、`lorentz` は予約済みです。
カタログの未使用項目は残せます。未知の項目、重複キー、不正な分割はnative初期化前に拒否します。
資源の不足見積もりは警告のみです。既存の一括直接保存／再開と保護付き失敗・再試行も使えます。

**同一プロセス内に異なる入力のモデルを複数置く部分は未実装です。**
現在のCハンドルは `rt_context` を切り替えますが、`model_new` は最初に読んだ1モデルを使います。
別モデルIDの登録や第2入力の独立初期化はまだありません。残る工程は次の4段階です。

1. `read_input`、並列初期化、共有XC・プロファイル設定から、不変のモデル入力・既定値・初期化データを分離。
2. 生存中の状態を変えずに別モデルを構築・選択し、step／復元を選択モデルのdt・容量・XC・形状条件へ接続。
3. 共通キャッシュ、I/O・タイマー、初期化子通信以外のMPI資源の所有権とモデル別終了処理を整理。SOC／+U／DCは別の範囲。
4. Pythonのモデルハンドルへ接続し、異なるセル・格子・k点・擬ポテンシャル・XC、終了順、独立比較、直接復元を検証。

波動関数、密度・XCペイロード、時刻、保護・復元候補、初期化時の子通信はすでに状態ごとに所有しています。
これらを最初から作り直す必要はありません。各工程の規模が違うため、項目数だけで完了率は出せません。
今回のPBE／TB-mBJ接続テストは既存のSi入力2種類による実装一致の確認です。
異元素の物性精度や同一プロセス内での異入力モデル対応を保証する試験ではありません。

# 後処理用の出力と実行場所

## プロセス／CPUの分離

TDDFTの状態ごとに専用プロセスと作業ディレクトリを持ちます。同一マシンの
FDTDプロセスと2つのTDDFTプロセスを接続し、通常の2 MPIランクのマルチスケール
実行と比較済みです。OSが各プロセスを実行するコアを選びます。コア固定や別ソケット
への配置を検証したわけではありません。1つの状態に同時にstepを呼ばないでください。
SBE/LDの状態も独立に作成します。同じ状態を複数スレッドで更新する契約はありません。

別ノード用の通信・配置APIは未実装です。将来も波動関数は担当プロセスに置き、
結合には時刻・ステップ・点ID・A/c・Jを交換する構成にします。場の更新は必要な
全点の応答がそろってから進め、遅延した電流を別のステップに混ぜません。
出力も各状態が別のディレクトリへ書く設計とし、複数プロセスによる同一ファイルへの
同時書き込みを避けます。BLAS/OpenMPスレッド数はプロセス数と合わせて設定します。

## 出力の入口

共通A駆動のバッチ計算では、保存先を明示します（省略時は従来どおり保存しません）。

```python
from salmonpy.crystal import respond
from salmonpy.outputs import OutputArchive

result = respond(initial, drive, output_directory='run-001')
archive = OutputArchive.open('run-001')
time = archive.array('time')
A = archive.array('A')             # SALMON Ac=A/c、原子単位
j = archive.array('current')       # 電荷電流密度 j（SALMONの物質電流は -j）
```

`time`, `A`, `current`は初期状態を含む端点サンプルです（2026-11-14 以前のアーカイブでは `J` = 物質電流）。TDDFTでは、正常にプロセスを
終了してファイルをフラッシュした後、ワーカーファイルも保存します。
手動ステップの場合は `state.finish('point-001')` で保存してその状態を閉じます。
単なる `state.close()` は保存しません。初期状態ローダーは別の計算に再利用できます。

保存先は新規ディレクトリでなければなりません。保存に失敗すると `OutputSaveError`
が発生し、`error.recovery_directory` に元のワーカーファイルを保持します。
その中の `run/` が生データです。可能であれば `accepted-response.npz` と
`recovery.json` も書きます。共通respondの後始末やオブジェクト破棄でも保持されます。
回収したら利用者が削除してください。これは一時領域なのでOSの清掃までは保証しません。

アーカイブは再開チェックポイントではありません。変更されていない
`data_for_restart/` の大きなGSデータは複製せず、ハッシュを記録します。
その他のファイルは保持します。初期ディレクトリに以前の出力が含まれている場合も
捨てず、索引の `content_matches_initial_bundle` で元ファイルと同じ内容かを示します。
入力と出力の混在を避けるため、通常はGSデータと入力だけの初期ディレクトリを使います。

## TDDFTで選べる出力

初期状態のSALMON入力の `&analysis` で選択します。この設定と保存間隔は
状態作成後に変更せず、別条件には新しい初期状態を用意します。

| 物理量 | 有効化・間隔 | 後処理 |
|---|---|---|
| A・電子電流・時刻 | 共通APIで常時記録 | `archive.array(...)` |
| 全エネルギーと内訳 | 本体の `_rt_energy.data` | `archive.table(...)` |
| 電子密度・初期密度との差 | `yn_out_dns_rt='y'`, `out_dns_rt_step` | `archive.cube(...)` |
| 微視的電流3成分 | `yn_out_micro_je='y'`, `out_micro_je_step` | CUBE |
| ELF | `yn_out_elf_rt='y'`, `out_elf_rt_step` | CUBE |
| 静電場3成分 | `yn_out_estatic_rt='y'`, `out_estatic_rt_step` | CUBE |
| k点・バンド・スピン別電流 | `yn_out_current_decomposed='y'`, `out_current_decomposed_step` | 専用ブロック読み込み |
| 入力・ログ・実行ファイルハッシュ・初期ファイルハッシュ | アーカイブ保存時 | `input.inp`, `manifest.json`, `native/` |

例（空間分布は8ステップごと、形式はCUBE）:

```text
&analysis
 format_voxel_data='cube'
 out_rt_energy_step=8
 yn_out_dns_rt='y', out_dns_rt_step=8
 yn_out_micro_je='y', out_micro_je_step=8
 yn_out_elf_rt='y', out_elf_rt_step=8
 yn_out_estatic_rt='y', out_estatic_rt_step=8
 yn_out_current_decomposed='y', out_current_decomposed_step=8
/
```

`out_rt_energy_step` は本体の固有エネルギー再計算の間隔で、エネルギーファイルを
その間隔だけに間引く指定ではありません。TB-mBJの出力エネルギーを変分的な
TB-mBJエネルギー汎関数の検証値と解釈しないでください。
静電場は本体のVh勾配の出力で、全Maxwell場そのものではありません。
微視的電流とセル平均Jでは非局所擬ポテンシャルの扱いも含め定義を確認してください。

```python
energy = archive.table('native/Si_rt_energy.data')
density = archive.cube('native/Si_dns_000008.cube')
rho = density.values              # (nx, ny, nz)
origin, axes = density.origin, density.axes
# 座標 = origin + [ix, iy, iz] @ axes
print(density.unit_note)          # 本体の単位コメントを保持。自動変換しない

bands = archive.decomposed_current('native/Si_current_decomposed.data')
t_band = bands.step * archive.metadata['dt']
J_total = bands.current.sum(axis=1)
# indices の列は k点・バンド・スピン。占有数・k点重みは本体ですでに乗算済み。
```

時刻・単位・符号を推測して書き換えません。ネイティブ表の列情報と元ヘッダーも
索引に残します。CUBEはSALMONの正の格子数を持つスカラー形式に対応します。
AVS/VTKは生ファイルを保持できますが、このモジュールに専用リーダーはありません。

RT波動関数チェックポイント、射影・unfoldingの出力は現在のワーカーでは未対応です。
`write_rt_wfn_k='y'` は明示的に拒否します。終了時の本体の自動FFT・pulse出力も
呼ばないため、保存したA/Jから解析します。収束済みGSの波動関数をRT最終波動関数と
誤って扱わないでください。

## LD・SBE・FDTDの選択保存

`save_arrays` は実数・複素数の選択配列を保存できます。単位と条件は呼び出し側が
明示し、ライブラリは物理的意味を推測しません。例えばLD:

```python
from salmonpy.outputs import save_arrays
save_arrays('ld-run', {'time':response.time, 'E':E,
                     'P':response.polarization, 'J':response.current},
            metadata={'units':{'time':'a.u.','E':'a.u.','P':'a.u.','J':'a.u.'},
                      'description':'LD response; model parameters recorded separately'})
```

- LD: E・P・J、周波数と複素誘電率を選択保存できます。
- SBE: 共通A/J時系列に加え、手動stepで選んだ時刻の `state.density` をコピーして
  密度行列・占有・コヒーレンス解析に使えます。全時刻の密度行列は自動保存しません。
- FDTD: バッチ実行の観測点E/H・材料格子・入力は既に実行ディレクトリに残ります。
  NumPy配列に整理して同じ保存関数で出せます。2点マルチスケール検証用のA/J履歴は
  検証成果物に保存済みですが、公開E/H単一ステップAPIはまだありません。
- DFT: 通常SALMONのエネルギー・密度・固有値・k点・遷移行列などを元の出力として
  保持します。GS計算・保存・再読み込みはDFT.mdに対応します。全形式の統一リーダーではありません。
- スペクトル: `salmonpy.analysis.spectrum` の複素結果と `dict(result.config)` を
  一緒に保存すると窓・フーリエ符号・規格化・単位を残せます。

容量の大きい分布・行列は保存間隔と対象を選び、後から復元できない情報を優先します。
ハッシュ検証は `OutputArchive.open` で既定有効です。利用者が意図的に編集した
アーカイブを読む場合のみ `verify=False` を指定できます。

## アーカイブの整合性検証

`OutputArchive.open` は `arrays.npz` の一意な索引登録を必須とし、配列名・形状・dtype・
有限値を実データと照合します。`verify=False` はファイルのハッシュ照合だけを省略し、
構造検証は省略しません。検証時には全配列を順番に読み込むため、大きいアーカイブでは
読み込み時間と最大配列に応じたメモリが必要です。
これは読み込み時の整合性確認であり、電子署名や読み込み後の外部編集を監視する機能ではありません。

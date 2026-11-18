# 標準TDDFTの同一モデルnativeハンドル

**更新:** Pythonの共有モデルと通常の複数点配置への接続も実装済みです。
[Python共有モデル](SHARED_TDDFT.md) を参照してください。既定の別状態経路は維持します。

[English](STANDARD_HANDLES.en.md) · Topic: `standard_tddft_handles` · Revision: `2026-11-18.1`

明示的にビルドした標準host-MPIランタイムが `sp_tddft_handle_protocol()==1` と
`sp_tddft_model_protocol()==1` を提供します。対応ヘッダーは
[native_adapters/standard_tddft_handles.h](native_adapters/standard_tddft_handles.h)です。
このnative C／C++入口を使うPython共有モデルと通常の複数点配置も利用できます。
明示的な共有指定がない場合は、既存の `create_state` と一群一点の経路を維持します。

同一プロセス／モデルでは、物質・セル・全k点・擬ポテンシャル・XC入力を共通不変にします。
周期・対称性なし・非スピン分極・固定イオンのPZ／TB-mBJ／Libxc-PZ/PBE、中点伝播・transverseが対象です。
**dt・容量・MPI分割・初期Aも共通**です。異種入力・対称群、状態別dt／容量、並行スレッドや
他のnativeソルバー初期化を同一プロセスで混ぜる対応ではありません。共有・検査済みGS／PP／入力を
新規作業ディレクトリへ用意し、接続中は変更しません。所有する出力はRT／電流とエネルギーの2本です。
共通名を使う追加解析出力は黙って省略・上書きせず拒否します。従来の単一状態の対応範囲は維持します。

| 入口（`sp_tddft_`を付ける） | 意味 |
|---|---|
| `handle_init(input,dt,capacity,A0,decomposition,&h,j)` | 一度だけ入力を読み、最初の状態を作る。旧initとの併用はしない。 |
| `model_new(&h,j)` | 元GSから独立状態を作る。全状態終了後も再生成でき、状態別のGSコピーは作らない。 |
| `handle_step(h,Aold,Anew,j)` | 従来カーネルの更新順で進める。Aoldはその状態の受理済み値と一致させる。 |
| `handle_info(h,&step,A,j)` | 受理済み時刻番号・A・電流を読む。時刻はstep×モデルdt。 |
| `handle_checkpoint(h,path,mode)` | 0保存・1新規状態へ直接読込・2候補準備・3確定・4取消。各ランク別ファイルをホストが検査する。 |
| `handle_transaction(h,mode)` | 0退避・1数値巻戻しと所有ストリーム終了・3再開・2退避解放。1と3の間にホストが所有ファイルだけを戻す。 |
| `handle_label(h,buffer,capacity)` | 出力prefix。共通base_directoryと `_rt.data`／`_rt_energy.data` を組み合わせる。末尾nullを含む容量が必要。 |
| `handle_close(h)` | その状態・復元候補・保護データ・初期化MPI子資源を解放。兄弟状態とホストMPIは維持する。 |

ハンドルが返すのは原子単位の**電荷電流密度j**（SALMON物質電流の負）です。
テキスト出力は従来のSALMONの規約を維持するので列定義を確認します。
AはCartesian A/c、dtは原子時間です。**旧** `sp_tddft_init/step` は物質電流Jのままで、
Pythonが従来どおり符号変換します。ハンドルの戻り値を二重に反転しないでください。
全k点のこの構成ではx／y／z・回転する外場を渡せます。

ホストがMPIを初期化し、`sp_attach(MPI_Comm_c2f(comm))` で接続します。群の全ランクが
同順序でメインスレッドから呼びます。tokenは正の不透明なint64で、物質IDやアドレスではありません。
再利用時は世代番号を更新し、古いtoken・ランク間で異なるtokenは数値状態を選ぶ前に拒否します。
label容量不足も群内で同じ結果を返します。全ハンドルを閉じてから `sp_detach` し、MPI終了はホストが行います。
一プロセス一接続寿命の契約です。importでコンパイル・取得はしません。

波動関数・密度・A／電流／時計・XC payload・保護用数値／出力情報・復元候補・保存配列の検査上限は
各contextが所有します。逐次呼出し用の選択ポインターは所有オブジェクトをコピーしません。
生成時だけ固有の出力名を選び、毎ステップで物質設定を入れ替えません。外側タイマーはモデルが所有し、
最後の状態終了で止め、空のモデルからの再生成時に再開します。通信要求・XC終了後に所有MPI子を解放し、
管理表は総生成回数でなく同時生存数に応じた枠を使います。

直接再開は同じバイナリ・不変入力・dt／容量・MPI分割で、新規状態へ行います。
受理済みステップは再計算しません。数値復元だけでは過去ログを統合しません。
同一性・ハッシュ検査と所有ファイルの退避はホスト側です。汎用保存形式ではありません。
数値・保護処理の既存の拒否保証は維持しますが、初期化・確保・I/O・MPI・プロセス障害全般の
再試行保証をハンドルだけで追加するものではありません。

開発用検査はwheelの `salmonpy/_runtime_source/tools/test_standard_tddft_handles_mpi.py`、
またはソースの `tools/` から実行します。

```sh
mpiexec -n 6 python tools/test_standard_tddft_handles_mpi.py \
  NEW_BUILD OLD_BUILD GS NEW_OUTPUT --group-size 2 --split orbital
```

3群で共有2状態と旧版の独立2状態を比較します。管理表の拡張、8／9ステップ保存からの新規復元、
同時保護／候補準備、兄弟終了、古いtokenと空モデルからの再生成も検査します。
物理的収束の保証ではありません。既存Python入口は[結合ガイド](COUPLED.md)、HSE固有の
Python共有モデルは[共有HSE](SHARED_HSE.md)を参照してください。

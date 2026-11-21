# 高速PR版HSEの計算処理ハンドル

[English](NATIVE_HANDLES.en.md)

明示的にビルドした呼び出し元プログラムが管理するMPI高速PR計算プログラムは `sp_fast_hse_handle_protocol()==1` と
`native_adapters/fast_hse_handles.h` のC／C++入口を提供します。同じMPI群・作業ディレクトリで、
**一つの不変な物質モデル**から独立した複数の電子状態を作る範囲です。
Pythonの旧 `create_state` は1プロセスにつき初期化1回のままです。
[open_model](SHARED_HSE.md) はモデルを共有し、状態別leaseで接続を保持します。
通常の混在サンプルのランク配置や、同一プロセス内の状態数を自動変更する機能ではありません。

呼び出し元プログラムがMPIを初期化し、`sp_attach(MPI_Comm_c2f(comm))` で接続します。
コミュニケータは最後まで維持してください。protocol照会以外は、全ランクが同じ順序で
メインスレッドから逐次呼びます。呼び出しはMPIの初期化・終了、作業ディレクトリ変更、
計算プログラムのコンパイルを行いません。

| 操作 | 意味 |
|---|---|
| `handle_init(path,dt,capacity,symmetry,&h,j)` | プロセスで一度だけ物質を初期化し、最初の状態を返す。旧singleton初期化との併用はしない。 |
| `model_new(dt,capacity,&h,j)` | model protocol1に対応する初期化済みモデルから、全状態終了後も再生成する。隠れたテンプレートは持たない。 |
| `handle_new(model,dt,capacity,&h,j)` | 生存するモデルの元GSから独立状態を作る。呼び出し元の励起・時間発展履歴はコピーしない。 |
| `handle_step(h,Aold,Anew,j,&energy)` | 選択した状態で従来の高速PR端点更新を実行する。 |
| `handle_info(h,&step,A,j)`、`handle_history(h,values)` | 受理済み端点と、ランク内の8要素の履歴診断値を返す。 |
| `handle_transaction(h,mode)` | 0退避、1数値巻き戻し・ログ閉鎖、2再接続・flush・退避解放、3再接続。再接続失敗時は所有情報を保持する。再接続前に呼び出し元プログラムがログ内容を戻す。 |
| `handle_checkpoint(h,path,mode)` | 0保存、1新規状態への直接読込、2候補作成、3確定、4取消。ランク別pathは異なってよい。従来の同一性・MPI分割・ファイル検証条件を守る。 |
| `handle_label(h,buffer,capacity)` | 呼び出し元プログラムのログ退避用の出力basenameをNULL終端で返す。`_rt.data`／`_rt_energy.data` を付ける。 |
| `handle_close(h)` | その数値・FFT・XC状態を解放する。他状態や呼び出し元プログラムが管理するMPI接続は終了しない。 |

表の名前には全て `sp_fast_hse_` を付けます。ハンドルは正の不透明なint64値であり、
アドレスや利用者の物質IDではありません。世代番号により、同じスロットが再利用されても
終了済みハンドルは無効のままです。管理領域は必要に応じて拡張し、確保エラーはMPI群内で揃えて返します。
statusは0成功、11ランク間の指定不一致、12無効・終了済みハンドルまたはMPI接続不一致、
13ラベル用バッファ不足、16管理領域の確保失敗です。他は従来のHSEの意味を維持します。
初期dt／capacityの不正・不一致は物質初期化の寿命消費前に拒否し、修正後に再試行できます。

dtは原子時間、AはCartesian **A/cの原子単位**、jは**電荷電流密度の原子単位**
（SALMON物質電流の負）、energyはHartree単位の全エネルギーです。
Aoldは受理済み端点と厳密に一致させます。初期化・復元後はinfoで照会してください。
x／y／zと回転・円偏光の入力を扱います。縮約時は実際の外場と保持対称性の整合が必要です。
HSEは高速PR版だけです。同じ対称性の複数状態にはshared_symmetry_protocol1が必要です。異なる対称群の混在には追加の分離が必要です。

[examples/native_fast_hse_handles.c](examples/native_fast_hse_handles.c) を、計算プログラムと
同じMPI・ツールチェーンで明示的にコンパイルします。

```sh
mpicc -std=c11 -I/path/to/salmonpy/native_adapters \
  /path/to/salmonpy/examples/native_fast_hse_handles.c \
  /path/to/libsalmonpy_fast_hse.so -lm -o native_handles
mpiexec -n 2 ./native_handles rt.inp
```

macOSは `.dylib` を指定します。指定MPI数用に準備した、初期外場ゼロ・全k点のRT／GS一式を
**新しい作業ディレクトリ**に置いて実行してください。例はdt・偏光が異なる2状態を進め、
一方を終了し、生存状態のモデルから新状態を作り、古いハンドルが無効のままであることを確認します。
部品接続の例であり、実材料の物理的収束を保証する計算ではありません。

物質入力・対称性・初期化・I/O・タイマー・MPIの所有には共通部分が残ります。
context_mpi_protocol1では初期化時の子を状態終了時に解放し、共通資源は `sp_detach` まで保持します。無制限な生成・終了の
繰り返しを保証する入口ではありません。全ハンドルを閉じてからdetachし、その後呼び出し元プログラムがMPIを終了します。
異なる物質モデル・任意のコミュニケータ変更・同時並行／スレッド呼び出しは対象外です。
ハンドルの安定性と保存形式の可搬性は別です。既存の未使用ワーク領域は独立初期化間で正規化されていません。
復元検証では受理済みの観測量・交換履歴・復元後の将来計算を比較し、独立状態間の保存ファイル全体の
バイト一致を物理的要件にはしません。

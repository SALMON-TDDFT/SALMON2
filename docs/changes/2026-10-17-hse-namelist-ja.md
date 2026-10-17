# HSE機能選択のnamelist化

統合commit d7084f4db を基点に、HSE/ACE/SRの実行時環境変数を
`&functional`へ移行した。新入力と旧環境変数の対応は
[機能選択の説明](../inputs/hse-feature-selection.md)を参照。
入力の既定値、MPI配布、大文字小文字の正規化、variables.log、
半径の単位換算と組合せ検査を追加。旧環境変数は読まない。
checkpointのモード番号と保存形式は維持する。

ユーザー指定によりSR wave batchは16固定。
LCFOキャッシュの固定容量制限、容量指定、超過時停止を削除した。
必要メモリ量は内部保持し、定例の容量出力は行わない。
実メモリ確保失敗とI/O不整合の検出は残す。

## 検証

- MPI/HSE/LibXC/ScaLAPACK、serial HSE、serial non-HSEの3ビルドが成功。
- 新namelist入力と旧環境変数無視を確認。変更前のbinaryでは新入力を拒否。
- 既存422/423/424の計9テストが合格。422で新入力の正規化とログを検査。
- 通常・シフト付きdense k交換oracleと2/4 MPI学習履歴が合格。
- 実入力でstrict/learnedの32step、SR direct/haloの2 MPI・4step、
  dense k再利用on/offの2 SCFを実行。学習予測開始も確認。
- SRのペア再利用有無・halo経路の出力最大差4.75e-21、
  dense k固有値最大差9.82e-15。不正な4組合せを入力段階で拒否。
- 読取りのみの独立レビューで具体的不具合なし。

ログ・入力・source SHAは `work/namelist-hse-20261017` に保存。
カーネル試験と既存回帰ログは `work/coding-rules-e798-20261017` に保存。
長時間計算・性能測定・pushは実施していない。
Fujitsu、GPU、EigenExaは未検証。restart試験をtestsuiteへ追加していない。

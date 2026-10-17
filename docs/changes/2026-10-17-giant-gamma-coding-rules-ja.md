# e798eef4a の Coding Rules 確認と修正

対象は e798eef4a00688cb8914f5512428d78faa281f34。巨大 Gamma 点の
共役ペア再利用、batch 8/16 の計算式、選択条件、許容誤差は変更しない。

## 修正

- 対象コミットの変更ファイルで不足していた各手続きの `implicit none` を追加。
- 132 文字を超える行を継続行に分割し、依存モジュールの import を必要な範囲に限定。
- 対応する既存 wrapper がある通信と時刻取得を SALMON の抽象化へ統一。
  未対応の整数配列通信、バッファ形状、整数種別の MPI 呼び出しは保持。
- 有限値検査をスカラー IEEE 検査に統一。
- MPI 専用 SR 実験経路を `USE_MPI` で保護。非 MPI では出力を初期化して
  未対応 status を返し、通常の serial HSE ビルドを維持。

## 検証

- GNU Fortran Release の MPI/HSE/LibXC/ScaLAPACK ビルド成功。
- 非 MPI の HSE 有効、HSE 無効の両ビルド成功。
- bounds check 付き SR モジュールで batch16 の 8 MPI 試験、ペア再利用の
  8/12 MPI 試験を完走。batch8/16、再利用 on/off、空 owner、非均等軌道分割、
  不正入力の既存 work 内 probe を使用。
- 既存 422/423/424 の preparation/run/verification 計 9 件が合格。
- 対象 Fortran の明示的 `implicit none` と行長検査、`git diff --check` が合格。

ログ、設定、source SHA は `work/coding-rules-e798-20261017` に保存。
新規 testsuite は追加していない。Fujitsu コンパイラ、GPU、EigenExa は今回未検証。
元 snapshot の provenance は歴史的記録として保持。

Distribution notes recorded 2026-11-24; current instructions are in SI_INTRO.md / SI_INTRO.en.md.

# Si入門例：基底状態から電流へ

[English](SI_INTRO.en.md)

固定Si入力でSCFを計算し、A/cで駆動した21時刻の電荷電流を期待値と比べます。
2026-11-23に、現行ソースから作った電子数付き応答の計算プログラムで再検証しました。
入力・擬ポテンシャル・期待値・許容誤差は以前と同じです。

## 必要な配布物

- OS・CPU・Python版に合う、更新済みSALMONpy wheel。Python依存はNumPyです。
- `provenance.json` の現在の `runtime_artifacts` に登録された外側ZIP。

現時点で新しく検証したのはApple Silicon Macで作成したローカル配布物です。公開Releaseにはまだ置いていません。配布者からwheelと対応ZIPを受け取ってください。保存先URLを推測しないでください。新ZIPがない環境では、まずwheelだけの [LD/SBEテスト](README.md)を使います。

旧 `runtime-2026.11.14` Releaseは現在の電子数付き応答と合いません。過去のハッシュと出典は `historical_runtime_artifacts` 等に保持し、現在の入口では受け付けません。Intel Mac・Linux用の新しい配布物については再検証が必要です。WindowsのDFT/TDDFT配布を保証しません。


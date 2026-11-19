# Solver source / 計算本体ソース

## 日本語

現行ソースは `runtime/sources/solver` の一系統です。高速PR版を基盤として、
通常XC・高速HSE・SALMONpyの状態分離と出力処理を統合しています。
共有LD/SBEカーネルはnative側とのバイト一致を検査します。
`source.json` と `source-trees/solver.json` が出典と全ファイルのSHA256を保持します。

通常の入口は `python -m salmonpy.runtime build NEW_BUILD --jobs N` です。
既定はUSE_HSE=OFF。高速HSEは必要なFFTW等を準備して明示的に有効化します。
`tools/prepare_salmon_runtime.py NEW_SOURCE` と旧HSE準備入口は同じファイルを
新しい外部ディレクトリへ検証コピーします。import時にはビルドしません。
禁止固有値ソルバーのソース・object監査を維持します。

旧standard・高速版基準ソース・overridesは `source-history/20261119` の
ハッシュ付きアーカイブとGit履歴に保存しました。現行ビルドは参照しません。
大きな履歴アーカイブはwheel/sdistに含めません。`--reference-only` は
開発チェックアウトでの明示的比較用です。過去の失敗バイナリは使用しません。
既存ランタイムの検証記録やNLCC力の認証は、この配置変更では更新されません。

## English

The single maintained tree is `runtime/sources/solver`, based on fast PR HSE,
with ordinary XC, state ownership and output changes integrated directly.
Shared LD/SBE kernels are checked byte-for-byte against native sources.
`source.json` and `source-trees/solver.json` retain provenance and file hashes.

Use `python -m salmonpy.runtime build NEW_BUILD --jobs N` for explicit builds.
The default is USE_HSE=OFF; fast HSE is explicitly enabled with FFTW dependencies.
`tools/prepare_salmon_runtime.py NEW_SOURCE` and the legacy HSE preparation entry
copy the same verified inventory to a fresh external directory. Import never builds.
Source and application-object prohibition audits remain mandatory.

Frozen standard, fast reference and overrides are preserved as hash-checked archives
in `source-history/20261119` and Git history. Current builds do not use them.
History archives are excluded from wheel/sdist; reference-only preparation is an
explicit developer-checkout comparison. Failed binaries remain evidence only.
Existing runtime validation receipts and NLCC force certification are unchanged.

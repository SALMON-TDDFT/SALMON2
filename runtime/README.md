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

## FFTW and offline fast-HSE builds / オフラインでの高速HSEビルド

FFTW is a dependency of the separately built fast-PR HSE runtime, not of the
lightweight LD/SBE wheel itself. The explicit HSE source build first searches
`FFTW_INSTALLDIR`, pkg-config and system paths for `fftw3.f03` and FFTW3. If no
usable installation is found, the current builder downloads/builds pinned
FFTW3.3.10. That fallback requires network access; it does not happen on import.
On offline nodes, provision FFTW beforehand and pass
`--cmake-arg=-DFFTW_INSTALLDIR=YOUR_FFTW_PREFIX`; check the configure log for
`HSE: using installed FFTW` before proceeding. Do not silently replace HSE with
an ordinary-XC runtime. Ordinary `USE_HSE=OFF` builds do not add this dependency.

FFTWは別途構築する高速PR版HSEランタイムの依存で、軽量LD/SBE wheel自体の依存ではありません。
明示的なHSEソースビルドは、まず`FFTW_INSTALLDIR`、pkg-config、システム内の`fftw3.f03`とFFTW3を探します。
利用可能なものがない場合は、現行ビルダーが固定版FFTW3.3.10をダウンロード・構築するため通信が必要です。
import時には行いません。通信できないノードでは事前にFFTWを用意し、
`--cmake-arg=-DFFTW_INSTALLDIR=YOUR_FFTW_PREFIX`を渡して、設定ログの
`HSE: using installed FFTW`を確認してください。HSEを通常XCへ黙って置き換えません。
通常の`USE_HSE=OFF`ビルドには、この依存を追加しません。

## Vendored Git inventory / 同梱ソースのGit登録

Keep bundled `.gitignore` bytes and source manifests unchanged when repairing
Git omissions. Root ignore exceptions do not override nested upstream ignore
files. After verifying bytes against the manifest, explicitly stage newly
imported manifest paths with `git add -f -- <exact paths>`. Never force-add a
whole working directory or generated build outputs. The source inventory test
checks Git tracking in development checkouts; also prepare/build a Git-only
export after a source-tree import so locally ignored files cannot mask omissions.

Git登録漏れの補修では、同梱`.gitignore`や目録のバイト列を変更しません。ルートの除外例外だけでは、
内側の元の除外規則に勝てません。目録とハッシュ一致を確認してから、新しく導入した目録上のパスだけを
`git add -f -- <対象パス>`で明示登録します。作業ディレクトリ全体や生成物を強制登録しないでください。
開発チェックアウトのソース目録テストはGitの追跡状態も検査します。ソース統合後はGitだけから取り出した
ツリーでも準備・ビルドし、ローカルの未追跡ファイルで欠落が隠れないことを確認します。

## Active GS implementation locations / 現行GS実装の所在

The maintained fast core defines the GS context in `src/gs/main_dft.f90` and
initial seed data in `src/gs/gs_seed.f90`. The old `embedded_gs`/`gs_model`/
`gs_storage` filenames describe the former tree; their inactive `#else` references
do not identify supported fallback implementations. `src/config.h.in` selects
`SALMONPY_FAST_CORE`. HSE interface adapters additionally live in
`salmonpy/native_adapters`, outside the solver-tree inventory. Do not infer a
missing active module from an old filename alone; verify the prepared build.

現行高速コアのGSコンテキストは`src/gs/main_dft.f90`、初期状態データは`src/gs/gs_seed.f90`にあります。
旧`embedded_gs`／`gs_model`／`gs_storage`というファイル名は旧構成のもので、残る非アクティブな`#else`参照は
対応済みの代替実装を意味しません。`src/config.h.in`が`SALMONPY_FAST_CORE`を選択します。
HSEインターフェースのアダプターは別途`salmonpy/native_adapters`にあり、solverツリー目録の外側です。
旧ファイル名だけから欠落と判断せず、準備後のビルド経路を確認してください。

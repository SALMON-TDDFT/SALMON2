# ハイブリッド実装の規約対応・追補

2026-09-29。develop-2.0.0のCODING_RULES.mdと照合した修正。

## 入力の変更

`hse_lcfo_wf_radius`は入力の長さ単位に従う。

| unit_system | 入力半径の単位 | 10 bohrと同じ半径 |
|---|---|---|
| `a.u.` | bohr | `hse_lcfo_wf_radius=10d0` |
| `A_eV_fs` | Å | `hse_lcfo_wf_radius=5.29177210903d0` |

読込み・broadcast後に内部bohrへ変換し、`variables.log`にもbohrで出力する。
半径0は全支持のまま。旧開発版では`A_eV_fs`でもこの変数だけ常にbohrだったため、
その形式の既存入力は半径をÅへ換算してから再利用する。

`rvv10_fft`と`hse_fft_layout`を既存の`string_lowercase`で正規化し、
`FFTW`、`AUTO`なども受理する。

## 構成と検証

- EXX・LCFO・rVV10の対象手続き・インターフェースへ明示的な`implicit none`を補う。
- LCFO実時間接続はdefault privateと明示publicを用い、依存を`only`で限定する。
- 単一手続きだけに必要なimportを局所化する。
- 既存の通信ラッパーが対応する集約・broadcast・ランク照会を共用する。
  対応するラッパーのない可変長転送やScaLAPACK自体は変更しない。
- 209文字のimport行を132文字以下の継続行に分割する。
- PBEh(40)、PBEh(40)+rVV10のGS→RT/Ehrenfest代表試験を番号付きCTestへ登録する。
- 富岳コンパイラ診断ツールのソース名・手続き名をEXX改名へ追従させる。

計算用キャッシュの生存期間や交換の数式、分散方法は変更しない。
実行中のスペクトル計算用バイナリには触れず、検証用ビルドを使用する。

## 検証記録

修正前に追加した入力試験2件が失敗することを確認した。
Å入力が未変換で残ること、および大文字のFFT選択が拒否されることを再現した。

最終結果:

- GNU Fortran 15、MPI/HSE/Libxc/ScaLAPACK有効の統合ビルド：成功。
- GNU Fortran 15、HSE有効・MPI/ScaLAPACK無効の統合ビルド：成功。
- 422–430の標準GS/RT/Ehrenfest（準備・計算・検証）：27/27通過。
  425/427はGamma・nspin=1・LCFO nstate=4に対応する4固有値を検証する。
  初回検証で二k点条件から流用した8固有値の誤った期待数が見つかり修正した。
- 単位・大文字入力、source ACE、fallback、PBE予備収束の回帰試験：22/22通過。
- MPI無効バイナリでも単位・大文字入力の2/2試験が通過。
- LCFO分散構築：MPI2/LAPACKとMPI4/ScaLAPACKで通過。
  129軌道ケースも含み、参照との誤差は約6e-16。
- Streamed seed：MPI2/4およびserialで通過。
  分散QR seed：MPI1/2/3/4で通過。非対応backendの拒否も確認。
- 本番serial通信を用いたtransport/cadence/sphereの単体試験：通過。
- 診断ツール：GNU compile-only全13ケース成功。
- 公式マニュアルの変更節：docutils構文警告なし。
- 独立レビュー：入力・MPI引数・公開範囲・試験依存に残る指摘なし。

検証ログは `work/rules-followup-*.log`。
富岳frtpx、アクセラレータ、他コンパイラは今回未実行。
診断ツールのGNU確認はfrtpx本体の動作確認を意味しない。

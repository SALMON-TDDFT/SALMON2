# 6edccd13e時点のCoding Rules再監査

e798eef4a以後の一連の変更をCODING_RULES.mdと照合した。
無関係なlegacy実装の全面整形は行わない。

## 発見・補正

- symmetry_validate_atoms_cartesianに明示的implicit noneを追加。
  atom/field検査の未使用IEEE importを削除。
- exx_k_orbitsの型・通信importを唯一の使用手続きへ限定。
  FFTW宣言に必要なiso_c_bindingはmodule scopeに保持。
- exx_spatialのint64 MPI reduction importを使用手続きへ限定し、USE_MPI保護を保持。
  communicationには対応するint64 wrapperがないため直接MPIは残す。
- 別リポジトリSALMON-DOCS（監査開始時3553bc6、clean）の
  source/input_keyword_list.rstを更新。9個のnamelist入力、環境変数切替の廃止、
  PBE先行の自動既定値、batch16、LCFO容量制限廃止、PBE0/PBEhの対称性再利用を反映。
  SALMON本体のdocsだけでは満たしていなかった公式マニュアル更新要件を補完した。

## 確認

- 対象変更ファイルの全手続きで明示的implicit none、Fortran行長132以下。
- 入力宣言、namelist/default、MPI配布、小文字化、yn検査、半径の単位換算、
  variables.logが整合。最終版では実行時環境変数の読取りなし。
- 新moduleのprivate/public、共有framework状態、配列寿命とMPI通信順序、
  checkpoint mode番号/保存形式の維持を確認。
- testsuite新規440の番号・構成・CPU/MPI登録と既存422の入力検査を確認。
  ユーザー指定に従いrestart試験をtestsuiteへ追加していない。
- MPI/HSE/LibXC/ScaLAPACK、serial HSE、serial non-HSEを再ビルドし成功。
  422/423/424/440の計12工程が合格。
- 読取りのみの独立レビューで追加の具体的不具合なし。
- 本体・マニュアルのgit diff --checkが合格。

検証ログはwork/coding-rules-final-20261017に保存。
Sphinx9では既存conf.pyの古いintersphinx_mappingがエラーになる。
tracked設定は変更せず、一時設定で当該mappingのみ無効化した上で、
変更前後のHTMLビルドが成功。既存46警告は同数で、新規警告なし。
この設定互換性問題は今回の入力記述変更以前から存在する。

Fujitsu、GPU、EigenExaは今回未検証であり、全構成の実行保証とはしない。
本体とSALMON-DOCSの補正はいずれも未commit・未push。

# Native k交換の自動選択整理（2026-10-19）

通常入力は3つの最適化スイッチを省略する。比較時だけ&functionalに
exx_k_algorithm=auto/reference/pair/orbitsを指定する（既定auto）。
autoは適格なスカラーに既存Hermitian pair、SOCと未対応条件には従来経路を使う。
SOC pairの富岳性能・RSSは未検証なので、自動的には優先しない。
物理条件・精度・SR打切り・履歴近似・診断設定を勝手に変更しない。

旧yn_exx_fft_orbits、yn_exx_k_pair_symmetry、yn_exx_spinor_pair_fftは互換入力として保持。
内部の未指定値を区別して一度だけ解決し、明示nも尊重する。
新しい手動指定と矛盾する旧指定、未知の値、未対応の手動指定は拒否する。
旧設定との優先関係・既定値はstateless exx_k_selectionに一本化した。
requested、実際のselected、選択理由をログへ記録する。

Coding Rulesのdefault/read/broadcast/lowercase/variables.log/値・互換性チェックを実装。
新規Fortran単位にimplicit none、局所use/only、変更行132文字以内。
標準testsuite445を追加し、18通りの解決規則と生成入力の省略チェックを登録。
既存440は手動orbits、441は自動scalar_pairの実選択も確認する。
prepare.pyは冗長な既定値を省略し、14ケースを新規生成してpreflight合格。
preflightは生成入力向けの簡易検証で、汎用namelist解析器ではない。

独立レビューで指摘されたroot判定、手動指定のtheory制約・実選択確認、
preflightの省略処理、比較値のNaN見逃しを修正。最終レビューに未解決の必須指摘なし。

GNU Fortran / MPI / OpenMP / HSE有効のout-of-sourceビルド成功。
最終binary SHA256 f440b276a0b8459e65aa08b20a0f279be9b19f2fdcfd04f690f010485e70356e。
標準CTest25/25、既存回帰17/17、生成入力チェック2/2合格。
新モジュール未実装時の失敗、preflightヘルパー追加前の失敗も確認した。

粗い矩形Si2のSOC/NLCC開発fixtureで、設定を全省略したGS3条件が残差<1e-8で収束。
x/y各20stepのRTは対称性ON/OFFとACE/直接交換の比較に合格。
偏光方向の対称性はOFF（x=nyy、y=yny）。auto-gs-rt.jsonを参照。
手動pair GSも収束、autoとのエネルギー差1.0e-8eV。
同じ保存済みGSからauto/reference/pair/旧pairを両偏光で各20step比較。
auto/referenceは一致、pair/旧pairとの差はx=3.715e-16、y=5.524e-15。
全出力の有限性・16列・ステップ数と実選択を検査。矛盾・未知・未対応条件5例を拒否。
manual-comparison.jsonを参照。WSe2の物理収束や性能の証明ではない。

HSE無効・MPI無効の追加ビルドで、既存noncollinear処理の純対称性ヘルパーが
ビルド対象にないことを検出。NOT USE_HSE時だけ対象へ加え、計算処理を変えず修正した。
この構成のビルド、および標準C2H2 GSとRT200stepが正常終了。
初回の検証2件はpython起動名が環境にないため失敗した。
計算済み出力に対しPython環境を指定して検証段階だけ再実行し2/2合格。
HSE有効構成のbinary SHAはこのビルド条件修正前後で同一。

SALMON-DOCS source/input_keyword_list.rstを更新：a9c06214b807b06108577a4fc2f5d2dc6ac6d312。
変更した節をdocutilsで警告なしに検証。サイト全体の生成は未実施。
Fujitsu再ビルド、GPU/OpenACCの実行、学習経路の新しい手動指定組み合わせ、
WSe2本計算の性能評価は今回実施していない。再起動形式は変更していない。
富岳の既存バイナリ、配置済み専用版、ジョブ入力・起動先・資源は一切変更していない。

主な再現手順：
- cmake --build build-local -j4
- ctest --test-dir build-local -R '44[0-5]_' --output-on-failure
- validate_symmetry.py --salmon build-local/salmon --pseudo validation/Si_PBE19.vps --output NEW --steps 20
- validate_k_algorithms.py --salmon build-local/salmon --reference NEW --output COMPARE_NEW

入力例はdocs/inputs/hse-feature-selection.md。通常はxcと物理条件だけを指定し、
比較時だけexx_k_algorithmを追加する。旧スイッチを新規入力へ追加する必要はない。

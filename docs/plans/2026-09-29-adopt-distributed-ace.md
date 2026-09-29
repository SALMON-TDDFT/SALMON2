# ACE分散経路の採用

ユーザーは実計算比較に基づきACE分散を採用した。

目的：native Γ点ではACEだけを分散し、MLWFは従来経路を使う。

1. exx_nativeのMLWF呼出しと逆変換を計測済みACE-only variantと一致させる。ACEのcomm_matrixは維持する。
2. MLWF分散モジュールと単体試験は比較研究用に残す。新しい入力項目は追加しない。
3. benchmarkの3経路生成を新しい既定に合わせる。既存の計測結果は変更しない。
4. ビルド、HSE06/PBE0/PBEh(40)のsource ACE・空間／軌道並列RTと対選別回帰を検証する。
5. 採用状態を実装ノートに反映する。今回pushは行わない。

## 実施結果

- 通常ソースのexx_nativeが実測済みace variantとバイト単位で一致。
- GNU MPI/ScaLAPACKビルド成功。
- test_source_ace.py: 7試験合格（HSE06/PBE0/PBEh(40)、空間・軌道分割、入力ガード）。
- test_pair_screen.py: 6試験合格（対省略、診断、SCF/RT、逆変換）。
- benchmark生成を再実行し、3経路すべてのexx_nativeが測定時ソースと一致することを確認。
- git diff --check合格。富岳・GPUでの追加検証は未実施。

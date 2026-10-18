# HSE06 NLCC・SOC ローカル検証

基準ソース: `0f1bddf3c25535860471e69412142d62f79d0fa4`。独立ブランチ: `codex/hse-nlcc-soc`。
富岳の既存バイナリ・ジョブ52445014・入力は変更していない。

固定核・周期系・非DCのnative HSE06に、固定コアNLCCと完全スピノル交換を実装。
CPU、均一kメッシュ、k方向のみのMPI分割、整数占有、Taylor4 RTに限定する。
再現入力と実行手順は `samples/hse_soc_nlcc/README.md` を参照。

## 検証結果

GNU Fortran/OpenMPI、Libxc、FFTW、OpenBLASでローカルビルド。
最終バイナリSHA256: `877adaa12e7e5e526b4a88c3a0c27363feb7fc98f13505a416648fdd08a43fbc`。

- 既存Au非ハイブリッドSOC参照試験501/502: GSと1000ステップRTとも合格（MPI 4ランク×OpenMP 2スレッド）。元の入力・参照値・検証スクリプトを使用。ログは `validation/legacy-au-final` に保持。
- 単体・MPI試験8件合格。NLCC密度/勾配微分、真空、スピン極限、非対角Fock、SU(2)、ACE、k分割と不正入力の集団拒否を確認。
- OpenMP 2スレッドの小セルGS/RT一式: PASS。数値は `local-results.json` に保存。
- SOC+NLCC GS残差: 9.63e-8、電子数8。PBE予備SCFからHSEへの切替も収束。
- RTのACE/直接計算の電流差: 3.46e-12 au。無電場エネルギー変動: 2.90e-12 Ha。
- 時間刻み誤差の減少比: 5.46。非共線XC行列の中点平均を含む。

試験コマンド:

```sh
# samples/hse_mlwf_reference 内
RUN_MPI_TESTS=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python3 -m unittest test_native_spin_nlcc test_native_exchange
# ソースルート内
python3 samples/hse_soc_nlcc/validate.py --salmon build-local/salmon --pseudo validation/Si_PBE19.vps --output validation/threaded-final --threads 2
```

生ログと実行入力は開発worktreeの `validation/threaded-final` に保持。

## 適用上の限界

Si 2原子・粗い8³格子・Gamma点の短時間実装検証であり、WSe2誘電関数や物性の収束検証ではない。
スカラーとj依存の擬ポテンシャル演算子の等価性は保証していない。
富岳コンパイラでの新規ビルド、W/Se完全相対論的擬ポテンシャルの選定、
WSe2の格子・真空・k点・時間刻み・伝播時間の収束確認は別途必要。
既存実行ファイルを置換せず、新しいソース/ビルド領域を使うこと。

# DG-Fragment RT-TDDFT クイックスタートガイド

**バージョン:** SALMON 2.2.2  
**実装日:** 2026年2月21日  
**状態:** 実装完了・HPC環境テスト待ち

---

## 概要

DG-Fragment法は、DC-LCFO基底状態計算で得られた局在化フラグメント基底を用いて、実時間時間依存密度汎関数理論（RT-TDDFT）計算を行う手法です。

### 主な特徴

- ✅ **実数基底関数**: DC-LCFOの実数波動関数を使用
- ✅ **適応的基底更新**: 強電場下で自動的に基底を拡張
- ✅ **ハロー交換**: フラグメント境界でのMPI通信（周期境界条件対応）
- ✅ **擬ポテンシャル完全対応**: 局所・非局所擬ポテンシャル両方
- ✅ **高次精度**: 4次精度有限差分ステンシル

---

## 基本的な使い方

### ステップ1: DC-LCFOによる基底状態計算

```fortran
&calculation
  theory = 'dft'
  yn_dc = 'y'
/

&control
  sysname = 'C2H2'
/

&system
  ! (通常のシステムパラメータ)
  iperiodic = 0
  al = 16.0, 16.0, 16.0
  nstate = 6
  nelec = 10
  nelem = 2
  natom = 4
/

&dg_fragment
  num_fragment = 2       ! フラグメント分割数
  nstate_frag = 10      ! フラグメントあたりの状態数
/

&pseudo
  ! (擬ポテンシャルファイル指定)
/

&rgrid
  dl = 0.25, 0.25, 0.25
/

&scf
  ncg = 500
/
```

**実行:**
```bash
mpirun -np 4 salmon < input_gs.inp
```

**出力:** `out_frag_*/frag_XXXXXX/wf_XXXXXX.bin` に基底関数が生成される

---

### ステップ2: DG-Fragmentによる実時間計算

```fortran
&calculation
  theory = 'tddft_pulse'
  yn_dg_fragment_rt = 'y'
  yn_conventional_from_dcdft = 'y'
/

&control
  sysname = 'C2H2'
/

&system
  ! (ステップ1と同じシステムパラメータ)
/

&tgrid
  nt = 10000
  dt = 0.08   ! fs
/

&dg_fragment
  num_fragment = 2      ! ステップ1と一致させる
  nstate_frag = 10      ! ステップ1と一致させる
  
  ! オプション: 適応的基底更新
  yn_adaptive_basis = 'y'
  basis_update_threshold = 0.1   ! eV
  
  ! オプション: 時間積分法
  time_integrator_dg_fragment = 'ssprk3'  ! 'ssprk3', 'rk4', 'aetrs'
  
  ! オプション: 基底更新手法
  yn_dc_cg_basis_update = 'y'  ! 'y': DC-CG法（推奨）, 'n': 対角化
/

&emfield
  ae_shape1 = 'Acos2'
  E_amplitude1 = 1.0e-3  ! a.u.
  I_wcm2_1 = 1.0e12      ! W/cm^2
  tw1 = 6.0              ! fs
  omega1 = 0.05696       ! a.u. (800 nm)
  epdir_re1 = 0.0, 0.0, 1.0
/

&analysis
  ! (通常の解析パラメータ)
/
```

**実行:**
```bash
mpirun -np 4 salmon < input_rt.inp
```

---

## 入力パラメータ詳細

### 必須パラメータ

| パラメータ | 説明 | 例 |
|-----------|------|-----|
| `yn_dg_fragment_rt` | DG-Fragment法の有効化 | `'y'` |
| `yn_conventional_from_dcdft` | DC-LCFO基底の使用 | `'y'` |
| `num_fragment` | フラグメント分割数（GS計算と一致） | `2` |
| `nstate_frag` | フラグメントあたりの状態数（GS計算と一致） | `10` |

### オプションパラメータ

| パラメータ | デフォルト | 説明 |
|-----------|----------|------|
| `yn_adaptive_basis` | `'n'` | 適応的基底更新の有効化 |
| `basis_update_threshold` | `0.1` eV | 基底更新判定のしきい値 |
| `time_integrator_dg_fragment` | `'ssprk3'` | 時間積分法 |
| `yn_dc_cg_basis_update` | `'n'` | DC-CG法による基底更新 |

---

## 出力ファイル

### 標準出力

- `rt_dg_fragment_energy.data`: エネルギーの時間発展
- `rt_dg_fragment_dipole.data`: 双極子モーメントの時間発展
- `rt_dg_fragment_current.data`: 電流密度の時間発展

### ログファイル

- `rt_dg_fragment_basis_update.log`: 基底更新イベント記録
- `salmon.log`: 通常のSALMONログ

---

## 実装済み機能

### ✅ 完了項目

1. **実数基底関数**: DC-LCFOの実数波動関数を直接使用（ゲージ回転不要）
2. **ハロー交換**: フラグメント境界でのMPI通信
   - 最大26近傍（3³-1）
   - 非ブロッキング通信（MPI_ISEND/IRECV）
   - 周期境界条件対応
3. **擬ポテンシャル統合**:
   - 局所擬ポテンシャル（Vpsl）: 全処理で考慮
   - 非局所擬ポテンシャル（V_NL）: プロジェクター法で実装
4. **運動エネルギー行列**: 4次精度有限差分、全グリッド領域対応
5. **適応的基底更新**:
   - ハミルトニアン変化量の監視
   - DC-CG法による基底拡張
   - 波動関数の新基底への射影
6. **入力パラメータ**: `yn_dc_cg_basis_update`で手法選択

### 🔧 マイナータスク（残存）

1. 非直方体フラグメントのグリッドマッピング最適化
2. スピン取り扱いの精緻化
3. HPC環境での実行テスト

---

## トラブルシューティング

### よくあるエラー

**1. フラグメント数が一致しない**
```
Error: num_fragment mismatch between GS and RT
```
**解決策**: 基底状態計算とRT計算で`num_fragment`と`nstate_frag`を一致させる

**2. 基底ファイルが見つからない**
```
Error: Cannot open wf_XXXXXX.bin
```
**解決策**: 
- `out_frag_*/`ディレクトリが存在するか確認
- 基底状態計算が正常終了したか確認
- 作業ディレクトリが正しいか確認

**3. メモリ不足**
```
Error: Cannot allocate phi_frag
```
**解決策**:
- `nstate_frag`を減らす
- MPIランク数を増やす
- より大きなメモリのノードを使用

**4. ハロー交換エラー**
```
Error: phi_frag not allocated with halo regions
```
**解決策**: コード実装エラー。`init_halo_communication`が正しく呼ばれているか確認

---

## 物理的背景

### フラグメント基底表現

波動関数をフラグメント基底で展開:

```
|ψ(t)⟩ = Σ_i c_i(t) |φ_i⟩
```

ここで`|φ_i⟩`はDC-LCFO計算で得られた局在化軌道。

### ハミルトニアン行列

時間依存ハミルトニアンの行列要素:

```
H_ij(t) = ⟨φ_i| T + V_local + V_NL + A(t)·p + A²(t)/2 |φ_j⟩
```

- **T**: 運動エネルギー演算子（-∇²/2）
- **V_local**: 局所ポテンシャル（Hartree + XC + 局所PP）
- **V_NL**: 非局所擬ポテンシャル
- **A(t)**: ベクトルポテンシャル（速度ゲージ）
- **p**: 運動量演算子

### 非局所擬ポテンシャル

Kleinman-Bylander分離可能形式:

```
V_NL = Σ_ilma |proj_ilma⟩ ⟨proj_ilma|

⟨φ_i|V_NL|φ_j⟩ = Σ_ilma ⟨φ_i|proj⟩ · V_coeff · ⟨proj|φ_j⟩
```

数値精度のため、`V_coeff = rinv_uvu`は一度だけ適用。

---

## 性能特性

### 計算コスト

| 操作 | スケーリング | 備考 |
|------|-------------|------|
| 係数の時間発展 | O(N_frag² × N_state²) | 行列積 |
| ハミルトニアン更新 | O(N_frag × N_state² × N_grid) | フラグメント局所 |
| ハロー交換 | O(N_boundary) | MPI通信 |
| 基底更新 | O(SCF反復数) | 必要時のみ |

### 推奨リソース

**小規模系（〜10原子）:**
- CPU: 4-8コア
- メモリ: 〜8 GB/ランク
- 時間: 数時間

**中規模系（〜50原子）:**
- CPU: 16-32コア
- メモリ: 〜16 GB/ランク
- 時間: 数十時間

---

## 実装詳細

### 重要なサブルーチン

| サブルーチン | 行番号 | 機能 |
|------------|--------|------|
| `init_halo_communication` | 388-478 | ハロー通信の初期化 |
| `exchange_phi_frag_halo` | 480-563 | ハロー領域の交換 |
| `calculate_hamiltonian_matrix` | 823-930 | ハミルトニアン行列計算 |
| `add_nonlocal_pp_matrix` | 1028-1133 | 非局所PP行列要素 |
| `apply_kinetic_to_basis` | 939-1026 | 運動エネルギー演算子適用 |
| `check_basis_update_trigger` | 1750-1860 | 基底更新判定 |
| `trigger_basis_update` | 1967-2110 | 基底更新実行 |

### 修正履歴

**2026-02-21**: 初期実装完了
- 実数基底関数（ゲージ回転除去）
- ハロー交換実装（周期境界対応）
- 擬ポテンシャル完全対応
- 運動エネルギー行列（4次精度）
- 適応的基底更新（DC-CG統合）
- 数値精度改善（rinv_uvu×1回）

---

## 次のステップ

1. **HPC環境でのビルド**
   - CMake設定の確認
   - MPI/OpenMPの最適化

2. **テストケース実行**
   - 小規模系での検証
   - 標準RT計算との比較

3. **性能プロファイリング**
   - ボトルネック特定
   - メモリ使用量の最適化

---

## 参考情報

### 関連ファイル

- **メインモジュール**: `src/rt/rt_dg_fragment.f90` (2742行)
- **入力パラメータ**: `src/io/salmon_global.f90`, `src/io/inputoutput.f90`
- **DC-LCFO基底**: `src/gs/dc/`
- **非局所PP**: `src/common/nonlocal_potential.f90`

### 詳細ドキュメント

より詳しい情報は以下を参照:
- `doc/DG_Fragment_Implementation_Notes.md` - 完全な実装ノート（英語）
- SALMONマニュアル - DC-LCFO手法の説明

---

## 連絡先

実装に関する質問:
- モジュール: `src/rt/rt_dg_fragment.f90`
- 実装日: 2026-02-21
- バージョン: SALMON 2.2.2

---

*以上*

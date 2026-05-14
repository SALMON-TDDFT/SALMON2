# DC から DG-RT への使い方

Status: Current

このメモは、DC 計算で作ったフラグメント基底を DG-RT の初期基底として使うための、現行の正規ルートをまとめたものです。対象は `yn_dc_for_dg='y'` による DG seed 直接出力です。

## 1. 基本方針

DC 側では、CG で得たフラグメント状態から DG-RT seed を直接出力します。このとき core 領域の overlap を対角化し、overcomplete になり得る成分を `lambda_cut` で除去します。除去後の同じ線形変換を buffer 領域にも適用し、RT 側はその cutoff 後の本数だけを物理基底として読みます。

正規ファイルは次の組です。

- `data_dcdft/fragments/<ifrag>/basis_functions_buffer.bin`
- `data_dcdft/fragments/<ifrag>/rgrid_index.bin`
- `data_dcdft/fragments/<ifrag>/wavefunctions.bin`
- `data_dcdft/total/dg_occupation.bin`

`basis_functions.bin` ではなく `basis_functions_buffer.bin` を使います。buffer 情報は RT 側の stencil 計算に必要なので、現在の DG-RT seed では必須です。

## 2. DC 入力の要点

DC は `theory='dft'` と `yn_dc='y'` で走らせます。DG-RT 用 seed を出すには `yn_dc_for_dg='y'` を指定します。

```fortran
&calculation
  theory = 'dft'
  yn_dc = 'y'
/

&dc
  num_fragment(1:3) = 2,2,2
  num_rgrid_buffer(1:3) = 4,4,4
  nproc_rgrid_tot(1:3) = 2,2,2
  nstate_frag = 64
  yn_dc_lcfo = 'n'
  yn_dc_lcfo_diag = 'n'
  yn_dc_for_dg = 'y'
/
```

`nstate_frag` は cutoff 前の候補状態数です。RT で使う基底本数ではありません。C64 の buffer=4, 2x2x2 fragment の確認では、DC 側 `nstate_frag=64` から cutoff 後に各 fragment 60 本になりました。`nstate_frag=32` は候補数として不足する場合があります。

`yn_dc_lcfo='y'` と `yn_dc_lcfo_diag='y'` は、通常の LCFO 対角化ルートです。`yn_dc_for_dg='y'` の場合は DG seed 直接出力が優先されるため、現在の DG-RT buffer seed では `yn_dc_lcfo='n'`, `yn_dc_lcfo_diag='n'` を推奨します。

DC の入力チェックでは `base_directory='./'` が前提です。富岳などで quota を避ける場合は、`base_directory` を変えるのではなく、ジョブの実行ディレクトリ自体を容量のある scratch/work 領域に置きます。

## 3. DC 実行後の確認

DC ログで次を確認します。

```sh
grep -E "DG seed basis cutoff|DG seed export|end SALMON|SIG|STOP|Error|FATAL" dc.log
```

正常例:

```text
DG seed basis cutoff: spin=1 min_nbasis=60 max_nbasis=60 n_mat=480
DG seed fragment basis file rows=60 original nstate_frag=64
DG seed export complete (non-LCFO): wavefunctions.bin + basis_functions_buffer.bin + rgrid_index.bin + dg_occupation.bin
end SALMON
```

`min_nbasis` と `max_nbasis` は cutoff 後の物理基底本数です。seed ファイルは cutoff 後の最大本数だけを書き、除去されたゼロ基底列は保存しません。RT 側は入力の `nstate_frag` ではなく、seed ファイル内の本数を優先して使います。

## 4. DG-RT 入力の要点

RT 側では `yn_dg_fragment_from_dcdft='y'` と `yn_dg_fragment_rt='y'` を指定します。DG fragment の分割数と buffer 幅は DC と一致させます。

```fortran
&calculation
  theory = 'tddft_response'
/

&parallel
  nproc_k = 1
  nproc_ob = 4
  nproc_rgrid(1:3) = 2,2,2
  process_allocation = 'orbital_sequential'
/

&dc
  yn_dg_fragment_from_dcdft = 'y'
  num_fragment(1:3) = 2,2,2
  num_rgrid_buffer(1:3) = 4,4,4
  nproc_rgrid_tot(1:3) = 2,2,2
  nstate_frag = 32
/

&propagation
  yn_dg_fragment_rt = 'y'
  time_integrator_dg_fragment = 'rk4'
  yn_fix_func = 'n'
  yn_plane_wave_basis = 'n'
  n_plane_waves_dg = 0
/
```

RT 入力の `nstate_frag` が seed ファイルと違う場合でも、現行コードはファイル側の本数を優先します。ログには次のような情報行が出ます。

```text
[INFO] nstate_frag differs: file=60 runtime=32 (using fragment-state count from file)
```

これは異常ではありません。実際の DG-RT 基底本数は `basis_functions_buffer.bin` と `wavefunctions.bin` の header で決まります。

フラグメント内軌道並列は `&parallel` の `nproc_ob` で指定します。DG-RT では既定で `dg_fragment_parallel_mode='orbital'` です。総 MPI 数は、基本的に

```text
np = num_fragment(1) * num_fragment(2) * num_fragment(3) * nproc_ob
```

に合わせます。

現行の orbital mode では、RT 側の親 real-space 並列は fragment 分割と一致させます。

```text
nproc_rgrid(1:3) = num_fragment(1:3)
process_allocation = 'orbital_sequential'
```

`nproc_ob` がフラグメント内の軌道並列数です。`process_allocation='grid_sequential'` のままだと、連続 rank が orbital rank ではなく real-space rank になり、H 構築で全 fragment のポテンシャルを global reduction する遅い経路に落ちます。現在のコードはこの誤配置を初期化時点でエラーにします。

RT-DG の初期化ログでは、`nproc_ob` を `nproc_rgrid` に折り込んではいけません。例えば `num_fragment=8,8,8`, `nproc_ob=4`, RT 側 `nproc_rgrid(1:3)=8,8,8` の場合、期待される表示は `nproc_rgrid : 8 8 8` と `nproc_ob : 4` です。`nproc_rgrid : 32 8 8` のように出る場合は、フラグメント内軌道並列が空間並列へ誤変換されています。

DC 計算時の `&dc nproc_rgrid_tot(1:3)` と DG-RT 計算時の `&parallel nproc_rgrid(1:3)` は一致させる必要はありません。DC 側の `nproc_rgrid_tot` は fragment SCF をどう空間分割して解くかの指定で、seed ファイルには cutoff 後の fragment box と buffer box が保存されます。DG-RT 側では、その seed を読み直したうえで `nproc_ob` によるフラグメント内軌道並列を使います。

## 5. 15 fs impulse 計算の指定

`dt=0.02` a.u. で 15 fs 程度を見る場合、必要ステップ数は約 31000 step です。

```fortran
&tgrid
  dt = 2.0d-2
  nt = 31006
/
```

誘電関数を見る impulse 計算では `theory='tddft_response'` を使い、`response.data` の `Re(eps_*)`, `Im(eps_*)` を確認します。

## 6. DG-RT 実行後の確認

RT ログで次を確認します。

```sh
grep -E "DG-Fragment RT initialization complete|Total states|end SALMON|SIG|STOP|Error|FATAL|NaN|Inf" rt.log
```

正常例:

```text
DG-Fragment RT initialization complete
  Total states: 480
end SALMON
```

`Total states` は cutoff 後の fragment 基底の総数です。C64 buffer=4, 2x2x2 fragment の確認では `60 * 8 = 480` でした。

密度や current の drift を見る場合は、通常ログだけでなく `*_dg_current_decomp.data` の raw charge や energy/current 列も併せて見ます。`rho%f` は全電子数へ正規化されるため、係数ノルムの破れは raw charge 側に先に出ます。

大規模並列では、DG density の実空間 owner への返送は full communicator の `Alltoallv` ではなく、peer ごとの sparse ring exchange を使います。富岳/Tofu で `tofu-mrq-overflow` が出る場合は、まず `nproc_ob` を下げて再現性を見てください。`nproc_ob` を下げると止まる場合、物理入力ではなく通信キュー圧が主因です。

## 7. 容量見積もり

`basis_functions_buffer.bin` は fragment ごとに buffer box 全体を保存します。概算は次です。

```text
bytes/fragment ~= product(nxyz_domain + 2*num_rgrid_buffer) * n_basis * 8
```

C64, 32^3 grid, 2x2x2 fragment, buffer=4, n_basis=60 では、およそ 6.8 MB/fragment、8 fragment 合計で約 54 MB です。64^3 grid や大きい fragment 数では急に大きくなるため、富岳では quota の小さいホームではなく scratch/work 側で実行してください。

`Disk quota exceeded` が seed 出力中に出た場合、計算の物理設定ではなく出力先容量の問題です。古い `data_dcdft/` や `output.*` を整理するか、容量のある作業ディレクトリで再実行します。

## 8. 旧ルートとの違い

通常の DC-LCFO ルートは `yn_dc_lcfo='y'`, `yn_dc_lcfo_diag='y'` で全系 LCFO 状態を作る経路です。この経路は現在の DG-RT buffer seed に必要な `basis_functions_buffer.bin` を生成しません。

現行の DC -> DG-RT では、次を正規ルートとします。

```fortran
yn_dc_for_dg = 'y'
yn_dc_lcfo = 'n'
yn_dc_lcfo_diag = 'n'
```

これにより、overcomplete 除去後の基底を core と buffer の両方へ同じ変換で出力し、DG-RT 側でそのまま読みます。

# DG Density State-Distribution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `calculate_density_from_fragments` の project ステージを状態分散並列化によって高速化する（n_pw=0 で ~32×、n_pw>0 で ~15× 通信削減）。

**Architecture:** n_pw=0 パスでは密度行列構築（zgemm→dsyrk+AllReduce）を状態分散化してブロックループ前に移動。n_pw>0 パスでは係数の upfront bcast と状態スライスによる状態ループへ切り替え、ブロックごとの comm_summation で rho を集約する。

**Tech Stack:** Fortran 90、BLAS (dsyrk/dgemm)、MPI (comm_summation/comm_bcast)

---

## ファイルマップ

| ファイル | 変更内容 |
|---|---|
| `src/rt/dg/rt_dg_density_reconstruct.f90` | タイミング追加・D 前処理・n_pw=0/n_pw>0 アルゴリズム変更（全 Task） |
| `src/rt/dg/rt_dg_fragment_types.f90` | subgroup_send フィールド削除（Task 5） |
| `src/rt/dg/rt_dg_fragment.f90` | subgroup_send allocate/deallocate 削除（Task 5） |
| `src/rt/dg/rt_dg_fragment_soi.f90` | subgroup_send allocate/deallocate 削除（Task 5） |

---

## 前提知識

コードの座標系：
- `dg_frag%coef(global_state_idx, orbital_idx, ispin)` — complex(8)、frag_root のみ保持
- `dg_frag%density_matrix_frag(nbf_max, nbf_max, nspin, ifrag_count)` — complex(8)、フラグメントごとの密度行列
- `distribute_project = (icomm_frag /= COMM_GROUP_NULL .and. isize_frag > 1)` — fragment 内分散フラグ
- `dg_frag%id_frag` — icomm_frag 内のランク番号（0 origin）
- `dg_frag%isize_frag` — icomm_frag のランク数
- `dg_frag%is_frag_root` — `id_frag == 0`

---

## Task 1: タイミング計測を未計時 3 箇所に追加

**目的:** project 合計と breakdown の差をゼロにし、各ボトルネックの寄与を可視化する。動作変更なし。

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

- [ ] **Step 1: 変数宣言を追加（line 36 付近の timing 変数ブロックに追記）**

```fortran
! 既存の宣言の後に追加
real(8) :: time_project_dmat_build
real(8) :: t_dmat0, t_dmat1
```

- [ ] **Step 2: 初期化を追加（line 68 付近の time_project_* = 0 ブロックに追記）**

```fortran
time_project_dmat_build = 0.0d0
```

- [ ] **Step 3: zgemm + bcast(D) の周囲に cpu_time を追加（line 547 付近）**

現在のコード（line 547-566）を以下に置き換える：

```fortran
              if (.not. dg_frag%density_matrix_frag_valid(ispin, i_local)) then
                call cpu_time(t_dmat0)
                ! nocc >> nbf: build D on root only, bcast D (O(nbf^2)) instead of coef (O(nbf*nocc))
                if (.not. distribute_project .or. dg_frag%is_frag_root) then
                  coef_occ_weighted(1:nbf, 1:nocc_spin) = (0.0d0, 0.0d0)
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
                  do io = 1, nocc_spin
                    do idx_local = 1, valid_basis_count
                      istate_frag = valid_basis_ids(idx_local)
                      coef_occ_weighted(istate_frag, io) = occ_scale * dg_frag%coef(basis_gid(istate_frag), io, ispin)
                    end do
                  end do
!$omp end parallel do
                  call zgemm('N', 'C', nbf, nbf, nocc_spin, (1.0d0, 0.0d0), coef_occ_weighted, nbf_max, &
                             coef_occ_weighted, nbf_max, (0.0d0, 0.0d0), dg_frag%density_matrix_frag(1, 1, ispin, i_local), nbf_max)
                end if
                if (distribute_project) then
                  call comm_bcast(dg_frag%density_matrix_frag(:, :, ispin, i_local), dg_frag%icomm_frag, 0)
                end if
                dg_frag%density_matrix_frag_valid(ispin, i_local) = .true.
                call cpu_time(t_dmat1)
                time_project_dmat_build = time_project_dmat_build + (t_dmat1 - t_dmat0)
              end if
```

- [ ] **Step 4: breakdown 出力行（line 714-716 付近）に dmat 項目を追加**

```fortran
        write(*,'(1x,a,7(a,1pe12.4))') "        density trace: project breakdown", &
          " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
          " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, &
          " over=", time_project_overhead, " dmat=", time_project_dmat_build
```

- [ ] **Step 5: ビルド**

```bash
cd /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build && make -j
```

Expected: コンパイルエラーなし

- [ ] **Step 6: 動作確認**

実行後の density trace 出力で `dmat=` が表示され、
`setup + psi + rho + grid + phi + over + dmat ≈ total project time` となることを確認。

- [ ] **Step 7: コミット**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "perf: add timing to untimed density matrix build regions"
```

---

## Task 2: n_pw=0 — D 構築をブロックループ前に移動（構造的リファクタ）

**目的:** D をフラグメントごとに一括計算し、ブロックループから切り離す。アルゴリズム変更なし。結果は bit-identical であること。

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

- [ ] **Step 1: ローカル変数 D_frag_re を宣言（line 45 付近の real(8), allocatable ブロックに追記）**

```fortran
    real(8), allocatable :: D_frag_re(:,:,:)  ! (nbf_max, nbf_max, nspin) per fragment
```

- [ ] **Step 2: アロケーションを fragment ループ直前に追加（line 340 の `call cpu_time(t_project0)` の直前）**

```fortran
    allocate(D_frag_re(nbf_max, nbf_max, system%nspin))
```

- [ ] **Step 3: fragment ループ内、ブロックループ（`do block_offset`）の直前に D プリパスを挿入**

`block_idx_global = 0` の直後、`do ifrag = ...` ループ本体の最初の `call cpu_time(t_setup0)` より前（line 358 付近）に以下を追加する。

```fortran
        ! --- D pre-pass: compute density matrix for all spins before block loop ---
        D_frag_re(:,:,:) = 0.0d0
        do ispin = 1, system%nspin
          nbf = dg_frag%n_basis(ifrag, ispin)
          nocc_spin = nocc_per_spin
          if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
            nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
          end if
          if (nbf <= 0 .or. nocc_spin <= 0) cycle
          valid_basis_count = 0
          do istate_frag = 1, nbf
            basis_gid(istate_frag) = dg_frag%index_basis(istate_frag, ifrag, ispin)
            if (basis_gid(istate_frag) < 1 .or. basis_gid(istate_frag) > dg_frag%n_mat_max) cycle
            valid_basis_count = valid_basis_count + 1
            valid_basis_ids(valid_basis_count) = istate_frag
          end do
          call cpu_time(t_dmat0)
          if (.not. distribute_project .or. dg_frag%is_frag_root) then
            coef_occ_weighted(1:nbf, 1:nocc_spin) = (0.0d0, 0.0d0)
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
            do io = 1, nocc_spin
              do idx_local = 1, valid_basis_count
                istate_frag = valid_basis_ids(idx_local)
                coef_occ_weighted(istate_frag, io) = occ_scale * dg_frag%coef(basis_gid(istate_frag), io, ispin)
              end do
            end do
!$omp end parallel do
            call zgemm('N', 'C', nbf, nbf, nocc_spin, (1.0d0, 0.0d0), coef_occ_weighted, nbf_max, &
                       coef_occ_weighted, nbf_max, (0.0d0, 0.0d0), dg_frag%density_matrix_frag(1, 1, ispin, i_local), nbf_max)
          end if
          if (distribute_project) then
            call comm_bcast(dg_frag%density_matrix_frag(:, :, ispin, i_local), dg_frag%icomm_frag, 0)
          end if
          dg_frag%density_matrix_frag_valid(ispin, i_local) = .true.
          do io = 1, nbf
            do istate_frag = 1, nbf
              D_frag_re(istate_frag, io, ispin) = real(dg_frag%density_matrix_frag(istate_frag, io, ispin, i_local), kind=8)
            end do
          end do
          call cpu_time(t_dmat1)
          time_project_dmat_build = time_project_dmat_build + (t_dmat1 - t_dmat0)
        end do
        ! --- end D pre-pass ---
```

- [ ] **Step 4: ブロックループ内の D 構築ブロック（line 547-566 付近）をスキップガードに置き換える**

Task 1 で計時を追加した `if (.not. dg_frag%density_matrix_frag_valid...) then ... end if` ブロック（line 547-566）を、以下に置き換える（D は既にプリパスで計算済みなので処理は不要だが、密度行列参照を D_frag_re に切り替える）：

```fortran
              ! D_frag_re(:,:,ispin) already computed in pre-pass above
```

また、line 574 の `density_mat_re` への代入（`density_mat_re = real(density_matrix_frag, ...)`）を D_frag_re から直接取得するよう変更する：

```fortran
!$omp parallel do private(io, istate_frag) schedule(static)
                do io = 1, nbf
!$omp simd
                  do istate_frag = 1, nbf
                    density_mat_re(istate_frag, io) = D_frag_re(istate_frag, io, ispin)
                  end do
                end do
!$omp end parallel do
```

- [ ] **Step 5: deallocate を fragment ループ後（project timing 終了付近）に追加**

```fortran
    if (allocated(D_frag_re)) deallocate(D_frag_re)
```

- [ ] **Step 6: ビルド**

```bash
cd build && make -j
```

Expected: コンパイルエラーなし

- [ ] **Step 7: 数値検証**

実行後に `total_charge` と `elec_num_scaled` を Task 1 のリファレンス値と比較。
`density trace: stage=after-normalize` の行の値が一致すること（差は丸め誤差の範囲 < 1e-10）。

また `dmat=` が依然として 2 秒程度であること（アルゴリズム変更なしのため）。

- [ ] **Step 8: コミット**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "refactor: extract density matrix build to pre-pass before block loop"
```

---

## Task 3: n_pw=0 — zgemm を dsyrk + 状態分散 AllReduce に置き換え

**目的:** D 構築を実数演算 (dsyrk) + icomm_frag 内状態分散に切り替える。isize_frag=4 で ~32× 高速化。

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

- [ ] **Step 1: 新規変数を宣言（line 45 付近に追加）**

```fortran
    real(8), allocatable :: coef_re_frag(:,:)   ! (nbf_max, nocc_spin) real coef for current fragment
    real(8), allocatable :: D_partial_re(:,:)    ! (nbf_max, nbf_max) partial D per rank
    integer :: io_s_frag, io_e_frag, nocc_loc, nocc_per_rank_loc
```

- [ ] **Step 2: アロケーションを D_frag_re と同じ場所に追加**

```fortran
    allocate(coef_re_frag(nbf_max, max(1, nocc_cache)))
    allocate(D_partial_re(nbf_max, nbf_max))
```

- [ ] **Step 3: Task 2 で追加した D プリパスの zgemm+bcast ブロックを dsyrk+AllReduce に置き換える**

プリパス内の以下のコード（Task 2 の Step 3 で追加したもの）：

```fortran
          ! OLD (Task 2 の構造)
          if (.not. distribute_project .or. dg_frag%is_frag_root) then
            coef_occ_weighted(1:nbf, 1:nocc_spin) = (0.0d0, 0.0d0)
            ! ... fill coef_occ_weighted ...
            call zgemm('N', 'C', nbf, nbf, nocc_spin, ...)
          end if
          if (distribute_project) then
            call comm_bcast(dg_frag%density_matrix_frag(:,:,ispin,i_local), dg_frag%icomm_frag, 0)
          end if
          dg_frag%density_matrix_frag_valid(ispin, i_local) = .true.
          D_frag_re(:,:,ispin) = real(dg_frag%density_matrix_frag(:,:,ispin,i_local))
```

を以下に置き換える（ブロック全体を差し替え）：

```fortran
          ! NEW: real dsyrk + state-distributed AllReduce
          ! Step 3a: root fills coef_re_frag, bcasts to all ranks in icomm_frag
          coef_re_frag(1:nbf, 1:nocc_spin) = 0.0d0
          if (.not. distribute_project .or. dg_frag%is_frag_root) then
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
            do io = 1, nocc_spin
              do idx_local = 1, valid_basis_count
                istate_frag = valid_basis_ids(idx_local)
                coef_re_frag(istate_frag, io) = real(dg_frag%coef(basis_gid(istate_frag), io, ispin), kind=8)
              end do
            end do
!$omp end parallel do
          end if
          if (distribute_project) then
            call comm_bcast(coef_re_frag(1:nbf, 1:nocc_spin), dg_frag%icomm_frag, 0)
          end if

          ! Step 3b: each rank computes dsyrk on its state slice
          nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
          io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
          io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
          nocc_loc = max(0, io_e_frag - io_s_frag + 1)

          D_partial_re(1:nbf_max, 1:nbf_max) = 0.0d0
          if (nocc_loc > 0 .and. nbf > 0) then
            ! D_partial = occ_factor * coef_re_frag[:,io_s:io_e]^T * coef_re_frag[:,io_s:io_e]
            ! upper triangle only
            ! 'N': C = alpha*A*A^T, A is (nbf x nocc_loc) stored as (nbf_max x nocc_loc), LDA=nbf_max>=nbf
            call dsyrk('U', 'N', nbf, nocc_loc, occ_factor, &
                       coef_re_frag(1, io_s_frag), nbf_max, &
                       0.0d0, D_partial_re, nbf_max)
          end if

          ! Step 3c: AllReduce partial D across icomm_frag
          if (distribute_project) then
            call comm_summation(D_partial_re(1:nbf_max, 1:nbf_max), &
                                D_frag_re(1:nbf_max, 1:nbf_max, ispin), &
                                nbf_max * nbf_max, dg_frag%icomm_frag)
          else
            D_frag_re(1:nbf_max, 1:nbf_max, ispin) = D_partial_re(1:nbf_max, 1:nbf_max)
          end if

          ! Step 3d: symmetrize (copy upper triangle to lower)
          do io = 1, nbf
            do istate_frag = io + 1, nbf
              D_frag_re(istate_frag, io, ispin) = D_frag_re(io, istate_frag, ispin)
            end do
          end do
          dg_frag%density_matrix_frag_valid(ispin, i_local) = .true.
```

- [ ] **Step 4: deallocate を追加（D_frag_re と同じ場所に）**

```fortran
    if (allocated(coef_re_frag)) deallocate(coef_re_frag)
    if (allocated(D_partial_re)) deallocate(D_partial_re)
```

- [ ] **Step 5: ビルド**

```bash
cd build && make -j
```

Expected: コンパイルエラーなし（dsyrk は BLAS に含まれており別途リンク不要）

- [ ] **Step 6: 数値検証**

```
期待: total_charge が Task 1 比較値と一致（差 < 1e-10）
期待: dmat= が 0.05 秒以下（Task 1 の ~2 秒から大幅改善）
```

isize_frag=1 でも同一結果が出ることを確認：
その場合 `io_s_frag=1, io_e_frag=nocc_spin, nocc_loc=nocc_spin` となり distribute_project=.false. で AllReduce なし → 完全等価。

- [ ] **Step 7: コミット**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "perf: replace zgemm with state-distributed dsyrk for density matrix build"
```

---

## Task 4: n_pw>0 — upfront coef bcast + 状態スライスループ + per-block Reduce

**目的:** n_pw>0 パスでグリッドブロック×バッチごとに発生していた comm_bcast を排除し、状態分散で各ランクが独立に rho を計算する。

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

- [ ] **Step 1: 新規変数宣言（n_pw>0 パス用）**

```fortran
    real(8), allocatable :: coef_re_full(:,:,:)  ! (nbf_max, nocc_spin, nspin) - upfront bcast coef
    real(8), allocatable :: coef_im_full(:,:,:)  ! (nbf_max, nocc_spin, nspin)
    real(8), allocatable :: rho_blk_partial(:)   ! (grid_block_size) - partial rho for state slice
    real(8), allocatable :: rho_blk_full(:)      ! (grid_block_size) - full rho after reduce
```

- [ ] **Step 2: Task 3 の D プリパスと同じループ内に n_pw>0 用 upfront bcast を追加する**

Task 3 の `do ispin = 1, system%nspin` プリパスループを包む `if (n_pw == 0)` / `else` 分岐を作成する。
`n_pw == 0` は Task 3 のコード、`n_pw > 0` は以下：

```fortran
        if (n_pw == 0) then
          ! ... Task 3 の D プリパス (既存) ...

        else
          ! n_pw > 0: upfront bcast of coef_re/im + coef_pw_full_cache
          if (.not. allocated(coef_re_full)) then
            allocate(coef_re_full(nbf_max, nocc_cache, system%nspin))
            allocate(coef_im_full(nbf_max, nocc_cache, system%nspin))
          end if
          coef_re_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          coef_im_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0

          do ispin = 1, system%nspin
            nocc_spin = nocc_per_spin
            if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
              nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
            end if
            nbf = dg_frag%n_basis(ifrag, ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = 0
            do istate_frag = 1, nbf
              basis_gid(istate_frag) = dg_frag%index_basis(istate_frag, ifrag, ispin)
              if (basis_gid(istate_frag) < 1 .or. basis_gid(istate_frag) > dg_frag%n_mat_max) cycle
              valid_basis_count = valid_basis_count + 1
              valid_basis_ids(valid_basis_count) = istate_frag
            end do
            if (.not. distribute_project .or. dg_frag%is_frag_root) then
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
              do io = 1, nocc_spin
                do idx_local = 1, valid_basis_count
                  istate_frag = valid_basis_ids(idx_local)
                  coef_re_full(istate_frag, io, ispin) = real(dg_frag%coef(basis_gid(istate_frag), io, ispin), kind=8)
                  coef_im_full(istate_frag, io, ispin) = aimag(dg_frag%coef(basis_gid(istate_frag), io, ispin))
                end do
              end do
!$omp end parallel do
            end if
            if (distribute_project) then
              call comm_bcast(coef_re_full(1:nbf, 1:nocc_spin, ispin), dg_frag%icomm_frag, 0)
              call comm_bcast(coef_im_full(1:nbf, 1:nocc_spin, ispin), dg_frag%icomm_frag, 0)
            end if
          end do

          ! coef_pw_full_cache は既存（refresh_pw_coef_cache で用意済み）
          ! 全ランクへ bcast（root のみ保持しているため）
          if (distribute_project .and. allocated(dg_frag%coef_pw_full_cache)) then
            call comm_bcast(dg_frag%coef_pw_full_cache(:, 1:nocc_cache, :), dg_frag%icomm_frag, 0)
          end if

          ! 状態範囲を計算（全 ispin 共通）
          nocc_per_rank_loc = (nocc_per_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
          io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
          io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_per_spin)
        end if
```

- [ ] **Step 3: ブロックループ内の n_pw>0 状態ループを状態スライス化する**

既存の n_pw>0 パス（`else` ブランチ、line 596 以降）を以下に変更する。

変更点：
1. `do io0 = 1, nocc_spin, state_block_size` → `do io0 = io_s_frag, io_e_frag, state_block_size`
2. `nbatch = min(state_block_size, nocc_spin - io0 + 1)` → `nbatch = min(state_block_size, io_e_frag - io0 + 1)`
3. `coef_blk_re/im` の fill を upfront の `coef_re_full/im_full` からのコピーに変更
4. `coef_pw_blk` の fill を `coef_pw_full_cache` から直接取得（bcast 削除）
5. `rho_blk_accum` を `rho_blk_partial` に変更し、ループ後に `comm_summation` で集約

```fortran
              ! n_pw > 0 path: state-distributed, no per-batch bcast
              if (.not. allocated(rho_blk_partial)) then
                allocate(rho_blk_partial(grid_block_size))
                allocate(rho_blk_full(grid_block_size))
              end if
              rho_blk_partial(1:npt_blk) = 0.0d0

              do io0 = io_s_frag, io_e_frag, state_block_size  ! ← state slice only
                nbatch = min(state_block_size, io_e_frag - io0 + 1)

                ! copy coef from local upfront buffer (no bcast)
                coef_blk_re(1:nbf, 1:nbatch) = coef_re_full(1:nbf, io0:io0+nbatch-1, ispin)
                coef_blk_im(1:nbf, 1:nbatch) = coef_im_full(1:nbf, io0:io0+nbatch-1, ispin)

                call cpu_time(t_psi0)
                call dgemm('N', 'N', npt_blk, nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                           coef_blk_re, nbf_max, 0.0d0, psi_blk_re, grid_block_size)
                call dgemm('N', 'N', npt_blk, nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                           coef_blk_im, nbf_max, 0.0d0, psi_blk_im, grid_block_size)
!$omp parallel do collapse(2) private(io, igrid) schedule(static)
                do io = 1, nbatch
                  do igrid = 1, npt_blk
                    psi_blk(igrid, io) = cmplx(psi_blk_re(igrid, io), psi_blk_im(igrid, io), kind=8)
                  end do
                end do
!$omp end parallel do
                do ipw0 = 1, n_pw, pw_block_size
                  npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
                  ! direct access from full_cache (no bcast)
                  coef_pw_blk(1:npw_blk, 1:nbatch) = dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)
                  psi_blk(1:npt_blk, 1:nbatch) = psi_blk(1:npt_blk, 1:nbatch) + &
                    matmul(phase_cache(1:npt_blk, ipw0:ipw0+npw_blk-1), coef_pw_blk(1:npw_blk, 1:nbatch))
                end do
                call cpu_time(t_psi1)
                time_project_psi = time_project_psi + (t_psi1 - t_psi0)

                call cpu_time(t_rho0)
!$omp parallel do private(io, igrid, rho_accum) schedule(static)
                do igrid = 1, npt_blk
                  rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                  do io = 1, nbatch
                    rho_accum = rho_accum + occ_factor * real(conjg(psi_blk(igrid, io)) * psi_blk(igrid, io), kind=8)
                  end do
                  rho_blk_partial(igrid) = rho_blk_partial(igrid) + rho_accum
                end do
!$omp end parallel do
                call cpu_time(t_rho1)
                time_project_rho = time_project_rho + (t_rho1 - t_rho0)
              end do  ! io0

              ! Reduce partial rho across icomm_frag → full rho on all ranks
              call cpu_time(t_rho0)
              if (distribute_project) then
                call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_full(1:npt_blk), &
                                    npt_blk, dg_frag%icomm_frag)
              else
                rho_blk_full(1:npt_blk) = rho_blk_partial(1:npt_blk)
              end if
              call cpu_time(t_rho1)
              time_project_rho_reduce = time_project_rho_reduce + (t_rho1 - t_rho0)

              ! write rho_blk_full to rho%f and rho_send (root only)
              if (.not. distribute_project .or. dg_frag%is_frag_root) then
                do idx_local = 1, local_grid_count
                  igrid = local_grid_ids(idx_local)
                  ixg = ixg_buf(igrid)
                  iyg = iyg_buf(igrid)
                  izg = izg_buf(igrid)
                  rho_contrib = rho_blk_full(igrid)
                  rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + rho_contrib
                  rho_s(ispin)%f(ixg, iyg, izg) = rho_s(ispin)%f(ixg, iyg, izg) + rho_contrib
                end do
                do idx_remote = 1, valid_remote_grid_count
                  igrid = valid_remote_grid_ids(idx_remote)
                  owner_rank = owner_buf(igrid)
                  slot = slot_buf(igrid)
                  rho_contrib = rho_blk_full(igrid)
                  rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                  rho_s_send(owner_rank, ispin)%f(slot, 1, 1) = rho_s_send(owner_rank, ispin)%f(slot, 1, 1) + rho_contrib
                end do
              end if
```

注意: 上記のコードは既存の send_pack ブロック（line 666-704 付近）を **丸ごと置き換える**。

- [ ] **Step 4: time_project_rho_reduce の宣言・初期化・出力を追加**

宣言（line 36 付近）：
```fortran
    real(8) :: time_project_rho_reduce
```

初期化（line 68 付近）：
```fortran
    time_project_rho_reduce = 0.0d0
```

breakdown 出力行に追加：
```fortran
        write(*,'(1x,a,8(a,1pe12.4))') "        density trace: project breakdown", &
          " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
          " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, &
          " over=", time_project_overhead, " dmat=", time_project_dmat_build, &
          " rho_red=", time_project_rho_reduce
```

- [ ] **Step 5: deallocate を追加**

```fortran
    if (allocated(coef_re_full)) deallocate(coef_re_full)
    if (allocated(coef_im_full)) deallocate(coef_im_full)
    if (allocated(rho_blk_partial)) deallocate(rho_blk_partial)
    if (allocated(rho_blk_full)) deallocate(rho_blk_full)
```

- [ ] **Step 6: ビルド**

```bash
cd build && make -j
```

- [ ] **Step 7: 数値検証（n_pw>0 ケース）**

n_pw>0 のテストケースで実行し、`total_charge` が旧来値と一致することを確認（差 < 1e-10）。
`setup=` が大幅に減少していること（旧来のバッチ bcast コストの消失）を確認。

- [ ] **Step 8: コミット**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "perf: n_pw>0 state-distributed loop with upfront coef bcast and per-block reduce"
```

---

## Task 5: 使用されなくなった subgroup_send インフラを削除

**目的:** n_pw>0 パスへの切り替えにより不要になった `density_subgroup_*`、`send_pack` 関連フィールドとロジックを除去してコードを整理する。

**注意:** このタスクは Task 4 が完全に動作することを確認した後に実施すること。

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_fragment_soi.f90`

- [ ] **Step 1: 削除対象フィールドを確認**

以下のフィールドが `send_pack` のみで使われていることを grep で確認する：

```bash
grep -rn "density_subgroup_send_count\|density_subgroup_send_slot_map\|density_subgroup_self_ix\|send_pack\|send_sum\|subgroup_send_offset\|density_block_first_offset\|density_block_step\|density_block_nblocks" \
  src/rt/dg/ --include="*.f90"
```

- [ ] **Step 2: rt_dg_density_reconstruct.f90 から send_pack 関連ロジックを削除**

削除対象（n_pw>0 パスへの切り替え後に使われなくなったもの）：
- `send_pack`, `send_sum` の allocate/deallocate とゼロ初期化
- `subgroup_send_offset` の allocate/deallocate
- `density_subgroup_send_count/slot_map` の初期化ブロック（line 176-241 付近）
- `density_block_nblocks/first_offset/step` の初期化ブロック（line 221-241 付近）
- line 721-761 の `comm_summation(send_pack, send_sum, ...)` ブロック全体

削除後、`if (distribute_project .and. allocated(send_sum)) then ... end if` ブロック全体を削除する。

- [ ] **Step 3: rt_dg_fragment_types.f90 からフィールドを削除**

```bash
grep -n "density_subgroup\|density_block_first\|density_block_step\|density_block_nblocks" \
  src/rt/dg/rt_dg_fragment_types.f90
```

確認後、該当する `logical/integer, allocatable` フィールド行を削除。

- [ ] **Step 4: rt_dg_fragment.f90 と rt_dg_fragment_soi.f90 から allocate/deallocate を削除**

```bash
grep -n "density_subgroup\|density_block_first\|density_block_step\|density_block_nblocks" \
  src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_fragment_soi.f90
```

確認後、該当行を削除。

- [ ] **Step 5: ビルド**

```bash
cd build && make -j
```

Expected: コンパイルエラーなし（削除したフィールドへの参照が残っていればエラーが出るので修正する）

- [ ] **Step 6: 数値・性能の最終確認**

```
期待: total_charge が Task 1 比較値と一致
期待: project breakdown の全項目の合計 ≈ project 合計（未計時ゼロ）
期待: n_pw=0 時の project 時間 < 0.2 s（旧 2.1 s）
```

- [ ] **Step 7: コミット**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90 \
        src/rt/dg/rt_dg_fragment_types.f90 \
        src/rt/dg/rt_dg_fragment.f90 \
        src/rt/dg/rt_dg_fragment_soi.f90
git commit -m "cleanup: remove unused subgroup_send infrastructure after state-distribution refactor"
```

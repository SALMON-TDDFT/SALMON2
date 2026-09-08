# 新セッション指示書: WPW S直交PW補空間

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## 新セッション開始文

```text
/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement/docs/plans/2026-07-23-wpw-s-orthogonal-pw-complement-session-handoff.md
を読み、executing-plansを使ってWPW S直交PW補空間計画のTask 1から開始してください。
専用worktreeとbranchを維持し、REDを確認してから実装してください。
```

## 作業場所

- worktree:
  `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`
- branch: `codex/wpw-s-orthogonal-complement`
- HEAD: `38f1c06` (`plan(dg-wpw): implement S-orthogonal PW complement`)
- 状態: clean。Task 1のコード編集はまだない。

root worktree
`/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG`にはユーザーの多数の未コミット差分と
stage2d runがある。root worktreeをcheckout、reset、clean、stash、削除してはならない。
実装・テスト・commitは必ず上記専用worktreeで行う。

## 必読資料

1. 設計書:
   `docs/plans/2026-07-23-wpw-s-orthogonal-pw-complement-design.md`
2. 実装計画:
   `docs/plans/2026-07-23-wpw-s-orthogonal-pw-complement.md`
3. これまでのoccupied-W診断:
   `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`
4. 直前のmetric-block比較:
   `docs/plans/2026-07-23-wpw-metric-correction.md`

## ユーザーが承認した物理方針

occupied-Wannierを主空間とし、PWをそのS-metric上の完全な補空間とする方式1を採用する。

\[
A=S_{WW}^{-1}S_{WP},\qquad P_\perp=P-WA,
\]

\[
T=\begin{pmatrix}I&-A\\0&I\end{pmatrix},\qquad
S'=T^\dagger ST,\qquad H'=T^\dagger HT.
\]

- `S'_{WP}=0`を作る。
- `H'_{WP}`は消さず、Hamiltonianによる物理的W--PW混成を残す。
- 初回比較ではPW列を削除せず、元のWPW spanを完全に保持する。
- 数値PW-complement rankとcutoff-low weightは診断だけに使う。
- 実空間成分別に変換せず、体積・非局所・DG界面を含む完成済みbounded H/S actionへ
  `T`と`T^dagger`を適用する。
- density/projector等の物理境界では
  `c_W=c'_W-Ac'_P`, `c_P=c'_P`へ戻す。
- density mixingは既存`dg_wpw_scf_mix`だけを使用し、Hamiltonianや変換をmixしない。

## 直前の証拠

Task 19のfragment-local `S_block^{-1}` correctionは、B=6でoccupied residualを
inner 32/96/160で約25.6/149.6/129.2倍改善したが、extra residualを
約1.56/1.14/4.73倍悪化させた。block dimension 138、condition `4.5069E+05`、
inner 160のextra correction ratioは`8.4597E+02`だった。

従って今回否定されたのは局所block-Jacobi近似だけであり、metric-aware correction一般や
occupied-W基底空間ではない。次はW--PW metric重複を同一spanの座標変換で除く。

## baseline verification

専用worktree作成直後に以下を実行し、全てPASS済み。

```bash
python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py
python3 tests/dg/test_wpw_generalized_algebra.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
```

結果:

- `PASS fixed-H WPW relaxation source contract`
- generalized algebra: `Ran 1 test ... OK`
- `PASS bounded matrix-free DG-DC production route contract`

新セッションでは最初に次を再確認する。

```bash
git status --short --branch
git log -3 --oneline
git diff --check
```

## Task 1開始点

実装計画Task 1をTDDで実施する。対象:

- Create: `src/common/dg_wpw_s_orthogonal_complement.f90`
- Modify: `src/common/CMakeLists.txt`
- Create: `tests/dg/test_dg_wpw_s_orthogonal_complement_mpi.f90`
- Create: `tests/dg/run_dg_wpw_s_orthogonal_complement_mpi.py`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`

最初に2-rank dense-oracle fixtureとsource contractを追加し、module/API不在によるREDを
実際に確認する。その後、分散`S_WW A=S_WP` solve、metric-only fingerprint、係数写像、
wrapped H/S action、collective stale rejectionを実装する。

Task 1 production実装ではglobal dense `S_WW`や`A`を単一rankへ集約しない。既存の
owner schedule、fetch/reduce、global Gramを利用する。小規模MPI fixtureだけはdense oracleを
使用してよい。Task 19のfragment-local inverseをprojection solveへ流用してはならない。

## lifecycle

変換の有効性はHamiltonian epochではなくmetric-only fingerprintで判断する。

- fingerprint対象: W/P IDs、ownership、metric convention、WW/WP/PP S entries。
- 対象外: H entries、interface lambda、Hamiltonian operator epoch。
- `set_dg_wpw_interface_lambda`およびH-only potential rebuildでは変換を再利用する。
- Sまたはbasis/layout IDが変われば全rankでstale拒否する。
- 既存bounded H/S action自身のepoch/layout検査は残す。

## 実行・レビュー手順

`executing-plans`とTDDに従う。

1. Task 1 RED。
2. Task 1実装とfocused verification。
3. 仕様適合レビュー。
4. コード品質レビュー。
5. Critical/Importantを全修正して再検証、commit。
6. Task 2のdefault-off fixed-H integrationへ進む。
7. Task 3でTask-16相当B=6、両preconditioner=`n`、history=`y`を8 rank実行する。
8. Task 4で結果をユーザーへ提示し、normal SCF/checkpoint/RTへの昇格は別途承認を得る。

既存dirty差分をcommitへ混ぜない。テストやbuildがsandboxのMPI socket制限で失敗した場合は、
同じ限定コマンドを承認付きで再実行する。

## commit履歴

- `9409c8b` design S-orthogonal PW complement
- `38f1c06` implementation plan and span-preserving refinement

root側で以前作成したローカル`main` merge worktree
`/private/tmp/SALMON2_RTDG-main-merge`はこの新実装より前の状態である。今回の実装先ではない。

# WPW S-Orthogonal PW Complement Design

## 目的

occupied-Wannier空間をWPW法の主空間とし、PWをそのS-metric上の厳密な補空間として扱う。
W--PW metric重複を除去しつつ、Hamiltonianによる物理的なW--PW混成、元のWPW span、
固有値問題、密度表現を保持する。

## 現状と仮説

現在の基底はoccupied-W行とselected-PW行をそのまま結合しているため、metricは

\[
S=\begin{pmatrix}S_{WW}&S_{WP}\\S_{PW}&S_{PP}\end{pmatrix}
\]

である。Task 19のfragment-local `S_block^{-1}`補正はoccupied residualを最大約150倍
改善した一方、extra residualを最終約4.73倍悪化させた。局所block conditionは
`4.51E+05`で、extra correctionは最大約846倍に増幅された。

次の仮説は、PWがoccupied-Wとmetric上で重複したまま追加されているため、補空間方向と
Wannier方向が補正・Ritz空間内で分離されず、小さいmetric modeとextra-state過増幅を
生んでいる、である。

## 採用方式

PW列をoccupied-W spanへS直交射影し、

\[
A=S_{WW}^{-1}S_{WP},\qquad P_\perp=P-WA
\]

と定義する。基底変換を

\[
B=[W\ P],\qquad
T=\begin{pmatrix}I&-A\\0&I\end{pmatrix},\qquad B'=BT=[W\ P_\perp]
\]

と書く。このとき

\[
S'=T^\dagger ST,\qquad H'=T^\dagger HT
\]

であり、有限精度内で`S'_{WP}=0`となる。`H'_{WP}`は消さない。従って固有状態は
HamiltonianによりWと補空間を混成でき、metric重複だけが除去される。

この変換は非特異な座標変換であり、元のWPW spanと一般化固有値を変えない。初回比較では
補空間PW列を削除しない。`S_perp`の数値rankとcutoff以下のweightは診断するが、実際の
方向選別は従来どおりreduced solveのmetric cutoffだけに任せる。これにより座標変換の効果と
基底縮小の効果を混同しない。

## 演算子構成

実空間PW値を個別に変更せず、体積、非局所、DG界面を含む完成済みbounded operatorへ
合同変換を適用する。これにより全Hamiltonian成分へ同じ変換が一貫して掛かる。

変換後係数`c'=(c_W',c_P')`から元基底係数への写像は

\[
c=Tc',\qquad c_W=c_W'-Ac_P',\qquad c_P=c_P'
\]

である。H/S actionは入力を`T`で元座標へ写し、既存bounded actionを実行し、出力を
`T^dagger`で補空間座標へ戻す。密度、projector、checkpointへ係数を渡すときも同じ
`c=Tc'`を使用する。

`A`は完成したreal-space metricからcollectiveに構築する。productionではglobal dense
`S_WW`を一rankへ集約せず、既存owner scheduleとglobal Gramを使う分散Hermitian solveで
`S_WW A=S_WP`を解く。小規模fixtureだけはdense oracleと比較する。solveではHermitian
defect、最小Rayleigh値、condition estimate、finite性、最終residualを検査する。単なる
fragment-local inverseではなく、occupied-W全空間に対する射影とする。

## 生命周期

補空間変換はmetricと基底layoutにのみ依存する。

- fixed-H反復中は固定する。
- `set_dg_wpw_interface_lambda`はHだけを変更するため、変換を再利用する。
- 外側SCFの密度・局所ポテンシャル更新でもSと基底layoutが不変なら再利用する。
- Sブロック、occupied-W/PW ID、ownershipから補空間専用`metric_fingerprint`を作る。
  これが変わればstaleとしてcollectiveに拒否し、再構築する。
- Hamiltonianだけのoperator epoch更新は変換を無効化しない。ただし通常のH/S action自身は
  従来どおりoperator epoch/layout検証を行う。
- density mixingは既存`dg_wpw_scf_mix`だけを使用し、変換やHamiltonianをmixしない。

## 初回比較条件

補空間化そのものの効果を分離するため、初回B=6 runでは新しいpreconditionerを併用しない。

- `yn_dg_wpw_preconditioner='n'`
- `yn_dg_wpw_metric_preconditioner='n'`
- search historyあり
- Task 16と同じnormalized occupied-W seed、basis、cutoff、publication gate
- fixed-Hの既存`interface_lambda` continuationを維持

## 診断と受入条件

実装fixtureで以下を検証する。

1. `||W^dagger S P_perp||`がscale相対`1E-11`以下。
2. 変換前後のprojector/span defectが`1E-11`以下。
3. numerical PW-complement rankとcutoff以下のmetric weightを記録するが、基底列は削除しない。
4. 変換前後の一般化固有値がrank保持ケースで`1E-10`相対以内。
5. 係数往復と密度再構築が`1E-11`以内。
6. interface lambda変更後も変換が有効。
7. epoch/fingerprintまたはmetric変更後はstale変換をcollectiveに拒否。
8. normal production route、fixed-H、continuationのcollective順序を維持。

物理比較ではinner 32/96/160のoccupied/extra residual、effective rank、metric spectrum、
Ritz boundary defect、最終`info`、publication状態をTask 16およびTask 19と比較する。

## 解釈

補空間化で改善すれば、occupied-W空間自体よりW--PW metric重複と座標条件が主要因だったと
判断する。改善しなければ、同一span内の単なる座標変換ではplateauを解消できず、
`H-epsilon S` correctionやstate別補正へ進む。PW列を実際に削る実験は基底空間変更を
伴うため、必要になった場合は別設計とする。

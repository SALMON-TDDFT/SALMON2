# Si 4×4×4 / 512原子: HSE06 DC-SCF → LCFO GS入力

`inputfile` と同梱の `Si_rps.dat` を同じ実行ディレクトリに置き、計算ノードで64 MPIランクで実行する。
OpenMPスレッド数はジョブの割当てに合わせて設定する。入力のMPI配置を変えずに総MPI数だけ変えない。
この一式はGS専用。RTの波動関数表現は実空間メッシュであり、LCFOは初期状態の再構成に用いる。

## 計算条件

- Siダイヤモンド構造の通常8原子セルを4×4×4に複製: 512原子、価電子2048個。
- 格子定数10.26 bohr、全体セル41.04 bohrの立方体、Γ点。
- 全体格子64³、格子間隔0.64125 bohr。原子座標はinputfile内に全て記載。
- DCコア4×4×4、各コア8原子・16³格子。
- 各軸の両側に8格子点（5.13 bohr）のバッファ。各周期フラグメントは64原子・32³格子。
- 全体Hartreeの空間分割4×4×4、各フラグメント内は1 MPI。合計64 MPI。
- 電子温度300 K、DCの共通化学ポテンシャル合わせ。
- フラグメント256状態、LCFOの全体系出力2048状態。後続RTの占有状態数は1024。
- フラグメントSCFでは局所化なし: yn_exx_dc_mlwf='n'、半径・ノルム切詰めなし、対スクリーニングoff。
- PBE予備SCFの閾値1e-4を3回連続で満たした後にHSE06へ移行。
- HSE06の遮蔽係数0.11 bohr^-1、交換率は実装の既定値25%。rVV10なし。
- SCFの最終閾値1e-8、最大1800反復、CG4回。Broyden係数0.3を初期設定とした。
- LCFO: CheFSI、フィルタ次数60、最大200サイクル、残差閾値1e-7。
  energy_cut=100 hartree、lambda_cut=1e-7は既存Si準備例を踏襲した初期設定。

## ScaLAPACKについて

ビルドにはScaLAPACKを有効にする。HSE06はΓ点でも複素軌道を使い、
`lcfo_eigensolver='scalapack'` と `'chefsi'` の両方に対応する。
今回は大きいLCFO行列から必要な2048状態を求めるため、ScaLAPACKを内部で使う
`lcfo_eigensolver='chefsi'` を指定している。
`&parallel` の `yn_scalapack='n'` は各フラグメント内の通常対角化の指定であり、
LCFO eigensolverの選択やScaLAPACKのビルド設定を無効化するものではない。

## 後続RTへ残すデータ

`yn_dc_lcfo='y'` と `yn_dc_lcfo_diag='y'` により、DC-SCF後にLCFO固有状態を出力する。
`write_gs_restart_data='no'` でもLCFO出力は行われる。通常GSのwfn再開ファイルを使う経路とは異なる。

GS終了後は **data_dcdftディレクトリ全体** を保持し、別の汎関数や計算結果と混ぜない。
各フラグメントのbasis_functions.bin、wavefunctions.bin、rgrid_index.bin等を含む。
後続RTの実行ディレクトリから `data_dcdft` として参照できるようにする。
後続RTは `yn_dc='n'`、`yn_conventional_from_dcdft='y'`、HSE06、同じ原子・セル・メッシュ・
フラグメント配置を使用する。GSの温度指定はRT入力へそのままコピーしない。
本一式にはRT入力・外場設定・ジョブ投入スクリプトは含めていない。

## 確認した範囲

原子数、座標の重複・セル内範囲、全64フラグメントのコア原子数8とバッファ付き原子数64、
入力変数のnamelist登録、擬ポテンシャルの価電子数4を静的に確認した。
Si_rps.datはtestsuites/pseudo/Si_rps.datと同一で、直近の通常Si Nk=4³比較と揃えている。
512原子のGS自体は未実行。SCF収束、格子間隔・バッファ幅・保持状態数の精度収束は実計算で確認する。

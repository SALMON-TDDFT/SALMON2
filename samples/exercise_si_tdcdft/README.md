# Si: experimental macroscopic TDCDFT

この例は、PZ-ALDA に巨視的な交換相関ベクトルポテンシャルを加える実装の出発点です。
**12³実空間格子・2³ k点は動作確認用です。励起子ピークの収束条件ではありません。**

## モデルと入力

`&functional` に以下を指定します。既定値は `tdcdft='none'`、各係数はゼロです。

```fortran
xc = 'PZ'
tdcdft = 'lrc'
tdcdft_alpha = 0.2
```

SALMON内部の `a=Axc/c` と電子数電流密度 `j`（電荷電流は `-j`）に対し、原子単位で

```
a'' + beta*a' + gamma*a = alpha*j
```

を伝播します。`lrc` は beta=gamma=0。`proca` は正の `tdcdft_restoring=gamma` を要求し、
`tdcdft_damping=beta` は任意の非負値です。alphaは原子単位の無次元係数で、betaとgammaの入力単位は
選択した時間単位の逆数、逆数の二乗です。`A_eV_fs` なら fs^-1、fs^-2。
`variables.log` には変換後の原子単位を記録します。

alpha=0.2 は Sun et al. のSiに使われたLRC値を開始点にしたものです。
例の **gamma=0.01 a.u. は安定化の動作確認用で、Siに合わせた物理パラメータではありません**。
LRCの長時間不安定性を、減衰を大きくするだけで解決したと解釈しないでください。
固定alphaによる計算は励起キャリア依存の遮蔽を自己無撞着に記述しません。

CPU、周期系、横応答、固定イオン、非スピン偏極、PZ、middlepoint、対称性縮約なしに限定しています。
理論の指定は `tddft_response` または `tddft_pulse`。GPU、Maxwell結合、MD、DC、jellium、
spin-orbit、aetrs、分散self-checkpoint、再開時の時間原点リセットは拒否します。

## 基底状態とimpulse

各計算を別ディレクトリで実行してください。以下はこのサンプルディレクトリからの例です。
実行ファイルのパスは環境に合わせて変更します。

```sh
export SALMON_EXE=/absolute/path/to/salmon
mkdir gs alda zero lrc proca
cp ../exercise_04_bulkSi_gs/Si_rps.dat gs/
(cd gs && "$SALMON_EXE" < ../Si_gs.inp > outputfile)
for run in alda zero lrc proca; do
  cp ../exercise_04_bulkSi_gs/Si_rps.dat "$run/"
  ln -s ../gs/data_for_restart "$run/restart"
  (cd "$run" && "$SALMON_EXE" < "../Si_impulse_${run}.inp" > outputfile)
done
python3 analyze.py alda/Si_rt.data --impulse 0.001 --output alda.csv
python3 analyze.py lrc/Si_rt.data --impulse 0.001 --output lrc.csv
```

`analyze.py` はPythonとNumPyが必要です。**入力のrt.dataは原子単位に限定**します。
z電流を既定とし、正のA/cステップで規格化します。
SALMONと同じ三次窓 `1-3(t/T)^2+2(t/T)^3`、追加のブロードニングなし、エネルギーシフトなしで
Re/Im epsilonを出力し、条件をJSONに保存します。`--peak-range 2.5 4.5` のように範囲を
指定すると最大値と積分強度を記録できますが、ピークの同定は収束後に行ってください。
本例の全観測時間は約23.2 fsです。周波数刻みを細かくしても有限時間の分解能は改善しません。

比較の順序：alpha=0、0.1、0.2、0.3、キック振幅半減、dt半減、k点増加、実空間格子細分化、
観測時間延長。k点や格子を変更したときは対応するGSも再計算します。
同じ観測時間・窓で比較し、SiのE1/E2の位置・高さ・積分強度を確認します。
E1の移動をそのまま束縛励起子の結合エネルギーと呼ぶことはできません。
PZのバンドギャップ誤差も別に評価する必要があります。

## レーザーと遅延プローブ

`Si_pump.inp` と `Si_pump_probe.inp` は同一GSから別々に実行するテンプレートです。
例は1.6 eV、全包絡長10 fsのAcos2、強度10^10 W/cm²。`tw1` は強度FWHMではありません。
`tdcdft='proca'` の係数は上記の未較正の例です。

Acos2の場合、プローブの閾値は `t1_start + tw1/2 + T1_T2`。本例は約20 fs
（826.827466原子時間）です。SALMONはこの閾値より**後**の格子点でAのステップを適用します。
時間原点の離散化誤差はdt収束で確認してください。

```sh
python3 analyze.py pump_probe/Si_rt.data --subtract pump/Si_rt.data \
  --impulse 0.001 --probe-time 826.827466 --output pumped.csv
```

ポンプ電流を差し引き、プローブ後の時刻を原点にした有限時間応答を求めます。
ポンプのみとプローブ付きの時間格子が異なる場合は解析を拒否します。
両計算のGS、係数、ポンプ、出力間隔は一致させてください。
負のプローブも実行する場合、差 `J(+probe)-J(-probe)` を取り、`--impulse` を振幅の二倍にします。
プローブ振幅半減と、ポンプなしの遅延プローブが通常の線形応答と一致することを確認してから、
強度・遅延依存性を解釈します。これらの長時間・非平衡スペクトルの物理検証は未完了です。

## 出力と再開

- `Si_rt.data` のA/Eは従来どおり古典場。電流はXC場を含む有効場で評価します。
- `Si_rt_xc.data` に時刻、Axc/cの3成分、Excの3成分を出力します。伝播に使うAは古典場とXC場の和。
- `checkpoint_interval` による再開には追加の `tdcdft.bin` が必要です。係数、dt、モードの変更を拒否し、
  XC場の2時刻と応答解析用の電流履歴を保持します。GS restartにはこのファイルは不要です。
- エネルギー出力にモデルXC場の保存エネルギーは追加していません。減衰を実験的な寿命と同一視しません。
- XC有効時は電流用の非局所擬ポテンシャル位相を端点で更新します。alpha=0のimpulseは従来結果と一致します。
  時間依存レーザーでは従来経路の中点位相との差があるため、比較基準には **XC有効・alpha=0** を使用します。

## 参考文献

- Sun et al., PRL 127, 077401 (2021), https://doi.org/10.1103/PhysRevLett.127.077401
- Williams and Ullrich, JCTC 21, 4753 (2025), https://arxiv.org/abs/2501.13290
- Dewhurst et al., PRB 111, L060302 (2025), https://arxiv.org/abs/2401.16140

本実装の係数は正規化した運動方程式の係数です。Kohn–Sham–Proca論文のa0/a2を直接入力するものではありません。

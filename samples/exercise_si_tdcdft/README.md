# Si: experimental macroscopic TDCDFT

この例は、PZ-ALDA に巨視的な交換相関ベクトルポテンシャルを加える実装の出発点です。
**初期解析は12³実空間格子・4³ k点で実施します。ユーザー指定によりk点収束走査は行いません。**

## モデルと入力

基準は [Williams–Ullrich, arXiv:2501.13290v2](https://arxiv.org/html/2501.13290v2)
の式(26)です。SALMON内部の `a=Axc/c` とセル平均電子数電流密度 `j` に対し

```
a'' + beta*a' + gamma*a = alpha*j
```

を伝播します。`tdcdft_alpha` がalpha、`tdcdft_damping` がbeta、
`tdcdft_restoring` がgammaです。alphaは無次元、betaは逆時間、gammaは逆時間の二乗。
`A_eV_fs` の場合はそれぞれfs^-1、fs^-2で入力します。
論文のAはHamiltonian内で直接運動量に加わり、SALMON内部のaに対応します。
ただし論文の電場記号は +dA/dt、SALMONのXC電場出力は `Exc=-da/dt` です。

Siの既存結果と結合強度を揃え、alpha=0.2を維持します。無質量の参照は

```fortran
xc = 'PZ'
tdcdft = 'lrc'
tdcdft_alpha = 0.2
```

です。安定化項の影響を調べるProca入力例（原子単位）は

```fortran
xc = 'PZ'
tdcdft = 'proca'
tdcdft_alpha = 0.2
tdcdft_damping = 0.0
tdcdft_restoring = 0.0001
```

です。**gamma=1e-4はSi用に検証済みの値ではなく、感度比較の開始値です。**
論文の2次元モデルで使われたgamma=0.01等をSiへ直接移植しません。
本例のbeta=0では自由振動の尺度sqrt(gamma)=0.01 a.u.、約0.272 eVです。
以前の例のgamma=0.01は約2.72 eVとなり、今回調べたい光学構造と近すぎるため変更しました。
これは結合した電子系の実際の余分な共鳴位置を保証する値ではありません。

同論文の式(27)では有効カーネルが周波数依存となるため、gammaは単なる数値的な
発散止めではありません。まずbeta=0、alpha固定でgamma=0（LRC）、1e-4、4e-4を比較し、
電流・XC場の長時間挙動と2–4 eVのスペクトル感度を確認します。最後の値の自由振動尺度は
約0.544 eVです。この比較は今後の検証項目で、安定性や最適値を確認した結果ではありません。

式(62)のgamma_thrは平均的なcounter forceから説明されますが、論文自体も
普遍的に予測できる閾値とはしていません。時刻ごとにXC力を打ち消す操作や、
密度からgammaを自動決定する処理は実装していません。既存出力から
`[a(t2)-a(t1)]/(t2-t1)` を複数の時間窓で調べれば平均のda/dtを診断できますが、
有限時間で小さいことだけではゼロ力定理や長時間安定性の証明にはなりません。
また、この論文は全電子モデルであり、SALMONの非局所擬ポテンシャルを含めた
厳密な力の検証とは区別します。

既定は `tdcdft='none'`、すべての係数ゼロです。正規化したProca入力ではgamma>0、
beta>=0を要求し、gamma=0の比較には `tdcdft='lrc'` を使います。
固定係数による計算は励起キャリア依存の遮蔽を自己無撞着に記述しません。

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

比較の順序：alpha=0、0.1、0.2、0.3、キック振幅半減、dt半減、実空間格子細分化、
観測時間延長。k点は4³に固定します。実空間格子を変更したときは対応するGSも再計算します。
同じ観測時間・窓で比較し、SiのE1/E2の位置・高さ・積分強度を確認します。
E1の移動をそのまま束縛励起子の結合エネルギーと呼ぶことはできません。
PZのバンドギャップ誤差も別に評価する必要があります。

## レーザーと遅延プローブ

`Si_pump.inp` と `Si_pump_probe.inp` は同一GSから別々に実行するテンプレートです。
例は1.6 eV、全包絡長10 fsのAcos2、強度10^10 W/cm²。`tw1` は強度FWHMではありません。
`tdcdft='proca'` の係数は上記と同じalpha=0.2、beta=0、gamma=1e-4です。
安定性とgamma感度を確認してからレーザー励起による変化を解釈してください。

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

## 補助的なDewhurst係数による入力

主入力は上記のalpha/beta/gammaです。前段で追加した `tdcdft_a2` と `tdcdft_a0` も
比較用の別表現として利用できます。alpha/restoringとの混在は拒否します。

a2=-85はDewhurst et al., PRB 111 L060302 (2025) の
[補足資料Table 1](https://journals.aps.org/prb/supplemental/10.1103/PhysRevB.111.L060302/SI.pdf)、
a0=-0.2は本文の弱束縛励起子の指定です。図1のa0=+0.25はLiFの安定化例で、Si共通値ではありません。
元のSiスペクトルには0.75 eVのscissor補正と0.22 eVのsmearingが使われていますが、
このサンプルには自動で適用していません。

論文のHamiltonianは(p-A/c)^2/2、SALMONは(p+A/c)^2/2です。
SALMON内部の `a=Axc/c` とセル平均電子数電流密度 `j` では

```
a2*a'' + a0*a = -4*pi*j
alpha = -4*pi/a2
restoring = a0/a2
```

に変換します。Siの値はalpha=0.1478396543、restoring=0.002352941176 a.u.に対応します。
`variables.log` にa2/a0と有効係数を記録します。
a2は無次元、a0は選択した時間単位の逆数の二乗です。`A_eV_fs` ならfs^-2なので、
論文の原子単位のa0=-0.2を入力する際は `-0.2/(0.02418884326505**2)` に換算します。

a2は有限の非ゼロ値で、正負とも指定可能。a0/a2>=0を要求します。
a0=0は質量項なしの比較用に許可します。正の比で得られる自由振動は、結合した電子系の
長時間安定性を保証するものではありません。離散時間刻みには `(a0/a2)*dt**2<4` を要求します。
`tdcdft_damping` は任意の非負の**正規化済み**減衰係数（逆時間）として併用できます。
論文のa1を直接入力するものではありません。

## 強レーザー用の瞬時遮蔽（試験モデル）

`Si_pump_instant.inp` は時間平均なしの推定を試す入力です。追加パラメータは原子単位限定。
`tdcdft_screening='none'` が既定で、従来のimpulse入力と固定alpha伝播は変わりません。
新モデルは `tddft_pulse`、非impulseの第一パルス、正規化alpha入力のみ対応します。

```
tdcdft_screening='instant'
tdcdft_screen_omega=0.05879892  ! representative angular frequency in atomic units
tdcdft_screen_reference=0.0   ! K0: placeholder; calibrate with weak laser response
tdcdft_screen_strength=1.0    ! s; zero is diagnostic-only, same propagation as fixed alpha
tdcdft_screen_floor=1e-8      ! threshold for sqrt(a.a + E.E/omega**2), atomic A/c units
```

横応答の古典場（XC場を含めない）を使い、a(t)=Ac_ext(t)-Ac_ext(0)、
E=-da/dt、P=-integral j dt とします。Pは台形積分で更新する3成分の状態で、
新たな時間履歴バッファはありません。既知の外場Eは中心差分、jとPは同じ端点で評価します。

```
K = (a.j - E.P) / (a.a + E.E/omega**2)
alpha_eff = alpha0 / [1 + s*4*pi*max(K-K0,0)/omega**2]
a_xc'' + beta*a_xc' + gamma*a_xc = alpha_eff*j
```

ベクトルの内積を使うスカラー縮約で、異方的な遮蔽テンソルではありません。
これは利用者の提案に基づく**現象論的な閉じ方**であり、Williams–Ullrich論文にある式ではありません。
実周波数のDrude誘電関数の逆数でも、有限qの電子–正孔相互作用から導出した式でもありません。
4*pi*Kの単位をプラズマ周波数の二乗に合わせていますが、分母を正に保つ形と強度sはモデル仮定です。
alphaを電流駆動項に掛ける定義を採用し、Exc=alpha(t)*Pを課す定義とは区別します。
後者から導けばalphaの時間微分が必要ですが、この実装はその式を採用していません。

最初は同じ周波数・包絡形状の弱レーザーでs=0の診断計算を行い、Kの基準を調べます。
単一K0ではバンド間応答の全時刻を再現できない場合があります。K0=0がSiの平衡応答という
意味ではありません。弱励起でalpha0を回復することは**校正して検証する条件**であり、
未校正のこの入力で自動的に保証されません。s=0では固定alphaと同じ伝播になります。

場のノルムがfloor以下ならKとalphaを保持し、Pの積分は継続します。
滑らかなパルス終端でも残留PがあるとKはfloor到達前に大きくなり得るため、
floorは単なる丸め誤差対策ではなくモデル結果に影響します。特にパルス後に保持されるalphaを
キャリア遮蔽の確定値と解釈してはいけません。floor、代表周波数、K0、sへの感度を調べてください。
広帯域・多色・遅延プローブでは単一周波数近似が弱くなります。プローブが推定値も変えるため、
過渡スペクトルの解釈は別途検証が必要です。精度改善や長時間安定性はまだ実証していません。

instantモードの `Si_rt_xc.data` は従来の7列に、alpha_eff、K、Pの3成分を追加した12列です。
追加列は原子単位。alpha_effはその時刻の電流から次のXC場を計算する際に使った値です。
checkpointはバージョン2でP,K,alphaと設定を保存します。旧バージョン1は固定alphaのみ再開可能。
instantの有無や係数変更での再開は拒否します。

## 分極と整合した瞬時遮蔽

`Si_pump_polarization.inp` の `tdcdft_screening='polarization'` は、上記の瞬時K推定と
alpha補正式を使いながら、XC場を `E_xc=alpha(t)*P` と定義する比較用モードです。
`instant` の `a_xc''=alpha(t)*j` とは異なります。beta=gamma=0を要求します。
既存のinstantモードは比較用に残してあります。

`a_xc'=-alpha*P` を二次精度の明示的な時間積分で進めます。
alphaの時間微分を後退差分で二次項に入れることで、`-alpha'*P` の寄与を含めます。
Pの台形積分と初期半ステップにより、alpha一定の極限では元の無減衰LRCと一致します。
出力E_xcはA_xcの中心差分なので、alpha変化中はalpha*Pとの差が離散化誤差として残ります。
alphaが一定になった後は、この差は丸め誤差の範囲で消えます。

時間精度は滑らかなalphaを仮定します。急なalpha変化や推定停止しきい値を横切る場合は、
時間刻みへの感度を検証してください。チェックポイントにはすでに保存するalphaとPを使い、
形式はバージョン2のままです。instant/polarization間で変更しての再開は拒否します。

これは残留電場の積分定数を避けるための**モデル定義の変更**であり、K推定のパルス終端感度や
遮蔽補正式そのものの物理的正しさを保証する修正ではありません。alphaの校正、しきい値感度、
強励起スペクトルの検証は引き続き必要です。外場終了時にA、J、Pをゼロへ強制する処理はありません。

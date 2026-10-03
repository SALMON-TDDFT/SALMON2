# 実験版の履歴数指定

```fortran
&functional
 exx_factor_history_frames=4
 exx_factor_exact_interval=12
 exx_factor_warmup_steps=192
/
```

framesは3:5、既定3。予測有効化は従来と同じSALMON_FACTOR_HISTORY=learned。厳密補正間隔12なら間の11stepを予測。warmupは補正間隔の整数倍、既存の10倍以上という入力制約を維持するが、4/5履歴では10倍ちょうどは8教師に届かず予測開始が遅れる。この比較では16倍で全条件8教師以上を確保。履歴数だけで予測区間を自動延長しない。

保持3/4/5は空間分布を3/4/5枚保持することで、学習期間3/4/5stepという意味ではない。学習係数は毎回の厳密教師で累積更新。warmup16intervalなら教師14/13/12組。予測モデル保存・GPU履歴拡張はこの実験では追加しない。

現段階はCPU実験機能。3/4/5履歴の比較継続中、10fsスペクトルの20%精度は未確定。GPUへの対応は本commitの対象外。

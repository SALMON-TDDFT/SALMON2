# RT比較サマリー（共通有効時間窓）

- 共通窓長 N_common = 600
- tは1列目、Jxは14列目（コメント行除外）
- HHG: Jxの平均除去 + Hann窓 + rFFT, h=omega/omega1 (omega1=0.05)
- 高調波比較は h=1..15 の整数格子へ補間して実施

> 注: full/dg/dgpw/dgpw_hi は run_fix23 を使用。

## 指標（FULL基準）
- FULL: RMSE=0.000000e+00, Pearson=1.000000, max|dJx|=0.000000e+00, HHG_L2=0.000000e+00
- DG: RMSE=5.286611e-06, Pearson=0.146538, max|dJx|=1.493190e-05, HHG_L2=3.425119e-01
- DGPW: RMSE=1.464467e-05, Pearson=-0.064944, max|dJx|=3.628650e-05, HHG_L2=3.190547e-01
- DGPWHI: RMSE=1.056206e-04, Pearson=-0.049928, max|dJx|=1.163976e-04, HHG_L2=7.333049e-01

## 上位5ピーク（h=1..15）
### FULL
- 1. h=11.0, intensity=1.000000e+00
- 2. h=10.0, intensity=9.869777e-01
- 3. h=12.0, intensity=9.754831e-01
- 4. h=13.0, intensity=9.509663e-01
- 5. h=9.0, intensity=9.319583e-01
### DG
- 1. h=15.0, intensity=1.000000e+00
- 2. h=14.0, intensity=9.746532e-01
- 3. h=13.0, intensity=9.493065e-01
- 4. h=12.0, intensity=9.239597e-01
- 5. h=11.0, intensity=8.986129e-01
### DGPW
- 1. h=10.0, intensity=1.000000e+00
- 2. h=11.0, intensity=9.847878e-01
- 3. h=9.0, intensity=9.478526e-01
- 4. h=12.0, intensity=9.093663e-01
- 5. h=8.0, intensity=8.957053e-01
### DGPWHI
- 1. h=11.0, intensity=1.000000e+00
- 2. h=10.0, intensity=9.976313e-01
- 3. h=12.0, intensity=9.186603e-01
- 4. h=9.0, intensity=9.016135e-01
- 5. h=13.0, intensity=8.373207e-01

## 近さランキング判定（FULLへの近さ）
- RMSE順: DG -> DGPW -> DGPWHI
- HHG L2順: DGPW -> DG -> DGPWHI

（判定は共通窓 N_common のみ）

# CPU学習ACE：履歴3固定

```fortran
&functional
 exx_factor_mode='learned'
 exx_factor_history_frames=3
 exx_factor_exact_interval=12
 exx_factor_n_teach=8
/
```

履歴は厳密な時間発展完了状態の因子3枚に固定。`exx_factor_history_frames`は既存入力との互換性のため残し、既定値・許可値とも3。4/5その他の指定は入力エラー。内部にも可変履歴数を保持せず、係数2個、学習行列2×2を使用する。

予測は`&functional`の`exx_factor_mode='learned'`で有効化する。環境変数による切替は廃止。厳密補正間隔12なら間の11stepを予測。warmupは `(H+N_teach-1)K` で自動計算する。H=3、N_teachの既定値8なので、この例のK=12では120stepとなる。`exx_factor_n_teach`で教師増分数を指定する。学習係数は厳密教師で累積更新し、指定数の教師と係数・整列の検査を満たしてから予測を許可する。旧`exx_factor_warmup_steps`は入力互換のため受理するが、自動計算値で置き換える。履歴3枚は学習期間3stepという意味ではない。

現段階はCPU実験機能。履歴4/5の比較は中止し、履歴3で予測区間と10fsスペクトルの20%精度を検証中。モデル保存・GPU対応はこの変更の対象外。履歴3への固定は予測区間を自動延長せず、全物質での精度も保証しない。

Review clarification (2026-10-03): full support (`exx_mlwf_radius=0`,
`exx_mlwf_norm_fraction=1`) permits an unshifted full k mesh with k-only
parallelism. Actual WF truncation remains restricted to Gamma; restart and
snapshot restrictions remain in effect. Internal prediction horizons use the
existing f/f² interpolation, which reproduces linear variation but has an
intermediate-time bias for quadratic variation, even with exact endpoint
coefficients. The time interpolation bias remains with the fixed three-frame model. The unit probe checks every horizon from 1 through 8.

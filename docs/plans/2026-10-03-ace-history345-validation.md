# CPU学習ACE公開版の検証記録

基点f8b01d74、SR-only full-support MLWF bypass版。学習ACE前提実装とk点毎履歴、frames3/4/5、exact_interval/warmup_steps namelistを含む。モデルは厳密accepted endpointだけを教師としPBE密度更新を保持。GPUとPTIM/Andersonは含めない。

- 標準GNU Fortran/MPI専用ビルド（既存object再利用）のsource/object/binarySHA保存済み。新moduleを含むリンク成功。
- 3/4/5履歴単体、OMP2、fcheck=allで合格。MPI2の整列・global内積合格。
- p3旧新回帰、horizon1:8で係数・予測max差3.16e-13。
- Si8MPI8OMP2の3/4/5履歴16stepは全rank正常。p3従来電流差2.4e-17、energy差5.7e-14Ha。
- 512step比較は実行中。間隔8/12は全rank正常、ACE拒否0。履歴追加の明確な改善なし。16と10fsは未確定。
- publish専用copyはimplicit none明記を追加、CPU unit CMake登録。公開copyの単体を再コンパイルし合格。
- 実行中の重い計算を重ねないため全ソースclean buildはこの公開turnでは未実施。GPU、全構成、保存モデル再利用の保証はしない。

学習ACEはSALMON_FACTOR_HISTORY=learnedで明示有効化する実験機能。デフォルトinterval8/warmup128/frames3は現行互換値、一般物質の推奨精度を意味しない。有限振幅比例性未確認。スペクトル20%基準は長時間検証で判定する。

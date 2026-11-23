# CIの重複整理（2026-11-23）

## 起動と検査範囲

| 入口 | 維持する検査 |
| --- | --- |
| push | Linuxのsdist→wheel・未修復ライブラリのハッシュ・導入後全テスト |
| 関連ソース等のpush/PR | Git追跡ファイルのみの配布元検査、MPIライブラリ構築、GS不要のC/Fortran照会・入力検査 |
| wheelのPR/手動 | 同一sdistから4環境×Python3.10～3.14を構築・修復・導入。全版でLD/SBE・公開API・AI用例の検査 |
| 上記のPython3.12 | 各OSで全テストとASE/plot依存を1回。別の依存最小環境ではcore検査 |
| runtime手動 | `scope=worker` はLinuxの短いワーカー検査、`scope=bundle` は3環境の配布・結合検査 |
| NLCC手動 | 検証用Libxcを使う固有の検査。設定の異なるビルドは安易に共用しない |

PR/手動wheelの全テスト実行は25回→4回です。残る版もcore検査を実行し、20種類のwheel構築は維持します。
84%は全テストの**呼び出し回数**の削減で、料金・Actions時間の削減率ではありません。
共通sdistの検査に失敗した場合は、後続の4環境を起動しません。

単独のtddft-worker.ymlを削除し、その検査はruntime-bundles.ymlへ統合しました。
短いworker選択では1回のLinuxビルドだけで済み、bundleの全環境を強制しません。
移動時に準備テストの古い二引数呼び出しも現行APIへ修正しました。
新しいCIポリシーテストは、対象OS/Python版、必須検査、空のテスト選択、手動限定の物理検査を確認します。
共有カーネル・ソース設定・監査／ビルド補助・ABI試験データの変更でも、必要なMPI検査が起動するようにしました。

同じworkflow・ブランチ／PRの古い自動実行は、新しい実行で取り消します。
明示的な手動実行は個別グループとし、取り消しません。既存のリモートジョブの停止・再実行は行っていません。
pushとPRは検査するGitの内容が異なり得るため両方残します。
ソース／オブジェクトの禁止固有値ソルバー監査、未修復wheelの検査、配布後の検査は省略しません。

## 確認範囲

変更前の同じコミット77a7f004について、fork元の
[LD wheels](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/45115812759) と
[MPI照会・入力検査](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/45115812843) は成功を確認しました。
ジョブ経過時間はそれぞれ約68秒・175秒で、前者は配布用の軽いライブラリ、後者はMPI用ソルバーを構築します。
独立privateリポジトリの履歴は接続権限で参照できなかったため、実行例はfork元の同SHAを参照しています。

新構成はactionlintとCIポリシーテスト5件が通過しました。cibuildwheelの実際の設定読込みで20識別子とfull/coreの振り分けを確認しました。
Macでは実際のsdist→wheelを構築し、未修復wheelのハッシュ検査が通過しました。導入後の全711件は710件成功・1件スキップ、coreは51件成功、NumPyだけの環境では49件成功・2件スキップでした。ランタイム準備テスト2件も成功しています。
これは新構成のGitHub上の実行成功ではありません。重いCIの手動起動・GS計算は行いません。
詳細記録はcalc2026/salmonpy-ci-consolidation-20261123へ保存します。

## 設定の一次資料

- [GitHubの同時実行制御](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/control-workflow-concurrency)
- [GitHubのYAMLアンカー](https://docs.github.com/en/actions/reference/workflows-and-actions/reusing-workflow-configurations)
- [cibuildwheelの版別設定](https://cibuildwheel.pypa.io/en/v3.4.1/configuration/)
- [cibuildwheel actionの入力](https://github.com/pypa/cibuildwheel/blob/v3.4.0/action.yml)

# 日本語ドキュメント完全ガイド

**SALMON v2.2.2 DG-Fragment RT実装 - 日本語リソースマップ**

---

## 📚 利用可能な日本語ドキュメント

### 1. **DG_FRAGMENT_RT_SESSION23_FINAL_JP.md** ← おすすめ出発点
- **内容**: セッション23完了レポート（日本語最初の選択肢）
- **対象**: プロジェクト管理者、実装完了状況確認
- **分量**: 約400行
- **主要セクション**:
  - エグゼクティブサマリー
  - セッション23実施内容
  - 物理検証ステータス
  - 技術実装詳細

### 2. **DG_FRAGMENT_TECHNICAL_SPEC_JP.md** ← 技術者向け
- **内容**: 完全技術仕様書（日本語）
- **対象**: コード実装者、デバッグが必要な場合
- **分量**: 約500行
- **主要セクション**:
  - システム概要
  - 主要パラメータ詳細
  - メモリアクセス修正（深掘り）
  - 実装詳細（関数インタフェース）
  - テスト・検証方法
  - トラブルシューティング

### 3. **英語ドキュメント参照**

#### 基本・全般向け
- **DG_FRAGMENT_RT_STATUS.md**
  - プロジェクト全体状況（英語）
  - クイックスタートガイド
  - パフォーマンス特性

#### 詳細レポート
- **PHYSICS_VALIDATION_REPORT.md**
  - 物理検証詳細
  - SIGBUS修正の根本原因分析
  - テスト設定と期待値

#### セッション記録
- **SESSION_23_SUMMARY.md**
  - セッション別実施内容
  - コード変更の詳細

#### ナビゲーション
- **DOCUMENTATION_INDEX.md**
  - すべてのドキュメントへのリンク
  - トピック別クイックナビ

---

## 🎯 用途別ドキュメント選択ガイド

### シナリオ1: 「今何ができるか知りたい」

**読むドキュメント** (優先順):
1. DG_FRAGMENT_RT_SESSION23_FINAL_JP.md → 「完了チェックリスト」セクション
2. DG_FRAGMENT_TECHNICAL_SPEC_JP.md → 「テスト・検証」セクション

**時間**: 10分

### シナリオ2: 「SIGBUSエラーの修正内容を理解したい」

**読むドキュメント**:
1. DG_FRAGMENT_TECHNICAL_SPEC_JP.md → 「メモリアクセス修正」セクション
2. PHYSICS_VALIDATION_REPORT.md → 「Root Cause Analysis」セクション（英語）
3. DG_FRAGMENT_TECHNICAL_SPEC_JP.md → 「実装詳細」セクション

**時間**: 20分

### シナリオ3: 「テストを実行したい」

**実施ステップ**:
```bash
cd /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test/H2

# テスト実行（自動化）
bash run_physics_validation.sh
```

**参考ドキュメント**:
- DG_FRAGMENT_RT_SESSION23_FINAL_JP.md → 「本番テスト実施方法」
- DG_FRAGMENT_TECHNICAL_SPEC_JP.md → 「テスト・検証」

**時間**: テスト実行 5-10分 + 分析 5分

### シナリオ4: 「問題が発生。何をすればいい?」

**トラブルシューティング**:
1. DG_FRAGMENT_TECHNICAL_SPEC_JP.md → 「トラブルシューティング」
2. DG_FRAGMENT_RT_SESSION23_FINAL_JP.md → 「技術知見」
3. 必要に応じて英語ドキュメント参照

**時間**: 20-30分

### シナリオ5: 「コードを修正・拡張したい」

**事前準備**:
1. DG_FRAGMENT_TECHNICAL_SPEC_JP.md → 「実装詳細」（完全読破）
2. コードを確認: src/rt/rt_dg_fragment.f90
3. 関連関数理解: read_fragment_basis_data, calculate_observables

**参考リソース**:
- DG_FRAGMENT_TECHNICAL_SPEC_JP.md → 「関数インタフェース」
- コード内のコメント
- 英語ドキュメント（詳細は英語に記載）

**時間**: コード理解 1-2時間

---

## 📖 各ドキュメント詳細説明

### DG_FRAGMENT_RT_SESSION23_FINAL_JP.md

**適した読者**:
- プロジェクト管理者
- 実装完了確認が必要な人
- 全体状況を把握したい人

**構成**:
```
1. エグゼクティブサマリー（5分読）
2. セッション23での実施内容（10分読）
3. 物理検証ステータス（5分読）
4. コード修正サマリー（10分読）
5. 本番テスト実施方法（5分読）
6. 重要技術知見（10分読）
```

**活用例**:
- クライアントへの進捗報告
- チームミーティングでの説明資料
- 全体的な完成度確認

---

### DG_FRAGMENT_TECHNICAL_SPEC_JP.md

**適した読者**:
- コード実装者・デバッガー
- 深い技術理解が必要な人
- 拡張・修正を計画している人

**構成**:
```
1. システム概要（画像付き流れ図）
2. 主要パラメータ（詳細説明 + 相互作用表）
3. メモリアクセス修正（修正前後コード比較）
4. 実装詳細（関数シグネチャ + 補足）
5. テスト・検証（実行方法 + 期待値）
6. トラブルシューティング（Q&A形式）
```

**活用例**:
- コード修正時の参考書
- 新しい機能追加時のガイド
- メモリ問題発生時の診断
- パフォーマンス最適化時の指標

---

## 🔄 ドキュメント間の関連性

```
DG_FRAGMENT_RT_SESSION23_FINAL_JP.md （全体概観）
  ↓
  ├─→ DG_FRAGMENT_TECHNICAL_SPEC_JP.md （詳細仕様）
  │     ├─→ コード内コメント
  │     └─→ 英語ドキュメント（詳細説明）
  │
  └─→ 実行スクリプト
     ├─→ run_physics_validation.sh
     └─→ physics_validation.py

実装時の参照順序:
  DG_FRAGMENT_RT_SESSION23_FINAL_JP.md
    ↓
  DG_FRAGMENT_TECHNICAL_SPEC_JP.md
    ↓
  実装ファイル: src/rt/rt_dg_fragment.f90
    ↓
  英語リファレンス: DG_FRAGMENT_RT_STATUS.md
```

---

## 💡 日本語ドキュメント活用Tips

### Tip 1: 素早く全体を理解する（15分）

1. DG_FRAGMENT_RT_SESSION23_FINAL_JP.md
   - エグゼクティブサマリー (2分)
   - マイルストーン達成表 (1分)
   - ステータステーブル (2分)

2. DG_FRAGMENT_TECHNICAL_SPEC_JP.md
   - システム概要の図 (3分)
   - パラメータテーブル (3分)
   - メモリ修正パターン (4分)

### Tip 2: 詳しく学ぶ（45分）

1. DG_FRAGMENT_RT_SESSION23_FINAL_JP.md （20分）
   - すべてのセクション読破
   
2. DG_FRAGMENT_TECHNICAL_SPEC_JP.md （25分）
   - すべてのセクション読破
   - コード例を検証

### Tip 3: 問題解決時に使う（10-20分）

DG_FRAGMENT_TECHNICAL_SPEC_JP.md の
「トラブルシューティング」セクションを確認
→ 問題パターン別対応が記載

### Tip 4: コード修正時に使う（参照用）

- DG_FRAGMENT_TECHNICAL_SPEC_JP.md の「実装詳細」
- 関数シグネチャとインタフェース
- コード例を確認
- 必要に応じて英語ドキュメント参照

---

## 📋 ドキュメント チェックリスト

### 読むべき順序（第1回）

- [ ] DG_FRAGMENT_RT_SESSION23_FINAL_JP.md 全体
- [ ] DG_FRAGMENT_TECHNICAL_SPEC_JP.md のシステム概要
- [ ] DG_FRAGMENT_TECHNICAL_SPEC_JP.md のメモリ修正
- [ ] 納得できるまで再読 ✓

### 実装作業前に確認

- [ ] DG_FRAGMENT_TECHNICAL_SPEC_JP.md の「実装詳細」セクション
- [ ] src/rt/rt_dg_fragment.f90 のコード確認
- [ ] 関連関数の理解確認

### テスト実行前に確認

- [ ] DG_FRAGMENT_RT_SESSION23_FINAL_JP.md の「本番テスト実施方法」
- [ ] run_physics_validation.sh スクリプト確認
- [ ] 期待出力の理解

### 問題発生時に確認

- [ ] DG_FRAGMENT_TECHNICAL_SPEC_JP.md の「トラブルシューティング」
- [ ] エラーメッセージの照合
- [ ] 対応方法の実施

---

## 🌐 英語ドキュメントとの使い分け

### 日本語が詳しい項目

- ✅ システム概要
- ✅ パラメータ説明
- ✅ メモリ修正の詳細
- ✅ テスト実施方法
- ✅ 初期設定手順

### 英語が詳しい項目

- 🔗 物理検証の詳細（PHYSICS_VALIDATION_REPORT.md）
- 🔗 パフォーマンス特性（DG_FRAGMENT_RT_STATUS.md）
- 🔗 長期拡張計画（DG_FRAGMENT_RT_STATUS.md）

### 推奨: 両言語を参照する

特に以下の場合:
1. 完全理解が必要な重要な修正
2. パフォーマンス問題対応時
3. 新しい機能追加時

---

## ✨ まとめ

### 最初の5分

```
DG_FRAGMENT_RT_SESSION23_FINAL_JP.md
    → 「エグゼクティブサマリー」を読む
    ✓ プロジェクトの現状がわかる
```

### 次の15分

```
DG_FRAGMENT_RT_SESSION23_FINAL_JP.md
    → 「セッション23での実施内容」を読む
    ✓ 何が修正されたか理解
```

### さらに詳しく（20分）

```
DG_FRAGMENT_TECHNICAL_SPEC_JP.md
    → 「メモリアクセス修正」セクション
    ✓ 修正の詳細を学ぶ
```

### テスト実行（10分）

```
bash run_physics_validation.sh
    ✓ 実装が正常か確認
```

---

## 🎓 学習パス（推奨）

**初級（初めて接する）**:
1. DG_FRAGMENT_RT_SESSION23_FINAL_JP.md （30分）
2. DG_FRAGMENT_TECHNICAL_SPEC_JP.md 「システム概要」 （10分）
3. テスト実行 （10分）

**中級（詳しく理解したい）**:
1. 初級のすべて
2. DG_FRAGMENT_TECHNICAL_SPEC_JP.md 「メモリアクセス修正」 （20分）
3. DG_FRAGMENT_TECHNICAL_SPEC_JP.md 「実装詳細」 （30分）
4. コード確認 （30分）

**上級（修正・拡張する）**:
1. 中級のすべて
2. DG_FRAGMENT_TECHNICAL_SPEC_JP.md 全セクション （60分）
3. 英語ドキュメント参照 （必要に応じて）
4. コード実装 （60分以上）

---

**最終更新**: 2026年2月23日  
**推奨読了時間総計**: 60-90分（基本習得）


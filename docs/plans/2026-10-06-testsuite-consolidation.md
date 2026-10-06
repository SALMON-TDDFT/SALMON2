# Testsuite consolidation implementation plan

> 2026-10-06：以下の`benchmarks`は当時のローカル性能測定です。マージ対象から除外し、測定コードはブランチ外へ保存しています。

Goal: Reduce duplicate PBEh40 integration cases and give retained test directories numeric identifiers.

1. Remove 429, 433, 434; retain DC PBEh40, DC/conventional PBEh40+rVV10, and one Ehrenfest case.
2. Move performance measurements to benchmarks; retain numerical unit tests under numbered directories at unchanged depth.
3. Update source/document/script references, check Python syntax and CMake fixture graph, and run registered mathematical tests.
4. Keep shared pseudo data and existing calculation outputs untouched. No commit/push requested.

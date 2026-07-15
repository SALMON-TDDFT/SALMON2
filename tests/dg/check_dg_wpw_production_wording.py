#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
paths = [ROOT / "docs/plans/2026-07-12-wannier-pw-exp-full-tddft.md",
         ROOT / "docs/plans/2026-07-12-wannier-pw-exp-full-tddft-design.md"]
text = "\n".join(path.read_text().lower() for path in paths)
for stale in ("first production target is the global-ownership",
              "only after the global path passes",
              "then plan distributed ownership as a separate milestone",
              "transferred to the fragment/distributed implementation"):
    assert stale not in text, f"stale global-first production wording: {stale}"
for required in ("windowed_kg", "legacy g-only", "small-system oracle", "matrix-free"):
    assert required in text, f"missing production/reference distinction: {required}"
print("PASS production wording separates scalable windowed_kg from G-only reference")

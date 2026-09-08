#!/usr/bin/env python3
"""Mutation contract for the non-self-referential Task 10 evidence manifest."""

from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VERIFIER = ROOT / "tests/dg/verify_task10_wf_reuse_evidence.py"
SPEC = importlib.util.spec_from_file_location("task10_evidence_verifier", VERIFIER)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


with tempfile.TemporaryDirectory(prefix="task10-manifest-contract-") as name:
    root = Path(name)
    scope = root / "raw"
    scope.mkdir()
    (scope / "a.log").write_text("complete raw output\n")
    nested = scope / "checkpoint"
    nested.mkdir()
    (nested / "rank-000000.bin").write_bytes(b"rank zero payload")
    individual = root / "analyzer.py"
    individual.write_text("print('analyzer')\n")
    manifest_path = root / "task10-manifest.json"

    manifest = MODULE.build_manifest(
        repository_head="abc123",
        scopes=[("raw", scope)],
        individual_files=[("analyzer", individual)],
        validation={"status": "PASS"},
        excluded_paths={manifest_path.resolve()},
    )
    MODULE.write_manifest(manifest_path, manifest)
    MODULE.verify_manifest(manifest_path, expected_repository_head="abc123")
    serialized = json.loads(manifest_path.read_text())
    assert all(
        item["path"] != str(manifest_path.resolve())
        for item in serialized["individual_files"]
    )
    assert all(
        entry["path"] != str(manifest_path.resolve())
        for scope_entry in serialized["scopes"]
        for entry in scope_entry["files"]
    )

    (scope / "a.log").write_text("X" * len("complete raw output\n"))
    try:
        MODULE.verify_manifest(manifest_path, expected_repository_head="abc123")
    except RuntimeError as error:
        assert "hash" in str(error).lower()
    else:
        raise AssertionError("manifest accepted mutated raw evidence")
    (scope / "a.log").write_text("complete raw output\n")

    (scope / "unlisted.log").write_text("late artifact\n")
    try:
        MODULE.verify_manifest(manifest_path, expected_repository_head="abc123")
    except RuntimeError as error:
        assert "inventory" in str(error).lower()
    else:
        raise AssertionError("manifest accepted a changed scope inventory")

print("Task10 non-self-referential evidence manifest contract: PASS")

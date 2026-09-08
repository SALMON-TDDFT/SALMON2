#!/usr/bin/env python3
"""Validate Task 10 raw receipts and create/verify an exact SHA-256 inventory."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import run_dg_fragment_wf_production_smoke as smoke
import run_dg_hybrid_si64_scdm_reuse as si64


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, *, label: str | None = None, root: Path | None = None) -> dict[str, Any]:
    absolute = path.resolve(strict=True)
    if not absolute.is_file():
        raise RuntimeError(f"manifest artifact is not a regular file: {absolute}")
    record: dict[str, Any] = {
        "path": str(absolute),
        "size": absolute.stat().st_size,
        "sha256": sha256(absolute),
    }
    if label is not None:
        record["label"] = label
    if root is not None:
        record["relative_path"] = absolute.relative_to(root.resolve(strict=True)).as_posix()
    return record


def scope_files(root: Path, excluded_paths: set[Path]) -> list[Path]:
    absolute_root = root.resolve(strict=True)
    if not absolute_root.is_dir():
        raise RuntimeError(f"manifest scope is not a directory: {absolute_root}")
    return sorted(
        (path for path in absolute_root.rglob("*")
         if path.is_file() and path.resolve() not in excluded_paths),
        key=lambda path: path.relative_to(absolute_root).as_posix(),
    )


def build_manifest(
    *,
    repository_head: str,
    scopes: list[tuple[str, Path]],
    individual_files: list[tuple[str, Path]],
    validation: dict[str, Any],
    excluded_paths: set[Path] | None = None,
) -> dict[str, Any]:
    """Build an exact inventory without ever including the manifest itself."""
    excluded = {path.resolve() for path in (excluded_paths or set())}
    scope_records = []
    recorded_paths: set[Path] = set()
    for label, root in scopes:
        absolute_root = root.resolve(strict=True)
        files = scope_files(absolute_root, excluded)
        records = []
        for path in files:
            absolute = path.resolve(strict=True)
            if absolute in recorded_paths:
                raise RuntimeError(f"artifact appears in overlapping manifest scopes: {absolute}")
            recorded_paths.add(absolute)
            records.append(file_record(absolute, root=absolute_root))
        scope_records.append({"label": label, "root": str(absolute_root), "files": records})
    individual_records = []
    for label, path in individual_files:
        absolute = path.resolve(strict=True)
        if absolute in excluded:
            raise RuntimeError(f"excluded manifest path requested as an artifact: {absolute}")
        if absolute in recorded_paths:
            raise RuntimeError(f"artifact appears in both a scope and individual inventory: {absolute}")
        recorded_paths.add(absolute)
        individual_records.append(file_record(absolute, label=label))
    return {
        "format": "SALMON_TASK10_SHA256_INVENTORY_V1",
        "hash_algorithm": "sha256",
        "repository_head": repository_head,
        "validation": validation,
        "excluded_paths": sorted(str(path) for path in excluded),
        "scopes": scope_records,
        "individual_files": individual_records,
        "artifact_count": sum(len(scope["files"]) for scope in scope_records)
        + len(individual_records),
    }


def write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def verify_record(record: dict[str, Any]) -> None:
    path = Path(record["path"])
    if not path.is_file():
        raise RuntimeError(f"manifest artifact is missing: {path}")
    if path.stat().st_size != record["size"]:
        raise RuntimeError(f"manifest artifact size changed: {path}")
    if sha256(path) != record["sha256"]:
        raise RuntimeError(f"manifest artifact hash changed: {path}")


def verify_manifest(path: Path, *, expected_repository_head: str) -> dict[str, Any]:
    manifest_path = path.resolve(strict=True)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("format") != "SALMON_TASK10_SHA256_INVENTORY_V1":
        raise RuntimeError("unknown Task10 evidence manifest format")
    if manifest.get("repository_head") != expected_repository_head:
        raise RuntimeError("Task10 manifest repository HEAD changed")
    if str(manifest_path) not in manifest.get("excluded_paths", []):
        raise RuntimeError("Task10 manifest does not explicitly exclude itself")
    excluded = {Path(value).resolve() for value in manifest["excluded_paths"]}
    counted = 0
    for scope in manifest["scopes"]:
        root = Path(scope["root"]).resolve(strict=True)
        expected_relatives = [record["relative_path"] for record in scope["files"]]
        actual_relatives = [
            item.resolve().relative_to(root).as_posix()
            for item in scope_files(root, excluded)
        ]
        if actual_relatives != expected_relatives:
            raise RuntimeError(f"manifest scope inventory changed: {scope['label']}")
        for record in scope["files"]:
            verify_record(record)
            counted += 1
    for record in manifest["individual_files"]:
        verify_record(record)
        counted += 1
    if counted != manifest.get("artifact_count"):
        raise RuntimeError("manifest artifact count changed")
    return manifest


def elapsed_from_files(run_dir: Path) -> float:
    return max(
        0.0,
        (run_dir / "overlapping_wannier_occupied.chk").stat().st_mtime
        - (run_dir / "inputfile").stat().st_mtime,
    )


def validate_raw_evidence(si64_result: Path, smoke_result: Path) -> dict[str, Any]:
    si64_clean_dir = si64_result / "01-clean-scdm"
    si64_reuse_dir = si64_result / "02-exact-reuse"
    si64_clean = si64.parse_run(
        si64_clean_dir.name, si64_clean_dir, elapsed_from_files(si64_clean_dir), None)
    si64_reuse = si64.parse_run(
        si64_reuse_dir.name, si64_reuse_dir, elapsed_from_files(si64_reuse_dir), None)
    si64_comparison = si64.compare_generation_and_reuse(si64_clean, si64_reuse)
    gauge = si64.fragment_wf_contract_receipts(
        si64_result / "fragment-wf-checkpoint", si64_clean["wf_publication_id"])

    smoke_dirs = [smoke_result / name for name in (
        "01-miss", "02-hit", "03-incomplete-recovery")]
    smoke_runs = [smoke.parse_evidence(
        directory.name, directory, elapsed_from_files(directory), None)
        for directory in smoke_dirs]
    seed_preloaded = smoke_runs[0]["dc_scf_skipped"]
    smoke.compare_smoke_runs(*smoke_runs, seed_preloaded)

    occupied = si64_comparison["occupied_state_numeric_comparison"]
    return {
        "status": "PASS",
        "si64": {
            "mpi_ranks": si64.EXPECTED_MPI_SIZE,
            "completion_receipts": [run["completion_receipts"] for run in (si64_clean, si64_reuse)],
            "scdm_rank_fragment_pairs": [[item["rank"], item["fragment_id"]] for item in gauge],
            "clean_wannier90_fragment_ids": [item["fragment_id"] for item in si64_clean["wannier90"]],
            "reuse_wannier90_fragment_ids": [item["fragment_id"] for item in si64_reuse["wannier90"]],
            "terminal_residuals": [si64_clean["terminal_residual"], si64_reuse["terminal_residual"]],
            "terminal_electron_defects": [
                si64_clean["terminal_electron_defect"], si64_reuse["terminal_electron_defect"]],
            "occupied_density_relative_frobenius_difference":
                occupied["occupied_density_relative_frobenius_difference"],
            "exact_fingerprint_fields": si64_comparison["exact_fingerprint_fields"],
        },
        "smoke": {
            "mpi_ranks": smoke.RANKS,
            "completion_receipts": [run["completion_receipts"] for run in smoke_runs],
            "wannier90_fragment_ids": [run["wannier_fragment_ids"] for run in smoke_runs],
            "projected_basis_fingerprints": [run["projected_basis_fingerprint"] for run in smoke_runs],
            "continuation_receipts_identical": smoke_runs[0]["schwarz_receipts"]
                == smoke_runs[1]["schwarz_receipts"] == smoke_runs[2]["schwarz_receipts"],
            "terminal_residuals": [run["terminal_residual"] for run in smoke_runs],
            "terminal_electron_defects": [run["electron_defect"] for run in smoke_runs],
        },
        "long_si64_rerun_required": False,
        "long_si64_rerun_reason": (
            "hardened preserved-Si64 analysis and current integrated-binary smoke both pass"
        ),
    }


def extract_si64_seed_directory(input_path: Path) -> Path:
    text = input_path.read_text()
    match = re.search(r"^\s*dg_dc_seed_directory\s*=\s*['\"]([^'\"]+)['\"]", text, re.I | re.M)
    if not match:
        raise RuntimeError("Si64 input does not identify its authoritative DC seed directory")
    return Path(match.group(1)).resolve(strict=True)


def git_head(repository: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repository, check=True,
        text=True, stdout=subprocess.PIPE).stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--si64-result-dir", type=Path, required=True)
    parser.add_argument("--smoke-result-dir", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--evidence-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    args = parser.parse_args()

    repository = args.repository.resolve(strict=True)
    si64_result = args.si64_result_dir.resolve(strict=True)
    smoke_result = args.smoke_result_dir.resolve(strict=True)
    build_dir = args.build_dir.resolve(strict=True)
    evidence_root = args.evidence_root.resolve(strict=True)
    manifest_path = args.manifest.resolve()
    head = git_head(repository)
    validation = validate_raw_evidence(si64_result, smoke_result)

    if args.write:
        si64_seed = extract_si64_seed_directory(si64_result / "01-clean-scdm/inputfile")
        scopes = [
            ("si64_clean_raw", si64_result / "01-clean-scdm"),
            ("si64_reuse_raw", si64_result / "02-exact-reuse"),
            ("si64_fragment_wf_publication", si64_result / "fragment-wf-checkpoint"),
            ("si64_authoritative_dc_seed", si64_seed),
            ("small_smoke_raw", smoke_result),
        ]
        individual_files = [
            ("si64_analyzer", repository / "tests/dg/run_dg_hybrid_si64_scdm_reuse.py"),
            ("smoke_analyzer", repository / "tests/dg/run_dg_fragment_wf_production_smoke.py"),
            ("analyzer_mutation_contract", repository / "tests/dg/check_dg_hybrid_si64_scdm_runner.py"),
            ("manifest_mutation_contract", repository / "tests/dg/check_task10_wf_reuse_manifest.py"),
            ("manifest_verifier", Path(__file__)),
            ("integrated_binary", build_dir / "salmon"),
            ("cmake_cache", build_dir / "CMakeCache.txt"),
            ("task10_build_log", evidence_root / "task10/build.stdout.txt"),
            ("task8_build_log", evidence_root / "task8-build.log"),
            ("task9_build_log", evidence_root / "task9/final-production-build.log"),
            ("status_before", evidence_root / "status-before.txt"),
            ("status_after", evidence_root / "status-after.txt"),
            ("dirty_manifest_before", evidence_root / "dirty-content-manifest-before.txt"),
            ("dirty_manifest_after", evidence_root / "dirty-content-manifest-after.txt"),
            ("preservation_manifest_before", evidence_root / "preservation-manifest-before.txt"),
            ("preservation_manifest_after", evidence_root / "preservation-manifest-after.txt"),
            ("task2_evidence_checksums", evidence_root / "task2-evidence-checksums.sha256"),
            ("task9_original_content_preservation", evidence_root / "task9/original-content-preservation.txt"),
            ("task9_original_status", evidence_root / "task9/original-status-current.txt"),
        ]
        manifest = build_manifest(
            repository_head=head,
            scopes=scopes,
            individual_files=individual_files,
            validation=validation,
            excluded_paths={manifest_path},
        )
        write_manifest(manifest_path, manifest)
    manifest = verify_manifest(manifest_path, expected_repository_head=head)
    print(json.dumps({
        "status": "PASS",
        "repository_head": head,
        "artifact_count": manifest["artifact_count"],
        "manifest": str(manifest_path),
        "validation": validation,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

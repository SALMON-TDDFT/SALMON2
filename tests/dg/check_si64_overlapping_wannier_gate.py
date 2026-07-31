#!/usr/bin/env python3
"""Validate the raw Si64 overlapping-Wannier ground-state evidence matrix."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import math
import re
import subprocess
import sys
from pathlib import Path


DECOMPOSITIONS = ("2x2x2",)
assert DECOMPOSITIONS == ("2x2x2",), "Task 9 accepts only the fixed 2x2x2 physical decomposition"
BOX_PROFILES = {
    "buffer5": (5, 5, 5),
    "buffer6": (6, 6, 6),
}
PRODUCTION_BOX_PROFILE = "buffer5"
WINDOWS = ("c192-sp48", "c256-sp48")
FORBIDDEN = (
    "[DG-DC-DIRECT]",
    "[DC-LCFO",
    "EigenExa",
    "[DG-WPW",
    "checkpoint data is printed",
    "writing restart data",
    "time evolution",
)
FLOAT_KEYS = {
    "runtime_s",
    "peak_memory_mib",
    "metric_min",
    "metric_max",
    "metric_condition",
    "occupied_inclusion_residual",
    "complete_shell_inclusion_residual",
    "center_defect",
    "spread_max",
    "tail_value_norm",
    "tail_gradient_norm",
    "s_hermiticity",
    "h_hermiticity",
    "t_hermiticity",
    "vnl_hermiticity",
    "vlocal_hermiticity",
    "density_residual",
    "unmixed_density_residual",
    "coefficient_residual",
    "s_orthogonality",
    "trace_charge",
    "integrated_charge",
    "total_energy",
}
INTEGER_KEYS = {
    "mpi_ranks",
    "omp_threads",
    "candidate_per_fragment",
    "target_per_fragment",
    "occupied_included",
    "occupied_required",
    "core_atoms_per_fragment",
    "global_target",
    "complete_shell_channels",
    "complete_shell_residual_rank",
    "bond_center_orbit_closed",
    "checkpoint_format_version",
    "matrix_owned_rows_max",
    "overlap_local_bytes_max",
    "hamiltonian_local_bytes_max",
}
TOLERANCES = {
    "center_defect": 1.0e-8,
    "occupied_inclusion_residual": 1.0e-8,
    "complete_shell_inclusion_residual": 1.0e-8,
    "s_hermiticity": 1.0e-10,
    "h_hermiticity": 1.0e-10,
    "t_hermiticity": 1.0e-10,
    "vnl_hermiticity": 1.0e-10,
    "vlocal_hermiticity": 1.0e-10,
    "density_residual": 1.0e-7,
    "unmixed_density_residual": 1.0e-7,
    "coefficient_residual": 1.0e-7,
    "s_orthogonality": 1.0e-6,
}
CONVERGENCE_KEYS = (
    "metric_min",
    "metric_max",
    "s_hermiticity",
    "h_hermiticity",
    "t_hermiticity",
    "vnl_hermiticity",
    "vlocal_hermiticity",
    "density_residual",
    "coefficient_residual",
    "integrated_charge",
    "total_energy",
)
TAG = re.compile(r"^\[OW-GS-EVIDENCE\]\s+([a-z0-9_]+)=\s*(\S+)\s*$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
DC_SCF = re.compile(
    r"DC\s+#SCF\s*=\s*(\d+)\s+Total Energy\s*=\s*([+\-0-9.EeDd]+)\s+diff\s*=\s*([+\-0-9.EeDd]+)"
)


class GateFailure(RuntimeError):
    pass


def fail(message: str) -> None:
    raise GateFailure(message)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_evidence(log_path: Path) -> tuple[dict[str, float], dict[str, int], dict[str, str]]:
    if not log_path.is_file():
        fail(f"missing raw log: {log_path}")
    text = log_path.read_text(errors="replace")
    for marker in FORBIDDEN:
        if marker.lower() in text.lower():
            fail(f"forbidden route marker {marker!r} in {log_path}")
    raw: dict[str, str] = {}
    for line in text.splitlines():
        match = TAG.match(line.strip())
        if match:
            key, value = match.groups()
            if key in raw:
                fail(f"duplicate evidence key {key!r} in {log_path}")
            raw[key] = value
    missing = (FLOAT_KEYS | INTEGER_KEYS | {"checkpoint_manifest_sha256"}) - raw.keys()
    if missing:
        fail(f"missing evidence keys in {log_path}: {', '.join(sorted(missing))}")
    floats = {key: float(raw[key].replace("D", "E").replace("d", "e")) for key in FLOAT_KEYS}
    integers = {key: int(raw[key]) for key in INTEGER_KEYS}
    if not all(math.isfinite(value) for value in floats.values()):
        fail(f"non-finite evidence in {log_path}")
    return floats, integers, raw


def validate_restart_log(log_path: Path) -> None:
    if not log_path.is_file():
        fail(f"missing restart log: {log_path}")
    text = log_path.read_text(errors="replace")
    lowered = text.lower()
    for marker in FORBIDDEN:
        if marker.lower() in lowered:
            fail(f"forbidden route marker {marker!r} in {log_path}")
    if text.count("[OW-GS] reused accepted route checkpoint") != 1:
        fail(f"{log_path}: checkpoint reuse must be reported exactly once")
    if "[ow-scf-diagnostic]" in lowered:
        fail(f"{log_path}: restart recomputed the overlapping-Wannier SCF")
    if "error stop" in lowered or "failed" in lowered:
        fail(f"{log_path}: restart contains an error or failed gate")


def validate_fixed_decomposition(variables: str, context: str) -> None:
    match = re.search(
        r"num_fragment\s*=\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)",
        variables,
        re.I,
    )
    if not match or tuple(map(int, match.groups())) != (2, 2, 2):
        fail(f"{context}: runtime physical decomposition is not 2x2x2")


def validate_row(root: Path, decomposition: str, box_profile: str, window: str) -> dict[str, float]:
    row_name = f"decomp-{decomposition}_box-{box_profile}_window-{window}"
    row = root / row_name
    metrics, integers, raw = parse_evidence(row / "run.log")
    variables = (row / "variables.log").read_text(errors="replace")
    validate_fixed_decomposition(variables, row_name)
    match = re.search(r"num_rgrid_buffer\s*=\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)", variables, re.I)
    if not match or tuple(map(int, match.groups())) != BOX_PROFILES[box_profile]:
        fail(f"{row_name}: runtime buffer vector does not match the physical-box profile")
    if integers["mpi_ranks"] != 8 or integers["omp_threads"] != 1:
        fail(f"{row_name}: requires exactly 8 MPI ranks and 1 OpenMP thread")
    candidate = int(window.split("-")[0][1:])
    if integers["candidate_per_fragment"] != candidate or integers["target_per_fragment"] != 48:
        fail(f"{row_name}: runtime window does not match matrix coordinate")
    if integers["core_atoms_per_fragment"] != 8 or integers["global_target"] != 384:
        fail(f"{row_name}: occupied plus complete-s+p target is not 8 atoms/48 local/384 global")
    if integers["checkpoint_format_version"] != 3:
        fail(f"{row_name}: row-owned checkpoint format V2 was not published")
    max_rows = integers["matrix_owned_rows_max"]
    if max_rows <= 0 or max_rows >= integers["global_target"]:
        fail(f"{row_name}: projected matrices are not distributed by owned rows")
    expected_local_bytes = max_rows * integers["global_target"] * 16
    if (
        integers["overlap_local_bytes_max"] != expected_local_bytes
        or integers["hamiltonian_local_bytes_max"] != expected_local_bytes
    ):
        fail(f"{row_name}: row-owned projected-matrix byte evidence is inconsistent")
    if (
        integers["complete_shell_channels"] != 32
        or integers["complete_shell_residual_rank"] != 32
        or integers["bond_center_orbit_closed"] != 1
    ):
        fail(f"{row_name}: complete shell or bond-center orbit closure is missing")
    if integers["occupied_required"] != 16:
        fail(f"{row_name}: occupied rank is not 16 per fragment")
    if integers["occupied_included"] != integers["occupied_required"]:
        fail(f"{row_name}: occupied inclusion is incomplete")
    if metrics["runtime_s"] <= 0.0 or metrics["peak_memory_mib"] <= 0.0:
        fail(f"{row_name}: runtime/memory evidence is invalid")
    if metrics["metric_min"] <= 0.0 or metrics["metric_max"] < metrics["metric_min"]:
        fail(f"{row_name}: overlap metric is not positive")
    if metrics["metric_condition"] > 1.0e10:
        fail(f"{row_name}: overlap metric is ill-conditioned")
    if metrics["tail_value_norm"] <= 0.0 or metrics["tail_gradient_norm"] <= 0.0:
        fail(f"{row_name}: retained buffer tails were not measured")
    for key, tolerance in TOLERANCES.items():
        if metrics[key] > tolerance:
            fail(f"{row_name}: {key}={metrics[key]:.6e} exceeds {tolerance:.6e}")
    for key in ("trace_charge", "integrated_charge"):
        if abs(metrics[key] - 256.0) > 1.0e-8:
            fail(f"{row_name}: {key} does not contain 256 electrons")
    manifest = row / "overlapping_wannier_gs.manifest"
    if not manifest.is_file() or digest(manifest) != raw["checkpoint_manifest_sha256"]:
        fail(f"{row_name}: checkpoint manifest hash mismatch")
    if manifest.read_bytes()[:32].rstrip(b" \x00") != b"SALMON_OW_GS_CHECKPOINT_V2":
        fail(f"{row_name}: checkpoint manifest is not V2")
    shards = sorted(row.glob("overlapping_wannier_gs.*.rank*"))
    if len(shards) != integers["mpi_ranks"]:
        fail(f"{row_name}: expected one V2 checkpoint shard per MPI rank")
    for shard in shards:
        if shard.read_bytes()[:32].rstrip(b" \x00") != b"SALMON_OW_GS_RANK_SHARD_V2":
            fail(f"{row_name}: checkpoint shard is not V2: {shard.name}")
    if not SHA256.fullmatch(raw["checkpoint_manifest_sha256"]):
        fail(f"{row_name}: malformed checkpoint manifest digest")
    validate_restart_log(row / "restart.log")
    return metrics


def relative_difference(a: float, b: float) -> float:
    return abs(a - b) / max(1.0, abs(a), abs(b))


def validate_reference(root: Path) -> None:
    reference = root / "normal-reference"
    log_path = reference / "run.log"
    variables_path = reference / "variables.log"
    if not log_path.is_file() or not variables_path.is_file():
        fail("normal reference raw log or variables.log is missing")
    log = log_path.read_text(errors="replace")
    variables = variables_path.read_text(errors="replace")
    validate_fixed_decomposition(variables, "normal reference")
    variables = variables.lower().replace(" ", "")
    trace = DC_SCF.findall(log)
    if not trace:
        fail("normal reference contains no DC SCF trace")
    final_difference = float(trace[-1][2].replace("D", "E").replace("d", "e"))
    if final_difference > 1.0e-9:
        fail(f"normal reference SCF did not reach 1e-9: {final_difference:.6e}")
    if "dc-lcfo" not in log.lower() or "eigenexa diag" not in log.lower():
        fail("normal reference lacks LCFO or EigenExa runtime markers")
    for assignment in ("yn_dg_dc_overlapping_wannier=n", "yn_dc_lcfo=y", "yn_eigenexa=y"):
        if assignment not in variables:
            fail(f"normal reference variables lack {assignment}")


def validate_row_storage_contract(repo: Path) -> None:
    checker = repo / "tests/dg/check_dg_overlapping_wannier_row_storage.py"
    result = subprocess.run(
        [sys.executable, str(checker)],
        cwd=repo,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0 or "row-storage source contract: PASS" not in result.stdout:
        fail(f"row-storage source contract failed:\n{result.stdout.rstrip()}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--skip-reference", action="store_true")
    parser.add_argument("--minimum-row", action="store_true")
    args = parser.parse_args()
    if args.minimum_row != args.skip_reference:
        parser.error("--minimum-row and --skip-reference must be used together")
    root = args.result_root.resolve()
    rows: dict[tuple[str, int, str], dict[str, float]] = {}
    try:
        if not root.is_dir():
            fail(f"result root does not exist: {root}")
        validate_row_storage_contract(args.repo.resolve())
        if args.minimum_row:
            validate_row(root, "2x2x2", "buffer6", "c192-sp48")
            print("Si64 overlapping-Wannier minimum row gate: PASS")
            return 0
        validate_reference(root)
        for coordinate in itertools.product(DECOMPOSITIONS, BOX_PROFILES, WINDOWS):
            rows[coordinate] = validate_row(root, *coordinate)
        for decomposition in DECOMPOSITIONS:
            for key in CONVERGENCE_KEYS:
                values = [
                    rows[(decomposition, PRODUCTION_BOX_PROFILE, window)][key]
                    for window in WINDOWS
                ]
                if relative_difference(values[0], values[1]) > 1.0e-4:
                    fail(f"{decomposition}: production candidate-window convergence failed for {key}")
    except (GateFailure, OSError, ValueError) as exc:
        print(f"Si64 overlapping-Wannier gate: FAIL: {exc}", file=sys.stderr)
        return 1
    print("Si64 overlapping-Wannier gate: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

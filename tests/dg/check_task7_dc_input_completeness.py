#!/usr/bin/env python3
"""Require complete &dc default/broadcast/variables.log coverage."""

from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "src/io/inputoutput.f90"
SOURCE = PATH.read_text(errors="replace").lower()


def preprocess(source: str) -> str:
    result = subprocess.run(
        ["cpp", "-P", "-traditional-cpp", f"-I{ROOT / 'gnu_makefiles'}", "-DUSE_MPI", "-DUSE_SCALAPACK", "-"],
        input=source,
        text=True,
        capture_output=True,
        cwd=ROOT,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout.lower()


def dc_members(source: str) -> list[str]:
    block = source[source.index("namelist/dc/") : source.index("namelist/unfolding/")]
    block = block.split("/dc/", 1)[1]
    return re.findall(r"\b[a-z][a-z0-9_]*\b", block.replace("&", " "))


def coverage_violations(source: str) -> dict[str, list[str]]:
    members = dc_members(source)
    assert len(members) == len(set(members)) == 102, (len(members), len(set(members)))
    cooked = preprocess(source)
    defaults = cooked[: cooked.index("!! == bcast for &calculation")]
    broadcasts = cooked[cooked.index("!! == bcast for &calculation") :]
    log_start = source.index("'dc', inml_dc")
    logs = source[log_start : source.index("'unfolding', inml_unfolding", log_start)]
    result: dict[str, list[str]] = {
        "missing_default": [],
        "duplicate_default": [],
        "missing_broadcast": [],
        "duplicate_broadcast": [],
        "missing_log": [],
        "duplicate_log": [],
    }
    for name in members:
        default_count = len(
            re.findall(rf"^\s*{re.escape(name)}(?:\s*\([^\n=]*\))?\s*=", defaults, re.M)
        )
        broadcast_count = len(
            re.findall(rf"call\s+comm_bcast\s*\(\s*{re.escape(name)}\b", broadcasts)
        )
        log_count = len(re.findall(rf"['\"]{re.escape(name)}['\"]", logs))
        if default_count == 0:
            result["missing_default"].append(name)
        elif default_count > 1:
            result["duplicate_default"].append(name)
        if broadcast_count == 0:
            result["missing_broadcast"].append(name)
        elif broadcast_count > 1:
            result["duplicate_broadcast"].append(name)
        if log_count == 0:
            result["missing_log"].append(name)
        elif log_count > 1:
            result["duplicate_log"].append(name)
    return result


coverage = coverage_violations(SOURCE)
assert not any(coverage.values()), coverage

# Independent semantic mutations prove all three dimensions enforce exactly
# one occurrence rather than mere token presence.
default_line = "    dg_dc_seed_mode = 'off'"
missing_default = coverage_violations(SOURCE.replace(default_line, "    removed_default = 'off'", 1))
assert missing_default["missing_default"] == ["dg_dc_seed_mode"], missing_default
duplicate_default = coverage_violations(SOURCE.replace(default_line, default_line + "\n" + default_line, 1))
assert duplicate_default["duplicate_default"] == ["dg_dc_seed_mode"], duplicate_default

broadcast_line = "    call comm_bcast(dg_dc_seed_mode, nproc_group_global)"
missing_broadcast = coverage_violations(SOURCE.replace(broadcast_line, "", 1))
assert missing_broadcast["missing_broadcast"] == ["dg_dc_seed_mode"], missing_broadcast
duplicate_broadcast = coverage_violations(
    SOURCE.replace(broadcast_line, broadcast_line + "\n" + broadcast_line, 1)
)
assert duplicate_broadcast["duplicate_broadcast"] == ["dg_dc_seed_mode"], duplicate_broadcast

missing_log = coverage_violations(SOURCE.replace("'dg_dc_seed_mode'", "'removed_dg_dc_seed_mode'", 1))
assert missing_log["missing_log"] == ["dg_dc_seed_mode"], missing_log
duplicate_log = coverage_violations(
    SOURCE.replace("'dg_dc_seed_mode'", "'dg_dc_seed_mode'//'dg_dc_seed_mode'", 1)
)
assert duplicate_log["duplicate_log"] == ["dg_dc_seed_mode"], duplicate_log

print("PASS all 102 &dc controls have one default, broadcast, and variables.log entry")

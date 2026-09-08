#!/usr/bin/env python3
"""Short Si8 production smoke for fragment-WF miss, hit, and incomplete recovery."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import struct
import subprocess
import time
from pathlib import Path


RANKS = 8
ELECTRON_TOLERANCE = 1.0e-8
UINT64_MASK = (1 << 64) - 1


def mix_occupied_hash(seed: int, value: int) -> int:
    seed &= UINT64_MASK
    rotated = ((seed << 9) | (seed >> 55)) & UINT64_MASK
    return rotated ^ (value & UINT64_MASK)


def validate_occupied_checkpoint(path: Path) -> int:
    """Read back and authenticate the native occupied-checkpoint stream."""
    payload = path.read_bytes()
    offset = 0

    def take(fmt: str) -> tuple:
        nonlocal offset
        size = struct.calcsize(fmt)
        if offset + size > len(payload):
            raise RuntimeError(f"truncated occupied checkpoint: {path}")
        values = struct.unpack_from(fmt, payload, offset)
        offset += size
        return values

    magic = payload[:16]
    offset = 16
    version, global_count, occupied_count = take("=iii")
    if magic != b"SALMON_DG_OCC02 " or version != 2:
        raise RuntimeError(f"unknown occupied checkpoint header: {path}")
    if global_count < 1 or occupied_count < 1:
        raise RuntimeError(f"invalid occupied checkpoint dimensions: {path}")
    fingerprints = take("=10Q")
    catalog, basis = fingerprints[:2]
    provenance = fingerprints[2:8]
    operator, state = fingerprints[8:]
    if any(value == 0 for value in fingerprints):
        raise RuntimeError(f"zero occupied checkpoint fingerprint: {path}")

    occupation_bits = take(f"={occupied_count}Q")
    eigenvalue_bits = take(f"={occupied_count}Q")
    receipt_bits = take("=5Q")
    occupations = struct.unpack(f"={occupied_count}d", struct.pack(f"={occupied_count}Q", *occupation_bits))
    eigenvalues = struct.unpack(f"={occupied_count}d", struct.pack(f"={occupied_count}Q", *eigenvalue_bits))
    receipts = struct.unpack("=5d", struct.pack("=5Q", *receipt_bits))
    if (not all(math.isfinite(value) and value >= 0.0 for value in occupations)
            or not all(math.isfinite(value) for value in eigenvalues)
            or not all(math.isfinite(value) and value >= 0.0 for value in receipts)):
        raise RuntimeError(f"nonfinite or negative occupied checkpoint payload: {path}")

    fingerprint = catalog
    for value in (basis, *provenance, operator, state, global_count, occupied_count):
        fingerprint = mix_occupied_hash(fingerprint, value)
    for occupation, eigenvalue in zip(occupation_bits, eigenvalue_bits):
        fingerprint = mix_occupied_hash(fingerprint, occupation)
        fingerprint = mix_occupied_hash(fingerprint, eigenvalue)
    for value in receipt_bits:
        fingerprint = mix_occupied_hash(fingerprint, value)
    coefficient_count = global_count * occupied_count
    coefficient_bits = take(f"={2 * coefficient_count}Q")
    coefficient_values = struct.unpack(
        f"={2 * coefficient_count}d", struct.pack(f"={2 * coefficient_count}Q", *coefficient_bits)
    )
    if not all(math.isfinite(value) for value in coefficient_values):
        raise RuntimeError(f"nonfinite occupied checkpoint coefficients: {path}")
    cursor = 0
    for row in range(1, global_count + 1):
        fingerprint = mix_occupied_hash(fingerprint, row)
        for _ in range(occupied_count):
            fingerprint = mix_occupied_hash(fingerprint, coefficient_bits[cursor])
            fingerprint = mix_occupied_hash(fingerprint, coefficient_bits[cursor + 1])
            cursor += 2
    if fingerprint == 0:
        fingerprint = 1
    stored_fingerprint, = take("=Q")
    if offset != len(payload) or stored_fingerprint != fingerprint:
        raise RuntimeError(f"corrupt occupied checkpoint fingerprint or extent: {path}")
    return stored_fingerprint if stored_fingerprint < (1 << 63) else stored_fingerprint - (1 << 64)


def replace_once(text: str, pattern: str, replacement: str) -> str:
    compiled = re.compile(pattern, re.IGNORECASE | re.MULTILINE)
    if len(compiled.findall(text)) != 1:
        raise RuntimeError(f"expected one input selector: {pattern}")
    return compiled.sub(replacement, text)


def render_input(template: str, dc_seed: Path, wf_checkpoint: Path) -> str:
    controls = f"""
  yn_dg_hybrid_scf = 'n'
  yn_dg_hybrid_continuation_scf = 'n'
  yn_dg_hybrid_divided_scf = 'y'
  dg_hybrid_divided_mixing = 'pulay'
  dg_dc_seed_mode = 'auto'
  dg_dc_seed_directory = '{dc_seed}'
  dg_fragment_wf_checkpoint_mode = 'auto'
  dg_fragment_wf_checkpoint_directory = '{wf_checkpoint}'
  dg_fragment_w90_initial_projection = 'scdm'
  dg_ow_localization_gradient_tolerance = 2.0d-2
  wannier_num_iter = 10000
  wannier_pw_cutoff = 0.5d0
"""
    rendered = replace_once(template, r"(?m)^\s*&dc\s*$", "&dc" + controls)
    rendered = replace_once(
        rendered,
        r"(?m)^(\s*yn_eigenexa\s*=\s*)['\"]\w+['\"]\s*$",
        r"\1'n'\n  yn_scalapack = 'y'",
    )
    return rendered


def one_match(pattern: str, text: str, label: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {label}, found {len(matches)}")
    return matches[0]


def parse_evidence(case: str, run_dir: Path, elapsed: float, return_code: int) -> dict:
    log_path = run_dir / "run.log"
    text = log_path.read_text(errors="replace")
    wf = one_match(
        r"\[DG-FRAGMENT-WF\]\s+mode=auto\s+checkpoint_hit=([TF])\s+"
        r"publication_id=(-?\d+)\s+reason=(.*)$",
        text,
        "fragment-WF receipt",
    )
    projected = one_match(
        r"\[DG-HYBRID-DIVIDED\]\s+projected_basis_fingerprint=(-?\d+)",
        text,
        "projected-basis fingerprint",
    )
    terminal = one_match(
        r"\[OW-GS\]\s+fixed-density/non-self-consistent divided WF\+PW LCFO solved once\s+"
        r"residual=\s*([0-9Ee+\-.]+)\s+orthogonality=\s*([0-9Ee+\-.]+)\s+"
        r"projector=\s*([0-9Ee+\-.]+)\s+electron_defect=\s*([0-9Ee+\-.]+)",
        text,
        "terminal LCFO receipt",
    )
    seed = one_match(
        r"\[DG-DC-SEED\]\s+mode=auto\s+publication_id=(-?\d+)\s+"
        r"scf_skipped=([TF])\s+mpi_size=(\d+)\s+mapping_fingerprint=(-?\d+)",
        text,
        "DC seed receipt",
    )
    schwarz = list(re.finditer(
        r"\[DG-HYBRID-SCHWARZ\].*?electron_defect=\s*([0-9Ee+\-.]+).*?"
        r"temperature=\s*([0-9Ee+\-.]+)",
        text,
        re.IGNORECASE,
    ))
    terminal_position = terminal.start()
    occupied = run_dir / "overlapping_wannier_occupied.chk"
    occupied_checkpoint_fingerprint = validate_occupied_checkpoint(occupied)
    evidence = {
        "case": case,
        "return_code": return_code,
        "elapsed_seconds": elapsed,
        "dc_publication_id": int(seed.group(1)),
        "dc_scf_skipped": seed.group(2) == "T",
        "dc_mpi_size": int(seed.group(3)),
        "mapping_fingerprint": int(seed.group(4)),
        "checkpoint_hit": wf.group(1) == "T",
        "wf_publication_id": int(wf.group(2)),
        "wf_reason": wf.group(3).strip(),
        "projected_basis_fingerprint": int(projected.group(1)),
        "terminal_residual": float(terminal.group(1)),
        "terminal_orthogonality": float(terminal.group(2)),
        "terminal_projector": float(terminal.group(3)),
        "electron_defect": float(terminal.group(4)),
        "schwarz_temperature_kelvin": float(schwarz[-1].group(2)) if schwarz else None,
        "terminal_lcfo_count": text.count(
            "[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once"
        ),
        "post_lcfo_density_updates": text[terminal_position:].count(
            "[DG-HYBRID-DIVIDED-SCF] iteration="
        ),
        "wannier_wout_count": len(list(run_dir.rglob("*.wout"))),
        "occupied_checkpoint": True,
        "occupied_checkpoint_fingerprint": occupied_checkpoint_fingerprint,
    }
    if return_code != 0:
        raise RuntimeError(f"{case} exited with {return_code}: {log_path}")
    if evidence["dc_mpi_size"] != RANKS:
        raise RuntimeError(f"{case} changed the rank/fragment contract")
    if evidence["terminal_lcfo_count"] != 1 or evidence["post_lcfo_density_updates"] != 0:
        raise RuntimeError(f"{case} did not perform exactly one terminal LCFO without density updates")
    if evidence["electron_defect"] > ELECTRON_TOLERANCE:
        raise RuntimeError(f"{case} did not reproduce 32 electrons at 300 K")
    if evidence["schwarz_temperature_kelvin"] != 300.0:
        raise RuntimeError(f"{case} changed the 300 K Schwarz temperature contract")
    if not evidence["occupied_checkpoint"]:
        raise RuntimeError(f"{case} did not publish a valid occupied checkpoint")
    return evidence


def run_case(binary: Path, run_dir: Path, input_text: str, fixture: Path, timeout: int) -> dict:
    run_dir.mkdir(parents=True)
    (run_dir / "inputfile").write_text(input_text)
    shutil.copy2(fixture / "atom.dat", run_dir / "atom.dat")
    root = Path(__file__).resolve().parents[2]
    shutil.copy2(root / "samples/exercise_04_bulkSi_gs/Si_rps.dat", run_dir / "Si_rps.dat")
    started = time.monotonic()
    with (run_dir / "inputfile").open("rb") as source, (run_dir / "run.log").open("wb") as log:
        completed = subprocess.run(
            [shutil.which("mpirun") or "mpirun", "-np", str(RANKS), str(binary)],
            cwd=run_dir,
            stdin=source,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={**os.environ, "OMP_NUM_THREADS": "1", "OMPI_MCA_rmaps_base_oversubscribe": "1"},
            timeout=timeout,
        )
    return parse_evidence(case=run_dir.name, run_dir=run_dir,
                          elapsed=time.monotonic() - started, return_code=completed.returncode)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path)
    parser.add_argument(
        "--seed-directory",
        type=Path,
        help="copy an existing authoritative DC seed into the fresh smoke directory",
    )
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    result = (args.result_dir or Path("/tmp") /
              f"dg-fragment-wf-smoke-{time.strftime('%Y%m%d-%H%M%S')}").resolve()
    if result.exists():
        raise RuntimeError(f"result directory must be fresh: {result}")
    result.mkdir(parents=True)
    fixture = root / "tests/dg/data/si8_overlapping_wannier"
    template = (fixture / "inputfile.in").read_text()
    dc_seed = result / "dc-seed"
    wf_checkpoint = result / "fragment-wf"
    seed_preloaded = args.seed_directory is not None
    if seed_preloaded:
        seed_source = args.seed_directory.resolve(strict=True)
        if not seed_source.is_dir():
            raise RuntimeError(f"DC seed source is not a directory: {seed_source}")
        shutil.copytree(seed_source, dc_seed)

    miss = run_case(binary, result / "01-miss", render_input(template, dc_seed, wf_checkpoint),
                    fixture, args.timeout)
    hit = run_case(binary, result / "02-hit", render_input(template, dc_seed, wf_checkpoint),
                   fixture, args.timeout)

    incomplete_checkpoint = result / "fragment-wf-incomplete"
    shutil.copytree(wf_checkpoint, incomplete_checkpoint)
    incomplete_manifest = incomplete_checkpoint / "dg_fragment_wf.manifest"
    incomplete_manifest.unlink()
    incomplete = run_case(
        binary,
        result / "03-incomplete-recovery",
        render_input(template, dc_seed, incomplete_checkpoint),
        fixture,
        args.timeout,
    )

    if miss["checkpoint_hit"] or not hit["checkpoint_hit"] or incomplete["checkpoint_hit"]:
        raise RuntimeError("fragment-WF miss/hit/incomplete auto decisions are incorrect")
    if miss["wf_publication_id"] != hit["wf_publication_id"]:
        raise RuntimeError("cache hit did not retain the published fragment-WF generation")
    if incomplete["wf_publication_id"] == miss["wf_publication_id"]:
        raise RuntimeError("incomplete publication was consumed instead of regenerated")
    fingerprints = {item["projected_basis_fingerprint"] for item in (miss, hit, incomplete)}
    if len(fingerprints) != 1:
        raise RuntimeError("projected basis changed across miss/hit/incomplete recovery")
    if miss["wannier_wout_count"] != RANKS or hit["wannier_wout_count"] != 0:
        raise RuntimeError("Wannier90 artifacts do not prove miss generation and hit bypass")
    if incomplete["wannier_wout_count"] != RANKS:
        raise RuntimeError("incomplete checkpoint recovery did not regenerate every fragment")
    expected_dc_skip = (True, True, True) if seed_preloaded else (False, True, True)
    if tuple(item["dc_scf_skipped"] for item in (miss, hit, incomplete)) != expected_dc_skip:
        raise RuntimeError("ordinary DC execution/reuse did not match the seed setup")
    if len({item["dc_publication_id"] for item in (miss, hit, incomplete)}) != 1:
        raise RuntimeError("ordinary-DC publication changed across the restart smoke")
    if len({item["mapping_fingerprint"] for item in (miss, hit, incomplete)}) != 1:
        raise RuntimeError("rank-fragment mapping changed across the restart smoke")

    evidence = {"mpi_ranks": RANKS, "omp_threads": 1, "seed_preloaded": seed_preloaded,
                "manifest_removed_for_incomplete_case": str(incomplete_manifest),
                "runs": [miss, hit, incomplete]}
    (result / "fragment_wf_smoke_evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

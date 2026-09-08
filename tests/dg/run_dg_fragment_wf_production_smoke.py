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
TERMINAL_SOLVER_TOLERANCE = 1.0e-10
UINT64_MASK = (1 << 64) - 1
FLOAT_TOKEN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
FLOAT_TOKEN_RE = re.compile(rf"{FLOAT_TOKEN}\Z", re.ASCII)


def parse_finite_float(
    token: str, label: str, *, minimum: float | None = None, maximum: float | None = None
) -> float:
    """Parse one complete decimal token and reject overflow/nonfinite values."""
    if FLOAT_TOKEN_RE.fullmatch(token) is None:
        raise RuntimeError(f"{label} is not a valid finite decimal scalar: {token!r}")
    value = float(token.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise RuntimeError(f"{label} is not finite: {token!r}")
    if minimum is not None and value < minimum:
        raise RuntimeError(f"{label} is below {minimum}: {value}")
    if maximum is not None and value > maximum:
        raise RuntimeError(f"{label} exceeds {maximum}: {value}")
    return value


def require_mpi_completion(text: str, expected_ranks: int, label: str) -> None:
    """Require the normal SALMON completion receipt from every MPI rank."""
    completion_count = len(re.findall(r"^\s*end SALMON\s*$", text, re.MULTILINE))
    if completion_count != expected_ranks:
        raise RuntimeError(
            f"{label} has incomplete process completion evidence: "
            f"expected {expected_ranks} end SALMON receipts, found {completion_count}"
        )


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


def continuation_receipts(text: str, case: str) -> list[dict[str, int | float | str]]:
    """Parse and validate the actual six DG continuation identity receipts."""
    pattern = re.compile(
        rf"\[DG-HYBRID-CONTINUATION\]\s+lambda=\s*({FLOAT_TOKEN})\s+"
        rf"diagnostic_state_lambda=\s*({FLOAT_TOKEN})\s+accepted_cg_steps=(\d+)\s+"
        rf"residual=\s*({FLOAT_TOKEN})\s+orthogonality_defect=\s*({FLOAT_TOKEN})\s+"
        rf"electron_defect=\s*({FLOAT_TOKEN})\s+rayleigh_energy_trace=\s*({FLOAT_TOKEN})\s+"
        rf"scaled_interface_action_norm=\s*({FLOAT_TOKEN})\s+measurement_status=(\w+)\s+"
        r"status=(\w+)\s+continuation_fingerprint=(-?\d+)",
        re.IGNORECASE,
    )
    receipts = [{
        "lambda": parse_finite_float(
            row.group(1), f"{case} continuation lambda", minimum=0.0, maximum=1.0),
        "diagnostic_state_lambda": parse_finite_float(
            row.group(2), f"{case} continuation diagnostic lambda", minimum=0.0,
            maximum=1.0),
        "accepted_cg_steps": int(row.group(3)),
        "residual": parse_finite_float(
            row.group(4), f"{case} continuation residual", minimum=0.0),
        "orthogonality_defect": parse_finite_float(
            row.group(5), f"{case} continuation orthogonality defect", minimum=0.0,
            maximum=TERMINAL_SOLVER_TOLERANCE),
        "electron_defect": parse_finite_float(
            row.group(6), f"{case} continuation electron defect", minimum=0.0,
            maximum=ELECTRON_TOLERANCE),
        "rayleigh_energy_trace": parse_finite_float(
            row.group(7), f"{case} continuation Rayleigh energy trace"),
        "scaled_interface_action_norm": parse_finite_float(
            row.group(8), f"{case} continuation scaled interface action norm", minimum=0.0),
        "measurement_status": row.group(9),
        "status": row.group(10),
        "continuation_fingerprint": int(row.group(11)),
    } for row in pattern.finditer(text)]
    valid_schedule = (
        len(receipts) == 6
        and all(math.isclose(item["lambda"], 0.2 * index, rel_tol=0.0, abs_tol=1.0e-14)
                and math.isclose(item["diagnostic_state_lambda"], item["lambda"],
                                 rel_tol=0.0, abs_tol=1.0e-14)
                for index, item in enumerate(receipts))
        and all(item["measurement_status"].lower() == "valid"
                and item["status"].lower() == "accepted"
                and 1 <= item["accepted_cg_steps"] <= 3
                and item["continuation_fingerprint"] != 0 for item in receipts)
    )
    if not valid_schedule:
        raise RuntimeError(f"{case} did not complete the exact six-record DG continuation contract")
    return receipts


def wannier90_fragment_ids(run_dir: Path) -> list[int]:
    """Return a unique fragment inventory for every Wannier90 output."""
    fragment_ids = []
    for path in run_dir.rglob("*.wout"):
        match = re.search(r"(?:^|/)(?:f|fragment-)(\d{6})(?:/|$)", path.as_posix())
        if not match:
            raise RuntimeError(f"cannot identify fragment from Wannier90 output: {path}")
        fragment_ids.append(int(match.group(1)))
    if len(fragment_ids) != len(set(fragment_ids)):
        raise RuntimeError("duplicate Wannier90 fragment output")
    return sorted(fragment_ids)


def parse_evidence(case: str, run_dir: Path, elapsed: float, return_code: int | None) -> dict:
    log_path = run_dir / "run.log"
    text = log_path.read_text(errors="replace")
    if return_code is not None and return_code != 0:
        raise RuntimeError(f"{case} exited with {return_code}: {log_path}")
    require_mpi_completion(text, RANKS, case)
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
        rf"residual=\s*({FLOAT_TOKEN})\s+orthogonality=\s*({FLOAT_TOKEN})\s+"
        rf"projector=\s*({FLOAT_TOKEN})\s+electron_defect=\s*({FLOAT_TOKEN})",
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
        rf"\[DG-HYBRID-SCHWARZ\]\s+epoch=(\d+)\s+neighbor_exchanges=(\d+)\s+"
        rf"accepted_cg_steps=(\d+)\s+common_extensions=(\d+)\s+"
        rf"residual=\s*({FLOAT_TOKEN})\s+electron_defect=\s*({FLOAT_TOKEN})\s+"
        rf"temperature=\s*({FLOAT_TOKEN})",
        text,
        re.IGNORECASE,
    ))
    continuation = continuation_receipts(text, case)
    terminal_position = terminal.start()
    occupied = run_dir / "overlapping_wannier_occupied.chk"
    occupied_checkpoint_fingerprint = validate_occupied_checkpoint(occupied)
    terminal_residual = parse_finite_float(
        terminal.group(1), f"{case} terminal residual", minimum=0.0,
        maximum=TERMINAL_SOLVER_TOLERANCE)
    terminal_orthogonality = parse_finite_float(
        terminal.group(2), f"{case} terminal orthogonality", minimum=0.0,
        maximum=TERMINAL_SOLVER_TOLERANCE)
    terminal_projector = parse_finite_float(
        terminal.group(3), f"{case} terminal projector", minimum=0.0,
        maximum=TERMINAL_SOLVER_TOLERANCE)
    terminal_electron_defect = parse_finite_float(
        terminal.group(4), f"{case} terminal electron defect", minimum=0.0,
        maximum=ELECTRON_TOLERANCE)
    schwarz_receipts = [{
        "epoch": int(row.group(1)),
        "neighbor_exchanges": int(row.group(2)),
        "accepted_cg_steps": int(row.group(3)),
        "common_extensions": int(row.group(4)),
        "residual": parse_finite_float(
            row.group(5), f"{case} Schwarz residual", minimum=0.0),
        "electron_defect": parse_finite_float(
            row.group(6), f"{case} Schwarz electron defect", minimum=0.0,
            maximum=ELECTRON_TOLERANCE),
        "temperature": parse_finite_float(
            row.group(7), f"{case} Schwarz temperature", minimum=0.0),
    } for row in schwarz]
    if ([item["epoch"] for item in schwarz_receipts] != list(range(1, 7))
            or any(item["neighbor_exchanges"] != RANKS - 1
                   or not 1 <= item["accepted_cg_steps"] <= 3
                   or item["common_extensions"] < 0 for item in schwarz_receipts)):
        raise RuntimeError(f"{case} did not complete the exact six-stage DG continuation schedule")
    if not math.isfinite(elapsed) or elapsed < 0.0:
        raise RuntimeError(f"{case} elapsed time is not finite and nonnegative")
    publication_values = tuple(int(value) for value in (
        seed.group(1), seed.group(4), wf.group(2), projected.group(1)))
    if any(value == 0 for value in publication_values):
        raise RuntimeError(f"{case} contains a zero publication or fingerprint identity")
    w90_fragment_ids = wannier90_fragment_ids(run_dir)
    evidence = {
        "case": case,
        "return_code": return_code,
        "completion_receipts": RANKS,
        "elapsed_seconds": elapsed,
        "dc_publication_id": int(seed.group(1)),
        "dc_scf_skipped": seed.group(2) == "T",
        "dc_mpi_size": int(seed.group(3)),
        "mapping_fingerprint": int(seed.group(4)),
        "checkpoint_hit": wf.group(1) == "T",
        "wf_publication_id": int(wf.group(2)),
        "wf_reason": wf.group(3).strip(),
        "projected_basis_fingerprint": int(projected.group(1)),
        "terminal_residual": terminal_residual,
        "terminal_orthogonality": terminal_orthogonality,
        "terminal_projector": terminal_projector,
        "electron_defect": terminal_electron_defect,
        "schwarz_receipts": schwarz_receipts,
        "continuation_receipts": continuation,
        "continuation_fingerprints": [
            item["continuation_fingerprint"] for item in continuation],
        "schwarz_temperature_kelvin": schwarz_receipts[-1]["temperature"] if schwarz_receipts else None,
        "terminal_lcfo_count": text.count(
            "[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once"
        ),
        "post_lcfo_density_updates": text[terminal_position:].count(
            "[DG-HYBRID-DIVIDED-SCF] iteration="
        ),
        "wannier_wout_count": len(w90_fragment_ids),
        "wannier_fragment_ids": w90_fragment_ids,
        "occupied_checkpoint": True,
        "occupied_checkpoint_fingerprint": occupied_checkpoint_fingerprint,
    }
    if evidence["dc_mpi_size"] != RANKS:
        raise RuntimeError(f"{case} changed the rank/fragment contract")
    if evidence["terminal_lcfo_count"] != 1 or evidence["post_lcfo_density_updates"] != 0:
        raise RuntimeError(f"{case} did not perform exactly one terminal LCFO without density updates")
    if evidence["schwarz_temperature_kelvin"] != 300.0:
        raise RuntimeError(f"{case} changed the 300 K Schwarz temperature contract")
    if not evidence["occupied_checkpoint"]:
        raise RuntimeError(f"{case} did not publish a valid occupied checkpoint")
    return evidence


def compare_smoke_runs(miss: dict, hit: dict, incomplete: dict, seed_preloaded: bool) -> None:
    """Enforce miss/hit/recovery identity and restart decisions."""
    if miss["checkpoint_hit"] or not hit["checkpoint_hit"] or incomplete["checkpoint_hit"]:
        raise RuntimeError("fragment-WF miss/hit/incomplete auto decisions are incorrect")
    if miss["wf_publication_id"] != hit["wf_publication_id"]:
        raise RuntimeError("cache hit did not retain the published fragment-WF generation")
    if incomplete["wf_publication_id"] == miss["wf_publication_id"]:
        raise RuntimeError("incomplete publication was consumed instead of regenerated")
    fingerprints = {item["projected_basis_fingerprint"] for item in (miss, hit, incomplete)}
    if len(fingerprints) != 1:
        raise RuntimeError("projected basis changed across miss/hit/incomplete recovery")
    if not (miss["continuation_receipts"] == hit["continuation_receipts"]
            == incomplete["continuation_receipts"]):
        raise RuntimeError("DG continuation identity receipts changed across miss/hit/incomplete recovery")
    if not (miss["schwarz_receipts"] == hit["schwarz_receipts"]
            == incomplete["schwarz_receipts"]):
        raise RuntimeError("DG Schwarz summaries changed across miss/hit/incomplete recovery")
    expected_fragment_ids = list(range(1, RANKS + 1))
    if miss["wannier_fragment_ids"] != expected_fragment_ids or hit["wannier_fragment_ids"]:
        raise RuntimeError("Wannier90 artifacts do not prove miss generation and hit bypass")
    if incomplete["wannier_fragment_ids"] != expected_fragment_ids:
        raise RuntimeError("incomplete checkpoint recovery did not regenerate every fragment")
    expected_dc_skip = (True, True, True) if seed_preloaded else (False, True, True)
    if tuple(item["dc_scf_skipped"] for item in (miss, hit, incomplete)) != expected_dc_skip:
        raise RuntimeError("ordinary DC execution/reuse did not match the seed setup")
    if len({item["dc_publication_id"] for item in (miss, hit, incomplete)}) != 1:
        raise RuntimeError("ordinary-DC publication changed across the restart smoke")
    if len({item["mapping_fingerprint"] for item in (miss, hit, incomplete)}) != 1:
        raise RuntimeError("rank-fragment mapping changed across the restart smoke")


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
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--result-dir", type=Path)
    parser.add_argument(
        "--seed-directory",
        type=Path,
        help="copy an existing authoritative DC seed into the fresh smoke directory",
    )
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--analyze-existing", action="store_true",
                        help="analyze an already completed smoke result without launching MPI")
    parser.add_argument("--evidence-output", type=Path,
                        help="explicit JSON output path; analyze-existing requires it outside result-dir")
    args = parser.parse_args()
    result = (args.result_dir or Path("/tmp") /
              f"dg-fragment-wf-smoke-{time.strftime('%Y%m%d-%H%M%S')}").resolve()
    if result.exists() and not args.analyze_existing:
        raise RuntimeError(f"result directory must be fresh: {result}")
    if not result.exists() and args.analyze_existing:
        raise RuntimeError(f"existing result directory is missing: {result}")
    if args.binary is None and not args.analyze_existing:
        raise RuntimeError("--binary is required unless --analyze-existing is selected")
    binary = args.binary.resolve(strict=True) if args.binary is not None else None
    result.mkdir(parents=True, exist_ok=args.analyze_existing)
    fixture = root / "tests/dg/data/si8_overlapping_wannier"
    template = (fixture / "inputfile.in").read_text()
    dc_seed = result / "dc-seed"
    wf_checkpoint = result / "fragment-wf"
    seed_preloaded = args.seed_directory is not None
    if seed_preloaded and not args.analyze_existing:
        seed_source = args.seed_directory.resolve(strict=True)
        if not seed_source.is_dir():
            raise RuntimeError(f"DC seed source is not a directory: {seed_source}")
        shutil.copytree(seed_source, dc_seed)

    incomplete_checkpoint = result / "fragment-wf-incomplete"
    incomplete_manifest = incomplete_checkpoint / "dg_fragment_wf.manifest"
    if args.analyze_existing:
        elapsed = lambda directory: max(
            0.0,
            (directory / "overlapping_wannier_occupied.chk").stat().st_mtime
            - (directory / "inputfile").stat().st_mtime,
        )
        miss_dir = result / "01-miss"
        hit_dir = result / "02-hit"
        incomplete_dir = result / "03-incomplete-recovery"
        miss = parse_evidence(miss_dir.name, miss_dir, elapsed(miss_dir), None)
        hit = parse_evidence(hit_dir.name, hit_dir, elapsed(hit_dir), None)
        incomplete = parse_evidence(
            incomplete_dir.name, incomplete_dir, elapsed(incomplete_dir), None)
        seed_preloaded = miss["dc_scf_skipped"]
    else:
        assert binary is not None
        miss = run_case(binary, result / "01-miss", render_input(template, dc_seed, wf_checkpoint),
                        fixture, args.timeout)
        hit = run_case(binary, result / "02-hit", render_input(template, dc_seed, wf_checkpoint),
                       fixture, args.timeout)
        shutil.copytree(wf_checkpoint, incomplete_checkpoint)
        incomplete_manifest.unlink()
        incomplete = run_case(
            binary,
            result / "03-incomplete-recovery",
            render_input(template, dc_seed, incomplete_checkpoint),
            fixture,
            args.timeout,
        )

    compare_smoke_runs(miss, hit, incomplete, seed_preloaded)

    evidence = {"mpi_ranks": RANKS, "omp_threads": 1, "seed_preloaded": seed_preloaded,
                "manifest_removed_for_incomplete_case": str(incomplete_manifest),
                "runs": [miss, hit, incomplete]}
    if args.analyze_existing:
        payload = {"status": "PASS", "result_dir": str(result), "evidence": evidence}
        if args.evidence_output is not None:
            output = args.evidence_output.resolve()
            try:
                output.relative_to(result)
            except ValueError:
                pass
            else:
                raise RuntimeError("analyze-existing evidence output must be outside result-dir")
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        output = (args.evidence_output or result / "fragment_wf_smoke_evidence.json").resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(evidence, indent=2) + "\n")
        print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

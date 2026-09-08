#!/usr/bin/env python3
"""Certify one clean Si64 SCDM fragment-WF build and its exact reuse."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import struct
import subprocess
import time
from pathlib import Path

import numpy as np

from run_dg_fragment_wf_production_smoke import (
    FLOAT_TOKEN,
    parse_finite_float,
    require_mpi_completion,
    validate_occupied_checkpoint,
)


EXPECTED_PUBLICATION_ID = 7047888166118007469
EXPECTED_MPI_SIZE = 8
EXPECTED_MAPPING_FINGERPRINT = 254086644876463474
DEFAULT_SEED_DIRECTORY = Path("/tmp/si64-task8-bounded-smoke-20260905/dc-seed")
RANDOM_W90_ITERATIONS = [1342, 1192, 2640, 3058, 1927, 1837, 779, 783]
ELECTRON_TOLERANCE = 1.0e-8
STATE_NUMERIC_TOLERANCE = 1.0e-9
TERMINAL_SOLVER_TOLERANCE = 1.0e-10
SEED_FILE_SHA256 = {
    "dg_dc_seed.manifest": "abf70a65c25a249c8575ed6d3550a7c54ba8f4946e93917046aed9bc96fc11aa",
    "rank00000000": "74bf3fa914879d99f024f64fd53c45a566b0501fefda01d652b0bf59274cba0e",
    "rank00000001": "947d294a881d1fbf5fd35c28173b011c7a8d77b22dfbdd1f08f444d08e4bc2f5",
    "rank00000002": "5490f01b353d283f16b2ba8e426e9b536cb7c959c3a1ef15091b5f06e15a6905",
    "rank00000003": "7515a2bf58325e69621726d70a6cb45620b87888efaa6e17611f6502169dcce4",
    "rank00000004": "0447a0836816cc845eea9ecf4e7a85980f7964cf7aa6f15a05d29012378c8561",
    "rank00000005": "dc12673cda8197fa7f5432e8f9a659301c670e8ecd55eaacd55ac3b983d12132",
    "rank00000006": "243ebe8748fcd906f73647c3d9513a7320e5021af393577d1576b476415f4d19",
    "rank00000007": "4f983b01725e0e2e2e9f69587688f39308699335c130f61f8668f5f9add580e7",
}


def validate_seed_directory(seed_directory: Path) -> dict[str, int]:
    """Reject before MPI startup unless the exact authoritative DC seed is present."""
    manifest = seed_directory / "dg_dc_seed.manifest"
    if not manifest.is_file():
        raise RuntimeError(f"missing reusable seed manifest: {manifest}")
    header = manifest.read_bytes()[:48]
    if len(header) != 48:
        raise RuntimeError("truncated reusable seed manifest")
    magic, version, mpi_size, publication_id = struct.unpack("=32siiq", header)
    if magic.rstrip() != b"SALMON_DG_DC_SEED_MANIFEST_V1":
        raise RuntimeError("unexpected reusable seed manifest magic")
    if (version, mpi_size, publication_id) != (1, EXPECTED_MPI_SIZE, EXPECTED_PUBLICATION_ID):
        raise RuntimeError("reusable seed manifest identity changed")
    publication_hex = f"{publication_id:016X}"
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
        return digest.hexdigest()

    if sha256(manifest) != SEED_FILE_SHA256["dg_dc_seed.manifest"]:
        raise RuntimeError("reusable seed manifest content changed")
    mapping = None
    for rank in range(EXPECTED_MPI_SIZE):
        shard = seed_directory / f"dg_dc_seed.{publication_hex}.rank{rank:08d}.shard"
        if shard.is_file():
            with shard.open("rb") as stream:
                shard_header = stream.read(80)
        else:
            shard_header = b""
        if len(shard_header) != 80:
            raise RuntimeError(f"missing or truncated reusable seed shard for rank {rank}")
        values = struct.unpack("=32siiiiqqqq", shard_header)
        shard_magic, shard_version, shard_mpi, shard_rank, fragment_id = values[:5]
        shard_publication, _, ownership, _ = values[5:]
        if (
            shard_magic.rstrip() != b"SALMON_DG_DC_SEED_SHARD_V1"
            or (shard_version, shard_mpi, shard_rank, fragment_id)
            != (1, EXPECTED_MPI_SIZE, rank, rank + 1)
            or shard_publication != publication_id
        ):
            raise RuntimeError(f"reusable seed rank-fragment contract changed at rank {rank}")
        if ownership != EXPECTED_MAPPING_FINGERPRINT:
            raise RuntimeError(f"reusable seed mapping fingerprint changed at rank {rank}")
        expected_sha256 = SEED_FILE_SHA256[f"rank{rank:08d}"]
        if sha256(shard) != expected_sha256:
            raise RuntimeError(f"reusable seed shard content changed at rank {rank}")
        mapping = ownership
    if mapping != EXPECTED_MAPPING_FINGERPRINT:
        raise RuntimeError("reusable seed mapping fingerprint changed")
    return {"publication_id": publication_id, "mpi_size": mpi_size, "mapping_fingerprint": mapping}


def replace_one(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"^(\s*{re.escape(name)}\s*=\s*)[^\n]+$", re.I | re.M)
    if len(pattern.findall(text)) != 1:
        raise RuntimeError(f"Si64 fixture must contain exactly one {name}")
    return pattern.sub(lambda match: match.group(1) + value, text)


def render_input(template: str, seed_directory: Path, wf_checkpoint: Path) -> str:
    rendered = template
    for name, value in {
        "yn_dg_hybrid_scf": "'n'",
        "yn_dg_hybrid_divided_scf": "'y'",
        "dg_hybrid_divided_mixing": "'pulay'",
        "dg_dc_seed_mode": "'read'",
        "dg_dc_seed_directory": f"'{seed_directory}'",
        "dg_fragment_wf_checkpoint_mode": "'auto'",
        "dg_fragment_wf_checkpoint_directory": f"'{wf_checkpoint}'",
        "dg_fragment_w90_initial_projection": "'scdm'",
    }.items():
        rendered = replace_one(rendered, name, value)
    return rendered


def one_match(pattern: str, text: str, label: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, text, re.I | re.M))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {label}, found {len(matches)}")
    return matches[0]


def wannier90_receipts(run_dir: Path) -> list[dict[str, int | float]]:
    receipts = []
    for path in run_dir.rglob("*.wout"):
        match = re.search(r"(?:^|/)(?:f|fragment-)(\d{6})(?:/|$)", path.as_posix())
        if not match:
            raise RuntimeError(f"cannot identify fragment from Wannier90 output: {path}")
        text = path.read_text(errors="replace")
        iterations = [int(value) for value in re.findall(r"^\s*(\d+).*?<--\s*CONV\s*$", text, re.M)]
        wall = re.findall(rf"Total Execution Time\s+({FLOAT_TOKEN})\s*\(sec\)", text)
        if not iterations or len(wall) != 1 or "Wannierisation convergence criteria satisfied" not in text:
            raise RuntimeError(f"incomplete or unconverged Wannier90 output: {path}")
        wall_seconds = parse_finite_float(
            wall[0], f"Wannier90 wall time for fragment {int(match.group(1))}", minimum=0.0)
        if max(iterations) < 1:
            raise RuntimeError(f"invalid Wannier90 iteration receipt: {path}")
        receipts.append({
            "fragment_id": int(match.group(1)),
            "iterations": max(iterations),
            "wall_seconds": wall_seconds,
        })
    fragment_ids = [int(item["fragment_id"]) for item in receipts]
    if len(fragment_ids) != len(set(fragment_ids)):
        raise RuntimeError("duplicate Wannier90 fragment receipt")
    return sorted(receipts, key=lambda item: int(item["fragment_id"]))


def fragment_wf_contract_receipts(
    checkpoint_directory: Path, publication_id: int
) -> list[dict[str, int | float | str]]:
    """Read the authenticated gauge/version identities stored in rank shards."""
    header_format = "=32s5iq10i16s13q"
    header_size = struct.calcsize(header_format)
    receipts = []
    publication_hex = f"{publication_id:016X}"
    pattern = f"dg_fragment_wf.publication-{publication_hex}.rank-*.bin"
    for path in sorted(checkpoint_directory.glob(pattern)):
        with path.open("rb") as stream:
            header = stream.read(header_size + 8)
        if len(header) != header_size + 8:
            raise RuntimeError(f"truncated fragment-WF shard header: {path}")
        values = struct.unpack(header_format, header[:header_size])
        seed_reconstruction_defect, = struct.unpack("=d", header[header_size:])
        if values[0].rstrip() != b"SALMON_DG_FRAGMENT_WF_SHARD_V1":
            raise RuntimeError(f"unexpected fragment-WF shard magic: {path}")
        contract = values[7:17]
        fingerprints = values[18:]
        if (values[6] != publication_id or fingerprints[0] != EXPECTED_MAPPING_FINGERPRINT
                or fingerprints[1] != EXPECTED_PUBLICATION_ID):
            raise RuntimeError(f"fragment-WF shard provenance changed: {path}")
        receipts.append({
            "rank": contract[2], "fragment_id": contract[3],
            "gauge_algorithm_version": contract[5],
            "candidate_rank": contract[6], "retained_rank": contract[7],
            "gauge_mode": values[17].rstrip().decode("ascii"),
            "selection_fingerprint": fingerprints[10],
            "gauge_fingerprint": fingerprints[11],
            "seed_reconstruction_defect": seed_reconstruction_defect,
        })
    if len(receipts) != EXPECTED_MPI_SIZE:
        raise RuntimeError("fragment-WF publication does not contain eight gauge contracts")
    rank_fragment_pairs = {
        (int(item["rank"]), int(item["fragment_id"])) for item in receipts
    }
    expected_pairs = {(rank, rank + 1) for rank in range(EXPECTED_MPI_SIZE)}
    if rank_fragment_pairs != expected_pairs:
        raise RuntimeError("fragment-WF rank-fragment gauge contract is not exact and unique")
    if any(item["rank"] + 1 != item["fragment_id"] or item["gauge_mode"] != "scdm"
           or item["gauge_algorithm_version"] != 1
           or item["candidate_rank"] != 400 or item["retained_rank"] != 400
           or item["selection_fingerprint"] == 0 or item["gauge_fingerprint"] == 0
           or not math.isfinite(item["seed_reconstruction_defect"])
           or item["seed_reconstruction_defect"] < 0.0
           or item["seed_reconstruction_defect"] > 1.0e-9 for item in receipts):
        raise RuntimeError("fragment-WF SCDM gauge contract changed")
    return receipts


def occupied_numeric_comparison(clean_path: Path, reuse_path: Path) -> dict[str, int | float | list[int]]:
    """Compare state payloads numerically; the integrity hash is deliberately bitwise."""
    def read(path: Path) -> tuple[tuple[int, int, int], tuple[int, ...], tuple[float, ...],
                                  tuple[float, ...], tuple[float, ...], memoryview]:
        payload = path.read_bytes()
        if len(payload) < 28 or payload[:16] != b"SALMON_DG_OCC02 ":
            raise RuntimeError(f"invalid occupied checkpoint header: {path}")
        version, global_count, occupied_count = struct.unpack_from("=iii", payload, 16)
        if version != 2 or global_count < 1 or occupied_count < 1:
            raise RuntimeError(f"invalid occupied checkpoint dimensions: {path}")
        offset = 28
        minimum_size = offset + 80 + 16 * occupied_count + 40 + 16 * global_count * occupied_count
        if len(payload) < minimum_size:
            raise RuntimeError(f"truncated occupied checkpoint numeric payload: {path}")
        fingerprints = struct.unpack_from("=10Q", payload, offset)
        offset += 80
        occupations = struct.unpack_from(f"={occupied_count}d", payload, offset)
        offset += 8 * occupied_count
        eigenvalues = struct.unpack_from(f"={occupied_count}d", payload, offset)
        offset += 8 * occupied_count
        receipts = struct.unpack_from("=5d", payload, offset)
        offset += 40
        coefficient_count = 2 * global_count * occupied_count
        coefficients = memoryview(payload)[offset:offset + 8 * coefficient_count].cast("d")
        if (not all(math.isfinite(value) and value >= 0.0 for value in occupations)
                or not all(math.isfinite(value) for value in eigenvalues)
                or not all(math.isfinite(value) and value >= 0.0 for value in receipts)
                or not all(math.isfinite(value) for value in coefficients)):
            raise RuntimeError(f"occupied checkpoint numeric payload is not finite: {path}")
        return ((version, global_count, occupied_count), fingerprints, occupations,
                eigenvalues, receipts, coefficients)

    clean = read(clean_path)
    reuse = read(reuse_path)
    if clean[0] != reuse[0]:
        raise RuntimeError("occupied checkpoint dimensions changed on reuse")
    # Fingerprints 0:9 describe catalog, basis, provenance, and operator.  Index
    # 9 is a bitwise hash of the floating state and is reported, not equated.
    if clean[1][:9] != reuse[1][:9]:
        raise RuntimeError("occupied checkpoint representation fingerprints changed on reuse")
    maximums = [
        max((abs(a - b) for a, b in zip(clean[index], reuse[index])), default=0.0)
        for index in (2, 3, 4)
    ]
    if not all(math.isfinite(value) for value in maximums):
        raise RuntimeError("occupied checkpoint spectral comparison is not finite")
    global_count, occupied_count = clean[0][1:]
    # The writer emits one contiguous row_values(:) vector per global row,
    # rather than writing coefficients(:,:) as one Fortran matrix record.
    # The stream is consequently row-major here.
    clean_coefficients = np.frombuffer(clean[5], dtype=np.float64).view(np.complex128).reshape(
        global_count, occupied_count)
    reuse_coefficients = np.frombuffer(reuse[5], dtype=np.float64).view(np.complex128).reshape(
        global_count, occupied_count)
    coefficient_difference = float(np.max(np.abs(clean_coefficients - reuse_coefficients), initial=0.0))
    if not math.isfinite(coefficient_difference):
        raise RuntimeError("occupied checkpoint coefficient comparison is not finite")

    # Compare rho = C diag(f) C^H directly, in bounded row blocks.  This is
    # invariant to phases and to rotations inside equally occupied degenerate
    # blocks, but rejects a change of occupied subspace or an illicit rotation
    # between states with different occupations.
    clean_occupations = np.asarray(clean[2])
    reuse_occupations = np.asarray(reuse[2])
    density_difference_squared = 0.0
    clean_density_squared = 0.0
    reuse_density_squared = 0.0
    with np.errstate(over="ignore", invalid="ignore"):
        for first in range(0, global_count, 128):
            last = min(first + 128, global_count)
            clean_rows = (clean_coefficients[first:last, :] * clean_occupations) @ clean_coefficients.conj().T
            reuse_rows = (reuse_coefficients[first:last, :] * reuse_occupations) @ reuse_coefficients.conj().T
            if not np.isfinite(clean_rows).all() or not np.isfinite(reuse_rows).all():
                raise RuntimeError("occupied checkpoint density comparison is not finite")
            density_difference_squared += float(np.vdot(clean_rows - reuse_rows, clean_rows - reuse_rows).real)
            clean_density_squared += float(np.vdot(clean_rows, clean_rows).real)
            reuse_density_squared += float(np.vdot(reuse_rows, reuse_rows).real)
    if not all(math.isfinite(value) and value >= 0.0 for value in (
            density_difference_squared, clean_density_squared, reuse_density_squared)):
        raise RuntimeError("occupied checkpoint density norm comparison is not finite")
    density_scale = max(math.sqrt(clean_density_squared), math.sqrt(reuse_density_squared), np.finfo(float).tiny)
    density_difference = math.sqrt(max(0.0, density_difference_squared)) / density_scale
    comparison = {
        "clean_state_fingerprint": clean[1][9],
        "reuse_state_fingerprint": reuse[1][9],
        "state_fingerprint_bitwise_equal": clean[1][9] == reuse[1][9],
        "occupation_max_abs_difference": maximums[0],
        "eigenvalue_max_abs_difference": maximums[1],
        "receipt_max_abs_difference": maximums[2],
        "coefficient_component_max_abs_difference": coefficient_difference,
        "occupied_density_relative_frobenius_difference": density_difference,
        "numeric_tolerance": STATE_NUMERIC_TOLERANCE,
        "exact_representation_fingerprints": list(clean[1][:9]),
    }
    # Eigenvectors may acquire independent phases or rotate within an exactly
    # degenerate, equally occupied block.  Spectra/occupations and certified
    # solver receipts are physical invariants; raw components remain diagnostic.
    if any(value > STATE_NUMERIC_TOLERANCE for value in maximums):
        raise RuntimeError("occupied checkpoint spectral state changed beyond numerical tolerance")
    if density_difference > STATE_NUMERIC_TOLERANCE:
        raise RuntimeError("occupied checkpoint density changed beyond numerical tolerance")
    return comparison


def continuation_receipts(text: str) -> list[dict[str, int | float | str]]:
    pattern = re.compile(
        rf"\[DG-HYBRID-CONTINUATION\]\s+lambda=\s*({FLOAT_TOKEN})\s+"
        rf"diagnostic_state_lambda=\s*({FLOAT_TOKEN})\s+accepted_cg_steps=(\d+)\s+"
        rf"residual=\s*({FLOAT_TOKEN})\s+orthogonality_defect=\s*({FLOAT_TOKEN})\s+"
        rf"electron_defect=\s*({FLOAT_TOKEN})\s+rayleigh_energy_trace=\s*({FLOAT_TOKEN})\s+"
        rf"scaled_interface_action_norm=\s*({FLOAT_TOKEN})\s+measurement_status=(\w+)\s+"
        r"status=(\w+)\s+continuation_fingerprint=(-?\d+)", re.I,
    )
    return [{
        "lambda": parse_finite_float(row.group(1), "continuation lambda", minimum=0.0, maximum=1.0),
        "diagnostic_state_lambda": parse_finite_float(
            row.group(2), "continuation diagnostic lambda", minimum=0.0, maximum=1.0),
        "accepted_cg_steps": int(row.group(3)),
        "residual": parse_finite_float(row.group(4), "continuation residual", minimum=0.0),
        "orthogonality_defect": parse_finite_float(
            row.group(5), "continuation orthogonality defect", minimum=0.0,
            maximum=TERMINAL_SOLVER_TOLERANCE),
        "electron_defect": parse_finite_float(
            row.group(6), "continuation electron defect", minimum=0.0,
            maximum=ELECTRON_TOLERANCE),
        "rayleigh_energy_trace": parse_finite_float(
            row.group(7), "continuation Rayleigh energy trace"),
        "scaled_interface_action_norm": parse_finite_float(
            row.group(8), "continuation scaled interface action norm", minimum=0.0),
        "measurement_status": row.group(9),
        "status": row.group(10),
        "continuation_fingerprint": int(row.group(11)),
    } for row in pattern.finditer(text)]


def parse_run(case: str, run_dir: Path, elapsed: float, return_code: int | None) -> dict:
    log_path = run_dir / "run.log"
    text = log_path.read_text(errors="replace")
    if return_code is not None and return_code != 0:
        raise RuntimeError(f"{case} exited with {return_code}: {log_path}")
    require_mpi_completion(text, EXPECTED_MPI_SIZE, case)
    seed = one_match(
        r"\[DG-DC-SEED\]\s+mode=read\s+publication_id=(-?\d+)\s+"
        r"scf_skipped=([TF])\s+mpi_size=(\d+)\s+mapping_fingerprint=(-?\d+)",
        text, "exact DC seed receipt")
    wf = one_match(
        r"\[DG-FRAGMENT-WF\]\s+mode=auto\s+checkpoint_hit=([TF])\s+"
        r"publication_id=(-?\d+)\s+reason=(.*)$", text, "fragment-WF receipt")
    projected = one_match(
        r"\[DG-HYBRID-DIVIDED\]\s+projected_basis_fingerprint=(-?\d+)",
        text, "projected-basis fingerprint")
    terminal = one_match(
        r"\[OW-GS\]\s+fixed-density/non-self-consistent divided WF\+PW LCFO solved once\s+"
        rf"residual=\s*({FLOAT_TOKEN})\s+orthogonality=\s*({FLOAT_TOKEN})\s+"
        rf"projector=\s*({FLOAT_TOKEN})\s+electron_defect=\s*({FLOAT_TOKEN})",
        text, "terminal fixed-density LCFO receipt")
    continuation = continuation_receipts(text)
    temperatures = [parse_finite_float(value, f"{case} Schwarz temperature", minimum=0.0)
                    for value in re.findall(
        rf"\[DG-HYBRID-SCHWARZ\].*?temperature=\s*({FLOAT_TOKEN})", text)]
    occupied = run_dir / "overlapping_wannier_occupied.chk"
    if not occupied.is_file():
        raise RuntimeError(f"missing occupied checkpoint: {occupied}")
    w90 = wannier90_receipts(run_dir)
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
    if not math.isfinite(elapsed) or elapsed < 0.0:
        raise RuntimeError(f"{case} elapsed time is not finite and nonnegative")
    identity_values = tuple(int(value) for value in (
        seed.group(1), seed.group(4), wf.group(2), projected.group(1)))
    if any(value == 0 for value in identity_values):
        raise RuntimeError(f"{case} contains a zero publication or fingerprint identity")
    evidence = {
        "case": case, "return_code": return_code, "elapsed_seconds": elapsed,
        "completion_receipts": EXPECTED_MPI_SIZE,
        "dc_publication_id": int(seed.group(1)), "dc_scf_skipped": seed.group(2) == "T",
        "dc_mpi_size": int(seed.group(3)), "mapping_fingerprint": int(seed.group(4)),
        "checkpoint_hit": wf.group(1) == "T", "wf_publication_id": int(wf.group(2)),
        "wf_reason": wf.group(3).strip(), "projected_basis_fingerprint": int(projected.group(1)),
        "wannier90": w90, "wannier90_iterations": [int(item["iterations"]) for item in w90],
        "random_wannier90_iterations": RANDOM_W90_ITERATIONS,
        "continuation": continuation,
        "continuation_fingerprints": [item["continuation_fingerprint"] for item in continuation],
        "terminal_residual": terminal_residual,
        "terminal_orthogonality": terminal_orthogonality,
        "terminal_projector": terminal_projector,
        "terminal_electron_defect": terminal_electron_defect,
        "schwarz_temperatures_kelvin": temperatures,
        "occupied_checkpoint_fingerprint": validate_occupied_checkpoint(occupied),
        "occupied_checkpoint_path": str(occupied),
        "ordinary_dc_scf_executed": "DC #SCF =" in text,
        "terminal_lcfo_count": text.count(
            "[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once"),
        "post_lcfo_density_updates": text[terminal.start():].count("[DG-HYBRID-DIVIDED-SCF] iteration="),
        "invalid_number_seen": re.search(
            r"(?<![A-Za-z])(?:nan|[+\-]?inf)(?![A-Za-z])", text, re.I) is not None,
    }
    exact_seed = (
        evidence["dc_publication_id"] == EXPECTED_PUBLICATION_ID and evidence["dc_scf_skipped"]
        and evidence["dc_mpi_size"] == EXPECTED_MPI_SIZE
        and evidence["mapping_fingerprint"] == EXPECTED_MAPPING_FINGERPRINT
        and not evidence["ordinary_dc_scf_executed"])
    valid_schedule = (
        len(continuation) == 6
        and all(math.isclose(item["lambda"], 0.2 * index, abs_tol=1.0e-14)
                for index, item in enumerate(continuation))
        and all(item["measurement_status"].lower() == "valid"
                and item["status"].lower() == "accepted"
                and 1 <= item["accepted_cg_steps"] <= 3 for item in continuation))
    if not exact_seed:
        raise RuntimeError(f"{case} did not reuse the exact authoritative DC seed")
    if not valid_schedule:
        raise RuntimeError(f"{case} did not complete the six-point fixed-density DG schedule")
    if not temperatures or any(value != 300.0 for value in temperatures):
        raise RuntimeError(f"{case} changed the 300 K occupation contract")
    if evidence["terminal_lcfo_count"] != 1 or evidence["post_lcfo_density_updates"] != 0:
        raise RuntimeError(f"{case} did not perform exactly one terminal LCFO without density updates")
    if evidence["invalid_number_seen"]:
        raise RuntimeError(f"{case} produced a nonfinite diagnostic")
    return evidence


def run_case(binary: Path, run_dir: Path, input_text: str, fixture: Path, timeout: int) -> dict:
    run_dir.mkdir(parents=True)
    (run_dir / "inputfile").write_text(input_text)
    shutil.copy2(fixture / "atom.dat", run_dir / "atom.dat")
    root = Path(__file__).resolve().parents[2]
    shutil.copy2(root / "samples/exercise_04_bulkSi_gs/Si_rps.dat", run_dir / "Si_rps.dat")
    print(f"starting {run_dir.name}: {run_dir}", flush=True)
    started = time.monotonic()
    with (run_dir / "inputfile").open("rb") as source, (run_dir / "run.log").open("wb") as log:
        completed = subprocess.run(
            [shutil.which("mpirun") or "mpirun", "-np", str(EXPECTED_MPI_SIZE), str(binary)],
            cwd=run_dir, stdin=source, stdout=log, stderr=subprocess.STDOUT,
            env={**os.environ, "OMP_NUM_THREADS": "1", "OMPI_MCA_rmaps_base_oversubscribe": "1"},
            timeout=timeout)
    elapsed = time.monotonic() - started
    evidence = parse_run(run_dir.name, run_dir, elapsed, completed.returncode)
    print(f"completed {run_dir.name} in {elapsed:.1f} s", flush=True)
    return evidence


def compare_generation_and_reuse(clean: dict, reuse: dict) -> dict:
    if clean["checkpoint_hit"] or not reuse["checkpoint_hit"]:
        raise RuntimeError("clean/reuse fragment-WF cache decisions are incorrect")
    expected_fragments = list(range(1, EXPECTED_MPI_SIZE + 1))
    if ([int(item["fragment_id"]) for item in clean["wannier90"]] != expected_fragments
            or reuse["wannier90"]):
        raise RuntimeError("Wannier90 artifacts do not prove eight clean builds and zero reuse calls")
    exact_fields = (
        "dc_publication_id", "mapping_fingerprint", "wf_publication_id",
        "projected_basis_fingerprint", "continuation_fingerprints")
    changed = [field for field in exact_fields if clean[field] != reuse[field]]
    if changed:
        raise RuntimeError("exact reuse changed deterministic fingerprints: " + ", ".join(changed))
    float_fields = (
        "terminal_residual", "terminal_orthogonality", "terminal_projector",
        "terminal_electron_defect")
    if any(not math.isfinite(float(run[field])) for run in (clean, reuse) for field in float_fields):
        raise RuntimeError("terminal observables are not finite")
    inconsistent = [field for field in float_fields if not math.isclose(
        clean[field], reuse[field], rel_tol=2.0e-10, abs_tol=1.0e-10)]
    if inconsistent:
        raise RuntimeError("reuse changed terminal observables: " + ", ".join(inconsistent))
    occupied_comparison = occupied_numeric_comparison(
        Path(clean["occupied_checkpoint_path"]), Path(reuse["occupied_checkpoint_path"]))
    iterations = clean["wannier90_iterations"]
    return {
        "exact_fingerprint_fields": list(exact_fields),
        "tolerance_consistent_float_fields": list(float_fields),
        "occupied_state_numeric_comparison": occupied_comparison,
        "scdm_vs_random_iteration_delta": [
            observed - baseline for observed, baseline in zip(iterations, RANDOM_W90_ITERATIONS)],
        "scdm_vs_random_total_iteration_ratio": sum(iterations) / sum(RANDOM_W90_ITERATIONS),
        "historical_interface_diagnostic": {
            "lambda_zero_residual": 1.7943383143e2,
            "lambda_one_residual": 2.0819639807e3,
            "interpretation": "interface term dominates residual growth; occupation unit bug is separate"}}


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--binary", type=Path, default=root / "build-hybrid-release/salmon")
    parser.add_argument("--seed-directory", type=Path, default=DEFAULT_SEED_DIRECTORY)
    parser.add_argument("--result-dir", type=Path)
    parser.add_argument("--timeout", type=int, default=14400)
    parser.add_argument("--analyze-existing", action="store_true",
                        help="analyze an already completed clean/reuse result without launching MPI")
    args = parser.parse_args()
    if args.mpi_ranks != EXPECTED_MPI_SIZE:
        raise RuntimeError("Si64 SCDM/reuse validation requires exactly 8 MPI ranks")
    result = (args.result_dir or Path("/tmp") /
              f"si64-task6-scdm-reuse-{time.strftime('%Y%m%d-%H%M%S')}").resolve()
    if result.exists() and not args.analyze_existing:
        raise RuntimeError(f"result directory must be fresh: {result}")
    if not result.exists() and args.analyze_existing:
        raise RuntimeError(f"existing result directory is missing: {result}")
    result.mkdir(parents=True, exist_ok=args.analyze_existing)
    wf_checkpoint = result / "fragment-wf-checkpoint"
    clean_dir = result / "01-clean-scdm"
    reuse_dir = result / "02-exact-reuse"
    if args.analyze_existing:
        elapsed = lambda directory: max(
            0.0,
            (directory / "overlapping_wannier_occupied.chk").stat().st_mtime
            - (directory / "inputfile").stat().st_mtime)
        clean = parse_run(clean_dir.name, clean_dir, elapsed(clean_dir), None)
        reuse = parse_run(reuse_dir.name, reuse_dir, elapsed(reuse_dir), None)
        expected_seed = {
            "publication_id": clean["dc_publication_id"],
            "mpi_size": clean["dc_mpi_size"],
            "mapping_fingerprint": clean["mapping_fingerprint"],
        }
    else:
        binary = args.binary.resolve(strict=True)
        seed_directory = args.seed_directory.resolve(strict=True)
        expected_seed = validate_seed_directory(seed_directory)
        fixture = root / "tests/dg/data/si64_overlapping_wannier_rt"
        template = (fixture / "input_hybrid_divided_lcfo.in").read_text()
        input_text = render_input(template, seed_directory, wf_checkpoint)
        clean = run_case(binary, clean_dir, input_text, fixture, args.timeout)
        reuse = run_case(binary, reuse_dir, input_text, fixture, args.timeout)
    comparison = compare_generation_and_reuse(clean, reuse)
    gauge_contracts = fragment_wf_contract_receipts(wf_checkpoint, clean["wf_publication_id"])
    evidence = {"mpi_ranks": EXPECTED_MPI_SIZE, "omp_threads": 1,
                "expected_seed": expected_seed, "wf_checkpoint_directory": str(wf_checkpoint),
                "fragment_wf_gauge_contracts": gauge_contracts,
                "runs": [clean, reuse], "comparison": comparison}
    (result / "si64_scdm_reuse_evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

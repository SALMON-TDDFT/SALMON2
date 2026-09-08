#!/usr/bin/env python3
"""Contract checks for the long Si64 SCDM generation/reuse runner."""

from __future__ import annotations

import importlib.util
import hashlib
import math
import struct
import tempfile
from pathlib import Path

import run_dg_fragment_wf_production_smoke as SMOKE


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "tests/dg/run_dg_hybrid_si64_scdm_reuse.py"
SPEC = importlib.util.spec_from_file_location("si64_scdm_runner", RUNNER)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


template = (ROOT / "tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in").read_text()
rendered = MODULE.render_input(template, Path("/seed"), Path("/wf"))
for token in (
    "dg_dc_seed_mode='read'",
    "dg_dc_seed_directory='/seed'",
    "dg_fragment_wf_checkpoint_mode='auto'",
    "dg_fragment_wf_checkpoint_directory='/wf'",
    "dg_fragment_w90_initial_projection='scdm'",
    "yn_dg_hybrid_divided_scf='y'",
):
    assert token.replace(" ", "").lower() in rendered.replace(" ", "").lower(), token

with tempfile.TemporaryDirectory(prefix="si64-scdm-runner-") as name:
    root = Path(name)
    seed = root / "seed"
    seed.mkdir()
    manifest = seed / "dg_dc_seed.manifest"
    manifest.write_bytes(struct.pack(
        "=32siiq", b"SALMON_DG_DC_SEED_MANIFEST_V1".ljust(32), 1,
        MODULE.EXPECTED_MPI_SIZE, MODULE.EXPECTED_PUBLICATION_ID,
    ))
    publication_hex = f"{MODULE.EXPECTED_PUBLICATION_ID:016X}"
    synthetic_hashes = {"dg_dc_seed.manifest": hashlib.sha256(manifest.read_bytes()).hexdigest()}
    for rank in range(MODULE.EXPECTED_MPI_SIZE):
        shard = seed / f"dg_dc_seed.{publication_hex}.rank{rank:08d}.shard"
        shard.write_bytes(struct.pack(
            "=32siiiiqqqq", b"SALMON_DG_DC_SEED_SHARD_V1".ljust(32), 1,
            MODULE.EXPECTED_MPI_SIZE, rank, rank + 1, MODULE.EXPECTED_PUBLICATION_ID,
            97, MODULE.EXPECTED_MAPPING_FINGERPRINT, 101,
        ))
        synthetic_hashes[f"rank{rank:08d}"] = hashlib.sha256(shard.read_bytes()).hexdigest()
    authoritative_hashes = MODULE.SEED_FILE_SHA256
    MODULE.SEED_FILE_SHA256 = synthetic_hashes
    assert MODULE.validate_seed_directory(seed)["mapping_fingerprint"] == MODULE.EXPECTED_MAPPING_FINGERPRINT
    damaged = seed / f"dg_dc_seed.{publication_hex}.rank00000003.shard"
    damaged.write_bytes(damaged.read_bytes() + b"damage")
    try:
        MODULE.validate_seed_directory(seed)
    except RuntimeError as error:
        assert "content changed" in str(error)
    else:
        raise AssertionError("damaged exact DC seed must be rejected")
    MODULE.SEED_FILE_SHA256 = authoritative_hashes

    for fragment, iteration, seconds in ((1, 11, 3.25), (8, 19, 7.5)):
        wout = root / f"dgfw-attempt-1/f{fragment:06d}/g00000001/w.wout"
        wout.parent.mkdir(parents=True)
        wout.write_text(
            f" {iteration:6d} -0.100E-09  1.000000  {seconds:8.2f}  <-- CONV\n"
            " <<< Wannierisation convergence criteria satisfied >>>\n"
            f" Total Execution Time {seconds:12.3f} (sec)\n"
        )
    receipts = MODULE.wannier90_receipts(root)
    assert receipts == [
        {"fragment_id": 1, "iterations": 11, "wall_seconds": 3.25},
        {"fragment_id": 8, "iterations": 19, "wall_seconds": 7.5},
    ]

    duplicate_wout = root / "dgfw-attempt-2/f000001/g00000002/w.wout"
    duplicate_wout.parent.mkdir(parents=True)
    duplicate_wout.write_text(
        "     12 -0.100E-09  1.000000      3.50  <-- CONV\n"
        " <<< Wannierisation convergence criteria satisfied >>>\n"
        " Total Execution Time        3.500 (sec)\n"
    )
    try:
        MODULE.wannier90_receipts(root)
    except RuntimeError as error:
        assert "duplicate" in str(error).lower()
    else:
        raise AssertionError("duplicate Wannier90 fragment receipt must be rejected")
    duplicate_wout.unlink()

    complete_log = "\n".join(["end SALMON"] * MODULE.EXPECTED_MPI_SIZE) + "\n"
    MODULE.require_mpi_completion(complete_log, MODULE.EXPECTED_MPI_SIZE, "complete fixture")
    try:
        MODULE.require_mpi_completion(
            "[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once "
            "residual=1.0E-13 orthogonality=1.0E-14 projector=1.0E-14 electron_defect=1.0E-10\n",
            MODULE.EXPECTED_MPI_SIZE,
            "truncated-after-terminal fixture",
        )
    except RuntimeError as error:
        assert "completion" in str(error).lower()
    else:
        raise AssertionError("log truncated after terminal receipt must be rejected")

    try:
        MODULE.parse_finite_float("1e999", "overflow fixture")
    except RuntimeError as error:
        assert "finite" in str(error).lower()
    else:
        raise AssertionError("overflow scalar must be rejected")

    parser_run = root / "parser-run"
    parser_run.mkdir()
    (parser_run / "overlapping_wannier_occupied.chk").write_bytes(b"checkpoint fixture")
    continuation_lines = "\n".join(
        "[DG-HYBRID-CONTINUATION] "
        f"lambda={0.2 * index:.1f} diagnostic_state_lambda={0.2 * index:.1f} "
        "accepted_cg_steps=3 residual=1.0E+02 orthogonality_defect=1.0E-14 "
        "electron_defect=1.0E-14 rayleigh_energy_trace=-1.0E+02 "
        f"scaled_interface_action_norm={float(index):.1f} measurement_status=valid "
        f"status=accepted continuation_fingerprint={101 + index}"
        for index in range(6)
    )
    parser_log = (
        f"[DG-DC-SEED] mode=read publication_id={MODULE.EXPECTED_PUBLICATION_ID} "
        f"scf_skipped=T mpi_size={MODULE.EXPECTED_MPI_SIZE} "
        f"mapping_fingerprint={MODULE.EXPECTED_MAPPING_FINGERPRINT}\n"
        "[DG-FRAGMENT-WF] mode=auto checkpoint_hit=T publication_id=303 reason=compatible checkpoint hit\n"
        "[DG-HYBRID-DIVIDED] projected_basis_fingerprint=404\n"
        + continuation_lines
        + "\n[DG-HYBRID-SCHWARZ] electron_defect=1.0E-14 temperature=3.0E+02\n"
        "[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once "
        "residual=2.0E-13 orthogonality=1.0E-14 projector=1.0E-14 electron_defect=1.0E-10\n"
        + complete_log
    )
    saved_validator = MODULE.validate_occupied_checkpoint
    MODULE.validate_occupied_checkpoint = lambda path: 505
    try:
        (parser_run / "run.log").write_text(parser_log)
        assert MODULE.parse_run("valid-parser-fixture", parser_run, 1.0, None)[
            "completion_receipts"
        ] == MODULE.EXPECTED_MPI_SIZE

        terminal_end = parser_log.index("\n", parser_log.index("[OW-GS]")) + 1
        (parser_run / "run.log").write_text(parser_log[:terminal_end])
        try:
            MODULE.parse_run("truncated-after-terminal", parser_run, 1.0, None)
        except RuntimeError as error:
            assert "completion" in str(error).lower()
        else:
            raise AssertionError("parse_run accepted a log truncated after terminal receipt")

        overflow_log = parser_log.replace(
            "residual=2.0E-13 orthogonality=", "residual=1e999 orthogonality=", 1
        )
        (parser_run / "run.log").write_text(overflow_log)
        try:
            MODULE.parse_run("overflow-terminal", parser_run, 1.0, None)
        except RuntimeError as error:
            assert "finite" in str(error).lower()
        else:
            raise AssertionError("parse_run accepted a finite-overflow terminal scalar")
    finally:
        MODULE.validate_occupied_checkpoint = saved_validator

    smoke_run = root / "smoke-parser-run"
    smoke_run.mkdir()
    (smoke_run / "overlapping_wannier_occupied.chk").write_bytes(b"checkpoint fixture")
    smoke_schwarz = "\n".join(
        "[DG-HYBRID-SCHWARZ] "
        f"epoch={epoch} neighbor_exchanges=7 accepted_cg_steps=3 "
        "common_extensions=0 residual=1.0E+02 electron_defect=1.0E-14 temperature=3.0E+02"
        for epoch in range(1, 7)
    )
    smoke_log = (
        "[DG-FRAGMENT-WF] mode=auto checkpoint_hit=T publication_id=606 reason=compatible checkpoint hit\n"
        "[DG-HYBRID-DIVIDED] projected_basis_fingerprint=707\n"
        + smoke_schwarz + "\n"
        "[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once "
        "residual=2.0E-13 orthogonality=1.0E-14 projector=1.0E-14 electron_defect=1.0E-10\n"
        "[DG-DC-SEED] mode=auto publication_id=808 scf_skipped=T mpi_size=8 "
        "mapping_fingerprint=909\n"
        + complete_log
    )
    saved_smoke_validator = SMOKE.validate_occupied_checkpoint
    SMOKE.validate_occupied_checkpoint = lambda path: 1001
    try:
        (smoke_run / "run.log").write_text(smoke_log)
        assert SMOKE.parse_evidence("valid-smoke-parser", smoke_run, 1.0, None)[
            "completion_receipts"
        ] == SMOKE.RANKS
        terminal_end = smoke_log.index("\n", smoke_log.index("[OW-GS]")) + 1
        (smoke_run / "run.log").write_text(smoke_log[:terminal_end])
        try:
            SMOKE.parse_evidence("truncated-smoke", smoke_run, 1.0, None)
        except RuntimeError as error:
            assert "completion" in str(error).lower()
        else:
            raise AssertionError("smoke parser accepted a log truncated after terminal receipt")
        (smoke_run / "run.log").write_text(smoke_log.replace(
            "residual=2.0E-13 orthogonality=", "residual=1e999 orthogonality=", 1
        ))
        try:
            SMOKE.parse_evidence("overflow-smoke", smoke_run, 1.0, None)
        except RuntimeError as error:
            assert "finite" in str(error).lower()
        else:
            raise AssertionError("smoke parser accepted a finite-overflow terminal scalar")
    finally:
        SMOKE.validate_occupied_checkpoint = saved_smoke_validator

    duplicate_smoke_wout_a = root / "smoke-w90-a/f000001/g00000001/w.wout"
    duplicate_smoke_wout_b = root / "smoke-w90-b/f000001/g00000002/w.wout"
    duplicate_smoke_wout_a.parent.mkdir(parents=True)
    duplicate_smoke_wout_b.parent.mkdir(parents=True)
    duplicate_smoke_wout_a.write_text("a")
    duplicate_smoke_wout_b.write_text("b")
    try:
        SMOKE.wannier90_fragment_ids(root)
    except RuntimeError as error:
        assert "duplicate" in str(error).lower()
    else:
        raise AssertionError("smoke parser accepted duplicate Wannier90 fragment IDs")

    def write_occupied(path: Path, state_fingerprint: int, delta: float, phase: float) -> None:
        path.write_bytes(
            b"SALMON_DG_OCC02 "
            + struct.pack("=iii", 2, 1, 1)
            + struct.pack("=10Q", 11, 12, 13, 14, 15, 16, 17, 18, 19, state_fingerprint)
            + struct.pack("=d", 2.0 + delta)
            + struct.pack("=d", -0.5 + delta)
            + struct.pack("=5d", *(delta for _ in range(5)))
            + struct.pack("=2d", phase * (1.0 + delta), delta)
        )

    clean = root / "clean.chk"
    reuse = root / "reuse.chk"
    write_occupied(clean, 101, 0.0, 1.0)
    write_occupied(reuse, 202, 1.0e-12, -1.0)
    comparison = MODULE.occupied_numeric_comparison(clean, reuse)
    assert not comparison["state_fingerprint_bitwise_equal"]
    assert comparison["coefficient_component_max_abs_difference"] > MODULE.STATE_NUMERIC_TOLERANCE
    assert comparison["occupied_density_relative_frobenius_difference"] < MODULE.STATE_NUMERIC_TOLERANCE

    nonfinite_spectral = root / "nonfinite-spectral.chk"
    write_occupied(nonfinite_spectral, 203, float("nan"), 1.0)
    try:
        MODULE.occupied_numeric_comparison(clean, nonfinite_spectral)
    except RuntimeError as error:
        assert "finite" in str(error).lower()
    else:
        raise AssertionError("nonfinite occupied spectrum must be rejected")

    nonfinite_density = root / "nonfinite-density.chk"
    write_occupied(nonfinite_density, 204, 0.0, float("inf"))
    try:
        MODULE.occupied_numeric_comparison(clean, nonfinite_density)
    except RuntimeError as error:
        assert "finite" in str(error).lower()
    else:
        raise AssertionError("nonfinite occupied-density coefficients must be rejected")

    def write_occupied_vector(path: Path, state_fingerprint: int, vector: tuple[complex, complex]) -> None:
        coefficient_values = tuple(value for coefficient in vector
                                   for value in (coefficient.real, coefficient.imag))
        path.write_bytes(
            b"SALMON_DG_OCC02 "
            + struct.pack("=iii", 2, 2, 1)
            + struct.pack("=10Q", 11, 12, 13, 14, 15, 16, 17, 18, 19, state_fingerprint)
            + struct.pack("=d", 2.0)
            + struct.pack("=d", -0.5)
            + struct.pack("=5d", *(0.0 for _ in range(5)))
            + struct.pack("=4d", *coefficient_values)
        )

    orthogonal = root / "orthogonal.chk"
    write_occupied_vector(clean, 301, (1.0 + 0.0j, 0.0 + 0.0j))
    write_occupied_vector(orthogonal, 302, (0.0 + 0.0j, 1.0 + 0.0j))
    try:
        MODULE.occupied_numeric_comparison(clean, orthogonal)
    except RuntimeError as error:
        assert "density" in str(error)
    else:
        raise AssertionError("different occupied subspaces must be rejected")

    def write_occupied_matrix(
        path: Path, state_fingerprint: int, matrix: tuple[tuple[complex, ...], ...]
    ) -> None:
        global_count = len(matrix)
        occupied_count = len(matrix[0])
        coefficient_values = tuple(
            value for row in matrix for coefficient in row
            for value in (coefficient.real, coefficient.imag)
        )
        path.write_bytes(
            b"SALMON_DG_OCC02 "
            + struct.pack("=iii", 2, global_count, occupied_count)
            + struct.pack("=10Q", 11, 12, 13, 14, 15, 16, 17, 18, 19, state_fingerprint)
            + struct.pack(f"={occupied_count}d", *(1.5 for _ in range(occupied_count)))
            + struct.pack(f"={occupied_count}d", *(-0.5 for _ in range(occupied_count)))
            + struct.pack("=5d", *(0.0 for _ in range(5)))
            + struct.pack(f"={2 * global_count * occupied_count}d", *coefficient_values)
        )

    clean_matrix = (
        (1.0 + 0.0j, 0.25 + 0.5j),
        (0.0 + 0.75j, 1.25 + 0.0j),
        (0.5 - 0.25j, -0.5 + 0.125j),
    )
    cosine = math.cos(0.37)
    sine = math.sin(0.37)
    rotated_matrix = tuple((
        (row[0] * cosine + row[1] * sine) * 1.0j,
        (-row[0] * sine + row[1] * cosine) * (-1.0 + 0.0j),
    ) for row in clean_matrix)
    rotated = root / "rotated.chk"
    write_occupied_matrix(clean, 401, clean_matrix)
    write_occupied_matrix(rotated, 402, rotated_matrix)
    rotated_comparison = MODULE.occupied_numeric_comparison(clean, rotated)
    assert rotated_comparison["coefficient_component_max_abs_difference"] > 1.0
    assert rotated_comparison["occupied_density_relative_frobenius_difference"] < 1.0e-14

    publication = 303
    checkpoint = root / "fragment-wf"
    checkpoint.mkdir()
    for rank in range(8):
        header = struct.pack(
            "=32s5iq10i16s13q",
            b"SALMON_DG_FRAGMENT_WF_SHARD_V1".ljust(32), 1, 4, 8, 8, 16, publication,
            1, 8, rank, rank + 1, 1, 1, 400, 400, 1, 400, b"scdm".ljust(16),
            MODULE.EXPECTED_MAPPING_FINGERPRINT, MODULE.EXPECTED_PUBLICATION_ID,
            21, 22, 23, 24, 25, 26, 27, 28, 100 + rank, 200 + rank, 29,
        )
        path = checkpoint / (
            f"dg_fragment_wf.publication-{publication:016X}.rank-{rank:06d}.bin"
        )
        path.write_bytes(header + struct.pack("=d", 1.0e-12))
    gauge = MODULE.fragment_wf_contract_receipts(checkpoint, publication)
    assert len(gauge) == 8 and all(item["seed_reconstruction_defect"] == 1.0e-12 for item in gauge)

    duplicate_rank_path = checkpoint / (
        f"dg_fragment_wf.publication-{publication:016X}.rank-000007.bin"
    )
    duplicate_rank_payload = bytearray(duplicate_rank_path.read_bytes())
    duplicate_rank_values = list(struct.unpack("=32s5iq10i16s13q", duplicate_rank_payload[:-8]))
    duplicate_rank_values[9] = 6
    duplicate_rank_values[10] = 7
    duplicate_rank_path.write_bytes(
        struct.pack("=32s5iq10i16s13q", *duplicate_rank_values) + duplicate_rank_payload[-8:]
    )
    try:
        MODULE.fragment_wf_contract_receipts(checkpoint, publication)
    except RuntimeError as error:
        assert "rank-fragment" in str(error).lower()
    else:
        raise AssertionError("duplicate rank-fragment gauge receipt must be rejected")

print("Si64 SCDM generation/reuse runner contract: PASS")

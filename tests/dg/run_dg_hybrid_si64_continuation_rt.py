#!/usr/bin/env python3
"""Run and certify the Si64 DG-continuation to zero-field hybrid-RT route."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import selectors
import signal
import shutil
import subprocess
import tempfile
import time
from collections.abc import Callable
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests/dg/data/si64_overlapping_wannier_rt"
FINAL_TOLERANCE = 1.0e-7
ELECTRON_TOLERANCE = 1.0e-8
SYMMETRY_TOLERANCE = 1.0e-10
RETAINED_BASIS_RECEIPT = "[HYBRID-RETAINED-BASIS-SYMMETRY]"
GS_RECEIPT = "[HYBRID-GS-ACCEPTANCE]"
LCFO_WINDOW_RECEIPT = "[HYBRID-LCFO-WINDOW]"
LCFO_SYMMETRY_RECEIPT = "[HYBRID-LCFO-SYMMETRY]"
RT_BASIS_RECEIPT = "[HYBRID-RT-BASIS]"
LEGACY_WARNING = "[HYBRID-SYMMETRY-COMPATIBILITY-WARNING]"
RECEIPT_FIELD = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=[ \t]*([^\s]+)")
HANDOFF_RECEIPT = "[HYBRID-RT-HANDOFF]"
STATIONARITY_RECEIPT = "[HYBRID-RT-STATIONARITY]"


def receipt_records(log: str, prefix: str) -> list[dict[str, str]]:
    return [dict(RECEIPT_FIELD.findall(line)) for line in log.splitlines() if prefix in line]


def printed_receipt_close(value: float, expected: float, *scale_values: float) -> bool:
    scale = max(1.0, abs(value), abs(expected), *(abs(item) for item in scale_values))
    return abs(value - expected) <= 2.0e-8 * scale


def certify(
    gs_log: str,
    rt_log: str,
    expected_steps: int = 1,
    expected_delta_e: float | None = None,
) -> dict[str, float | int]:
    if expected_steps < 1:
        raise ValueError("expected RT step count must be positive")
    retained_matches = receipt_records(gs_log, RETAINED_BASIS_RECEIPT)
    gs_matches = receipt_records(gs_log, GS_RECEIPT)
    window_matches = receipt_records(gs_log, LCFO_WINDOW_RECEIPT)
    symmetry_matches = receipt_records(gs_log, LCFO_SYMMETRY_RECEIPT)
    rt_basis_matches = receipt_records(gs_log, RT_BASIS_RECEIPT)
    handoff_matches = receipt_records(rt_log, HANDOFF_RECEIPT)
    stationarity_matches = receipt_records(rt_log, STATIONARITY_RECEIPT)
    if not retained_matches:
        raise ValueError("missing retained-basis symmetry diagnostic")
    if not gs_matches:
        raise ValueError("missing complete lambda-one DG continuation receipt")
    if not window_matches:
        raise ValueError("missing certified LCFO energy-window receipt")
    if not symmetry_matches:
        raise ValueError("missing final LCFO physical-symmetry receipt")
    if not rt_basis_matches:
        raise ValueError("missing certified RT-basis receipt")
    if not handoff_matches:
        raise ValueError("missing exact GS-to-RT payload handoff receipt")
    if not stationarity_matches:
        raise ValueError("missing zero-field hybrid RT stationarity receipt")
    retained = retained_matches[-1]
    try:
        retained_closed = int(retained["closed"])
        retained_defect = float(retained["defect"])
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete retained-basis symmetry diagnostic") from error
    if retained_closed not in (0, 1) or not math.isfinite(retained_defect) or retained_defect < 0.0:
        raise ValueError("invalid retained-basis symmetry diagnostic")
    gs = gs_matches[-1]
    try:
        flags = tuple(int(gs[name]) for name in ("seed_identity", "lambda_zero", "lambda_one", "final_refresh"))
        requested_target_rank = int(gs["requested_target_rank"])
        extended_target_rank = int(gs["extended_target_rank"])
        residuals = tuple(float(gs[name]) for name in (
            "r_h", "r_rho", "r_t", "r_s", "electron", "symmetry", "real_space",
        ))
        physical_symmetry_defects = tuple(float(gs[name]) for name in (
            "occupied_defect", "target_defect", "target_energy_defect", "density_defect",
        ))
        full_operator_defect = float(gs["full_operator_defect"])
        gs_fingerprint = int(gs["payload_fingerprint"])
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete or invalid lambda-one DG continuation receipt") from error
    if flags != (1, 1, 1, 1):
        raise ValueError("continuation seed, lambda-zero, lambda-one, or final refresh was not accepted")
    if requested_target_rank < 1 or extended_target_rank < requested_target_rank:
        raise ValueError("invalid requested or degeneracy-extended symmetry target rank")
    if any(not math.isfinite(value) or value < 0.0 for value in residuals):
        raise ValueError("invalid final DG continuation residual")
    limits = (FINAL_TOLERANCE,) * 4 + (ELECTRON_TOLERANCE, SYMMETRY_TOLERANCE, FINAL_TOLERANCE)
    if any(value > limit for value, limit in zip(residuals, limits)):
        raise ValueError("final DG continuation residual exceeds its production tolerance")
    if any(not math.isfinite(value) or value < 0.0 or value > SYMMETRY_TOLERANCE
           for value in physical_symmetry_defects):
        raise ValueError("post-LCFO physical symmetry defect exceeds its production tolerance")
    if not math.isfinite(full_operator_defect) or full_operator_defect < 0.0:
        raise ValueError("invalid full retained-Hamiltonian symmetry diagnostic")
    window = window_matches[-1]
    try:
        window_mode = window["mode"]
        compatibility = int(window["compatibility"])
        delta_e = float(window["delta_e"])
        e_homo = float(window["homo"])
        requested_cutoff = float(window["requested_cutoff"])
        requested_rank = int(window["requested_rank"])
        certified_cutoff = float(window["certified_cutoff"])
        certified_rank = int(window["certified_rank"])
        boundary_rank = int(window["boundary_rank"])
        extension_states = int(window["extension_states"])
        extension_energy = float(window["extension_energy"])
        proof_state = int(window["proof_state"])
        proof_energy = float(window["proof_energy"])
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete or invalid certified LCFO energy-window receipt") from error
    window_energies = (delta_e, e_homo, requested_cutoff, certified_cutoff, extension_energy, proof_energy)
    if (window_mode not in ("explicit", "legacy_dynamic") or compatibility not in (0, 1)
            or any(not math.isfinite(value) for value in window_energies)
            or proof_state not in (0, 1)
            or requested_rank < 1 or certified_rank < requested_rank
            or boundary_rank != certified_rank
            or extension_states != certified_rank - requested_rank
            or extension_energy < 0.0
            or not printed_receipt_close(
                extension_energy, max(0.0, certified_cutoff - requested_cutoff),
                certified_cutoff, requested_cutoff,
            )):
        raise ValueError("invalid certified LCFO energy-window receipt")
    if requested_target_rank != requested_rank or extended_target_rank != certified_rank:
        raise ValueError("legacy GS rank receipt disagrees with the certified LCFO window")
    if delta_e == -1.0:
        if window_mode != "legacy_dynamic" or compatibility != 1 or LEGACY_WARNING not in gs_log:
            raise ValueError("legacy dynamic-rank receipt or compatibility warning is missing")
    else:
        if (delta_e < 0.0 or window_mode != "explicit" or compatibility != 0
                or not printed_receipt_close(requested_cutoff, e_homo + delta_e, e_homo, delta_e)
                or proof_state != 1):
            raise ValueError("invalid explicit LCFO energy-window receipt")
    if expected_delta_e is not None:
        if (not math.isfinite(expected_delta_e) or expected_delta_e < 0.0
                or not printed_receipt_close(delta_e, expected_delta_e)):
            raise ValueError("LCFO receipt disagrees with the rendered symmetry energy window")
    symmetry = symmetry_matches[-1]
    try:
        lcfo_symmetry_defects = tuple(float(symmetry[name]) for name in ("occupied", "target", "energy", "density"))
        worst_operation = int(symmetry["worst_operation"])
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete or invalid final LCFO physical-symmetry receipt") from error
    if (worst_operation < 1 or any(not math.isfinite(value) or value < 0.0 or value > SYMMETRY_TOLERANCE
                                   for value in lcfo_symmetry_defects)):
        raise ValueError("final LCFO physical-symmetry defect exceeds its production tolerance")
    rt_basis = rt_basis_matches[-1]
    try:
        construction_rank = int(rt_basis["construction_rank"])
        basis_certified_rank = int(rt_basis["certified_rank"])
        certified_rt_rank = int(rt_basis["rt_rank"])
        localization_spread = float(rt_basis["localization_spread"])
        embedding_fingerprint = int(rt_basis["embedding_fingerprint"])
        operator_covariance = float(rt_basis["operator_covariance"])
        checkpoint_version = int(rt_basis["checkpoint_version"])
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete or invalid certified RT-basis receipt") from error
    if (construction_rank < certified_rank or basis_certified_rank != certified_rank
            or certified_rt_rank != certified_rank):
        raise ValueError("published RT rank differs from the certified LCFO rank")
    if proof_state != int(certified_rank < construction_rank):
        raise ValueError("LCFO proof-state receipt disagrees with the solved construction spectrum")
    proof_boundary = max(requested_cutoff, certified_cutoff) if window_mode == "explicit" else certified_cutoff
    proof_is_strictly_above_boundary = (
        proof_energy > proof_boundary
        and not printed_receipt_close(
            proof_energy, proof_boundary, requested_cutoff, certified_cutoff,
        )
    )
    if ((proof_state == 1 and not proof_is_strictly_above_boundary)
            or (proof_state == 0 and proof_energy != 0.0)):
        raise ValueError("LCFO proof energy does not bracket the certified space")
    if (checkpoint_version != 3 or embedding_fingerprint == 0
            or not math.isfinite(localization_spread) or localization_spread < 0.0
            or not math.isfinite(operator_covariance) or operator_covariance < 0.0
            or operator_covariance > SYMMETRY_TOLERANCE):
        raise ValueError("invalid certified RT-basis publication receipt")
    handoff = handoff_matches[-1]
    try:
        rt_fingerprint = int(handoff["payload_fingerprint"])
        handoff_lcfo_residual = float(handoff["lcfo_residual"])
        handoff_metric_defect = float(handoff["metric_defect"])
        handoff_projector_defect = float(handoff["projector_defect"])
        occupied_rank = int(handoff["occupied_rank"])
        rt_extents = tuple(int(handoff[name]) for name in (
            "certified_rank", "state_rank", "metric_rank", "operator_rank",
            "basis_rank", "coefficient_rows",
        ))
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete or invalid certified RT handoff receipt") from error
    if gs_fingerprint == 0 or rt_fingerprint != gs_fingerprint:
        raise ValueError("GS and RT payload identities differ")
    if occupied_rank < 1 or rt_extents[0] < occupied_rank or rt_extents[0] != certified_rank or \
            any(extent != construction_rank for extent in rt_extents[1:]):
        raise ValueError("RT certification or construction-space extent is inconsistent")
    if any(not math.isfinite(value) or value < 0.0 or value > FINAL_TOLERANCE for value in (
            handoff_lcfo_residual, handoff_metric_defect, handoff_projector_defect)):
        raise ValueError("invalid measured LCFO/metric/projector handoff invariant")
    stationarity_limits = (FINAL_TOLERANCE, FINAL_TOLERANCE, FINAL_TOLERANCE,
                           ELECTRON_TOLERANCE, FINAL_TOLERANCE)
    stationarity_values: list[tuple[float, ...]] = []
    stationarity_steps: list[int] = []
    for record in stationarity_matches:
        try:
            stationarity_steps.append(int(record["step"]))
            values = tuple(float(record[name]) for name in (
                "density", "energy", "projector", "electron", "h_residual",
            ))
        except (KeyError, ValueError) as error:
            raise ValueError("incomplete or invalid zero-field stationarity receipt") from error
        if any(not math.isfinite(value) or value < 0.0 for value in values):
            raise ValueError("invalid zero-field stationarity receipt")
        if any(value > limit for value, limit in zip(values, stationarity_limits)):
            raise ValueError("zero-field stationarity drift exceeds its production tolerance")
        stationarity_values.append(values)
    if stationarity_steps != list(range(1, expected_steps + 1)):
        raise ValueError("zero-field stationarity receipts do not cover every requested RT step")
    stationarity = tuple(max(values[index] for values in stationarity_values) for index in range(5))
    return {
        "payload_fingerprint": gs_fingerprint,
        "retained_basis_closed": retained_closed,
        "retained_basis_defect": retained_defect,
        "requested_target_rank": requested_target_rank,
        "extended_target_rank": extended_target_rank,
        "requested_rank": requested_rank,
        "certified_rank": certified_rank,
        "boundary_rank": boundary_rank,
        "construction_rank": construction_rank,
        "certified_rt_rank": certified_rt_rank,
        "delta_e": delta_e,
        "homo": e_homo,
        "requested_cutoff": requested_cutoff,
        "certified_cutoff": certified_cutoff,
        "extension_energy": extension_energy,
        "proof_state": proof_state,
        "proof_energy": proof_energy,
        "lcfo_occupied_defect": lcfo_symmetry_defects[0],
        "lcfo_target_defect": lcfo_symmetry_defects[1],
        "lcfo_energy_defect": lcfo_symmetry_defects[2],
        "lcfo_density_defect": lcfo_symmetry_defects[3],
        "worst_operation": worst_operation,
        "localization_spread": localization_spread,
        "checkpoint_version": checkpoint_version,
        "embedding_fingerprint": embedding_fingerprint,
        "operator_covariance": operator_covariance,
        "legacy_dynamic_rank": compatibility,
        "r_h": residuals[0],
        "r_rho": residuals[1],
        "r_t": residuals[2],
        "r_s": residuals[3],
        "occupied_defect": physical_symmetry_defects[0],
        "target_defect": physical_symmetry_defects[1],
        "target_energy_defect": physical_symmetry_defects[2],
        "density_defect": physical_symmetry_defects[3],
        "full_operator_defect": full_operator_defect,
        "rt_lcfo_residual": handoff_lcfo_residual,
        "rt_metric_defect": handoff_metric_defect,
        "rt_projector_defect": handoff_projector_defect,
        "rt_certified_rank": rt_extents[0],
        "rt_occupied_rank": occupied_rank,
        "rt_state_rank": rt_extents[1],
        "rt_metric_rank": rt_extents[2],
        "rt_operator_rank": rt_extents[3],
        "rt_basis_rank": rt_extents[4],
        "rt_coefficient_rows": rt_extents[5],
        "stationarity_steps": expected_steps,
        "density_drift": stationarity[0],
        "energy_drift": stationarity[1],
        "projector_drift": stationarity[2],
        "electron_drift": stationarity[3],
        "hamiltonian_residual": stationarity[4],
    }


def parser_self_test() -> None:
    complete_gs = (
        "[HYBRID-RETAINED-BASIS-SYMMETRY] closed=0 defect=2e-2\n"
        "[HYBRID-LCFO-WINDOW] mode=explicit compatibility=0 delta_e=0.5 homo=-0.2 "
        "requested_cutoff=0.3 requested_rank=3 certified_cutoff=0.8 certified_rank=4 "
        "boundary_rank=4 extension_states=1 extension_energy=0.5 proof_state=1 proof_energy=1.1\n"
        "[HYBRID-LCFO-SYMMETRY] occupied=1e-12 target=2e-12 energy=3e-12 "
        "density=4e-12 worst_operation=2\n"
        "[HYBRID-RT-BASIS] construction_rank=7 certified_rank=4 rt_rank=4 localization_spread=1.25 "
        "embedding_fingerprint=271828 operator_covariance=5e-12 checkpoint_version=3\n"
        "[HYBRID-GS-ACCEPTANCE] seed_identity=1 lambda_zero=1 lambda_one=1 "
        "final_refresh=1 requested_target_rank=3 extended_target_rank=4 "
        "r_h=1e-10 r_rho=2e-10 r_t=3e-10 r_s=4e-10 "
        "electron=0e0 symmetry=4e-12 real_space=6e-10 "
        "occupied_defect=1e-12 target_defect=2e-12 target_energy_defect=3e-12 "
        "density_defect=4e-12 full_operator_defect=2e-2 payload_fingerprint=314159\n"
    )
    complete_rt = (
        "[HYBRID-RT-HANDOFF] payload_fingerprint=314159 lcfo_residual=2e-12 metric_defect=2e-12 "
        "projector_defect=2e-12 occupied_rank=2 certified_rank=4 state_rank=7 metric_rank=7 operator_rank=7 "
        "basis_rank=7 coefficient_rows=7\n"
        "[HYBRID-RT-STATIONARITY] step=1 density=1e-12 energy=2e-12 projector=3e-12 "
        "electron=4e-12 h_residual=5e-12\n"
    )
    evidence = certify(complete_gs, complete_rt)
    assert evidence["requested_rank"] == 3
    assert evidence["certified_rank"] == evidence["certified_rt_rank"] == 4
    assert evidence["boundary_rank"] == 4
    assert evidence["construction_rank"] == 7
    assert evidence["delta_e"] == 0.5
    assert evidence["requested_cutoff"] == 0.3
    assert evidence["certified_cutoff"] == 0.8
    assert evidence["proof_energy"] == 1.1
    assert evidence["lcfo_target_defect"] == 2.0e-12
    assert evidence["lcfo_density_defect"] == 4.0e-12
    assert evidence["checkpoint_version"] == 3
    assert evidence["embedding_fingerprint"] == 271828
    assert evidence["operator_covariance"] == 5.0e-12
    assert evidence["rt_certified_rank"] == evidence["certified_rank"]
    assert evidence["rt_certified_rank"] >= evidence["rt_occupied_rank"]
    assert evidence["rt_state_rank"] == evidence["construction_rank"]
    assert evidence["rt_metric_rank"] == evidence["construction_rank"]
    assert evidence["rt_operator_rank"] == evidence["construction_rank"]
    assert evidence["rt_basis_rank"] == evidence["construction_rank"]
    assert evidence["rt_coefficient_rows"] == evidence["construction_rank"]
    fixed_width_evidence = certify(
        complete_gs.replace("=", "=  "),
        complete_rt.replace("=", "=  "),
    )
    assert fixed_width_evidence == evidence
    two_step_rt = complete_rt + (
        "[HYBRID-RT-STATIONARITY] step=2 density=2e-12 energy=3e-12 projector=4e-12 "
        "electron=5e-12 h_residual=6e-12\n"
    )
    two_step_evidence = certify(complete_gs, two_step_rt, expected_steps=2)
    assert two_step_evidence["stationarity_steps"] == 2
    assert two_step_evidence["density_drift"] == 2.0e-12
    for incomplete_rt in (
        complete_rt,
        two_step_rt.replace("step=1", "step=2", 1),
        two_step_rt.replace("step=2 density=2e-12", "step=3 density=2e-12"),
        two_step_rt.replace("step=1 density=1e-12", "step=1 density=1e-2"),
        two_step_rt.replace("step=1 density=1e-12", "step=1 density=nan"),
    ):
        try:
            certify(complete_gs, incomplete_rt, expected_steps=2)
        except ValueError:
            continue
        raise AssertionError("incomplete or non-stationary RT step sequence was accepted")
    certify(complete_gs, complete_rt, expected_delta_e=0.5)
    try:
        certify(complete_gs, complete_rt, expected_delta_e=0.25)
    except ValueError:
        pass
    else:
        raise AssertionError("LCFO receipt disagreed with the rendered energy window but was accepted")
    between_levels_gs = complete_gs.replace(
        "[HYBRID-LCFO-WINDOW] mode=explicit compatibility=0 delta_e=0.5 homo=-0.2 "
        "requested_cutoff=0.3 requested_rank=3 certified_cutoff=0.8 certified_rank=4 "
        "boundary_rank=4 extension_states=1 extension_energy=0.5 proof_state=1 proof_energy=1.1\n",
        "[HYBRID-LCFO-WINDOW] mode=explicit compatibility=0 delta_e=0.234567894 homo=0.123456784 "
        "requested_cutoff=0.358024679 requested_rank=3 certified_cutoff=0.2 certified_rank=3 "
        "boundary_rank=3 extension_states=0 extension_energy=0 proof_state=1 proof_energy=0.8\n",
    ).replace(
        "construction_rank=7 certified_rank=4 rt_rank=4",
        "construction_rank=7 certified_rank=3 rt_rank=3",
    ).replace("extended_target_rank=4", "extended_target_rank=3")
    between_levels_rt = complete_rt.replace(
        "certified_rank=4 state_rank=7 metric_rank=7 operator_rank=7 basis_rank=7 coefficient_rows=7",
        "certified_rank=3 state_rank=7 metric_rank=7 operator_rank=7 basis_rank=7 coefficient_rows=7",
    )
    between_levels = certify(between_levels_gs, between_levels_rt)
    assert between_levels["certified_cutoff"] < between_levels["requested_cutoff"]
    assert between_levels["requested_rank"] == between_levels["certified_rank"] == 3
    proof_at_print_tolerance = between_levels_gs.replace(
        "proof_energy=0.8", "proof_energy=0.358024680",
    )
    try:
        certify(proof_at_print_tolerance, between_levels_rt)
    except ValueError:
        pass
    else:
        raise AssertionError(
            "proof energy indistinguishable from the larger requested/certified cutoff was accepted"
        )
    incomplete_logs = [
        "[OW-GS-DIAGNOSTIC] one_shot_generalized_eigenexa residual=0\n",
        "[OW-GS] fully refreshed DG continuation converged lambda=1\n",
        complete_gs.replace("lambda_zero=1", "lambda_zero=0"),
        complete_gs.replace("r_h=1e-10", "r_h=1e-2"),
        complete_gs.replace("[HYBRID-RETAINED-BASIS-SYMMETRY]", "[MISSING-RETAINED-BASIS-SYMMETRY]"),
        complete_gs.replace(" requested_target_rank=3", ""),
        complete_gs.replace(" extended_target_rank=4", ""),
        complete_gs.replace(" occupied_defect=1e-12", ""),
        complete_gs.replace(" target_defect=2e-12", ""),
        complete_gs.replace(" target_energy_defect=3e-12", ""),
        complete_gs.replace(" density_defect=4e-12", ""),
        complete_gs.replace(" full_operator_defect=2e-2", ""),
        complete_gs.replace("target_defect=2e-12", "target_defect=1e-2"),
        complete_gs.replace("target_defect=2e-12", "target_defect=nan"),
        complete_gs.replace("density_defect=4e-12", "density_defect=inf"),
        complete_gs.replace("extended_target_rank=4", "extended_target_rank=2"),
        complete_gs.replace("requested_target_rank=3", "requested_target_rank=2"),
        complete_gs.replace("extended_target_rank=4", "extended_target_rank=5"),
        complete_gs.replace("[HYBRID-LCFO-WINDOW]", "[MISSING-HYBRID-LCFO-WINDOW]"),
        complete_gs.replace(" requested_cutoff=0.3", ""),
        complete_gs.replace(" certified_rank=4", "", 1),
        complete_gs.replace(" boundary_rank=4", ""),
        complete_gs.replace(" extension_energy=0.5", ""),
        complete_gs.replace(" proof_state=1", ""),
        complete_gs.replace(" proof_energy=1.1", ""),
        complete_gs.replace("certified_rank=4", "certified_rank=2", 1),
        complete_gs.replace("boundary_rank=4", "boundary_rank=5"),
        complete_gs.replace("proof_energy=1.1", "proof_energy=0.7"),
        complete_gs.replace("[HYBRID-LCFO-SYMMETRY]", "[MISSING-HYBRID-LCFO-SYMMETRY]"),
        complete_gs.replace(" occupied=1e-12", ""),
        complete_gs.replace("energy=3e-12", "energy=1e-2"),
        complete_gs.replace("density=4e-12", "density=nan"),
        complete_gs.replace("[HYBRID-RT-BASIS]", "[MISSING-HYBRID-RT-BASIS]"),
        complete_gs.replace(" checkpoint_version=3", ""),
        complete_gs.replace("checkpoint_version=3", "checkpoint_version=2"),
        complete_gs.replace(" construction_rank=7", ""),
        complete_gs.replace(" rt_rank=4", ""),
        complete_gs.replace("construction_rank=7", "construction_rank=3"),
        complete_gs.replace("rt_rank=4", "rt_rank=7"),
        complete_gs.replace("embedding_fingerprint=271828", "embedding_fingerprint=0"),
        complete_gs.replace("operator_covariance=5e-12", "operator_covariance=1e-2"),
        complete_gs.replace("construction_rank=7 certified_rank=4 rt_rank=4",
                            "construction_rank=7 certified_rank=5 rt_rank=4"),
    ]
    for incomplete in incomplete_logs:
        try:
            certify(incomplete, complete_rt)
        except ValueError:
            continue
        raise AssertionError("incomplete or old one-shot evidence was accepted")
    legacy_gs = complete_gs.replace(
        "[HYBRID-LCFO-WINDOW] mode=explicit compatibility=0 delta_e=0.5 homo=-0.2 "
        "requested_cutoff=0.3 requested_rank=3 certified_cutoff=0.8 certified_rank=4 "
        "boundary_rank=4 extension_states=1 extension_energy=0.5 proof_state=1 proof_energy=1.1\n",
        "[HYBRID-SYMMETRY-COMPATIBILITY-WARNING] energy_window=-1 uses dynamic requested-rank selection\n"
        "[HYBRID-LCFO-WINDOW] mode=legacy_dynamic compatibility=1 delta_e=-1 homo=-0.2 "
        "requested_cutoff=0.3 requested_rank=3 certified_cutoff=0.8 certified_rank=4 "
        "boundary_rank=4 extension_states=1 extension_energy=0.5 proof_state=1 proof_energy=1.1\n",
    )
    certify(legacy_gs, complete_rt)
    invalid_legacy_logs = [
        legacy_gs.replace(
            "[HYBRID-SYMMETRY-COMPATIBILITY-WARNING] energy_window=-1 uses dynamic requested-rank selection\n", ""
        ),
        legacy_gs.replace("mode=legacy_dynamic", "mode=explicit"),
        legacy_gs.replace("compatibility=1", "compatibility=0"),
        complete_gs.replace("delta_e=0.5", "delta_e=-1"),
    ]
    for invalid in invalid_legacy_logs:
        try:
            certify(invalid, complete_rt)
        except ValueError:
            continue
        raise AssertionError("incomplete legacy dynamic-rank evidence was accepted")
    invalid_rt_logs = [
        complete_rt.replace("projector_defect=2e-12", "projector_defect=1e-2"),
        complete_rt.replace("lcfo_residual=2e-12", "lcfo_residual=1e-2"),
        complete_rt.replace("metric_defect=2e-12", "metric_defect=-1e-2"),
        complete_rt.replace("metric_defect=2e-12", "metric_defect=nan"),
        complete_rt.replace(" certified_rank=4", ""),
        complete_rt.replace("occupied_rank=2", "occupied_rank=5"),
        complete_rt.replace("state_rank=7", "state_rank=8"),
        complete_rt.replace("metric_rank=7", "metric_rank=8"),
        complete_rt.replace("operator_rank=7", "operator_rank=8"),
        complete_rt.replace("basis_rank=7", "basis_rank=8"),
        complete_rt.replace("coefficient_rows=7", "coefficient_rows=8"),
    ]
    for invalid in invalid_rt_logs:
        try:
            certify(complete_gs, invalid)
        except ValueError:
            continue
        raise AssertionError("invalid measured RT handoff evidence was accepted")


def production_rt_receipt_self_test() -> None:
    state_source = (ROOT / "src/rt/dg/rt_dg_hybrid_initialization_v5.f90").read_text().lower()
    state_type = state_source[
        state_source.index("type,public::s_rt_dg_hybrid_state"):
        state_source.index("end type s_rt_dg_hybrid_state")
    ]
    for field in ("startup_orbital_residual", "startup_metric_defect", "startup_projector_defect"):
        assert field in state_type, f"RT state does not retain measured {field}"
    assert "build_certified_rt_state" not in state_source

    rt_source = (ROOT / "src/rt/main_tddft.f90").read_text().lower()
    marker = rt_source.index("[hybrid-rt-handoff]")
    receipt = rt_source[marker:marker + 1400]
    for field in (
        "certified_rank=", "state_rank=", "metric_rank=", "operator_rank=", "basis_rank=",
        "coefficient_rows=", "occupied_rank=", "lcfo_residual=", "metric_defect=", "projector_defect=",
    ):
        assert field in receipt, f"production RT handoff receipt omits {field}"
    for forbidden in ("operator_symmetry=", "operation_count=", "nonidentity_count="):
        assert forbidden not in receipt, f"fabricated v5 symmetry evidence remains: {forbidden}"


def certify_seed_receipt(log: str, ranks: int, skipped: bool) -> dict[str, int]:
    records = receipt_records(log, "[DG-DC-SEED]")
    if len(records) != 1:
        raise ValueError("expected exactly one canonical DG DC seed receipt")
    record = records[0]
    try:
        mode = record["mode"]
        publication_id = int(record["publication_id"])
        receipt_skipped = {"T": True, "F": False}[record["scf_skipped"].upper()]
        mpi_size = int(record["mpi_size"])
        mapping_fingerprint = int(record["mapping_fingerprint"])
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete or invalid DG DC seed receipt") from error
    if (mode != "auto" or publication_id == 0 or mapping_fingerprint == 0
            or mpi_size != ranks or receipt_skipped != skipped):
        raise ValueError("DG DC seed receipt disagrees with the requested reuse route")
    has_conventional_scf = "DC #SCF =" in log
    if skipped == has_conventional_scf:
        raise ValueError("DG DC seed SCF skip receipt disagrees with conventional SCF output")
    if not skipped:
        lines = log.splitlines()
        seed_line = next(index for index, line in enumerate(lines) if "[DG-DC-SEED]" in line)
        wait_lines = [
            index for index, line in enumerate(lines) if "[DG-DC-SEED-WAIT]" in line
        ]
        wait_records = receipt_records(log, "[DG-DC-SEED-WAIT]")
        if not wait_records:
            raise ValueError("cold DG DC seed publication lacks mixed-charge convergence wait evidence")
        if any(index >= seed_line for index in wait_lines):
            raise ValueError("DG DC seed wait evidence was emitted after the publication receipt")
        try:
            wait_values = [
                (float(wait["error"]), float(wait["tolerance"])) for wait in wait_records
            ]
        except (KeyError, ValueError) as error:
            raise ValueError("incomplete or invalid DG DC seed wait evidence") from error
        if any(
            not math.isfinite(electron_error)
            or not math.isfinite(tolerance)
            or tolerance < 0.0
            or electron_error <= tolerance
            for electron_error, tolerance in wait_values
        ):
            raise ValueError("DG DC seed wait evidence does not prove a blocked convergence exit")
    return {
        "publication_id": publication_id,
        "mapping_fingerprint": mapping_fingerprint,
    }


def certify_seed_workflow(cold_log: str, reused_logs: list[str], ranks: int) -> dict[str, int]:
    if not reused_logs:
        raise ValueError("DC seed workflow has no same-layout reuse evidence")
    cold = certify_seed_receipt(cold_log, ranks, skipped=False)
    for log in reused_logs:
        reused = certify_seed_receipt(log, ranks, skipped=True)
        if reused != cold:
            raise ValueError("DC seed publication or rank-to-fragment mapping changed during reuse")
    return cold


def certify_seed_rejection(log: str, returncode: int, expected_cause: str) -> None:
    if returncode == 0:
        raise ValueError("incompatible DG DC seed exited successfully")
    if "DC #SCF =" in log:
        raise ValueError("incompatible DG DC seed reached conventional SCF")
    if "[DG-DC-SEED]" in log:
        raise ValueError("incompatible DG DC seed emitted an acceptance receipt")
    records = receipt_records(log, "[DG-DC-SEED-ERROR]")
    if len(records) != 1:
        raise ValueError("DG DC seed mismatch did not emit exactly one tagged rejection")
    record = records[0]
    if (record.get("status") != "invalid" or record.get("cause") != expected_cause
            or not record.get("detail")):
        raise ValueError("DG DC seed rejection lacks the expected cause-specific receipt")
    if "ERROR STOP DG DC seed is absent, invalid, or incompatible" not in log:
        raise ValueError("DG DC seed mismatch did not terminate at the pre-DC seed gate")


def render_input(template: str, settings: dict[str, str]) -> str:
    rendered = template
    for name, value in settings.items():
        assignment = re.compile(
            rf"^(?P<prefix>[ \t]*{re.escape(name)}[ \t]*=[ \t]*).*$",
            re.IGNORECASE | re.MULTILINE,
        )
        matches = list(assignment.finditer(rendered))
        if len(matches) != 1:
            raise ValueError(f"expected one Si64 input assignment for {name}, found {len(matches)}")
        rendered = assignment.sub(lambda match: match.group("prefix") + value, rendered, count=1)
    return rendered


def fortran_real_assignment(text: str, name: str) -> float:
    number = r"[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+\-]?\d+)?"
    matches = re.findall(
        rf"^[ \t]*{re.escape(name)}[ \t]*=[ \t]*({number})[ \t]*(?:!.*)?$",
        text,
        re.IGNORECASE | re.MULTILINE,
    )
    if len(matches) != 1:
        raise ValueError(f"expected one real Si64 input assignment for {name}, found {len(matches)}")
    value = float(matches[0].replace("d", "e").replace("D", "E"))
    if not math.isfinite(value):
        raise ValueError(f"Si64 input assignment for {name} is not finite")
    return value


def build_run_environment(parent: dict[str, str]) -> dict[str, str]:
    environment = {
        name: value for name, value in parent.items()
        if not name.startswith("SALMON_")
    }
    environment["OMP_NUM_THREADS"] = "1"
    environment["GFORTRAN_UNBUFFERED_PRECONNECTED"] = "y"
    environment["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    return environment


def certify_zero_field_rt_input(text: str) -> int:
    shape_matches = re.findall(
        r"^[ \t]*ae_shape1[ \t]*=[ \t]*['\"]([^'\"]+)['\"][ \t]*(?:!.*)?$",
        text,
        re.IGNORECASE | re.MULTILINE,
    )
    nt_matches = re.findall(
        r"^[ \t]*nt[ \t]*=[ \t]*([+\-]?\d+)[ \t]*(?:!.*)?$",
        text,
        re.IGNORECASE | re.MULTILINE,
    )
    if len(shape_matches) != 1 or shape_matches[0].strip().lower() != "none":
        raise ValueError("Si64 RT acceptance requires exactly one ae_shape1='none' assignment")
    if len(nt_matches) != 1:
        raise ValueError("Si64 RT acceptance requires exactly one integer nt assignment")
    steps = int(nt_matches[0])
    if steps < 1:
        raise ValueError("Si64 RT acceptance requires a positive nt")
    return steps


def certify_production_gs_receipts(
    log: str,
    expected_pw_cutoff: float,
) -> dict[str, float | int | str]:
    group_matches = re.findall(
        r"\[OW-GS-DIAGNOSTIC\]\s+affine_generator_proof\s+"
        r"group_order=(\d+)\s+generator_count=(\d+)",
        log,
    )
    localization_matches = receipt_records(log, "[HYBRID-WF-LOCALIZATION]")
    cutoff_matches = receipt_records(log, "[HYBRID-PW-CUTOFF]")
    if not group_matches or not localization_matches or not cutoff_matches:
        raise ValueError("missing localization-first reciprocal-action evidence")
    group_order, generator_count = (int(value) for value in group_matches[-1])
    if group_order <= 1 or generator_count <= 0 or generator_count >= group_order:
        raise ValueError("Si64 calculation did not exercise a nonidentity crystal action")
    localization = localization_matches[-1]
    try:
        symmetry_constraint = localization["symmetry_constraint"]
        raw_rank = int(localization["raw_rank"])
        retained_rank = int(localization["retained_rank"])
        localization_iterations = int(localization["iterations"])
        spreads = tuple(float(localization[name]) for name in (
            "spread_min", "spread_max", "spread_mean", "spread_total", "unitarity",
        ))
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete unconstrained WF localization receipt") from error
    if (symmetry_constraint != "off" or raw_rank < 1 or retained_rank != raw_rank
            or localization_iterations < 0 or any(not math.isfinite(value) or value < 0.0 for value in spreads)
            or spreads[-1] > SYMMETRY_TOLERANCE):
        raise ValueError("invalid unconstrained WF localization receipt")
    cutoff = cutoff_matches[-1]
    try:
        requested_cutoff = float(cutoff["requested"])
        effective_cutoff = float(cutoff["effective"])
        shell_added = int(cutoff["shell_added"])
        orbit_added = int(cutoff["orbit_added"])
    except (KeyError, ValueError) as error:
        raise ValueError("incomplete reciprocal PW cutoff receipt") from error
    if (not math.isfinite(expected_pw_cutoff) or expected_pw_cutoff <= 0.0
            or not math.isfinite(requested_cutoff) or not math.isfinite(effective_cutoff)
            or requested_cutoff <= 0.0 or effective_cutoff < 0.0
            or not printed_receipt_close(requested_cutoff, expected_pw_cutoff)
            or shell_added < 0 or orbit_added < 0):
        raise ValueError("invalid reciprocal PW cutoff completion receipt")
    return {
        "group_order": group_order,
        "generator_count": generator_count,
        "wf_symmetry_constraint": symmetry_constraint,
        "wf_raw_rank": raw_rank,
        "wf_localization_iterations": localization_iterations,
        "pw_requested_cutoff": requested_cutoff,
        "pw_effective_cutoff": effective_cutoff,
        "pw_shell_added": shell_added,
        "pw_orbit_added": orbit_added,
    }


def seed_workflow_self_test() -> None:
    fixture = (FIXTURES / "input_hybrid_dg_continuation.in").read_text().lower()
    assert re.search(r"\bdg_dc_seed_mode\s*=\s*['\"]auto['\"]", fixture), (
        "Si64 continuation fixture must cold-publish and then reuse with seed mode auto"
    )
    assert re.search(r"\bdg_dc_seed_directory\s*=\s*['\"]\.\./dc-seed['\"]", fixture), (
        "Si64 continuation fixture must name the runner's shared sibling seed directory"
    )

    cold = (
        "DC #SCF = 1 Total Energy = -1 diff = 1e-3\n"
        "[DG-DC-SEED-WAIT] iteration=85 mixed_electrons=255.9999991 error=8.6e-7 tolerance=1e-8\n"
        "[DG-DC-SEED-WAIT] iteration=86 mixed_electrons=255.9999995 error=5.0e-7 tolerance=1e-8\n"
        "[DG-DC-SEED] mode=auto publication_id=17 scf_skipped=F mpi_size=8 mapping_fingerprint=271828\n"
    )
    reused = (
        "[DG-DC-SEED] mode=auto publication_id=17 scf_skipped=T mpi_size=8 mapping_fingerprint=271828\n"
    )
    evidence = certify_seed_workflow(cold, [reused, reused, reused, reused], 8)
    assert evidence == {"publication_id": 17, "mapping_fingerprint": 271828}
    cold_without_wait = (
        "DC #SCF = 1 Total Energy = -1 diff = 1e-3\n"
        "[DG-DC-SEED] mode=auto publication_id=17 scf_skipped=F mpi_size=8 mapping_fingerprint=271828\n"
    )
    cold_wait_after_receipt = cold.replace(
        "[DG-DC-SEED-WAIT] iteration=86 mixed_electrons=255.9999995 error=5.0e-7 tolerance=1e-8\n",
        "",
    ) + "[DG-DC-SEED-WAIT] iteration=86 mixed_electrons=255.9999995 error=5.0e-7 tolerance=1e-8\n"
    invalid_cold_logs = (
        cold_without_wait,
        cold_wait_after_receipt,
        cold.replace("error=5.0e-7", "error=1e-8"),
        cold.replace("error=5.0e-7", "error=nan"),
        cold.replace("tolerance=1e-8", "tolerance=inf", 1),
    )
    for invalid in invalid_cold_logs:
        try:
            certify_seed_workflow(invalid, [reused], 8)
        except ValueError:
            continue
        raise AssertionError("missing, late, or invalid DG DC seed wait evidence was accepted")
    invalid_reuse = (
        reused + "DC #SCF = 1 Total Energy = -1 diff = 1e-3\n",
        reused.replace("publication_id=17", "publication_id=19"),
        reused.replace("mapping_fingerprint=271828", "mapping_fingerprint=314159"),
        reused.replace("mpi_size=8", "mpi_size=4"),
        reused.replace("scf_skipped=T", "scf_skipped=F"),
        reused.replace("mode=auto", "mode=off"),
        reused.replace(" publication_id=17", ""),
    )
    for invalid in invalid_reuse:
        try:
            certify_seed_workflow(cold, [invalid], 8)
        except ValueError:
            continue
        raise AssertionError("invalid or non-reused DC seed evidence was accepted")

    rejected = (
        "[DG-DC-SEED-ERROR] status=invalid cause=rank_fragment_mapping "
        "detail=rank_fragment_mismatch\n"
        "ERROR STOP DG DC seed is absent, invalid, or incompatible\n"
    )
    certify_seed_rejection(rejected, 1, "rank_fragment_mapping")
    try:
        certify_seed_rejection(
            "ERROR STOP invalid preserved total-system topology for DG DC seed\n",
            1,
            "rank_fragment_mapping",
        )
    except ValueError:
        pass
    else:
        raise AssertionError("irrelevant preserved-topology failure was accepted as a seed rejection")
    for invalid_log, returncode in (
        (rejected, 0),
        (rejected + "DC #SCF = 1\n", 1),
        ("ERROR STOP DG DC seed is absent, invalid, or incompatible\n", 1),
        ("[DG-DC-SEED-ERROR]\nERROR STOP DG DC seed is absent, invalid, or incompatible\n", 1),
        ("other failure\n", 1),
    ):
        try:
            certify_seed_rejection(invalid_log, returncode, "rank_fragment_mapping")
        except ValueError:
            continue
        raise AssertionError("unclear or late DC seed mismatch was accepted")


def rendering_self_test() -> None:
    template = (
        "&dc\n wannier_num_iter=10000\n wannier_pw_cutoff=0.5d0\n"
        " dg_dc_seed_directory='../dc-seed'\n/\n"
    )
    rendered = render_input(template, {
        "wannier_num_iter": "9999",
        "wannier_pw_cutoff": "0.55d0",
        "dg_dc_seed_directory": "'/absolute/dc-seed'",
    })
    for expected in (
        "wannier_num_iter=9999",
        "wannier_pw_cutoff=0.55d0",
        "dg_dc_seed_directory='/absolute/dc-seed'",
    ):
        assert expected in re.sub(r"\s+", "", rendered)
    assert fortran_real_assignment(rendered, "wannier_pw_cutoff") == 0.55
    for invalid_template, setting in ((template, "missing_control"), (template + "wannier_num_iter=4\n", "wannier_num_iter")):
        try:
            render_input(invalid_template, {setting: "1"})
        except ValueError:
            continue
        raise AssertionError("ambiguous or missing Si64 input assignment was accepted")


def rt_output_layout_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="si64-hybrid-rt-layout-") as temporary:
        output = Path(temporary)
        gs_case = output / "certified-gs"
        gs_case.mkdir()
        for name, payload in (
            ("atom.dat", b"atoms"),
            ("Si_rps.dat", b"pseudo"),
            ("hybrid_dg_ground_state.chk", b"checkpoint"),
        ):
            (gs_case / name).write_bytes(payload)
        (gs_case / "variables.log").write_text("gs evidence\n")
        (gs_case / "overlapping_wannier_mlwf.wout").write_text("gs localization evidence\n")
        rt_case, rt_input = prepare_rt_case(output, "certified-rt", gs_case, "&calculation\n/\n")
        assert rt_case != gs_case and rt_input == rt_case / "input-rt.in"
        assert (rt_case / "hybrid_dg_ground_state.chk").read_bytes() == b"checkpoint"
        assert not (rt_case / "overlapping_wannier_mlwf.wout").exists()
        (rt_case / "variables.log").write_text("rt evidence\n")
        assert (gs_case / "variables.log").read_text() == "gs evidence\n"


def zero_field_input_self_test() -> None:
    fixture = (FIXTURES / "input_hybrid_dg_zero_field_rt.in").read_text()
    assert certify_zero_field_rt_input(fixture) == 10
    invalid_inputs = (
        fixture.replace("ae_shape1='none'", "ae_shape1='impulse'"),
        fixture.replace("ae_shape1='none'", ""),
        fixture + "ae_shape1='none'\n",
        fixture.replace("nt=10", "nt=0"),
        fixture.replace("nt=10", "nt=ten"),
    )
    for invalid in invalid_inputs:
        try:
            certify_zero_field_rt_input(invalid)
        except ValueError:
            continue
        raise AssertionError("non-zero-field or invalid RT input was accepted")


def production_gs_receipt_self_test() -> None:
    complete = (
        "[OW-GS-DIAGNOSTIC] affine_generator_proof group_order=1536 generator_count=5\n"
        "[HYBRID-WF-LOCALIZATION] symmetry_constraint=off raw_rank=37 retained_rank=37 "
        "iterations=10 spread_min=1 spread_max=2 spread_mean=1.5 spread_total=55.5 unitarity=1e-13\n"
        "[HYBRID-PW-CUTOFF] requested=0.5 effective=0.4687864262645563 "
        "shell_added=0 orbit_added=0\n"
    )
    evidence = certify_production_gs_receipts(complete, expected_pw_cutoff=0.5)
    assert evidence["group_order"] == 1536 and evidence["generator_count"] == 5
    assert evidence["wf_symmetry_constraint"] == "off"
    invalid_logs = (
        complete.replace("group_order=1536", "group_order=1"),
        complete.replace("generator_count=5", "generator_count=0"),
        complete.replace("symmetry_constraint=off", "symmetry_constraint=on"),
        complete.replace("[HYBRID-PW-CUTOFF]", "[MISSING-HYBRID-PW-CUTOFF]"),
        complete.replace("effective=0.4687864262645563", "effective=-0.1"),
    )
    for invalid in invalid_logs:
        try:
            certify_production_gs_receipts(invalid, expected_pw_cutoff=0.5)
        except ValueError:
            continue
        raise AssertionError("incomplete localization-first reciprocal evidence was accepted")
    try:
        certify_production_gs_receipts(complete, expected_pw_cutoff=0.55)
    except ValueError:
        pass
    else:
        raise AssertionError("PW receipt disagreed with the rendered cutoff but was accepted")


def seed_tree_digest_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="si64-seed-digest-") as temporary:
        seed = Path(temporary) / "seed"
        seed.mkdir()
        (seed / "payload.bin").write_bytes(b"payload")
        baseline = seed_tree_digest(seed)
        assert baseline == seed_tree_digest(seed)
        (seed / "broken-link").symlink_to(seed / "absent")
        try:
            seed_tree_digest(seed)
        except ValueError:
            pass
        else:
            raise AssertionError("broken symbolic link in the committed seed was ignored")
    assert certify_seed_directory_path(Path("/tmp/si64-seed")) == "/tmp/si64-seed"
    for invalid in (
        Path("/tmp/si64-'seed"),
        Path("/tmp/si64-\nseed"),
        Path("/tmp/si64-\tseed"),
        Path("/" + "x" * 241),
    ):
        try:
            certify_seed_directory_path(invalid)
        except ValueError:
            continue
        raise AssertionError("unsafe DG DC seed input path was accepted")


def run_environment_self_test() -> None:
    environment = build_run_environment({
        "PATH": "/bin",
        "SALMON_DG_W90_REPLAY_DIRECTORY": "/outside/replay",
        "SALMON_DG_VARIATIONAL_PAYLOAD_CAPTURE": "/outside/payload",
        "SALMON_POISSON_KERNEL_MODE": "hse_sr",
    })
    assert environment["PATH"] == "/bin"
    assert environment["OMP_NUM_THREADS"] == "1"
    assert environment["GFORTRAN_UNBUFFERED_PRECONNECTED"] == "y"
    assert environment["OMPI_MCA_rmaps_base_oversubscribe"] == "1"
    assert "SALMON_DG_W90_REPLAY_DIRECTORY" not in environment
    assert "SALMON_DG_VARIATIONAL_PAYLOAD_CAPTURE" not in environment
    assert not any(name.startswith("SALMON_") for name in environment)


def prepare_case(output: Path, name: str, input_text: str) -> tuple[Path, Path]:
    case = output / name
    if case.exists():
        raise ValueError(f"case directory must be fresh: {case}")
    case.mkdir(parents=True)
    shutil.copy2(FIXTURES / "atom.dat", case / "atom.dat")
    shutil.copy2(ROOT / "samples/exercise_04_bulkSi_gs/Si_rps.dat", case / "Si_rps.dat")
    input_path = case / "input-gs.in"
    input_path.write_text(input_text)
    return case, input_path


def prepare_rt_case(
    output: Path,
    name: str,
    gs_case: Path,
    input_text: str,
) -> tuple[Path, Path]:
    case = output / name
    if case.exists():
        raise ValueError(f"RT case directory must be fresh: {case}")
    required = ("atom.dat", "Si_rps.dat", "hybrid_dg_ground_state.chk")
    for filename in required:
        source = gs_case / filename
        if source.is_symlink() or not source.is_file():
            raise ValueError(f"certified GS handoff file is missing or unsafe: {source}")
    case.mkdir(parents=True)
    for filename in required:
        shutil.copy2(gs_case / filename, case / filename)
    input_path = case / "input-rt.in"
    input_path.write_text(input_text)
    return case, input_path


def process_session_members(session_id: int) -> set[int] | None:
    try:
        listing = subprocess.run(
            ["ps", "-axo", "pid="],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except OSError:
        return None
    if listing.returncode != 0:
        return None
    members: set[int] = set()
    parsed_any = False
    for line in listing.stdout.splitlines():
        try:
            process_id = int(line.strip())
        except ValueError:
            continue
        parsed_any = True
        try:
            process_session = os.getsid(process_id)
        except (ProcessLookupError, PermissionError):
            continue
        if process_session == session_id:
            members.add(process_id)
    return members if parsed_any else None


def process_group_is_live(group_id: int) -> bool:
    try:
        os.killpg(group_id, 0)
    except ProcessLookupError:
        return False
    return True


def process_session_status(session_id: int) -> tuple[bool, bool]:
    members = process_session_members(session_id)
    if members is not None:
        return bool(members), True
    return process_group_is_live(session_id), False


def process_session_is_live(session_id: int) -> bool:
    return process_session_status(session_id)[0]


def wait_for_process_session_exit(session_id: int, timeout_seconds: float) -> tuple[bool, bool]:
    deadline = time.monotonic() + timeout_seconds
    inspection_verified = True
    while time.monotonic() < deadline:
        live, inspected = process_session_status(session_id)
        inspection_verified = inspection_verified and inspected
        if not live:
            return True, inspection_verified
        time.sleep(0.05)
    live, inspected = process_session_status(session_id)
    return not live, inspection_verified and inspected


def signal_process_session(session_id: int, signal_number: int) -> bool:
    members = process_session_members(session_id)
    if members is not None:
        for process_id in sorted(members - {session_id}) + ([session_id] if session_id in members else []):
            try:
                if os.getsid(process_id) != session_id:
                    continue
                os.kill(process_id, signal_number)
            except (ProcessLookupError, PermissionError):
                pass
        return True
    try:
        os.killpg(session_id, signal_number)
    except ProcessLookupError:
        pass
    return False


def stop_process(
    process: subprocess.Popen[bytes],
    require_session_inspection: bool = False,
) -> tuple[int, bytes]:
    session_id = process.pid
    inspection_verified = signal_process_session(session_id, signal.SIGTERM)
    try:
        remaining, _ = process.communicate(timeout=2.0)
    except subprocess.TimeoutExpired:
        inspection_verified = signal_process_session(session_id, signal.SIGKILL) and inspection_verified
        remaining, _ = process.communicate(timeout=5.0)
    exited, inspected = wait_for_process_session_exit(session_id, 0.5)
    inspection_verified = inspection_verified and inspected
    if not exited:
        inspection_verified = signal_process_session(session_id, signal.SIGKILL) and inspection_verified
        exited, inspected = wait_for_process_session_exit(session_id, 2.0)
        inspection_verified = inspection_verified and inspected
        if not exited:
            members = process_session_members(session_id)
            raise RuntimeError(f"failed to reap launcher session: {members}")
    if require_session_inspection and not inspection_verified:
        raise RuntimeError("launcher session cleanup could not be verified because process inspection failed")
    return process.returncode, remaining or b""


def cleanup_after_failure(
    process: subprocess.Popen[bytes],
    original_error: BaseException,
    record_output: Callable[[bytes], object] | None = None,
    require_session_inspection: bool = False,
) -> None:
    try:
        _, tail = stop_process(process, require_session_inspection=require_session_inspection)
        if record_output is not None:
            record_output(tail)
    except BaseException as cleanup_error:
        original_error.add_note(
            f"launcher cleanup also failed: {type(cleanup_error).__name__}: {cleanup_error}"
        )


def reject_unverified_session_exit(process: subprocess.Popen[bytes]) -> None:
    error = RuntimeError("launcher exit could not be verified because process inspection failed")
    cleanup_after_failure(process, error, require_session_inspection=True)
    raise error


def process_lifecycle_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="si64-hybrid-runner-") as temporary:
        root = Path(temporary)
        separate_group = process_session_members(os.getsid(0)) is not None
        launcher = root / "fake-mpiexec"
        launcher.write_text(
            "#!/usr/bin/env python3\n"
            "import os, subprocess, sys, time\n"
            "from pathlib import Path\n"
            "separate_group = os.environ['FAKE_SEPARATE_GROUP'] == '1'\n"
            "child = subprocess.Popen(\n"
            "    [sys.executable, '-c', \"import signal,time; "
            "signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(5)\"],\n"
            "    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,\n"
            "    preexec_fn=os.setpgrp if separate_group else None)\n"
            "Path(os.environ['FAKE_PARENT_PID']).write_text(str(os.getpid()))\n"
            "Path(os.environ['FAKE_CHILD_PID']).write_text(str(child.pid))\n"
            "if os.environ['FAKE_MODE'] == 'partial_receipt':\n"
            "    os.write(1, b'[DG-DC-')\n"
            "    time.sleep(0.05)\n"
            "    os.write(1, b'SEED] mode=auto publication_id=17 scf_skipped=T mpi_size=1 mapping_finger')\n"
            "    time.sleep(0.20)\n"
            "    Path(os.environ['FAKE_LINE_COMPLETE']).write_text('yes')\n"
            "    os.write(1, b'print=271828\\n')\n"
            "if os.environ['FAKE_MODE'] != 'orphan_exit':\n"
            "    time.sleep(5)\n"
        )
        launcher.chmod(0o755)
        input_path = root / "input.in"
        input_path.write_text("")

        def pid_is_ours(path: Path, session_id: int) -> bool:
            if not path.exists():
                return False
            try:
                process_id = int(path.read_text())
                return os.getsid(process_id) == session_id
            except (ProcessLookupError, PermissionError, ValueError):
                return False

        def require_reaped(*pid_paths: Path) -> None:
            session_id = int(pid_paths[0].read_text())
            deadline = time.monotonic() + 2.0
            while time.monotonic() < deadline and any(pid_is_ours(path, session_id) for path in pid_paths):
                time.sleep(0.01)
            live = [int(path.read_text()) for path in pid_paths if pid_is_ours(path, session_id)]
            assert not live, f"runner left its launcher session alive: pids={live}, session={session_id}"

        partial_case = root / "partial"
        partial_case.mkdir()
        parent_pid = root / "partial-parent.pid"
        child_pid = root / "partial-child.pid"
        line_complete = root / "line-complete"
        partial_env = {
            **os.environ,
            "FAKE_MODE": "partial_receipt",
            "FAKE_PARENT_PID": str(parent_pid),
            "FAKE_CHILD_PID": str(child_pid),
            "FAKE_LINE_COMPLETE": str(line_complete),
            "FAKE_SEPARATE_GROUP": "1" if separate_group else "0",
        }
        log = run_until_seed_receipt(
            str(launcher), Path("/unused"), partial_case, input_path, 1, partial_env, 2.0,
        )
        assert line_complete.exists(), "runner stopped on a partial receipt line"
        certify_seed_receipt(log, 1, skipped=True)
        require_reaped(parent_pid, child_pid)

        timeout_case = root / "timeout"
        timeout_case.mkdir()
        timeout_parent_pid = root / "timeout-parent.pid"
        timeout_child_pid = root / "timeout-child.pid"
        timeout_env = {
            **os.environ,
            "FAKE_MODE": "timeout",
            "FAKE_PARENT_PID": str(timeout_parent_pid),
            "FAKE_CHILD_PID": str(timeout_child_pid),
            "FAKE_LINE_COMPLETE": str(root / "unused-line-complete"),
            "FAKE_SEPARATE_GROUP": "1" if separate_group else "0",
        }
        try:
            run_to_completion(
                str(launcher), Path("/unused"), timeout_case, input_path, "timeout.log", 1,
                timeout_env, 0.2,
            )
        except TimeoutError:
            pass
        else:
            raise AssertionError("fake timeout launcher completed successfully")
        require_reaped(timeout_parent_pid, timeout_child_pid)

        orphan_case = root / "orphan"
        orphan_case.mkdir()
        orphan_parent_pid = root / "orphan-parent.pid"
        orphan_child_pid = root / "orphan-child.pid"
        orphan_env = {
            **os.environ,
            "FAKE_MODE": "orphan_exit",
            "FAKE_PARENT_PID": str(orphan_parent_pid),
            "FAKE_CHILD_PID": str(orphan_child_pid),
            "FAKE_LINE_COMPLETE": str(root / "unused-orphan-line"),
            "FAKE_SEPARATE_GROUP": "1" if separate_group else "0",
        }
        returncode, _ = run_to_completion(
            str(launcher), Path("/unused"), orphan_case, input_path, "orphan.log", 1,
            orphan_env, 2.0,
        )
        assert returncode == 0
        require_reaped(orphan_parent_pid, orphan_child_pid)


def exception_cleanup_self_test() -> None:
    global stop_process
    real_stop_process = stop_process

    def failing_stop_process(_process: object, **_kwargs: object) -> tuple[int, bytes]:
        raise RuntimeError("sentinel cleanup failure")

    stop_process = failing_stop_process  # type: ignore[assignment]
    original = KeyboardInterrupt("sentinel original interrupt")
    try:
        try:
            raise original
        except BaseException as error:
            cleanup_after_failure(None, error)  # type: ignore[arg-type]
            raise
    except KeyboardInterrupt as caught:
        assert caught is original
        assert any("sentinel cleanup failure" in note for note in getattr(caught, "__notes__", ()))
        try:
            reject_unverified_session_exit(None)  # type: ignore[arg-type]
        except RuntimeError as unverified:
            assert str(unverified) == "launcher exit could not be verified because process inspection failed"
            assert any("sentinel cleanup failure" in note for note in getattr(unverified, "__notes__", ()))
        else:
            raise AssertionError("unverified normal launcher exit was accepted")
    finally:
        stop_process = real_stop_process


def run_until_seed_receipt(
    launcher: str,
    binary: Path,
    case: Path,
    input_path: Path,
    ranks: int,
    environment: dict[str, str],
    timeout_seconds: float,
    require_session_inspection: bool = False,
) -> str:
    log_path = case / "gs.log"
    print(f"START {case.name}: run through the canonical DC seed receipt", flush=True)
    with input_path.open("rb") as source, log_path.open("wb") as log:
        process = subprocess.Popen(
            [launcher, "-n", str(ranks), str(binary)],
            cwd=case,
            stdin=source,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=environment,
            bufsize=0,
            start_new_session=True,
        )
        assert process.stdout is not None
        selector = selectors.DefaultSelector()
        selector.register(process.stdout, selectors.EVENT_READ)
        deadline = time.monotonic() + timeout_seconds
        receipt_seen = False
        pending_line = b""
        returncode: int | None = None

        def record_output(chunk: bytes) -> bool:
            nonlocal pending_line
            if not chunk:
                return False
            log.write(chunk)
            log.flush()
            pending_line += chunk
            found = False
            while b"\n" in pending_line:
                line, pending_line = pending_line.split(b"\n", 1)
                if b"[DG-DC-SEED]" in line:
                    found = True
            return found

        try:
            while True:
                remaining_time = deadline - time.monotonic()
                if remaining_time <= 0.0:
                    raise TimeoutError(f"{case.name} timed out before the DG DC seed receipt")
                events = selector.select(timeout=min(1.0, remaining_time))
                if events:
                    chunk = os.read(process.stdout.fileno(), 1024 * 1024)
                    if record_output(chunk):
                        receipt_seen = True
                        returncode, tail = stop_process(
                            process, require_session_inspection=require_session_inspection,
                        )
                        record_output(tail)
                        break
                polled = process.poll()
                if polled is not None:
                    returncode, tail = stop_process(
                        process, require_session_inspection=require_session_inspection,
                    )
                    if record_output(tail):
                        receipt_seen = True
                    break
        except BaseException as error:
            cleanup_after_failure(
                process,
                error,
                record_output,
                require_session_inspection=require_session_inspection,
            )
            raise
        finally:
            selector.close()
    status = {
        "intentional_stop_after_seed_receipt": receipt_seen,
        "returncode": returncode,
        "ranks": ranks,
    }
    (case / "runner-status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    log_text = log_path.read_text(errors="replace")
    if not receipt_seen:
        raise RuntimeError(
            f"{case.name} exited with {returncode} before publishing or accepting the DC seed; "
            f"see {log_path}"
        )
    print(f"PASS {case.name}: canonical DC seed receipt observed", flush=True)
    return log_text


def run_to_completion(
    launcher: str,
    binary: Path,
    case: Path,
    input_path: Path,
    log_name: str,
    ranks: int,
    environment: dict[str, str],
    timeout_seconds: float,
    require_session_inspection: bool = False,
) -> tuple[int, str]:
    log_path = case / log_name
    print(f"START {case.name}/{log_name} on {ranks} ranks", flush=True)
    with input_path.open("rb") as source, log_path.open("wb") as log:
        process = subprocess.Popen(
            [launcher, "-n", str(ranks), str(binary)],
            cwd=case,
            stdin=source,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
            start_new_session=True,
        )
        try:
            returncode = process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired as wait_error:
            error = TimeoutError(f"{case.name}/{log_name} timed out; output is preserved")
            cleanup_after_failure(
                process, error, require_session_inspection=require_session_inspection,
            )
            raise error from wait_error
        except BaseException as error:
            cleanup_after_failure(
                process, error, require_session_inspection=require_session_inspection,
            )
            raise
        session_live, inspected = process_session_status(process.pid)
        if require_session_inspection and not inspected:
            reject_unverified_session_exit(process)
        if session_live:
            stop_process(process, require_session_inspection=require_session_inspection)
    (case / f"{log_name}.status.json").write_text(json.dumps({
        "returncode": returncode,
        "ranks": ranks,
    }, indent=2, sort_keys=True) + "\n")
    text = log_path.read_text(errors="replace")
    print(f"DONE {case.name}/{log_name}: exit={returncode}", flush=True)
    return returncode, text


def seed_tree_digest(directory: Path) -> str:
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError(f"committed DG DC seed directory is missing: {directory}")
    digest = hashlib.sha256()
    entries = list(directory.rglob("*"))
    for path in entries:
        if path.is_symlink():
            raise ValueError(f"DG DC seed must not contain symbolic links: {path}")
        if not path.is_file() and not path.is_dir():
            raise ValueError(f"DG DC seed contains an unsupported entry: {path}")
    files = [path for path in entries if path.is_file()]
    if not files:
        raise ValueError("committed DG DC seed directory is empty")
    digest.update(b"SALMON-DG-DC-SEED-TREE-v1\0")
    for path in sorted(entries, key=lambda item: item.relative_to(directory).as_posix()):
        relative = path.relative_to(directory).as_posix().encode()
        kind = b"F" if path.is_file() else b"D"
        size = path.stat().st_size if path.is_file() else 0
        digest.update(kind)
        digest.update(len(relative).to_bytes(8, "little"))
        digest.update(relative)
        digest.update(size.to_bytes(8, "little"))
        if path.is_file():
            with path.open("rb") as stream:
                while chunk := stream.read(1024 * 1024):
                    digest.update(chunk)
    return digest.hexdigest()


def certify_seed_directory_path(directory: Path) -> str:
    value = str(directory)
    if not directory.is_absolute():
        raise ValueError("DG DC seed directory must be absolute")
    if "'" in value:
        raise ValueError("DG DC seed directory cannot be represented in a quoted Fortran input value")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise ValueError("DG DC seed directory must not contain control characters")
    if len(value) > 240:
        raise ValueError("absolute DG DC seed path does not fit the production input field")
    return value


def require_unchanged_seed(seed_directory: Path, expected_digest: str, case_name: str) -> None:
    actual = seed_tree_digest(seed_directory)
    if actual != expected_digest:
        raise ValueError(f"{case_name} changed the committed conventional-DC seed")


def run(binary: Path, output: Path, ranks: int, timeout_seconds: float) -> None:
    if output.exists():
        raise ValueError(f"output directory must be fresh: {output}")
    if ranks != 8:
        raise ValueError("the certified Si64 rank-to-fragment fixture requires exactly 8 ranks")
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0.0:
        raise ValueError("timeout must be positive and finite")
    launcher = shutil.which("mpiexec")
    if launcher is None:
        raise RuntimeError("mpiexec is required for the Si64 continuation/RT verification")
    if process_session_members(os.getsid(0)) is None:
        raise RuntimeError("process-session inspection is required for safe Si64 MPI cleanup")
    output.mkdir(parents=True)
    seed_directory = (output / "dc-seed").resolve()
    seed_directory_input = certify_seed_directory_path(seed_directory)
    template = (FIXTURES / "input_hybrid_dg_continuation.in").read_text()
    base_input = render_input(template, {
        "dg_dc_seed_directory": f"'{seed_directory_input}'",
    })
    requested_pw_cutoff = fortran_real_assignment(base_input, "wannier_pw_cutoff")
    requested_symmetry_window = fortran_real_assignment(base_input, "dg_hybrid_symmetry_energy_window")
    env = build_run_environment(dict(os.environ))

    cold_case, cold_input = prepare_case(output, "cold-auto", base_input)
    cold_log = run_until_seed_receipt(
        launcher, binary, cold_case, cold_input, ranks, env, timeout_seconds,
        require_session_inspection=True,
    )
    cold_receipt = certify_seed_receipt(cold_log, ranks, skipped=False)
    seed_digest = seed_tree_digest(seed_directory)

    reused_logs: list[str] = []
    reuse_case, reuse_input = prepare_case(output, "reuse-auto", base_input)
    reuse_log = run_until_seed_receipt(
        launcher, binary, reuse_case, reuse_input, ranks, env, timeout_seconds,
        require_session_inspection=True,
    )
    certify_seed_receipt(reuse_log, ranks, skipped=True)
    reused_logs.append(reuse_log)
    require_unchanged_seed(seed_directory, seed_digest, reuse_case.name)

    changed_inputs = (
        ("localization-change", {"wannier_num_iter": "9999"}),
        ("pw-cutoff-change", {"wannier_pw_cutoff": "0.55d0"}),
        ("symmetry-window-change", {"dg_hybrid_symmetry_energy_window": "0.25d0"}),
    )
    for case_name, settings in changed_inputs:
        case, input_path = prepare_case(output, case_name, render_input(base_input, settings))
        log = run_until_seed_receipt(
            launcher, binary, case, input_path, ranks, env, timeout_seconds,
            require_session_inspection=True,
        )
        certify_seed_receipt(log, ranks, skipped=True)
        reused_logs.append(log)
        require_unchanged_seed(seed_directory, seed_digest, case_name)

    rank_count_input = render_input(base_input, {
        "nproc_rgrid_tot": "2,2,1",
    })
    rank_case, rank_input_path = prepare_case(output, "rank-count-rejected", rank_count_input)
    rank_returncode, rank_log = run_to_completion(
        launcher, binary, rank_case, rank_input_path, "gs.log", 4, env,
        min(timeout_seconds, 1800.0),
        require_session_inspection=True,
    )
    certify_seed_rejection(rank_log, rank_returncode, "mpi_rank_count")
    require_unchanged_seed(seed_directory, seed_digest, rank_case.name)

    mapping_input = render_input(base_input, {
        "num_fragment": "2,2,1",
        "num_rgrid_buffer": "6,6,0",
        "nproc_rgrid": "2,1,1",
    })
    mapping_case, mapping_input_path = prepare_case(output, "rank-fragment-rejected", mapping_input)
    mapping_returncode, mapping_log = run_to_completion(
        launcher, binary, mapping_case, mapping_input_path, "gs.log", ranks, env,
        min(timeout_seconds, 1800.0),
        require_session_inspection=True,
    )
    certify_seed_rejection(mapping_log, mapping_returncode, "rank_fragment_mapping")
    require_unchanged_seed(seed_directory, seed_digest, mapping_case.name)

    production_case, production_input = prepare_case(output, "certified-gs-rt", base_input)
    production_returncode, production_gs_log = run_to_completion(
        launcher, binary, production_case, production_input, "gs.log", ranks, env, timeout_seconds,
        require_session_inspection=True,
    )
    if production_returncode:
        raise RuntimeError(
            f"Si64 seed-reused DG continuation failed with exit {production_returncode}; "
            f"see {production_case / 'gs.log'}"
        )
    certify_seed_receipt(production_gs_log, ranks, skipped=True)
    reused_logs.append(production_gs_log)
    require_unchanged_seed(seed_directory, seed_digest, production_case.name)
    route_evidence = certify_production_gs_receipts(
        production_gs_log, expected_pw_cutoff=requested_pw_cutoff,
    )

    rt_input_text = (FIXTURES / "input_hybrid_dg_zero_field_rt.in").read_text()
    rt_steps = certify_zero_field_rt_input(rt_input_text)
    rt_case, rt_input = prepare_rt_case(
        output, "certified-rt", production_case, rt_input_text,
    )
    rt_returncode, rt_log = run_to_completion(
        launcher, binary, rt_case, rt_input, "rt.log", ranks, env, timeout_seconds,
        require_session_inspection=True,
    )
    if rt_returncode:
        raise RuntimeError(
            f"Si64 zero-field certified-space hybrid RT failed with exit {rt_returncode}; "
            f"see {rt_case / 'rt.log'}"
        )
    acceptance = certify(
        production_gs_log,
        rt_log,
        expected_steps=rt_steps,
        expected_delta_e=requested_symmetry_window,
    )
    workflow = certify_seed_workflow(cold_log, reused_logs, ranks)
    if workflow != cold_receipt:
        raise ValueError("final DC seed workflow receipt changed after certification")
    require_unchanged_seed(seed_directory, seed_digest, "zero-field RT")

    summary = {
        "runner_contract_version": 1,
        "ranks": ranks,
        "seed_sha256": seed_digest,
        "seed": workflow,
        "localization_first": route_evidence,
        "certified_rt": acceptance,
        "output_layout": {
            "gs_case": production_case.name,
            "rt_case": rt_case.name,
            "separate_directories": True,
        },
        "rejections": {
            "rank_count_returncode": rank_returncode,
            "rank_fragment_returncode": mapping_returncode,
        },
    }
    (output / "acceptance.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", nargs="?", type=Path)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--timeout-seconds", type=float, default=14400.0)
    parser.add_argument("--parser-only", action="store_true")
    args = parser.parse_args()
    parser_self_test()
    production_rt_receipt_self_test()
    seed_workflow_self_test()
    rendering_self_test()
    rt_output_layout_self_test()
    zero_field_input_self_test()
    production_gs_receipt_self_test()
    seed_tree_digest_self_test()
    run_environment_self_test()
    process_lifecycle_self_test()
    exception_cleanup_self_test()
    if args.parser_only:
        print("PASS Si64 DG continuation RT acceptance parser")
        return
    if args.binary is None or args.output is None:
        parser.error("binary and output are required unless --parser-only is used")
    run(args.binary.resolve(strict=True), args.output.resolve(), args.ranks, args.timeout_seconds)
    print(f"PASS Si64 DG continuation to zero-field RT on {args.ranks} ranks")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run and certify the Si64 DG-continuation to zero-field hybrid-RT route."""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import os
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
RECEIPT_FIELD = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")
HANDOFF_RECEIPT = re.compile(
    r"\[HYBRID-RT-HANDOFF\].*payload_fingerprint=(-?\d+).*"
    r"operator_symmetry=([+\-\d.Ee]+).*projector_symmetry=([+\-\d.Ee]+)"
)
STATIONARITY_RECEIPT = re.compile(
    r"\[HYBRID-RT-STATIONARITY\].*density=([+\-\d.Ee]+).*"
    r"energy=([+\-\d.Ee]+).*projector=([+\-\d.Ee]+).*"
    r"electron=([+\-\d.Ee]+).*h_residual=([+\-\d.Ee]+)"
)


def receipt_records(log: str, prefix: str) -> list[dict[str, str]]:
    return [dict(RECEIPT_FIELD.findall(line)) for line in log.splitlines() if prefix in line]


def printed_receipt_close(value: float, expected: float, *scale_values: float) -> bool:
    scale = max(1.0, abs(value), abs(expected), *(abs(item) for item in scale_values))
    return abs(value - expected) <= 2.0e-8 * scale


def certify(gs_log: str, rt_log: str) -> dict[str, float | int]:
    retained_matches = receipt_records(gs_log, RETAINED_BASIS_RECEIPT)
    gs_matches = receipt_records(gs_log, GS_RECEIPT)
    window_matches = receipt_records(gs_log, LCFO_WINDOW_RECEIPT)
    symmetry_matches = receipt_records(gs_log, LCFO_SYMMETRY_RECEIPT)
    rt_basis_matches = receipt_records(gs_log, RT_BASIS_RECEIPT)
    handoff_matches = list(HANDOFF_RECEIPT.finditer(rt_log))
    stationarity_matches = list(STATIONARITY_RECEIPT.finditer(rt_log))
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
    if ((proof_state == 1 and proof_energy <= certified_cutoff)
            or (proof_state == 0 and proof_energy != 0.0)):
        raise ValueError("LCFO proof energy does not bracket the certified space")
    if (checkpoint_version != 3 or embedding_fingerprint == 0
            or not math.isfinite(localization_spread) or localization_spread < 0.0
            or not math.isfinite(operator_covariance) or operator_covariance < 0.0
            or operator_covariance > SYMMETRY_TOLERANCE):
        raise ValueError("invalid certified RT-basis publication receipt")
    handoff = handoff_matches[-1]
    rt_fingerprint = int(handoff.group(1))
    if gs_fingerprint == 0 or rt_fingerprint != gs_fingerprint:
        raise ValueError("GS and RT payload identities differ")
    operator_symmetry = float(handoff.group(2))
    projector_symmetry = float(handoff.group(3))
    if not math.isfinite(operator_symmetry) or operator_symmetry < 0.0:
        raise ValueError("invalid zero-field full-operator symmetry diagnostic")
    if (not math.isfinite(projector_symmetry) or projector_symmetry < 0.0
            or projector_symmetry > SYMMETRY_TOLERANCE):
        raise ValueError("invalid zero-field occupied-projector symmetry receipt")
    stationarity = tuple(float(value) for value in stationarity_matches[-1].groups())
    if any(not value >= 0.0 for value in stationarity):
        raise ValueError("invalid zero-field stationarity receipt")
    stationarity_limits = (FINAL_TOLERANCE, FINAL_TOLERANCE, FINAL_TOLERANCE,
                           ELECTRON_TOLERANCE, FINAL_TOLERANCE)
    if any(value > limit for value, limit in zip(stationarity, stationarity_limits)):
        raise ValueError("zero-field stationarity drift exceeds its production tolerance")
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
        "rt_operator_symmetry_defect": operator_symmetry,
        "rt_projector_symmetry_defect": projector_symmetry,
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
        "[HYBRID-RT-HANDOFF] payload_fingerprint=314159 operator_symmetry=2e-2 projector_symmetry=2e-12\n"
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
    between_levels = certify(between_levels_gs, complete_rt)
    assert between_levels["certified_cutoff"] < between_levels["requested_cutoff"]
    assert between_levels["requested_rank"] == between_levels["certified_rank"] == 3
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
        complete_rt.replace("projector_symmetry=2e-12", "projector_symmetry=1e-2"),
        complete_rt.replace("operator_symmetry=2e-2", "operator_symmetry=-1e-2"),
        complete_rt.replace("operator_symmetry=2e-2", "operator_symmetry=nan"),
    ]
    for invalid in invalid_rt_logs:
        try:
            certify(complete_gs, invalid)
        except ValueError:
            continue
        raise AssertionError("invalid RT symmetry evidence was accepted")


def run(binary: Path, output: Path, ranks: int) -> None:
    if output.exists():
        raise ValueError(f"output directory must be fresh: {output}")
    output.mkdir(parents=True)
    shutil.copy2(FIXTURES / "atom.dat", output / "atom.dat")
    shutil.copy2(ROOT / "samples/exercise_04_bulkSi_gs/Si_rps.dat", output / "Si_rps.dat")
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    gs_log = output / "gs.log"
    with (FIXTURES / "input_hybrid_dg_continuation.in").open("rb") as source, gs_log.open("wb") as log:
        result = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(binary)],
                                cwd=output, stdin=source, stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode:
        raise RuntimeError(f"Si64 DG continuation failed with exit {result.returncode}")
    rt_log = output / "rt.log"
    with (FIXTURES / "input_hybrid_dg_zero_field_rt.in").open("rb") as source, rt_log.open("wb") as log:
        result = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(binary)],
                                cwd=output, stdin=source, stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode:
        raise RuntimeError(f"Si64 zero-field hybrid RT failed with exit {result.returncode}")
    certify(gs_log.read_text(errors="replace"), rt_log.read_text(errors="replace"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", nargs="?", type=Path)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--parser-only", action="store_true")
    args = parser.parse_args()
    parser_self_test()
    if args.parser_only:
        print("PASS Si64 DG continuation RT acceptance parser")
        return
    if args.binary is None or args.output is None:
        parser.error("binary and output are required unless --parser-only is used")
    run(args.binary.resolve(strict=True), args.output.resolve(), args.ranks)
    print(f"PASS Si64 DG continuation to zero-field RT on {args.ranks} ranks")


if __name__ == "__main__":
    main()

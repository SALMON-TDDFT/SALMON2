#!/usr/bin/env python3
"""Run and certify the Si64 DG-continuation to zero-field hybrid-RT route."""

from __future__ import annotations

import argparse
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
GS_RECEIPT = re.compile(
    r"\[HYBRID-GS-ACCEPTANCE\].*seed_identity=(\d).*lambda_zero=(\d).*"
    r"lambda_one=(\d).*final_refresh=(\d).*r_h=([+\-\d.Ee]+).*"
    r"r_rho=([+\-\d.Ee]+).*r_t=([+\-\d.Ee]+).*r_s=([+\-\d.Ee]+).*"
    r"electron=([+\-\d.Ee]+).*symmetry=([+\-\d.Ee]+).*"
    r"real_space=([+\-\d.Ee]+).*payload_fingerprint=(-?\d+)"
)
HANDOFF_RECEIPT = re.compile(
    r"\[HYBRID-RT-HANDOFF\].*payload_fingerprint=(-?\d+).*"
    r"operator_symmetry=([+\-\d.Ee]+).*projector_symmetry=([+\-\d.Ee]+)"
)
STATIONARITY_RECEIPT = re.compile(
    r"\[HYBRID-RT-STATIONARITY\].*density=([+\-\d.Ee]+).*"
    r"energy=([+\-\d.Ee]+).*projector=([+\-\d.Ee]+).*"
    r"electron=([+\-\d.Ee]+).*h_residual=([+\-\d.Ee]+)"
)


def certify(gs_log: str, rt_log: str) -> dict[str, float | int]:
    gs_matches = list(GS_RECEIPT.finditer(gs_log))
    handoff_matches = list(HANDOFF_RECEIPT.finditer(rt_log))
    stationarity_matches = list(STATIONARITY_RECEIPT.finditer(rt_log))
    if not gs_matches:
        raise ValueError("missing complete lambda-one DG continuation receipt")
    if not handoff_matches:
        raise ValueError("missing exact GS-to-RT payload handoff receipt")
    if not stationarity_matches:
        raise ValueError("missing zero-field hybrid RT stationarity receipt")
    gs = gs_matches[-1]
    if tuple(int(gs.group(i)) for i in range(1, 5)) != (1, 1, 1, 1):
        raise ValueError("continuation seed, lambda-zero, lambda-one, or final refresh was not accepted")
    residuals = tuple(float(gs.group(i)) for i in range(5, 12))
    if any(not value >= 0.0 for value in residuals):
        raise ValueError("invalid final DG continuation residual")
    limits = (FINAL_TOLERANCE,) * 4 + (ELECTRON_TOLERANCE, SYMMETRY_TOLERANCE, FINAL_TOLERANCE)
    if any(value > limit for value, limit in zip(residuals, limits)):
        raise ValueError("final DG continuation residual exceeds its production tolerance")
    gs_fingerprint = int(gs.group(12))
    handoff = handoff_matches[-1]
    rt_fingerprint = int(handoff.group(1))
    if gs_fingerprint == 0 or rt_fingerprint != gs_fingerprint:
        raise ValueError("GS and RT payload identities differ")
    handoff_defects = (float(handoff.group(2)), float(handoff.group(3)))
    if any(value < 0.0 or value > SYMMETRY_TOLERANCE for value in handoff_defects):
        raise ValueError("invalid zero-field handoff covariance receipt")
    stationarity = tuple(float(value) for value in stationarity_matches[-1].groups())
    if any(not value >= 0.0 for value in stationarity):
        raise ValueError("invalid zero-field stationarity receipt")
    stationarity_limits = (FINAL_TOLERANCE, FINAL_TOLERANCE, FINAL_TOLERANCE,
                           ELECTRON_TOLERANCE, FINAL_TOLERANCE)
    if any(value > limit for value, limit in zip(stationarity, stationarity_limits)):
        raise ValueError("zero-field stationarity drift exceeds its production tolerance")
    return {
        "payload_fingerprint": gs_fingerprint,
        "r_h": residuals[0],
        "r_rho": residuals[1],
        "r_t": residuals[2],
        "r_s": residuals[3],
        "density_drift": stationarity[0],
        "energy_drift": stationarity[1],
        "projector_drift": stationarity[2],
        "electron_drift": stationarity[3],
        "hamiltonian_residual": stationarity[4],
    }


def parser_self_test() -> None:
    complete_gs = (
        "[HYBRID-GS-ACCEPTANCE] seed_identity=1 lambda_zero=1 lambda_one=1 "
        "final_refresh=1 r_h=1e-10 r_rho=2e-10 r_t=3e-10 r_s=4e-10 "
        "electron=0e0 symmetry=5e-11 real_space=6e-10 payload_fingerprint=314159\n"
    )
    complete_rt = (
        "[HYBRID-RT-HANDOFF] payload_fingerprint=314159 operator_symmetry=1e-12 projector_symmetry=2e-12\n"
        "[HYBRID-RT-STATIONARITY] step=1 density=1e-12 energy=2e-12 projector=3e-12 "
        "electron=4e-12 h_residual=5e-12\n"
    )
    certify(complete_gs, complete_rt)
    incomplete_logs = [
        "[OW-GS-DIAGNOSTIC] one_shot_generalized_eigenexa residual=0\n",
        "[OW-GS] fully refreshed DG continuation converged lambda=1\n",
        complete_gs.replace("lambda_zero=1", "lambda_zero=0"),
        complete_gs.replace("r_h=1e-10", "r_h=1e-2"),
    ]
    for incomplete in incomplete_logs:
        try:
            certify(incomplete, complete_rt)
        except ValueError:
            continue
        raise AssertionError("incomplete or old one-shot evidence was accepted")


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

#!/usr/bin/env python3
"""Focused RED/GREEN fixture for harmonic peak-versus-dip classification."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


checker_path = Path(__file__).with_name("check_si64_overlapping_wannier_response_hhg.py")
spec = importlib.util.spec_from_file_location("si64_response_checker", checker_path)
assert spec is not None and spec.loader is not None
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)

runner_path = Path(__file__).with_name("run_si64_overlapping_wannier_response_hhg.py")
runner_spec = importlib.util.spec_from_file_location("si64_response_runner", runner_path)
assert runner_spec is not None and runner_spec.loader is not None
runner = importlib.util.module_from_spec(runner_spec)
runner_spec.loader.exec_module(runner)

assert runner.LASER_CYCLES == 10.0
assert runner.POST_PULSE_CYCLES == 2.0
assert runner.LASER_DT_AU == 2.0
assert runner.LASER_REFERENCE_DT_AU == 1.0
assert runner.HHG_WINDOW == "hann"
runner_source = Path(runner.__file__).read_text()
checker_source = checker_path.read_text()
assert "exact_site_group_order" not in runner_source
assert "local_exact_group_orders" not in checker_source
assert "--displaced-root" not in checker_source
assert runner.laser_sample_count(runner.LASER_DT_AU) * runner.LASER_DT_AU == \
    runner.laser_sample_count(runner.LASER_REFERENCE_DT_AU) * runner.LASER_REFERENCE_DT_AU

spectrum = np.zeros((9, 4))
spectrum[:, 2] = np.arange(9, dtype=float)
spectrum[:, 3] = 1.0
spectrum[1:4, 3] = [4.0, 0.2, 3.0]
spectrum[3:6, 3] = [3.0, 8.0, 2.0]

dip = checker.harmonic_local_morphology(spectrum, 2)
peak = checker.harmonic_local_morphology(spectrum, 4)
assert dip["kind"] == "dip" and dip["neighbor_ratio"] < 0.1
assert peak["kind"] == "peak" and peak["neighbor_ratio"] > 2.0
print("PASS long-pulse harmonic morphology contract")

#!/usr/bin/env python3
import importlib.util
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools" / "expand_conventional_restart.py"
spec = importlib.util.spec_from_file_location("expand_conventional_restart", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

rng = np.random.default_rng(7)
npoint, nold, nnew = 48, 3, 8
hvol = 0.25
occupied, _ = np.linalg.qr(rng.standard_normal((npoint, nold)))
occupied /= np.sqrt(hvol)
expanded = module.expand_orbitals(occupied, nnew, hvol, seed=11)

assert np.array_equal(expanded[:, :nold], occupied)
overlap = hvol * expanded.T @ expanded
assert np.max(np.abs(overlap - np.eye(nnew))) < 1.0e-12

with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "info.bin"
    module.write_info(path, nk=1, nstate=8, iteration=2000, nprocs=8, real_orbital=True)
    assert module.read_info(path) == (1, 8, 2000, 8, True)

print("PASS conventional restart expansion preserves occupied WFs and orthonormalizes added states")

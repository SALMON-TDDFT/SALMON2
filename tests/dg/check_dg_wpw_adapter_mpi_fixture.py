#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
fixture = ROOT / "tests/dg/test_dg_wpw_adapter_mpi.f90"
runner = ROOT / "tests/dg/run_dg_wpw_adapter_mpi_fixture.py"
assert fixture.exists(), "missing end-to-end RT callback MPI fixture"
assert runner.exists(), "missing executable MPI fixture harness"
text = fixture.read_text().lower()
for token in ("apply_h_wpw_callback", "apply_s_wpw_callback", "mpi_allreduce",
              "dense_h", "dense_s", "ww_nl_blocks", "mixed_h", "mixed_s"):
    assert token in text, f"adapter MPI fixture omits {token}"
print("PASS end-to-end RT callback MPI fixture contract")

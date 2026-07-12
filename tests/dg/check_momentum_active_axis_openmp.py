#!/usr/bin/env python3
from pathlib import Path
import re


source = (Path(__file__).resolve().parents[2] / "src/rt/dg/rt_dg_fragment_ops.f90").read_text()
start = source.index("  subroutine apply_momentum_blocks(")
end = source.index("  end subroutine apply_momentum_blocks", start)
body = source[start:end]

assert "firstprivate(active_dir_count, active_dirs)" in body
assert re.search(r"!\$omp&\s+private\([^\n]*active_dir_count", body) is None
print("PASS momentum active-axis metadata survives OpenMP entry")

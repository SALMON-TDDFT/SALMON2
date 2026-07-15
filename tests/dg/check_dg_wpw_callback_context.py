#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
common = (ROOT / "src/common/dg_wpw_matrix_free_operator.f90").read_text().lower()
adapter = (ROOT / "src/rt/dg/rt_dg_wpw_matrix_free_adapter.f90").read_text().lower()

for token in ("class(*), intent(inout) :: context", "xw_owned", "xp_owned", "yw_owned", "yp_owned"):
    assert token in common, f"common callback cannot carry split WPW state: {token}"
for token in ("s_rt_dg_wpw_operator_context", "bind_rt_dg_wpw_operator_context",
              "context_bound", "apply_h_wpw_callback", "apply_s_wpw_callback", "select type"):
    assert token in adapter, f"RT adapter lacks usable callback context bridge: {token}"
assert "if (.not. ctx%context_bound)" in adapter, "callbacks accept a context that was not bound"
callback_part = common.split("subroutine global_gram_batch", 1)[0]
assert "n_local_rows" not in callback_part, "flattened callback contract still hides W/P ownership"
print("PASS common callback has a usable RT WPW context bridge")

#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
text = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text().lower()
assert "use dg_wpw_checkpoint" in text, "GS production path does not use the shared WPW checkpoint"
start = text.index("subroutine publish_wpw_production_checkpoint")
end = text.index("end subroutine publish_wpw_production_checkpoint", start)
body = text[start:end]
for token in (
    "wpw_qw", "wpw_qp", "wpw_eigenvalues", "wpw_occupations",
    "wpw_volume_accumulator%potentials", "bounded_operator%layout_fingerprint",
    "bounded_operator%operator_epoch", "bounded_operator%ww_h", "bounded_operator%ww_s",
    "bounded_operator%wp_h", "bounded_operator%wp_s", "bounded_operator%pp_h",
    "bounded_operator%pp_s", "write_dg_wpw_checkpoint",
):
    assert token in body, f"checkpoint publication omits {token}"
converged = re.search(r"if\s*\(\s*state_info==0\.and\.wpw_scf_state%converged\s*\)\s*then(.*?)endif", text, re.S)
assert converged and "publish_wpw_production_checkpoint" in converged.group(1), \
    "SCF convergence returns without publishing the WPW checkpoint"
assert "if(state_info/=0)return" in converged.group(1).replace(" ", ""), \
    "checkpoint publication failure does not revoke SCF success"
print("PASS converged WPW SCF publishes complete shared checkpoint before cleanup")

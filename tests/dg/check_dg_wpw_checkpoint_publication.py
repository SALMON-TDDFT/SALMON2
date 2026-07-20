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
for token in (
    "fixed_h_mode", "density_carrying_fragment_seed", "metric_residual",
    "captured_norm", "projection_rank", "projection_charge",
    "final_interface_lambda", "tolerance_profile", "frozen_layout_fingerprint",
    "frozen_ww_provenance_fingerprint",
):
    assert token in body, f"fixed-H checkpoint metadata omits {token}"
assert "validate_wpw_frozen_h_state" in body, "checkpoint publication lacks exact frozen-state gate"
assert "interface_lambda" in body and "1d0" in body, "checkpoint publication is not gated at lambda one"
assert "status='delete'" in body, "publication attempt must revoke any stale manifest before rank writes"
assert "read_dg_wpw_checkpoint" in body, "rank checkpoint must pass read-back before manifest publication"
assert "expected_fingerprint" in body, "checkpoint read-back must validate layout identity"
assert "wpw_metric_diagnostic_only" in body, "diagnostic metric continuation is not publication-gated"
assert ".not.wpw_metric_diagnostic_only" in body, "diagnostic-only route can publish a checkpoint"
converged = re.search(r"if\s*\(\s*state_info==0\.and\.wpw_scf_state%converged\s*\)\s*then(.*?)endif", text, re.S)
assert converged and "publish_wpw_production_checkpoint" in converged.group(1), \
    "SCF convergence returns without publishing the WPW checkpoint"
assert "if(state_info/=0)return" in converged.group(1).replace(" ", ""), \
    "checkpoint publication failure does not revoke SCF success"
print("PASS converged WPW SCF publishes complete shared checkpoint before cleanup")

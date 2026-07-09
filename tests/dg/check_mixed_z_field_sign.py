from pathlib import Path
import re


src = Path("src/rt/dg/rt_dg_integrator_expdiag.f90").read_text()

match = re.search(r"-\s*E_(?:mid|use)\(\d\)\s*\*\s*dg_frag%mixed_wannier_bpw_z", src)
if match:
    line = src[: match.start()].count("\n") + 1
    raise SystemExit(
        "mixed-Z length-gauge field coupling must use +E*Z to match "
        f"SALMON Vbox += E*r; found negative sign at line {line}"
    )

match = re.search(r"-\s*E_use\(\d\)\s*\*\s*dg_frag%global_wannier_position", src)
if match:
    line = src[: match.start()].count("\n") + 1
    raise SystemExit(
        "global-Wannier length-gauge field coupling must use +E*R to match "
        f"SALMON Vbox += E*r; found negative sign at line {line}"
    )

match = re.search(r"-\s*E_mid\(\d\)\s*\*\s*dg_frag%full_h_seed_xi", src)
if match:
    line = src[: match.start()].count("\n") + 1
    raise SystemExit(
        "full-DG seed length-gauge field coupling must use +E*R to match "
        f"SALMON Vbox += E*r; found negative sign at line {line}"
    )

print("DG length-gauge field coupling signs match +E*R")

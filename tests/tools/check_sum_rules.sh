#!/bin/bash
# check_sum_rules.sh -- verify analytic sum rules from a *_info.data file
# Usage: check_sum_rules.sh <path/to/sysname_info.data>
# Exit code: 0 if all pass, 1 if any fail

set -euo pipefail

INFO="${1:-}"
if [[ -z "$INFO" ]]; then
  echo "ERROR: usage: $0 <path/to/sysname_info.data>"
  exit 1
fi
if [[ ! -f "$INFO" ]]; then
  echo "ERROR: file not found: $INFO"
  exit 1
fi

echo "=== Sum Rule Check: $(basename "$INFO") ==="

python3 - "$INFO" <<'PY'
import re
import sys

fname = sys.argv[1]
lines = open(fname).read()

def extract(pat):
    m = re.search(pat + r'\s*=\s*([-+]?\d+\.?\d*(?:[eE][+-]?\d+)?)', lines)
    if not m:
        print(f"  NOT FOUND: {pat}")
        return None
    return float(m.group(1))

tol_ha = 1e-6
tol_rat = 1e-3
results = []

val_kin_vir = extract(r'Tr\(kin\)\*V \+ 2E_kin')
if val_kin_vir is not None:
    ok = abs(val_kin_vir) < tol_ha
    status = "PASS" if ok else "FAIL"
    print(f"  Tr(kin)*V + 2*E_kin = {val_kin_vir:+.2e} Ha  [{status}]")
    results.append(ok)
else:
    print("  Tr(kin)*V + 2*E_kin: line not found in info.data -- skipping")

val_har_vir = extract(r'Tr\(har\)\*V - E_h')
if val_har_vir is not None:
    ok = abs(val_har_vir) < tol_ha
    status = "PASS" if ok else "FAIL"
    print(f"  Tr(har)*V - E_h     = {val_har_vir:+.2e} Ha  [{status}]")
    results.append(ok)
else:
    print("  Tr(har)*V - E_h: line not found -- skipping")

p_lr_grad = extract(r'P_loc_lr_grad \[GPa\]')
p_lr_diag = extract(r'P_loc_lr_diag \[GPa\]')
if p_lr_grad is not None and p_lr_diag is not None and abs(p_lr_grad) > 1.0:
    ratio = p_lr_diag / p_lr_grad
    ok = abs(ratio - (-1.5)) < tol_rat
    status = "PASS" if ok else "FAIL"
    print(f"  P_lr_diag/P_lr_grad = {ratio:+.6f}  (expect -1.5) [{status}]")
    results.append(ok)
else:
    print("  P_lr ratio: lines not found or P_lr_grad too small -- skipping")

if not results:
    print("  No checks ran (is yn_out_stress='y' in the run?)")
    sys.exit(1)

all_pass = all(results)
n_pass = sum(results)
print(f"=== {n_pass}/{len(results)} sum rules PASS ===")
sys.exit(0 if all_pass else 1)
PY

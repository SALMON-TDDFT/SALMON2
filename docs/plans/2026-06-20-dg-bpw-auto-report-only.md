# DG-BPW Auto Report-Only Diagnostics

## Purpose

`dg_bpw_auto` is the first step toward automatic DG-Wannier+BPW basis selection.  The first-stage implementation is deliberately conservative: it diagnoses the current BPW complement but does not select a new propagation basis.

In the first-stage implementation, `dg_bpw_auto='y'` does not change the propagation basis; it only generates diagnostic reports and recommendation labels.

## Inputs

The relevant propagation namelist inputs are:

```text
dg_bpw_auto = 'y' / 'n'
dg_bpw_auto_accuracy = 1.0d-3
dg_bpw_auto_min_n = -1
dg_bpw_auto_max_n = -1
dg_bpw_auto_report = 'y' / 'n'
```

When `dg_bpw_auto='n'`, SALMON follows the existing manual BPW path.  When `dg_bpw_auto='y'`, SALMON still follows the same propagation path and writes `dg_bpw_auto_report.dat` on the reporting rank.

## Report Contents

The report includes:

- candidate BPW count
- complete shell selection and selected `G2`/energy maximum
- selected raw and Lowdin-filtered BPW counts
- singular value minimum, maximum, and condition number
- projectability boundary diagnostics
- hermiticity of the mixed position matrix
- commutator diagnostic `C_comm_sum`
- f-sum proxy values
- warnings and a machine-readable `SUMMARY` line

The `SUMMARY` line is intended for scripts:

```text
SUMMARY candidate_n=... selected_raw_n=... selected_perp_n=... shell_ecut=... sigma_min=... sigma_max=... cond=... proj_gap=... warnings=... recommendation=...
```

## Recommendation Labels

The report may label a basis as:

```text
recommended
acceptable
risky
rejected
```

`recommended` means numerically healthy in the current diagnostic criteria, not automatically selected for production propagation.

`rejected` is reserved for severe numerical issues such as zero selected BPW states, extreme condition number, vanishing singular values, or explicit min/max count violations.  Small projectability gaps are treated as `risky`, not as fatal.

## Regression Check

The lightweight regression check is:

```bash
tools/check_dg_bpw_auto_report_only.sh CASE_DIR SALMON_BIN
```

It verifies only the first-stage contract:

```text
auto=n/y polarization outputs are identical
auto=y writes dg_bpw_auto_report.dat
the report contains SUMMARY
the CSV helper emits a recommendation column
```

It intentionally does not check spectra or f-sum convergence, since those are physics diagnostics and are less suitable for stable CI.
